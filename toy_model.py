"""Toy model: объяснение показателя alpha в lambda2(giant) ~ n^(-alpha).

Гипотеза H1: если интегрированная плотность малых собственных значений
нормализованного лапласиана ведёт себя как F(lambda) ~ lambda^beta, то
lambda2 ~ n^(-1/beta), т.е. alpha = 1/beta.

Скрипт для каждой модели семейства распределений степеней:
  1. строит граф (min degree 2 для всех, кроме ER);
  2. считает lambda2 giant;
  3. считает K наименьших собственных чисел giant, фитирует beta, alpha_pred=1/beta;
  4. собирает статистику длин цепочек степени 2.

Запуск:
    python toy_model.py --quick
    python toy_model.py --n 10000 --K 150 --out results_toy_model
"""

from __future__ import annotations

import argparse
import csv
import math
import os

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

from src.generators import configuration_graph

# известные alpha (по core/giant) из основных прогонов --- для сверки
KNOWN_ALPHA = {
    "unif_D3": 0.066, "unif_D4": 0.070, "unif_D5": 0.061, "unif_D6": 0.081,
    "unif_D8": 0.069, "unif_D12": 0.076, "unif_D20": 0.064,
    "pl_g2.5_D3": 0.096, "pl_g2.5_D4": 0.101, "pl_g2.5_D5": 0.155,
    "pl_g2.5_D8": 0.146, "pl_g2.5_D15": 0.143, "pl_g2.5_D50": 0.164,
    "two23_p0.6": 0.072,
    "er_c3": 0.089,
}

DEFAULT_MODELS = [
    {"name": "unif_D4", "kind": "unif", "D": 4},
    {"name": "unif_D8", "kind": "unif", "D": 8},
    {"name": "unif_D20", "kind": "unif", "D": 20},
    {"name": "pl_g2.5_D5", "kind": "pl", "gamma": 2.5, "D": 5},
    {"name": "pl_g2.5_D15", "kind": "pl", "gamma": 2.5, "D": 15},
    {"name": "pl_g2.5_D50", "kind": "pl", "gamma": 2.5, "D": 50},
    {"name": "two23_p0.6", "kind": "two23", "p2": 0.6},
    {"name": "er_c3", "kind": "er", "c": 3},
]

QUICK_MODELS = [
    {"name": "unif_D8", "kind": "unif", "D": 8},
    {"name": "pl_g2.5_D50", "kind": "pl", "gamma": 2.5, "D": 50},
    {"name": "er_c3", "kind": "er", "c": 3},
]

FIELDS = [
    "model", "kind", "n", "m", "avg_deg", "p2", "q", "giant_n", "lambda2",
    "n_eigs", "beta", "alpha_pred", "alpha_known",
    "n_chains", "chain_mean", "chain_max", "chain_p_gt10",
]


def make_degrees(spec, n, rng):
    if spec["kind"] == "unif":
        return rng.choice(np.arange(2, spec["D"] + 1), size=n).astype(np.int64)
    if spec["kind"] == "pl":
        ks = np.arange(2, spec["D"] + 1)
        w = ks.astype(float) ** (-spec["gamma"])
        return rng.choice(ks, size=n, p=w / w.sum()).astype(np.int64)
    if spec["kind"] == "two23":
        p2 = spec["p2"]
        return rng.choice([2, 3], size=n, p=[p2, 1 - p2]).astype(np.int64)
    if spec["kind"] == "expo":
        ks = np.arange(2, spec.get("Dcap", 200) + 1)
        w = np.exp(-spec["theta"] * ks)
        return rng.choice(ks, size=n, p=w / w.sum()).astype(np.int64)
    raise ValueError(spec["kind"])


def make_graph(spec, n, seed):
    if spec["kind"] == "er":
        return nx.fast_gnp_random_graph(n, spec["c"] / max(n - 1, 1), seed=seed)
    rng = np.random.default_rng(seed)
    deg = make_degrees(spec, n, rng)
    return configuration_graph(deg, seed=seed)


def giant(G):
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()


def low_eigvals(G, k):
    nodes = list(G.nodes())
    n = len(nodes)
    if n <= 1:
        return np.array([])
    A = nx.to_scipy_sparse_array(G, nodelist=nodes, dtype=float, format="csr")
    deg = np.asarray(A.sum(axis=1)).ravel()
    inv = np.zeros_like(deg)
    nz = deg > 0
    inv[nz] = 1.0 / np.sqrt(deg[nz])
    M = 0.5 * (sp.identity(n, format="csr") + sp.diags(inv) @ A @ sp.diags(inv))
    kk = min(k, n - 1)
    if kk < 1:
        return np.array([])
    mu, _ = eigsh(M, k=kk, which="LA")
    return np.sort(2.0 * (1.0 - mu))


def density_beta(eigs, n, lam_min=1e-9):
    lam = np.sort(eigs[eigs > lam_min])
    if len(lam) < 10:
        return float("nan")
    i = np.arange(1, len(lam) + 1)
    lx, ly = np.log(lam), np.log(i / n)
    return float(np.polyfit(lx, ly, 1)[0])


def chain_lengths(G):
    d2 = [v for v in G.nodes() if G.degree(v) == 2]
    sub = G.subgraph(d2)
    return [len(c) for c in nx.connected_components(sub)]


def analyze(spec, n, seed, K):
    g = make_graph(spec, n, seed)
    gc = giant(g)
    nv = g.number_of_nodes()
    deg = [d for _, d in g.degree()]
    dbar = float(np.mean(deg)) if deg else float("nan")
    p2 = (sum(1 for d in deg if d == 2) / nv) if nv else float("nan")
    eigs = low_eigvals(gc, K)
    lam2 = float(eigs[1]) if len(eigs) > 1 else float("nan")
    beta = density_beta(eigs, gc.number_of_nodes())
    cl = chain_lengths(gc)
    row = {
        "model": spec["name"], "kind": spec["kind"], "n": nv, "m": g.number_of_edges(),
        "avg_deg": dbar, "p2": p2, "q": (2 * p2 / dbar) if dbar else float("nan"),
        "giant_n": gc.number_of_nodes(), "lambda2": lam2,
        "n_eigs": len(eigs), "beta": beta,
        "alpha_pred": (1.0 / beta) if beta and beta == beta and beta > 0 else float("nan"),
        "alpha_known": KNOWN_ALPHA.get(spec["name"], ""),
        "n_chains": len(cl),
        "chain_mean": float(np.mean(cl)) if cl else 0.0,
        "chain_max": int(max(cl)) if cl else 0,
        "chain_p_gt10": (sum(1 for x in cl if x > 10) / len(cl)) if cl else 0.0,
    }
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--K", type=int, default=150)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default="results_toy_model")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    models = QUICK_MODELS if args.quick else DEFAULT_MODELS
    n = 2000 if args.quick else args.n
    K = min(args.K, max(10, (n // 2) - 1))
    reps = 2 if args.quick else args.reps

    rows = []
    for spec in models:
        for rep in range(reps):
            rows.append(analyze(spec, n, seed=1000 + 7 * rep, K=K))
    path = os.path.join(args.out, "toy_metrics.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    print(f"n={n} K={K} reps={reps} -> {path}\n")
    hdr = f"{'model':>14} {'dbar':>6} {'p2':>6} {'q':>6} {'beta':>7} {'a_pred':>7} {'a_known':>8} {'chain_max':>9} {'chain_p>10':>10}"
    print(hdr)
    print("-" * len(hdr))
    by = {}
    for r in rows:
        by.setdefault(r["model"], []).append(r)
    for name in [m["name"] for m in models]:
        rs = by[name]
        f = lambda k: float(np.mean([x[k] for x in rs]))
        kl = rs[0]["alpha_known"]
        print(f"{name:>14} {f('avg_deg'):>6.2f} {f('p2'):>6.3f} {f('q'):>6.3f} "
              f"{f('beta'):>7.3f} {f('alpha_pred'):>7.3f} {str(kl):>8} "
              f"{int(f('chain_max')):>9} {f('chain_p_gt10'):>10.3f}")
    print("\nИнтерпретация: если alpha_pred (=1/beta) совпадает с a_known --- H1 верна.")


if __name__ == "__main__":
    main()
