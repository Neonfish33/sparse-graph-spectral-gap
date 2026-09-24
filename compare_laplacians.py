"""Нормализованный vs комбинаторный лапласиан: скейлинг lambda2(giant) по n.

Проверяем разницу: 2008 (Samukhin et al.) даёт для КОМБИНАТОРНОГО лапласиана
lambda2(N) ~ (ln N)^-2 при q_m<=2. У нас — НОРМАЛИЗОВАННЫЙ. Считаем оба на одних
и тех же сохранённых графах (results_v2/graphs, pl_dmin2) и сравниваем форму
спада: степенной n^-alpha vs A (ln n)^-beta.

Запуск:
    python compare_laplacians.py
"""

from __future__ import annotations

import os
import numpy as np
import networkx as nx
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from src.generators import load_graph
from src.diagnostics import smallest_eigenpairs

GRAPHS = os.path.join("results_v2", "graphs")
GAMMAS = [2.1, 2.5, 3.0]
NS = [100, 200, 500, 1000, 3000, 10000]
REPS = 10
DMIN = 2


def giant(G):
    return G.subgraph(max(nx.connected_components(G), key=len)).copy()


def norm_l2(G):
    vals, _, _ = smallest_eigenpairs(G, k=2)
    return float(max(vals[1], 0.0)) if len(vals) > 1 else float("nan")


def comb_l2(G):
    nodes = list(G.nodes())
    if len(nodes) < 3:
        return float("nan")
    A = nx.to_scipy_sparse_array(G, nodelist=nodes, dtype=float, format="csr")
    deg = np.asarray(A.sum(axis=1)).ravel()
    L = diags(deg) - A
    try:
        vals = np.sort(eigsh(L, k=2, which="SA", return_eigenvectors=False))
        return float(max(vals[1], 0.0))
    except Exception:
        return float("nan")


def fit_power(ns, ys):
    ns, ys = np.asarray(ns, float), np.asarray(ys, float)
    m = np.isfinite(ys) & (ys > 0)
    if m.sum() < 3:
        return float("nan"), float("nan")
    lx, ly = np.log(ns[m]), np.log(ys[m])
    s, b = np.polyfit(lx, ly, 1)
    r2 = 1 - np.sum((ly - (s * lx + b)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return -s, r2


def fit_logpow(ns, ys):
    ns, ys = np.asarray(ns, float), np.asarray(ys, float)
    m = np.isfinite(ys) & (ys > 0)
    if m.sum() < 3:
        return float("nan"), float("nan")
    lx, ly = np.log(np.log(ns[m])), np.log(ys[m])
    s, b = np.polyfit(lx, ly, 1)
    r2 = 1 - np.sum((ly - (s * lx + b)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return -s, r2


def main():
    by = {}
    total = len(GAMMAS) * len(NS) * REPS
    done = 0
    for g in GAMMAS:
        for n in NS:
            for r in range(REPS):
                p = os.path.join(GRAPHS, f"pl_dmin2_g{g}_n{n}_d{DMIN}_r{r}.npz")
                if not os.path.exists(p):
                    continue
                G, _ = load_graph(p)
                gc = giant(G)
                by.setdefault((g, n), []).append((norm_l2(gc), comb_l2(gc)))
                done += 1
                if done % 30 == 0:
                    print(f"[cmp] {done}/{total}", flush=True)

    print("\n=== lambda2(giant): нормализованный vs комбинаторный ===")
    print(f"{'gamma':>6} {'n':>7} {'norm':>10} {'comb':>10}")
    for g in GAMMAS:
        for n in NS:
            v = by.get((g, n))
            if not v:
                continue
            nm = np.nanmean([x[0] for x in v])
            cm = np.nanmean([x[1] for x in v])
            print(f"{g:>6} {n:>7} {nm:>10.5f} {cm:>10.5f}")

    print("\n=== фиты (n=100..10^4) ===")
    print(f"{'gamma':>6} {'obs':>5} {'alpha(power)':>13} {'R2':>6} | "
          f"{'beta(log)':>10} {'R2':>6}")
    for g in GAMMAS:
        ns = [n for n in NS if (g, n) in by]
        norm = [np.nanmean([x[0] for x in by[(g, n)]]) for n in ns]
        comb = [np.nanmean([x[1] for x in by[(g, n)]]) for n in ns]
        for name, ys in (("norm", norm), ("comb", comb)):
            a, r2a = fit_power(ns, ys)
            b, r2b = fit_logpow(ns, ys)
            print(f"{g:>6} {name:>5} {a:>13.4f} {r2a:>6.3f} | {b:>10.4f} {r2b:>6.3f}")


if __name__ == "__main__":
    main()
