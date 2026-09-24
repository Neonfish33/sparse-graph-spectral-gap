"""Independent agent: chains of degree-2 internal vertices vs lambda2.

Read-only on inputs. Writes only into results_agents/chains/.

Definition of a chain (maximal path whose internal vertices have degree 2):
  On the giant component H of a graph, let D2 = {v : deg_H(v) == 2}.
  The induced subgraph H[D2] decomposes into connected components; each
  component is a path (or, degenerately, a cycle).  A path component with
  k degree-2 vertices corresponds to the maximal path
        x0 - a1 - a2 - ... - ak - x1
  whose internal vertices a_i have degree 2 and whose endpoints x0,x1 are
  NOT degree 2 (branch deg>=3 or leaf deg==1).  We use as chain length
  L = k + 1   (number of EDGES of the full maximal path), and also record k.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import networkx as nx
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from src.generators import load_graph  # noqa: E402
from src.diagnostics import smallest_eigenpairs  # noqa: E402

GRAPHS = os.path.join(ROOT, "results_v2", "graphs")
DIAG = os.path.join(ROOT, "results_merged", "diagnostics.csv")
OUT = os.path.join(ROOT, "results_agents", "chains")
os.makedirs(OUT, exist_ok=True)


# --------------------------------------------------------------------------
# chain extraction
# --------------------------------------------------------------------------
def extract_chains(G: nx.Graph, giant_only: bool = True):
    """Return (stats_dict, lengths_L ndarray, lengths_k ndarray, bridge_flags)."""
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return _empty_chain_stats(), np.array([]), np.array([]), np.array([])

    if giant_only:
        H = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    else:
        H = G

    deg = dict(H.degree())
    D2 = [v for v, d in deg.items() if d == 2]
    sub = H.subgraph(D2)
    bridges = set(map(frozenset, nx.bridges(H))) if H.number_of_edges() else set()

    lens_L, lens_k, bridge_flags = [], [], []
    for comp in nx.connected_components(sub):
        comp = list(comp)
        cnt = len(comp)
        e_int = sub.subgraph(comp).number_of_edges()
        # path (cnt-1 internal edges) vs cycle (cnt internal edges)
        is_cycle = (e_int == cnt)
        # attachment edges from comp to non-D2 vertices
        attach = []
        cset = set(comp)
        for u in comp:
            for w in H.neighbors(u):
                if w not in cset:
                    attach.append(frozenset((u, w)))
        # every edge of the maximal path
        internal = [frozenset(e) for e in sub.subgraph(comp).edges()]
        chain_edges = internal + attach
        if is_cycle:
            # degenerate 2-regular giant: treat component as one cyclic chain
            lens_k.append(cnt)
            lens_L.append(cnt)
        else:
            k = cnt
            lens_k.append(k)
            lens_L.append(k + 1)
        all_bridge = len(chain_edges) > 0 and all(e in bridges for e in chain_edges)
        bridge_flags.append(bool(all_bridge))

    lens_L = np.asarray(lens_L, dtype=float)
    lens_k = np.asarray(lens_k, dtype=float)
    bridge_flags = np.asarray(bridge_flags, dtype=bool)
    stats = _summarize(lens_L, lens_k, bridge_flags)
    return stats, lens_L, lens_k, bridge_flags


def _empty_chain_stats():
    return {
        "n_chain": 0, "mean_L": float("nan"), "median_L": float("nan"),
        "p90_L": float("nan"), "p95_L": float("nan"), "p99_L": float("nan"),
        "Lmax": float("nan"), "mean_k": float("nan"), "kmax": float("nan"),
        "mean_L_weighted": float("nan"),
        "frac_bridge_chains": float("nan"), "frac_2ecc_chains": float("nan"),
        "frac_bridge_lenweighted": float("nan"),
        "n_deg2": 0, "giant_n": 0,
    }


def _summarize(lens_L, lens_k, bridge_flags):
    n = len(lens_L)
    if n == 0:
        return _empty_chain_stats()
    w = lens_L / lens_L.sum()
    return {
        "n_chain": int(n),
        "mean_L": float(np.mean(lens_L)),
        "median_L": float(np.median(lens_L)),
        "p90_L": float(np.percentile(lens_L, 90)),
        "p95_L": float(np.percentile(lens_L, 95)),
        "p99_L": float(np.percentile(lens_L, 99)),
        "Lmax": float(np.max(lens_L)),
        "mean_k": float(np.mean(lens_k)),
        "kmax": float(np.max(lens_k)),
        "mean_L_weighted": float(np.sum(w * lens_L)),
        "frac_bridge_chains": float(np.mean(bridge_flags)),
        "frac_2ecc_chains": float(np.mean(~bridge_flags)),
        "frac_bridge_lenweighted": float(np.sum(w * bridge_flags)),
        "n_deg2": int(lens_k.sum()),
        "giant_n": 0,
    }


# --------------------------------------------------------------------------
# jobs
# --------------------------------------------------------------------------
def job_pl(spec):
    (graph_id, lam2, n_target, gamma, dmin, rep, dbar, N2) = spec
    path = os.path.join(GRAPHS, graph_id + ".npz")
    if not os.path.exists(path):
        return None
    G, _ = load_graph(path)
    stats, lens_L, lens_k, bf = extract_chains(G, giant_only=True)
    giant_n = G.subgraph(max(nx.connected_components(G), key=len)).number_of_nodes() \
        if G.number_of_nodes() else 0
    stats["giant_n"] = int(giant_n)
    q = 2.0 * N2 / G.number_of_nodes() / dbar if dbar else float("nan")
    rec = {
        "source": "pl", "graph_id": graph_id, "ensemble": f"pl_dmin{dmin}",
        "gamma": float(gamma), "d_min": int(dmin), "n_target": int(n_target),
        "rep": int(rep), "n_actual": G.number_of_nodes(),
        "lam2_giant": float(lam2), "n_chain": stats["n_chain"],
        "mean_L": stats["mean_L"], "median_L": stats["median_L"],
        "p90_L": stats["p90_L"], "p95_L": stats["p95_L"], "p99_L": stats["p99_L"],
        "Lmax": stats["Lmax"], "mean_k": stats["mean_k"], "kmax": stats["kmax"],
        "mean_L_weighted": stats["mean_L_weighted"],
        "frac_bridge_chains": stats["frac_bridge_chains"],
        "frac_2ecc_chains": stats["frac_2ecc_chains"],
        "frac_bridge_lenweighted": stats["frac_bridge_lenweighted"],
        "n_deg2": stats["n_deg2"], "giant_n": stats["giant_n"],
        "dbar": dbar, "N2": N2, "q": q,
    }
    return rec


# --------------------------------------------------------------------------
# toy / bounded models (regenerated, own seeds)
# --------------------------------------------------------------------------
TOY_MODELS = [
    {"name": "pl_g2.5_D3", "kind": "pl", "gamma": 2.5, "D": 3},
    {"name": "pl_g2.5_D4", "kind": "pl", "gamma": 2.5, "D": 4},
    {"name": "pl_g2.5_D5", "kind": "pl", "gamma": 2.5, "D": 5},
    {"name": "pl_g2.5_D8", "kind": "pl", "gamma": 2.5, "D": 8},
    {"name": "two23_p0.4", "kind": "two23", "p2": 0.4},
    {"name": "two23_p0.6", "kind": "two23", "p2": 0.6},
    {"name": "two23_p0.8", "kind": "two23", "p2": 0.8},
    {"name": "b234_p0.5", "kind": "b234", "p2": 0.5},
    {"name": "b234_p0.7", "kind": "b234", "p2": 0.7},
    {"name": "unif_D3", "kind": "unif", "D": 3},
    {"name": "unif_D4", "kind": "unif", "D": 4},
    {"name": "unif_D5", "kind": "unif", "D": 5},
    {"name": "unif_D8", "kind": "unif", "D": 8},
]


def make_degrees(spec, n, rng):
    if spec["kind"] == "pl":
        ks = np.arange(2, spec["D"] + 1)
        w = ks.astype(float) ** (-spec["gamma"])
        return rng.choice(ks, size=n, p=w / w.sum()).astype(np.int64)
    if spec["kind"] == "two23":
        p2 = spec["p2"]
        return rng.choice([2, 3], size=n, p=[p2, 1 - p2]).astype(np.int64)
    if spec["kind"] == "b234":
        p2 = spec["p2"]
        rem = (1 - p2) / 2
        return rng.choice([2, 3, 4], size=n, p=[p2, rem, rem]).astype(np.int64)
    if spec["kind"] == "unif":
        ks = np.arange(2, spec["D"] + 1)
        return rng.choice(ks, size=n).astype(np.int64)
    raise ValueError(spec["kind"])


def job_toy(task):
    spec, n, rep, seed = task
    from src.generators import configuration_graph
    rng = np.random.default_rng(seed)
    deg = make_degrees(spec, n, rng)
    G = configuration_graph(deg, seed=seed)
    if G.number_of_nodes() == 0:
        return None
    giant = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    vals, _, _ = smallest_eigenpairs(giant, k=2)
    lam2 = float(max(0.0, vals[1])) if len(vals) > 1 else float("nan")
    stats, lens_L, lens_k, bf = extract_chains(G, giant_only=True)
    stats["giant_n"] = giant.number_of_nodes()
    m = G.number_of_edges()
    N2 = sum(1 for _, d in G.degree() if d == 2)
    dbar = 2.0 * m / G.number_of_nodes()
    return {
        "source": "toy", "graph_id": f"{spec['name']}_n{n}_r{rep}",
        "ensemble": spec["name"], "gamma": spec.get("gamma", ""),
        "d_min": 2, "n_target": n, "rep": rep, "n_actual": G.number_of_nodes(),
        "lam2_giant": lam2, "n_chain": stats["n_chain"],
        "mean_L": stats["mean_L"], "median_L": stats["median_L"],
        "p90_L": stats["p90_L"], "p95_L": stats["p95_L"], "p99_L": stats["p99_L"],
        "Lmax": stats["Lmax"], "mean_k": stats["mean_k"], "kmax": stats["kmax"],
        "mean_L_weighted": stats["mean_L_weighted"],
        "frac_bridge_chains": stats["frac_bridge_chains"],
        "frac_2ecc_chains": stats["frac_2ecc_chains"],
        "frac_bridge_lenweighted": stats["frac_bridge_lenweighted"],
        "n_deg2": stats["n_deg2"], "giant_n": stats["giant_n"],
        "dbar": dbar, "N2": N2, "q": 2.0 * N2 / G.number_of_nodes() / dbar if dbar else float("nan"),
    }


FIELDS = [
    "source", "graph_id", "ensemble", "gamma", "d_min", "n_target", "rep",
    "n_actual", "lam2_giant", "n_chain", "mean_L", "median_L", "p90_L",
    "p95_L", "p99_L", "Lmax", "mean_k", "kmax", "mean_L_weighted",
    "frac_bridge_chains", "frac_2ecc_chains", "frac_bridge_lenweighted",
    "n_deg2", "giant_n", "dbar", "N2", "q",
]


def run():
    t0 = time.time()
    # ---- select PL graphs ----
    rows = list(csv.DictReader(open(DIAG, encoding="utf-8")))
    want_ens = {"pl_dmin2", "pl_dmin3", "pl_dmin1"}
    gammas = {"2.5", "3.0"}
    ntargets = {"1000", "3000", "10000"}
    reps_max = 20
    pl_specs = []
    for r in rows:
        if r["ensemble"] not in want_ens:
            continue
        if r["γ"] not in gammas:
            continue
        if r["n_target"] not in ntargets:
            continue
        if int(r["rep"]) >= reps_max:
            continue
        pl_specs.append((
            r["graph_id"], float(r["λ₂_giant"]), int(r["n_target"]),
            float(r["γ"]), int(r["d_min"]), int(r["rep"]),
            float(r["d̄"]), float(r["N₂"]),
        ))
    # dedupe by graph_id
    seen = set()
    pl_specs = [s for s in pl_specs if not (s[0] in seen or seen.add(s[0]))]

    toy_specs = []
    for spec in TOY_MODELS:
        for n in (1000, 10000):
            for rep in range(20):
                # deterministic, process-independent seed (own seeds, NOT the
                # ones used in results_toy/diagnostics.csv -> internally
                # consistent chains+lambda2 per regenerated graph)
                seed = (sum(ord(c) for c in spec["name"]) * 1000 + n + rep * 97) % (2 ** 31)
                toy_specs.append((spec, n, rep, seed))

    print(f"[chains] PL graphs={len(pl_specs)}  toy graphs={len(toy_specs)}", flush=True)
    out_rows = []
    with ProcessPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(job_pl, s) for s in pl_specs]
        futs += [ex.submit(job_toy, s) for s in toy_specs]
        done = 0
        for fut in as_completed(futs):
            rec = fut.result()
            done += 1
            if rec is not None:
                out_rows.append(rec)
            if done % 25 == 0:
                print(f"[chains] done={done}/{len(futs)} {time.time()-t0:6.1f}s", flush=True)

    with open(os.path.join(OUT, "chains_raw.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for rec in out_rows:
            w.writerow({k: rec.get(k, "") for k in FIELDS})
    print(f"[chains] wrote {len(out_rows)} rows -> chains_raw.csv  ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    run()
