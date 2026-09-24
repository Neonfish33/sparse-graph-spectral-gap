"""Matching library: minimal spectral metrics for ER vs power-law 2-core.

Read-only w.r.t. the rest of the project (imports src.* but does not modify it).
Writes only inside results_agents/matcher/.
"""
from __future__ import annotations

import os
import sys

import networkx as nx
import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from src.diagnostics import smallest_eigenpairs
from src.generators import power_law_graph


def gen_er(n: int, c: float, seed: int) -> nx.Graph:
    p = c / max(n - 1, 1)
    return nx.fast_gnp_random_graph(n, p, seed=seed)


def gen_pl(n: int, gamma: float, d_max: int, seed: int) -> nx.Graph:
    return power_law_graph(n, gamma, d_min=2, d_max=d_max, seed=seed)


def _lam2(G: nx.Graph) -> float:
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return float("nan")
    if nx.number_connected_components(G) > 1:
        return 0.0
    vals, _, _ = smallest_eigenpairs(G, k=3)
    if len(vals) > 1:
        return max(0.0, float(vals[1]))
    return 0.0


def metrics(G: nx.Graph) -> dict:
    n = G.number_of_nodes()
    m = G.number_of_edges()
    degs = [d for _, d in G.degree()] if n else []
    avg = (2.0 * m / n) if n else float("nan")
    N2 = sum(1 for d in degs if d == 2)
    P2 = (N2 / n) if n else float("nan")
    out = {
        "n": n,
        "m": m,
        "avg_deg": avg,
        "N2": N2,
        "P2": P2,
        "q": (2.0 * P2 / avg) if avg else float("nan"),
        "giant_frac": float("nan"),
        "giant_n": 0,
        "giant_lambda2": float("nan"),
        "core_frac": float("nan"),
        "core_n": 0,
        "core_m": 0,
        "d_bar_core": float("nan"),
        "core_lambda2": float("nan"),
        "core_N2": 0,
        "core_P2": float("nan"),
        "q_core": float("nan"),
    }
    if n == 0:
        return out

    giant_nodes = max(nx.connected_components(G), key=len)
    giant = G.subgraph(giant_nodes).copy()
    out["giant_frac"] = len(giant_nodes) / n
    out["giant_n"] = len(giant_nodes)
    out["giant_lambda2"] = _lam2(giant)

    core = nx.k_core(G, 2)
    cn = core.number_of_nodes()
    cm = core.number_of_edges()
    out["core_frac"] = cn / n
    out["core_n"] = cn
    out["core_m"] = cm
    if cn:
        out["d_bar_core"] = 2.0 * cm / cn
        cdegs = [d for _, d in core.degree()]
        cN2 = sum(1 for d in cdegs if d == 2)
        out["core_N2"] = cN2
        out["core_P2"] = cN2 / cn
        out["q_core"] = (2.0 * (cN2 / cn) / out["d_bar_core"]) if out["d_bar_core"] else float("nan")
    out["core_lambda2"] = _lam2(core)
    return out


def fit_alpha(ns, ys):
    ns = np.asarray(ns, float)
    ys = np.asarray(ys, float)
    mask = np.isfinite(ns) & np.isfinite(ys) & (ns > 0) & (ys > 0)
    if mask.sum() < 3:
        return float("nan"), float("nan")
    lx, ly = np.log(ns[mask]), np.log(ys[mask])
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return -slope, r2
