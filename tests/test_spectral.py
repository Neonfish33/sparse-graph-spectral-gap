"""Тесты генераторов и спектральных диагностик.

Запуск:
    python tests/test_spectral.py        # как скрипт
    pytest tests/test_spectral.py        # через pytest
"""

from __future__ import annotations

import os
import sys

import networkx as nx
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.diagnostics import _degrees_and_adj, _scan, analyze_graph
from src.generators import (
    ensemble_configs,
    generate,
    planted_power_law_graph,
    power_law_degrees,
)


def _brute_phi(degrees, adj, S):
    Sset = set(S)
    vol = sum(degrees[v] for v in S)
    cut = sum(1 for v in S for u in adj[v] if u not in Sset)
    return cut / vol if vol else float("inf")


def test_power_law_degrees_basic():
    deg = power_law_degrees(1000, gamma=2.5, d_min=1, seed=42)
    assert deg.sum() % 2 == 0, "сумма степеней должна быть чётной"
    assert deg.min() >= 1
    assert deg.max() < 1000


def test_sweep_cut_matches_bruteforce():
    """Регрессия: инкрементальный подсчёт cut должен совпадать с brute force."""
    for G in (nx.path_graph(10), nx.gnp_random_graph(60, 0.1, seed=1)):
        deg, adj = _degrees_and_adj(G)
        f2 = np.random.default_rng(0).standard_normal(len(deg))
        order = np.argsort(f2)

        vol_V = deg.sum()
        for scan_order in (order, order[::-1]):
            cands = []
            for k in range(len(scan_order)):
                S = [int(v) for v in scan_order[: k + 1]]
                if sum(deg[v] for v in S) <= 0.5 * vol_V:
                    cands.append(_brute_phi(deg, adj, S))
            true_min = min(cands)
            best = _scan(scan_order, deg, adj, vol_V, eps=0.05)
            assert abs(best["global"][0] - true_min) < 1e-9, (
                best["global"][0],
                true_min,
            )


def test_disconnected_lambda2_is_zero():
    G = nx.disjoint_union(nx.path_graph(20), nx.path_graph(20))
    row = analyze_graph(G, k=3)
    assert row["full_lambda2"] == 0.0
    assert row["full_lambda3"] == 0.0


def test_cheeger_upper_bound():
    """lambda_2 <= 2 phi для giant component.

    Sweep cut даёт конкретный разрез, поэтому его conductance phi_star есть
    верхняя оценка истинного минимума phi(G) <= phi_star, откуда
    lambda_2 <= 2 phi(G) <= 2 phi_star. Обратная граница Cheeger относится к
    истинному phi(G) и по sweep-оценке не проверяется.
    """
    checked = 0
    for spec in ensemble_configs([600], [2.5], d_mins=[1, 2, 3]):
        for rep in range(3):
            g, _ = generate(spec, seed=100 + rep)
            if g.number_of_nodes() == 0:
                continue
            row = analyze_graph(g, k=3)
            lg = row["giant_lambda2"]
            phi = row["giant_phi_star"]
            if lg == lg and phi == phi and phi > 0 and not np.isinf(phi):
                assert lg <= 2 * phi + 1e-6
                checked += 1
    assert checked > 0


def test_planted_mu_control():
    """Доля межблочных концов рёбер должна быть близка к mu."""
    for mu in (0.05, 0.1, 0.2):
        g, labels = planted_power_law_graph(2000, 2.5, d_min=1, mu=mu, seed=11)
        total = 2 * g.number_of_edges()
        cross = sum(2 for u, v in g.edges() if labels[u] != labels[v])
        assert abs(cross / total - mu) < 0.05, (mu, cross / total)


def test_log_delta_is_nan_or_finite():
    for spec in ensemble_configs([700], [2.5], d_mins=[1, 2, 3]):
        for rep in range(2):
            g, _ = generate(spec, seed=200 + rep)
            if g.number_of_nodes() == 0:
                continue
            row = analyze_graph(g, k=3)
            v = row["full_log_Delta"]
            assert v != v or abs(v) < 50


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
