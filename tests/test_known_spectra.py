"""Тесты корректности λ₂ на графах с известным спектром + IPR и φ-профиль.

Известные собственные значения нормализованного лапласиана:
  * K_n          : λ₂ = n/(n-1)
  * C_n (цикл)   : λ₂ = 1 - cos(2π/n)
  * P_n (путь)   : λ₂ = 1 - cos(π/(n-1))
  * K_{a,b}      : λ₂ = 1  (a, b >= 2)
  * дизъюнктное объединение : λ₂ = 0

Запуск:
    python tests/test_known_spectra.py
    pytest tests/test_known_spectra.py
"""

from __future__ import annotations

import os
import sys

import networkx as nx
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.diagnostics import (
    DENSE_THRESHOLD,
    _degrees_and_adj,
    ipr,
    ipr_pi,
    smallest_eigenpairs,
    sweep_curve,
)


def _lambda2(G):
    vals, _, _ = smallest_eigenpairs(G, k=3)
    return float(vals[1]) if len(vals) > 1 else 0.0


def _all_eigs(G):
    vals, _, _ = smallest_eigenpairs(G, k=G.number_of_nodes())
    return np.sort(vals)


def test_complete_graph():
    for n in (5, 20, 100):
        assert abs(_lambda2(nx.complete_graph(n)) - n / (n - 1)) < 1e-9


def test_cycle():
    for n in (10, 50, 200):
        expected = 1.0 - np.cos(2.0 * np.pi / n)
        assert abs(_lambda2(nx.cycle_graph(n)) - expected) < 1e-9


def test_path():
    for n in (3, 10, 50):
        expected = 1.0 - np.cos(np.pi / (n - 1))
        assert abs(_lambda2(nx.path_graph(n)) - expected) < 1e-9


def test_complete_bipartite():
    for a, b in ((2, 3), (5, 5), (20, 30)):
        assert abs(_lambda2(nx.complete_bipartite_graph(a, b)) - 1.0) < 1e-9


def test_disjoint_union_lambda2_zero():
    G = nx.disjoint_union(nx.path_graph(20), nx.path_graph(20))
    vals, _, _ = smallest_eigenpairs(G, k=3)
    assert abs(vals[0]) < 1e-12
    assert abs(vals[1]) < 1e-12


def test_spectrum_in_range():
    for G in (nx.cycle_graph(30), nx.path_graph(30), nx.complete_graph(15)):
        vals = _all_eigs(G)
        assert vals.min() > -1e-9
        assert vals.max() < 2.0 + 1e-9


def test_sparse_path_large_bipartite():
    """n > DENSE_THRESHOLD: проверяем sparse eigsh через трюк M = I - L/2."""
    n = DENSE_THRESHOLD + 500
    G = nx.complete_bipartite_graph(n // 2, n // 2)
    assert G.number_of_nodes() > DENSE_THRESHOLD
    assert abs(_lambda2(G) - 1.0) < 1e-6


def test_sparse_path_cycle_small_gap():
    """Большой цикл: очень малый зазор λ₂ ≈ 2π²/n², проверяем относительно."""
    n = 2000
    G = nx.cycle_graph(n)
    expected = 1.0 - np.cos(2.0 * np.pi / n)
    got = _lambda2(G)
    assert abs(got - expected) / expected < 0.05


def test_ipr_delocalized():
    n = 100
    v = np.ones(n) / np.sqrt(n)
    assert abs(ipr(v) - 1.0 / n) < 1e-12


def test_ipr_localized():
    v = np.zeros(50)
    v[7] = 1.0
    assert abs(ipr(v) - 1.0) < 1e-12
    assert abs(ipr_pi(v, np.ones(50)) - 1.0) < 1e-12


def test_ipr_pi_uniform_degrees_matches_ipr():
    n = 100
    v = np.ones(n) / np.sqrt(n)
    degrees = np.full(n, 3.0)
    assert abs(ipr_pi(v, degrees) - ipr(v)) < 1e-12


def test_ipr_localization_ordering():
    """v₂ у графа с висячей цепочкой локализован сильнее, чем у цикла.

    Эталон -- цикл (невырожденное λ₂, делокализованный вектор, IPR ~ 1/n).
    K_n как эталон не годится: там λ₂ вырождено с кратностью n-1, и вектор
    из eigenspace произволен.
    """
    g = nx.disjoint_union(nx.complete_graph(150), nx.path_graph(60))
    g.add_edge(0, 150)
    _, vecs, _ = smallest_eigenpairs(g, k=3)
    v2 = vecs[:, 1]
    assert float(np.sum(v2[150:] ** 2)) > 0.9  # масса на висячей цепочке
    _, vecs_c, _ = smallest_eigenpairs(nx.cycle_graph(210), k=3)
    v2_c = vecs_c[:, 1]
    assert ipr(v2) > ipr(v2_c)


def test_phi_curve_shape():
    G = nx.cycle_graph(60)
    deg, adj = _degrees_and_adj(G)
    _, vecs, _ = smallest_eigenpairs(G, k=3)
    curves = sweep_curve(deg, adj, vecs[:, 1])
    for side in ("asc", "desc"):
        rho, phi = curves[side]
        assert len(rho) == len(phi) == G.number_of_nodes()
        assert rho[0] > 0
        assert abs(rho[-1] - 1.0) < 1e-12
        assert abs(phi[-1]) < 1e-12  # S = V даёт нулевой cut
        assert np.all(phi >= -1e-12)


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
