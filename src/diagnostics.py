"""Спектральные и комбинаторные диагностики (разделы 6, 7, 13 конспекта).

Считаем для графа G и его подграфов (giant component, 2-core):
  * lambda2, lambda3 нормализованного лапласиана L = I - D^{-1/2} A D^{-1/2};
  * phi_periph, phi_mod, rho_periph, rho_mod, rho_star через spectral sweep cut;
  * Delta = phi_mod / phi_periph, rho_2 (спектральный cut).

Трюк для устойчивого вычисления малых собственных значений:
  L и M = I - L/2 = (I + D^{-1/2} A D^{-1/2}) / 2 имеют спектры
  lambda in [0,2]  <->  mu = 1 - lambda/2 in [0,1].
  Малые lambda = большие mu, поэтому берём eigsh(M, which='LA').
"""

from __future__ import annotations

import math

import networkx as nx
import numpy as np

try:
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh

    HAVE_SCIPY = True
except ImportError:  # pragma: no cover - fallback для окружений без scipy
    HAVE_SCIPY = False

DENSE_THRESHOLD = 1500
INF = float("inf")


def _degrees_and_adj(G: nx.Graph):
    nodes = list(G.nodes())
    index = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    deg = np.zeros(n, dtype=float)
    adj = [[] for _ in range(n)]
    for u, v in G.edges():
        iu, iv = index[u], index[v]
        if iu == iv:
            continue
        deg[iu] += 1.0
        deg[iv] += 1.0
        adj[iu].append(iv)
        adj[iv].append(iu)
    return deg, adj


def _normalized_laplacian_dense(G: nx.Graph):
    nodes = list(G.nodes())
    n = len(nodes)
    A = nx.to_numpy_array(G, nodelist=nodes, weight=None, dtype=float)
    deg = A.sum(axis=1)
    inv = np.zeros_like(deg)
    nz = deg > 0
    inv[nz] = 1.0 / np.sqrt(deg[nz])
    D = np.diag(inv)
    L = np.eye(n) - D @ A @ D
    return L, deg


def ipr(v) -> float:
    """Inverse participation ratio: IPR(v) = sum v_i^4 / (sum v_i^2)^2.

    Делокализованный вектор на n узлах: IPR ~ 1/n. Локализованный на k
    узлах: IPR ~ 1/k. У собственных векторов из eigh/eigsh sum v_i^2 = 1.
    """
    v = np.asarray(v, dtype=float)
    s2 = float(np.sum(v ** 2))
    if s2 <= 0.0:
        return float("nan")
    return float(np.sum(v ** 4) / s2 ** 2)


def ipr_pi(v, degrees) -> float:
    """pi-взвешенный IPR: локализация относительно стационарной меры pi_i = d_i/vol.

    Считается IPR распределения w_i = pi_i v_i^2 / sum_j pi_j v_j^2:
        IPR_pi(v) = sum_i w_i^2.
    Для делокализованного вектора даёт sum_i pi_i^2 (то же, что IPR при
    однородных степенях), для локализованного на узле --- 1.
    """
    v = np.asarray(v, dtype=float)
    d = np.asarray(degrees, dtype=float)
    vol = float(d.sum())
    if vol <= 0.0 or v.size == 0:
        return float("nan")
    w = (d / vol) * v ** 2
    s = float(w.sum())
    if s <= 0.0:
        return float("nan")
    return float(np.sum(w ** 2) / s ** 2)


def smallest_eigenpairs(G: nx.Graph, k: int = 3):
    """Возвращает (eigenvalues по возрастанию, eigenvectors, degrees)."""
    nodes = list(G.nodes())
    n = len(nodes)
    if n == 0:
        return np.array([]), np.empty((0, 0)), np.array([])

    if n <= DENSE_THRESHOLD or not HAVE_SCIPY:
        L, deg = _normalized_laplacian_dense(G)
        vals, vecs = np.linalg.eigh(L)
        kk = min(k, n)
        return vals[:kk], vecs[:, :kk], deg

    A = nx.to_scipy_sparse_array(G, nodelist=nodes, weight=None, dtype=float, format="csr")
    deg = np.asarray(A.sum(axis=1)).ravel()
    inv = np.zeros_like(deg)
    nz = deg > 0
    inv[nz] = 1.0 / np.sqrt(deg[nz])
    D = sp.diags(inv)
    M = 0.5 * (sp.identity(n, format="csr") + D @ A @ D)
    kk = min(k, n - 1)
    if kk < 1:
        return np.array([0.0]), np.ones((n, 1)), deg
    mu, vecs = eigsh(M, k=kk, which="LA")
    lam = 2.0 * (1.0 - mu)
    order = np.argsort(lam)
    return lam[order], vecs[:, order], deg


def _scan(order, degrees, adj, vol_V, eps):
    """Один проход sweep cut. Возвращает лучшие (phi, rho) по двум бакетам."""
    n = len(degrees)
    in_S = np.zeros(n, dtype=bool)
    vol_S = 0.0
    cut = 0
    best = {
        "periph": (INF, None),
        "mod": (INF, None),
        "global": (INF, None),
    }
    for v in order:
        v = int(v)
        in_S[v] = True
        vol_S += degrees[v]
        for u in adj[v]:
            if in_S[u]:
                cut -= 1
            else:
                cut += 1
        if vol_S <= 0:
            continue
        phi = cut / vol_S
        # Ограничиваемся vol(S) <= vol(V)/2: для S с vol > vol(V)/2
        # дополнение имеет тот же conductance, а S = V даёт тривиальный phi = 0.
        if vol_S <= 0.5 * vol_V:
            if phi < best["global"][0]:
                best["global"] = (phi, vol_S / vol_V)
            if vol_S <= eps * vol_V:
                if phi < best["periph"][0]:
                    best["periph"] = (phi, vol_S / vol_V)
            elif phi < best["mod"][0]:
                best["mod"] = (phi, vol_S / vol_V)
    return best


def sweep_curve(degrees, adj, f2):
    """Полный профиль conductance phi(rho) вдоль sweep cut.

    Возвращает dict с двумя направлениями:
      "asc"  -- (rho, phi) для возрастающего порядка f2,
      "desc" -- (rho, phi) для убывающего порядка f2.
    Здесь rho = vol(S_k)/vol(V), phi = |cut(S_k)|/vol(S_k). Последняя точка
    rho = 1 даёт тривиальный phi = 0. Профиль используется для визуализации
    и для определения режимов post hoc вместо фиксированного eps.
    """
    degrees = np.asarray(degrees, dtype=float)
    vol_V = float(degrees.sum())
    n = len(degrees)
    empty = (np.array([]), np.array([]))
    if n == 0 or vol_V <= 0:
        return {"asc": empty, "desc": empty}

    order = np.argsort(f2)
    out = {}
    for name, scan_order in (("asc", order), ("desc", order[::-1])):
        in_S = np.zeros(n, dtype=bool)
        vol_S = 0.0
        cut = 0
        rhos, phis = [], []
        for v in scan_order:
            v = int(v)
            in_S[v] = True
            vol_S += degrees[v]
            for u in adj[v]:
                if in_S[u]:
                    cut -= 1
                else:
                    cut += 1
            if vol_S > 0:
                rhos.append(vol_S / vol_V)
                phis.append(cut / vol_S)
        out[name] = (np.array(rhos), np.array(phis))
    return out


def graph_sweep_curve(G: nx.Graph):
    """φ-профиль второго собственного вектора для giant component графа G.

    Возвращает результат sweep_curve (dict с 'asc'/'desc') или None, если
    граф пуст/без рёбер. Для несвязного G профиль считается по giant
    component --- там, где λ₂ содержательно.
    """
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return None
    giant = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    vals, vecs, _ = smallest_eigenpairs(giant, k=3)
    if vecs.shape[1] < 2:
        return None
    f2 = vecs[:, 1]
    deg, adj = _degrees_and_adj(giant)
    return sweep_curve(deg, adj, f2)


def sweep_diagnostics(G: nx.Graph, f2, degrees, adj, eps: float = 0.05):
    """phi_periph/phi_mod/rho_*/Delta по второму собственному вектору f2.

    Сканируем в обоих направлениях (низкая и высокая сторона f2), чтобы не
    потерять малый cut на стороне с высокими значениями.
    """
    n = len(degrees)
    vol_V = degrees.sum()
    if n == 0 or vol_V <= 0:
        return _empty_sweep()

    order = np.argsort(f2)
    best = {"periph": (INF, None), "mod": (INF, None), "global": (INF, None)}
    for scan_order in (order, order[::-1]):
        b = _scan(scan_order, degrees, adj, vol_V, eps)
        for key in best:
            if b[key][0] < best[key][0]:
                best[key] = b[key]

    phi_periph, rho_periph = best["periph"]
    phi_mod, rho_mod = best["mod"]
    phi_star, rho_star = best["global"]

    # log Delta = log(phi_mod / phi_periph); nan, если conductance не найден
    # (нулевой или бесконечный), чтобы не путать с содержательными значениями.
    if (
        phi_periph <= 0.0
        or phi_mod <= 0.0
        or math.isinf(phi_periph)
        or math.isinf(phi_mod)
    ):
        log_Delta = float("nan")
    else:
        log_Delta = math.log(phi_mod) - math.log(phi_periph)

    if len(f2) == 0 or vol_V <= 0:
        rho_2 = float("nan")
    else:
        vol_pos = float(degrees[f2 > 0].sum())
        vol_neg = vol_V - vol_pos
        rho_2 = min(vol_pos, vol_neg) / vol_V

    return {
        "phi_periph": phi_periph,
        "phi_mod": phi_mod,
        "phi_star": phi_star,
        "rho_periph": rho_periph,
        "rho_mod": rho_mod,
        "rho_star": rho_star,
        "rho_2": rho_2,
        "log_Delta": log_Delta,
    }


def _empty_sweep():
    return {
        "phi_periph": INF,
        "phi_mod": INF,
        "phi_star": INF,
        "rho_periph": None,
        "rho_mod": None,
        "rho_star": None,
        "rho_2": None,
        "log_Delta": INF,
    }


def _graph_metrics(G: nx.Graph, prefix: str, eps: float = 0.05, k: int = 3):
    n = G.number_of_nodes()
    m = G.number_of_edges()
    out = {
        f"{prefix}_n": n,
        f"{prefix}_m": m,
        f"{prefix}_lambda2": float("nan"),
        f"{prefix}_lambda3": float("nan"),
        f"{prefix}_gap23": float("nan"),
        f"{prefix}_phi_periph": float("nan"),
        f"{prefix}_phi_mod": float("nan"),
        f"{prefix}_phi_star": float("nan"),
        f"{prefix}_rho_star": float("nan"),
        f"{prefix}_rho_2": float("nan"),
        f"{prefix}_log_Delta": float("nan"),
        f"{prefix}_ipr2": float("nan"),
        f"{prefix}_ipr2_pi": float("nan"),
        f"{prefix}_pr2": float("nan"),
    }
    if n == 0 or m == 0:
        return out

    # Для несвязного графа lambda_2 = 0 тривиально (вырожденное собственное
    # значение с кратностью = числу компонент), а spectral sweep cut теряет
    # смысл: f_2 локализован на одной компоненте. Диагностики щели для такого
    # графа неинтерпретируемы --- их считаем отдельно для giant/core.
    if nx.number_connected_components(G) > 1:
        out[f"{prefix}_lambda2"] = 0.0
        out[f"{prefix}_lambda3"] = 0.0
        out[f"{prefix}_gap23"] = 0.0
        return out

    vals, vecs, deg = smallest_eigenpairs(G, k=k)
    lam2 = max(0.0, float(vals[1])) if len(vals) > 1 else 0.0
    lam3 = max(0.0, float(vals[2])) if len(vals) > 2 else float("nan")
    out[f"{prefix}_lambda2"] = lam2
    out[f"{prefix}_lambda3"] = lam3
    out[f"{prefix}_gap23"] = lam3 - lam2 if lam3 == lam3 else float("nan")

    f2 = vecs[:, 1] if vecs.shape[1] > 1 else vecs[:, 0]
    deg_full, adj = _degrees_and_adj(G)
    sw = sweep_diagnostics(G, f2, deg_full, adj, eps=eps)
    out[f"{prefix}_phi_periph"] = sw["phi_periph"]
    out[f"{prefix}_phi_mod"] = sw["phi_mod"]
    out[f"{prefix}_phi_star"] = sw["phi_star"]
    out[f"{prefix}_rho_star"] = sw["rho_star"] if sw["rho_star"] is not None else float("nan")
    out[f"{prefix}_rho_2"] = sw["rho_2"] if sw["rho_2"] is not None else float("nan")
    out[f"{prefix}_log_Delta"] = sw["log_Delta"]

    # Локализация второго собственного вектора (делокализованный -> 1/n).
    ipr2 = ipr(f2)
    out[f"{prefix}_ipr2"] = ipr2
    out[f"{prefix}_ipr2_pi"] = ipr_pi(f2, deg_full)
    out[f"{prefix}_pr2"] = (1.0 / ipr2) if ipr2 and ipr2 > 0 else float("nan")
    return out


def analyze_graph(G: nx.Graph, eps: float = 0.05, k: int = 3) -> dict:
    """Полный набор диагностик для G, giant component и 2-core."""
    row: dict = {}
    n = G.number_of_nodes()
    m = G.number_of_edges()
    degs = [d for _, d in G.degree()] if n else []
    row["n"] = n
    row["m"] = m
    row["avg_deg"] = (2.0 * m / n) if n else float("nan")
    row["d_max"] = max(degs) if degs else 0
    row["d_min_actual"] = min(degs) if degs else 0
    row["n_components"] = nx.number_connected_components(G) if n else 0
    row["N1"] = sum(1 for d in degs if d == 1)
    row["N2"] = sum(1 for d in degs if d == 2)

    row.update(_graph_metrics(G, "full", eps=eps, k=k))

    if n:
        giant_nodes = max(nx.connected_components(G), key=len)
        giant = G.subgraph(giant_nodes).copy()
        row["giant_frac"] = len(giant_nodes) / n
    else:
        giant = nx.Graph()
        row["giant_frac"] = float("nan")
    row.update(_graph_metrics(giant, "giant", eps=eps, k=k))

    core = nx.k_core(G, 2) if n else nx.Graph()
    row["core_frac"] = core.number_of_nodes() / n if n else float("nan")
    row.update(_graph_metrics(core, "core", eps=eps, k=k))

    return row


if __name__ == "__main__":
    import sys

    g = nx.path_graph(200)
    print(analyze_graph(g, k=3))
