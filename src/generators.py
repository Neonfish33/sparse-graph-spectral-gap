"""Генерация графов для четырёх ансамблей проекта.

Ансамбли:
  1. power-law, d_min = 1  -> периферия
  2. power-law, d_min = 2  -> цепочки / 2-core
  3. power-law, d_min = 3  -> экспандероподобный core
  4. power-law + planted communities -> настоящая модульность

Все генераторы возвращают простой граф (без петель и кратных рёбер):
используется стёртая конфигурационная модель.
"""

from __future__ import annotations

import networkx as nx
import numpy as np


def default_dmax(n: int, gamma: float) -> int:
    """Естественный cutoff d_max ~ n^{1/(gamma-1)} (раздел 11 конспекта)."""
    if gamma <= 1:
        raise ValueError("gamma must be > 1")
    dmax = int(round(n ** (1.0 / (gamma - 1.0))))
    return max(2, min(dmax, n - 1))


def power_law_degrees(
    n: int,
    gamma: float,
    d_min: int = 1,
    d_max: int | None = None,
    seed: int | None = None,
) -> np.ndarray:
    """Степени i.i.d. из P(D=k) ~ k^{-gamma}, k in [d_min, d_max].

    Сумма приводится к чётной (необходимое условие для configuration model).
    """
    if d_min < 1:
        raise ValueError("d_min must be >= 1")
    if d_max is None:
        d_max = default_dmax(n, gamma)
    d_max = min(d_max, n - 1)
    if d_min > d_max:
        raise ValueError(f"d_min={d_min} > d_max={d_max}")

    rng = np.random.default_rng(seed)
    ks = np.arange(d_min, d_max + 1)
    weights = ks.astype(float) ** (-gamma)
    probs = weights / weights.sum()
    deg = rng.choice(ks, size=n, p=probs).astype(np.int64)

    if deg.sum() % 2 == 1:
        idx = rng.integers(n)
        if deg[idx] < d_max:
            deg[idx] += 1
        elif deg[idx] > d_min:
            deg[idx] -= 1
        else:
            deg[idx] += 1
    return deg


def _erase_configuration_model(multigraph: nx.MultiGraph) -> nx.Graph:
    """Стираем петли и кратные рёбра (erased configuration model)."""
    g = nx.Graph(multigraph)
    g.remove_edges_from(list(nx.selfloop_edges(g)))
    g.remove_nodes_from(list(nx.isolates(g)))
    return g


def configuration_graph(degrees, seed: int | None = None) -> nx.Graph:
    """Стёртая конфигурационная модель по заданной последовательности степеней."""
    degrees = [int(d) for d in degrees]
    if sum(degrees) % 2 == 1:
        degrees[0] += 1
    mg = nx.configuration_model(degrees, seed=seed)
    mg = nx.MultiGraph(mg)
    return _erase_configuration_model(mg)


def power_law_graph(
    n: int,
    gamma: float,
    d_min: int = 1,
    d_max: int | None = None,
    seed: int | None = None,
) -> nx.Graph:
    """Ансамбли 1-3: чистая power-law конфигурационная модель."""
    deg = power_law_degrees(n, gamma, d_min=d_min, d_max=d_max, seed=seed)
    g = configuration_graph(deg, seed=seed)
    g.graph["d_min"] = d_min
    g.graph["gamma"] = gamma
    return g


def _block_stub_match(
    degrees: np.ndarray,
    groups: np.ndarray,
    mu: float,
    rng: np.random.Generator,
) -> nx.Graph:
    """Stub-matching с блочным предпочтением.

    Каждый stub используется ровно один раз => степени сохраняются точно.
    С вероятностью (1 - mu) партнёр берётся из той же группы, иначе из другой.
    """
    stubs = np.repeat(np.arange(len(degrees)), degrees)
    rng.shuffle(stubs)
    alive = np.ones(len(stubs), dtype=bool)
    pos_by_stub = stubs

    edges = []
    order = rng.permutation(len(stubs))
    for i in order:
        if not alive[i]:
            continue
        u = int(pos_by_stub[i])
        gu = groups[u]
        alive[i] = False
        want_cross = rng.random() < mu
        partner = -1
        for _ in range(64):
            j = int(rng.integers(len(stubs)))
            if not alive[j]:
                continue
            w = int(pos_by_stub[j])
            if w == u:
                continue
            cross = groups[w] != gu
            if cross == want_cross:
                partner = j
                break
        if partner < 0:
            for _ in range(64):
                j = int(rng.integers(len(stubs)))
                if alive[j]:
                    partner = j
                    break
        if partner < 0:
            continue
        alive[partner] = False
        w = int(pos_by_stub[partner])
        if w != u:
            edges.append((u, w))

    g = nx.Graph()
    g.add_nodes_from(range(len(degrees)))
    g.add_edges_from(edges)
    g.remove_edges_from(list(nx.selfloop_edges(g)))
    g.remove_nodes_from(list(nx.isolates(g)))
    return g


def planted_power_law_graph(
    n: int,
    gamma: float,
    d_min: int = 1,
    n_communities: int = 2,
    mu: float = 0.1,
    d_max: int | None = None,
    seed: int | None = None,
) -> tuple[nx.Graph, np.ndarray]:
    """Ансамбль 4: power-law степени + planted блочная структура.

    mu — целевая доля межблочных концов рёбер (fraction of inter-community
    edge endpoints). Чем меньше mu, тем сильнее модульность.
    Возвращает (graph, community_labels). Реальные степени сохраняются.
    """
    rng = np.random.default_rng(seed)
    deg = power_law_degrees(n, gamma, d_min=d_min, d_max=d_max, seed=seed)
    labels = rng.integers(0, n_communities, size=n)
    g = _block_stub_match(deg, labels, mu=mu, rng=rng)
    g.graph["d_min"] = d_min
    g.graph["gamma"] = gamma
    g.graph["mu"] = mu
    g.graph["n_communities"] = n_communities
    return g, labels


def ensemble_configs(
    n_values,
    gammas,
    d_mins=(1, 2, 3),
    n_communities=2,
    mus=(0.1,),
    d_max=None,
    planted_d_mins=(1,),
):
    """Итерируемая спецификация экспериментов (раздел 12 конспекта).

    d_max=None означает естественный cutoff d_max ~ n^{1/(gamma-1)};
    d_max=k фиксирует максимальную степень (устраняет конфаунд с n).
    mus -- список значений доли межблочных концов для planted-ансамбля.
    planted_d_mins -- значения d_min для planted (d_min>=2 убирает листья
    и делает граф связным, изолируя модульный сигнал от периферии).
    """
    for n in n_values:
        for gamma in gammas:
            for d_min in d_mins:
                yield {
                    "ensemble": f"pl_dmin{d_min}",
                    "kind": "power_law",
                    "n": n,
                    "gamma": gamma,
                    "d_min": d_min,
                    "d_max": d_max,
                }
            for d_min in planted_d_mins:
                for mu in mus:
                    yield {
                        "ensemble": f"planted_dmin{d_min}_mu{mu}",
                        "kind": "planted",
                        "n": n,
                        "gamma": gamma,
                        "d_min": d_min,
                        "n_communities": n_communities,
                        "mu": mu,
                        "d_max": d_max,
                    }


def generate(spec: dict, seed: int) -> tuple[nx.Graph, np.ndarray | None]:
    """Генерирует граф по спецификации из ensemble_configs."""
    n = spec["n"]
    gamma = spec["gamma"]
    d_max = spec.get("d_max")
    if spec["kind"] == "power_law":
        return (
            power_law_graph(n, gamma, d_min=spec["d_min"], d_max=d_max, seed=seed),
            None,
        )
    if spec["kind"] == "planted":
        return planted_power_law_graph(
            n,
            gamma,
            d_min=spec["d_min"],
            n_communities=spec.get("n_communities", 2),
            mu=spec.get("mu", 0.1),
            d_max=d_max,
            seed=seed,
        )
    raise ValueError(f"unknown kind: {spec['kind']}")


def graph_id(spec: dict, rep: int) -> str:
    """Детерминированный идентификатор графа (имя файла без расширения)."""
    ens = str(spec["ensemble"])
    parts = [ens, f"g{spec['gamma']}", f"n{spec['n']}", f"d{spec['d_min']}"]
    mu = spec.get("mu")
    if mu not in (None, "") and f"mu{mu}" not in ens:
        parts.append(f"mu{mu}")
    parts.append(f"r{rep}")
    return "_".join(parts)


def save_graph(g: nx.Graph, path: str, labels=None) -> None:
    """Сохраняет граф как сжатый .npz (edges + n + опционально labels).

    Узлы перенумеровываются в 0..n-1, чтобы изолированные/удалённые вершины
    не создавали дыр; для planted сохраняются метки сообществ.
    """
    nodes = list(g.nodes())
    index = {v: i for i, v in enumerate(nodes)}
    edges = np.array([[index[u], index[v]] for u, v in g.edges()], dtype=np.int32)
    if edges.size == 0:
        edges = np.zeros((0, 2), dtype=np.int32)
    arrays = {
        "edges": edges,
        "n": np.array([len(nodes)], dtype=np.int32),
    }
    if labels is not None:
        lab = np.asarray(labels)
        arrays["labels"] = np.array([lab[v] for v in nodes], dtype=np.int32)
    np.savez_compressed(path, **arrays)


def load_graph(path: str) -> tuple[nx.Graph, np.ndarray | None]:
    """Загружает граф, сохранённый save_graph."""
    data = np.load(path, allow_pickle=False)
    n = int(data["n"][0])
    g = nx.Graph()
    g.add_nodes_from(range(n))
    edges = data["edges"]
    if edges.size:
        g.add_edges_from(map(tuple, edges.tolist()))
    labels = data["labels"] if "labels" in data else None
    return g, labels


if __name__ == "__main__":
    for spec in ensemble_configs([1000], [2.5]):
        g, labels = generate(spec, seed=0)
        print(spec["ensemble"], "n", g.number_of_nodes(), "m", g.number_of_edges())
