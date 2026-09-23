"""Этап 1.5: ER-контроль.

Вопрос: скейлинг lambda2 ~ n^(-alpha) --- следствие тяжёлого хвоста степеней
или общий эффект разреженных графов с min degree <= 2?

Генерируем Erdos-Renyi G(n, c/n) (Poisson-степени, нет тяжёлого хвоста) для
c in {2,3,4} и считаем lambda2 нормализованного лапласиана для полного графа,
giant component и 2-core --- как в основных экспериментах.

Если ER giant/2-core тоже дают n^(-alpha) с alpha ~ 0.15 -> эффект generic
(min degree <= 2), а не power-law.  Если спадает иначе -> хвост важен.

Запуск:
    python run_er_control.py --jobs 6 --out results_er
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import networkx as nx
import numpy as np

from src.diagnostics import analyze_graph

N_VALUES = [100, 200, 500, 1000, 3000, 10000, 30000, 100000]
C_VALUES = [2, 3, 4]
FIELDS = [
    "model", "c", "n_target", "rep", "seed", "n", "m", "avg_deg", "n_components",
    "N1", "N2", "giant_frac", "giant_n", "giant_lambda2", "giant_phi_star",
    "giant_ipr2", "core_frac", "core_n", "core_lambda2", "core_phi_star", "core_ipr2",
]


def _one(task):
    c, n, rep, seed = task
    p = c / max(n - 1, 1)
    g = nx.fast_gnp_random_graph(n, p, seed=seed)
    row = analyze_graph(g, k=3)
    keep = {k: row.get(k, float("nan")) for k in FIELDS if k in row}
    keep.update({"model": "ER", "c": c, "n_target": n, "rep": rep, "seed": seed})
    return {k: keep.get(k, float("nan")) for k in FIELDS}


def fit_alpha(ns, ys):
    ns = np.asarray(ns, float)
    ys = np.asarray(ys, float)
    m = np.isfinite(ns) & np.isfinite(ys) & (ns > 0) & (ys > 0)
    if m.sum() < 3:
        return float("nan"), float("nan")
    lx, ly = np.log(ns[m]), np.log(ys[m])
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    r2 = 1 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return -slope, r2


def summarize(rows, outdir):
    lines = ["=== ER-контроль: alpha в lambda2 ~ n^(-alpha) ===",
             f"{'c':>3} {'obs':>6} {'alpha':>8} {'R2':>6}  (n точек)"]
    summary = []
    for c in sorted({r["c"] for r in rows}):
        for obs, col in (("giant", "giant_lambda2"), ("core", "core_lambda2")):
            by_n = {}
            for r in rows:
                if r["c"] != c:
                    continue
                by_n.setdefault(r["n_target"], []).append(float(r[col]))
            ns = sorted(by_n)
            ys = [np.mean(by_n[n]) for n in ns]
            a, r2 = fit_alpha(ns, ys)
            lines.append(f"{c:>3} {obs:>6} {a:>8.4f} {r2:>6.3f}  ({len(ns)})")
            summary.append({"c": c, "obs": obs, "alpha": a, "r2": r2})
    print("\n".join(lines))
    with open(os.path.join(outdir, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(outdir, "scaling_by_c.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["c", "obs", "alpha", "r2"])
        w.writeheader()
        w.writerows(summary)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--out", default="results_er")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    tasks = [
        (c, n, rep, 1000 * c + n + rep)
        for c in C_VALUES
        for n in N_VALUES
        for rep in range(args.reps)
    ]
    print(f"[er] всего графов: {len(tasks)} (c={C_VALUES}, n={N_VALUES}, reps={args.reps})")

    csv_path = os.path.join(args.out, "diagnostics.csv")
    rows = []
    t0 = time.time()
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        done = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(_one, t) for t in tasks]
            for fut in as_completed(futs):
                row = fut.result()
                rows.append(row)
                w.writerow(row)
                fh.flush()
                done += 1
                if done % max(1, len(tasks) // 20) == 0 or done == len(tasks):
                    print(f"[er] done={done}/{len(tasks)}  {time.time()-t0:6.1f}s")
    print(f"[er] rows={len(rows)} -> {csv_path}")
    summarize(rows, args.out)
    print(f"[er] выходные файлы в {args.out}/")


if __name__ == "__main__":
    main()
