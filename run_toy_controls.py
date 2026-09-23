"""Этап 3a: контроли toy-модели --- хвост (A) vs доля degree-2 (B).

Гипотеза A: спад λ₂ специфично усиливается тяжёлым хвостом степеней.
Гипотеза B: всё определяется долей вершин степени 2 (структурный фактор).

Сравниваем сети с min degree 2 и разными распределениями степеней:
  * pl_g<γ>_D<D>  --- power-law на {2..D} (хвост + обрезка);
  * two23_p<p2>   --- носитель {2,3}, доля degree-2 = p2 (без хвоста);
  * b234_p<p2>    --- носитель {2,3,4}, доля degree-2 = p2 (без хвоста).

Метрика-кандидат «структурного» параметра: q = 2·P(2)/d̄ (вероятность, что
ребро ведёт в вершину степени 2).

Запуск:
    python run_toy_controls.py --jobs 6 --reps 20 --out results_toy
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from src.diagnostics import analyze_graph
from src.generators import configuration_graph

N_VALUES = [1000, 10000, 100000]
MODELS = [
    {"name": "pl_g2.5_D3", "kind": "pl", "gamma": 2.5, "D": 3},
    {"name": "pl_g2.5_D4", "kind": "pl", "gamma": 2.5, "D": 4},
    {"name": "pl_g2.5_D5", "kind": "pl", "gamma": 2.5, "D": 5},
    {"name": "pl_g2.5_D8", "kind": "pl", "gamma": 2.5, "D": 8},
    {"name": "pl_g2.5_D15", "kind": "pl", "gamma": 2.5, "D": 15},
    {"name": "pl_g2.5_D50", "kind": "pl", "gamma": 2.5, "D": 50},
    {"name": "two23_p0.4", "kind": "two23", "p2": 0.4},
    {"name": "two23_p0.6", "kind": "two23", "p2": 0.6},
    {"name": "two23_p0.8", "kind": "two23", "p2": 0.8},
    {"name": "b234_p0.5", "kind": "b234", "p2": 0.5},
    {"name": "b234_p0.7", "kind": "b234", "p2": 0.7},
    {"name": "unif_D3", "kind": "unif", "D": 3},
    {"name": "unif_D4", "kind": "unif", "D": 4},
    {"name": "unif_D5", "kind": "unif", "D": 5},
    {"name": "unif_D6", "kind": "unif", "D": 6},
    {"name": "unif_D8", "kind": "unif", "D": 8},
    {"name": "unif_D12", "kind": "unif", "D": 12},
    {"name": "unif_D20", "kind": "unif", "D": 20},
]
FIELDS = [
    "model", "kind", "gamma", "D", "p2_target", "n_target", "rep", "n", "m",
    "avg_deg", "N1", "N2", "p2_emp", "giant_frac", "giant_lambda2",
    "core_frac", "core_n", "core_lambda2",
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


def _one(task):
    spec, n, rep, seed = task
    rng = np.random.default_rng(seed)
    deg = make_degrees(spec, n, rng)
    g = configuration_graph(deg, seed=seed)
    row = analyze_graph(g, k=3)
    nv = row.get("n", 0)
    out = {
        "model": spec["name"], "kind": spec["kind"],
        "gamma": spec.get("gamma", ""), "D": spec.get("D", ""),
        "p2_target": spec.get("p2", ""), "n_target": n, "rep": rep,
        "n": nv, "m": row.get("m"), "avg_deg": row.get("avg_deg"),
        "N1": row.get("N1"), "N2": row.get("N2"),
        "p2_emp": (row.get("N2", 0) / nv) if nv else float("nan"),
        "giant_frac": row.get("giant_frac"),
        "giant_lambda2": row.get("giant_lambda2"),
        "core_frac": row.get("core_frac"), "core_n": row.get("core_n"),
        "core_lambda2": row.get("core_lambda2"),
    }
    return out


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
    by_model = {}
    for r in rows:
        by_model.setdefault(r["model"], []).append(r)
    lines = ["=== toy-контроли: alpha в lambda2(n) ~ n^-alpha ===",
             f"{'model':>14} {'dbar':>6} {'p2':>6} {'q=2p2/dbar':>11} "
             f"{'a_giant':>8} {'R2g':>6} {'a_core':>8} {'R2c':>6}"]
    out = []
    for name in [m["name"] for m in MODELS]:
        rs = by_model.get(name, [])
        if not rs:
            continue
        big = [r for r in rs if float(r["n_target"]) >= 10000]
        dbar = float(np.mean([float(r["avg_deg"]) for r in big]))
        p2 = float(np.mean([float(r["p2_emp"]) for r in big]))
        q = 2 * p2 / dbar if dbar else float("nan")
        ns = sorted({float(r["n_target"]) for r in rs})
        ag, rg = fit_alpha(ns, [np.mean([float(r["giant_lambda2"]) for r in rs if float(r["n_target"]) == n]) for n in ns])
        ac, rc = fit_alpha(ns, [np.mean([float(r["core_lambda2"]) for r in rs if float(r["n_target"]) == n]) for n in ns])
        lines.append(f"{name:>14} {dbar:>6.2f} {p2:>6.3f} {q:>11.3f} "
                     f"{ag:>8.4f} {rg:>6.3f} {ac:>8.4f} {rc:>6.3f}")
        out.append({"model": name, "dbar": dbar, "p2": p2, "q": q,
                    "alpha_giant": ag, "r2_giant": rg, "alpha_core": ac, "r2_core": rc})
    txt = "\n".join(lines)
    print("\n" + txt)
    with open(os.path.join(outdir, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    with open(os.path.join(outdir, "alpha_vs_p2.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--out", default="results_toy")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    csv_path = os.path.join(args.out, "diagnostics.csv")

    all_tasks = [
        (spec, n, rep, abs(hash((spec["name"], n, rep))) % (2 ** 31))
        for spec in MODELS
        for n in N_VALUES
        for rep in range(args.reps)
    ]
    existing = set()
    if args.resume and os.path.exists(csv_path):
        with open(csv_path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                existing.add((r["model"], int(float(r["n_target"])), int(r["rep"])))
    tasks = [t for t in all_tasks if (t[0]["name"], t[1], t[2]) not in existing]
    print(f"[toy] моделей={len(MODELS)} n={N_VALUES} reps={args.reps} "
          f"всего={len(all_tasks)} к выполнению={len(tasks)} (пропущено {len(all_tasks)-len(tasks)})")

    t0 = time.time()
    with open(csv_path, "a" if existing else "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if not existing:
            w.writeheader()
        done = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(_one, t) for t in tasks]
            for fut in as_completed(futs):
                w.writerow(fut.result())
                fh.flush()
                done += 1
                step = max(1, len(tasks) // 20)
                if done % step == 0 or done == len(tasks):
                    print(f"[toy] done={done}/{len(tasks)}  {time.time()-t0:6.1f}s", flush=True)

    with open(csv_path, encoding="utf-8") as fh:
        all_rows = list(csv.DictReader(fh))
    print(f"[toy] всего строк={len(all_rows)} -> {csv_path}")
    summarize(all_rows, args.out)
    print(f"[toy] выходные файлы в {args.out}/")


if __name__ == "__main__":
    main()
