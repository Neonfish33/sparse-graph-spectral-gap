"""Оркестрация численных экспериментов (раздел 12 конспекта).

Примеры:
    python run_experiments.py --quick
    python run_experiments.py --n 1000 3000 --reps 20 --out results/
    python run_experiments.py --n 1000 3000 10000 --gammas 2.1 2.5 3.0 \
        --dmin 1 2 3 --reps 50 --out results/
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.diagnostics import analyze_graph, graph_sweep_curve
from src.generators import ensemble_configs, generate, graph_id, save_graph

SIMPLE_LABELS = {
    "avg_deg": "d̄",
    "d_min_actual": "d_min_obs",
    "d_max": "d_max",
    "n_components": "n_comp",
    "N1": "N₁",
    "N2": "N₂",
    "gamma": "γ",
    "mu": "μ",
    "giant_frac": "frac_giant",
    "core_frac": "frac_2core",
}

METRIC_LABELS = {
    "n": "n",
    "m": "m",
    "lambda2": "λ₂",
    "lambda3": "λ₃",
    "gap23": "(λ₃−λ₂)",
    "phi_periph": "φ_periph",
    "phi_mod": "φ_mod",
    "phi_star": "φ*",
    "rho_star": "ρ*",
    "rho_2": "ρ₂",
    "log_Delta": "logΔ",
    "ipr2": "IPR₂",
    "ipr2_pi": "IPR₂^π",
    "pr2": "PR₂",
}


def pretty_name(key: str) -> str:
    """ASCII-ключ диагностики -> математическое обозначение для вывода."""
    if key in SIMPLE_LABELS:
        return SIMPLE_LABELS[key]
    for prefix in ("full", "giant", "core"):
        p = prefix + "_"
        if key.startswith(p) and key[len(p):] in METRIC_LABELS:
            return f"{METRIC_LABELS[key[len(p):]]}_{prefix}"
    return key


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, nargs="+", default=[1000, 3000, 10000])
    p.add_argument("--gammas", type=float, nargs="+", default=[2.1, 2.5, 3.0])
    p.add_argument("--dmin", type=int, nargs="+", default=[1, 2, 3])
    p.add_argument("--reps", type=int, default=10)
    p.add_argument("--q", type=int, default=2, help="число planted communities")
    p.add_argument(
        "--mu", type=float, nargs="+", default=[0.1],
        help="доля межблочных концов (можно несколько значений)",
    )
    p.add_argument(
        "--dmax", type=int, default=None,
        help="фиксировать d_max (по умолчанию естественный n^{1/(gamma-1)})",
    )
    p.add_argument("--eps", type=float, default=0.05, help="порог периферийного объёма")
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--out", type=str, default="results")
    p.add_argument(
        "--planted-dmin", type=int, nargs="+", default=[1],
        help="d_min для planted-ансамбля (>=2 убирает листья)",
    )
    p.add_argument("--no-planted", action="store_true",
                   help="не генерировать planted-ансамбль")
    p.add_argument("--jobs", type=int, default=1, help="число процессов (1 = последовательно)")
    p.add_argument(
        "--save-curves", action="store_true",
        help="сохранить φ-профили (giant component, rep=0) в phi_curves.csv",
    )
    p.add_argument(
        "--save-graphs", action="store_true",
        help="сохранить каждый граф в {out}/graphs/*.npz (воспроизводимость)",
    )
    p.add_argument(
        "--resume", action="store_true",
        help="продолжить: пропустить уже готовые графы (по graph_id в CSV)",
    )
    p.add_argument("--quick", action="store_true")
    return p.parse_args(argv)


CURVE_MAX_POINTS = 500


def _curve_rows(curve, spec, rep):
    """Строки φ-профиля с даунсэмплингом до ~CURVE_MAX_POINTS на сторону."""
    rows = []
    if not curve:
        return rows
    for side in ("asc", "desc"):
        rho, phi = curve[side]
        m = len(rho)
        if m > CURVE_MAX_POINTS:
            step = int(math.ceil(m / CURVE_MAX_POINTS))
            idx = list(range(0, m, step))
            if idx[-1] != m - 1:
                idx.append(m - 1)
        else:
            idx = range(m)
        for i in idx:
            rows.append(
                {
                    "ensemble": spec["ensemble"],
                    "gamma": spec["gamma"],
                    "d_min": spec["d_min"],
                    "mu": spec.get("mu", ""),
                    "n": spec["n"],
                    "rep": rep,
                    "side": side,
                    "rho": float(rho[i]),
                    "phi": float(phi[i]),
                }
            )
    return rows


def _run_one(task):
    """Одна реализация: генерация + диагностики (+ φ-профиль, + сохранение графа)."""
    spec, rep, seed, eps, want_curve, graphs_dir = task
    g, labels = generate(spec, seed)
    row = {
        "ensemble": spec["ensemble"],
        "gamma": spec["gamma"],
        "d_min": spec["d_min"],
        "mu": spec.get("mu", ""),
        "n_target": spec["n"],
        "rep": rep,
        "seed": seed,
    }
    gid = graph_id(spec, rep)
    row["graph_id"] = gid
    if graphs_dir:
        path = os.path.join(graphs_dir, gid + ".npz")
        save_graph(g, path, labels)
        row["graph_file"] = gid + ".npz"
    row.update(analyze_graph(g, eps=eps))
    curves = _curve_rows(graph_sweep_curve(g), spec, rep) if want_curve else []
    return row, curves


def mean_or_nan(values):
    vals = [v for v in values if v == v and not math.isinf(v)]
    return statistics.fmean(vals) if vals else float("nan")


def scaling_exponent(rows, ensemble, field):
    """Оценка alpha и R^2 в log(lambda2) = -alpha log(n) + C.

    Регрессия по средним значениям на каждом n. Возвращает
    (alpha, intercept, R^2, n_points).
    """
    by_n = {}
    for r in rows:
        if r["ensemble"] != ensemble:
            continue
        v = r[field]
        if v != v or v <= 0 or math.isinf(v):
            continue
        by_n.setdefault(r.get("n_target", r["n"]), []).append(v)
    ns = sorted(by_n)
    if len(ns) < 2:
        return float("nan"), float("nan"), float("nan"), 0
    x = np.log(np.array(ns, dtype=float))
    y = np.log(np.array([mean_or_nan(by_n[n]) for n in ns], dtype=float))
    mask = np.isfinite(y)
    if mask.sum() < 2:
        return float("nan"), float("nan"), float("nan"), int(mask.sum())
    xm, ym = x[mask], y[mask]
    slope, intercept = np.polyfit(xm, ym, 1)
    pred = slope * xm + intercept
    ss_res = float(np.sum((ym - pred) ** 2))
    ss_tot = float(np.sum((ym - ym.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return -slope, intercept, r2, int(mask.sum())


def main(argv=None):
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    args = parse_args(argv)
    if args.quick:
        args.n = [400, 800]
        args.gammas = [2.5]
        args.reps = 3

    os.makedirs(args.out, exist_ok=True)
    csv_path = os.path.join(args.out, "diagnostics.csv")

    specs = list(
        ensemble_configs(
            args.n,
            args.gammas,
            d_mins=args.dmin,
            n_communities=args.q,
            mus=args.mu,
            d_max=args.dmax,
            planted_d_mins=[] if args.no_planted else args.planted_dmin,
        )
    )
    graphs_dir = os.path.join(args.out, "graphs") if args.save_graphs else None
    if graphs_dir:
        os.makedirs(graphs_dir, exist_ok=True)

    all_tasks = [
        (spec, rep, args.seed + rep, args.eps, args.save_curves and rep == 0, graphs_dir)
        for spec in specs
        for rep in range(args.reps)
    ]

    completed = set()
    header_fields = None
    if args.resume and os.path.exists(csv_path):
        with open(csv_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            header_fields = list(reader.fieldnames or [])
            for r in reader:
                gid = r.get("graph_id")
                if gid:
                    completed.add(gid)
        if "graph_id" not in header_fields:
            print("[run] WARNING: нет колонки graph_id -> resume невозможен, начинаю заново")
            completed, header_fields = set(), None
        else:
            print(f"[run] resume: уже готово {len(completed)} графов")

    tasks = [t for t in all_tasks if graph_id(t[0], t[1]) not in completed]
    total = len(tasks)
    print(f"[run] specs={len(specs)} reps={args.reps} всего={len(all_tasks)} "
          f"к выполнению={total} (пропущено {len(all_tasks) - total})")
    print(f"[run] output -> {csv_path}")

    rows = []
    state = {"fieldnames": None, "writer": None, "curve_writer": None, "done": 0}
    t0 = time.time()
    curve_path = os.path.join(args.out, "phi_curves.csv") if args.save_curves else None
    curve_append = bool(args.resume and curve_path and os.path.exists(curve_path))
    if curve_path:
        curve_fh = open(
            curve_path, "a" if curve_append else "w", newline="", encoding="utf-8"
        )
        if curve_append:
            with open(curve_path, encoding="utf-8") as cf:
                try:
                    cf_fields = next(csv.reader(cf))
                except StopIteration:
                    cf_fields = None
            if cf_fields:
                state["curve_writer"] = csv.DictWriter(curve_fh, fieldnames=cf_fields)
    else:
        curve_fh = None

    def emit(row, curves, fh):
        rows.append(row)
        if state["writer"] is None:
            state["fieldnames"] = list(row.keys())
            state["writer"] = csv.DictWriter(
                fh, fieldnames=[pretty_name(k) for k in state["fieldnames"]]
            )
            state["writer"].writeheader()
        state["writer"].writerow({pretty_name(k): v for k, v in row.items()})
        if curve_fh is not None and curves:
            if state["curve_writer"] is None:
                state["curve_writer"] = csv.DictWriter(
                    curve_fh, fieldnames=list(curves[0].keys())
                )
                state["curve_writer"].writeheader()
            for cr in curves:
                state["curve_writer"].writerow(cr)
            curve_fh.flush()
        state["done"] += 1
        step = max(1, total // 20)
        if total and (state["done"] % step == 0 or state["done"] == total):
            print(
                f"[run] done={state['done']}/{total}  "
                f"{time.time() - t0:6.1f}s"
            )

    csv_mode = "a" if header_fields is not None else "w"
    with open(csv_path, csv_mode, newline="", encoding="utf-8") as fh:
        if header_fields is not None:
            state["writer"] = csv.DictWriter(fh, fieldnames=header_fields)
        if args.jobs and args.jobs > 1:
            with ProcessPoolExecutor(max_workers=args.jobs) as ex:
                futures = [ex.submit(_run_one, t) for t in tasks]
                for fut in as_completed(futures):
                    row, curves = fut.result()
                    emit(row, curves, fh)
        else:
            for t in tasks:
                row, curves = _run_one(t)
                emit(row, curves, fh)
        fh.flush()
    if curve_fh is not None:
        curve_fh.close()

    ensembles = sorted({r["ensemble"] for r in rows})

    print("\n[summary] среднее λ₂ по ансамблям (G / G_giant / G_2core):")
    header = f"{'ensemble':>12} {'γ':>6} {'d_min':>6} " \
             f"{'λ₂(G)':>10} {'λ₂(G_giant)':>12} {'λ₂(G_2core)':>12}"
    print(header)
    print("-" * 64)
    for ens in ensembles:
        sub = [r for r in rows if r["ensemble"] == ens]
        for gamma in sorted({r["gamma"] for r in sub}):
            s = [r for r in sub if r["gamma"] == gamma]
            if not s:
                continue
            print(
                f"{ens:>12} {gamma:>6} {s[0]['d_min']:>6} "
                f"{mean_or_nan([r['full_lambda2'] for r in s]):>10.4f} "
                f"{mean_or_nan([r['giant_lambda2'] for r in s]):>12.4f} "
                f"{mean_or_nan([r['core_lambda2'] for r in s]):>12.4f}"
            )

    print("\n[scaling] λ₂(n) ~ n^(−α), регрессия log λ₂ vs log n (R², n точек):")
    header2 = f"{'ensemble':>12} {'γ':>6} {'d_min':>6} " \
              f"{'α(G)':>7} {'R²(G)':>6} " \
              f"{'α(giant)':>9} {'R²(giant)':>9} " \
              f"{'α(core)':>8} {'R²(core)':>8}"
    print(header2)
    print("-" * len(header2))
    for ens in ensembles:
        sub = [r for r in rows if r["ensemble"] == ens]
        for gamma in sorted({r["gamma"] for r in sub}):
            s = [r for r in sub if r["gamma"] == gamma]
            if not s:
                continue
            af, _, rf, nf = scaling_exponent(s, ens, "full_lambda2")
            ag, _, rg, ng = scaling_exponent(s, ens, "giant_lambda2")
            ac, _, rc, nc = scaling_exponent(s, ens, "core_lambda2")
            print(
                f"{ens:>12} {gamma:>6} {s[0]['d_min']:>6} "
                f"{af:>7.3f} {rf:>6.3f} "
                f"{ag:>9.3f} {rg:>9.3f} "
                f"{ac:>8.3f} {rc:>8.3f}"
            )
    print("  (α по <2 точкам n не оценивается -> nan)")

    print(f"\n[done] rows={len(rows)} -> {csv_path}")


if __name__ == "__main__":
    main()
