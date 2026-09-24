"""ERT<->power-law ensemble matcher.

Stage `grid`  : coarse full-factorial metric grid (ER: c x n ; PL: D x gamma x n).
Stage `final` : high-rep runs at the matched pairs chosen from the grid.

Only writes into results_agents/matcher/.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mlib import gen_er, gen_pl, metrics, fit_alpha  # noqa: E402

OUT = HERE

N_GRID = [1000, 3000, 10000, 30000]
C_GRID = [2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0]
D_GRID = [10, 20, 30, 50, 80, 120]
GAMMAS = [2.1, 2.5, 3.0]

MCOL = ["n", "m", "avg_deg", "N2", "P2", "q", "giant_frac", "giant_n",
        "giant_lambda2", "core_frac", "core_n", "core_m", "d_bar_core",
        "core_lambda2", "core_P2", "q_core"]


def _seed(*parts) -> int:
    x = 1469598103934665603
    for p in parts:
        x = (x ^ (int(p) & 0xFFFFFFFF)) * 1099511628211
        x &= 0xFFFFFFFFFFFFFFFF
    return int(x & 0x7FFFFFFF)


def _task_er(args):
    n, c, rep = args
    g = gen_er(n, float(c), _seed(1, n, int(round(c * 1000)), rep))
    r = metrics(g)
    r.update({"model": "ER", "n_target": n, "c": c, "gamma": "", "D": "", "rep": rep})
    return r


def _task_pl(args):
    n, gamma, D, rep = args
    g = gen_pl(n, float(gamma), int(D), _seed(2, n, int(round(float(gamma) * 1000)), D, rep))
    r = metrics(g)
    r.update({"model": "PL", "n_target": n, "c": "", "gamma": gamma, "D": D, "rep": rep})
    return r


FIELDS = ["model", "n_target", "c", "gamma", "D", "rep"] + MCOL


def _run(tasks, worker, path, jobs):
    rows = []
    t0 = time.time()
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        done = 0
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(worker, t) for t in tasks]
            for fut in as_completed(futs):
                row = fut.result()
                rows.append(row)
                w.writerow({k: row.get(k, "") for k in FIELDS})
                fh.flush()
                done += 1
                if done % max(1, len(tasks) // 10) == 0 or done == len(tasks):
                    print(f"  done={done}/{len(tasks)} {time.time()-t0:6.1f}s", flush=True)
    print(f"  -> {path} ({len(rows)} rows)", flush=True)
    return rows


def stage_grid(jobs, reps):
    er_tasks = [(n, c, rep) for n in N_GRID for c in C_GRID for rep in range(reps)]
    pl_tasks = [(n, g, D, rep) for n in N_GRID for g in GAMMAS for D in D_GRID for rep in range(reps)]
    print(f"[grid] ER={len(er_tasks)} PL={len(pl_tasks)} tasks, jobs={jobs}", flush=True)
    _run(er_tasks, _task_er, os.path.join(OUT, "grid_er.csv"), jobs)
    _run(pl_tasks, _task_pl, os.path.join(OUT, "grid_pl.csv"), jobs)


def read_csv(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in MCOL:
            try:
                r[k] = float(r[k])
            except (ValueError, TypeError):
                r[k] = float("nan")
        r["n_target"] = int(float(r["n_target"]))
        r["rep"] = int(float(r["rep"]))
        if r.get("c"):
            r["c"] = float(r["c"])
        if r.get("D"):
            r["D"] = int(float(r["D"]))
        if r.get("gamma"):
            r["gamma"] = float(r["gamma"])
    return rows


def mean_by(rows, keyfns, col):
    by = {}
    for r in rows:
        k = keyfns(r)
        by.setdefault(k, []).append(r[col])
    return {k: float(np.nanmean(v)) for k, v in by.items()}


def choose_matches(er_rows, pl_rows):
    """Return dict gamma -> dict n -> (c, D, pl_dbar, pl_q, er_dbar_curve...)."""
    er_n = mean_by(er_rows, lambda r: (r["n_target"], r["c"]), "d_bar_core")
    er_q = mean_by(er_rows, lambda r: (r["n_target"], r["c"]), "q")
    pl_d = mean_by(pl_rows, lambda r: (r["n_target"], r["gamma"], r["D"]), "d_bar_core")
    pl_q = mean_by(pl_rows, lambda r: (r["n_target"], r["gamma"], r["D"]), "q")

    matches = {}
    for n in N_GRID:
        cs = sorted(c for (nn, c) in er_n if nn == n)
        yc = np.array([er_n[(n, c)] for c in cs])
        yq = np.array([er_q[(n, c)] for c in cs])
        for g in GAMMAS:
            best = None
            for D in D_GRID:
                if (n, g, D) not in pl_d:
                    continue
                dpl = pl_d[(n, g, D)]
                qpl = pl_q[(n, g, D)]
                if dpl < yc.min() or dpl > yc.max():
                    cstar = cs[int(np.argmin(np.abs(yc - dpl)))] if dpl < yc.min() else cs[-1]
                    der = min(yc.max(), max(yc.min(), dpl))
                else:
                    cstar = float(np.interp(dpl, yc, cs))
                    der = float(np.interp(cstar, cs, yc))
                qer = float(np.interp(cstar, cs, yq))
                cost = abs(qer - qpl)
                cand = (abs(dpl - der), cost, cstar, D, dpl, qpl, qer, der)
                if best is None or cand[:2] < best[:2]:
                    best = cand
            matches[(n, g)] = best
    return matches, (er_n, er_q, pl_d, pl_q)


def _task_final_er(args):
    n, c, gamma, D, rep = args
    g = gen_er(n, float(c), _seed(11, n, int(round(float(c) * 1000)), rep))
    r = metrics(g)
    r.update({"model": "ER", "n_target": n, "c": c, "gamma": gamma, "D": D, "rep": rep})
    return r


def _task_final_pl(args):
    n, c, gamma, D, rep = args
    g = gen_pl(n, float(gamma), int(D), _seed(12, n, int(round(float(gamma) * 1000)), D, rep))
    r = metrics(g)
    r.update({"model": "PL", "n_target": n, "c": c, "gamma": gamma, "D": D, "rep": rep})
    return r


def stage_final(jobs, reps, matches):
    er_tasks, pl_tasks = [], []
    info = {}
    for g in GAMMAS:
        ref = matches[(10000, g)]
        cstar, Dstar = float(ref[2]), int(ref[3])
        info[g] = (cstar, Dstar)
        print(f"[final] gamma={g} c*={cstar:.4f} D*={Dstar}", flush=True)
        for n in N_GRID:
            for rep in range(reps):
                er_tasks.append((n, cstar, g, Dstar, rep))
                pl_tasks.append((n, cstar, g, Dstar, rep))
    print(f"[final] ER={len(er_tasks)} PL={len(pl_tasks)} jobs={jobs}", flush=True)
    _run(er_tasks, _task_final_er, os.path.join(OUT, "final_er.csv"), jobs)
    _run(pl_tasks, _task_final_pl, os.path.join(OUT, "final_pl.csv"), jobs)
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["grid", "final", "all"], default="all")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--reps-grid", type=int, default=5)
    ap.add_argument("--reps-final", type=int, default=12)
    args = ap.parse_args()

    if args.stage in ("grid", "all"):
        stage_grid(args.jobs, args.reps_grid)
    if args.stage in ("final", "all"):
        er_rows = read_csv(os.path.join(OUT, "grid_er.csv"))
        pl_rows = read_csv(os.path.join(OUT, "grid_pl.csv"))
        matches, _ = choose_matches(er_rows, pl_rows)
        stage_final(args.jobs, args.reps_final, matches)


if __name__ == "__main__":
    main()
