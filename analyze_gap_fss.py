"""Направление A: finite-size scaling и флуктуации спектральной щели.

A1 (q_m>2, d_min=3): lambda2(n) -> lambda_c > 0. Фит lambda2 = lambda_c + C n^-theta
    (и альтернатива lambda_c + C/log n). Оценка края lambda_c и показателя theta.
A2 (q_m=2, d_min=2): распределение lambda2 по реализациям; стабильность рескейла
    n^alpha * lambda2; относительные флуктуации; корреляция с локализацией (IPR).

Запуск:
    python analyze_gap_fss.py
"""

from __future__ import annotations

import collections
import csv
import math
import os

import numpy as np
from scipy.optimize import curve_fit

CSV = os.path.join("results_merged", "diagnostics.csv")
L2 = "\u03bb\u2082_giant"       # lambda2_giant
IPR = "IPR\u2082_giant"         # IPR2_giant


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def fit_offpow(ns, ys):
    def f(n, C0, C, th):
        return C0 + C * np.asarray(n, float) ** (-th)
    try:
        p, _ = curve_fit(f, ns, ys, p0=[0.5 * min(ys), ys[-1], 0.1],
                         bounds=([0, 1e-9, 1e-3], [np.inf, np.inf, 8]), maxfev=40000)
        return p
    except Exception:
        return (float("nan"),) * 3


def fit_offlog(ns, ys):
    def f(n, C0, C):
        return C0 + C / np.log(np.asarray(n, float))
    try:
        p, _ = curve_fit(f, ns, ys, p0=[0.5 * min(ys), ys[-1]], maxfev=40000)
        return p
    except Exception:
        return (float("nan"),) * 2


def main():
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(CSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["ensemble"] in ("pl_dmin2", "pl_dmin3"):
                by[(r["ensemble"], r["\u03b3"])][fnum(r["n_target"])].append(
                    (fnum(r[L2]), fnum(r[IPR])))
    for (ens, gam) in sorted(by):
        d = by[(ens, gam)]
        ns = np.array(sorted(d), float)
        means = np.array([np.nanmean([x[0] for x in d[n]]) for n in ns])
        stds = np.array([np.nanstd([x[0] for x in d[n]]) for n in ns])
        q05 = np.array([np.nanpercentile([x[0] for x in d[n]], 5) for n in ns])
        q95 = np.array([np.nanpercentile([x[0] for x in d[n]], 95) for n in ns])
        ipr = np.array([np.nanmean([x[1] for x in d[n]]) for n in ns])
        print(f"\n=== {ens} gamma={gam} (reps/n={[len(d[n]) for n in ns][0]}) ===")
        print(f"{'n':>8} {'mean':>9} {'std':>8} {'std/mean':>9} "
              f"{'q05':>8} {'q95':>8} {'q95/q05':>8} {'mean*IPR':>9} {'n*mean':>10}")
        for i, n in enumerate(ns):
            r = stds[i] / means[i] if means[i] else float("nan")
            q = q95[i] / q05[i] if q05[i] else float("nan")
            print(f"{n:>8.0f} {means[i]:>9.5f} {stds[i]:>8.5f} {r:>9.3f} "
                  f"{q05[i]:>8.5f} {q95[i]:>8.5f} {q:>8.3f} "
                  f"{means[i]*ipr[i]:>9.4f} {ns[i]*means[i]:>10.4f}")
        if ens == "pl_dmin3":
            C0, C, th = fit_offpow(ns, means)
            a0, a1 = fit_offlog(ns, means)
            print(f"  A1 fit  lambda_c + C n^-theta : lambda_c={C0:.4f}  C={C:.4f}  theta={th:.3f}")
            print(f"  A1 fit  lambda_c + C/log n    : lambda_c={a0:.4f}  C={a1:.4f}")
        if ens == "pl_dmin2":
            a = -np.polyfit(np.log(ns), np.log(means), 1)[0]
            print(f"  A2 mean ~ n^-alpha, alpha={a:.4f}")
            print(f"  A2 rescaled n^alpha*lambda2: "
                  f"{np.round(ns**a * means, 4).tolist()}")
    print("\nИнтерпретация:")
    print("  A1: если lambda_c>0 заметно, theta показывает скорость выхода на плато.")
    print("  A2: если n^alpha*lambda2 ~ const и std/mean ~ const --- масштабная инвариантность.")


if __name__ == "__main__":
    main()
