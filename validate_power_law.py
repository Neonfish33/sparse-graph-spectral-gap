"""Этап 0: является ли спад lambda2(n) степенным законом?

Источник: results_merged/diagnostics.csv (сырые реализации).
Для каждого (d_min, gamma) и наблюдаемой lambda2_giant:
  * сравнение моделей через AICc:
      M1 power law   : lambda2 = A * n^(-alpha)              (k=2)
      M2 inverse log : lambda2 = A / log n                   (k=1)
      M3 exp-log     : lambda2 = A * exp(-c * (log n)^beta)  (k=3)
  * bootstrap CI для alpha (ресэмпл реализаций ВНУТРИ каждого n, B=1000);
  * локальный alpha по декадам n.

Запуск:
    python validate_power_law.py
"""

from __future__ import annotations

import csv
import math
import os
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit

CSV = os.path.join("results_merged", "diagnostics.csv")
OBS = "\u03bb\u2082_giant"          # lambda2_giant
DMINS = [2.0, 3.0]
GAMMAS = [2.1, 2.5, 3.0]
ENS = {"pl_dmin2": 2.0, "pl_dmin3": 3.0}   # только чистые power-law ансамбли
B = 1000
RNG = np.random.default_rng(0)


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def gather():
    g = defaultdict(lambda: defaultdict(list))
    with open(CSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            ens = r["ensemble"]
            if ens not in ENS:
                continue
            d = ENS[ens]
            gam, n, v = fnum(r["\u03b3"]), fnum(r["n_target"]), fnum(r[OBS])
            if gam in GAMMAS and v == v:
                g[(d, gam)][n].append(v)
    return g


def fit_power_law(ns, ys):
    lx, ly = np.log(ns), np.log(ys)
    slope, intercept = np.polyfit(lx, ly, 1)
    return -slope, math.exp(intercept)


def ssr_power_law(ns, ys):
    a, A = fit_power_law(ns, ys)
    pred = A * ns ** (-a)
    return float(np.sum((ys - pred) ** 2)), a


def ssr_inv_log(ns, ys):
    x = 1.0 / np.log(ns)
    A = float(np.sum(ys * x) / np.sum(x * x))
    return float(np.sum((ys - A * x) ** 2)), A


def ssr_exp_log(ns, ys):
    def f(n, A, c, b):
        return A * np.exp(-c * np.log(n) ** b)
    try:
        p, _ = curve_fit(f, ns, ys, p0=[ys[0], 0.1, 1.0],
                         bounds=([0, 0, 0.1], [np.inf, 10, 5]), maxfev=40000)
        return float(np.sum((ys - f(ns, *p)) ** 2)), p
    except Exception:
        return float("inf"), None


def aicc(ssr, k, n):
    if not np.isfinite(ssr) or ssr <= 0 or n - k - 1 <= 0:
        return float("inf")
    return n * math.log(ssr / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)


def bootstrap_alpha(by_n, ns, B=B):
    alphas = []
    for _ in range(B):
        means = []
        for n in ns:
            v = np.asarray(by_n[n])
            means.append(float(v[RNG.integers(0, len(v), len(v))].mean()))
        try:
            a, _ = fit_power_law(ns, np.asarray(means))
            alphas.append(a)
        except Exception:
            pass
    if len(alphas) < 50:
        return (float("nan"), float("nan"))
    return tuple(np.percentile(alphas, [2.5, 97.5]))


def main():
    g = gather()
    for (d, gam) in sorted(g):
        by_n = g[(d, gam)]
        ns = np.asarray(sorted(by_n), dtype=float)
        ys = np.asarray([np.mean(by_n[n]) for n in ns])
        if len(ns) < 4:
            continue
        ssr1, sigma_lin = ssr_power_law(ns, ys)
        ssr2, _ = ssr_inv_log(ns, ys)
        ssr3, p3 = ssr_exp_log(ns, ys)
        n = len(ns)
        aic = {"power_law": aicc(ssr1, 2, n),
               "1/log": aicc(ssr2, 1, n),
               "exp-log": aicc(ssr3, 3, n)}
        best = min(aic, key=aic.get)
        a, A = fit_power_law(ns, ys)
        ss_tot = float(np.sum((ys - ys.mean()) ** 2))
        r2 = 1 - ssr1 / ss_tot
        ci = bootstrap_alpha(by_n, ns)
        print(f"\n=== d_min={d:.0f}, \u03b3={gam} (точек n: {n}) ===")
        print(f"  power law (log-fit): \u03b1 = {a:.4f}, A = {A:.4f}, R\u00b2 = {r2:.5f}")
        print(f"  bootstrap 95% CI \u03b1: [{ci[0]:.4f}, {ci[1]:.4f}]")
        for k in ("power_law", "1/log", "exp-log"):
            print(f"    AICc[{k:>9}] = {aic[k]:8.2f}")
        print(f"  WINNER: {best}   (M3 params: {None if p3 is None else np.round(p3,3)})")
        print("  local \u03b1 по декадам:")
        for i in range(n - 1):
            la = -math.log(ys[i + 1] / ys[i]) / math.log(ns[i + 1] / ns[i])
            print(f"    n {ns[i]:.0f} \u2192 {ns[i+1]:.0f}: \u03b1_local = {la:.4f}")


if __name__ == "__main__":
    main()
