"""Асимптотическая проверка: добавляем n=3e5,1e6 к основному датасету (до 1e5)
и смотрим, меняется ли форма спада lambda2(giant).

Сравниваем на полном диапазоне n=100..1e6:
  M1: A n^-alpha          (степенной)
  M2: A (ln n)^-beta      (логарифмический, предсказание 2008 при q_m<=2)
плюс локальные наклоны по декадам.

Запуск: python compare_asymptotics.py
"""

from __future__ import annotations

import csv
import os
import numpy as np

MAIN = os.path.join("results_merged", "diagnostics.csv")
ASYM = os.path.join("results_asymptotic", "diagnostics.csv")
L2 = "\u03bb\u2082_giant"
GAMMAS = [2.5, 3.0]
ENS = ["pl_dmin2", "pl_dmin3"]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_means():
    by = {}
    for path in (MAIN, ASYM):
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                ens = r["ensemble"]
                if ens not in ENS:
                    continue
                g = fnum(r["\u03b3"])
                if g not in GAMMAS:
                    continue
                n = fnum(r["n_target"])
                v = fnum(r[L2])
                if v == v:
                    by.setdefault((ens, g), {}).setdefault(n, []).append(v)
    return by


def fit_power(ns, ys):
    lx, ly = np.log(ns), np.log(ys)
    s, b = np.polyfit(lx, ly, 1)
    r2 = 1 - np.sum((ly - (s * lx + b)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return -s, r2


def fit_logpow(ns, ys):
    lx, ly = np.log(np.log(ns)), np.log(ys)
    s, b = np.polyfit(lx, ly, 1)
    r2 = 1 - np.sum((ly - (s * lx + b)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return -s, r2


def main():
    by = load_means()
    for (ens, g) in sorted(by):
        d = by[(ens, g)]
        ns = np.array(sorted(d), float)
        ys = np.array([np.mean(d[n]) for n in ns])
        print(f"\n=== {ens} gamma={g} ===")
        print("n / mean lambda2_giant:")
        print("  " + "  ".join(f"{int(n)}:{y:.5f}" for n, y in zip(ns, ys)))
        a, r2a = fit_power(ns, ys)
        b, r2b = fit_logpow(ns, ys)
        print(f"  full (100..1e6): power alpha={a:.4f} R2={r2a:.4f} | "
              f"log beta={b:.4f} R2={r2b:.4f}")
        # до 1e5 vs с 1e6
        m = ns <= 1e5
        a1, r1 = fit_power(ns[m], ys[m])
        a2, r2 = fit_power(ns, ys)
        print(f"  alpha(<=1e5)={a1:.4f} R2={r1:.4f}  ->  alpha(<=1e6)={a2:.4f} R2={r2:.4f}")
        print("  local slopes (log-log, per step):")
        for i in range(len(ns) - 1):
            la = -np.log(ys[i + 1] / ys[i]) / np.log(ns[i + 1] / ns[i])
            print(f"    {int(ns[i])}->{int(ns[i+1])}: {la:.4f}")


if __name__ == "__main__":
    main()
