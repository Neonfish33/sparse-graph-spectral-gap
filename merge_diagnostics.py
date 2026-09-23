"""Слияние diagnostics.csv из нескольких прогонов в один для анализа α(n).

Нужен потому, что большой прогон (`results_large`) добавил только 2 новых
значения n (3×10⁴, 10⁵), а `fit_powerlaw` требует ≥3 точек по n. Малые n
(n = 100…10⁴) берутся из `results_v2`.

Запуск:
    python merge_diagnostics.py
    python analyze_results.py --csv results_merged/diagnostics.csv --out results_merged/analysis
"""

from __future__ import annotations

import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
SOURCES = [
    os.path.join(BASE, "results_v2", "diagnostics.csv"),
    os.path.join(BASE, "results_large", "diagnostics.csv"),
]
OUT = os.path.join(BASE, "results_merged")
TARGETS = [100, 200, 500, 1000, 3000, 10000, 30000, 100000]


def nearest_target(n: int) -> int:
    return min(TARGETS, key=lambda t: abs(t - n))


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    with open(SOURCES[-1], encoding="utf-8") as fh:
        fields = csv.DictReader(fh).fieldnames

    rows = []
    filled = 0
    for src in SOURCES:
        with open(src, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                row = {k: r.get(k, "") for k in fields}
                if not row.get("n_target"):
                    n = int(float(r["n"]))
                    tgt = nearest_target(n)
                    if abs(tgt - n) > max(5, 0.05 * n):
                        sys.stderr.write(f"[warn] n={n} далеко от цели {tgt}\n")
                    row["n_target"] = str(tgt)
                    filled += 1
                rows.append(row)

    path = os.path.join(OUT, "diagnostics.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"merged -> {path}")
    print(f"  всего строк: {len(rows)} | восстановлено n_target: {filled}")


if __name__ == "__main__":
    main()
