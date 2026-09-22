"""Анализ результатов run_experiments.py.

Строит:
  * scaling_by_n.csv    -- α и R² в λ₂ ~ n^(−α) для G / G_giant / G_2core;
  * scaling_by_phi.csv  -- α и R² в λ₂ ~ (φ*)^α (контролирующий conductance);
  * phase.csv           -- по-спецификационные средние (logΔ, ρ*, IPR₂) + режим;
  * summary.txt         -- читаемые таблицы;
  * PNG-фигуры (если есть matplotlib): фазовая диаграмма, λ₂ vs φ*, λ₂(G) vs λ₂(core).

Запуск:
    python analyze_results.py --csv results/diagnostics.csv --out analysis/
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAVE_MPL = True
except Exception:  # pragma: no cover
    HAVE_MPL = False


# --- имена колонок (Unicode-заголовки из run_experiments.pretty_name) ---
C_ENS = "ensemble"
C_GAMMA = "γ"
C_DMIN = "d_min"
C_MU = "μ"
C_N = "n"
C_NTARGET = "n_target"
C_FRAC_CORE = "frac_2core"

L2_FULL = "λ₂_full"
L2_GIANT = "λ₂_giant"
L2_CORE = "λ₂_core"
PHI_GIANT = "φ*_giant"
PHI_BY_OBS = {"full": "φ*_full", "giant": "φ*_giant", "core": "φ*_core"}
LOG_DELTA = "logΔ_giant"
RHO_STAR = "ρ*_giant"
IPR_GIANT = "IPR₂_giant"
IPR_PI_GIANT = "IPR₂^π_giant"

# пороги классификатора режимов (эвристика, см. summary)
RHO_LOCAL = 0.25
IPR_EXCESS_LOCAL = 1.5
LOG_DELTA_MODULAR = 0.5


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return x


def load(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return [{k: fnum(v) for k, v in r.items()} for r in rows]


def is_num(x):
    return isinstance(x, float) and x == x and not math.isinf(x)


def _n(r):
    """Целевое n из спецификации (а не фактическое число вершин)."""
    v = r.get(C_NTARGET)
    return v if v is not None else r[C_N]


def fit_powerlaw(xs, ys):
    """Регрессия log y = slope·log x + c. Возвращает (slope, intercept, R², n)."""
    x = np.asarray([v for v in xs], dtype=float)
    y = np.asarray([v for v in ys], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if mask.sum() < 3:
        return float("nan"), float("nan"), float("nan"), int(mask.sum())
    lx, ly = np.log(x[mask]), np.log(y[mask])
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return slope, intercept, r2, int(mask.sum())


def mean(values):
    vals = [v for v in values if is_num(v)]
    return float(np.mean(vals)) if vals else float("nan")


def group_by(rows, keyfunc):
    g = defaultdict(list)
    for r in rows:
        g[keyfunc(r)].append(r)
    return g


def spec_key(r):
    return (r[C_ENS], r[C_GAMMA], r[C_DMIN], r[C_MU])


def classify(d_min, log_delta, rho, ipr_excess):
    """Эвристическая классификация механизма (документирована в summary)."""
    if not (is_num(log_delta) and is_num(rho)):
        return "unknown"
    if d_min >= 2:
        return "core-no-periphery"
    localized = is_num(ipr_excess) and ipr_excess >= IPR_EXCESS_LOCAL
    if rho < RHO_LOCAL:
        return "peripheral" if (log_delta >= 0 and localized) else "local-mixed"
    return "modular" if log_delta < LOG_DELTA_MODULAR else "modular-mixed"


def scaling_by_n(rows):
    out = []
    keys = sorted({(r[C_ENS], r[C_GAMMA], r[C_DMIN]) for r in rows})
    for ens, gamma, dmin in keys:
        sub = [r for r in rows if r[C_ENS] == ens and r[C_GAMMA] == gamma]
        by_n = defaultdict(list)
        for r in sub:
            by_n[_n(r)].append(r)
        ns = sorted(by_n)
        rec = {"ensemble": ens, "gamma": gamma, "d_min": dmin, "n_points": len(ns)}
        for name, col in (("full", L2_FULL), ("giant", L2_GIANT), ("core", L2_CORE)):
            xs = [n for n in ns]
            ys = [mean([r[col] for r in by_n[n]]) for n in ns]
            slope, _, r2, npts = fit_powerlaw(xs, ys)
            rec[f"alpha_{name}"] = -slope if slope == slope else float("nan")
            rec[f"r2_{name}"] = r2
        out.append(rec)
    return out


def scaling_by_phi(rows):
    """λ₂ ~ (φ*)^α: пулированная регрессия по всем реализациям спецификации."""
    out = []
    keys = sorted({(r[C_ENS], r[C_GAMMA], r[C_DMIN]) for r in rows})
    for ens, gamma, dmin in keys:
        sub = [r for r in rows if r[C_ENS] == ens and r[C_GAMMA] == gamma]
        for name, ycol in (("full", L2_FULL), ("giant", L2_GIANT), ("core", L2_CORE)):
            xcol = PHI_BY_OBS[name]
            xs = [r[xcol] for r in sub]
            ys = [r[ycol] for r in sub]
            slope, _, r2, npts = fit_powerlaw(xs, ys)
            out.append(
                {
                    "ensemble": ens,
                    "gamma": gamma,
                    "d_min": dmin,
                    "observable": name,
                    "alpha_phi": slope,
                    "r2_phi": r2,
                    "n_points": npts,
                    "in_cheeger_range": (is_num(slope) and 1.0 - 0.25 <= slope <= 2.0 + 0.25),
                }
            )
    return out


def phase_table(rows):
    out = []
    for (ens, gamma, dmin, mu), sub in sorted(
        group_by(rows, spec_key).items(), key=lambda kv: str(kv[0])
    ):
        for n in sorted({_n(r) for r in sub}):
            s = [r for r in sub if _n(r) == n]
            ld = mean([r[LOG_DELTA] for r in s])
            rho = mean([r[RHO_STAR] for r in s])
            ipr = mean([r[IPR_GIANT] for r in s])
            ipr_excess = ipr * n if is_num(ipr) else float("nan")
            out.append(
                {
                    "ensemble": ens,
                    "gamma": gamma,
                    "d_min": dmin,
                    "mu": mu,
                    "n": n,
                    "log_Delta": ld,
                    "rho_star": rho,
                    "ipr2_giant": ipr,
                    "ipr_excess": ipr_excess,
                    "frac_2core": mean([r[C_FRAC_CORE] for r in s]),
                    "regime": classify(dmin, ld, rho, ipr_excess),
                }
            )
    return out


def core_degree_comparison(rows):
    """d_min=2 vs d_min=3 при фиксированных (gamma, n): роль вершин степени 2.

    Оба ансамбля без листьев (core-no-periphery), поэтому различие в lambda2
    и alpha(n) не связано с периферией --- это чистый эффект малых степеней.
    """
    out = []
    gammas = sorted({r[C_GAMMA] for r in rows})
    ns = sorted({_n(r) for r in rows})
    for gamma in gammas:
        for n in ns:
            s2 = [r for r in rows if r[C_ENS] == "pl_dmin2"
                  and r[C_GAMMA] == gamma and _n(r) == n]
            s3 = [r for r in rows if r[C_ENS] == "pl_dmin3"
                  and r[C_GAMMA] == gamma and _n(r) == n]
            if not s2 or not s3:
                continue
            l2 = mean([r[L2_GIANT] for r in s2])
            l3 = mean([r[L2_GIANT] for r in s3])
            out.append(
                {
                    "gamma": gamma,
                    "n": n,
                    "lam2_dmin2": l2,
                    "lam2_dmin3": l3,
                    "ratio_d3_over_d2": (l3 / l2) if is_num(l2) and l2 > 0 else float("nan"),
                    "phi_dmin2": mean([r[PHI_GIANT] for r in s2]),
                    "phi_dmin3": mean([r[PHI_GIANT] for r in s3]),
                    "ipr_dmin2": mean([r[IPR_GIANT] for r in s2]),
                    "ipr_dmin3": mean([r[IPR_GIANT] for r in s3]),
                }
            )
    return out


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def _marker(ens):
    if "planted" in str(ens):
        return "D"
    if "dmin3" in str(ens):
        return "^"
    if "dmin2" in str(ens):
        return "s"
    return "o"


def make_plots(rows, phase, outdir):
    if not HAVE_MPL:
        print("[analyze] matplotlib недоступен -> фигуры пропущены")
        return []

    made = []
    # 1. фазовая диаграмма (log Δ, ρ*), цвет -- IPR₂·n
    fig, ax = plt.subplots(figsize=(7, 5))
    for ens in sorted({str(r[C_ENS]) for r in rows}):
        pts = [p for p in phase if p["ensemble"] == ens]
        if not pts:
            continue
        xs = [p["log_Delta"] for p in pts]
        ys = [p["rho_star"] for p in pts]
        cs = [p["ipr_excess"] for p in pts]
        sc = ax.scatter(xs, ys, c=cs, cmap="viridis", vmin=1, vmax=6,
                        marker=_marker(ens), s=45, alpha=0.8, label=ens)
    ax.set_xlabel("log Δ")
    ax.set_ylabel("ρ*")
    ax.set_title("Фазовая диаграмма (цвет: IPR₂·n, 1 = делокализовано)")
    ax.axhline(RHO_LOCAL, color="gray", ls="--", lw=0.8)
    ax.legend(fontsize=7, ncol=2)
    p = os.path.join(outdir, "phase_diagram.png")
    fig.tight_layout()
    fig.savefig(p, dpi=130)
    plt.close(fig)
    made.append(p)

    # 2. λ₂ vs φ* (log-log)
    fig, ax = plt.subplots(figsize=(6, 5))
    for ens in sorted({str(r[C_ENS]) for r in rows}):
        xs = [r[PHI_GIANT] for r in rows if str(r[C_ENS]) == ens]
        ys = [r[L2_GIANT] for r in rows if str(r[C_ENS]) == ens]
        pts = [(x, y) for x, y in zip(xs, ys) if is_num(x) and is_num(y) and x > 0 and y > 0]
        if not pts:
            continue
        px, py = zip(*pts)
        ax.scatter(px, py, s=10, alpha=0.5, marker=_marker(ens), label=ens)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("φ*  (giant)")
    ax.set_ylabel("λ₂ (giant)")
    ax.set_title("λ₂ против контролирующего conductance")
    ax.legend(fontsize=7, ncol=2)
    p = os.path.join(outdir, "lambda2_vs_phi.png")
    fig.tight_layout()
    fig.savefig(p, dpi=130)
    plt.close(fig)
    made.append(p)

    # 3. λ₂(G) vs λ₂(core)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    xs = [r[L2_FULL] for r in rows if is_num(r[L2_FULL]) and is_num(r[L2_CORE])]
    ys = [r[L2_CORE] for r in rows if is_num(r[L2_FULL]) and is_num(r[L2_CORE])]
    ax.scatter(xs, ys, s=10, alpha=0.5)
    lim = max([1e-3] + xs + ys) * 1.1
    ax.plot([0, lim], [0, lim], "r--", lw=0.8)
    ax.set_xlabel("λ₂(G)")
    ax.set_ylabel("λ₂(G_2core)")
    ax.set_title("Вклад удаления периферии")
    p = os.path.join(outdir, "lambda2_full_vs_core.png")
    fig.tight_layout()
    fig.savefig(p, dpi=130)
    plt.close(fig)
    made.append(p)
    return made


def main(argv=None):
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", default="results/diagnostics.csv")
    ap.add_argument("--out", default="analysis")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    rows = load(args.csv)
    print(f"[analyze] загружено {len(rows)} строк из {args.csv}")

    sn = scaling_by_n(rows)
    sph = scaling_by_phi(rows)
    ph = phase_table(rows)
    cc = core_degree_comparison(rows)
    write_csv(os.path.join(args.out, "scaling_by_n.csv"), sn)
    write_csv(os.path.join(args.out, "scaling_by_phi.csv"), sph)
    write_csv(os.path.join(args.out, "phase.csv"), ph)
    write_csv(os.path.join(args.out, "core_degree_comparison.csv"), cc)

    lines = []
    lines.append("=== α в λ₂ ~ n^(−α)  (R², точек n) ===")
    lines.append(f"{'ensemble':>14} {'γ':>5} {'dmin':>5} "
                 f"{'α(G)':>7}{'R²':>6} {'α(giant)':>9}{'R²':>6} "
                 f"{'α(core)':>8}{'R²':>6}")
    for r in sn:
        lines.append(
            f"{r['ensemble']:>14} {r['gamma']:>5} {r['d_min']:>5} "
            f"{r['alpha_full']:>7.3f}{r['r2_full']:>6.2f} "
            f"{r['alpha_giant']:>9.3f}{r['r2_giant']:>6.2f} "
            f"{r['alpha_core']:>8.3f}{r['r2_core']:>6.2f}"
        )

    lines.append("")
    lines.append("=== α в λ₂ ~ (φ*)^α  (giant; Cheeger ⇒ ожидаем α∈[1,2]) ===")
    lines.append(f"{'ensemble':>14} {'γ':>5} {'obs':>6} {'α_φ':>8}{'R²':>6} "
                 f"{'в [1,2]':>8}")
    for r in sph:
        lines.append(
            f"{r['ensemble']:>14} {r['gamma']:>5} {r['observable']:>6} "
            f"{r['alpha_phi']:>8.3f}{r['r2_phi']:>6.2f} "
            f"{str(r['in_cheeger_range']):>8}"
        )

    if cc:
        lines.append("")
        lines.append("=== d_min=2 vs d_min=3 (оба без листьев): эффект степени 2 ===")
        lines.append(f"{'γ':>5} {'n':>6} {'λ₂(dmin2)':>10} {'λ₂(dmin3)':>10} "
                     f"{'ratio':>6} {'φ*(d2)':>8} {'φ*(d3)':>8}")
        for r in cc:
            lines.append(
                f"{r['gamma']:>5} {r['n']:>6} {r['lam2_dmin2']:>10.4f} "
                f"{r['lam2_dmin3']:>10.4f} {r['ratio_d3_over_d2']:>6.2f} "
                f"{r['phi_dmin2']:>8.3f} {r['phi_dmin3']:>8.3f}"
            )

    lines.append("")
    lines.append("=== режимы (эвристика) ===")
    lines.append(f"peripheral: d_min=1, ρ*<{RHO_LOCAL}, logΔ≥0, IPR₂·n≥{IPR_EXCESS_LOCAL}")
    lines.append(f"modular:    ρ*≥{RHO_LOCAL}, logΔ<{LOG_DELTA_MODULAR}")
    counts = defaultdict(int)
    for p in ph:
        counts[p["regime"]] += 1
    for k in sorted(counts):
        lines.append(f"  {k:>18}: {counts[k]}")

    summary = "\n".join(lines)
    print(summary)
    with open(os.path.join(args.out, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write(summary + "\n")

    if not args.no_plot:
        made = make_plots(rows, ph, args.out)
        for p in made:
            print(f"[analyze] figure -> {p}")

    print(f"[analyze] выходные файлы в {args.out}/")


if __name__ == "__main__":
    main()
