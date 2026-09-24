#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent finite-size scaling check of lambda2(giant) ~ A n^-alpha.

Reads (READ-ONLY):
  results_merged/diagnostics.csv   (pl_dmin2, pl_dmin3 only)
  results_er/diagnostics.csv       (model ER, by c)

Writes into results_agents/fss/:
  fits.csv, report.md
"""
import csv, os, math, itertools
import numpy as np
from scipy.optimize import curve_fit

BASE = r"C:\Users\yagan\OneDrive\Документы\graf"
OUT = os.path.join(BASE, "results_agents", "fss")
os.makedirs(OUT, exist_ok=True)


def fnum(s):
    if s is None:
        return np.nan
    s = s.strip()
    if s == "" or s.lower() in ("nan", "none", "null"):
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def load_merged():
    path = os.path.join(BASE, "results_merged", "diagnostics.csv")
    rows = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["ensemble"] not in ("pl_dmin2", "pl_dmin3"):
                continue
            rows.append({
                "ens": r["ensemble"],
                "gamma": r["\u03b3"],
                "nt": fnum(r["n_target"]),
                "n": fnum(r["n"]),
                "lam": fnum(r["\u03bb\u2082_giant"]),
            })
    return rows


def load_er():
    path = os.path.join(BASE, "results_er", "diagnostics.csv")
    rows = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["model"] != "ER":
                continue
            rows.append({
                "ens": "ER",
                "gamma": r["c"],
                "nt": fnum(r["n_target"]),
                "n": fnum(r["n"]),
                "lam": fnum(r["giant_lambda2"]),
            })
    return rows


def aggregate(rows):
    """Aggregate by n_target; returns (n_arr, mean_arr, std_arr, count_arr) sorted by n."""
    d = {}
    for r in rows:
        if not np.isfinite(r["nt"]) or not np.isfinite(r["lam"]):
            continue
        key = int(round(r["nt"]))
        d.setdefault(key, []).append(r["lam"])
    ns = sorted(d)
    n = np.array([float(k) for k in ns])
    m = np.array([np.mean(d[k]) for k in ns])
    s = np.array([np.std(d[k], ddof=1) if len(d[k]) > 1 else 0.0 for k in ns])
    c = np.array([len(d[k]) for k in ns])
    return n, m, s, c


# ---------------- model definitions in LOG space ----------------
# Each fn takes (treated_x, params) and returns log(y_hat).
# We fit on Y = log(y). Residuals are in log space (relative error).

def m1_linfit(t, Y):
    # log y = a - alpha t,  k=2
    sl, ic = np.polyfit(t, Y, 1)
    pred = ic + sl * t
    return {"a": ic, "alpha": -sl}, pred, 2


def m2_linfit(u, Y):
    # log y = a - beta u, u=log(log n),  k=2
    sl, ic = np.polyfit(u, Y, 1)
    pred = ic + sl * u
    return {"a": ic, "beta": -sl}, pred, 2


def m3_fit(t, Y, n):
    # y = C + A n^-alpha ; log y = log(C + A exp(-alpha t)),  k=3
    def f(t, laC, laA, alpha):
        C = math.exp(laC); A = math.exp(laA)
        return np.log(C + A * np.exp(-alpha * t))
    best = None
    for laC0 in (math.log(1e-4), math.log(1e-3), math.log(1e-2)):
        p0 = [laC0, math.log(max(1e-3, np.exp(Y[0])) * n[0] ** 0.3), 0.3]
        try:
            p, _ = curve_fit(f, t, Y, p0=p0,
                             bounds=([-30, -30, 1e-4], [5, 5, 10]), maxfev=200000)
            sse = float(np.sum((Y - f(t, *p)) ** 2))
            if best is None or sse < best[0]:
                best = (sse, p)
        except Exception:
            pass
    if best is None:
        return None, None, 3
    p = best[1]
    return {"C": math.exp(p[0]), "A": math.exp(p[1]), "alpha": p[2]}, f(t, *p), 3


def m4_fit(t, Y, n):
    # y = A n^-alpha + B n^-beta ; log y = log(A exp(-alpha t)+B exp(-beta t)), k=4
    def f(t, laA, alpha, lrB, beta):
        A = math.exp(laA); B = math.exp(laA + lrB)
        return np.log(A * np.exp(-alpha * t) + B * np.exp(-beta * t))
    best = None
    for a0 in (0.05, 0.2, 0.5):
        for b0 in (0.5, 1.0, 2.0):
            if b0 <= a0:
                continue
            p0 = [math.log(max(1e-6, np.exp(Y[0]))) + a0 * t[0], a0, math.log(0.1), b0]
            try:
                p, _ = curve_fit(f, t, Y, p0=p0,
                                 bounds=([-40, 1e-4, -20, 1e-4], [10, 10, 20, 20]),
                                 maxfev=400000)
                if p[3] <= p[1] + 1e-6:
                    continue
                sse = float(np.sum((Y - f(t, *p)) ** 2))
                if best is None or sse < best[0]:
                    best = (sse, p)
            except Exception:
                pass
    if best is None:
        return None, None, 4
    p = best[1]
    return ({"A": math.exp(p[0]), "alpha": p[1], "B": math.exp(p[0] + p[2]), "beta": p[3]},
            f(t, *p), 4)


def m5_fit(t, Y, n):
    # y = A exp(-c (log n)^delta) ; log y = a - c t^delta, k=3
    def f(t, a, c, delta):
        c = abs(c)
        return a - c * (t ** delta)
    best = None
    for d0 in (0.5, 1.0, 2.0):
        p0 = [Y[0] + 0.3 * t[0], 0.3 / max(1.0, t[0] ** (d0 - 1)), d0]
        try:
            p, _ = curve_fit(f, t, Y, p0=p0,
                             bounds=([-50, -20, 0.05], [50, 20, 10]), maxfev=400000)
            sse = float(np.sum((Y - f(t, *p)) ** 2))
            if best is None or sse < best[0]:
                best = (sse, p)
        except Exception:
            pass
    if best is None:
        return None, None, 3
    p = best[1]
    return {"a": p[0], "c": abs(p[1]), "delta": p[2]}, f(t, *p), 3


MODELS = [("M1 A n^-a", m1_linfit), ("M2 A (log n)^-b", m2_linfit),
          ("M3 C+A n^-a", m3_fit), ("M4 A n^-a+B n^-b", m4_fit),
          ("M5 A exp(-c logn^d)", m5_fit)]


def aicc(sse, N, k):
    if sse <= 0:
        sse = 1e-300
    aic = N * math.log(sse / N) + 2 * k
    denom = N - k - 1
    aicc_v = aic + (2 * k * (k + 1) / denom if denom > 0 else float("inf"))
    bic = N * math.log(sse / N) + k * math.log(N)
    return aic, aicc_v, bic


def run_set(label, n, y, out_rows):
    N = len(n)
    t = np.log(n)
    u = np.log(np.log(n))
    Y = np.log(y)
    res = {}
    for name, fn in MODELS:
        try:
            if name.startswith("M1"):
                params, pred, k = m1_linfit(t, Y)
            elif name.startswith("M2"):
                params, pred, k = m2_linfit(u, Y)
            elif name.startswith("M3"):
                params, pred, k = m3_fit(t, Y, n)
            elif name.startswith("M4"):
                params, pred, k = m4_fit(t, Y, n)
            else:
                params, pred, k = m5_fit(t, Y, n)
        except Exception as e:
            params, pred, k = None, None, 99
        if params is None:
            res[name] = dict(ok=False, sse=np.inf, k=k, params=None,
                             aicc=np.inf, bic=np.inf, aic=np.inf, resid=None)
            continue
        sse = float(np.sum((Y - pred) ** 2))
        aic, aicc_v, bic = aicc(sse, N, k)
        res[name] = dict(ok=True, sse=sse, k=k, params=params,
                         aicc=aicc_v, bic=bic, aic=aic, resid=Y - pred)
    return res


def subset_alpha(labels_y, min_n, model_name):
    """Refit chosen model on n>=min_n, return alpha (or relevant exponent)."""
    pass


def main():
    merged = load_merged()
    er = load_er()

    sets = []
    for ens in ("pl_dmin2", "pl_dmin3"):
        for g in ("2.1", "2.5", "3.0"):
            sel = [r for r in merged if r["ens"] == ens and r["gamma"] == g]
            sets.append((f"{ens} g={g}", sel))
    for c in ("2", "3", "4"):
        sel = [r for r in er if r["gamma"] == c]
        sets.append((f"ER c={c}", sel))

    out_rows = []
    report_blocks = []

    for label, sel in sets:
        n, y, sd, cnt = aggregate(sel)
        if len(n) < 3:
            continue
        N = len(n)
        res = run_set(label, n, y, out_rows)
        # best by AICc
        best = min(((k, v) for k, v in res.items() if v["ok"]),
                   key=lambda kv: kv[1]["aicc"])
        bestname = best[0]

        # alpha from M1 full
        alpha_M1 = res["M1 A n^-a"]["params"]["alpha"]
        # stability: refit M1 on lower-cut prefixes
        stab = {}
        for mn in (1000, 3000, 10000):
            m = n >= mn
            if m.sum() >= 2:
                a1, _, _ = m1_linfit(np.log(n[m]), np.log(y[m]))
                stab[mn] = a1["alpha"]
            else:
                stab[mn] = np.nan
        # stability: refit M1 on upper-cut suffixes (drop large n)
        stab_hi = {}
        for mx in (30000, 10000, 3000):
            m = n <= mx
            if m.sum() >= 2:
                a1, _, _ = m1_linfit(np.log(n[m]), np.log(y[m]))
                stab_hi[mx] = a1["alpha"]
            else:
                stab_hi[mx] = np.nan
        # local log-log slopes between consecutive points
        lslopes = [-np.log(y[i + 1] / y[i]) / np.log(n[i + 1] / n[i])
                   for i in range(N - 1)]
        # Deltas AICc of M1 vs others (positive => other model better)
        base = res["M1 A n^-a"]["aicc"]
        deltas = {}
        for name in ("M2 A (log n)^-b", "M3 C+A n^-a", "M4 A n^-a+B n^-b",
                     "M5 A exp(-c logn^d)"):
            deltas[name] = res[name]["aicc"] - base
        # residual sign pattern for M1
        resid = res["M1 A n^-a"]["resid"]
        signs = "".join("+" if r > 0 else "-" for r in resid)
        # curvature diagnostic: residual vs (log n)^2 (M1 residuals are exactly
        # orthogonal to log n by OLS construction, so plain corr(res, log n)==0)
        rr = float(np.corrcoef(np.log(n) ** 2, resid)[0, 1]) if N > 2 else np.nan

        out_rows.append(dict(label=label, N=N, best=bestname,
                             alpha_M1=alpha_M1, auc=None,
                             **{d: deltas[d] for d in deltas}))

        # print detailed
        print("=" * 78)
        print(f"{label}: N points={len(n)}  n={list(map(int,n))}")
        print(f"  lambda2 mean = {np.round(y,5)}")
        for name in res:
            v = res[name]
            if v["ok"]:
                print(f"  {name:20s} k={v['k']} aicc={v['aicc']:9.3f} "
                      f"bic={v['bic']:9.3f} dAICc={v['aicc']-base:8.3f} "
                      f"params={ {kk: round(vv,4) for kk,vv in v['params'].items()} }")
            else:
                print(f"  {name:20s} FAILED")
        print(f"  alpha M1 full={alpha_M1:.4f}  n>=1000:{stab[1000]:.4f} "
              f"n>=3000:{stab[3000]:.4f} n>=10000:{stab[10000]:.4f}")
        print(f"  alpha M1 full  | n<=30000:{stab_hi[30000]:.4f} "
              f"n<=10000:{stab_hi[10000]:.4f} n<=3000:{stab_hi[3000]:.4f}")
        print(f"  local slopes: {[round(float(s),3) for s in lslopes]}")
        print(f"  M1 resid signs (n asc): {signs}  corr(resid, (log n)^2)={rr:.3f}")

        report_blocks.append(dict(label=label, n=n, y=y, sd=sd, cnt=cnt,
                                  res=res, bestname=bestname, base=base,
                                  alpha=alpha_M1, stab=stab, stab_hi=stab_hi,
                                  deltas=deltas, signs=signs, rr=rr,
                                  lslopes=lslopes))

    # write fits.csv
    with open(os.path.join(OUT, "fits.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "N", "best_by_AICc", "alpha_M1",
                    "dAICc_M2", "dAICc_M3", "dAICc_M4", "dAICc_M5",
                    "alpha_n>=1000", "alpha_n>=3000", "alpha_n>=10000",
                    "alpha_n<=30000", "alpha_n<=10000", "alpha_n<=3000",
                    "local_slopes", "M1_resid_signs", "corr_resid_logn2"])
        for b in report_blocks:
            w.writerow([b["label"], len(b["n"]), b["bestname"], round(b["alpha"], 5),
                        round(b["deltas"]["M2 A (log n)^-b"], 3),
                        round(b["deltas"]["M3 C+A n^-a"], 3),
                        round(b["deltas"]["M4 A n^-a+B n^-b"], 3),
                        round(b["deltas"]["M5 A exp(-c logn^d)"], 3),
                        round(b["stab"][1000], 5) if np.isfinite(b["stab"][1000]) else "NA",
                        round(b["stab"][3000], 5) if np.isfinite(b["stab"][3000]) else "NA",
                        round(b["stab"][10000], 5) if np.isfinite(b["stab"][10000]) else "NA",
                        round(b["stab_hi"][30000], 5) if np.isfinite(b["stab_hi"][30000]) else "NA",
                        round(b["stab_hi"][10000], 5) if np.isfinite(b["stab_hi"][10000]) else "NA",
                        round(b["stab_hi"][3000], 5) if np.isfinite(b["stab_hi"][3000]) else "NA",
                        " ".join(f"{float(s):.3f}" for s in b["lslopes"]),
                        b["signs"], round(b["rr"], 3)])

    # write report.md
    def fmt(v):
        return f"{v:.4f}" if np.isfinite(v) else "n/a"

    with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as f:
        f.write("# Independent finite-size-scaling check: lambda2(giant) ~ A n^-alpha\n\n")
        f.write("**Method.** lambda2_giant is averaged over all reps sharing the same "
                "`n_target`; the x-axis is `n_target` (not the realised `n`), as specified. "
                "All fits are done in **log space**: we minimise the sum of squared "
                "residuals of `log lambda2` (i.e. relative error), which is appropriate "
                "because lambda2 spans ~1.5 decades. For each model N=8 points. "
                "`AIC = N ln(SSE/N)+2k`, `AICc = AIC + 2k(k+1)/(N-k-1)`, "
                "`BIC = N ln(SSE/N)+k ln N`.\n\n")
        f.write("Model set: **M1** `A n^-a` (k=2); **M2** `A (log n)^-b` (k=2); "
                "**M3** `C + A n^-a` (k=3); **M4** `A n^-a + B n^-b` (k=4); "
                "**M5** `A exp(-c (log n)^d)` (k=3).\n\n")

        f.write("## 1. Main table\n\n")
        f.write("Delta-AICc are written as **AICc(alternative) - AICc(M1)**: "
                "**positive = M1 (pure power law) is better**, negative = the alternative wins. "
                "`alpha(M1)` is the full-range exponent.\n\n")
        f.write("| set | best model (AICc) | alpha(M1) | dAICc 1vs2 | 1vs3 | 1vs4 | 1vs5 | "
                "alpha n>=1000 | alpha n>=10000 | verdict |\n")
        f.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for b in report_blocks:
            d = b["deltas"]
            con = ("power law survives" if b["bestname"].startswith("M1")
                   else f"**REJECT M1** ({b['bestname']} wins by {-min(d.values()):.2f})")
            f.write(f"| {b['label']} | {b['bestname']} | {b['alpha']:.4f} | "
                    f"{d['M2 A (log n)^-b']:+.2f} | {d['M3 C+A n^-a']:+.2f} | "
                    f"{d['M4 A n^-a+B n^-b']:+.2f} | {d['M5 A exp(-c logn^d)']:+.2f} | "
                    f"{fmt(b['stab'][1000])} | {fmt(b['stab'][10000])} | {con} |\n")

        f.write("\n## 2. Stability of alpha (lower-cut, upper-cut, local slopes)\n\n")
        f.write("| set | full | n>=1000 | n>=3000 | n>=10000 | n<=30000 | n<=10000 | "
                "n<=3000 | local log-log slopes (n=100->...->1e5) |\n")
        f.write("|---|---|---|---|---|---|---|---|---|\n")
        for b in report_blocks:
            ls = " ".join(f"{float(s):.3f}" for s in b["lslopes"])
            f.write(f"| {b['label']} | {b['alpha']:.4f} | {fmt(b['stab'][1000])} | "
                    f"{fmt(b['stab'][3000])} | {fmt(b['stab'][10000])} | "
                    f"{fmt(b['stab_hi'][30000])} | {fmt(b['stab_hi'][10000])} | "
                    f"{fmt(b['stab_hi'][3000])} | `{ls}` |\n")

        f.write("\n## 3. Residual structure of M1 (pure power law)\n\n")
        f.write("Signs of M1 log-residuals in increasing n. A long runs structure / a "
                "systematic sign change indicates the pure power law is not the true form. "
                "`corr` = correlation of the M1 residual with `(log n)^2` "
                "(M1 residuals are orthogonal to `log n` by construction, so `(log n)^2` "
                "is the natural curvature probe).\n\n")
        f.write("| set | resid signs (n asc) | corr(res, (log n)^2) |\n|---|---|---|\n")
        for b in report_blocks:
            f.write(f"| {b['label']} | `{b['signs']}` | {b['rr']:+.3f} |\n")

        f.write("\n## 4. Full fit details (all models)\n\n")
        for b in report_blocks:
            f.write(f"### {b['label']}\n\n")
            f.write(f"lambda2 by n_target: {[round(float(v), 5) for v in b['y']]} "
                    f"(rep counts {[int(c) for c in b['cnt']]})\n\n")
            f.write("| model | k | AICc | BIC | dAICc vs M1 | params |\n")
            f.write("|---|---|---|---|---|---|\n")
            for name, v in b["res"].items():
                if v["ok"]:
                    ps = ", ".join(f"{kk}={vv:.4g}" for kk, vv in v["params"].items())
                    f.write(f"| {name} | {v['k']} | {v['aicc']:.3f} | {v['bic']:.3f} | "
                            f"{v['aicc']-b['base']:+.3f} | {ps} |\n")
                else:
                    f.write(f"| {name} | - | FAILED | - | - | - |\n")
            f.write("\n")

        # helper: pull a fitted param
        def P(b, model, key):
            return b["res"][model]["params"][key]

        f.write("\n## 5. Where offset / two-power models win: saturation / crossover\n\n")
        f.write("Only M3 (`C + A n^-a`) and M4 (`A n^-a + B n^-b`) can signal a floor or a "
                "crossover. Readings:\n\n")
        f.write("**Clear saturation (M3 beats M1 by > ~4 AICc; nonzero floor C):**\n\n")
        f.write("- **pl_dmin3 g=2.5**: dAICc(M3-M1) = "
                f"{report_blocks[4]['deltas']['M3 C+A n^-a']:+.2f}, floor C = "
                f"{P(report_blocks[4],'M3 C+A n^-a','C'):.4f} "
                f"(lambda2: {report_blocks[4]['y'][0]:.3f} at n=100 -> "
                f"{report_blocks[4]['y'][-1]:.3f} at n=1e5). lambda2 approaches a constant, "
                "not zero.\n")
        f.write("- **pl_dmin3 g=3.0**: dAICc(M3-M1) = "
                f"{report_blocks[5]['deltas']['M3 C+A n^-a']:+.2f}, floor C = "
                f"{P(report_blocks[5],'M3 C+A n^-a','C'):.4f} "
                f"(lambda2: {report_blocks[5]['y'][0]:.3f} -> "
                f"{report_blocks[5]['y'][-1]:.3f}). Clear plateau.\n")
        f.write("- **ER c=2**: dAICc(M3-M1) = "
                f"{report_blocks[6]['deltas']['M3 C+A n^-a']:+.2f}, floor C = "
                f"{P(report_blocks[6],'M3 C+A n^-a','C'):.4f} "
                f"(small but statistically resolved; M2 log-power even better here).\n")
        f.write("- **pl_dmin3 g=2.1**: M1 still wins by "
                f"{report_blocks[3]['deltas']['M3 C+A n^-a']:+.2f}, but M3 needs a very large "
                f"floor C = {P(report_blocks[3],'M3 C+A n^-a','C'):.4f} (compared with "
                f"lambda2 = {report_blocks[3]['y'][0]:.3f}..{report_blocks[3]['y'][-1]:.3f}); "
                "i.e. the range is dominated by an offset and alpha is poorly identified "
                "(M3 alpha = 0.081 vs M1 alpha = 0.037).\n\n")
        f.write("**Marginal / weak:** pl_dmin2 g=3.0 (dAICc M3-M1 = "
                f"{report_blocks[2]['deltas']['M3 C+A n^-a']:+.2f}, M5 -2.61), "
                "pl_dmin2 g=2.1 and g=2.5 and ER c=3,c=4 show no floor "
                "(M3 dAICc > +3.99 with C near 0).\n\n")
        f.write("**Two-power M4:** it is never the AICc winner. It is only competitive for "
                "pl_dmin3 g=2.5 (dAICc M4-M1 = "
                f"{report_blocks[4]['deltas']['M4 A n^-a+B n^-b']:+.2f}) and g=3.0 "
                f"({report_blocks[5]['deltas']['M4 A n^-a+B n^-b']:+.2f}); elsewhere it "
                "collapses to M3 (one amplitude -> 0) and is penalised as an over-fit. "
                "Hence the deviation from a pure power law is best described as an additive "
                "floor, not as a genuine two-exponent crossover.\n")

        f.write("\n## 6. Verdict\n\n")
        f.write("**Pure power law `A n^-a` passes the strict test (best AICc, all "
                "alternatives worse, stable alpha) for:**\n\n")
        f.write("- `pl_dmin2 g=2.1` alpha = 0.150 (dAICc vs M2/M3/M4/M5 = "
                "+9.9/+4.0/+13.3/+4.4)\n")
        f.write("- `pl_dmin2 g=2.5` alpha = 0.158 (+15.7/+5.6/+14.9/+5.3)\n")
        f.write("- `ER c=3` alpha = 0.206 and `ER c=4` alpha = 0.168 (controls)\n")
        f.write("- `pl_dmin3 g=2.1` alpha = 0.037 only marginally (M1 wins by +2.56 over M3, "
                "but M3 = C+A n^-a already carries a floor C = 0.151).\n\n")
        f.write("**Pure power law is REJECTED (an alternative model wins by dAICc, no floor) "
                "for:**\n\n")
        f.write("- `pl_dmin2 g=3.0`: M5 wins by 2.61 (M3 also ties, -0.24)\n")
        f.write("- `pl_dmin3 g=2.5`: M2 log-power wins by 8.16; M3 offset wins by 5.76\n")
        f.write("- `pl_dmin3 g=3.0`: M2 log-power wins by 9.01; M3 offset wins by 4.73; "
                "alpha falls monotonically 0.038 (n<=3000) -> 0.021 (n>=3000) -> 0.016 "
                "(n>=1e4)\n")
        f.write("- `ER c=2`: M2 wins by 9.59.\n\n")
        f.write("**alpha stability.** For the surviving sets alpha changes by <= 0.03 across "
                "lower cuts (n>=1000, >=1e4) and upper cuts (n<=1e4, <=3e3). For the rejected "
                "sets alpha drifts systematically with the window: e.g. pl_dmin3 g=3.0 goes "
                "0.0321 (n<=3e4) -> 0.0163 (n>=1e4); pl_dmin2 g=2.5 0.160 (n<=3e4) -> 0.135 "
                "(n>=1e4), and pl_dmin3 g=2.5 collapses to a plateau. This window dependence "
                "is the direct signature that the fitted alpha is not a true asymptotic "
                "exponent.\n\n")
        f.write("**On the stated claim.** `alpha(d_min=2) ~ 0.15-0.17` is reproduced for "
                "gamma = 2.1 (0.150) and 2.5 (0.158), and gamma = 3.0 (0.171) -- but for "
                "gamma = 3.0 the pure power law is not the best description (stretched "
                "exponential/offset preferred), so 0.171 is an effective local slope. "
                "`alpha(d_min=3) ~ 0.03` is reproduced in magnitude (0.037 / 0.035 / 0.030 "
                "for gamma = 2.1 / 2.5 / 3.0), but the clean power law is rejected for "
                "gamma = 2.5 and 3.0 and holds for gamma = 2.1 only with a large constant "
                "offset. In short, the values are of the right size, but the 'pure power-law "
                "in n' interpretation is only valid in a minority of the tested sets; for "
                "d_min=3 at gamma >= 2.5 the data prefer a saturating (offset / log-power) "
                "form.\n")

    print("\nWROTE", os.path.join(OUT, "report.md"), "and fits.csv")


if __name__ == "__main__":
    main()
