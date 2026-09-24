"""Devil's-advocate numeric falsification tests on existing CSVs (read-only)."""
import csv, math, os, sys
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MERGED = os.path.join(ROOT, "results_merged", "diagnostics.csv")
ER = os.path.join(ROOT, "results_er", "diagnostics.csv")
TOY = os.path.join(ROOT, "results_toy", "diagnostics.csv")
OUT = os.path.join(HERE, "report_data.txt")

L2G = "\u03bb\u2082_giant"
L2C = "\u03bb\u2082_core"
L2F = "\u03bb\u2082_full"
GAM = "\u03b3"
DBAR = "d\u0304"

buf = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    buf.append(s)

def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")

def load(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

def fit(ns, ys):
    ns = np.asarray(ns, float); ys = np.asarray(ys, float)
    m = np.isfinite(ns) & np.isfinite(ys) & (ns > 0) & (ys > 0)
    if m.sum() < 2:
        return float("nan"), float("nan"), int(m.sum())
    lx, ly = np.log(ns[m]), np.log(ys[m])
    sl, ic = np.polyfit(lx, ly, 1)
    pred = sl * lx + ic
    sst = np.sum((ly - ly.mean()) ** 2)
    r2 = 1 - np.sum((ly - pred) ** 2) / sst if sst > 0 else float("nan")
    return -sl, r2, int(m.sum())

def gmean(v):
    v = [x for x in v if x == x and x > 0]
    return math.exp(sum(math.log(x) for x in v) / len(v)) if v else float("nan")

# ---------------- load merged ----------------
rows = load(MERGED)
P("loaded merged:", len(rows))

# group: (ensemble,gamma,dmin) -> n_target -> list rows
G = defaultdict(lambda: defaultdict(list))
for r in rows:
    g = fnum(r[GAM]); dm = fnum(r["d_min"]); nt = fnum(r["n_target"])
    G[(r["ensemble"], g, dm)][nt].append(r)

def series(sub, col):
    ns = sorted(sub)
    return ns, [np.mean([fnum(r[col]) for r in sub[n] if fnum(r[col]) == fnum(r[col])]) for n in ns]

def geo_series(sub, col):
    ns = sorted(sub)
    return ns, [gmean([fnum(r[col]) for r in sub[n]]) for n in ns]

P("\n================ A. PL ensembles: alpha(giant/core/full) & range stability ================")
for key in sorted(G, key=str):
    ens, g, dm = key
    if not ens.startswith("pl_"):
        continue
    sub = G[key]
    ns, yg = series(sub, L2G); _, yc = series(sub, L2C); _, yf = series(sub, L2F)
    P(f"\n-- {ens} gamma={g} d_min={dm}  (n points: {[int(x) for x in ns]})")
    for col, ys in (("giant", yg), ("core", yc), ("full", yf)):
        a, r2, k = fit(ns, ys)
        P(f"   alpha_{col:5s} full-range = {a:.4f}  R2={r2:.4f}  (pts={k})   means={[round(y,4) for y in ys]}")
    for thr in (1000, 3000, 10000):
        idx = [i for i, n in enumerate(ns) if n >= thr]
        if len(idx) >= 2:
            a, r2, k = fit([ns[i] for i in idx], [yg[i] for i in idx])
            ac, rc, _ = fit([ns[i] for i in idx], [yc[i] for i in idx])
            P(f"   n>={thr:6d}: alpha_giant={a:.4f} R2={r2:.4f} | alpha_core={ac:.4f} R2={rc:.4f} (pts={k})")
    # local alpha giant
    loc = []
    for i in range(len(ns) - 1):
        if yg[i] > 0 and yg[i + 1] > 0:
            loc.append(-math.log(yg[i + 1] / yg[i]) / math.log(ns[i + 1] / ns[i]))
        else:
            loc.append(float("nan"))
    P("   local alpha_giant per step:", [round(x, 4) for x in loc])
    # geomean and pooled-realization fits
    _, ygg = geo_series(sub, L2G)
    ag, rg, _ = fit(ns, ygg)
    P(f"   alpha_giant on GEOMEAN = {ag:.4f} R2={rg:.4f}")
    allx = [fnum(r["n_target"]) for n in sub for r in sub[n] if fnum(r[L2G]) == fnum(r[L2G]) and fnum(r[L2G]) > 0]
    ally = [fnum(r[L2G]) for n in sub for r in sub[n] if fnum(r[L2G]) == fnum(r[L2G]) and fnum(r[L2G]) > 0]
    ap, rp, kp = fit(allx, ally)
    P(f"   alpha_giant POOLED all realizations = {ap:.4f} R2={rp:.4f} (N={kp})")

P("\n================ B. jitter n_actual vs n_target, and cutoff d_max ================")
P(f"{'ensemble':16s} {'gamma':>5} {'nt':>7} {'n_mean':>9} {'n_sd':>7} {'ratio':>7} {'dmax_mean':>9} {'N1/...':>7} {'frac2core':>9} {'core_dbar':>9}")
for key in sorted(G, key=str):
    ens, g, dm = key
    if not ens.startswith("pl_"):
        continue
    for nt in sorted(G[key]):
        sub = G[key][nt]
        nn = [fnum(r["n"]) for r in sub]
        dmax = [fnum(r["d_max"]) for r in sub]
        n1 = [fnum(r["N\u2081"]) for r in sub]
        p2c = [fnum(r["frac_2core"]) for r in sub]
        ncor = [fnum(r["n_core"]) for r in sub]
        mcor = [fnum(r["m_core"]) for r in sub]
        cd = [2 * m / n for m, n in zip(mcor, ncor) if n and n > 0]
        P(f"{ens:16s} {g:5.1f} {int(nt):7d} {np.mean(nn):9.1f} {np.std(nn):7.2f} "
          f"{np.mean(nn)/nt:7.4f} {np.mean(dmax):9.2f} {np.mean(n1):7.2f} "
          f"{np.mean(p2c):9.4f} {np.mean(cd):9.3f}")

P("\n================ C. ER control ================")
er = load(ER)
EG = defaultdict(lambda: defaultdict(list))
for r in er:
    EG[fnum(r["c"])][fnum(r["n_target"])].append(r)
P(f"{'c':>4} {'obs':>6} {'alpha_full':>10} {'R2':>7} {'a>=1000':>9} {'a>=3000':>9} {'a>=10000':>9}  means")
for c in sorted(EG):
    for obs, col in (("giant", "giant_lambda2"), ("core", "core_lambda2")):
        sub = EG[c]
        ns = sorted(sub)
        ys = [np.mean([fnum(r[col]) for r in sub[n] if fnum(r[col]) == fnum(r[col])]) for n in ns]
        a, r2, _ = fit(ns, ys)
        res = []
        for thr in (1000, 3000, 10000):
            idx = [i for i, n in enumerate(ns) if n >= thr]
            aa, _, _ = fit([ns[i] for i in idx], [ys[i] for i in idx]) if len(idx) >= 2 else (float("nan"), 0, 0)
            res.append(aa)
        P(f"{c:4.0f} {obs:>6} {a:10.4f} {r2:7.3f} {res[0]:9.4f} {res[1]:9.4f} {res[2]:9.4f}  {[round(y,4) for y in ys]}")

P("\n================ D. PL vs ER matched on core mean degree ================")
# ER core mean degree: not stored; approximate via giant? use whole-graph c (ER core denser).
P("PL core mean degree (2*m_core/n_core) and giant frac per gamma, all n:")
for key in sorted(G, key=str):
    ens, g, dm = key
    if not ens.startswith("pl_"):
        continue
    allcd = []
    for nt in G[key]:
        sub = G[key][nt]
        for r in sub:
            nc = fnum(r["n_core"]); mc = fnum(r["m_core"])
            if nc and nc > 0:
                allcd.append(2 * mc / nc)
    P(f"  {ens} gamma={g}: core_dbar mean={np.mean(allcd):.3f}  (min {np.min(allcd):.3f}, max {np.max(allcd):.3f})")

P("\n================ E. leaves / low-degree in erased config model ================")
P("N1 (observed degree-1 count) and d_min_obs distributions (they should be 0 and 2 for pl_dmin2):")
for key in sorted(G, key=str):
    ens, g, dm = key
    if not ens.startswith("pl_"):
        continue
    n1all = [fnum(r["N\u2081"]) for nt in G[key] for r in G[key][nt]]
    dmobs = [fnum(r["d_min_obs"]) for nt in G[key] for r in G[key][nt]]
    n2all = [fnum(r["N\u2082"]) for nt in G[key] for r in G[key][nt]]
    nall = [fnum(r["n"]) for nt in G[key] for r in G[key][nt]]
    fracN1 = np.mean([a / b for a, b in zip(n1all, nall) if b]) if nall else float("nan")
    P(f"  {ens:16s} g={g}: mean N1={np.mean(n1all):8.2f}  mean N1/n={fracN1:.5f}  "
      f"d_min_obs: min={np.min(dmobs):.0f} mean={np.mean(dmobs):.3f}  "
      f"mean N2/n={np.mean([a/b for a,b in zip(n2all,nall) if b]):.4f}")

P("\n================ F. lambda2==0 / disconnected frequency ================")
for lab, path, gc, cc in (("merged", MERGED, L2G, L2C), ("er", ER, "giant_lambda2", "core_lambda2")):
    dat = rows if lab == "merged" else er
    zc = sum(1 for r in dat if fnum(r[cc]) == 0)
    zg = sum(1 for r in dat if fnum(r[gc]) == 0)
    nan_c = sum(1 for r in dat if fnum(r[cc]) != fnum(r[cc]))
    nan_g = sum(1 for r in dat if fnum(r[gc]) != fnum(r[gc]))
    P(f"  {lab}: core_l2==0 {zc}/{len(dat)} (nan {nan_c}) | giant_l2==0 {zg} (nan {nan_g})")
if "n_comp" in rows[0]:
    multi = sum(1 for r in rows if fnum(r["n_comp"]) > 1)
    P(f"  merged n_comp>1: {multi}/{len(rows)}  (full graph disconnected)")

P("\n================ G. mean vs median vs heavy-tail (per-realization spread) ================")
for key in sorted(G, key=str):
    ens, g, dm = key
    if ens != "pl_dmin2":
        continue
    P(f"\n  {ens} gamma={g} lambda2_giant per n: mean / median / geomean / sd / max/mean")
    ns = sorted(G[key])
    for n in ns:
        v = np.array([fnum(r[L2G]) for r in G[key][n] if fnum(r[L2G]) == fnum(r[L2G])])
        P(f"    n={int(n):7d}: mean={v.mean():.5f} med={np.median(v):.5f} geo={gmean(v):.5f} "
          f"sd={v.std():.5f} max/mean={v.max()/v.mean():.2f} frac<head/>")

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(buf) + "\n")
print("\nwritten", OUT)
