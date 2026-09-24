"""Devil's advocate round 2: merge consistency, model selection, bootstrap CI, solver check."""
import csv, math, os, sys
from collections import defaultdict
import numpy as np
from scipy.optimize import curve_fit

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MERGED = os.path.join(ROOT, "results_merged", "diagnostics.csv")
ER = os.path.join(ROOT, "results_er", "diagnostics.csv")
OUT = os.path.join(HERE, "report_data2.txt")

L2G = "\u03bb\u2082_giant"; L2C = "\u03bb\u2082_core"; GAM = "\u03b3"
buf = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); buf.append(s)
def fnum(x):
    try: return float(x)
    except: return float("nan")
def load(p):
    with open(p, encoding="utf-8") as f: return list(csv.DictReader(f))

def fit(ns, ys):
    ns = np.asarray(ns, float); ys = np.asarray(ys, float)
    m = np.isfinite(ns) & np.isfinite(ys) & (ns > 0) & (ys > 0)
    lx, ly = np.log(ns[m]), np.log(ys[m])
    sl, ic = np.polyfit(lx, ly, 1)
    pred = sl * lx + ic
    sst = np.sum((ly - ly.mean()) ** 2)
    r2 = 1 - np.sum((ly - pred) ** 2) / sst if sst > 0 else float("nan")
    return -sl, ic, r2

def aicc(ssr, k, n):
    if not np.isfinite(ssr) or ssr <= 0 or n - k - 1 <= 0: return float("inf")
    return n * math.log(ssr / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)

def ssr_pl(ns, ys):
    a, ic, _ = fit(ns, ys); pred = np.exp(ic) * ns ** (-a)
    return float(np.sum((ys - pred) ** 2))
def ssr_explog(ns, ys):
    def f(n, A, c, b): return A * np.exp(-c * np.log(n) ** b)
    try:
        p, _ = curve_fit(f, ns, ys, p0=[ys[0], 0.1, 1.0],
                         bounds=([0, 0, 0.1], [np.inf, 10, 5]), maxfev=40000)
        return float(np.sum((ys - f(ns, *p)) ** 2)), p
    except Exception: return float("inf"), None

rows = load(MERGED); er = load(ER)
G = defaultdict(lambda: defaultdict(list))
for r in rows:
    G[(r["ensemble"], fnum(r[GAM]))][fnum(r["n_target"])].append(r)
EG = defaultdict(lambda: defaultdict(list))
for r in er:
    EG[fnum(r["c"])][fnum(r["n_target"])].append(r)

P("================ H. v2-only (n<=1e4) vs full (incl n=3e4,1e5) & extrapolation ================")
for ens in ("pl_dmin2", "pl_dmin3"):
    for g in sorted({k[1] for k in G if k[0] == ens}):
        sub = G[(ens, g)]
        ns = sorted(sub)
        gm = lambda c: [np.mean([fnum(r[c]) for r in sub[n] if fnum(r[c]) == fnum(r[c])]) for n in ns]
        yg = gm(L2G)
        small = [i for i, n in enumerate(ns) if n <= 10000]
        a_s, ic_s, r_s = fit([ns[i] for i in small], [yg[i] for i in small])
        a_f, ic_f, r_f = fit(ns, yg)
        P(f" {ens} g={g}: alpha v2-only(<=1e4)={a_s:.4f} R2={r_s:.4f} | full={a_f:.4f} R2={r_f:.4f}")
        for i, n in enumerate(ns):
            if n > 10000:
                pred = math.exp(ic_s) * n ** (-a_s)
                P(f"      n={int(n):7d}: observed={yg[i]:.5f} extrap_v2={pred:.5f} ratio_obs/pred={yg[i]/pred:.4f}")

P("\n================ I. model selection: power law vs stretched-exp (AICc) ================")
P(f"{'series':28s} {'alpha_pl':>9} {'AICc_pl':>9} {'AICc_explog':>11} {'winner':>10} beta3")
for ens in ("pl_dmin2", "pl_dmin3"):
    for g in sorted({k[1] for k in G if k[0] == ens}):
        sub = G[(ens, g)]; ns = sorted(sub)
        yg = [np.mean([fnum(r[L2G]) for r in sub[n] if fnum(r[L2G]) == fnum(r[L2G])]) for n in ns]
        n = len(ns)
        a = fit(ns, yg)[0]
        s1 = ssr_pl(np.asarray(ns, float), np.asarray(yg)); s3, p3 = ssr_explog(np.asarray(ns, float), np.asarray(yg))
        P(f"{ens+' g'+str(g):28s} {a:9.4f} {aicc(s1,2,n):9.2f} {aicc(s3,3,n):11.2f} "
          f"{'explog' if aicc(s3,3,n)<aicc(s1,2,n) else 'power':>10} {None if p3 is None else round(float(p3[2]),3)}")
for c in sorted(EG):
    for obs, col in (("giant", "giant_lambda2"), ("core", "core_lambda2")):
        sub = EG[c]; ns = sorted(sub)
        ys = [np.mean([fnum(r[col]) for r in sub[n] if fnum(r[col]) == fnum(r[col])]) for n in ns]
        n = len(ns); a = fit(ns, ys)[0]
        s1 = ssr_pl(np.asarray(ns, float), np.asarray(ys)); s3, p3 = ssr_explog(np.asarray(ns, float), np.asarray(ys))
        P(f"{'ER c='+str(int(c))+' '+obs:28s} {a:9.4f} {aicc(s1,2,n):9.2f} {aicc(s3,3,n):11.2f} "
          f"{'explog' if aicc(s3,3,n)<aicc(s1,2,n) else 'power':>10} {None if p3 is None else round(float(p3[2]),3)}")

P("\n================ J. bootstrap 95% CI of alpha (resample reps within n, B=800) ================")
rng = np.random.default_rng(0)
def boot(by_n, ns, col, B=800):
    al = []
    for _ in range(B):
        means = []
        for n in ns:
            v = np.array([fnum(r[col]) for r in by_n[n] if fnum(r[col]) == fnum(r[col])])
            means.append(v[rng.integers(0, len(v), len(v))].mean())
        a = fit(ns, means)[0]
        if a == a: al.append(a)
    return np.percentile(al, [2.5, 97.5])
for ens in ("pl_dmin2", "pl_dmin3"):
    for g in sorted({k[1] for k in G if k[0] == ens}):
        sub = G[(ens, g)]; ns = sorted(sub)
        ci = boot(sub, ns, L2G)
        P(f" {ens} g={g} alpha_giant = {fit(ns,[np.mean([fnum(r[L2G]) for r in sub[n] if fnum(r[L2G])==fnum(r[L2G])]) for n in ns])[0]:.4f} CI={ci.round(4)}")
for c in sorted(EG):
    for obs, col in (("giant", "giant_lambda2"), ("core", "core_lambda2")):
        sub = EG[c]; ns = sorted(sub)
        ci = boot(sub, ns, col)
        P(f" ER c={int(c)} {obs:5s} alpha = {fit(ns,[np.mean([fnum(r[col]) for r in sub[n] if fnum(r[col])==fnum(r[col])]) for n in ns])[0]:.4f} CI={ci.round(4)}")

with open(OUT, "w", encoding="utf-8") as f: f.write("\n".join(buf) + "\n")
print("written", OUT)
