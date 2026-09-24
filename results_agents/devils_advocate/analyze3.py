"""Round 3: jackknife, IPR/gap localization, explicit alpha ratios."""
import csv, math, os, sys
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MERGED = os.path.join(ROOT, "results_merged", "diagnostics.csv")
ER = os.path.join(ROOT, "results_er", "diagnostics.csv")
OUT = os.path.join(HERE, "report_data3.txt")
L2G="\u03bb\u2082_giant"; L2C="\u03bb\u2082_core"; GAM="\u03b3"
buf=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); buf.append(s)
def fnum(x):
    try: return float(x)
    except: return float("nan")
def load(p):
    with open(p,encoding="utf-8") as f: return list(csv.DictReader(f))
def fit(ns,ys):
    ns=np.asarray(ns,float); ys=np.asarray(ys,float)
    m=np.isfinite(ns)&np.isfinite(ys)&(ns>0)&(ys>0)
    lx,ly=np.log(ns[m]),np.log(ys[m]); sl,ic=np.polyfit(lx,ly,1)
    pred=sl*lx+ic; sst=np.sum((ly-ly.mean())**2)
    return -sl, (1-np.sum((ly-pred)**2)/sst if sst>0 else float("nan"))
rows=load(MERGED); er=load(ER)
G=defaultdict(lambda: defaultdict(list))
for r in rows: G[(r["ensemble"],fnum(r[GAM]))][fnum(r["n_target"])].append(r)
EG=defaultdict(lambda: defaultdict(list))
for r in er: EG[fnum(r["c"])][fnum(r["n_target"])].append(r)
def mean_col(sub,n,col): return np.mean([fnum(r[col]) for r in sub[n] if fnum(r[col])==fnum(r[col])])

P("================ K. jackknife alpha (drop one n) ================")
for ens in ("pl_dmin2","pl_dmin3"):
    for g in sorted({k[1] for k in G if k[0]==ens}):
        sub=G[(ens,g)]; ns=sorted(sub); ys=[mean_col(sub,n,L2G) for n in ns]
        full=fit(ns,ys)[0]; al=[]
        for i in range(len(ns)):
            keep=[j for j in range(len(ns)) if j!=i]
            al.append(fit([ns[j] for j in keep],[ys[j] for j in keep])[0])
        P(f" {ens} g={g}: full={full:.4f} jack[min={min(al):.4f} max={max(al):.4f}] "
          f"drop100={al[0]:.4f} drop1e5={al[-1]:.4f}")
for c in sorted(EG):
    for obs,col in (("giant","giant_lambda2"),("core","core_lambda2")):
        sub=EG[c]; ns=sorted(sub); ys=[mean_col(sub,n,col) for n in ns]
        full=fit(ns,ys)[0]; al=[]
        for i in range(len(ns)):
            keep=[j for j in range(len(ns)) if j!=i]
            al.append(fit([ns[j] for j in keep],[ys[j] for j in keep])[0])
        P(f" ER c={int(c)} {obs:5s}: full={full:.4f} jack[min={min(al):.4f} max={max(al):.4f}]")

P("\n================ L. localization diagnostics (IPR2_giant * n, gap23) ================")
for ens in ("pl_dmin2","pl_dmin3"):
    for g in sorted({k[1] for k in G if k[0]==ens}):
        sub=G[(ens,g)]
        P(f" {ens} g={g}:")
        for n in sorted(sub):
            ipr=[fnum(r["IPR\u2082_giant"]) for r in sub[n] if fnum(r["IPR\u2082_giant"])==fnum(r["IPR\u2082_giant"])]
            gap=[fnum(r["(\u03bb\u2083\u2212\u03bb\u2082)_giant"]) for r in sub[n] if fnum(r["(\u03bb\u2083\u2212\u03bb\u2082)_giant"])==fnum(r["(\u03bb\u2083\u2212\u03bb\u2082)_giant"])]
            frac_small=np.mean([1 if x<1e-3 else 0 for x in gap]) if gap else float("nan")
            P(f"   n={int(n):7d}: IPR2*n mean={np.mean([x*n for x in ipr]):.3f}  gap23 mean={np.mean(gap):.5f}  frac(gap23<1e-3)={frac_small:.3f}")

P("\n================ M. explicit alpha ratios at ~matched core dbar ================")
# PL dmin2 core dbar (mean over n) and alpha_core
P("PL dmin2: alpha_giant / alpha_core, core dbar (~n-independent):")
for g in sorted({k[1] for k in G if k[0]=="pl_dmin2"}):
    sub=G[("pl_dmin2",g)]; ns=sorted(sub)
    ag=fit(ns,[mean_col(sub,n,L2G) for n in ns])[0]
    ac=fit(ns,[mean_col(sub,n,L2C) for n in ns])[0]
    cd=np.mean([2*fnum(r["m_core"])/fnum(r["n_core"]) for n in sub for r in sub[n] if fnum(r["n_core"])])
    P(f"  g={g}: core_dbar={cd:.2f} alpha_giant={ag:.4f} alpha_core={ac:.4f}")
P("ER: alpha_giant vs alpha_core (core dbar from LITERATURE: c2=2.68 c3=3.44 c4=4.26):")
cdmap={2:2.68,3:3.44,4:4.26}
for c in sorted(EG):
    sub=EG[c]; ns=sorted(sub)
    ag=fit(ns,[mean_col(sub,n,"giant_lambda2") for n in ns])[0]
    ac=fit(ns,[mean_col(sub,n,"core_lambda2") for n in ns])[0]
    P(f"  c={int(c)} core_dbar~{cdmap[int(c)]}: alpha_giant={ag:.4f} alpha_core={ac:.4f}")

with open(OUT,"w",encoding="utf-8") as f: f.write("\n".join(buf)+"\n")
print("written",OUT)
