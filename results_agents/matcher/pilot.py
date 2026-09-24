import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from mlib import gen_er, gen_pl, metrics

n = 3000
reps = 3
cs = [1.5, 1.8, 2.0, 2.3, 2.6, 3.0, 3.5, 4.0, 4.5, 5.0]
print("=== ER (full-graph q) ===")
print(f"{'c':>5} {'d_bar_core':>11} {'q':>8} {'q_core':>8} {'gfrac':>7} {'cfrac':>7}")
for c in cs:
    rows = [metrics(gen_er(n, c, 1000 + i)) for i in range(reps)]
    dc = np.nanmean([r['d_bar_core'] for r in rows])
    q = np.nanmean([r['q'] for r in rows])
    qc = np.nanmean([r['q_core'] for r in rows])
    gf = np.nanmean([r['giant_frac'] for r in rows])
    cf = np.nanmean([r['core_frac'] for r in rows])
    print(f"{c:>5.2f} {dc:>11.3f} {q:>8.3f} {qc:>8.3f} {gf:>7.3f} {cf:>7.3f}")

print("=== PL d_min=2 (D grid), gamma ===")
for g in [2.1, 2.5, 3.0]:
    print(f"-- gamma={g} --")
    print(f"{'D':>5} {'d_bar_core':>11} {'q':>8} {'gfrac':>7} {'cfrac':>7}")
    for D in [10, 15, 20, 30, 40, 50, 70, 100]:
        rows = [metrics(gen_pl(n, g, D, 2000 + i)) for i in range(reps)]
        dc = np.nanmean([r['d_bar_core'] for r in rows])
        q = np.nanmean([r['q'] for r in rows])
        gf = np.nanmean([r['giant_frac'] for r in rows])
        cf = np.nanmean([r['core_frac'] for r in rows])
        print(f"{D:>5d} {dc:>11.3f} {q:>8.3f} {gf:>7.3f} {cf:>7.3f}")
