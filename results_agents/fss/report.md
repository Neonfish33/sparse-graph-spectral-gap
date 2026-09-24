# Independent finite-size-scaling check: lambda2(giant) ~ A n^-alpha

**Method.** lambda2_giant is averaged over all reps sharing the same `n_target`; the x-axis is `n_target` (not the realised `n`), as specified. All fits are done in **log space**: we minimise the sum of squared residuals of `log lambda2` (i.e. relative error), which is appropriate because lambda2 spans ~1.5 decades. For each model N=8 points. `AIC = N ln(SSE/N)+2k`, `AICc = AIC + 2k(k+1)/(N-k-1)`, `BIC = N ln(SSE/N)+k ln N`.

Model set: **M1** `A n^-a` (k=2); **M2** `A (log n)^-b` (k=2); **M3** `C + A n^-a` (k=3); **M4** `A n^-a + B n^-b` (k=4); **M5** `A exp(-c (log n)^d)` (k=3).

## 1. Main table

Delta-AICc are written as **AICc(alternative) - AICc(M1)**: **positive = M1 (pure power law) is better**, negative = the alternative wins. `alpha(M1)` is the full-range exponent.

| set | best model (AICc) | alpha(M1) | dAICc 1vs2 | 1vs3 | 1vs4 | 1vs5 | alpha n>=1000 | alpha n>=10000 | verdict |
|---|---|---|---|---|---|---|---|---|---|
| pl_dmin2 g=2.1 | M1 A n^-a | 0.1496 | +9.89 | +3.99 | +13.32 | +4.40 | 0.1400 | 0.1344 | power law survives |
| pl_dmin2 g=2.5 | M1 A n^-a | 0.1580 | +15.73 | +5.60 | +14.93 | +5.26 | 0.1615 | 0.1353 | power law survives |
| pl_dmin2 g=3.0 | M5 A exp(-c logn^d) | 0.1708 | +3.25 | -0.24 | +5.02 | -2.61 | 0.1584 | 0.1690 | **REJECT M1** (M5 A exp(-c logn^d) wins by 2.61) |
| pl_dmin3 g=2.1 | M1 A n^-a | 0.0374 | +5.92 | +2.56 | +11.90 | +3.00 | 0.0347 | 0.0297 | power law survives |
| pl_dmin3 g=2.5 | M2 A (log n)^-b | 0.0345 | -8.16 | -5.76 | +0.04 | -7.48 | 0.0308 | 0.0323 | **REJECT M1** (M2 A (log n)^-b wins by 8.16) |
| pl_dmin3 g=3.0 | M2 A (log n)^-b | 0.0295 | -9.01 | -4.73 | +1.96 | -2.98 | 0.0250 | 0.0163 | **REJECT M1** (M2 A (log n)^-b wins by 9.01) |
| ER c=2 | M2 A (log n)^-b | 0.2697 | -9.59 | -5.75 | +2.49 | -6.89 | 0.2423 | 0.2279 | **REJECT M1** (M2 A (log n)^-b wins by 9.59) |
| ER c=3 | M1 A n^-a | 0.2059 | +2.68 | +5.47 | +6.57 | +5.15 | 0.2047 | 0.2087 | power law survives |
| ER c=4 | M1 A n^-a | 0.1675 | +8.40 | +4.75 | +14.08 | +5.25 | 0.1582 | 0.1492 | power law survives |

## 2. Stability of alpha (lower-cut, upper-cut, local slopes)

| set | full | n>=1000 | n>=3000 | n>=10000 | n<=30000 | n<=10000 | n<=3000 | local log-log slopes (n=100->...->1e5) |
|---|---|---|---|---|---|---|---|---|
| pl_dmin2 g=2.1 | 0.1496 | 0.1400 | 0.1413 | 0.1344 | 0.1520 | 0.1566 | 0.1568 | `0.066 0.184 0.213 0.127 0.160 0.122 0.146` |
| pl_dmin2 g=2.5 | 0.1580 | 0.1615 | 0.1563 | 0.1353 | 0.1608 | 0.1597 | 0.1465 | `0.140 0.125 0.161 0.162 0.201 0.140 0.132` |
| pl_dmin2 g=3.0 | 0.1708 | 0.1584 | 0.1533 | 0.1690 | 0.1727 | 0.1802 | 0.1956 | `0.210 0.200 0.180 0.196 0.123 0.156 0.180` |
| pl_dmin3 g=2.1 | 0.0374 | 0.0347 | 0.0359 | 0.0297 | 0.0389 | 0.0394 | 0.0383 | `0.036 0.037 0.058 0.023 0.049 0.033 0.027` |
| pl_dmin3 g=2.5 | 0.0345 | 0.0308 | 0.0299 | 0.0323 | 0.0359 | 0.0369 | 0.0413 | `0.048 0.045 0.037 0.037 0.023 0.039 0.026` |
| pl_dmin3 g=3.0 | 0.0295 | 0.0250 | 0.0209 | 0.0163 | 0.0321 | 0.0359 | 0.0380 | `0.079 0.040 0.009 0.038 0.033 0.011 0.021` |
| ER c=2 | 0.2697 | 0.2423 | 0.2234 | 0.2279 | 0.2858 | 0.2925 | 0.3221 | `0.333 0.398 0.208 0.328 0.189 0.305 0.160` |
| ER c=3 | 0.2059 | 0.2047 | 0.2129 | 0.2087 | 0.2060 | 0.2067 | 0.2109 | `0.561 0.073 0.221 0.165 0.222 0.210 0.207` |
| ER c=4 | 0.1675 | 0.1582 | 0.1588 | 0.1492 | 0.1740 | 0.1698 | 0.1716 | `0.057 0.191 0.248 0.146 0.162 0.206 0.099` |

## 3. Residual structure of M1 (pure power law)

Signs of M1 log-residuals in increasing n. A long runs structure / a systematic sign change indicates the pure power law is not the true form. `corr` = correlation of the M1 residual with `(log n)^2` (M1 residuals are orthogonal to `log n` by construction, so `(log n)^2` is the natural curvature probe).

| set | resid signs (n asc) | corr(res, (log n)^2) |
|---|---|---|
| pl_dmin2 g=2.1 | `-++---++` | +0.049 |
| pl_dmin2 g=2.5 | `--+++--+` | -0.017 |
| pl_dmin2 g=3.0 | `++---+++` | +0.090 |
| pl_dmin3 g=2.1 | `+++-+--+` | +0.067 |
| pl_dmin3 g=2.5 | `++---+-+` | +0.103 |
| pl_dmin3 g=3.0 | `+-----++` | +0.099 |
| ER c=2 | `++---+-+` | +0.106 |
| ER c=3 | `+---++++` | +0.022 |
| ER c=4 | `-++----+` | +0.032 |

## 4. Full fit details (all models)

### pl_dmin2 g=2.1

lambda2 by n_target: [0.19097, 0.18244, 0.15411, 0.13298, 0.11566, 0.09535, 0.0834, 0.06999] (rep counts [50, 50, 50, 50, 50, 50, 100, 100])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -54.304 | -56.545 | +0.000 | a=-0.9513, alpha=0.1496 |
| M2 A (log n)^-b | 2 | -44.417 | -46.658 | +9.887 | a=0.1449, beta=1.126 |
| M3 C+A n^-a | 3 | -50.314 | -56.076 | +3.990 | C=0.01829, A=0.406, alpha=0.1788 |
| M4 A n^-a+B n^-b | 4 | -40.980 | -53.996 | +13.324 | A=0.01833, alpha=0.0001, B=0.4059, beta=0.1789 |
| M5 A exp(-c logn^d) | 3 | -49.902 | -55.664 | +4.402 | a=-0.6846, c=0.2743, delta=0.8083 |

### pl_dmin2 g=2.5

lambda2 by n_target: [0.11907, 0.10808, 0.09635, 0.08615, 0.07212, 0.05659, 0.04855, 0.04144] (rep counts [50, 50, 50, 50, 50, 50, 100, 100])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -55.616 | -57.857 | +0.000 | a=-1.381, alpha=0.158 |
| M2 A (log n)^-b | 2 | -39.886 | -42.127 | +15.730 | a=-0.235, beta=1.184 |
| M3 C+A n^-a | 3 | -50.016 | -55.778 | +5.600 | C=2.612e-11, A=0.2513, alpha=0.158 |
| M4 A n^-a+B n^-b | 4 | -40.683 | -53.698 | +14.933 | A=0.2513, alpha=0.158, B=4.32e-06, beta=8.814 |
| M5 A exp(-c logn^d) | 3 | -50.354 | -56.116 | +5.262 | a=-1.479, c=0.1201, delta=1.091 |

### pl_dmin2 g=3.0

lambda2 by n_target: [0.07062, 0.06105, 0.05083, 0.04487, 0.03617, 0.03119, 0.02628, 0.02115] (rep counts [50, 50, 50, 50, 50, 50, 100, 100])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -51.752 | -53.993 | +0.000 | a=-1.902, alpha=0.1708 |
| M2 A (log n)^-b | 2 | -48.501 | -50.742 | +3.251 | a=-0.6399, beta=1.291 |
| M3 C+A n^-a | 3 | -51.990 | -57.752 | -0.238 | C=0.008505, A=0.1717, alpha=0.2236 |
| M4 A n^-a+B n^-b | 4 | -46.734 | -59.750 | +5.018 | A=0.1247, alpha=0.1531, B=0.3817, beta=0.8085 |
| M5 A exp(-c logn^d) | 3 | -54.357 | -60.119 | -2.606 | a=-0.9153, c=0.732, delta=0.5665 |

### pl_dmin3 g=2.1

lambda2 by n_target: [0.32054, 0.31271, 0.30216, 0.29019, 0.28291, 0.26682, 0.25739, 0.24916] (rep counts [50, 50, 50, 50, 50, 50, 100, 100])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -74.057 | -76.298 | +0.000 | a=-0.9679, alpha=0.03741 |
| M2 A (log n)^-b | 2 | -68.132 | -70.373 | +5.925 | a=-0.6926, beta=0.2822 |
| M3 C+A n^-a | 3 | -71.492 | -77.254 | +2.564 | C=0.1514, A=0.2474, alpha=0.08135 |
| M4 A n^-a+B n^-b | 4 | -62.158 | -75.174 | +11.898 | A=0.1519, alpha=0.0001, B=0.2469, beta=0.08146 |
| M5 A exp(-c logn^d) | 3 | -71.061 | -76.822 | +2.996 | a=-0.8387, c=0.1039, delta=0.6857 |

### pl_dmin3 g=2.5

lambda2 by n_target: [0.26435, 0.25575, 0.24545, 0.23921, 0.22967, 0.22349, 0.21409, 0.20745] (rep counts [50, 50, 50, 50, 50, 50, 100, 100])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -71.526 | -73.767 | +0.000 | a=-1.184, alpha=0.03453 |
| M2 A (log n)^-b | 2 | -79.685 | -81.926 | -8.159 | a=-0.9263, beta=0.2621 |
| M3 C+A n^-a | 3 | -77.289 | -83.050 | -5.762 | C=0.1668, A=0.172, alpha=0.1245 |
| M4 A n^-a+B n^-b | 4 | -71.489 | -84.504 | +0.037 | A=0.2883, alpha=0.02869, B=0.219, beta=0.6328 |
| M5 A exp(-c logn^d) | 3 | -79.009 | -84.770 | -7.482 | a=-0.5947, c=0.4604, delta=0.3078 |

### pl_dmin3 g=3.0

lambda2 by n_target: [0.20868, 0.19754, 0.19043, 0.18924, 0.18153, 0.17451, 0.17241, 0.16812] (rep counts [50, 50, 50, 50, 50, 50, 100, 100])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -61.225 | -63.466 | +0.000 | a=-1.46, alpha=0.02946 |
| M2 A (log n)^-b | 2 | -70.236 | -72.477 | -9.011 | a=-1.234, beta=0.2265 |
| M3 C+A n^-a | 3 | -65.952 | -71.713 | -4.727 | C=0.1587, A=0.1388, alpha=0.2295 |
| M4 A n^-a+B n^-b | 4 | -59.270 | -72.285 | +1.955 | A=0.2202, alpha=0.02389, B=13.62, beta=1.541 |
| M5 A exp(-c logn^d) | 3 | -64.202 | -69.963 | -2.977 | a=2.843, c=4.098, delta=0.05 |

### ER c=2

lambda2 by n_target: [0.02907, 0.02308, 0.01603, 0.01388, 0.00968, 0.00771, 0.00551, 0.00455] (rep counts [20, 20, 20, 20, 20, 20, 20, 20])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -36.957 | -39.198 | +0.000 | a=-2.384, alpha=0.2697 |
| M2 A (log n)^-b | 2 | -46.544 | -48.785 | -9.587 | a=-0.37, beta=2.049 |
| M3 C+A n^-a | 3 | -42.709 | -48.471 | -5.752 | C=0.002577, A=0.1462, alpha=0.3739 |
| M4 A n^-a+B n^-b | 4 | -34.469 | -47.485 | +2.488 | A=0.03033, alpha=0.1723, B=0.2069, beta=0.563 |
| M5 A exp(-c logn^d) | 3 | -43.847 | -49.609 | -6.890 | a=3.741, c=4.96, delta=0.2506 |

### ER c=3

lambda2 by n_target: [0.09286, 0.06292, 0.05886, 0.05048, 0.0421, 0.03224, 0.02558, 0.01994] (rep counts [20, 20, 20, 20, 20, 20, 20, 20])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -37.874 | -40.116 | +0.000 | a=-1.546, alpha=0.2059 |
| M2 A (log n)^-b | 2 | -35.196 | -37.437 | +2.678 | a=-0.03239, beta=1.552 |
| M3 C+A n^-a | 3 | -32.409 | -38.170 | +5.466 | C=0.003168, A=0.2257, alpha=0.2244 |
| M4 A n^-a+B n^-b | 4 | -31.301 | -44.317 | +6.573 | A=0.1865, alpha=0.1917, B=9.048e+07, beta=4.885 |
| M5 A exp(-c logn^d) | 3 | -32.728 | -38.490 | +5.147 | a=-0.9959, c=0.4776, delta=0.7379 |

### ER c=4

lambda2 by n_target: [0.15662, 0.15057, 0.12634, 0.10636, 0.09063, 0.07454, 0.05946, 0.05277] (rep counts [20, 20, 20, 20, 20, 20, 20, 20])

| model | k | AICc | BIC | dAICc vs M1 | params |
|---|---|---|---|---|---|
| M1 A n^-a | 2 | -48.738 | -50.979 | +0.000 | a=-1.053, alpha=0.1675 |
| M2 A (log n)^-b | 2 | -40.333 | -42.574 | +8.405 | a=0.1708, beta=1.259 |
| M3 C+A n^-a | 3 | -43.986 | -49.748 | +4.751 | C=0.012, A=0.3706, alpha=0.1943 |
| M4 A n^-a+B n^-b | 4 | -34.653 | -47.668 | +14.085 | A=0.01202, alpha=0.0001, B=0.3706, beta=0.1943 |
| M5 A exp(-c logn^d) | 3 | -43.493 | -49.254 | +5.245 | a=-0.8595, c=0.2542, delta=0.8665 |


## 5. Where offset / two-power models win: saturation / crossover

Only M3 (`C + A n^-a`) and M4 (`A n^-a + B n^-b`) can signal a floor or a crossover. Readings:

**Clear saturation (M3 beats M1 by > ~4 AICc; nonzero floor C):**

- **pl_dmin3 g=2.5**: dAICc(M3-M1) = -5.76, floor C = 0.1668 (lambda2: 0.264 at n=100 -> 0.207 at n=1e5). lambda2 approaches a constant, not zero.
- **pl_dmin3 g=3.0**: dAICc(M3-M1) = -4.73, floor C = 0.1587 (lambda2: 0.209 -> 0.168). Clear plateau.
- **ER c=2**: dAICc(M3-M1) = -5.75, floor C = 0.0026 (small but statistically resolved; M2 log-power even better here).
- **pl_dmin3 g=2.1**: M1 still wins by +2.56, but M3 needs a very large floor C = 0.1514 (compared with lambda2 = 0.321..0.249); i.e. the range is dominated by an offset and alpha is poorly identified (M3 alpha = 0.081 vs M1 alpha = 0.037).

**Marginal / weak:** pl_dmin2 g=3.0 (dAICc M3-M1 = -0.24, M5 -2.61), pl_dmin2 g=2.1 and g=2.5 and ER c=3,c=4 show no floor (M3 dAICc > +3.99 with C near 0).

**Two-power M4:** it is never the AICc winner. It is only competitive for pl_dmin3 g=2.5 (dAICc M4-M1 = +0.04) and g=3.0 (+1.96); elsewhere it collapses to M3 (one amplitude -> 0) and is penalised as an over-fit. Hence the deviation from a pure power law is best described as an additive floor, not as a genuine two-exponent crossover.

## 6. Verdict

**Pure power law `A n^-a` passes the strict test (best AICc, all alternatives worse, stable alpha) for:**

- `pl_dmin2 g=2.1` alpha = 0.150 (dAICc vs M2/M3/M4/M5 = +9.9/+4.0/+13.3/+4.4)
- `pl_dmin2 g=2.5` alpha = 0.158 (+15.7/+5.6/+14.9/+5.3)
- `ER c=3` alpha = 0.206 and `ER c=4` alpha = 0.168 (controls)
- `pl_dmin3 g=2.1` alpha = 0.037 only marginally (M1 wins by +2.56 over M3, but M3 = C+A n^-a already carries a floor C = 0.151).

**Pure power law is REJECTED (an alternative model wins by dAICc, no floor) for:**

- `pl_dmin2 g=3.0`: M5 wins by 2.61 (M3 also ties, -0.24)
- `pl_dmin3 g=2.5`: M2 log-power wins by 8.16; M3 offset wins by 5.76
- `pl_dmin3 g=3.0`: M2 log-power wins by 9.01; M3 offset wins by 4.73; alpha falls monotonically 0.038 (n<=3000) -> 0.021 (n>=3000) -> 0.016 (n>=1e4)
- `ER c=2`: M2 wins by 9.59.

**alpha stability.** For the surviving sets alpha changes by <= 0.03 across lower cuts (n>=1000, >=1e4) and upper cuts (n<=1e4, <=3e3). For the rejected sets alpha drifts systematically with the window: e.g. pl_dmin3 g=3.0 goes 0.0321 (n<=3e4) -> 0.0163 (n>=1e4); pl_dmin2 g=2.5 0.160 (n<=3e4) -> 0.135 (n>=1e4), and pl_dmin3 g=2.5 collapses to a plateau. This window dependence is the direct signature that the fitted alpha is not a true asymptotic exponent.

**On the stated claim.** `alpha(d_min=2) ~ 0.15-0.17` is reproduced for gamma = 2.1 (0.150) and 2.5 (0.158), and gamma = 3.0 (0.171) -- but for gamma = 3.0 the pure power law is not the best description (stretched exponential/offset preferred), so 0.171 is an effective local slope. `alpha(d_min=3) ~ 0.03` is reproduced in magnitude (0.037 / 0.035 / 0.030 for gamma = 2.1 / 2.5 / 3.0), but the clean power law is rejected for gamma = 2.5 and 3.0 and holds for gamma = 2.1 only with a large constant offset. In short, the values are of the right size, but the 'pure power-law in n' interpretation is only valid in a minority of the tested sets; for d_min=3 at gamma >= 2.5 the data prefer a saturating (offset / log-power) form.
