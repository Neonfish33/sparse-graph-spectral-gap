# Full project review — results, conclusions, what is closed

> **English translation.** The original and authoritative language is Russian:
> see [`REVIEW.md`](REVIEW.md).

A single document for re-checking. Everything we did, found, and why we stopped.
References to files/numbers are given so it can be verified independently.

Final date: 2026-09-23. Repository: https://github.com/Neonfish33/sparse-graph-spectral-gap

---

## 0. TL;DR

- We studied the **spectral gap λ₂ of the normalized Laplacian** of sparse random
  graphs with heavy-tailed degrees (theory + numerical analysis).
- We obtained clean scalings empirically (`λ₂ ~ n^{−α}`, etc.).
- **The literature check showed the core is already published (2003–2008).** Our
  result is a **reproduction** (minimum degree sets the spectral edge; the tail is
  secondary; the finite-size λ₂(N) is also known).
- **No substantive novelty.** Project closed. The repo is kept as a tool/benchmark.
- We considered ~6 alternative directions — **all also occupied** (see §8).

---

## 1. Original problem

Obtain the asymptotics of the second eigenvalue `λ₂` of the normalized Laplacian
`L = I − D^{−1/2} A D^{−1/2}` for sparse random graphs with `P(D=k) ∼ C k^{−γ}` in
the regime `d̄_min = O(1)`, not covered by Chung–Lu–Vu theory.

Ensembles: power-law (`d_min = 1, 2, 3`) and `planted` (modularity).

---

## 2. What was built

### Code (`src/`, root)
- `src/generators.py` — degree sequences, configuration model, 4 ensembles.
- `src/diagnostics.py` — `λ₂` (full/giant/2-core), sweep cut, `φ_periph`, `φ_mod`,
  `φ*`, `Δ`, `ρ*`, IPR₂, φ-profile.
- `run_experiments.py` — orchestration (flags `--jobs/--mu/--dmax/--resume/--save-graphs/...`).
- `analyze_results.py` — α(n), α(φ), phase diagram, regimes.
- `analyze_gap_fss.py` — finite-size scaling of the gap (`λ_c + C n^{−θ}`) and distributions.
- `validate_power_law.py` — AICc: power law vs `1/log n` vs exp-log; bootstrap α.
- `run_er_control.py` — ER control; `run_toy_controls.py` — tail vs degree-2 fraction.
- `toy_model.py`, `merge_diagnostics.py`, `export_graph.py`.
- `tests/` — 19 tests (spectra, sweep cut vs brute force, Cheeger).

### Data
| Directory | Content | In git? |
|---|---|---|
| `results_v2/` | 13500 graphs, n=100…10⁴, γ=2.1/2.5/3.0, d_min=1/2/3, planted | ❌ (graphs 231 MB) |
| `results_large/` | n=3×10⁴ and 10⁵, d_min=2/3, reps=100 | plots ✅ |
| `results_merged/` | v2+large, **14700 rows**, n=100…10⁵ | `diagnostics.csv` ✅ |
| `results_er/`, `results_toy/` | ER control, toy controls | ✅ |
| `results_agents/` | independent-agent reports | ✅ |

---

## 3. Empirical results (numbers)

### 3.1. α in `λ₂(giant) ~ n^{−α}` (over n = 100…10⁵, 8 points)
| ensemble | γ=2.1 | γ=2.5 | γ=3.0 |
|---|---|---|---|
| `pl_dmin1` | 0.268 | 0.401 | 1.100 |
| `pl_dmin2` | **0.150** (R²=1.00) | **0.158** (R²=1.00) | **0.171** (R²=1.00) |
| `pl_dmin3` | 0.037 | 0.035 | 0.029 |

- For d_min=2 the decay is **power-law** over 3 decades, no saturation (outcome (b)).
- Local α on the last decade (3×10⁴→10⁵): 0.145 / 0.133 / 0.183 — no systematic slowing.

### 3.2. Gap `λ₂(d3)/λ₂(d2)`
| γ | n=10³ | 10⁴ | 3×10⁴ | 10⁵ |
|---|---|---|---|---|
| 2.1 | 2.18 | 2.80 | 3.09 | **3.56** |
| 2.5 | 2.78 | 3.95 | 4.41 | **5.01** |
| 3.0 | 4.22 | 5.59 | 6.56 | **7.95** |

### 3.3. Cheeger
`α_φ` (λ₂ ~ (φ*)^α) ∈ [1, 2] for d2: 1.41 / 1.45 / 1.52; for d3 there is no power law.

### 3.4. Localization (key)
`IPR₂(giant) ≈ const` (0.07–0.19) as n grows ⇒ **λ₂ is localized on ~6–11 vertices**,
the support is a chain of degree-2 vertices (verified by dense recomputation by the
devil's-advocate agent).

### 3.5. Degree-2 fraction
`P(2) = N₂/n` for d2: 0.426 / 0.520 / 0.619 (γ=2.1/2.5/3.0); for d3: ≈0.003.
`q = 2·P(2)/d̄`.

### 3.6. Combinatorial vs normalized Laplacian (`compare_laplacians.py`)
On the **same** saved graphs (pl_dmin2, n=100…10⁴, 10 reps):

| γ | observable | α (power) | R² | β (log) | R² |
|---|---|---|---|---|---|
| 2.1 | norm | 0.151 | 0.963 | 1.01 | 0.964 |
| 2.1 | **comb** | 0.186 | 0.969 | 1.24 | 0.980 |
| 2.5 | norm | 0.151 | 0.890 | 1.00 | 0.885 |
| 2.5 | comb | 0.183 | 0.895 | 1.22 | 0.902 |
| 3.0 | norm | 0.192 | 0.976 | 1.28 | 0.984 |
| 3.0 | comb | 0.212 | 0.977 | 1.42 | 0.982 |

**Conclusions:**
- The combinatorial Laplacian decays **the same way** (power-law-like) as the
  normalized one ⇒ "normalized vs combinatorial" is **not** a gap; it is one phenomenon.
- Power law and `(ln N)^{−2}` are **indistinguishable** for n≤10⁴ (R² differ by ~0.001–0.01).
- ⇒ the 2008 asymptotics `(ln N)^{−2}` is **not visible** on our range (or it sets in at n≫10⁵).

### 3.7. Asymptotic check up to n=10⁶ (`compare_asymptotics.py`)
Added n=3×10⁵ and 10⁶ (pl_dmin2/3, γ=2.5/3.0, reps=5) in `results_asymptotic/`,
merged with the main dataset (n=100…10⁵). Full range n=100…10⁶.

| ensemble/γ | α(≤10⁵) | α(≤10⁶) | R²(power) | R²(log) | better |
|---|---|---|---|---|---|
| d2 γ=2.5 | 0.158 | 0.147 | 0.989 | 0.978 | power |
| d2 γ=3.0 | 0.171 | 0.157 | 0.990 | 0.995 | **log** |
| d3 γ=2.5 | 0.035 | 0.032 | 0.987 | 0.994 | **log** |
| d3 γ=3.0 | 0.030 | 0.026 | 0.958 | 0.991 | **log** |

Local log-log slopes at large n:
- d2 γ=2.5: 10⁵→3×10⁵ = **0.015**, 3×10⁵→10⁶ = **0.242** (noisy: reps=5);
- d2 γ=3.0: 0.108 → 0.086 (decreasing);
- d3 γ=2.5/3.0: ~0.02–0.04 → ~0.007–0.021 (plateau).

**Conclusion:** extending to n=10⁶, α **decreases** for all sets, and the logarithmic
form `A(ln n)^{−β}` (β≈1.2–1.33) becomes preferred in **3 of 4** cases. This is
**consistent** with the 2008 asymptotics `(ln N)^{−2}` — the decay **slows** toward a
logarithm. Caveat: β<2 (not yet the exact `(ln N)^{−2}`), and at n=10⁶ there are only
5 realizations → noisy; the crossover is very slow.

**Bottom line:** the core of 2008 (min degree → edge) is **confirmed**; the finite α is
effective, the asymptotic is logarithmic. **No novelty**, but the picture is closed.

### 3.8. Why the "normalized vs combinatorial" difference does NOT yield a result

**Observation.** The combinatorial λ₂ decays slightly faster than the normalized one
(α_comb ≈ 0.19–0.21 vs α_norm ≈ 0.15–0.19, n=100…10⁴). The ratio `λ_comb/λ_norm`
**decreases slowly**: γ=2.1 2.48→2.07; γ=2.5 2.39→2.07; γ=3.0 2.23→2.04.

**Interpretation.** This is a **finite-size normalization effect**, not a qualitative
difference:
- both decay (both are a min-degree-≤2 phenomenon); the `q_m>2` criterion is common;
- the ratio decreases slowly (not `~ n^α`), i.e. a rescaling, not a new law;
- asymptotics up to n=10⁶: α **decreases**; the log form `A(ln n)^{−β}` (β≈1.2–1.33)
  is preferred in **3 of 4** cases (exception: d2 γ=2.5 — noisy, reps=5).

**Consequence.** The difference lies in the **effective finite-size exponent** due to
normalization; the qualitative conclusion of Samukhin 2008 (the role of min degree)
applies to the normalized Laplacian as well. A strict derivation for the normalized case
remains open (§9.5), but the expected result is a **confirmation, not novelty**.

**Bottom line:** no novelty — **confirmed experimentally** (n up to 10⁶), not assumed.

---

## 4. Verification by independent agents

### 4.1. Devil's-advocate agent (`results_agents/devils_advocate/report.md`)
Refuted (low severity): n jitter, eigensolver switch, v2+large merge, disconnectedness.
Serious:
- **(#5/#6) PL vs ER comparison — gluing different objects.** We compared
  `λ₂(giant)` of PL with `λ₂(core)` of ER. Like-for-like:
  - **giant vs giant:** PL 0.15–0.17 is **smaller** than ER 0.17–0.27 → PL *slower*;
  - **core vs core:** 1.5–2.1×.
  ⇒ "PL is 2× faster than ER" is an **artifact**.
- **(#10) λ₂ is a localized defect**, not a global gap (IPR·n grows linearly).
- **(#11) The causal "tail" claim is confounded** by the degree-2 fraction (the d2→d3
  transition coincides with P(2) vanishing).
- **(#7/#8) R²≈1.0 is an averaging effect**; AICc for d2 γ3.0 and d3 selects exp-log.

### 4.2. FSS agent (`results_agents/fss/report.md`)
- A pure power law **passes** for: d2 γ=2.1 (0.150), γ=2.5 (0.158), ER c=3 (0.206),
  c=4 (0.168), d3 γ=2.1 (marginally).
- **Rejected** for: d2 γ=3.0 (`exp(−c(log n)^β)`), d3 γ=2.5 and 3.0 (`(log n)^{−β}`
  / plateau), ER c=2.
- **Saturation** (nonzero floor C≈0.16) for d3 γ=2.5, 3.0.
- α drifts with the n-range for the "rejected" sets ⇒ effective, not asymptotic.

---

## 5. Controls

### 5.1. ER control (`run_er_control.py`)
| model | core d̄ | α(core) | α(giant) |
|---|---|---|---|
| ER c=2 | 2.68 | 0.138 | 0.270 |
| ER c=3 | 3.44 | 0.089 | 0.206 |
| ER c=4 | 4.26 | 0.072 | 0.168 |

Conclusion: **the qualitative decay is generic** (ER decays too) ⇒ not power-law-specific.

### 5.2. Toy controls (`run_toy_controls.py`)
| family | tail | α_core |
|---|---|---|
| uniform {2..D}, D=3..20 | compact | **0.06–0.08 (flat)** |
| power-law {2..D}, D=3..50 | heavy | 0.096 → 0.164 |
| {2,3}, p2=0.4..0.8 | — | 0.045 → 0.092 |
| {2,3,4} | — | 0.094 / 0.098 |

Conclusion: the "hub" hypothesis (a jump at max degree ≥5) is **refuted**; but α is
**not** a function of (d̄, P(2), max degree, branching ρ). The tail shape matters — but
this was not confirmed unambiguously (see §4.1 #11).

---

## 6. Literature: why there is no novelty (details in `LITERATURE.md`)

🔴 **Key paper — Samukhin, Dorogovtsev, Mendes, PRE 77, 036115 (2008)**
(arXiv:0706.1176), *"Are scale-free architectures really important?"*:
- the lower spectral edge `λ_c` is set by the **minimum degree** `q_m`;
- **`λ_c > 0` for `q_m > 2`, `λ_c = 0` for `q_m ≤ 2`** ⇒ our d2/d3 dichotomy;
- "high-degree part … is not that essential" ⇒ our "the tail is secondary" conclusion;
- **finite size:** `λ₂(N) ∼ (ln N)^{−2}` for `q_m ≤ 2` and
  `λ₂(N) − λ_c ∼ (ln ln N)^{−2}` for `q_m > 2` ⇒ our "power law" is asymptotically logarithmic;
- the normalized Laplacian was anticipated ("main conclusions are also valid if the
  escape rate from a vertex is fixed").

Additionally:
- **DGMS, PRE 68, 046109 (2003)** — spectra of locally tree-like networks, `ρ(λ) ∼ |λ|^{1−2γ}`, role of low degrees.
- **ISCCSP 2012** (10.1109/isccsp.2012.6217860) — zero gap of the normalized Laplacian of the ER giant.
- **Hata & Nakao, Sci Rep 7, 1121 (2017)** — localization of Laplacian eigenvectors (dense regime).

**Conclusion:** the core, the finite size, and the qualitative localization are already
published (2003–2008).

---

## 7. What was NOT confirmed (refuted interpretations)

1. "Power-law is ~2× faster than ER" — artifact of an unmatched comparison (giant vs core).
2. The causal "heavy tail accelerates the decay" — confounded with the degree-2 fraction.
3. "A clean universal power law" — true only for d2 γ=2.1, 2.5.
4. "λ₂ is a global spectral gap" — in fact a localized defect.
5. My prediction "saturation (a)" — not confirmed (outcome (b) on the accessible range).
6. The "hub" hypothesis (max degree) — refuted.
7. The first explanation "bounded degrees ⇒ expander" — fails for d_min≤2.

---

## 8. Considered pivots — all occupied

| Direction | Status |
|---|---|
| Hypergraph SBM / thresholds | occupied (Pal–Zhu 2021, Stephan–Zhu 2022/24, Fernandez–Stephan–Zhu 2026) |
| Delocalization of Laplacian eigenvectors | occupied (Huang–Landon 2015, Metz–Neri PRL 2021, García-Mata et al.) |
| Temporal × spectral | small field; anchors 2013 (Masuda) + 2024–2026; **one gap — the non-Hermitian branch** |
| Distributed learning / mesh + TinyML | occupied (RockNet 2025, AIRPLAN 2026, DFL topology) |
| Biology (connectomes/Hi-C/GRN/PPI) | occupied; competition from systems/ML groups |
| Gacha/loot boxes | occupied by psychology; a niche may exist in the **math of pity mechanics** |

---

## 9. What potentially remains open (honest, no guarantees)

1. **Non-Hermitian spectral theory of time-respecting walks** — the query
   "non-Hermitian + temporal network" gives ~0 papers. Requires deep immersion.
2. **Spectral dimension** of time-respecting walks (`p(t) ~ t^{−d_s/2}`).
3. **Localization of temporal modes** (eigenvectors of temporal operators).
4. **Math of pity systems in gacha** (Markov chains/renewal/hitting times) — possibly unoccupied.
5. **Strict derivation for the normalized Laplacian.** Samukhin 2008 worked with the
   **combinatorial** Laplacian; the transfer to the normalized one is **asserted**
   ("main conclusions are also valid if the escape rate is fixed") but **not derived**
   — their scalar resolvent technique does not carry over to the non-symmetric operator
   `D^{−1}A` (non-normal, size-biased walk). Our numerical data (§3.6) **support** the
   transfer. A strict theorem is **open**. Assessment: technically hard; the expected
   result is a confirmation (not a discovery); without a specialist co-author **not recommended**.

All are narrow and risky; none guarantees a "big journal" without collaboration.
The only "own" unclosed door is item 5; but it leads to confirmation, not novelty.

---

## 10. How to re-check yourself

### Key numbers
- `results_merged/analysis/scaling_by_n.csv` — α(n) for all ensembles.
- `results_merged/analysis/core_degree_comparison.csv` — `ratio_d3_over_d2` vs n.
- `results_merged/analysis/scaling_by_phi.csv` — α_φ (Cheeger).
- `results_agents/devils_advocate/report.md`, `results_agents/fss/report.md` — critique.

### Commands
```bash
python analyze_results.py --csv results_merged/diagnostics.csv --out results_merged/analysis
python validate_power_law.py          # AICc + bootstrap α
python analyze_gap_fss.py             # finite-size λ_c, θ, distributions
python -m pytest tests/ -q            # tests
```

### What to re-check first (if looking for something missed)
1. **Does 2008 really cover the normalized Laplacian?** (they used the combinatorial
   one; they only asserted the transfer). Open PRE 77, 036115 and check the sections
   on "fixed escape rate".
2. **Is α=0.15 finite-size?** (their `(ln N)^{−2}`): needs n≫10⁵.
3. **The non-Hermitian temporal branch** — is it really empty (search wider:
   pseudospectra, non-normal, exceptional points).

---

## 11. Bottom line

- **Technically:** a full, reproducible spectral pipeline + data + tests was built.
- **Scientifically:** no substantive novelty — a known result (2003–2008) was reproduced.
- **The right call:** close the topic (done) and avoid overcrowded niches.
- **What remains:** either a narrow non-Hermitian/spectral temporal topic (risky, needs
  collaboration), or an artifact/engineering project, or another domain (e.g., the math of gacha).

---

*If, while re-checking, you find something we missed — send it: we will look at it point by point.*
