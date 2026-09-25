# Project changelog

> **English translation.** The original and authoritative language is Russian:
> see [`CHANGELOG.md`](CHANGELOG.md).

Project: the spectral gap of the normalized Laplacian in sparse random graphs with
heavy-tailed degrees.

---

## 2026-09-23 — Literature check: the core reproduces PRE 2008; pivot

Reading the key works showed that the central qualitative result is **already known**:

- **Samukhin, Dorogovtsev, Mendes**, *Are scale-free architectures really important?*,
  **PRE 77, 036115 (2008)**: the lower spectral edge is set by the **minimum degree**;
  `λ_c>0` for `q_m>2`, `λ_c=0` for `q_m≤2`; the degree tail is irrelevant. ⇒ **our
  d_min=2/3 dichotomy and the "tail secondary" conclusion reproduce a known result.**
- **DGMS**, *Spectra of complex networks*, **PRE 68, 046109 (2003)**: spectra of locally
  tree-like networks, `ρ(λ)∼|λ|^{1−2γ}`, role of low degrees.
- **Hata & Nakao**, **Sci Rep 7, 1121 (2017)**: localization of Laplacian eigenvectors
  (degree–eigenvalue correspondence) but in the **dense** regime `⟨k⟩≫1`.

**Bottom line:** the "new law" is retracted. The project pivots to a question **beyond
2003–2008**: finite-size scaling and fluctuations of the gap at the `q_m=2→3` transition
(direction A). The literature check for A found **no direct match**
(finite-size scaling of the spectral gap at the min-degree transition, preliminarily,
is unpublished).

The previous results (v2, large, ER, toy controls) retain value as a reproduction /
numerical detail of the known result **and as data for A**.

---

## 2026-09-23 — Stage 3a: toy controls (tail vs degree-2 fraction)

`run_toy_controls.py`, `results_toy/` (gitignored). Families: `pl{γ=2.5,D}`,
`uniform {2..D}`, `{2,3}` and `{2,3,4}`; n = 10³,10⁴,10⁵, reps = 20.

| family | tail | α_core |
|---|---|---|
| uniform {2..D}, D=3..20 | compact | **0.06–0.08 (flat)** over D, d̄, q |
| power-law {2..D}, D=3..50 | heavy | 0.096 → **0.164** |
| {2,3}, p2=0.4..0.8 | — | 0.045 → 0.092 |
| {2,3,4} | — | 0.094 / 0.098 |

**Conclusions:**
- The "hub" hypothesis (a jump at max degree ≥5) is **refuted**: uniform α≈0.07 for D up
  to 20, including d̄=11.
- α is **not** a function of (d̄, P(2), q, max degree, branching ρ). Uniform D=20
  (d̄=11, ρ≈13) gives α=0.064, while ER c=2 (d̄=2.7) gives α=0.138.
- The **tail shape** matters: a heavy/unbounded tail amplifies α ~2× versus compact.
  The degree-2 fraction (hypothesis B) alone does not explain it (q varies by 40×, α is flat).
- The correct structural parameter was not found → a task for the toy model.

### Toy model prepared for tomorrow

- `TOY_MODEL.md` — design: H1 (α = 1/β, `F(λ)~λ^β`), H2 (degree-2 chain lengths),
  H3 (backbone), exp-tail crossover; success criterion.
- `toy_model.py` — skeleton (low-eig density, β, chain lengths, models unif/pl/two23/er/expo), CLI `--quick`.
- `run_toy_model.bat` — launch.
- Quick run: H1 in its naive form (fitting β over all K) **did not match** α —
  needs refinement (fit only near 0, larger n).

---

## 2026-09-23 — Stages 0/1/1.5: law validation, literature search, ER control

### Stage 0 — `validate_power_law.py` (new)

Check on `results_merged/diagnostics.csv`: AICc (power law vs `A/log n` vs
`A·exp(−c·(log n)^β)`), bootstrap CI for α (resampling realizations within n), local α
per decade. **Important:** filter by ensemble (`pl_dmin*`), otherwise the planted
ensembles spoil α (the first version gave 0.117 instead of 0.150).

Results (λ₂_giant):

| d_min | γ | α (log-fit) | R² | bootstrap 95% CI | AICc winner |
|---|---|---|---|---|---|
| 2 | 2.1 | 0.1496 | 0.993 | [0.1433, 0.1561] | power_law |
| 2 | 2.5 | 0.1580 | 0.996 | [0.1493, 0.1664] | power_law |
| 2 | 3.0 | 0.1708 | 0.994 | [0.1612, 0.1796] | **exp-log** (ΔAICc≈10) |
| 3 | 2.1 | 0.0374 | 0.994 | [0.0348, 0.0399] | power_law |
| 3 | 2.5 | 0.0345 | 0.990 | [0.0316, 0.0373] | exp-log |
| 3 | 3.0 | 0.0295 | 0.953 | [0.0262, 0.0326] | exp-log |

Verdict: consistent with a power law (α matched `analyze_results`), α grows with γ
(the 2.1 and 3.0 CIs do not overlap), no saturation. **But** for γ=3.0 (d2) and for d3,
`exp(−c·(log n)^β)` (β≈0.5, slower than any power law) fits better — the "pure power
law" is not unambiguous.

### Stage 1 — literature search (details in `LITERATURE.md`)

arXiv API + Crossref + OpenAlex. No exact scaling found. **Critically:** a precedent was
found — the normalized Laplacian giant of **sparse ER** has a zero gap as n→∞
(ISCCSP 2012, 10.1109/isccsp.2012.6217860). So the *qualitative* decay is not new.

### Stage 1.5 — ER control (`run_er_control.py`, `results_er/`)

ER `G(n, c/n)`, c ∈ {2,3,4}, n = 100…10⁵, reps = 20 → 480 graphs (~10 min, `--jobs 6`).
`results_er/` is on disk (gitignored).

| model | core d̄ | α(λ₂_core) | α(λ₂_giant) |
|---|---|---|---|
| ER c=2 | 2.68 | 0.138 | 0.270 |
| ER c=3 | 3.44 | 0.089 | 0.206 |
| ER c=4 | 4.26 | 0.072 | 0.168 |
| pl_dmin2 γ=3.0 | 3.08 | 0.171 | — |
| pl_dmin2 γ=2.5 | 3.87 | 0.158 | — |
| pl_dmin2 γ=2.1 | 5.06 | 0.150 | — |
| pl_dmin3 | 4.8–7.3 | 0.029–0.037 | — |

**Verdict:** (1) the qualitative decay is **generic** (ER also gives n^(−α)), not
power-law; (2) at equal core d̄ the power-law α is **~2× larger** than ER (≈0.16 vs
≈0.09 at d̄≈3.5) → the tail/degree structure adds; (3) **confound**: pl_dmin2 has a large
degree-2 fraction, the ER core does not; the effect may come from the degree-2 fraction
rather than the tail. The d2→d3 transition is sharp (0.15→0.03).

### Next

Bounded-support controls: **M_b** degrees ∈ {2,3} (vary the fraction p₂, no tail) and
**M_a** degrees ∈ {2..D} (vary D = d_max test). Separate "tail" from "degree-2 fraction".
This is the numerical toy model.

---

## 2026-09-23 — v2 full run, diagnostic system, fixes

### Data (not in git, on disk)

| Directory | Content |
|---|---|
| `results_v2/` | **13500 graphs**: n = 100…10000, γ = 2.1/2.5/3.0, d_min = 1/2/3, planted d_min = 1/2/3, μ = 0.05/0.1/0.2/0.3, d_max = 50, reps = 50 |
| `results_v2/graphs/` | all graphs, `.npz` (13500 files, **230.7 MB**), nodes 0..n−1, labels for planted |
| `results_v2/analysis/` | `core_degree_comparison.csv`, `scaling_by_n.csv`, `scaling_by_phi.csv`, `phase.csv`, `summary.txt`, 3 figures |
| `results_fixed/` | previous run (6300 graphs, planted only d_min=1) |
| `results_large/` | **large-n test**: n = 3×10⁴ and 10⁵, γ = 2.1/2.5/3.0, d_min = 2/3, no planted, reps = 100 → 1200 graphs; `phi_curves.csv` present |
| `results_merged/` | `results_v2` + `results_large` (**14700 rows**, n = 100…10⁵), `n_target` restored for v2; `analysis/` final tables |

### Code

- `src/diagnostics.py` — added `ipr`, `ipr_pi` (π-weighted IPR), `sweep_curve`,
  `graph_sweep_curve`; IPR₂/PR₂ in the metrics; **fixed the `φ* = 0` bug** (global-minimum
  update at `S = V`); correct handling of disconnected graphs (λ₂ = 0, sweep → nan);
  nan logic for `log Δ`.
- `src/generators.py` — `planted_d_mins`, `graph_id`, `save_graph`, `load_graph`.
- `run_experiments.py` — flags `--jobs`, `--mu` (several), `--dmax`, `--save-curves`,
  `--save-graphs`, `--resume`, `--no-planted`, `--planted-dmin`; column `n_target`;
  fixed the **duplicate `d_min` column** (observed → `d_min_obs`).
- `analyze_results.py` — α(n) and α(φ) with R² **within (ensemble, γ)** (no pooling);
  `core_degree_comparison`; phase diagram; `n_target`.
- `export_graph.py` — export `.npz` → edgelist / CSV / GraphML / GML / JSON.
- `run_large.bat` — launch the large-n test (3×10⁴, 10⁵).
- `tests/test_spectral.py` (6 tests), `tests/test_known_spectra.py` (13) — all pass.

### Key empirical results (finite n)

1. **Three regimes** by α in `λ₂(giant) ~ n^(−α)`:
   d_min=1 (periphery) α = 0.27…0.97 (grows with γ); d_min=2 α ≈ 0.14…0.19;
   d_min=3 α ≈ 0.03…0.04 (λ₂ ≈ const).
2. **Matched d_min=2 vs d_min=3** (both without leaves), the gap `λ₂(d3)/λ₂(d2)`
   grows with n: γ=2.1 1.68→3.56; γ=2.5 2.22→5.01; γ=3.0 2.96→7.95 (n=100→10⁵,
   see the "Large run" section below).
3. **Cheeger:** α_φ ∈ [1,2] for periphery, d_min=2 and planted; for d_min=3 there is no
   power law (λ₂ ≈ const).
4. **Planted:** the modular signal is visible in the core (λ₂(core) grows with μ); for
   d_min ≥ 2 planted graphs are connected.

### Fixed bugs

- `φ* = 0` and `ρ* = 1` due to the trivial `S = V` (spoiled all previous φ*/ρ*).
- Duplicate `d_min` column in the CSV.
- α was computed from the actual n (jitter ±units) → from the target `n_target`.
- The previously reported "α_φ = 1.90 for core" is a **pooling artifact** (within groups 1.19–1.59).
- The false "sweep cut fix" from an external review was **rejected** (my code matched brute force).

### Large run n = 3×10⁴, 10⁵ — RESULT (outcome (b))

The run finished 2026-09-23 (~2.5 h, `--jobs 6`), 1200/1200 graphs.

**Verdict: the power-law decay of λ₂(d_min=2) continues, no saturation. Outcome (b),
my prediction (a) was not confirmed.**

α in λ₂(giant) ~ n^(−α) over n = 100…10⁵ (8 points, `results_merged/analysis/scaling_by_n.csv`):

| ensemble | γ=2.1 | γ=2.5 | γ=3.0 |
|---|---|---|---|
| pl_dmin2 | 0.150 (R²=1.00) | 0.158 (R²=1.00) | 0.171 (R²=1.00) |
| pl_dmin3 | 0.037 (R²=0.99) | 0.035 (R²=0.99) | 0.029 (R²=0.96) |

R² = 1.00 over three decades ⇒ the 10⁵ points lie perfectly on a straight line, no curvature.
The local α(λ₂ d2) on the 3×10⁴→10⁵ decade: 0.145 / 0.133 / 0.183 vs the global
0.150 / 0.158 / 0.171 — no systematic slowing.

`ratio = λ₂(d3)/λ₂(d2)` (n = 100…10⁵):

| γ | 10³ | 10⁴ | 3×10⁴ | 10⁵ |
|---|---|---|---|---|
| 2.1 | 2.18 | 2.80 | 3.09 | **3.56** |
| 2.5 | 2.78 | 3.95 | 4.41 | **5.01** |
| 3.0 | 4.22 | 5.59 | 6.56 | **7.95** |

The prediction was: γ=2.1 → 3.0–3.6 (actual 3.56 ✓); γ=2.5 → 4.2–5.0 (actual 5.01,
slightly above); γ=3.0 → 5.8–7.0 (actual 7.95, above; a pure law would give 7.88).

Cheeger preserved: α_φ(pl_dmin2) = 1.41 / 1.45 / 1.52 (all in [1,2]).

**Meaning:** the argument "bounded degrees (d_max=50) ⇒ expander ⇒ λ₂ = Θ(1)" does
**not** work for heavy-tail + min degree 2. The decay n^(−0.15…−0.17) is a stable result,
not a finite-size effect. This is the central question for the toy model.

### Next steps

1. Toy model: explain the **power-law** decay at bounded degrees (not saturation) —
   non-backtracking / Alon–Boppana does not give `n^(−α)` here; find the mechanism
   (probably the growth of degree-2 chain lengths or the degree tail).
2. Restore `n_target` for v2 more reliably (currently rounding to the nearest target);
   ideally regenerate v2 with the current code.
3. Results section in `paper_en.tex` + `paper.tex` with outcome (b); honest limitations.
   *(outdated: the papers were removed, the project is closed — see `LITERATURE.md`)*
4. (Optional) blind validation, real networks.

---

## Early history (brief)

- 2026-09-21: initial code (generators/diagnostics/run_experiments), the `README.md`
  conspectus, `paper.tex`/`paper.pdf`, the English version `paper_en.tex`/`paper_en.pdf`.
  **Later (2026-09-23) the papers were removed** — their main claim turned out to be a
  reproduction of a known result (PRE 2008), see `LITERATURE.md`.
- The external archive `files.zip` with "fixes" was reviewed: 2 of 3 critical "bugs"
  were false, but a real `φ* = 0` bug was found.
