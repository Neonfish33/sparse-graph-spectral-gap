# Toy model — plan for tomorrow

> **English translation.** The original and authoritative language is Russian:
> see [`TOY_MODEL.md`](TOY_MODEL.md).

Goal: **explain the exponent α** in `λ₂(giant) ~ n^(−α)` through the parameters of the
degree distribution — find the correct structural parameter and derive α.

---

## 0. Hard constraints from the experiments (do not rediscover!)

| family | tail | α_core |
|---|---|---|
| uniform {2..D} | compact | **0.06–0.08 (flat)** over D=3..20, d̄=2.5..11, q=0.40..0.01 |
| power-law {2..D} | heavy | 0.096 (D=3) → **0.164** (D=50) |
| ER G(n,c/n) | Poisson | 0.072 (c=4) → **0.138** (c=2) |
| ER giant (with leaves) | Poisson | 0.168 (c=4) → 0.270 (c=2) |

**Established:** α is **not** a function of (d̄, P(2), max degree, branching
ρ=E[D(D−1)]/E[D]). Uniform D=20 (d̄=11, ρ≈13) gives α=0.064, while ER c=2 (d̄=2.7)
gives α=0.138. ⇒ The **tail shape** matters (unboundedness/heaviness), not the moments.
The "hub" hypothesis (a jump at max degree ≥5) is **refuted**.

## 1. Hypotheses

- **H1 (spectral density).** `λ₂` is the edge of the finite-size tail of the density of
  small eigenvalues. If the integrated density near 0 behaves as `F(λ) ~ λ^β`, then
  `λ₂ ~ n^(−1/β)`, i.e. **α = 1/β**. Here β is the "spectral dimension". Numerically
  checkable: measure β and compare α with 1/β.
- **H2 (degree-2 chains).** α is determined by the length distribution of the maximal
  chains of degree-2 vertices (extremal length ℓ_max, tail P(ℓ>L)).
- **H3 (backbone).** After suppressing degree-2 chains one obtains a multigraph with
  resistances (≈ chain lengths); its spectral gap is related to α.

## 2. Numerical experiments (run by `toy_model.py`)

1. **H1:** for each model at fixed n=10⁴ measure the K=150 smallest eigenvalues of the
   giant, fit `log F(λ)` vs `log λ` → β, then α_pred=1/β. Compare with the known α
   (unif≈0.07, pl≈0.16, ER≈0.09–0.14).
2. **H2:** distribution of degree-2 chain lengths: mean, max, P(ℓ>L). Heavy tail for PL,
   geometric for uniform?
3. **H3:** suppressed multigraph with weights = chain lengths; λ₂ and its dependence.
4. **Crossover:** exponential tail `P(k) ∝ e^{−θk}` on [2,∞), vary θ: as θ→∞ compact
   (α≈0.07), as θ→0 heavy (α↑). Find the transition point.

## 3. Analysis (goal)

Derive β (or the ℓ-statistics) from the generating function of the degree distribution /
the local weak limit (unimodular Galton–Watson tree) and obtain `α = F({p_k})`.
Landmarks in the literature: Hermon–Li–Yao–Zhang (min degree ≥3 ⇒ Θ(1)),
Chung–Lu–Vu (d̄_min ≫ ln²n), "spectral dimension" (Avrachenkov et al.).

## 4. Success criterion

A single function φ(distribution) is found such that `α ≈ φ` simultaneously for uniform
(0.07), power-law (0.10–0.16) and ER (0.07–0.14). Then the toy model explains both the
generic baseline (B) and the tail amplification (A).

## 5. Files

- `toy_model.py` — skeleton (H1+H2, stubs for H3/H4), CLI.
- `run_toy_model.bat` — quick launch.
- Results → `results_toy_model/`.

Run (quick check):
```
python toy_model.py --quick
python toy_model.py --n 10000 --K 150 --out results_toy_model
```
