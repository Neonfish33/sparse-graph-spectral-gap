# Literature check for novelty

> **English translation.** The original and authoritative language is Russian:
> see [`LITERATURE.md`](LITERATURE.md).

Question: has anyone published the observed scaling
`λ₂(normalized Laplacian, giant) ~ n^(−α)` with α ≈ 0.15…0.17
for a power-law configuration model with **min degree 2** (d_max = 50)?

Method: arXiv API (full-text queries) + Crossref. Semantic Scholar returned 429.
Date: 2026-09-23.

---

## 🔴 UPDATE (after reading the key works): the core is ALREADY PUBLISHED

- **Samukhin, Dorogovtsev, Mendes**, *Laplacian spectra of complex networks and
  random walks on them: **Are scale-free architectures really important?***,
  **PRE 77, 036115 (2008)**, arXiv:0706.1176. Establishes: the lower spectral edge
  `λ_c` is set by the **minimum degree q_m**; "high-degree part … is not that
  essential"; **λ_c > 0 for q_m > 2, λ_c = 0 for q_m ≤ 2**. ⇒ our d_min=2/3 dichotomy
  and the "tail is secondary" conclusion are a **reproduction of a known result**.
- **Dorogovtsev, Goltsev, Mendes, Samukhin**, *Spectra of complex networks*,
  **PRE 68, 046109 (2003)**: exact equations for the spectra of locally tree-like
  networks; `ρ(λ) ∼ |λ|^{1−2γ}`; a separate analysis of the role of low degrees.
- **Hata & Nakao**, *Localization of Laplacian eigenvectors on random networks*,
  **Sci Rep 7, 1121 (2017)**: localization on nodes with similar degrees
  (degree–eigenvalue correspondence), but at `⟨k⟩≫1` (dense regime).

- **The same 2008 paper also gives the finite size λ₂(N)** (verbatim, Sec. I/Table 1):
  `λ₂(N) ∼ (ln N)^{−2}` for `q_m ≤ 2` and `λ₂(N) − λ_c ∼ (ln ln N)^{−2}` for
  `q_m > 2`. ⇒ our "power law `n^{−α}`" is an **effective** exponent for n ≤ 10⁵; the
  true asymptotic is **logarithmic**. The "finite-size scaling" direction (A) is also
  closed. They anticipated the normalized case ("main conclusions are also valid if the
  escape rate from a vertex is fixed").

**FINAL VERDICT (2026-09-23): the topic is closed.** The core (min degree → gap; the
tail secondary), the finite size, and the qualitative localization picture were
published in 2003–2008. No new result can be obtained in this niche. The project is
stopped. All data/code obtained are kept as a reproduction and numerical detail of the
known result; if desired — as a starting point for another topic.

The old notes below are kept for history.

---

## Conclusion (outdated, see the update above)

**No exact match** (precisely `n^(−α)`, power-law config, min degree 2) **was found**
— preliminarily novel. **But** the very fact that `λ₂ → 0` in the sparse regime is
already known (see the ER precedent), so novelty could only be in the **rate α and its
dependence on the degree distribution**, not in the phenomenon itself.

## What is already known (baseline) — IMPORTANT for novelty

- **ISCCSP 2012**, DOI 10.1109/isccsp.2012.6217860: the normalized Laplacian giant of
  a **sparse ER** `G(n, d/n)` (d>1) has a **zero gap** as n→∞.
  ⇒ the qualitative decay is **not** power-law-specific, but a property of sparse
  graphs with min degree ≤ 2 (leaves + degree-2 chains).
- **Coja-Oghlan–Lanka**: the 2-core has a gap separated from 0.
- **Hermon–Li–Yao–Zhang**: for min degree ≥ 3 the relaxation time is Θ(1).
  ⇒ the d2→d3 transition = "switching off" the degree-2 structure.

## What could be new

If, at **equal average degree**, power-law decays faster than ER and this is not
explained by the degree-2 fraction — that would be the heavy-tail contribution. It is
tested by the ER control and the bounded-support controls (see below).

## Nearest works

| Work | What it proves | Relation to us |
|---|---|---|
| **Hermon, Li, Yao, Zhang**, arXiv:2105.11585, *Ann. Prob.* 50(5) 2022 | config model with degrees in **[3, d̄]**, spectral gap / relaxation time | **key**: requires min degree ≥3 — exactly our "const" regime; d_min=2 is outside |
| **Chung–Lu–Vu** (semicircle) | concentration when **d̄_min ≫ ln²n** | our regime d̄_min=O(1) — outside |
| **Wang, Li, Huang**, arXiv:2412.02087 (2024) | semicircle for the normalized Laplacian of the config model, "mild assumptions" | about the bulk spectrum, not the gap; heavy-tail+deg2 seems uncovered |
| **Coste**, arXiv:1708.00530 (2017) | spectral gap of sparse random **digraphs** by degree sequence (trace method) | methodology |
| **Cai, Caputo, Perarnau, Quattropani**, arXiv:2104.08389 | directed CM, heavy-tailed in-degrees, mixing, power-law stationary | heavy-tail context |
| **Fernley**, arXiv:2211.13537 | voter model on scale-free, mixing, τ∈(2,3] | "ultrasmall world" context |
| **Lu & Peng**, arXiv:1204.6207 | Laplacian concentration requires δ ≫ ln n | why sparse is hard |
| **Ducatez & Rivier**, EJP 2025, 10.1214/25-ejp1303 | gap of the (combinatorial) Laplacian of ER and embedded trees | methodology for ER; does not contradict the "zero gap" of the normalized Laplacian |
| **Avrachenkov et al.**, arXiv:1910.08869 | smallest nonzero eigenvalue → 0, power-law density near 0, spectral dimension | similar phenomenology |

## What the search did NOT cover (risks)

- Citation graph: who cites Coja-Oghlan–Lanka, Hermon et al., Bordenave — most reliable,
  but the arXiv API does not provide it.
- Google Scholar / Web of Science / MathSciNet.
- "Spectral dimension" scale-free queries might have missed an already-known exponent.

## Our ER control (`run_er_control.py`, `results_er/`, 2026-09-23)

ER `G(n, c/n)`, c ∈ {2,3,4}, n = 100…10⁵, reps = 20. Comparison by the average degree
of the 2-core (the ER core is denser than the whole graph):

| model | core d̄ | α(λ₂_core) |
|---|---|---|
| ER c=2 | 2.68 | 0.138 |
| ER c=3 | 3.44 | 0.089 |
| ER c=4 | 4.26 | 0.072 |
| pl_dmin2 γ=3.0 | 3.08 | 0.171 |
| pl_dmin2 γ=2.5 | 3.87 | 0.158 |
| pl_dmin2 γ=2.1 | 5.06 | 0.150 |

**Conclusion:** the qualitative decay is generic (ER too), but at equal d̄ the power-law
α is **~2× larger**. Confound: the degree-2 fraction is higher for pl. Bounded-support
controls are needed to separate "tail" from "degree-2 fraction".

## What is needed for a reliable novelty check

1. **M_b/M_a controls** (bounded support) — key (see CHANGELOG).
2. Targeted "who cites" search for Hermon–Li–Yao–Zhang / ISCCSP 2012.
3. Google Scholar / MathSciNet (not covered).
4. `d_max` test: vary the tail length at fixed n.
