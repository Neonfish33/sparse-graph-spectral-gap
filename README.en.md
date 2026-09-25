# sparse-graph-spectral-gap

> **Note:** the primary language of this project is **Russian** — the full detailed
> notes are in [`README.md`](README.md). This file is the concise English landing page.

**Status: closed project — no novel result.** A reproducible toolkit / benchmark
for the second eigenvalue `λ₂` of the **normalized Laplacian** of sparse random
graphs with heavy-tailed degrees. It **reproduces a known result** and is kept as
a numerical / benchmark artifact.

> This work empirically reproduces and numerically details a **known** result:
> **Samukhin, Dorogovtsev, Mendes**, *Laplacian spectra of complex networks and
> random walks on them: Are scale-free architectures really important?*,
> **Phys. Rev. E 77, 036115 (2008)**. For sparse random graphs, the lower edge of
> the Laplacian spectrum is determined by the **minimum degree** (`λ_c > 0` iff
> `q_min > 2`; `λ_c = 0` iff `q_min ≤ 2`), and the heavy tail is secondary. The
> finite-size behavior is also known (`λ₂(N) ∼ (ln N)^{−2}` for `q_min ≤ 2`).
> We confirm the picture numerically up to `n = 10⁶`. No new result was obtained.
> See [`LITERATURE.md`](LITERATURE.md) and [`REVIEW.md`](REVIEW.md).

## What's inside

- `src/generators.py` — degree sequences, configuration model, four ensembles.
- `src/diagnostics.py` — `λ₂` (full / giant / 2-core), sweep cut, conductance
  (`φ_periph`, `φ_mod`, `φ*`), `Δ`, `ρ*`, IPR, φ-curve.
- Scripts: `run_experiments.py`, `analyze_results.py`, `analyze_gap_fss.py`,
  `validate_power_law.py`, `run_er_control.py`, `run_toy_controls.py`,
  `toy_model.py`, `compare_laplacians.py`, `compare_asymptotics.py`,
  `merge_diagnostics.py`, `export_graph.py`.
- `tests/` — 19 tests (spectra, sweep cut vs brute force, Cheeger, IPR).

### Data
Included in git (~13 MB): `results_*/diagnostics.csv` and `results_*/analysis/**`
(so the analysis is reproducible without regenerating graphs). The large graphs
(`results_*/graphs/`, 231 MB / 13500 `.npz`) are gitignored and can be regenerated
deterministically from the seeds stored in the CSVs.

## Quick start

```bash
pip install -r requirements.txt
python -m pytest tests/ -q                     # run tests
python run_experiments.py --quick              # smoke test
python analyze_results.py --csv results_merged/diagnostics.csv --out results_merged/analysis
python validate_power_law.py                   # AICc: power law vs log vs exp-log
python compare_laplacians.py                   # normalized vs combinatorial λ₂
```

## Documentation

Each document exists in English (`*.en.md`) and Russian. The **Russian originals are
authoritative**; the English ones are translations.

- [`REVIEW.en.md`](REVIEW.en.md) / [`REVIEW.md`](REVIEW.md) — full summary, findings, open questions.
- [`LITERATURE.en.md`](LITERATURE.en.md) / [`LITERATURE.md`](LITERATURE.md) — literature check; why there is no novelty.
- [`CHANGELOG.en.md`](CHANGELOG.en.md) / [`CHANGELOG.md`](CHANGELOG.md) — change history.
- [`TOY_MODEL.en.md`](TOY_MODEL.en.md) / [`TOY_MODEL.md`](TOY_MODEL.md) — toy-model design.
- [`README.md`](README.md) — detailed Russian project notes (the primary README).

## License

MIT — see [`LICENSE`](LICENSE).
