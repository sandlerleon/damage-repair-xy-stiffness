# Damage–repair kinetics and phase stiffness in a site-diluted XY model

Leon Sandler, Independent Researcher — sandler.leon@gmail.com
ORCID [0009-0007-4584-808X](https://orcid.org/0009-0007-4584-808X)

This repository holds the kinetic Monte Carlo study behind the manuscript *"Damage–repair kinetics and phase stiffness
in a site-diluted XY model: implications for irradiated superconductors"* (prepared for **Superconductor Science and
Technology**; preprint [10.5281/zenodo.23173924](https://doi.org/10.5281/zenodo.23173924), manuscript v4). Every number and every figure
of the manuscript is produced by the scripts here, from the raw run outputs that are also here.

## The model

A 2D XY lattice (phases θᵢ, couplings weakened to g₀ = 0.12 on damaged sites, T = 0.35) in which sites are damaged at a
rate λᵢ = λ₀[1 + β Dⁿ/(Kⁿ + Dⁿ)] that rises with the local phase disorder D (h = 6, K = 0.5, β = 20, λ₀ = 0.02) and
repaired at rate R = ρλ₀. It is a phenomenological model of competing displacement damage and recovery in an irradiated
superconductor. It has no pinning and is **not calibrated** to any conductor.

## What the Monte Carlo shows

* **The mean-field reduction predicts bistability, the lattice shows none in the ranges tested.** Re-derived for the damaged
  fraction, the reduction predicts a bistable window. On the lattice, two starts (ordered and 95 % damaged) agree to within
  0.0027 in the coherence over 50 conditions (L = 32–64, up to 10 seeds); 9 intervals exclude zero, none survives a Holm
  correction (smallest adjusted p = 0.19). Quasi-static ramps close. The reason is that the closure μ̄ = μ(q) fails: the
  kill multiplier the lattice applies is a smooth function of the damaged fraction.
* **Field-free stiffness.** Measured as the *twist response* Υ_tw of the non-equilibrium steady state with the winding number
  recorded (φ₀ = 0.3π/L), repair that preserves the phase symmetry restores a stiffness above a threshold ratio ρ: a finite-size estimate
  between 7.00 and 7.75 for repair that leaves the phase unchanged (`keep`, L up to 96) and between 5.50 and 6.00 for repair to the
  phase of the neighbors (`neighbor`, L up to 64), against ρ ≈ 2.5 without state-dependent damage. These are brackets from
  operational cutoffs, not thermodynamic-limit bounds. The twist response is validated (`study.py validation`) against the equilibrium
  formula on clean and diluted lattices, imposed windings, and an independent static code. The equilibrium helicity-modulus formula underestimates Υ_tw by up to about 30 % and goes negative where the
  stiffness vanishes.
* **Is ρ enough?** Increasingly so as damage becomes slow compared with phase relaxation: the three repair rules differ by 0.12 in U at
  the finite damage rates tested and agree within 0.012 at the slowest damage scanned (the mathematical limit is not claimed).
* **Static dilution.** Quenched random dilution loses stiffness at a removed fraction f ≈ 0.25 ± 0.01 (sizes 16–96, finite-size
  extrapolation and its sensitivity reported); compact damage clusters tolerate more.
* **Link to superconductors (illustrative).** With T_BKT(f)/T_BKT(0) mapped to Tc/Tc0 and a REBCO Tc decrease per dpa taken from the
  literature, an effective conversion of about 10 model sites per dpa results (no further physical interpretation); an Arrhenius annealing analysis is
  compatible with a broad activation spectrum, whose width the two reported recoveries do not determine. No calibration is claimed.

## Layout

```
code/xy_engine.py     the engine: Metropolis rotor sweeps + synchronous damage/repair (reset0 | keep | neighbor), twist, winding,
                      C(r), quenched and clustered dilution; update rules reproduce MyUncle bit for bit
code/test_engine.py   15 checks (run before any result)
code/study.py         all studies (steady, converge, branches, ramps, fieldfree, stiffness, stiff_L*, timescale*, sensitivity,
                      tauphi, closure0*, fss, cluster, tcurve)
code/analyze*.py      raw runs -> results/results.json (analyze_e.py: the twist-response validation)
code/meanfield.py     the mean-field reduction (windows, beta_c)
code/figures.py       Figures 1–8 -> figures/
results/*_raw.json    raw outputs of every run (results/run_log*.txt = logs)
manuscript/           manuscript v1–v4 (v4: language-edited; merge_rubriq.py), cover letter, and the build scripts that read results.json
tools/                Zenodo deposit scripts (token read from ZENODO_TOKEN, never stored)
```

`results/stiffness_phi0.1_deprecated_raw.json` is a first stiffness run with a twist angle that was too large for the larger
lattices (φ₀ must be ≪ π/L); it is kept for the record and is not used.

## Reproducing

```bash
pip install -r requirements.txt
cd code
python test_engine.py
python study.py steady converge branches ramps fieldfree stiffness   # the study list is in study.py's docstring; --quick for a smoke test
python analyze.py && python analyze_b.py && python analyze_c.py && python analyze_d.py
python figures.py fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8
```

The full programme takes several days on one workstation. Each run is seeded from `SeedSequence([MASTER, study, ...])`, so a
single run can be repeated alone. Intervals are 95 % percentile bootstrap intervals over independent runs or disorder
realizations; autocorrelation times use Sokal's automatic window.

## Citation

Software: v1.2.1 (manuscript v4) [10.5281/zenodo.23173923](https://doi.org/10.5281/zenodo.23173923); v1.2.0 (manuscript v3) [10.5281/zenodo.23173811](https://doi.org/10.5281/zenodo.23173811); v1.1.0 (manuscript v2) [10.5281/zenodo.23149157](https://doi.org/10.5281/zenodo.23149157); v1.0.0 (manuscript v1) [10.5281/zenodo.23138195](https://doi.org/10.5281/zenodo.23138195); concept DOI, always the latest, [10.5281/zenodo.23138194](https://doi.org/10.5281/zenodo.23138194). Update rules of the underlying framework:
[MyUncle](https://github.com/sandlerleon/MyUncle) ([10.5281/zenodo.21223569](https://doi.org/10.5281/zenodo.21223569)).

MIT licence.
