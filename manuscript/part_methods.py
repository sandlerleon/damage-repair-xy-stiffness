# -*- coding: utf-8 -*-
# Section 2: model and numerical methodology (exec'd by build_manuscript.py)
i_, j_ = V("i"), V("j")
TH, LM, RHO, MU, BETA = V("θ"), V("λ"), V("ρ"), V("μ"), V("β")
EQ_N, MI_N, PL_N = T(" = "), T(" − "), T(" + ")
sub0 = lambda sym, s: SUB(V(sym), T(s))
ST = R["steady"]
LAM0 = ST["lam0"]

HD("2. Model and numerical methodology")
HD("2.1 Hamiltonian and phase dynamics", 2)
P("An *L* × *L* square lattice with periodic boundaries carries a phase θ_{i} ∈ (−π, π] on each site, representing the phase of the local "
  "superconducting order parameter, and an integrity *g*_{i}. A site is intact (*g* = 1) or damaged (*g* = *g*_{0} = 0.12); a damaged "
  "site is one whose pair amplitude has been suppressed by displacement damage (pair breaking), so that it couples only weakly to its "
  "neighbors. The energy is")
EQ(V("H") + EQ_N + T("−") + sub0("J", "0") + T(" ") + T("∑") + T("⟨") + i_ + j_ + T("⟩") + T(" ") + Vs("g", "i") + Vs("g", "j") + T("cos") + DELIM(Vs("θ", "j") + MI_N + Vs("θ", "i") + MI_N + Vs("φ", "ij"))
   + MI_N + sub0("h", "0") + T(" ") + T("∑") + SUB(V("g"), i_) + T("cos") + T(" ") + SUB(TH, i_) + T(","), 1)
P("with *J*_{0} = 1 and *k*_{B} = 1 setting the units of energy and temperature, and *T* = 0.35 throughout unless stated. The twist φ_{ij} is "
  "zero except for bonds in the +*x* direction, where it takes the uniform value φ (a vector potential; Section 2.4). The term "
  "*h*_{0} breaks the global U(1) phase symmetry explicitly. Most results below are for *h*_{0} = 0, so that coherence is spontaneous; "
  "*h*_{0} = 0.4 is the field of the superconducting preset of MyUncle [[sandler_myuncle]] and is used only for the baseline scans. Phases are updated by checkerboard "
  "Metropolis sweeps with uniform proposals on (−π, π] at fixed integrity, two sweeps (*relax* = 2) between successive kinetic updates.")

HD("2.2 Damage and repair", 2)
P("After the phase sweeps each intact site is damaged, with probability per step")
EQ(SUB(LM, i_) + EQ_N + sub0("λ", "0") + T(" ") + MU + DELIM(SUB(V("D"), i_)) + T(",") + T("  ") + MU + DELIM(V("D")) + EQ_N + T("1") + PL_N + BETA
   + FRAC(SUP(V("D"), V("h")), SUP(V("K"), V("h")) + PL_N + SUP(V("D"), V("h"))) + T(","), 2)
P("where D_{i} = 1 − |¼ Σ_{j} exp(iθ_{j})| ∈ [0, 1], the sum running over the four nearest neighbors of site *i*, is the local phase disorder "
  "(the neighbor phases are measured in the gauge of the twist), *h* = 6 and *K* = 0.5. Damage is therefore faster where phase order is already disrupted; β = 0 is state-independent damage and "
  "β = 20 is the baseline feedback. Every damaged site is then repaired with probability *R* in the same step, so a site damaged in a step "
  "can be repaired in it. Repair restores *g* = 1 and treats the phase of the repaired site by one of three rules: *reset0*, "
  "θ → 0 (the rule of the MyUncle superconducting preset); *keep*, the phase is left unchanged and the "
  "dynamics relaxes it; or *neighbor*, the phase is set to that of the coupling-weighted neighboring phasors, so that the repaired site "
  "rejoins the coherent region around it. The first rule selects a phase that is common to all sites and therefore acts, at rate *R*, as a "
  "symmetry-breaking field even when *h*_{0} = 0; the other two preserve the U(1) symmetry (this is tested, Section 2.4). The control "
  "ratio is ρ = *R*/λ_{0}, with λ_{0} = 0.02, and one kinetic step is the unit of time. All rates are probabilities per step and must not "
  "exceed 1: for λ_{0} = 0.02 the largest kill probability is λ_{0}(1 + β) = %s < 1. In the scan over λ_{0} (Section 3.4) the value "
  "λ_{0} = 0.08 gives λ_{0}(1 + β) = %s > 1, so the kill probability saturates at 1 where *D*_{i} ≳ 0.53, and the corresponding points are "
  "flagged." % (f3(LAM0 * 21, 2), f3(0.08 * 21, 2)))

HD("2.3 Mean-field reduction", 2)
P("The reduction is derived from the update just defined. Let *q* be the damaged fraction (*g* = *g*_{0}), *p* = 1 − *q* the intact fraction, "
  "and μ̄ = ⟨μ(*D*_{i})⟩ the mean kill multiplier over the intact sites. Because damage is applied first and repair second, "
  "one step maps *q* to")
EQ(V("q") + T("′") + EQ_N + DELIM(T("1") + MI_N + V("R")) + DELIM(V("q") + PL_N + DELIM(T("1") + MI_N + V("q")) + sub0("λ", "0") + V("μ̄")), 3)
P("For λ_{0}, *R* ≪ 1 and τ = λ_{0}*t* this is the rate equation")
EQ(FRAC(T("d") + V("q"), T("d") + V("τ")) + EQ_N + V("μ̄") + DELIM(T("1") + MI_N + V("q")) + MI_N + RHO + V("q") + T("."), 4)
P("The variable is the *damaged* fraction. For β = 0, μ̄ = 1 and *q*^{*} = 1/(1 + ρ); the intact fraction is ρ/(1 + ρ), and it is the "
  "latter that increases with ρ like the order *U*; the damaged fraction 1/(1 + ρ) must not be used as a reference for *U*, which is "
  "not a comparison of like with like. The only approximation is the closure of μ̄, in which the mean over the intact sites of the "
  "local-disorder function is replaced by its value at the damaged fraction, μ̄ = μ(*q*) (the annealed identification *D*_{i} → *q*, "
  "motivated by the slaving of the phase order to the intact fraction, *U* ≈ 1 − *q*). Steady states then satisfy")
EQ(RHO + EQ_N + FRAC(MU + DELIM(V("q")) + DELIM(T("1") + MI_N + V("q")), V("q")) + T(","), 5)
WC = ST["windows"]
P("which is S-shaped, and hence bistable, for cooperativity *h* > 1 and gain above a threshold β_{c}(*h*, *K*) (β_{c} = %s at *h* = 6, *K* = 0.5; "
  "there is no bistability for *h* ≤ 1). For (*h*, *K*, β) = (6, 0.5, 20) the window of equation (5) is %s < ρ < %s, and %s < ρ < %s for β = 5. "
  "Equation (5) is the continuous-time limit; the map (3) that the code implements gives, from its fixed points ρ = (1 − *q*)μ/[*q* + (1 − *q*)λ_{0}μ], "
  "the windows %s < ρ < %s (β = 20) and %s < ρ < %s (β = 5). The shift is appreciable because λ_{0}μ reaches %s per step at high disorder, "
  "where the continuous-time limit is poor. The closure is tested against the lattice in Section 3.1."
  % (f3(ST["beta_c_h6"], 2), f3(WC["20"]["continuous"][0], 2), f3(WC["20"]["continuous"][1], 2), f3(WC["5"]["continuous"][0], 2),
     f3(WC["5"]["continuous"][1], 2), f3(WC["20"]["discrete"][0], 2), f3(WC["20"]["discrete"][1], 2), f3(WC["5"]["discrete"][0], 2),
     f3(WC["5"]["discrete"][1], 2), f3(LAM0 * 21, 2)))

HD("2.4 Observables, protocols and uncertainty", 2)
P("**Order and damage.** The coherence is *U* = |⟨exp(iθ)⟩| over all sites; the damaged fraction *q*, the kill and repair currents "
  "(events per site and step), and the mean local disorder and kill multiplier of the intact sites are recorded. "
  "**Phase stiffness.** The helicity modulus of an equilibrium system is the second derivative of the free energy with respect to a "
  "uniform twist; in a non-equilibrium state with fluctuating damage there is no free energy and the equilibrium fluctuation formula "
  "cannot be assumed. We therefore measure the *twist response*: the linear response of the stationary bond current to a uniform "
  "vector potential φ applied along *x*,")
EQ(sub0("Υ", "tw") + EQ_N + T("−") + FRAC(DELIM(Vs("j", "x") + T("(+") + sub0("φ", "0") + T(")") + MI_N + Vs("j", "x") + T("(−") + sub0("φ", "0") + T(")")), T("2") + sub0("φ", "0")) + T(","), 6)
P("with *j*_{x}(φ) = ⟨*N*^{−1} Σ_{⟨ij⟩∥x} *J*_{ij} sin(θ_{j} − θ_{i} − φ)⟩ the stationary mean current per site, φ_{0} = 0.3π/*L* (%s at *L* = 32; see below), and the "
  "same random-number stream for both signs. It requires no detailed balance and measures what a superconductor measures: the supercurrent "
  "that the damaged-and-repairing medium sustains per unit imposed phase gradient (it is proportional to the superfluid density, and to "
  "*d*/λ_{L}^{2} for a film of thickness *d* and London depth λ_{L}). In equilibrium it equals the helicity modulus. For comparison we also "
  "evaluate the equilibrium fluctuation formula on the same stationary runs," % f3(json.load(open(os.path.join(RESDIR, "stiffness_raw.json")))["meta"]["phi0"], 2))
EQ(sub0("Υ", "eq") + EQ_N + FRAC(T("1"), V("N")) + T("⟨") + T("∑") + T(" ") + sub0("J", "ij") + T("cos") + T(" ") + V("Δ") + TH + T("⟩") + MI_N + FRAC(V("N"), V("T")) + DELIM(T("⟨") + SUP(Vs("j", "x"), T("2")) + T("⟩") + MI_N + SUP(T("⟨") + Vs("j", "x") + T("⟩"), T("2"))) + T("."), 7)
P("Phase correlations are *C*(*r*) = ⟨cos(θ_{i} − θ_{i+*r*})⟩ over all pairs at separation *r* along *x* and *y*, and, separately, over pairs of "
  "intact sites.")
P("**Protocols.** Unless stated, runs start from the intact state (and, for tests of bistability, also from a state with 95% of the sites "
  "damaged and random phases), equilibrate for 3000 steps (2000 in the baseline scan; both are over 30 times the longest autocorrelation "
  "time measured) and are measured for 3000–4000 steps, sampled every 5–10 steps. Each point is the mean over 6–10 independent seeds; "
  "the error is the standard error over seeds, and intervals are 95% percentile bootstrap intervals (2000 resamples) over seeds or disorder "
  "realizations [[efron1979]]. Integrated autocorrelation times use the automatic window of Sokal [[madras1988]] and are reported, but "
  "never used in place of the seed-to-seed spread. Comparisons of two starting states over many conditions are corrected for multiple "
  "comparisons by Holm’s procedure. Every run is seeded from a SeedSequence of a master seed and a code naming the study, parameters, "
  "size, replicate and initial state, so any run can be repeated alone. The engine reproduces the update rules of the MyUncle framework "
  "[[sandler_myuncle]] bit for bit for the same seed (the superconducting preset, three parameter sets), and its test suite also "
  "verifies the discrete two-state fixed point at β = 0, the repair rules, the gauge consistency of the twisted lattice, "
  "the agreement of Υ_{tw} with the equilibrium formula on a clean lattice, and that with *h*_{0} = 0 the two symmetric repair rules "
  "leave the global phase undetermined while *reset0* selects phase 0 (|⟨exp(iψ)⟩| = 1.00).")
P("**Static dilution.** For the field-free transition of a lattice with quenched random site removal (removed sites have *g* = 0), the helicity "
  "modulus is computed from equation (7), which is legitimate in equilibrium, and compared with the Nelson–Kosterlitz value 2*T*/π [[nk1977]]. "
  "The crossing fraction *f*_{KT}(*L*) is located by linear interpolation of the realization-averaged Υ(*f*) and extrapolated in "
  "1/ln²*L* and by a Weber–Minnhagen form [[wm1988]]. The uncertainty of each crossing (a bootstrap over disorder realizations) is "
  "distinguished throughout from the scatter of the crossings over sizes, which is not an error estimate and no substitute for the "
  "extrapolation.")
HD("2.5 Validation of the twist response", 2)
VAL = R["validation"]
_vc, _vd, _vl, _vw = VAL["clean"], VAL["diluted"], VAL["diluted_long"], VAL["wound"]
_Ls = (16, 32, 64)
_cs = ("0.1", "0.3", "0.6", "1", "1.5")
_rowsv = [["*L*"] + ["φ_{0} = %s π/*L*" % c for c in _cs]]
for _L in _Ls:
    _rowsv.append([str(_L)] + [pm(_vc["%d|%s" % (_L, c)]["ups_tw"], 3) for c in _cs])
_dmax = max(abs(_vc["%d|0.3" % L]["ups_tw"][0] - _vc["%d|0.3" % L]["ups_eq"][0]) for L in _Ls)
spm = lambda v: ("%+.3f ± %.3f" % (v[0], v[1])).replace("-0.000", "0.000").replace("+0.000", "0.000")
_slips = sum(v["n_flagged"] for v in _vc.values())
P("The twist response is the central measurement of this paper, and the equilibrium formula that it replaces fails in the driven state, so it was "
  "validated in control simulations with the same engine and observables, with damage and repair switched off (λ_{0} = *R* = 0), at *T* = %s and "
  "*h*_{0} = 0. **Normalization and linearity.** On the clean lattice, where equation (7) is exact, the twist response at φ_{0} = 0.3π/*L* is "
  "%s, %s and %s for *L* = 16, 32 and 64 (8 seeds), equal to the equilibrium formula evaluated on the same runs to within %s and close to the first-order "
  "spin-wave value 1 − *T*/4 = %s; the normalization 2φ_{0} and the sign convention are therefore right, and Υ_{tw} does not depend on *L*. It "
  "varies little with the twist amplitude (@T:valid@): the decrease with φ_{0} is the anharmonic (cosine) nonlinearity of the response, and is "
  "visible only for the smallest lattice and φ_{0} > 0.6π/*L*. No clean-lattice run changed its winding number (%d slips in %d pairs); phase slips matter in states "
  "with vortices, and an earlier run of the damage–repair states with a fixed φ_{0} = 0.1, which exceeds π/*L* for *L* ≥ 32, was contaminated by "
  "them and was discarded in favor of φ_{0} = 0.3π/*L*."
  % (f3(VAL["T"], 2), pm(_vc["16|0.3"]["ups_tw"], 3), pm(_vc["32|0.3"]["ups_tw"], 3), pm(_vc["64|0.3"]["ups_tw"], 3), f3(max(_dmax, 0.001), 3), f3(VAL["harmonic_clean"], 4),
     _slips, sum(v["n"] for v in _vc.values())))
TAB(_rowsv, "Twist-response stiffness Υ_{tw} (mean ± s.e., 8 seeds) of the clean lattice (no damage or repair, *T* = 0.35) against the twist angle φ_{0} in units of π/*L*.",
    widths=[0.6, 1.2, 1.2, 1.2, 1.2, 1.2], size=8, label="valid")
P("**Windings.** A frozen global winding number *W* carries the current Υ(2π*W*/*L* − φ), so a pair of runs with different *W* gives a wrong response unless it is "
  "corrected by the integer difference. This was tested by preparing the two runs of each pair with imposed windings (*L* = 32, %d seeds): with *W* = 0 in both, "
  "Υ_{tw} = %s; with *W* = 1 in both, %s (a state with *W* = 1 is twisted by 2π/*L* per row, which lowers the response by about 2%%, as the cosine nonlinearity, cos(2π/*L*) = %s, suggests); with *W* = 1 in one run and "
  "0 in the other, the uncorrected value is %s and the corrected value, which is the one used throughout the damage–repair data (Section 3.3), is %s, within 1%% of the clean value. "
  "**Independent equilibrium code.** On static randomly diluted lattices (*L* = 32, *f* = 0.1 and 0.2, 16 coupling maps, runs of 12000 steps after 3000) the twist response "
  "of the dynamic engine, on the same coupling maps, is %s and %s, against %s and %s from the equilibrium helicity modulus computed by the separate static Metropolis "
  "code used in Section 3.5 (paired differences %s and %s) and %s and %s from the fluctuation formula on the same runs. With runs a quarter as long (8 maps) the response "
  "at *f* = 0.2 was lower than the fluctuation formula by %s, a sampling effect that disappears in the longer runs, which is a reminder that the stiffness of near-critical "
  "states needs long runs; the damage–repair data use 6000 measured steps and their intervals are over seeds. The twist response and the equilibrium formula agree in "
  "equilibrium, as they should, and differ in the driven damage–repair state, which is the result of Section 3.3."
  % (_vw["w00"]["n"], pm(_vw["w00"]["ups_tw"], 3), pm(_vw["w11"]["ups_tw"], 3), f3(math.cos(2 * math.pi / 32), 3), pm(_vw["w10"]["ups_uncorrected"], 2),
     pm(_vw["w10"]["ups_tw"], 3), pm(_vl["0.1"]["ups_tw"], 3), pm(_vl["0.2"]["ups_tw"], 3), pm(_vl["0.1"]["ups_eq_static"], 3), pm(_vl["0.2"]["ups_eq_static"], 3),
     spm(_vl["0.1"]["diff_tw_minus_static"]), spm(_vl["0.2"]["diff_tw_minus_static"]), pm(_vl["0.1"]["ups_eq_same_run"], 3), pm(_vl["0.2"]["ups_eq_same_run"], 3),
     pm((-_vd["0.2"]["diff_tw_minus_eq_same_run"][0], _vd["0.2"]["diff_tw_minus_eq_same_run"][1]), 3)))
