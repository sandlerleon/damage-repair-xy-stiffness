# -*- coding: utf-8 -*-
# Section 3.4: timescales and sensitivity (exec'd by build_manuscript.py)
TS = R["timescale"]
CT = TS["cells"]
CO = TS["collapse"]
SL = TS["slow_limit"]
TP = R["tauphi"]
TSS = R["timescale_s"]
SN = R["sensitivity"]


def U_(mode, rho, lam, rel):
    return CT["%s|%g|%g|%d" % (mode, rho, lam, rel)]["U"][0]


HD("3.4 Timescales: is the repair-to-damage ratio enough?", 2)
tphi_h = [TP["0.4|%s|%g" % (m, r)]["tau"][0] for m in ("keep", "reset0") for r in (1.0, 2.5, 5.0, 10.0)]
tphi_f = [TP["0|keep|%g" % r]["tau"][0] for r in (5.0, 10.0)]
P("Changing the damage probability λ_{0} at fixed ρ changes *U* appreciably, so ρ does not collapse the lattice data at a fixed number of phase "
  "sweeps per step. We decompose this dependence (Figure 6, @T:slow@ and @T:sens@). Three effects are involved. First, the kinetic update is discrete: at fixed ρ "
  "the per-step probabilities *R* = ρλ_{0} and λ_{0}μ grow with λ_{0}, and the hazard-equivalent ratio ln(1 − *R*)/ln(1 − λ_{0}) exceeds ρ by a factor "
  "%s at λ_{0} = 0.04 and %s at λ_{0} = 0.08 (where, in addition, λ_{0}(1 + β) = %s > 1 and the kill probability saturates). Second, the rule "
  "θ → 0 injects order at rate *R* (Section 3.3). Third, there is a genuine competition between kinetics and phase relaxation, measured by "
  "Ω = λ_{0}/*relax*, the number of damage events per site per phase sweep, which is proportional to the ratio of the phase-relaxation "
  "time to the damage time. The relaxation time of the phase field in a frozen damage landscape, measured as the integrated autocorrelation time of *U*, is "
  "%s–%s sweeps in the field (*h*_{0} = 0.4) and %s–%s sweeps for the global phase of the coherent field-free state."
  % (f3(CO["keep"]["rho_hazard_over_rho_at_lam0.04"], 2), f3(CO["keep"]["rho_hazard_over_rho_at_lam0.08"], 2), f3(0.08 * 21, 2),
     f3(min(tphi_h), 1), f3(max(tphi_h), 1), f3(min(tphi_f), 0), f3(max(tphi_f), 0)))
P("We therefore scanned λ_{0} from 0.005 to 0.08 and the number of phase sweeps per kinetic step from 1 to 8 (Ω from %s to %s) for the three repair "
  "rules at ρ = 2.5 and 5 (*h*_{0} = 0.4, β = 20, *L* = 32, 8 seeds, 3000 + 1500 steps). At fixed *relax*, raising λ_{0} from 0.005 to 0.04 at ρ = 5 raises *U* "
  "from %s to %s for *reset0* and from %s to %s for *neighbor*, but leaves it almost unchanged for *keep* (%s to %s); at fixed λ_{0}, "
  "increasing *relax* lowers *U* for *reset0* and *neighbor* and raises it for *keep*. The dependence of *U* on (ρ, λ_{0}, *relax*) is "
  "resolved by plotting *U* against the damaged fraction *q* that each run reaches (Figure 6c). For *neighbor*, *U* is a function of *q* "
  "alone: the rms deviation of the 32 unsaturated points from a single cubic curve *U*(*q*) is %s, with no trend in Ω (%s per decade). Its "
  "dependence on λ_{0} at fixed ρ is therefore carried entirely by *q*, that is by the discrete-time effect above. For *reset0* the deviation is "
  "%s with a small trend (%s per decade). For *keep* it is %s, and *U* at a given *q* rises by %s per decade as Ω falls (correlation %s with log Ω): a "
  "repaired site retains the random phase it had while damaged and needs sweeps to rejoin the coherent background, so that rapid kinetics, "
  "relative to the relaxation, keep it incoherent; a site repaired to the phase of its neighbors has nothing to relax."
  % (sci(0.005 / 8, 1), "%.2f" % 0.04, f3(U_("reset0", 5.0, 0.005, 2), 3), f3(U_("reset0", 5.0, 0.04, 2), 3),
     f3(U_("neighbor", 5.0, 0.005, 2), 3), f3(U_("neighbor", 5.0, 0.04, 2), 3), f3(U_("keep", 5.0, 0.005, 2), 3), f3(U_("keep", 5.0, 0.04, 2), 3),
     f3(CO["neighbor"]["master_q_rms"], 3), "%+.4f" % CO["neighbor"]["resid_slope_per_decade_omega"], f3(CO["reset0"]["master_q_rms"], 3),
     "%+.4f" % CO["reset0"]["resid_slope_per_decade_omega"], f3(CO["keep"]["master_q_rms"], 3), f3(abs(CO["keep"]["resid_slope_per_decade_omega"]), 3),
     f3(CO["keep"]["resid_vs_logomega_corr"], 2)))
rows = [["repair rule", "ρ = 2.5, Ω = 10^{−2}", "ρ = 2.5, Ω = 6 × 10^{−4}", "ρ = 5, Ω = 10^{−2}", "ρ = 5, Ω = 6 × 10^{−4}"]]
for mode, nm in (("reset0", "reset0 (θ → 0)"), ("keep", "keep"), ("neighbor", "neighbor")):
    rows.append([nm] + [pm(SL[r][mode][k]) for r in ("2.5", "5") for k in ("U_base", "U_slow")])
rows.append(["spread over the rules", f3(SL["2.5"]["spread_base"]), f3(SL["2.5"]["spread_slow"]), f3(SL["5"]["spread_base"]), f3(SL["5"]["spread_slow"])])
TAB(rows, "Coherence *U* (mean ± s.e., 8 seeds, *L* = 32, β = 20, *h*_{0} = 0.4) at the baseline kinetics (λ_{0} = 0.02, two phase sweeps per "
          "step, Ω = 0.01) and in the slowest-damage state scanned (λ_{0} = 0.005, eight sweeps, Ω = 6 × 10^{−4}), for the three repair rules.",
    widths=[1.5, 1.2, 1.4, 1.2, 1.4], label="slow")
P("The central result is the last row of @T:slow@. At the baseline kinetics the three repair rules differ in *U* by %s at ρ = 2.5 and %s at ρ = 5; at the "
  "slowest damage scanned they differ by %s and %s, a tenfold convergence (Figure 6d, in which the spread falls roughly in proportion to Ω), and "
  "*U* at ρ = 5 is %s, %s and %s for *reset0*, *keep* and *neighbor*. The numerical evidence therefore supports the repair-to-damage ratio becoming the controlling parameter as Ω falls: "
  "over the range scanned the dependence on the repair rule and on λ_{0} separately shrinks roughly in proportion to Ω, and at the slowest damage it is small. "
  "The mathematical limit Ω → 0 lies beyond the range tested, and we do not claim it; the dependences found at finite Ω are the finite-rate corrections of the three effects above. Two consequences follow. The values at the baseline "
  "kinetics, such as @T:steady@, are values at a finite Ω = 0.01 and lie up to %s above the slowest-damage value at ρ = 5 for *reset0*; and the equilibrium "
  "fluctuation formula for the stiffness, which fails at finite Ω (Section 3.3), should approach the twist response as Ω falls. Over the range scanned it does: for "
  "the field-free *keep* state at ρ = 10 and λ_{0} = 0.005 the ratio Υ_{eq}/Υ_{tw} is %s, %s and %s at Ω = %s, %s and %s (Figure 6f), the "
  "stiffness itself changing by only %s over the same range. In a conductor the damage and the annealing proceed over hours to years, while the "
  "relaxation of the superconducting phase is expected to be many orders of magnitude faster (we do not quantify it here), so that experiments "
  "should lie at much smaller Ω than scanned, where ρ would be the controlling quantity if the trend continues; that extrapolation is an assumption."
  % (f3(SL["2.5"]["spread_base"], 2), f3(SL["5"]["spread_base"], 2), f3(SL["2.5"]["spread_slow"], 3), f3(SL["5"]["spread_slow"], 3),
     f3(SL["5"]["reset0"]["U_slow"][0], 2), f3(SL["5"]["keep"]["U_slow"][0], 2), f3(SL["5"]["neighbor"]["U_slow"][0], 2),
     f3(SL["5"]["reset0"]["U_base"][0] - SL["5"]["reset0"]["U_slow"][0], 2),
     f3(TSS["0.005|1"]["ups_eq"][0] / TSS["0.005|1"]["ups_tw"][0], 2), f3(TSS["0.005|2"]["ups_eq"][0] / TSS["0.005|2"]["ups_tw"][0], 2),
     f3(TSS["0.005|8"]["ups_eq"][0] / TSS["0.005|8"]["ups_tw"][0], 2), sci(TSS["0.005|1"]["omega"]), sci(TSS["0.005|2"]["omega"], 1),
     sci(TSS["0.005|8"]["omega"], 1), "%d%%" % round(100 * abs(TSS["0.005|8"]["ups_tw"][0] / TSS["0.005|1"]["ups_tw"][0] - 1))))
FIG("Fig6_timescale.png",
    "Is ρ enough? (a, b) *U* against Ω = λ_{0}/*relax* for the three repair rules (*h*_{0} = 0.4, β = 20, *L* = 32; every combination of λ_{0} = "
    "0.005–0.08 and *relax* = 1–8; open symbols: λ_{0} = 0.08, whose kill probability saturates) at ρ = 5 and 2.5. (c) The same points against the "
    "damaged fraction *q* reached. (d) Spread of *U* across the three rules at each (λ_{0}, *relax*). (e) Field-free stiffness Υ_{tw} (filled) and the "
    "equilibrium formula Υ_{eq} (open) at ρ = 10 for λ_{0} = 0.005, 0.02 and 0.08 (*relax* = 1, 2, 8 along each line). (f) Υ_{eq}/Υ_{tw} for "
    "λ_{0} ≤ 0.02.")
# ---- sensitivity
P("@T:sens@ gives the sensitivity of *U* at the baseline kinetics to the other parameters, one at a time, with 10 seeds (6 at *L* = 96). *U* is "
  "independent of *L* from 24 to 96 within error, and depends strongly on temperature, on the field *h*_{0} and the coupling *g*_{0} of damaged sites, and "
  "on β (@T:sens@).")
rows = [["parameter", "values", "*U* at ρ = 2.5", "*U* at ρ = 5"]]
for name, vals, lab in (("T", [0.25, 0.35, 0.5, 0.7], "temperature *T*"), ("g0", [0.0, 0.03, 0.12, 0.3, 0.5], "damaged coupling *g*_{0} (*g*_{0} = 0: *U* over all sites)"),
                        ("h0", [0.0, 0.2, 0.4, 0.8], "field *h*_{0} (θ → 0 repair)"), ("beta", [0.0, 5.0, 10.0, 20.0, 30.0], "feedback β"),
                        ("L", [24, 32, 48, 64, 96], "lattice size *L*")):
    rows.append([lab, ", ".join("%g" % v for v in vals)] +
                [" / ".join(f3(SN["%s|%g|%g" % (name, v, rho)]["U"][0], 3) for v in vals) for rho in (2.5, 5.0)])
TAB(rows, "Sensitivity of the steady-state *U* (baseline λ_{0} = 0.02, two sweeps per step, repair θ → 0, *h*_{0} = 0.4 unless varied; mean over 10 seeds "
          "(6 for *L* = 96); the standard errors are at most %s). Values are listed in the order of the parameter values." %
    f3(max(v["U"][1] for v in SN.values()), 3), widths=[2.0, 1.4, 1.6, 1.6], size=8, label="sens")
