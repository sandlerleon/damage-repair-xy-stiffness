# -*- coding: utf-8 -*-
# Section 3.2: tests for bistability and hysteresis (exec'd by build_manuscript.py between 3.1 and 3.3)
BR = R["branches"]
BS = BR["summary"]
RP = R["ramps"]
CL0 = R["closure0"]
LAM = ST["lam0"]
WC_ = ST["windows"]

HD("3.2 Tests for bistability and hysteresis", 2)
qd = max(abs(v["q_diff"]) for k, v in BR.items() if k != "summary")
# closure-model branch separation at rho = 6 for beta = 20, in damaged fraction
qq = np.linspace(1e-4, 1 - 1e-4, 200000)
F6 = MF.rho_ode(qq, 20.0) - 6.0
roots = qq[np.where(np.sign(F6[:-1]) != np.sign(F6[1:]))[0]]
exc_txt = "; ".join("*L* = %s, β = %s, ρ = %s: %+.4f [%+.4f, %+.4f]" % (e["key"].split("|")[0], e["key"].split("|")[1], e["key"].split("|")[2], e["diff"], e["ci95"][0], e["ci95"][1])
                    for e in BS["exceptions"])
nsig = sum(1 for e in BS["exceptions"] if e["holm_p"] < 0.05)
P("Section 2.3 predicts, for β = 20, a bistable window %s < ρ < %s in the map the code implements (%s < ρ < %s in continuous time), and %s < ρ < %s "
  "for β = 5. We tested this on the lattice with the original repair rule and *h*_{0} = 0.4 (Figure 2). From an intact start and from a start with "
  "95%% of the sites damaged and random phases, at 11 values of ρ from 2 to 12, for *L* = 32 and 48 (β = 5 and 20) and *L* = 64 (β = 20), "
  "with 6–10 seeds per condition and 3000 + 3000 steps, the two starts reach the same state everywhere. Over the %d conditions the largest "
  "difference in *U* is %s and in the damaged fraction %s, whereas the closure model, equation (5), has stable states at damaged fractions "
  "differing by %s at ρ = 6. Nine of the %d bootstrap intervals of the difference exclude zero; they are (difference [95%% interval]) %s. "
  "All are at most %s in absolute value, they have no common sign, and none is significant after correction for the %d comparisons by Holm’s "
  "procedure (smallest adjusted *p* = %s). With only 6–10 runs per start percentile bootstrap intervals are somewhat too narrow, so a few "
  "exclusions among %d intervals are expected; they are not evidence of a residual dependence on the initial state."
  % (f3(WC_["20"]["discrete"][0], 2), f3(WC_["20"]["discrete"][1], 2), f3(WC_["20"]["continuous"][0], 2), f3(WC_["20"]["continuous"][1], 2),
     f3(WC_["5"]["discrete"][0], 2), f3(WC_["5"]["discrete"][1], 2), BS["n_tests"], f3(BS["max_abs_diff"], 4), f3(qd, 4),
     f3(roots[-1] - roots[0], 2) if len(roots) >= 3 else "—", BS["n_tests"], exc_txt, f3(BS["max_abs_diff"], 4), BS["n_tests"],
     f3(BS["min_holm_p"], 2), BS["n_tests"]))
dw = sorted(RP["dwells"], key=float)
fast, slow = RP["dwells"][dw[0]], RP["dwells"][dw[-1]]
P("The second test is a quasi-static ramp. At *L* = 48 the ratio ρ is stepped from 16 down to 1 and back in 31 steps, with a dwell of %s to %s "
  "steps per value (8 runs per rate). A hysteresis loop is present at fast sweeps, with signed area ∫(*U*_{↓} − *U*_{↑})dρ = %s ± %s at "
  "%s steps per value, but it closes as the sweep slows: %s ± %s at %s steps, and over the slowest dwells (%s steps and more, 24 runs) the area "
  "is %s (95%% bootstrap interval %s to %s), consistent with zero (Figure 2c,d). The loop at finite sweep rate is therefore the lag of a "
  "single-valued steady state behind a moving control parameter, not metastability; a bistable system keeps a finite loop as the sweep "
  "slows [[rao1990]]. We state the result as it is: within the tested ranges (*L* ≤ 64, β = 5 and 20, runs of 3000 steps, up to %s steps per "
  "ramp value), **no hysteresis or dependence on the initial state is resolved**. This does not prove that none exists: stronger feedback "
  "(β > 30), larger lattices or much longer runs were not explored."
  % (dw[0], dw[-1], f3(fast["area"][0], 2), f3(fast["area"][1], 2), dw[0], f3(slow["area"][0], 3), f3(slow["area"][1], 3), dw[-1], RP["plateau"]["dwell_min"],
     f3(RP["plateau"]["area"], 3), f3(RP["plateau"]["area_ci95"][0], 3), f3(RP["plateau"]["area_ci95"][1], 3), dw[-1]))
FIG("Fig2_hysteresis_tests.png",
    "Tests for bistability and hysteresis (*h*_{0} = 0.4, original repair rule). (a) Stationary *U* from an intact (filled) and a 95%%-damaged "
    "(open) start for β = 20 and *L* = 32, 48, 64; shaded: the bistable windows of the closure model in continuous time (light) and in the "
    "implemented discrete-time map (dark). (b) Difference between the two starts with 95%% bootstrap intervals, β = 20 (filled) and 5 (open). "
    "(c) Quasi-static ramps of ρ down (solid) and up (dashed), *L* = 48, at 25, 200 and 1600 steps per value. (d) Signed loop area against "
    "dwell, 95%% CI over 8 runs.")
# ---- why: the closure, tested against the lattice, with an independent prediction from beta = 0 lattices
rows = sorted({float(k.split("|")[2]) for k in CL0 if k.startswith("0.4|20|")})
tab = [["ρ", "*q* (lattice)", "μ̄ applied by the lattice", "μ(*q*), closure", "μ̄ from β = 0 lattices at the same *q*"]]
q0 = np.array([CL0["0.4|0|%g" % r]["q"][0] for r in sorted({float(k.split("|")[2]) for k in CL0 if k.startswith("0.4|0|")})])
m0 = np.array([CL0["0.4|0|%g" % r]["mu_ref_kill"][0] for r in sorted({float(k.split("|")[2]) for k in CL0 if k.startswith("0.4|0|")})])
o0 = np.argsort(q0)
for r in (2.0, 3.0, 4.0, 5.0, 7.0, 10.0):
    if r not in rows:
        continue
    c = CL0["0.4|20|%g" % r]
    tab.append([f3(r, 1), f3(c["q"][0]), f3(c["mu_kill"][0], 2), f3(c["mu_closure_q"], 2), f3(float(np.interp(c["q"][0], q0[o0], m0[o0])), 2)])
idc = max(abs(CL0["0.4|20|%g" % r]["rho_identity"] / r - 1) for r in rows)
# reduced steady-state relation with the multiplier measured on beta = 0 lattices (damage independent of phase): unique root?
qgrid = np.linspace(max(q0.min(), 0.01), min(q0.max(), 0.98), 3000)
mu_q = np.interp(qgrid, q0[o0], m0[o0])
rho_pred = MF.rho_map(qgrid, 20.0, LAM, mubar=mu_q)
mono = bool(np.all(np.diff(rho_pred) < 0))
q_pred = {r: float(np.interp(r, rho_pred[::-1], qgrid[::-1])) for r in (3.0, 5.0)}
P("The reduction fails to describe the lattice because its closure fails, not because of the kinetics. Because the kill step is where the "
  "closure enters, we recorded in every kill step the mean multiplier μ̄ that the lattice actually applies; with it the map (3) is satisfied "
  "by the stationary lattice to within %s in ρ (a check of the map, equation (3), itself). @T:closure@ and Figure 1c show that μ̄ differs from "
  "the closure μ(*q*) in both directions: at ρ ≤ 4 the lattice kills far more slowly than the closure assumes, because the intact sites see a mean "
  "local disorder well below the damaged fraction, and at ρ ≥ 5 it kills faster, because the few intact sites in disordered pockets are killed at nearly "
  "the full rate even in a mostly coherent lattice. The effective feedback is flattened, and a flattened μ̄(*q*) does not bend the steady-state "
  "relation into an S shape. The test is not circular: lattices with state-independent damage (β = 0) give, at each damaged fraction, the "
  "multiplier that the β = 20 feedback would apply to them (last column of @T:closure@); inserted in the map (3), this multiplier gives a "
  "steady-state relation ρ(*q*) that is %s (a single steady state for every ρ) and predicts *q* = %s and %s at ρ = 3 and 5, against the "
  "measured %s and %s for β = 20 (the difference is the effect of the correlation between damage and disorder, which the β = 0 lattices lack)."
  % ("%.1f%%" % (100 * idc), "monotonic" if mono else "not monotonic", f3(q_pred[3.0]), f3(q_pred[5.0]), f3(CL0["0.4|20|3"]["q"][0]), f3(CL0["0.4|20|5"]["q"][0])))
TAB(tab, "Kill multiplier in the steady-state lattice (*h*_{0} = 0.4, β = 20, *L* = 32; repair θ → 0): the value μ̄ applied by the lattice in its kill "
         "step, the closure value μ(*q*) at the measured damaged fraction *q*, and the multiplier that a β = 20 feedback would apply to "
         "β = 0 lattices at the same *q*.", widths=[0.5, 1.0, 1.4, 1.2, 1.9], label="closure")
