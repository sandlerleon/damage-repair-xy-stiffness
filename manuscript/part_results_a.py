# -*- coding: utf-8 -*-
# Section 3.1 and 3.3 (exec'd by build_manuscript.py)
ST = R["steady"]
FF = R["fieldfree"]
SW = R["stiffness"]
NK = 2 * 0.35 / math.pi


def logcross(xs, ys, level, decreasing=False):
    """Crossing of y(x) with `level`, interpolated linearly in log x (NaN if none)."""
    for i in range(len(xs) - 1):
        a, b = ys[i] - level, ys[i + 1] - level
        if a * b <= 0 and ys[i] != ys[i + 1]:
            t = (level - ys[i]) / (ys[i + 1] - ys[i])
            return float(math.exp(math.log(xs[i]) + t * (math.log(xs[i + 1]) - math.log(xs[i]))))
    return float("nan")


def sw_series(mode, beta, key="ups_tw"):
    rs = [r for r in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 40.0] if SW.get("32|%g|%s|%g|intact" % (beta, mode, r))]
    return rs, [SW["32|%g|%s|%g|intact" % (beta, mode, r)][key][0] for r in rs]


def eta_series(mode, beta, L=32):
    rs = [r for r in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 40.0]
          if FF.get("%d|%g|%s|%g" % (L, beta, mode, r)) and FF["%d|%g|%s|%g" % (L, beta, mode, r)]["intact"]["eta_int"]]
    return rs, [FF["%d|%g|%s|%g" % (L, beta, mode, r)]["intact"]["eta_int"] for r in rs]


# ---------------------------------------------------------------------------------- 3.1
HD("3. Results")
HD("3.1 Steady-state phase coherence", 2)
rhos_all = sorted({float(k.split("|")[2]) for k in ST["long"]})
dU = max(abs(ST["long"]["32|%g|%g" % (b, r)]["U"][0] - ST["long"]["64|%g|%g" % (b, r)]["U"][0]) for b in (0, 5, 20) for r in rhos_all)
sv = ST["short_vs_long"]
dmax_short = max(abs(c["diff"]) for c in sv)
Jmax = max(v["J_rel_imbalance"] for v in ST["long"].values())
tau_lo = min(v["tau_U"] for v in ST["long"].values())
tau_hi = max(v["tau_U"] for v in ST["long"].values())
CV = R["converge"]
P("With the original repair rule (θ → 0) and a weak field (*h*_{0} = 0.4), the phase coherence *U* rises smoothly with ρ and saturates at "
  "%s for ρ = 40 (Figure 1a, @T:steady@); the saturation is below unity because of thermal phase fluctuations. State-dependent damage "
  "(β > 0) depresses *U* at small ρ, from %s at ρ = 1 for β = 0 to %s and %s for β = 5 and 20 (*L* = 64), and its effect has disappeared "
  "by ρ ≈ 10. The results are independent of lattice size: the largest difference between *L* = 32 and *L* = 64 over all 33 conditions is "
  "%s. The state is reached quickly: runs started intact and 95%% damaged enter the stationary band of *U* within %s and %s steps at ρ = 2.5 and 5 "
  "(late-time difference %s and %s), a short protocol of 200 equilibration and 200 measurement steps agrees with the long runs to within "
  "%s (15 conditions), the kill and repair currents balance to a relative %s, and the integrated autocorrelation time of *U* is %s–%s steps. "
  "The dashed line in Figure 1a is the intact fraction ρ/(1 + ρ) of the rate equation; it is shown for orientation only, since *U* "
  "is not an intact fraction (damaged sites remain weakly coupled and thermal fluctuations reduce *U* even for an intact lattice)."
  % (f3(ST["long"]["64|0|40"]["U"][0]), f3(ST["long"]["64|0|1"]["U"][0]), f3(ST["long"]["64|5|1"]["U"][0]), f3(ST["long"]["64|20|1"]["U"][0]),
     f3(dU), CV["2.5|t_enter_intact"], CV["2.5|t_enter_damaged"], f3(abs(CV["2.5|late_gap"])), f3(abs(CV["5|late_gap"])), f3(dmax_short),
     sci(Jmax), "%.0f" % tau_lo, "%.0f" % tau_hi))
FIG("Fig1_steady_state.png",
    "Steady-state order and the mean-field reduction (*h*_{0} = 0.4, repair θ → 0). (a) *U* against ρ for three feedback strengths β, "
    "*L* = 32 (filled) and 64 (open), mean ± 95%% CI over 10 and 6 seeds (the bars are smaller than the symbols); dashed: the intact fraction "
    "ρ/(1 + ρ) of the rate equation. (b) Damaged fraction *q* of the lattice (*L* = 64) and of the closure model, equation (5) (lines; dotted "
    "where unstable). (c) Mean kill multiplier of the intact sites in the lattice and in the closure μ(*q*).")
rows = [["β", "ρ = 1", "ρ = 2.5", "ρ = 5", "ρ = 10", "ρ = 40"]]
for b in (0, 5, 20):
    rows.append([str(b)] + [pm(ST["long"]["64|%d|%g" % (b, r)]["U"]) for r in (1.0, 2.5, 5.0, 10.0, 40.0)])
TAB(rows, "Steady-state coherence *U* (mean ± s.e., 6 seeds, *L* = 64, *h*_{0} = 0.4, repair θ → 0). The corresponding values for *L* = 32 "
          "(10 seeds) differ by at most %s." % f3(dU), widths=[0.5, 1.1, 1.1, 1.1, 1.1, 1.1], label="steady")

# ---------------------------------------------------------------------------------- 3.3
HD("3.3 Phase stiffness and spatial correlations without a field", 2)
rk, yk = sw_series("keep", 20)
rn, yn = sw_series("neighbor", 20)
r0, y0 = sw_series("keep", 0)
rc_k, rc_n, rc_0 = logcross(rk, yk, NK), logcross(rn, yn, NK), logcross(r0, y0, NK)
ek = eta_series("keep", 20)
en = eta_series("neighbor", 20)
e0 = eta_series("keep", 0)
re_k, re_n, re_0 = (logcross(*e, 0.25) for e in (ek, en, e0))
_gaps = [(k, v) for k, v in FF.items() if "gap_U" in v]
_exc = [(k, v) for k, v in _gaps if not (v["gap_ci95"][0] <= 0 <= v["gap_ci95"][1])]
_trapk = max(_gaps, key=lambda kv: abs(kv[1]["gap_U"]))[0]
_trapU = [x for x in FF[_trapk]["damaged"]["U_per_seed"] if x < 0.5]
_trapD = FF[_trapk]["damaged"]["U_per_seed"]
assert len(_trapU) == 1, "the trapped-run description assumes exactly one trapped run"
P("The field *h*_{0} aligns the phases by itself, and the rule θ → 0 selects a phase that is common to all sites, so for *h*_{0} = 0 it still orders "
  "the lattice at rate *R* (Section 2.2; the global phase is pinned to 0, |⟨exp(iψ)⟩| = 1.00 in the test suite). To ask whether repair restores "
  "*spontaneous* phase coherence and stiffness we therefore set *h*_{0} = 0 and use the two repair rules that preserve the U(1) symmetry, *keep* "
  "and *neighbor*, with *reset0* shown for comparison. Figure 3a shows *U*(ρ) from an intact and from a 95%%-damaged start. For both symmetric "
  "rules there is a narrow range of ρ in which the lattice goes from an incoherent state (*U* at the finite-size floor of about %s for random "
  "phases at *L* = 32) to a coherent one, and the position of this range is set by the feedback and by the repair rule: with state-independent damage "
  "(β = 0) the stiffness sets in at ρ = %s, with the baseline feedback (β = 20) at ρ = %s for *keep* and ρ = %s for *neighbor* (*L* = 32, Figure 4a; "
  "finite-size estimates of the thresholds are given below). The "
  "reason is the one built into the model: damage is faster where phase order is already disrupted, so a coherent region must be repaired "
  "faster than a disordered one is damaged, and a repaired site that rejoins the coherence around it (*neighbor*) helps more than one that is "
  "left to relax (*keep*). With *reset0* there is no transition: the resets themselves align the phases, and *U* rises smoothly (Figure 3a). "
  "The initial state makes almost no difference: of the %d conditions compared (rules, β, ρ, *L*), the 95%% interval of the difference of *U* "
  "between the two starts excludes zero in %d, and in %d of these the difference is at most %s (Figure 3c). The remaining one is "
  "state-independent damage (β = 0, *keep*) at ρ = 5, *L* = 32, where one of the ten runs started from random phases stayed incoherent for the whole run "
  "(*U* = %s, against %s–%s for the other nine and for all ten intact starts) at the same damaged fraction. We read it as a kinetic trap, a "
  "persistent configuration of phase defects inherited from the random initial phases, in a regime in which the coherent state is the stable one, "
  "and not as bistability; it is the only such run among the 640 two-start runs of this study, and the stiffness measured from damaged starts "
  "agrees with that from ordered starts where we tested it (below)." % (
      f3(1 / math.sqrt(1024) * math.sqrt(math.pi / 4), 2), f3(rc_0, 1), f3(rc_k, 1), f3(rc_n, 1), len(_gaps), len(_exc), len(_exc) - 1, f3(max(abs(v["gap_U"]) for k, v in _exc if k != _trapk), 3), f3(min(_trapU), 3),
      f3(min(x for x in _trapD if x > 0.5), 3), f3(max(_trapD), 3)))
FIG("Fig3_field_free.png",
    "Field-free steady states (*h*_{0} = 0) from two initial states. (a) *U* against ρ for the repair rules *keep* and *neighbor* (β = 20), "
    "*reset0* (β = 20) and *keep* with state-independent damage (β = 0); *L* = 32, 10 seeds; filled: intact start, open: 95%%-damaged start; "
    "triangles: *L* = 48. (b) Damaged fraction. (c) Difference of *U* between the two starts with 95%% bootstrap intervals.")
P("Coherence in magnitude is not the same as phase rigidity, so we measure the stiffness. Figure 4a shows the twist-response stiffness "
  "Υ_{tw}(ρ) of equation (6) in the ordered start state (*L* = 32; the twist φ_{0} = 0.3π/*L* = %s is well below π/*L*, beyond which a "
  "periodic lattice changes its winding number, and the integer winding of every run is recorded; none of the %d seed pairs differs in "
  "global winding). The stiffness is zero within error (|Υ_{tw}| ≤ %s) at ρ ≤ 6 for *keep* and ρ ≤ 5 for *neighbor*, and rises through "
  "the Nelson–Kosterlitz value 2*T*/π = %s, used only as a reference because the equilibrium criterion has not been justified here, at "
  "ρ = %s (*keep*), %s (*neighbor*) and %s (β = 0), at *L* = 32. Repair therefore restores genuine phase stiffness, with a threshold in ρ. "
  "The spatial correlations (Figure 5) give the same threshold: for ρ below it *C*(*r*) decays exponentially with a correlation length of "
  "a few lattice spacings, and above it *C*(*r*) decays algebraically, ~ *r*^{−η}, with η = %s ± 0.01 at ρ = 40 and η = 1/4 (the "
  "Berezinskii–Kosterlitz–Thouless value at the transition) crossed at ρ = %s (*keep*), %s (*neighbor*) and %s (β = 0). In the algebraic "
  "phase the exponent measured from *C*(*r*) agrees with the spin-wave relation η = *T*/(2πΥ_{tw}) (Figure 5c) to within %s over ρ ≥ 10, which "
  "supports Υ_{tw} as the stiffness of the damaged-and-repairing medium even though the state is not an equilibrium state."
  % (f3(R["stiffness"]["32|20|keep|40|intact"]["phi0"], 4), sum(v["n"] for k, v in SW.items() if k.endswith("|keep|40|intact")),
     f3(max(abs(y) for r, y in zip(rk, yk) if r <= 6), 3), f3(NK), f3(rc_k, 1), f3(rc_n, 1), f3(rc_0, 1),
     f3(FF["32|20|keep|40"]["intact"]["eta_int"], 2), f3(re_k, 1), f3(re_n, 1), f3(re_0, 1),
     "%d%%" % round(100 * max(abs(0.35 / (2 * math.pi * SW["32|20|%s|%g|intact" % (m, r)]["ups_tw"][0]) / FF["32|20|%s|%g" % (m, r)]["intact"]["eta_int"] - 1)
                              for m in ("keep", "neighbor") for r in (10.0, 12.0, 15.0, 20.0, 40.0)))))
FIG("Fig4_stiffness.png",
    "Phase stiffness. (a) Twist-response stiffness Υ_{tw}, equation (6), against ρ for the symmetric repair rules (*h*_{0} = 0, *L* = 32, ordered "
    "start, 10 seeds, 95%% CI); the dashed line is the Nelson–Kosterlitz value 2*T*/π, shown as a reference only. (b) The equilibrium fluctuation "
    "formula, equation (7), evaluated on the same runs, against Υ_{tw}; the line is equality. (c) Υ_{tw} against *U*.")
FIG("Fig5_correlations.png",
    "Phase correlations *C*(*r*) of intact-site pairs in the field-free steady state (*L* = 32, β = 20) for (a) *keep* and (b) *neighbor* at "
    "ρ from 3 to 20. (c) The exponent η of an algebraic fit (2 ≤ *r* ≤ *L*/4) and the spin-wave value *T*/(2πΥ_{tw}); dotted: η = 1/4.")
rat = {(m, r): SW["32|20|%s|%g|intact" % (m, r)]["ratio_eq_tw"] for m in ("keep", "neighbor") for r in (6.0, 7.0, 8.0, 10.0, 15.0, 20.0)
       if SW["32|20|%s|%g|intact" % (m, r)]["ratio_eq_tw"] is not None}
neg = min(SW["32|20|keep|%g|intact" % r]["ups_eq"][0] for r in (5.0, 6.0, 7.0))
P("The equilibrium fluctuation formula, equation (7), cannot be applied here without checking, and the check fails near the transition "
  "(Figure 4b). It underestimates the twist response by %s at ρ = 8 for *keep* (ratio %s), by %s at ρ = 6 for *neighbor* and by %s and %s at ρ = 10 "
  "and 15 for *keep*; far from the transition it recovers, to within %s at ρ = 40. In the incoherent regime, where the twist response vanishes, "
  "it is negative (down to %s), which no stiffness can be. The cause is that the formula attributes all fluctuations of the bond current "
  "to the phases at a frozen coupling landscape, whereas here the couplings change with time as sites are damaged and repaired; the "
  "dependence of the discrepancy on the ratio of the damage and phase-relaxation times is examined in Section 3.4. Consequently the "
  "Nelson–Kosterlitz criterion Υ = 2*T*/π, which is an equilibrium statement, is used in this paper only for the static quenched "
  "dilution of Section 3.5, and the dynamic results are expressed through Υ_{tw}, η and the transition ratio ρ." %
  ("%d%%" % round(100 * (1 - rat[("keep", 8.0)])), f3(rat[("keep", 8.0)], 2), "%d%%" % round(100 * (1 - rat[("neighbor", 6.0)])),
   "%d%%" % round(100 * (1 - rat[("keep", 10.0)])), "%d%%" % round(100 * (1 - rat[("keep", 15.0)])),
   "%d%%" % round(100 * (1 - SW["32|20|keep|40|intact"]["ratio_eq_tw"])) if SW["32|20|keep|40|intact"]["ratio_eq_tw"] < 0.995 else "1%",
   f3(neg, 2)))
