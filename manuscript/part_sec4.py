# -*- coding: utf-8 -*-
# Section 4: relationship to irradiated superconductors (exec'd by build_manuscript.py)
MP = R["mapping"]
AF = MP["anneal_flat"]
AG = MP["anneal"]
a13, g13 = AF["%g" % 1e13], AG["%g" % 1e13]

HD("4. Relationship to irradiated superconductors")
HD("4.1 Experimental evidence the model is compared with", 2)
P("The experiments that bear on the model are those on REBCO coated conductors irradiated with fast neutrons and then annealed. "
  "(i) The transition temperature decreases linearly with the calculated displacement damage, at the same rate for neutrons and for "
  "helium ions; for neutron-irradiated tapes the decrease is from about 90 K to 81 K and from about 94 K to 84 K at 3.7–3.8 mdpa, a "
  "fractional decrease of about %s [[adams2023]]. (ii) The critical current first rises with fluence, passes a temperature-dependent "
  "maximum and falls below its pristine value at fluences of the order of 10^{22} m^{−2} (at 50 K in most tapes), and the n-value decreases even where "
  "the critical current is enhanced [[prokopec2015,fischer2018,fischer2019,unterrainer2022]]; this initial rise is due to added pinning and is outside "
  "the present model. (iii) The degradation of *T*_{c} and of the critical current density are closely related, attributed to a loss "
  "of superfluid density, and annealing for 12 h after a 5 °C min^{−1} ramp recovers about 25%% of the *T*_{c} decrease at 150 °C and about "
  "60%% at 400 °C, from near-linear increases of *T*_{c} with annealing temperature in several tapes (slopes of about 1.8–3.3 K per 100 °C); in an oxygen-poor "
  "environment oxygen loss limits the anneal to about 220 °C, so the higher-temperature recovery refers to anneals in oxygen [[unterrainer2022]]. These numbers are quoted from the cited papers; the "
  "model is compared with them only in the limited sense below." % f3(MP["dTc_over_Tc"], 2))
HD("4.2 From a damaged fraction to displacement damage", 2)
P("The superfluid density, which the experiments identify as what the damage removes, is the physical counterpart of the stiffness Υ "
  "(for a film of thickness *d* the phase stiffness is proportional to *d*/λ_{L}^{2}, with λ_{L} the London depth [[beasley1979,tinkham2004]]). The "
  "relative stiffness Υ_{tw}(ρ)/Υ_{tw}(∞) of the damage–repair steady state is therefore the model’s counterpart of the relative superfluid density and "
  "fixes λ_{L}/λ_{L}(0) = [Υ(∞)/Υ(ρ)]^{1/2}, a quantity that is measurable by penetration-depth experiments but that we do not compare with data "
  "(@T:rel@). We do not equate *U*, or Υ, with a critical current: the model has neither pinning nor a transport current.")
rows = [["ρ", "*q* (keep)", "Υ_{tw}/Υ_{tw}(40)", "λ_{L}/λ_{L}(0)", "*q* (neighbor)", "Υ_{tw}/Υ_{tw}(40)", "λ_{L}/λ_{L}(0)"]]
for r in (8.0, 10.0, 12.0, 15.0, 20.0):
    cells = [f3(r, 0)]
    for mode in ("keep", "neighbor"):
        e = SW["32|20|%s|%g|intact" % (mode, r)]
        e0 = SW["32|20|%s|40|intact" % mode]
        rel = e["ups_tw"][0] / e0["ups_tw"][0]
        cells += [f3(e["q"][0], 3), f3(rel, 3), f3(1 / math.sqrt(rel), 3)]
    rows.append(cells)
TAB(rows, "Relative stiffness and the implied London-depth ratio in the field-free steady state (*L* = 32, β = 20, ordered start); the "
          "reference is the stiffness at ρ = 40, where the damaged fraction is 0.5%.", widths=[0.5, 0.9, 1.1, 1.0, 1.0, 1.1, 1.0], label="rel")
P("To relate the damaged fraction *f* of static damage to displacement damage we need the dependence of the transition temperature on *f*, and "
  "here, unlike for the coherence *U*, the model has a precise answer: the BKT temperature *T*_{BKT}(*f*) of the diluted lattice, from the crossing "
  "of Υ(*T*) with 2*T*/π (Section 2.4) extrapolated in 1/ln²*L*, falls from the clean value %s [[hasenbusch2005]] to %s (95%% interval %s–%s) at "
  "*f* = 0.05 and to %s at *f* = 0.10 (Figure 8a). If the *T*_{c} of the conductor is taken to follow *T*_{BKT}(*f*), a fractional decrease "
  "of %s at 3.7–3.8 mdpa [[adams2023]] corresponds to *f* = %s, that is to *a* = %s model sites made non-superconducting per dpa (range %s–%s "
  "across the interval estimates and the two tapes). This number is an effective conversion factor, not a measured or derived physical quantity: it depends on what a model site is taken to represent "
  "and on the assumption that *T*_{c} follows *T*_{BKT}(*f*), and we attach no further physical interpretation to it. "
  "The assumption that *T*_{c} follows *T*_{BKT}(*f*) of a two-dimensional XY model is an assumption and not a result; REBCO is a layered "
  "three-dimensional superconductor, and the proportionality of *T*_{c} to the phase stiffness is a relation that the "
  "experiments suggest [[unterrainer2022]] but that the model does not derive."
  % (f3(MP["T_BKT_clean"], 3), f3(R["tcurve"]["0.05"]["extrap"]["T_BKT"], 2), f3(R["tcurve"]["0.05"]["extrap"]["ci95"][0], 2),
     f3(R["tcurve"]["0.05"]["extrap"]["ci95"][1], 2), f3(R["tcurve"]["0.1"]["extrap"]["T_BKT"], 2), f3(MP["dTc_over_Tc"], 2),
     f3(MP["f_at_dpa"], 3), f3(MP["sites_per_dpa"], 0), f3(MP["sites_per_dpa_range"][0], 0), f3(MP["sites_per_dpa_range"][1], 0)))
P("**Rates, probabilities and time.** The model’s probabilities are per kinetic step, and the physical duration of a step enters the mapping. With "
  "*ṅ* the rate at which a given site is made non-superconducting by displacement damage, *ṅ* = *a* × (dpa rate), and *ν*_{a} exp(−*E*/*k*_{B}*T*_{a}) "
  "the thermally activated recovery rate with attempt frequency *ν*_{a} and activation energy *E*, the per-step probabilities for a step of "
  "duration Δ*t* are λ_{0} = 1 − exp(−*ṅ*Δ*t*) and *R* = 1 − exp[−*ν*_{a}exp(−*E*/*k*_{B}*T*_{a})Δ*t*]; a rate cannot be used as a probability. Their ratio "
  "ρ ≈ *ν*_{a}exp(−*E*/*k*_{B}*T*_{a})/*ṅ* is independent of Δ*t* only when both probabilities are small, which "
  "is the condition under which the rate equation (4) holds; in the simulations the hazard-equivalent ratio ln(1 − *R*)/ln(1 − λ_{0}) exceeds "
  "*R*/λ_{0} by a factor %s at λ_{0} = 0.04 and %s at λ_{0} = 0.08 (Section 3.4). Accumulated damage (dpa) and damage rate (dpa per second) "
  "enter differently: with no recovery the damaged fraction grows as *f* = *a*·dpa and the static dilution results of Section 3.5 apply; "
  "with simultaneous recovery the steady state is set by ρ."
  % (f3(R["timescale"]["collapse"]["keep"]["rho_hazard_over_rho_at_lam0.04"], 2), f3(R["timescale"]["collapse"]["keep"]["rho_hazard_over_rho_at_lam0.08"], 2)))
HD("4.3 Annealing", 2)
P("The rule θ → 0, the instant restoration of coupling and a single repair probability for all damaged sites are not physical descriptions of annealing. "
  "The first is replaced here by the *keep* and *neighbor* rules (Section 3.3), and we turn to the recovery kinetics. Consider, illustratively, 12 h isothermal "
  "anneals [[unterrainer2022]] of a tape damaged to the state above (*f* = %s), with a distribution of activation energies *g*(*E*) among "
  "the damaged sites, so that the fraction of sites still damaged after the anneal is ∫*g*(*E*)exp[−*ν*_{a}*t* exp(−*E*/*k*_{B}*T*_{a})]d*E* and the recovered fraction of the *T*_{c} decrease follows from "
  "*T*_{BKT}(*f*). A single activation energy gives a recovery that rises from 0 to 100%% over a few tens of kelvin, which is not what is "
  "reported (about 25%% at 150 °C and 60%% at 400 °C, the latter in oxygen, increasing roughly linearly with the annealing temperature). A Gaussian spectrum "
  "reproduces both reported values for *ν*_{a} = 10^{13} s^{−1} with a mean of %s eV and a standard deviation of %s eV, and a flat spectrum "
  "does so with a width of %s eV (from %s to %s eV); the fit of two parameters to two numbers is not a validation and does not determine the width uniquely (the authors of the experiment infer only "
  "qualitatively a very broad distribution [[unterrainer2022]]), and it gives a testable "
  "consequence only through its shape: a recovery of %s at 275 °C for both, and a recovery of %s already at 25 °C for the flat spectrum. Varying "
  "*ν*_{a} from 10^{11} to 10^{15} s^{−1} shifts the energies by about %s eV and leaves the shape unchanged. What survives any such choice is only the "
  "qualitative statement that recovery growing linearly over 150–400 °C in a 12 h anneal requires an activation spectrum much broader than a single process "
  "(the model fits give widths of the order of an electronvolt, which should be read as illustrative values of these two assumed forms), if the process is the independent thermally activated repair of independent damaged sites assumed here."
  % (f3(MP["f_at_dpa"], 3), f3(g13["E0_eV"], 2), f3(g13["sigma_eV"], 2), f3(a13["width_eV"], 1), f3(a13["Emin_eV"], 2), f3(a13["Emax_eV"], 2),
     "%d%%" % round(100 * g13["recovery_at_275C"]), "%d%%" % round(100 * a13["recovery_25C"]),
     f3(AG["%g" % 1e15]["E0_eV"] - AG["%g" % 1e11]["E0_eV"], 1)))
FIG("Fig8_mapping_annealing.png",
    "From damage to experiment. (a) *T*_{c}/*T*_{c0} = *T*_{BKT}(*f*)/*T*_{BKT}(0) from the helicity-modulus crossings (points, L = 16–48 extrapolated in 1/ln²*L*; "
    "*T*_{BKT}(0) from [[hasenbusch2005]]) and the interpolation used; the square marks the reported fractional decrease of *T*_{c} at 3.7–3.8 mdpa [[adams2023]], "
    "which fixes the number of model sites per dpa on the upper axis. (b) Recovered fraction of the *T*_{c} decrease after 12 h anneals for a Gaussian "
    "spectrum of activation energies, fitted to the two reported recoveries [[unterrainer2022]] (points) for three attempt frequencies.")
HD("4.4 What a calibration and a validation would require", 2)
P("A defensible calibration needs more than two annealing numbers. It requires (i) damage and anneal data on the *same* samples for *T*_{c}, "
  "for a quantity proportional to the stiffness (the penetration depth) and for the critical current, at several fluences; (ii) isochronal, not only "
  "isothermal, anneals, so that a spectrum of activation energies can be determined rather than assumed; (iii) a held-out set, for example a different "
  "fluence or a different tape, on which the mapping of Section 4.2 and the spectrum of Section 4.3 are tested; (iv) an extension of the model by pinning, "
  "to address the critical current, and by interlayer coupling, since the *T*_{c} of the conductor is not that of a single plane; and (v) spatially "
  "correlated damage. We have not attempted any of these, and the model remains uncalibrated.")
