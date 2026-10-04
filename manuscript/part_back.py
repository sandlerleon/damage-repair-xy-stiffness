# -*- coding: utf-8 -*-
# Sections 5 and 6, declarations and references (exec'd by build_manuscript.py)
_TK, _TN = R["thresholds"]["keep"], R["thresholds"]["neighbor"]
_FSS = R["fss"]["0.35"]["extrap_all"]
HD("5. Discussion and limitations")
P("**What the results say.** The mean-field reduction is an exact description of the kinetics only when its closure holds, and for the lattice it does "
  "not: the multiplier that the lattice applies to the damage rate is a smooth function of the damaged fraction, nearly the same whether the "
  "damage depends on phase order (β = 20) or not (β = 0), and far from the sharp sigmoid of the closure (Figure 1c). The bistability that the reduction "
  "predicts is therefore not found: two starting states converge, ramps close as they slow, and the interval exceptions are small, of mixed sign "
  "and not significant after correction. What survives in the lattice is a threshold, not a bistability. With the field removed and a repair that "
  "preserves the phase symmetry, the lattice has a threshold in ρ for spontaneous coherence, with phase correlations that change from exponential to "
  "algebraic and a stiffness, measured as a twist response, that is zero below it and large above it; state-dependent damage raises the threshold by a "
  "factor of about three over state-independent damage, because disordered regions are damaged faster, and a repair that restores the phase of the "
  "surroundings lowers it by a quarter. The ratio ρ is a sufficient control parameter in the limit in which the damage is slow compared with the "
  "relaxation of the phase, which we expect to be the limit of irradiation and annealing experiments, and in that limit the repair rule becomes irrelevant.")
P("**What it does not say.** The equilibrium Nelson–Kosterlitz criterion and the equilibrium helicity formula are not applied to the dynamic states, "
  "because the formula fails there (Section 3.3); the static dilution is the only place where they are used. The threshold in the thermodynamic limit is "
  "bracketed, not determined: for *keep* between ρ = %s and %s, for *neighbor* between %s and %s (@T:thr@), from sizes up to *L* = 96 and 64. "
  "Whether the stiffness jumps at the threshold, as it does for a BKT transition, or rises continuously over an interval narrower than our grid is "
  "not resolved. The absence of hysteresis refers to the sizes, run lengths and feedback strengths tested. No entropy-production result is claimed: "
  "the entropy production of a closed-form two-state cycle is a property of a framework, not a measurement of the simulated dynamics." % (f3(_TK["rho_vanished_max"], 2) if _TK["rho_vanished_max"] else "—", f3(_TK["rho_finite_min"], 2),
                                                    f3(_TN["rho_vanished_max"], 2) if _TN["rho_vanished_max"] else "—", f3(_TN["rho_finite_min"], 2)))
P("**Limitations of the model.** The lattice is two-dimensional and the damage is removal of sites, whereas displacement cascades produce extended "
  "defects, which matter (Section 3.5), some of which improve pinning; repair restores the coupling of a site at once and independently of its "
  "neighbors; there are no vortices as dynamical objects with pinning, no applied field or transport current, so the critical current and "
  "its initial rise with fluence are out of reach; the layered three-dimensional structure of REBCO is absent, so that *T*_{c} is not the BKT "
  "temperature of a plane; the activation spectrum of the annealing is assumed. Lattices are at most 96 × 96 (128 × 128 near the crossing of the static "
  "analysis), runs are at most a few thousand steps (1600 per value in the ramps), and 6–10 seeds or realizations were used. The mapping to displacement "
  "damage and annealing is illustrative and is not a calibration.")
HD("6. Conclusions")
P("A minimal damage–repair XY model gives three results that bear on the phase stiffness of irradiated superconductors, each stated with its limits. "
  "First, the mean-field bistability that the reduction predicts is not realized on the lattice in the ranges tested, because the closure of the reduction "
  "fails, and no hysteresis or dependence on the initial state is resolved. Second, repair restores spontaneous phase coherence and a "
  "twist-response stiffness above a threshold ratio of repair to damage that depends on the feedback and on the repair rule, the equilibrium "
  "fluctuation formula underestimates that stiffness by up to about %d%% near the threshold, and the stiffness, not the magnitude *U* of the order, is the "
  "quantity that distinguishes the coherent from the incoherent state. Third, the repair-to-damage ratio is a sufficient control parameter only in the limit of "
  "slow damage relative to phase relaxation, in which the repair rule no longer matters, and the finite-rate values differ from the slow-damage values by "
  "up to %s in *U*. The static dilution crossing, *f*_{KT} = %s ± 0.01 at *T* = 0.35, and the clustered-damage results complete the picture; "
  "an illustrative mapping to displacement damage and annealing gives about ten model sites per dpa and an activation spectrum of about 1 eV or more in "
  "width, and a calibration, which would need data that we did not have, remains to be done." %
  (_eq_under(), f3(R["timescale"]["slow_limit"]["5"]["reset0"]["U_base"][0] - R["timescale"]["slow_limit"]["5"]["reset0"]["U_slow"][0], 2),
   f3(round(_FSS["f_inf"], 2), 2)))
HD("Declarations")
P("**Funding.** This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors. "
  "**Competing interests.** The author declares no competing interests. **Author contributions (CRediT).** Leon Sandler: conceptualization, "
  "methodology, software, formal analysis, investigation, visualization, writing – original draft, writing – review and editing. "
  "**Declaration of AI-assisted technology.** During the preparation of this work the author used Claude (Anthropic) to write and test the simulation "
  "and analysis code to the author’s specification, to check the reference list against Crossref records and to assist with language and typesetting. "
  "The author reviewed and edited the content and takes full responsibility for it.")
P("**Data and code availability.** The engine, its test suite, the study and analysis scripts, all raw run outputs, the figures and the manuscript build "
  "script are available at %s, archived at https://doi.org/%s. This manuscript is available as a preprint at https://doi.org/%s. The engine reproduces "
  "the update rules of the MyUncle framework [[sandler_myuncle]] bit for bit. Experimental values quoted in Section 4 are taken from the cited papers "
  "and no experimental data were measured or fitted other than the numbers stated there." % (REPO, SW_DOI, PP_DOI))
HD("References")
for k in CITE:
    q = doc.add_paragraph()
    H.add_rich(q, "[%d] %s" % (CITE.index(k) + 1, RF.entry(k, DB)), size=10)
    q.paragraph_format.line_spacing = 1.15
    q.paragraph_format.space_after = Pt(3)
