# -*- coding: utf-8 -*-
# Abstract and keywords (exec'd by build_manuscript.py right after the title block)
_TH = R["thresholds"]
_SWk = R["stiffness"]
_NKv = 2 * 0.35 / math.pi


def _brk(mode):
    t = _TH[mode]
    lo, hi = t["rho_vanished_max"], t["rho_finite_min"]
    if lo is None:                                                # only a lower bound on the size scan: use the finite-size crossing instead
        return "ρ ≈ %s" % f3(hi, 2)
    return "ρ between %s and %s" % (f3(lo, 2), f3(hi, 2))


def _eq_under():
    r = [1 - _SWk["32|20|%s|%g|intact" % (m, rho)]["ratio_eq_tw"] for m in ("keep", "neighbor") for rho in (6.0, 7.0, 8.0, 10.0)
         if _SWk["32|20|%s|%g|intact" % (m, rho)]["ratio_eq_tw"] is not None and _SWk["32|20|%s|%g|intact" % (m, rho)]["ups_tw"][0] > 0.2]
    return int(5 * round(100 * max(r) / 5))


_SL = R["timescale"]["slow_limit"]
_rho0 = None
_ABS = (
    "Neutron damage lowers the transition temperature of REBCO conductors roughly linearly with displacement damage, and annealing recovers part "
    "of it, implicating phase stiffness. We model this competition as a two-dimensional XY lattice of the "
    "superconducting phase whose sites are damaged at a rate rising with local phase disorder and repaired stochastically, simulated by kinetic "
    "Monte Carlo. A mean-field reduction predicts bistability; the lattice shows no hysteresis or dependence "
    "on the initial state in the ranges tested (*L* ≤ 64), because the closure of the reduction fails. At zero field, symmetry-preserving repair restores "
    "spontaneous coherence and a twist-response stiffness (validated on clean and diluted lattices) above an apparent threshold in the repair-to-damage "
    "ratio, bracketed by finite-size scans (*L* ≤ 96) at %s (%s with neighbor-phase repair), against ρ ≈ 2.5 without feedback; the equilibrium formula "
    "underestimates the driven stiffness by up to %d%%. Numerically, ρ controls the state ever better as damage slows relative to phase relaxation. "
    "Finite-size analysis of random dilution estimates the stiffness loss at *f* ≈ %s ± 0.01; clustered damage tolerates more. An illustrative effective mapping gives about ten model "
    "sites per dpa. The model is uncalibrated."
    % (_brk("keep"), _brk("neighbor"), _eq_under(), f3(round(R["fss"]["0.35"]["extrap_all"]["f_inf"], 2), 2)))
_nw = len(re.sub(r"\*", "", _ABS).split())
assert _nw <= 200, "abstract is %d words" % _nw
HD("Abstract")
P(_ABS)
P("**Keywords:** superconductivity; kinetic Monte Carlo; XY model; phase stiffness; irradiation damage; annealing; REBCO; fusion magnets", align="left")
print("abstract words:", _nw)
