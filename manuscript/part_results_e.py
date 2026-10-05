# -*- coding: utf-8 -*-
# Section 3.3 (continued): finite-size dependence of the dynamic stiffness and the threshold (exec'd by build_manuscript.py)
TH_ = R["thresholds"]
S1_, S2_ = R["stiff_L"], R["stiff_L2"]


def yL(mode, L, rho):
    """Twist-response stiffness [mean, se] for (mode, L, rho) from the size scans (None if not measured)."""
    v = S2_.get("%s|%d|%g" % (mode, L, rho))
    if v is None and mode == "keep" and rho in (5.0, 7.0, 10.0, 20.0):
        v = S1_.get("%d|%g" % (L, rho))
    return v["ups_tw"] if v else None


rhos_k = sorted({float(k.split("|")[2]) for k in S2_ if k.startswith("keep|")} | {5.0, 7.0, 10.0, 20.0})
rhos_n = sorted({float(k.split("|")[2]) for k in S2_ if k.startswith("neighbor|")})
rows = [["repair rule", "ρ", "*L* = 24", "32", "48", "64", "96", "status"]]
for mode, rr in (("keep", rhos_k), ("neighbor", rhos_n)):
    for rho in rr:
        cells = [yL(mode, L, rho) for L in (24, 32, 48, 64, 96)]
        rows.append([mode if rho == rr[0] else "", f3(rho, 2)] + [pm(c, 3) if c else "–" for c in cells] +
                    [TH_[mode]["rows"].get("%g" % rho, {}).get("status", "–")])
bk = (TH_["keep"]["rho_vanished_max"], TH_["keep"]["rho_finite_min"])
bn = (TH_["neighbor"]["rho_vanished_max"], TH_["neighbor"]["rho_finite_min"])
y_k = TH_["keep"]["Y_at_lowest_finite_rho_Lmax"]
y_n = TH_["neighbor"]["Y_at_lowest_finite_rho_Lmax"]
SDm = R["stiffness_dmg"]
_dif = []
for k, v in SDm.items():
    o = SW[k.replace("damaged", "intact")]
    _dif.append((k, v["ups_tw"][0] - o["ups_tw"][0], v["ups_tw"][1]))
_nfl, _ntot = sum(v["n_flagged"] for v in SDm.values()), sum(v["n"] for v in SDm.values())
_big = max(_dif, key=lambda d: abs(d[1]))
P("The stiffness does not depend on the initial state either. A start with 95%% of the sites damaged and random phases can leave windings and "
  "vortices in the lattice, so for such states the twist response was measured with the winding number of every run recorded and the response "
  "corrected for a global integer winding when the two runs of a seed pair differ in it (%d of the %d seed pairs). In the %d conditions tested "
  "(*keep*, ρ = 8–20; *neighbor*, ρ = 7–20; *L* = 32, 10 seeds) the twist response from the damaged start agrees with that from the ordered start to "
  "within %s, the largest difference being %s ± %s at %s, ρ = %s." %
  (_nfl, _ntot, len(SDm), f3(abs(_big[1]), 3), "%+.3f" % _big[1], f3(_big[2], 3), _big[0].split("|")[2], _big[0].split("|")[3]))
P("Is the threshold a finite-size effect? Deep in the coherent state the stiffness does not depend on the size of the lattice (@T:thr@): at ρ = 10 "
  "and 20 the twist response of *keep* is %s–%s and %s–%s for *L* = 24 to 96, so Υ_{tw} is a bulk quantity. Near the threshold it depends on *L*. For *keep* at "
  "ρ = 7 it falls with size (%s at *L* = 24, %s at 96) and is zero within error at the largest size, whereas at ρ = 8 it is finite and nearly "
  "size-independent (%s at *L* = 32, %s at 96, a decrease of %s that has stopped by *L* = 64). We call the stiffness *vanished* at a "
  "ρ if the twist response at the largest size is below 0.05 (with at least three sizes) and *finite* if it exceeds 0.15 there and has not fallen below "
  "0.7 of its value at the smallest size; the finite-size estimate of the threshold is bracketed by the largest vanished and the smallest finite ρ. The cutoffs 0.05, 0.15 and "
  "0.7 are operational choices, not properties of the model, and another choice would move the brackets. For the sizes studied the estimate lies "
  "between ρ = %s and %s for *keep* (*L* up to 96) and between %s and %s for *neighbor* (*L* up to 64) (@T:thr@); these brackets are consistent with a "
  "threshold, and are not bounds on its value in the thermodynamic limit. At the lowest ρ at "
  "which the stiffness is finite its value is %s ± %s for *keep* (*L* = %s) and %s ± %s for *neighbor* (*L* = %s), larger than 2*T*/π = %s. "
  "For *keep* the approach to the threshold is visible in the size dependence (@T:thr@): at ρ = 7.25 and 7.5 the stiffness still falls with *L* "
  "(%s → %s and %s → %s from *L* = 32 to 96), while at ρ = 7.75 it is size-independent at %s ± %s (*L* = 96), so the finite-size estimate of the threshold lies in the upper "
  "part of the bracket, probably between 7.5 and 7.75 for *L* ≤ 96. A stiffness that decreases with size through 2*T*/π and settles above it just beyond the "
  "threshold is what a BKT-type continuous onset looks like; our data are compatible with it and do not exclude a discontinuity narrower than "
  "the step of 0.25 in ρ. No hysteresis between the two starts is resolved (Figure 3c)."
  % (f3(min(yL("keep", L, 10.0)[0] for L in (24, 32, 48, 64, 96)), 2), f3(max(yL("keep", L, 10.0)[0] for L in (24, 32, 48, 64, 96)), 2),
     f3(min(yL("keep", L, 20.0)[0] for L in (24, 32, 48, 64, 96)), 2), f3(max(yL("keep", L, 20.0)[0] for L in (24, 32, 48, 64, 96)), 2),
     f3(yL("keep", 24, 7.0)[0], 2), f3(yL("keep", 96, 7.0)[0], 3), f3(yL("keep", 32, 8.0)[0], 2), f3(yL("keep", 96, 8.0)[0], 2),
     "%d%%" % round(100 * (1 - yL("keep", 96, 8.0)[0] / yL("keep", 32, 8.0)[0])),
     f3(bk[0], 2), f3(bk[1], 2), f3(bn[0], 2) if bn[0] else "(not reached)", f3(bn[1], 2),
     f3(y_k[0], 2), f3(y_k[1], 2), TH_["keep"]["Lmax_at_lowest_finite"], f3(y_n[0], 2), f3(y_n[1], 2), TH_["neighbor"]["Lmax_at_lowest_finite"], f3(NK),
     f3(yL("keep", 32, 7.25)[0], 2), f3(yL("keep", 96, 7.25)[0], 2), f3(yL("keep", 32, 7.5)[0], 2), f3(yL("keep", 96, 7.5)[0], 2),
     f3(yL("keep", 96, 7.75)[0], 2), f3(yL("keep", 96, 7.75)[1], 2)))
TAB(rows, "Twist-response stiffness Υ_{tw} (mean ± s.e.) of the field-free steady state (β = 20, *T* = 0.35, ordered start, 6–10 seeds) against the lattice "
          "size *L* across the threshold, with the status assigned by the operational criterion in the text (not a determination of the thermodynamic limit).", widths=[0.9, 0.5, 0.95, 0.95, 0.95, 0.95, 0.95, 0.8], size=7.5, label="thr")
