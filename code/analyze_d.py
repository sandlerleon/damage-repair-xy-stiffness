# -*- coding: utf-8 -*-
"""Fourth part of the analysis: from damage fraction to displacement damage, and from annealing temperature to recovery.

The experimental numbers used here are those reported in refs [Adams et al 2023] and [Unterrainer et al 2022] (see the
manuscript): neutron-irradiated REBCO tapes lose about 10% of Tc (90 -> 81 K and 94 -> 84 K) at 3.7-3.8 mdpa; annealing for
12 h recovers about 25% of the Tc decrease at 150 C and about 60% at 400 C.

Assumptions (all explicit in the manuscript): (i) Tc is proportional to the BKT temperature of the model, Tc/Tc0 =
T_BKT(f)/T_BKT(0); (ii) the damage is random point damage (the clustered-damage study shows how much that matters);
(iii) the activation energies of the recovery processes are Gaussian-distributed, with attempt frequency nu.
"""
import json
import math

import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

import analyze as A

KB_EV = 8.617333262e-5
EXP = {"dTc_over_Tc": [9.0 / 90.0, 10.0 / 94.0], "mdpa": [3.7, 3.8], "recovery": {150.0: 0.25, 400.0: 0.60}, "anneal_h": 12.0}


T_BKT_CLEAN = 0.8929           # Hasenbusch 2005, high-precision value for the clean square-lattice XY model


def tau_points(res_tc, res_fss, which="mid"):
    """(f, T_BKT) points: the extrapolated crossings of the T scan (which = mid, lo or hi of the 95% interval), the clean
    value from the literature, and the extrapolated f_KT(T) of the fixed-T analysis."""
    pts = {0.0: T_BKT_CLEAN}
    for k, v in res_tc.items():
        if k == "T_grid" or float(k) == 0.0 or not v.get("extrap"):
            continue
        e = v["extrap"]
        pts[float(k)] = {"mid": e["T_BKT"], "lo": e["ci95"][0] if e.get("ci95") else e["T_BKT"], "hi": e["ci95"][1] if e.get("ci95") else e["T_BKT"]}[which]
    for T, v in res_fss.items():
        e = v.get("extrap_all")
        if e and float(T) < 0.4:
            pts[float(e["f_inf"])] = float(T)
    fs = sorted(pts)
    Ts = [pts[f] for f in fs]
    for i in range(1, len(Ts)):                    # enforce a non-increasing T_BKT(f) for the interpolation
        Ts[i] = min(Ts[i], Ts[i - 1])
    return fs, Ts


def sec_mapping(out):
    R = A.load_out()
    fs, Ts = tau_points(R["tcurve"], R.get("fss", {}), "mid")
    T0 = T_BKT_CLEAN
    tau = PchipInterpolator(fs, np.array(Ts) / T0)
    fgrid = np.linspace(fs[0], fs[-1], 4000)
    tg = tau(fgrid)

    def inv(curve, x):
        return float(fgrid[int(np.argmin(np.abs(curve - x)))])
    taus = {}
    for w in ("lo", "mid", "hi"):
        f_, T_ = tau_points(R["tcurve"], R.get("fss", {}), w)
        taus[w] = PchipInterpolator(f_, np.array(T_) / T0)(fgrid)
    small = [f for f in fs if 0 < f <= 0.1]
    slope = float(np.polyfit([0] + small, [1.0] + [float(tau(f)) for f in small], 1)[0])
    dec_list = EXP["dTc_over_Tc"]
    dec = float(np.mean(dec_list))
    dpa = float(np.mean(EXP["mdpa"])) * 1e-3
    f_star = inv(tg, 1.0 - dec)
    sites_per_dpa = f_star / dpa
    rng = []
    for w in ("lo", "mid", "hi"):
        for x in dec_list:
            for m in EXP["mdpa"]:
                rng.append(inv(taus[w], 1.0 - x) / (m * 1e-3))
    res = {"T_BKT_clean": T0, "T_BKT_clean_source": "Hasenbusch 2005", "points_f": fs, "points_T": Ts, "slope_dtau_df": slope,
           "dTc_over_Tc": dec, "dpa": dpa, "f_at_dpa": f_star, "sites_per_dpa": sites_per_dpa,
           "sites_per_dpa_range": [float(min(rng)), float(max(rng))]}
    # annealing: Gaussian distribution of activation energies, attempt frequency nu, hold time t
    t_s = EXP["anneal_h"] * 3600.0
    f0 = f_star

    def recovered(Ta_C, E0, sig, nu):
        Tk = Ta_C + 273.15
        E = np.linspace(max(E0 - 6 * sig, 0.01), E0 + 6 * sig, 2000)
        g = np.exp(-0.5 * ((E - E0) / sig) ** 2)
        g /= trapezoid(g, E)
        surv = trapezoid(g * np.exp(-nu * t_s * np.exp(-E / (KB_EV * Tk))), E)
        f_after = f0 * surv
        return float((tau(min(f_after, fs[-1])) - tau(f0)) / (1.0 - tau(f0)))

    out_an = {}
    for nu in (1e11, 1e13, 1e15):
        def resid(x):
            return [recovered(150.0, x[0], x[1], nu) - 0.25, recovered(400.0, x[0], x[1], nu) - 0.60]
        sol = least_squares(resid, [1.6 + 0.05 * math.log10(nu / 1e13), 0.3], bounds=([0.5, 0.02], [4.0, 1.5]))
        T_C = np.linspace(25, 500, 96)
        curve = [recovered(t, sol.x[0], sol.x[1], nu) for t in T_C]
        lin = np.polyfit([150.0, 400.0], [0.25, 0.60], 1)
        mid = np.arange(150.0, 401.0, 25.0)
        dev = [recovered(t, sol.x[0], sol.x[1], nu) - np.polyval(lin, t) for t in mid]
        out_an["%g" % nu] = {"E0_eV": float(sol.x[0]), "sigma_eV": float(sol.x[1]), "cost": float(sol.cost),
                             "T_C": T_C.tolist(), "recovery": curve, "max_dev_from_line_150_400": float(np.max(np.abs(dev))),
                             "recovery_at_275C": recovered(275.0, sol.x[0], sol.x[1], nu)}
    res["anneal"] = out_an
    # flat spectrum of activation energies between Emin and Emax: a site with E below E*(T) = kT ln(nu t) has annealed
    # (the step approximation to exp(-nu t exp(-E/kT)), exact to within a few kT); recovery of Tc follows from tau(f)
    flat = {}
    for nu in (1e11, 1e13, 1e15):
        def estar(Ta_C):
            return KB_EV * (Ta_C + 273.15) * math.log(nu * t_s)

        def rec_flat(Ta_C, Emin, Emax):
            phi = min(max((estar(Ta_C) - Emin) / (Emax - Emin), 0.0), 1.0)
            return float((tau(f0 * (1 - phi)) - tau(f0)) / (1.0 - tau(f0)))
        sol = least_squares(lambda x: [rec_flat(150.0, x[0], x[1]) - 0.25, rec_flat(400.0, x[0], x[1]) - 0.60], [1.0, 3.0],
                            bounds=([0.0, 0.5], [5.0, 8.0]))
        flat["%g" % nu] = {"Emin_eV": float(sol.x[0]), "Emax_eV": float(sol.x[1]), "width_eV": float(sol.x[1] - sol.x[0]),
                           "Estar_150C": estar(150.0), "Estar_400C": estar(400.0),
                           "recovery_25C": rec_flat(25.0, *sol.x), "recovery_275C": rec_flat(275.0, *sol.x)}
    res["anneal_flat"] = flat
    res["anneal_assumptions"] = {"hold_h": EXP["anneal_h"], "recovery_data": {"150C": 0.25, "400C": 0.60}}
    # tau(f) on a grid for the figure
    res["tau_f"] = {"f": fgrid[::40].tolist(), "tau": tg[::40].tolist()}
    out["mapping"] = res


A.SECTIONS["mapping"] = sec_mapping

if __name__ == "__main__":
    import sys
    import analyze_b  # noqa: F401
    import analyze_c  # noqa: F401
    out = A.load_out()
    for s in sys.argv[1:]:
        A.SECTIONS[s](out)
        print("analysed", s)
    A.save_out(out)
