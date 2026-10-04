# -*- coding: utf-8 -*-
"""Second part of the analysis: bistability tests, field-free stiffness, timescale competition, sensitivity."""
import math

import numpy as np
from scipy.stats import ttest_ind

import analyze as A
import meanfield as MF
from analyze import B, RNG, boot_ci, group, holm, load, mean_se, rho_of


# ================================================================ two-start tests (h0 = 0.4, original repair)
def sec_branches(out):
    D = load("branches")["runs"]
    keys = sorted({(r["params"]["L"], r["params"]["beta"], rho_of(r)) for r in D})
    res, tests = {}, []
    for (L, beta, rho) in keys:
        rr = [r for r in D if r["params"]["L"] == L and r["params"]["beta"] == beta and rho_of(r) == rho]
        a = np.array([r["U"]["mean"] for r in rr if r["kw"]["start"] == "intact"])
        b = np.array([r["U"]["mean"] for r in rr if r["kw"]["start"] == "damaged"])
        qa = np.array([r["q"]["mean"] for r in rr if r["kw"]["start"] == "intact"])
        qb = np.array([r["q"]["mean"] for r in rr if r["kw"]["start"] == "damaged"])
        diff = float(a.mean() - b.mean())
        boots = [a[RNG.integers(0, len(a), len(a))].mean() - b[RNG.integers(0, len(b), len(b))].mean() for _ in range(B)]
        p = float(ttest_ind(a, b, equal_var=False).pvalue)
        res["%d|%g|%g" % (L, beta, rho)] = {
            "n": [len(a), len(b)], "U_intact": mean_se(a), "U_damaged": mean_se(b), "diff": diff,
            "diff_ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            "q_diff": float(qa.mean() - qb.mean()), "welch_p": p, "seed_sd": float(np.concatenate([a, b]).std(ddof=1))}
        tests.append(("%d|%g|%g" % (L, beta, rho), p))
    hp = holm([p for _, p in tests])
    exc = []
    for (k, p), h in zip(tests, hp):
        res[k]["holm_p"] = h
        lo, hi = res[k]["diff_ci95"]
        if not (lo <= 0 <= hi):
            exc.append({"key": k, "diff": res[k]["diff"], "ci95": res[k]["diff_ci95"], "welch_p": p, "holm_p": h})
    res["summary"] = {"n_tests": len(tests), "max_abs_diff": max(abs(v["diff"]) for k, v in res.items() if k != "summary" and isinstance(v, dict) and "diff" in v),
                      "n_ci_exclude_zero": len(exc), "exceptions": exc, "min_holm_p": min(hp),
                      "max_seed_sd": max(v["seed_sd"] for k, v in res.items() if k != "summary" and isinstance(v, dict) and "seed_sd" in v)}
    out["branches"] = res


# ================================================================ quasi-static ramps
def sec_ramps(out):
    D = load("ramps")["runs"]
    res = {"dwells": {}}
    for dwell in sorted({r["dwell"] for r in D}):
        rr = [r for r in D if r["dwell"] == dwell]
        rho = np.array(rr[0]["rho_grid"])
        o = np.argsort(rho)
        areas = np.array([np.trapezoid((np.array(r["down"]) - np.array(r["up"]))[o], rho[o]) for r in rr])
        dn = np.mean([np.array(r["down"])[o] for r in rr], axis=0)
        up = np.mean([np.array(r["up"])[o] for r in rr], axis=0)
        res["dwells"][str(dwell)] = {"n": len(rr), "area": mean_se(areas), "area_ci95": boot_ci(areas), "down": dn.tolist(),
                                     "up": up.tolist(), "rho": rho[o].tolist(), "max_abs_gap": float(np.max(np.abs(dn - up)))}
    slow = [r for r in D if r["dwell"] >= 400]
    rho = np.array(D[0]["rho_grid"])
    o = np.argsort(rho)
    pl = np.array([np.trapezoid((np.array(r["down"]) - np.array(r["up"]))[o], rho[o]) for r in slow])
    res["plateau"] = {"dwell_min": 400, "n": len(pl), "area": float(pl.mean()), "area_ci95": boot_ci(pl)}
    out["ramps"] = res


# ================================================================ helicity-modulus estimators for the dynamic runs
def upsilon_eq(st, N, T):
    """Equilibrium fluctuation formula, averaged over the x and y bonds (valid in equilibrium only)."""
    return 0.5 * ((st["c_x"] - N / T * st["var_s_x"]) + (st["c_y"] - N / T * st["var_s_y"]))


def eta_from_C(C, rmin=2, rmax=None):
    """Exponent of an algebraic fit C(r) ~ r^(-eta) over rmin <= r <= rmax (or None if C is not positive there)."""
    C = np.asarray(C, float)
    r = np.arange(1, len(C) + 1)
    rmax = rmax or len(C)
    m = (r >= rmin) & (r <= rmax) & (C > 1e-3)
    if m.sum() < 4:
        return None
    p = np.polyfit(np.log(r[m]), np.log(C[m]), 1)
    return float(-p[0])


def sec_fieldfree(out):
    D = load("fieldfree")["runs"]
    res = {}
    keys = sorted({(r["params"]["L"], r["params"]["beta"], r["params"]["repair"], rho_of(r)) for r in D})
    for (L, beta, mode, rho) in keys:
        N = L * L
        T = D[0]["params"]["T"]
        ent = {}
        for start in ("intact", "damaged"):
            rr = [r for r in D if r["params"]["L"] == L and r["params"]["beta"] == beta and r["params"]["repair"] == mode
                  and rho_of(r) == rho and r["kw"]["start"] == start]
            if not rr:
                continue
            ups = [upsilon_eq(r["stiff"], N, T) for r in rr]
            ent[start] = {"n": len(rr), "U": mean_se([r["U"]["mean"] for r in rr]), "Ui": mean_se([r["Ui"]["mean"] for r in rr]),
                          "q": mean_se([r["q"]["mean"] for r in rr]), "D": mean_se([r["D"]["mean"] for r in rr]),
                          "mu": mean_se([r["mu"]["mean"] for r in rr]), "ups_eq": mean_se(ups),
                          "tau_U": float(np.median([r["U"]["tau"] for r in rr])),
                          "C_all": np.mean([r["C_all"] for r in rr], axis=0).tolist(),
                          "C_int": np.mean([r["C_int"] for r in rr], axis=0).tolist(),
                          "eta_int": eta_from_C(np.mean([r["C_int"] for r in rr], axis=0), 2, L // 4),
                          "J_rel_imb": float(abs(np.mean([r["J_kill"] for r in rr]) - np.mean([r["J_rep"] for r in rr])) / max(np.mean([r["J_kill"] for r in rr]), 1e-9)),
                          "U_per_seed": [r["U"]["mean"] for r in rr], "q_per_seed": [r["q"]["mean"] for r in rr]}
        if "intact" in ent and "damaged" in ent:
            a, b = np.array(ent["intact"]["U_per_seed"]), np.array(ent["damaged"]["U_per_seed"])
            ent["gap_U"] = float(a.mean() - b.mean())
            ent["gap_ci95"] = [float(np.percentile([a[RNG.integers(0, len(a), len(a))].mean() - b[RNG.integers(0, len(b), len(b))].mean() for _ in range(B)], q)) for q in (2.5, 97.5)]
        res["%d|%g|%s|%g" % (L, beta, mode, rho)] = ent
    out["fieldfree"] = res


def twist_response(r, L):
    """Upsilon_tw from the +/- twist runs of one seed: -(j+ - j-)/(2 phi0). Thermally excited vortex pairs polarize in the twist and
    screen the response; that polarization is part of the physical response and is kept. What is not part of it is a *global*
    winding number frozen into a run (an integer: the state of winding W carries the current j = Upsilon (2 pi W/L - phi)), so
    the response is corrected for the integer parts W = round(mean winding) of the two runs when they differ.
    Returns (response, response without any correction, flag), the flag marking seeds whose two runs differ in global winding."""
    ph = r["phi0"]
    sp, sm = r["plus"]["stiff"], r["minus"]["stiff"]
    unc = -(sp["s_x"] - sm["s_x"]) / (2 * ph)
    dW = round(sp["w_mean"]) - round(sm["w_mean"])
    cor = -(sp["s_x"] - sm["s_x"]) / (2 * ph - 2 * math.pi * dW / L)
    return float(cor), float(unc), bool(dW != 0)


def _stiff_group(rr, N, L, T):
    tw = [twist_response(r, L) for r in rr]
    cor = np.array([t[0] for t in tw])
    unc = np.array([t[1] for t in tw])
    flag = np.array([t[2] for t in tw])
    eq = [0.5 * (upsilon_eq(r["plus"]["stiff"], N, T) + upsilon_eq(r["minus"]["stiff"], N, T)) for r in rr]
    return {"n": len(rr), "ups_tw": mean_se(cor), "ups_tw_ci95": boot_ci(cor), "ups_tw_uncorrected": mean_se(unc), "ups_eq": mean_se(eq),
            "n_flagged": int(flag.sum()), "U": mean_se([0.5 * (r["plus"]["U"]["mean"] + r["minus"]["U"]["mean"]) for r in rr]),
            "q": mean_se([0.5 * (r["plus"]["q"]["mean"] + r["minus"]["q"]["mean"]) for r in rr]),
            "ratio_eq_tw": float(np.mean(eq) / cor.mean()) if abs(cor.mean()) > 0.05 else None,
            "tau_s": float(np.median([r["plus"]["stiff"]["tau_s"] for r in rr])), "tw_per_seed": cor.tolist(),
            "phi0": rr[0]["phi0"], "w_mean_abs": float(np.mean([abs(r["plus"]["stiff"]["w_mean"]) for r in rr] + [abs(r["minus"]["stiff"]["w_mean"]) for r in rr])),
            "w_rowstd": float(np.mean([r["plus"]["stiff"]["w_rowstd"] for r in rr]))}


def sec_stiffness(out):
    D = load("stiffness")["runs"]
    res = {}
    keys = sorted({(r["params"]["L"], r["params"]["beta"], r["params"]["repair"], rho_of(r)) for r in D})
    for (L, beta, mode, rho) in keys:
        rr = [r for r in D if r["params"]["L"] == L and r["params"]["beta"] == beta and r["params"]["repair"] == mode and rho_of(r) == rho]
        res["%d|%g|%s|%g|intact" % (L, beta, mode, rho)] = _stiff_group(rr, L * L, L, D[0]["params"]["T"])
    out["stiffness"] = res


def sec_stiffness_dmg(out):
    D = load("stiffness_dmg")["runs"]
    res = {}
    keys = sorted({(r["params"]["L"], r["params"]["beta"], r["params"]["repair"], rho_of(r)) for r in D})
    for (L, beta, mode, rho) in keys:
        rr = [r for r in D if r["params"]["L"] == L and r["params"]["beta"] == beta and r["params"]["repair"] == mode and rho_of(r) == rho]
        res["%d|%g|%s|%g|damaged" % (L, beta, mode, rho)] = _stiff_group(rr, L * L, L, D[0]["params"]["T"])
    out["stiffness_dmg"] = res


def sec_stiff_L(out):
    D = load("stiff_L")["runs"]
    res = {}
    for (L, rho) in sorted({(r["params"]["L"], rho_of(r)) for r in D}):
        rr = [r for r in D if r["params"]["L"] == L and rho_of(r) == rho]
        res["%d|%g" % (L, rho)] = _stiff_group(rr, L * L, L, rr[0]["params"]["T"])
    out["stiff_L"] = res


# ================================================================ timescale competition
def sec_stiff_L2(out):
    import os
    D = load("stiff_L2")["runs"]
    if os.path.exists(os.path.join(A.RES, "stiff_L3_raw.json")):          # the finer scan across the two thresholds
        D = D + load("stiff_L3")["runs"]
    res = {}
    for (mode, L, rho) in sorted({(r["params"]["repair"], r["params"]["L"], rho_of(r)) for r in D}):
        rr = [r for r in D if r["params"]["repair"] == mode and r["params"]["L"] == L and rho_of(r) == rho]
        g = _stiff_group(rr, L * L, L, rr[0]["params"]["T"])
        g["tw_ci95"] = boot_ci(g["tw_per_seed"])
        res["%s|%d|%g" % (mode, L, rho)] = g
    out["stiff_L2"] = res


def rho_hazard(R_, lam0):
    """Ratio of the continuous-time hazard rates equivalent to the per-step probabilities: ln(1 - R)/ln(1 - lam0). It equals
    R/lam0 only when both probabilities are small."""
    return math.log(1.0 - R_) / math.log(1.0 - lam0)


def sec_timescale(out):
    """Is rho enough? U against (rho, lam0, relax) for the three repair rules. Omega = lam0/relax is the number of damage events
    per site per phase sweep. Three questions: (i) how much of the lam0-dependence at fixed rho is the discrete update (R not
    small); (ii) how much is the phase-reset injection of the original rule; (iii) what is left."""
    D = load("timescale")["runs"]
    res = {"cells": {}}
    keys = sorted({(r["params"]["repair"], rho_of(r), r["params"]["lam0"], r["params"]["relax"]) for r in D})
    for (mode, rho, lam, rel) in keys:
        rr = [r for r in D if r["params"]["repair"] == mode and rho_of(r) == rho and r["params"]["lam0"] == lam and r["params"]["relax"] == rel]
        R_ = rr[0]["params"]["R"]
        res["cells"]["%s|%g|%g|%d" % (mode, rho, lam, rel)] = {
            "n": len(rr), "U": mean_se([r["U"]["mean"] for r in rr]), "q": mean_se([r["q"]["mean"] for r in rr]), "omega": lam / rel,
            "R": R_, "rho_hazard": rho_hazard(R_, lam), "max_kill_prob": lam * 21.0,
            "saturated": bool(lam * 21.0 > 1.0)}
    coll = {}
    for mode in ("reset0", "keep", "neighbor"):
        cells = [(k.split("|"), v) for k, v in res["cells"].items() if k.startswith(mode + "|")]
        q = np.array([v["q"][0] for kk, v in cells])
        U = np.array([v["U"][0] for kk, v in cells])
        om = np.array([v["omega"] for kk, v in cells])
        rho = np.array([float(kk[1]) for kk, v in cells])
        lam = np.array([float(kk[2]) for kk, v in cells])
        rel = np.array([int(kk[3]) for kk, v in cells])
        ok = lam < 0.07                                   # exclude the saturated lam0 = 0.08 cells from the fits
        out_m = {}
        # (a) spread of U at fixed nominal rho over lam0 (relax = 2) and over relax (lam0 = 0.02)
        for r_ in (2.5, 5.0):
            m1 = (rho == r_) & (rel == 2)
            m2 = (rho == r_) & (np.abs(lam - 0.02) < 1e-9)
            out_m["range_lam0_at_rho%g" % r_] = [float(U[m1].min()), float(U[m1].max())]
            out_m["range_relax_at_rho%g" % r_] = [float(U[m2].min()), float(U[m2].max())]
        # (b) master curve U(q): polynomial of degree 3 over all unsaturated cells, residual against log Omega
        c = np.polyfit(q[ok], U[ok], 3)
        resid = U[ok] - np.polyval(c, q[ok])
        out_m["master_q_rms"] = float(np.sqrt(np.mean(resid ** 2)))
        out_m["master_q_max"] = float(np.max(np.abs(resid)))
        out_m["resid_vs_logomega_corr"] = float(np.corrcoef(np.log10(om[ok]), resid)[0, 1])
        out_m["resid_slope_per_decade_omega"] = float(np.polyfit(np.log10(om[ok]), resid, 1)[0])
        # (c) at fixed nominal rho and fixed lam0 (only relax varies): U against relax, with q
        small = {}
        for r_ in (2.5, 5.0):
            for l_ in (0.005, 0.01, 0.02):
                m = (rho == r_) & (np.abs(lam - l_) < 1e-9)
                small["%g|%g" % (r_, l_)] = {"relax": rel[m].tolist(), "U": U[m].tolist(), "q": q[m].tolist()}
        out_m["fixed_lam0"] = small
        # (d) hazard-equivalent ratio against nominal
        out_m["rho_hazard_over_rho_at_lam0.04"] = float(np.mean([v["rho_hazard"] / float(kk[1]) for kk, v in cells if abs(float(kk[2]) - 0.04) < 1e-9]))
        out_m["rho_hazard_over_rho_at_lam0.08"] = float(np.mean([v["rho_hazard"] / float(kk[1]) for kk, v in cells if abs(float(kk[2]) - 0.08) < 1e-9]))
        coll[mode] = out_m
    res["collapse"] = coll
    # the slow-damage limit: the smallest Omega in the grid (lam0 = 0.005, relax = 8) against the baseline (lam0 = 0.02, relax = 2)
    lim = {}
    for rho in (2.5, 5.0):
        ent = {"omega_slow": 0.005 / 8, "omega_base": 0.02 / 2}
        for mode in ("reset0", "keep", "neighbor"):
            cs, cb = res["cells"]["%s|%g|0.005|8" % (mode, rho)], res["cells"]["%s|%g|0.02|2" % (mode, rho)]
            ent[mode] = {"U_slow": cs["U"], "q_slow": cs["q"], "U_base": cb["U"], "q_base": cb["q"]}
        Us = [ent[m]["U_slow"][0] for m in ("reset0", "keep", "neighbor")]
        Ub = [ent[m]["U_base"][0] for m in ("reset0", "keep", "neighbor")]
        ent["spread_slow"] = float(max(Us) - min(Us))
        ent["spread_base"] = float(max(Ub) - min(Ub))
        # convergence along relax at lam0 = 0.005 for keep (increments halve as Omega halves?)
        ent["keep_U_vs_relax_lam0.005"] = [res["cells"]["keep|%g|0.005|%d" % (rho, r)]["U"][0] for r in (1, 2, 4, 8)]
        lim["%g" % rho] = ent
    res["slow_limit"] = lim
    out["timescale"] = res


def sec_timescale_s(out):
    D = load("timescale_s")["runs"]
    res = {}
    for (lam, rel) in sorted({(r["params"]["lam0"], r["params"]["relax"]) for r in D}):
        rr = [r for r in D if r["params"]["lam0"] == lam and r["params"]["relax"] == rel]
        g = _stiff_group(rr, rr[0]["params"]["L"] ** 2, rr[0]["params"]["L"], rr[0]["params"]["T"])
        g["omega"] = lam / rel
        res["%g|%d" % (lam, rel)] = g
    out["timescale_s"] = res


def sec_sensitivity(out):
    D = load("sensitivity")["runs"]
    res = {}
    for name, v, rho in sorted({(r["tag"], r["params"][r["tag"]], rho_of(r)) for r in D}):
        rr = [r for r in D if r["tag"] == name and r["params"][name] == v and rho_of(r) == rho]
        res["%s|%g|%g" % (name, v, rho)] = {"n": len(rr), "U": mean_se([r["U"]["mean"] for r in rr]),
                                            "Ui": mean_se([r["Ui"]["mean"] for r in rr]), "q": mean_se([r["q"]["mean"] for r in rr])}
    out["sensitivity"] = res


def sec_tauphi(out):
    D = load("tauphi")["runs"]
    res = {}
    for (h0, mode, rho) in sorted({(r["params"]["h0"], r["params"]["repair"], rho_of(r)) for r in D}):
        rr = [r for r in D if r["params"]["h0"] == h0 and r["params"]["repair"] == mode and rho_of(r) == rho]
        res["%g|%s|%g" % (h0, mode, rho)] = {"n": len(rr), "tau": mean_se([r["tau_sweeps"] for r in rr]), "U": mean_se([r["U"] for r in rr])}
    out["tauphi"] = res


def sec_closure0(out):
    """The kill multiplier the lattice applies (measured in the kill step) against the closure mu(q), and the independent
    prediction obtained from lattices in which damage is independent of phase (beta = 0)."""
    D = load("closure0")["runs"]
    import os
    if os.path.exists(os.path.join(A.RES, "closure0b_raw.json")):          # beta = 0 at small rho (high damaged fractions)
        D = D + load("closure0b")["runs"]
    lam = D[0]["params"]["lam0"]
    res = {}
    for (h0, beta, rho) in sorted({(r["params"]["h0"], r["params"]["beta"], rho_of(r)) for r in D}):
        rr = [r for r in D if r["params"]["h0"] == h0 and r["params"]["beta"] == beta and rho_of(r) == rho]
        q, mk, mr = (np.mean([r[k] for r in rr]) for k in ("q_pre", "mu_kill", "mu_ref_kill"))
        res["%g|%g|%g" % (h0, beta, rho)] = {
            "n": len(rr), "q": mean_se([r["q_pre"] for r in rr]), "mu_kill": mean_se([r["mu_kill"] for r in rr]),
            "mu_ref_kill": mean_se([r["mu_ref_kill"] for r in rr]), "U": mean_se([r["U"]["mean"] for r in rr]),
            "D_sampled": mean_se([r["D"]["mean"] for r in rr]),
            "mu_closure_q": float(MF.mu(q, beta)) if beta > 0 else 1.0,
            "rho_identity": float(MF.rho_map(q, beta, lam, mubar=mk))}
    out["closure0"] = res


def sec_thresholds(out):
    """Thermodynamic-limit bracket of the stiffness threshold from the size scans (field-free, beta = 20, ordered start, T = 0.35).
    Criterion, stated in the manuscript: at a given rho the stiffness has *vanished* if the twist response at the largest size
    is below 0.05 (with at least three sizes); it is *finite* if it exceeds 0.15 at the largest size and does not fall
    by more than a factor 0.7 from the smallest to the largest size. The bracket is (largest vanishing rho, smallest finite rho)."""
    S1, S2 = out["stiff_L"], out["stiff_L2"]
    NK = 2 * 0.35 / math.pi
    res = {"NK": NK, "criterion": {"vanished": "Y(Lmax) < 0.05 (at least three sizes)", "finite": "Y(Lmax) > 0.15 and Y(Lmax) >= 0.7 Y(Lmin)"}}
    for mode in ("keep", "neighbor"):
        table = {}
        for key, v in S2.items():
            m, L, rho = key.split("|")
            if m == mode:
                table.setdefault(float(rho), {})[int(L)] = v["ups_tw"]
        if mode == "keep":
            for key, v in S1.items():
                L, rho = key.split("|")
                if float(rho) in (5.0, 7.0, 10.0, 20.0):
                    table.setdefault(float(rho), {}).setdefault(int(L), v["ups_tw"])
        rows, vanished, finite = {}, [], []
        for rho in sorted(table):
            d = table[rho]
            Ls = sorted(d)
            lo, hi = d[Ls[0]][0], d[Ls[-1]][0]
            st = "vanished" if (hi < 0.05 and len(Ls) >= 3) else ("finite" if (hi > 0.15 and hi >= 0.7 * lo) else "undecided")
            rows["%g" % rho] = {"L": {str(L): d[L] for L in Ls}, "status": st}
            (vanished if st == "vanished" else finite if st == "finite" else []).append(rho)
        res[mode] = {"rows": rows, "rho_vanished_max": max(vanished) if vanished else None, "rho_finite_min": min(finite) if finite else None}
        if finite:
            r0 = min(finite)
            Ls = sorted(table[r0])
            res[mode]["Y_at_lowest_finite_rho_Lmax"] = table[r0][Ls[-1]]
            res[mode]["Lmax_at_lowest_finite"] = Ls[-1]
    out["thresholds"] = res


A.SECTIONS.update({"branches": sec_branches, "ramps": sec_ramps, "fieldfree": sec_fieldfree, "stiffness": sec_stiffness,
                   "stiff_L": sec_stiff_L, "thresholds": sec_thresholds, "stiff_L2": sec_stiff_L2, "stiffness_dmg": sec_stiffness_dmg, "timescale": sec_timescale, "timescale_s": sec_timescale_s, "sensitivity": sec_sensitivity,
                   "tauphi": sec_tauphi, "closure0": sec_closure0})

if __name__ == "__main__":
    import sys
    out = A.load_out()
    for s in sys.argv[1:]:
        A.SECTIONS[s](out)
        print("analysed", s)
    A.save_out(out)
