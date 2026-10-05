# -*- coding: utf-8 -*-
"""Fifth part of the analysis: control simulations for the twist-response stiffness (study validation)."""
import math

import numpy as np

import analyze as A
import analyze_b as AB


def _tw(r, L):
    """Twist response of one seed pair, as in analyze_b.twist_response, from the validation records."""
    ph = r["phi0"]
    sp, sm = r["plus"], r["minus"]
    dW = round(sp["w_mean"]) - round(sm["w_mean"])
    unc = -(sp["s_x"] - sm["s_x"]) / (2 * ph)
    cor = -(sp["s_x"] - sm["s_x"]) / (2 * ph - 2 * math.pi * dW / L)
    return float(cor), float(unc), bool(dW != 0)


def sec_validation(out):
    D = A.load("validation")["runs"]
    T = A.load("validation")["meta"]["T"]
    res = {"T": T, "harmonic_clean": 1.0 - T / 4.0}
    clean = {}
    for L in sorted({r["params"]["L"] for r in D if r["kind"] == "twist_clean"}):
        for tag in sorted({r["tag"] for r in D if r["kind"] == "twist_clean" and r["params"]["L"] == L}, key=lambda t: float(t[1:])):
            rr = [r for r in D if r["kind"] == "twist_clean" and r["params"]["L"] == L and r["tag"] == tag]
            tw = [_tw(r, L) for r in rr]
            eq = [0.5 * (AB.upsilon_eq(r["plus"], L * L, T) + AB.upsilon_eq(r["minus"], L * L, T)) for r in rr]
            clean["%d|%s" % (L, tag[1:])] = {"n": len(rr), "c": float(tag[1:]), "phi0": rr[0]["phi0"], "ups_tw": A.mean_se([t[0] for t in tw]),
                                             "ups_eq": A.mean_se(eq), "n_flagged": int(sum(t[2] for t in tw)),
                                             "w_abs_max": float(max(abs(r["plus"]["w_mean"]) for r in rr))}
    res["clean"] = clean
    def dil_block(D, T):
        dil = {}
        for tag in sorted({r["tag"] for r in D if r["kind"] == "twist_diluted"}):
            rr = [r for r in D if r["kind"] == "twist_diluted" and r["tag"] == tag]
            st = [r for r in D if r["kind"] == "static_eq" and r["tag"] == tag]
            L = rr[0]["params"]["L"]
            tw = [_tw(r, L)[0] for r in rr]
            eq = [0.5 * (AB.upsilon_eq(r["plus"], L * L, T) + AB.upsilon_eq(r["minus"], L * L, T)) for r in rr]
            sd = {tuple(r["seed"]): r["ups_eq_static"] for r in st}      # same coupling maps (same seed code) in both codes
            dd = [t - sd[tuple(r["seed"])] for t, r in zip(tw, rr)]
            de = [t - e for t, e in zip(tw, eq)]
            dil[tag[1:]] = {"n": len(rr), "ups_tw": A.mean_se(tw), "ups_eq_same_run": A.mean_se(eq), "ups_eq_static": A.mean_se(list(sd.values())),
                            "diff_tw_minus_static": A.mean_se(dd), "diff_tw_minus_eq_same_run": A.mean_se(de)}
        return dil
    res["diluted"] = dil_block(D, T)
    import os
    if os.path.exists(os.path.join(A.RES, "validation_long_raw.json")):
        res["diluted_long"] = dil_block(A.load("validation_long")["runs"], T)
    wd = {}
    for tag in ("w00", "w11", "w10"):
        rr = [r for r in D if r["kind"] == "twist_wound" and r["tag"] == tag]
        tw = [_tw(r, 32) for r in rr]
        wd[tag] = {"n": len(rr), "ups_tw": A.mean_se([t[0] for t in tw]), "ups_uncorrected": A.mean_se([t[1] for t in tw]),
                   "n_flagged": int(sum(t[2] for t in tw)), "w_plus": float(np.mean([r["plus"]["w_mean"] for r in rr])),
                   "w_minus": float(np.mean([r["minus"]["w_mean"] for r in rr]))}
    res["wound"] = wd
    out["validation"] = res


A.SECTIONS["validation"] = sec_validation

if __name__ == "__main__":
    import sys
    import analyze_b  # noqa: F401
    import analyze_c  # noqa: F401
    import analyze_d  # noqa: F401
    out = A.load_out()
    for s in sys.argv[1:]:
        A.SECTIONS[s](out)
        print("analysed", s)
    A.save_out(out)
