# -*- coding: utf-8 -*-
"""Reduce the raw runs to the numbers the manuscript reports (results/results.json).

    python analyze.py <section> [<section> ...]       sections: steady converge branches ramps fieldfree stiffness stiff_L
                                                       timescale sensitivity tauphi fss cluster tcurve closure0 anneal mapping

Statistics: every interval is a 95% percentile bootstrap interval (B = 2000) over independent runs (seeds or disorder
realizations); standard errors are seed standard deviations over sqrt(n). Autocorrelation times are Sokal's automatic
window (W >= 6 tau_int) and are reported, but the error of a point is always the spread over independent seeds.
"""
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit

import meanfield as MF

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
RNG = np.random.default_rng(20261005)
B = 2000
OUTP = os.path.join(RES, "results.json")


def load(study):
    with open(os.path.join(RES, "%s_raw.json" % study)) as f:
        return json.load(f)


def load_out():
    return json.load(open(OUTP)) if os.path.exists(OUTP) else {}


def save_out(out):
    with open(OUTP, "w") as f:
        json.dump(out, f, indent=1)


def rho_of(r):
    return round(r["params"]["R"] / r["params"]["lam0"], 3)


def mean_se(x):
    x = np.asarray(x, float)
    n = len(x)
    return float(np.mean(x)), float(np.std(x, ddof=1) / math.sqrt(n)) if n > 1 else float("nan")


def boot_ci(x, stat=np.mean, b=B):
    x = np.asarray(x, float)
    idx = RNG.integers(0, len(x), (b, len(x)))
    v = np.array([stat(x[i]) for i in idx])
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def group(runs, keyf, valf):
    d = defaultdict(list)
    for r in runs:
        d[keyf(r)].append(valf(r))
    return d


def holm(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    out, run = [0.0] * m, 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (m - k) * pvals[i]))
        out[i] = run
    return out


# ================================================================ steady, convergence, closure (h0 = 0.4, original repair)
def sec_steady(out):
    D = load("steady")["runs"]
    lam = D[0]["params"]["lam0"]
    res = {"long": {}, "short": {}}
    for tag in ("long", "short"):
        runs = [r for r in D if r["tag"] == tag]
        keys = sorted({(r["params"]["L"], r["params"]["beta"], rho_of(r)) for r in runs})
        for (L, beta, rho) in keys:
            rr = [r for r in runs if r["params"]["L"] == L and r["params"]["beta"] == beta and rho_of(r) == rho]
            U, q, Dd, mu = [[r[k]["mean"] for r in rr] for k in ("U", "q", "D", "mu")]
            Jk, Jr = [r["J_kill"] for r in rr], [r["J_rep"] for r in rr]
            res[tag]["%d|%g|%g" % (L, beta, rho)] = {
                "n": len(rr), "U": mean_se(U), "q": mean_se(q), "D": mean_se(Dd), "mu": mean_se(mu),
                "U_sd": float(np.std(U, ddof=1)), "tau_U": float(np.median([r["U"]["tau"] for r in rr])),
                "J_kill": float(np.mean(Jk)), "J_rep": float(np.mean(Jr)),
                "J_rel_imbalance": float(abs(np.mean(Jk) - np.mean(Jr)) / np.mean(Jk))}
    # short protocol against long, L = 32
    cmp_ = []
    for k, v in res["short"].items():
        if k in res["long"]:
            cmp_.append({"key": k, "short": v["U"][0], "long": res["long"][k]["U"][0], "long_se": res["long"][k]["U"][1],
                         "diff": v["U"][0] - res["long"][k]["U"][0]})
    res["short_vs_long"] = cmp_
    # intact-fraction references at beta = 0: continuous time rho/(1+rho) and the implemented discrete map
    res["lam0"] = lam
    # closure test: lattice kill multiplier against the closure mu(q) and mu(D)
    cl = {}
    for beta in (5.0, 20.0):
        rows = []
        for key, v in res["long"].items():
            L, b, rho = key.split("|")
            if int(L) != 64 or float(b) != beta:
                continue
            q, mub, Dm = v["q"][0], v["mu"][0], v["D"][0]
            rows.append({"rho": float(rho), "q": q, "D_intact": Dm, "mu_lattice": mub,
                         "mu_closure_q": float(MF.mu(q, beta)), "mu_at_mean_D": float(MF.mu(Dm, beta)),
                         "rho_map_check": float(MF.rho_map(q, beta, lam, mubar=mub))})
        rows.sort(key=lambda r: r["rho"])
        cl[str(int(beta))] = rows
    res["closure"] = cl
    res["windows"] = {str(int(b)): {"continuous": MF.window(b), "discrete": MF.window(b, lam0=lam)} for b in (5.0, 20.0)}
    res["beta_c_h6"] = MF.beta_c()
    out["steady"] = res


def sec_converge(out):
    R = load("converge")["runs"]
    res = {}
    for rho in sorted({rho_of(r) for r in R}):
        for start in ("intact", "damaged"):
            a = np.array([r["series"]["U"] for r in R if rho_of(r) == rho and r["kw"]["start"] == start])
            ev = R[0]["series"]["every"]
            res["%g|%s" % (rho, start)] = {"t": (np.arange(a.shape[1]) * ev).tolist(), "mean": a.mean(0).tolist(),
                                           "se": (a.std(0, ddof=1) / math.sqrt(len(a))).tolist(), "n": len(a)}
        t = np.array(res["%g|intact" % rho]["t"])
        mi, md = np.array(res["%g|intact" % rho]["mean"]), np.array(res["%g|damaged" % rho]["mean"])
        late = float(np.mean(np.concatenate([mi[t >= 4000], md[t >= 4000]])))
        late_gap = float(np.mean(mi[t >= 4000]) - np.mean(md[t >= 4000]))
        # first time each start comes within 0.02 of the late-time value and stays within 0.03 afterwards
        def enter(m):
            for k in range(len(m)):
                if abs(m[k] - late) < 0.02 and (np.abs(m[k:] - late) < 0.03).all():
                    return int(t[k])
        res["%g|U_late" % rho] = late
        res["%g|late_gap" % rho] = late_gap
        res["%g|t_enter_intact" % rho] = enter(mi)
        res["%g|t_enter_damaged" % rho] = enter(md)
    out["converge"] = res


SECTIONS = {"steady": sec_steady, "converge": sec_converge}

if __name__ == "__main__":
    out = load_out()
    for s in sys.argv[1:]:
        SECTIONS[s](out)
        print("analysed", s)
    save_out(out)
