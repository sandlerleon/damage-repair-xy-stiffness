# -*- coding: utf-8 -*-
"""Third part of the analysis: quenched-dilution finite-size scaling, clustered damage, the BKT line T_BKT(f) and the
mapping of a damaged fraction to displacement damage and of anneal temperature to recovery."""
import math

import numpy as np

import analyze as A
from analyze import B, RNG, load, mean_se

C_GRID = np.linspace(-2, 30, 641)


def crossing(x, y, level):
    """First crossing of y(x) with `level` by linear interpolation (NaN if none)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    for i in range(len(x) - 1):
        if (y[i] - level) * (y[i + 1] - level) <= 0 and y[i] != y[i + 1]:
            return float(x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return float("nan")


def lnL2_fit(Ls, vals):
    X = 1.0 / np.log(np.asarray(Ls, float)) ** 2
    c = np.polyfit(X, vals, 1)
    return float(c[1]), float(c[0])


def collect(runs, T, clusters=0):
    """Y[L][f] -> array of per-realization helicity moduli."""
    Y = {}
    for r in runs:
        if abs(r["T"] - T) > 1e-9 or r.get("clusters", 0) != clusters:
            continue
        Y.setdefault(r["L"], {}).setdefault(round(r["f"], 4), []).append(r["Y"])
    return {L: {f: np.array(v) for f, v in d.items()} for L, d in Y.items()}


def crossings_by_size(Y, T):
    NK = 2 * T / math.pi
    res = {}
    for L in sorted(Y):
        fs = sorted(Y[L])
        if len(fs) < 3:
            continue
        mean = np.array([Y[L][f].mean() for f in fs])
        fk = crossing(fs, mean, NK)
        boots = []
        for _ in range(B):
            m = np.array([Y[L][f][RNG.integers(0, len(Y[L][f]), len(Y[L][f]))].mean() for f in fs])
            c = crossing(fs, m, NK)
            if np.isfinite(c):
                boots.append(c)
        res[L] = {"f_KT": fk, "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))] if boots else None,
                  "n_real": int(min(len(Y[L][f]) for f in fs)), "n_f": len(fs), "f": fs, "Y": mean.tolist(),
                  "Y_se": [float(Y[L][f].std(ddof=1) / math.sqrt(len(Y[L][f]))) for f in fs]}
    return res


def extrapolate(Y, T, per_size, Lmin=16):
    NK = 2 * T / math.pi
    Ls = [L for L in sorted(per_size) if L >= Lmin and np.isfinite(per_size[L]["f_KT"])]
    if len(Ls) < 3:
        return None
    fk = np.array([per_size[L]["f_KT"] for L in Ls])
    f_inf, slope = lnL2_fit(Ls, fk)
    boots = []
    for _ in range(B // 2):
        fb = []
        for L in Ls:
            fs = sorted(Y[L])
            m = np.array([Y[L][f][RNG.integers(0, len(Y[L][f]), len(Y[L][f]))].mean() for f in fs])
            fb.append(crossing(fs, m, NK))
        if np.all(np.isfinite(fb)):
            boots.append(lnL2_fit(Ls, fb)[0])
    return {"sizes": Ls, "f_inf": f_inf, "slope": slope,
            "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))] if boots else None,
            "scatter_sd_over_sizes": float(np.std(fk, ddof=1)), "mean_over_sizes": float(np.mean(fk))}


def weber_minnhagen(Y, T, Lmin=24):
    """Upsilon_L(f) = (2T/pi)(1 + 1/(2 ln L + C)) at the true f_KT, with C free."""
    NK = 2 * T / math.pi
    Ls = [L for L in sorted(Y) if L >= Lmin and len(Y[L]) >= 6]
    if len(Ls) < 3:
        return None
    fs_common = sorted(set.intersection(*[set(Y[L]) for L in Ls]))
    grid = np.arange(fs_common[0], fs_common[-1] + 1e-9, 0.0025)
    lnL = np.log(np.array(Ls, float))
    pred = NK * (1 + 1 / (2 * lnL[None, :] + C_GRID[:, None]))

    def best(Ym):
        chis = []
        for f in grid:
            y = np.array([np.interp(f, fs_common, [Ym[L][g].mean() if not isinstance(Ym[L][g], float) else Ym[L][g] for g in fs_common]) for L in Ls])
            s = np.array([max(np.interp(f, fs_common, [Y[L][g].std(ddof=1) / math.sqrt(len(Y[L][g])) for g in fs_common]), 1e-4) for L in Ls])
            c = np.sum(((y[None, :] - pred) / s[None, :]) ** 2, axis=1)
            chis.append((float(c.min()), float(C_GRID[int(np.argmin(c))])))
        k = int(np.argmin([c[0] for c in chis]))
        return float(grid[k]), chis[k][1], chis[k][0]
    f0, C0, chi = best(Y)
    boots = []
    for _ in range(150):
        Yb = {L: {g: Y[L][g][RNG.integers(0, len(Y[L][g]), len(Y[L][g]))] for g in fs_common} for L in Ls}
        boots.append(best(Yb)[0])
    return {"f_KT": f0, "C": C0, "chi2": chi, "dof": len(Ls) - 1, "sizes": Ls, "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]}


def sec_fss(out):
    runs = load("fss")["runs"]
    res = {}
    for T in sorted({r["T"] for r in runs}):
        Y = collect(runs, T)
        per = crossings_by_size(Y, T)
        ent = {"T": T, "NK": 2 * T / math.pi, "per_L": {str(L): v for L, v in per.items()}}
        ent["extrap_all"] = extrapolate(Y, T, per)
        ent["sensitivity_lnL2"] = {str(Lmin): extrapolate(Y, T, per, Lmin) for Lmin in (16, 24, 32, 48)}
        if T > 0.3:
            ent["weber_minnhagen"] = {str(Lmin): weber_minnhagen(Y, T, Lmin) for Lmin in (24, 32)}
        # binder cumulant and magnetisation moments for the figure
        bind = {}
        for L in sorted(Y):
            rr = [r for r in runs if abs(r["T"] - T) < 1e-9 and r["L"] == L and r.get("clusters", 0) == 0]
            fsx = sorted({round(r["f"], 4) for r in rr})
            bind[str(L)] = {"f": fsx, "U4": [1 - np.mean([r["M4"] for r in rr if round(r["f"], 4) == f]) /
                                             (3 * np.mean([r["M2"] for r in rr if round(r["f"], 4) == f]) ** 2) for f in fsx],
                            "tau_Y_max": float(max(r["tau_Y"] for r in rr))}
        ent["binder"] = bind
        res["%g" % T] = ent
    out["fss"] = res


def sec_cluster(out):
    runs = load("cluster")["runs"]
    res = {}
    T = 0.35
    NK = 2 * T / math.pi
    for rad in sorted({r["clusters"] for r in runs}):
        for L in sorted({r["L"] for r in runs if r["clusters"] == rad}):
            Y = collect([r for r in runs if r["L"] == L], T, clusters=rad)
            per = crossings_by_size(Y, T)
            if L in per:
                res["%d|%d" % (rad, L)] = per[L]
    out["cluster"] = res


def sec_tcurve(out):
    runs = load("tcurve")["runs"]
    res = {}
    fs = sorted({round(r["f"], 4) for r in runs})
    Ts = sorted({r["T"] for r in runs})
    Ls = sorted({r["L"] for r in runs})
    for f in fs:
        perL = {}
        for L in Ls:
            Y = {T: np.array([r["Y"] for r in runs if r["L"] == L and abs(r["f"] - f) < 1e-9 and abs(r["T"] - T) < 1e-9]) for T in Ts}
            Tl = crossing_T(Ts, [Y[T].mean() for T in Ts])
            boots = []
            for _ in range(B // 2):
                m = [Y[T][RNG.integers(0, len(Y[T]), len(Y[T]))].mean() for T in Ts]
                c = crossing_T(Ts, m)
                if np.isfinite(c):
                    boots.append(c)
            perL[L] = {"T_L": Tl, "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))] if boots else None,
                       "Y": [float(Y[T].mean()) for T in Ts]}
        valid = [L for L in Ls if np.isfinite(perL[L]["T_L"])]
        ext = None
        if len(valid) >= 3:
            vals = [perL[L]["T_L"] for L in valid]
            Tinf, slope = lnL2_fit(valid, vals)
            boots = []
            for _ in range(B // 4):
                vb = []
                for L in valid:
                    Yb = [np.array([r["Y"] for r in runs if r["L"] == L and abs(r["f"] - f) < 1e-9 and abs(r["T"] - T) < 1e-9]) for T in Ts]
                    m = [y[RNG.integers(0, len(y), len(y))].mean() for y in Yb]
                    vb.append(crossing_T(Ts, m))
                if np.all(np.isfinite(vb)):
                    boots.append(lnL2_fit(valid, vb)[0])
            ext = {"T_BKT": Tinf, "slope": slope, "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))] if boots else None,
                   "sizes": valid}
        res["%g" % f] = {"per_L": {str(L): v for L, v in perL.items()}, "extrap": ext}
    res["T_grid"] = Ts
    out["tcurve"] = res


def crossing_T(Ts, Y):
    """Crossing of Upsilon(T) with 2T/pi: Upsilon decreases with T, 2T/pi increases."""
    Ts, Y = np.asarray(Ts, float), np.asarray(Y, float)
    g = Y - 2 * Ts / math.pi
    for i in range(len(Ts) - 1):
        if g[i] >= 0 > g[i + 1]:
            return float(Ts[i] + g[i] * (Ts[i + 1] - Ts[i]) / (g[i] - g[i + 1]))
    return float("nan")


A.SECTIONS.update({"fss": sec_fss, "cluster": sec_cluster, "tcurve": sec_tcurve})

if __name__ == "__main__":
    import sys
    import analyze_b  # noqa: F401  (registers the sections of part B)
    out = A.load_out()
    for s in sys.argv[1:]:
        A.SECTIONS[s](out)
        print("analysed", s)
    A.save_out(out)
