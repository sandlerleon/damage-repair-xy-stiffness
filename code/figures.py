# -*- coding: utf-8 -*-
"""Figures of the manuscript, from results/results.json (analyze*.py) and the raw runs.

    python figures.py <fig> [<fig> ...]        fig1 fig2 fig3 fig4 fig5 fig6 fig7 fig8
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import meanfield as MF

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8.5,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.fontsize": 7, "lines.linewidth": 1.2})
BETACOL = {0: "#1a7a4c", 5: "#e08e0b", 20: "#c0392b"}
MODECOL = {"reset0": "#7f7f7f", "keep": "#08519c", "neighbor": "#c0392b"}


def R():
    return json.load(open(os.path.join(RES, "results.json")))


def tag(ax, t):
    ax.text(-0.14, 1.06, t, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def fig1():
    """Steady-state order, damaged fraction against the corrected reduction, and the closure test."""
    S = R()["steady"]
    lam = S["lam0"]
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.5))
    fig.subplots_adjust(wspace=0.34)
    rhos = sorted({float(k.split("|")[2]) for k in S["long"]})
    ax = axs[0]
    for beta in (0, 5, 20):
        for L, mk in ((32, "o"), (64, "s")):
            y = [S["long"]["%d|%d|%g" % (L, beta, r)]["U"][0] for r in rhos]
            e = [S["long"]["%d|%d|%g" % (L, beta, r)]["U"][1] for r in rhos]
            ax.errorbar(rhos, y, yerr=np.array(e) * 1.96, fmt=mk + "-", ms=3.2, lw=0.8, color=BETACOL[beta], mfc=BETACOL[beta] if L == 32 else "white",
                        label=(r"$\beta$ = %d" % beta) if L == 32 else None)
    ax.plot([], [], "ko-", ms=3, mfc="k", lw=0.8, label="L = 32 (filled)")
    ax.plot([], [], "ks-", ms=3, mfc="white", lw=0.8, label="L = 64 (open)")
    rr = np.geomspace(0.9, 45, 200)
    ax.plot(rr, rr / (1 + rr), "k--", lw=0.9, label=r"intact fraction $\rho/(1+\rho)$ (rate equation)")
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"repair-to-damage ratio $\rho = R/\lambda_0$")
    ax.set_ylabel(r"phase coherence $U$")
    ax.set_title(r"Steady-state order ($h_0$ = 0.4, repair $\theta \to 0$)", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=5.8, loc="lower right", ncol=1, handlelength=1.6)
    tag(ax, "A")
    ax = axs[1]
    for beta in (0, 5, 20):
        q = [S["long"]["64|%d|%g" % (beta, r)]["q"][0] for r in rhos]
        ax.plot(rhos, q, "o", ms=3.2, color=BETACOL[beta], label=r"lattice, $\beta$ = %d" % beta)
        qq = np.linspace(1e-4, 1 - 1e-4, 4000)
        rho_c = MF.rho_ode(qq, beta) if beta > 0 else (1 - qq) / qq
        if beta == 0:
            ax.plot(rr, 1 / (1 + rr), "-", color=BETACOL[beta], lw=0.9, alpha=0.8)
        else:
            win = S["windows"][str(beta)]["continuous"]
            # stable branches of the closure model
            ok = np.ones_like(qq, bool)
            dr = np.gradient(rho_c, qq)
            ax.plot(rho_c[dr < 0], qq[dr < 0], "-", color=BETACOL[beta], lw=0.9, alpha=0.8)
            ax.plot(rho_c[dr >= 0], qq[dr >= 0], ":", color=BETACOL[beta], lw=0.9, alpha=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.9, 45)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"damaged fraction $q$")
    ax.set_title("Damaged fraction and the closure model (lines)", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    tag(ax, "B")
    ax = axs[2]
    C0 = R()["closure0"]
    qq = np.linspace(0.02, 0.98, 400)
    ax.plot(qq, MF.mu(qq, 20), "k-", lw=1.1, label=r"closure $\mu(q)$, $\beta$ = 20")
    r20 = sorted(float(k.split("|")[2]) for k in C0 if k.startswith("0.4|20|"))
    r0 = sorted(float(k.split("|")[2]) for k in C0 if k.startswith("0.4|0|"))
    ax.plot([C0["0.4|20|%g" % r]["q"][0] for r in r20], [C0["0.4|20|%g" % r]["mu_kill"][0] for r in r20], "o", color=BETACOL[20], ms=3.6,
            label=r"lattice, $\beta$ = 20 (multiplier applied in the kill step)")
    ax.plot([C0["0.4|0|%g" % r]["q"][0] for r in r0], [C0["0.4|0|%g" % r]["mu_ref_kill"][0] for r in r0], "s", mfc="none", color=BETACOL[0], ms=3.6,
            label=r"$\beta$ = 0 lattices, multiplier a $\beta$ = 20 feedback would apply")
    ax.set_xlabel(r"damaged fraction $q$")
    ax.set_ylabel(r"mean kill multiplier of intact sites $\bar\mu$")
    ax.set_title("Kill multiplier: closure against lattice", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=5.8, loc="upper left")
    tag(ax, "C")
    save(fig, "Fig1_steady_state.png")


def _ff(R_, L, beta, mode, rho):
    return R_["fieldfree"].get("%d|%g|%s|%g" % (L, beta, mode, rho))


RHO_FF = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 40.0]


def fig3():
    """Field-free (h0 = 0) steady states from two starting states, for the three repair rules."""
    Rr = R()
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.5))
    fig.subplots_adjust(wspace=0.34)
    series = [("keep", 20, "keep", r"keep, $\beta$ = 20"), ("neighbor", 20, "neighbor", r"neighbor, $\beta$ = 20"),
              ("reset0", 20, "reset0", r"reset0, $\beta$ = 20"), ("keep", 0, "keep0", r"keep, $\beta$ = 0")]
    cols = {"keep": MODECOL["keep"], "neighbor": MODECOL["neighbor"], "reset0": MODECOL["reset0"], "keep0": "#1a7a4c"}
    ax = axs[0]
    for mode, beta, key, lab in series:
        rs = [r for r in RHO_FF if _ff(Rr, 32, beta, mode, r)]
        for start, mf in (("intact", True), ("damaged", False)):
            y = [_ff(Rr, 32, beta, mode, r)[start]["U"][0] for r in rs]
            ax.plot(rs, y, "o-" if mf else "s--", ms=3.4 if mf else 3.8, lw=0.8, color=cols[key], mfc=cols[key] if mf else "none",
                    label=lab if mf else None)
        if key in ("keep", "neighbor"):
            rs48 = [r for r in RHO_FF if _ff(Rr, 48, beta, mode, r)]
            ax.plot(rs48, [_ff(Rr, 48, beta, mode, r)["intact"]["U"][0] for r in rs48], "^", ms=4.5, color=cols[key], mfc="none", mew=0.9)
    ax.plot([], [], "k^", ms=4.5, mfc="none", label="L = 48")
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"coherence $U$")
    ax.set_title(r"Field-free ($h_0$ = 0), L = 32; filled: intact start, open: damaged start", loc="left", fontsize=7.2)
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    tag(ax, "A")
    ax = axs[1]
    rr = np.geomspace(0.9, 45, 200)
    for mode, beta, key, lab in series:
        rs = [r for r in RHO_FF if _ff(Rr, 32, beta, mode, r)]
        ax.plot(rs, [_ff(Rr, 32, beta, mode, r)["intact"]["q"][0] for r in rs], "o-", ms=3.4, lw=0.8, color=cols[key], label=lab)
    ax.plot(rr, 1 / (1 + rr), "k:", lw=0.8, label=r"$\beta$ = 0 rate equation")
    ax.set_xscale("log")
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"damaged fraction $q$")
    ax.set_title("Damaged fraction", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    tag(ax, "B")
    ax = axs[2]
    for mode, beta, key, lab in series:
        for L, mk in ((32, "o"), (48, "^")):
            rs = [r for r in RHO_FF if _ff(Rr, L, beta, mode, r) and "gap_U" in _ff(Rr, L, beta, mode, r)]
            if not rs:
                continue
            g = np.array([_ff(Rr, L, beta, mode, r)["gap_U"] for r in rs])
            lo = g - np.array([_ff(Rr, L, beta, mode, r)["gap_ci95"][0] for r in rs])
            hi = np.array([_ff(Rr, L, beta, mode, r)["gap_ci95"][1] for r in rs]) - g
            ax.errorbar(np.array(rs) * (1 + 0.02 * (L == 48)), g, yerr=[lo, hi], fmt=mk, ms=3, lw=0.7, capsize=1.5, color=cols[key], mfc=cols[key] if L == 32 else "none")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xscale("log")
    ax.set_ylim(-0.12, 0.12)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$U_{\rm intact\ start} - U_{\rm damaged\ start}$")
    ax.set_title("No dependence on the initial state is resolved (95% CI)", loc="left", fontsize=7.5)
    tag(ax, "C")
    save(fig, "Fig3_field_free.png")


def fig5():
    """Phase correlations at low, intermediate and high repair-to-damage ratio."""
    Rr = R()
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.5))
    fig.subplots_adjust(wspace=0.34)
    for ax, mode, ttl in ((axs[0], "keep", "keep, L = 32, intact pairs"), (axs[1], "neighbor", "neighbor, L = 32, intact pairs")):
        cm = plt.cm.viridis
        for k, rho in enumerate([3.0, 5.0, 6.0, 7.0, 8.0, 10.0, 20.0]):
            e = _ff(Rr, 32, 20, mode, rho)
            if not e:
                continue
            C = np.array(e["intact"]["C_int"])
            r = np.arange(1, len(C) + 1)
            ax.plot(r, np.clip(C, 1e-3, None), "o-", ms=2.6, lw=0.8, color=cm(k / 6.5), label=r"$\rho$ = %g" % rho)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(2e-3, 1.1)
        ax.set_xlabel("separation r")
        ax.set_ylabel(r"$C(r)$")
        ax.set_title(ttl, loc="left", fontsize=8)
        ax.legend(frameon=False, fontsize=6, ncol=2)
    tag(axs[0], "A")
    tag(axs[1], "B")
    ax = axs[2]
    st = Rr.get("stiffness", {})
    for mode in ("keep", "neighbor"):
        rs = [r for r in RHO_FF if _ff(Rr, 32, 20, mode, r)]
        eta = [_ff(Rr, 32, 20, mode, r)["intact"]["eta_int"] for r in rs]
        ax.plot([r for r, e in zip(rs, eta) if e], [e for e in eta if e], "o-", ms=3.4, lw=0.8, color=MODECOL[mode], label=r"from $C(r)$, %s" % mode)
        tw = []
        for r in rs:
            e = st.get("32|20|%s|%g|intact" % (mode, r))
            tw.append(0.35 / (2 * math.pi * e["ups_tw"][0]) if e and e["ups_tw"][0] > 0.03 else None)
        ax.plot([r for r, t in zip(rs, tw) if t], [t for t in tw if t], "s--", ms=3, lw=0.8, color=MODECOL[mode], mfc="none",
                label=r"$T/(2\pi\Upsilon_{\rm tw})$, %s" % mode)
    ax.axhline(0.25, color="k", ls=":", lw=0.8)
    ax.text(40, 0.27, r"$\eta = 1/4$", ha="right", fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.04, 5)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"exponent $\eta$ of $C(r) \sim r^{-\eta}$")
    ax.set_title("Algebraic exponent against stiffness", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    tag(ax, "C")
    save(fig, "Fig5_correlations.png")


def fig4():
    """Phase stiffness from the twist response, against the equilibrium formula and against U."""
    Rr = R()
    St = Rr["stiffness"]
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.5))
    fig.subplots_adjust(wspace=0.36)
    NK = 2 * 0.35 / math.pi
    conds = [("keep", 20, "keep", r"keep, $\beta$ = 20"), ("neighbor", 20, "neighbor", r"neighbor, $\beta$ = 20"), ("keep", 0, "keep0", r"keep, $\beta$ = 0")]
    cols = {"keep": MODECOL["keep"], "neighbor": MODECOL["neighbor"], "keep0": "#1a7a4c"}
    ax = axs[0]
    for mode, beta, key, lab in conds:
        rs = [r for r in RHO_FF if St.get("32|%g|%s|%g|intact" % (beta, mode, r))]
        y = np.array([St["32|%g|%s|%g|intact" % (beta, mode, r)]["ups_tw"][0] for r in rs])
        e = np.array([St["32|%g|%s|%g|intact" % (beta, mode, r)]["ups_tw"][1] for r in rs])
        ax.errorbar(rs, y, yerr=1.96 * e, fmt="o-", ms=3.4, lw=0.9, capsize=1.5, color=cols[key], label=lab)
    ax.axhline(NK, color="k", ls="--", lw=0.8)
    ax.text(1.05, NK + 0.015, r"$2T/\pi$ (equilibrium reference)", fontsize=6.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"twist-response stiffness $\Upsilon_{\rm tw}$")
    ax.set_title("Stiffness, L = 32 (ordered start)", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    tag(ax, "A")
    ax = axs[1]
    for mode, beta, key, lab in conds:
        rs = [r for r in RHO_FF if St.get("32|%g|%s|%g|intact" % (beta, mode, r))]
        tw = np.array([St["32|%g|%s|%g|intact" % (beta, mode, r)]["ups_tw"][0] for r in rs])
        eq = np.array([St["32|%g|%s|%g|intact" % (beta, mode, r)]["ups_eq"][0] for r in rs])
        ax.plot(tw, eq, "o", ms=3.6, color=cols[key], label=lab)
    ax.plot([-0.2, 1], [-0.2, 1], "k-", lw=0.7)
    ax.set_xlim(-0.05, 0.95)
    ax.set_ylim(-0.3, 0.95)
    ax.set_xlabel(r"$\Upsilon_{\rm tw}$ (twist response)")
    ax.set_ylabel(r"$\Upsilon_{\rm eq}$ (equilibrium fluctuation formula)")
    ax.set_title("The equilibrium formula underestimates the response", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    tag(ax, "B")
    ax = axs[2]
    for mode, beta, key, lab in conds:
        rs = [r for r in RHO_FF if St.get("32|%g|%s|%g|intact" % (beta, mode, r))]
        tw = np.array([St["32|%g|%s|%g|intact" % (beta, mode, r)]["ups_tw"][0] for r in rs])
        U = np.array([St["32|%g|%s|%g|intact" % (beta, mode, r)]["U"][0] for r in rs])
        ax.plot(U, tw, "o-", ms=3.4, lw=0.7, color=cols[key], label=lab)
    ax.set_xlabel(r"coherence $U$")
    ax.set_ylabel(r"$\Upsilon_{\rm tw}$")
    ax.set_title(r"$U$ stays near the finite-size floor while $\Upsilon_{\rm tw}$ is zero", loc="left", fontsize=7.5)
    ax.axvline(0.03, color="#888", ls=":", lw=0.8)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")
    tag(ax, "C")
    save(fig, "Fig4_stiffness.png")


def fig2():
    """Tests for bistability and hysteresis (h0 = 0.4, original repair)."""
    Rr = R()
    Bd, Rp = Rr["branches"], Rr["ramps"]
    wins = Rr["steady"]["windows"]["20"]
    fig, axs = plt.subplots(2, 2, figsize=(8.6, 6.6))
    fig.subplots_adjust(hspace=0.42, wspace=0.32)
    ax = axs[0, 0]
    ax.axvspan(wins["continuous"][0], wins["continuous"][1], color="#fdf1c9", zorder=0)
    ax.axvspan(wins["discrete"][0], wins["discrete"][1], color="#f6c98a", alpha=0.45, zorder=0)
    cols = {32: "#9ecae1", 48: "#3182bd", 64: "#08306b"}
    for L in (32, 48, 64):
        rs = sorted(float(k.split("|")[2]) for k in Bd if k != "summary" and k.startswith("%d|20|" % L))
        ax.plot(rs, [Bd["%d|20|%g" % (L, r)]["U_intact"][0] for r in rs], "o-", ms=3.4, lw=0.8, color=cols[L], label="L = %d, intact start" % L)
        ax.plot(rs, [Bd["%d|20|%g" % (L, r)]["U_damaged"][0] for r in rs], "s", ms=4.5, mfc="none", color=cols[L], label="L = %d, damaged start" % L if L == 64 else None)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"stationary $U$")
    ax.set_title(r"Two starts, $\beta$ = 20 (shaded: closure-model windows)", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    tag(ax, "A")
    ax = axs[0, 1]
    for beta, mk in ((5, "o"), (20, "s")):
        for L in (32, 48, 64):
            rs = sorted(float(k.split("|")[2]) for k in Bd if k != "summary" and k.startswith("%d|%d|" % (L, beta)))
            if not rs:
                continue
            d = np.array([Bd["%d|%d|%g" % (L, beta, r)]["diff"] for r in rs])
            lo = d - np.array([Bd["%d|%d|%g" % (L, beta, r)]["diff_ci95"][0] for r in rs])
            hi = np.array([Bd["%d|%d|%g" % (L, beta, r)]["diff_ci95"][1] for r in rs]) - d
            ax.errorbar(np.array(rs) * (1 + 0.012 * (L // 16 - 2)), d, yerr=[lo, hi], fmt=mk, ms=3, lw=0.6, capsize=1.2, color=cols[L],
                        mfc=cols[L] if beta == 20 else "none")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylim(-0.03, 0.03)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$U_{\rm intact} - U_{\rm damaged}$")
    ax.set_title(r"Difference between the two starts (filled: $\beta$ = 20, open: 5)", loc="left", fontsize=7.5)
    tag(ax, "B")
    ax = axs[1, 0]
    for dwell, c in (("25", "#c0392b"), ("200", "#e08e0b"), ("1600", "#1a7a4c")):
        d = Rp["dwells"].get(dwell)
        if not d:
            continue
        ax.plot(d["rho"], d["down"], "-", color=c, lw=1.1, label="%s steps per $\\rho$" % dwell)
        ax.plot(d["rho"], d["up"], "--", color=c, lw=1.1)
    ax.axvspan(wins["continuous"][0], wins["continuous"][1], color="#fdf1c9", zorder=0)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$U$")
    ax.set_title("Quasi-static ramps, L = 48 (solid: down, dashed: up)", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    tag(ax, "C")
    ax = axs[1, 1]
    dw = sorted(Rp["dwells"], key=float)
    x = [float(d) for d in dw]
    y = np.array([Rp["dwells"][d]["area"][0] for d in dw])
    lo = y - np.array([Rp["dwells"][d]["area_ci95"][0] for d in dw])
    hi = np.array([Rp["dwells"][d]["area_ci95"][1] for d in dw]) - y
    ax.errorbar(x, y, yerr=[lo, hi], fmt="o-", ms=3.6, lw=0.9, capsize=2, color="#08519c")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("steps per $\\rho$ value")
    ax.set_ylabel(r"signed loop area $\int (U_\downarrow - U_\uparrow)\,d\rho$")
    ax.set_title("Loop area against sweep rate (95% CI, 8 runs)", loc="left", fontsize=8)
    tag(ax, "D")
    save(fig, "Fig2_hysteresis_tests.png")


def fig6():
    """Is rho enough? Dependence on the damage rate relative to phase relaxation, Omega = lam0/relax."""
    Rr = R()
    C = Rr["timescale"]["cells"]
    fig, axs = plt.subplots(2, 3, figsize=(11.2, 6.4))
    fig.subplots_adjust(hspace=0.42, wspace=0.34)
    for k, rho in enumerate((5.0, 2.5)):
        ax = axs[0, k]
        for mode in ("reset0", "keep", "neighbor"):
            cells = [(v["omega"], v["U"][0], v["U"][1], v["saturated"]) for key, v in C.items() if key.startswith("%s|%g|" % (mode, rho))]
            ok = [c for c in cells if not c[3]]
            sat = [c for c in cells if c[3]]
            ax.errorbar([c[0] for c in ok], [c[1] for c in ok], yerr=[1.96 * c[2] for c in ok], fmt="o", ms=3.4, color=MODECOL[mode], capsize=0,
                        label=mode)
            ax.plot([c[0] for c in sat], [c[1] for c in sat], "o", ms=3.4, mfc="none", color=MODECOL[mode])
        ax.set_xscale("log")
        ax.set_xlabel(r"$\Omega = \lambda_0/{\rm relax}$")
        ax.set_ylabel(r"$U$")
        ax.set_title(r"$\rho$ = %g, $h_0$ = 0.4 (open: $\lambda_0$ = 0.08)" % rho, loc="left", fontsize=8)
        if k == 0:
            ax.legend(frameon=False, fontsize=6.5, title="repair rule", title_fontsize=6.5)
        tag(ax, "AB"[k])
    ax = axs[0, 2]
    for mode in ("reset0", "keep", "neighbor"):
        pts = [(v["q"][0], v["U"][0], v["omega"]) for key, v in C.items() if key.startswith(mode + "|") and not v["saturated"]]
        sc = ax.scatter([p_[0] for p_ in pts], [p_[1] for p_ in pts], s=10, c=MODECOL[mode], alpha=0.8, label=mode)
    ax.set_xlabel(r"damaged fraction $q$")
    ax.set_ylabel(r"$U$")
    ax.set_title(r"$U$ against $q$ (all $\lambda_0 \leq$ 0.04, relax, $\rho$)", loc="left", fontsize=8)
    tag(ax, "C")
    ax = axs[1, 0]
    for rho, mk in ((5.0, "o"), (2.5, "s")):
        xs, ys = [], []
        for lam in (0.005, 0.01, 0.02, 0.04):
            for rel in (1, 2, 4, 8):
                us = [C["%s|%g|%g|%d" % (m, rho, lam, rel)]["U"][0] for m in ("reset0", "keep", "neighbor")]
                xs.append(lam / rel)
                ys.append(max(us) - min(us))
        ax.plot(xs, ys, mk, ms=4, color="#08519c", mfc="#08519c" if rho == 5 else "none", label=r"$\rho$ = %g" % rho)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Omega = \lambda_0/{\rm relax}$")
    ax.set_ylabel("spread of $U$ across the three repair rules")
    ax.set_title("The repair rule stops mattering as $\\Omega \\to 0$", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    tag(ax, "D")
    S_ = Rr["timescale_s"]
    ax = axs[1, 1]
    lamc = {0.005: "#9ecae1", 0.02: "#3182bd", 0.08: "#08306b"}
    for lam, c in lamc.items():
        ks = sorted([k for k in S_ if abs(float(k.split("|")[0]) - lam) < 1e-9], key=lambda k: S_[k]["omega"])
        om = [S_[k]["omega"] for k in ks]
        ax.errorbar(om, [S_[k]["ups_tw"][0] for k in ks], yerr=[1.96 * S_[k]["ups_tw"][1] for k in ks], fmt="o-", ms=4, lw=0.8, color=c, capsize=1.5,
                    label=r"$\lambda_0$ = %g" % lam)
        ax.plot(om, [S_[k]["ups_eq"][0] for k in ks], "s--", ms=4, lw=0.8, color=c, mfc="none")
    ax.plot([], [], "ko", ms=4, label=r"$\Upsilon_{\rm tw}$ (filled)")
    ax.plot([], [], "ks", ms=4, mfc="none", label=r"$\Upsilon_{\rm eq}$ (open)")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Omega$ (relax = 1, 2, 8 along each line)")
    ax.set_ylabel("stiffness")
    ax.set_title(r"Field-free, keep, $\rho$ = 10, $L$ = 32", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6, ncol=2, loc="lower right")
    tag(ax, "E")
    ax = axs[1, 2]
    for key, v in S_.items():
        lam, rel = key.split("|")
        if float(lam) < 0.07:
            ax.plot(v["omega"], v["ups_eq"][0] / v["ups_tw"][0], "o", ms=4.5, color="#c0392b")
    ax.axhline(1, color="k", lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Omega$")
    ax.set_ylabel(r"$\Upsilon_{\rm eq}/\Upsilon_{\rm tw}$")
    ax.set_title(r"The equilibrium formula converges as $\Omega \to 0$ ($\lambda_0 \leq$ 0.02)", loc="left", fontsize=7.5)
    tag(ax, "F")
    save(fig, "Fig6_timescale.png")


def fig7():
    """Quenched dilution: helicity-modulus crossings, extrapolation and clustered damage."""
    Rr = R()
    F = Rr["fss"]["0.35"]
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.6))
    fig.subplots_adjust(wspace=0.34)
    NK = F["NK"]
    cols = {"16": "#c6dbef", "24": "#9ecae1", "32": "#6baed6", "48": "#3182bd", "64": "#08519c", "96": "#08306b", "128": "#000000"}
    ax = axs[0]
    for L, v in F["per_L"].items():
        ax.errorbar(v["f"], v["Y"], yerr=1.96 * np.array(v["Y_se"]), fmt="o-", ms=2.4, lw=0.7, capsize=1, color=cols[L], label="L = %s" % L)
    ax.axhline(NK, color="k", ls="--", lw=0.8)
    ax.set_xlim(0.1, 0.45)
    ax.set_xlabel("removed fraction f")
    ax.set_ylabel(r"helicity modulus $\Upsilon$")
    ax.set_title("T = 0.35, random dilution (95% CI)", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=5.8, ncol=2)
    tag(ax, "A")
    ax = axs[1]
    Ls = [int(L) for L in F["per_L"] if F["per_L"][L]["f_KT"] == F["per_L"][L]["f_KT"] and F["per_L"][L]["ci95"]]
    X = 1 / np.log(np.array(Ls, float)) ** 2
    y = np.array([F["per_L"][str(L)]["f_KT"] for L in Ls])
    lo = y - np.array([F["per_L"][str(L)]["ci95"][0] for L in Ls])
    hi = np.array([F["per_L"][str(L)]["ci95"][1] for L in Ls]) - y
    ax.errorbar(X, y, yerr=[lo, hi], fmt="o", ms=4, capsize=2.5, color="#08519c", label="crossing and its 95% CI")
    e = F["extrap_all"]
    xx = np.linspace(0, X.max() * 1.05, 20)
    ax.plot(xx, e["f_inf"] + e["slope"] * xx, "-", color="#08519c", lw=0.9)
    ax.errorbar([0], [e["f_inf"]], yerr=[[e["f_inf"] - e["ci95"][0]], [e["ci95"][1] - e["f_inf"]]], fmt="s", ms=5, color="#c0392b", capsize=3,
                label=r"$L \to \infty$ (all sizes)")
    for k, (Lmin, mk) in enumerate((("24", "^"), ("32", "v"), ("48", "d"))):
        se = F["sensitivity_lnL2"].get(Lmin)
        if se:
            ax.errorbar([0.003 * (k + 1)], [se["f_inf"]], yerr=[[se["f_inf"] - se["ci95"][0]], [se["ci95"][1] - se["f_inf"]]], fmt=mk, ms=4, capsize=2,
                        color="#e08e0b", label=r"$L \geq$ 24, 32, 48 (triangles)" if k == 0 else None)
    ax.axhspan(e["mean_over_sizes"] - e["scatter_sd_over_sizes"], e["mean_over_sizes"] + e["scatter_sd_over_sizes"], color="#ddd", alpha=0.5, zorder=0)
    ax.set_xlabel(r"$1/\ln^2 L$")
    ax.set_ylabel(r"$f_{KT}(L)$")
    ax.set_title("Crossings and extrapolation (grey: scatter over sizes)", loc="left", fontsize=7.5)
    ax.legend(frameon=False, fontsize=6, loc="lower left")
    tag(ax, "B")
    ax = axs[2]
    Cl = Rr.get("cluster", {})
    for rad, c in ((0, "#08519c"), (1, "#1a7a4c"), (2, "#e08e0b"), (3, "#c0392b")):
        v = Cl.get("%d|48" % rad)
        if v:
            ax.errorbar(v["f"], v["Y"], yerr=1.96 * np.array(v["Y_se"]), fmt="o-", ms=2.6, lw=0.8, capsize=1, color=c, label="cluster radius %d" % rad)
    ax.axhline(NK, color="k", ls="--", lw=0.8)
    ax.set_xlabel("removed fraction f")
    ax.set_ylabel(r"$\Upsilon$")
    ax.set_title("Random against clustered damage, L = 48", loc="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    tag(ax, "C")
    save(fig, "Fig7_dilution.png")


def fig8():
    """From damaged fraction to displacement damage, and annealing."""
    Rr = R()
    M = Rr["mapping"]
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.6))
    fig.subplots_adjust(wspace=0.32)
    ax = axs[0]
    ax.plot(M["tau_f"]["f"], M["tau_f"]["tau"], "-", color="#08519c", lw=1.2, label=r"$T_{BKT}(f)/T_{BKT}(0)$ (interpolated)")
    ax.plot(M["points_f"], np.array(M["points_T"]) / M["T_BKT_clean"], "o", color="#08519c", ms=4, label="helicity-modulus crossings")
    ax.axhline(1 - M["dTc_over_Tc"], color="#c0392b", ls=":", lw=0.9)
    ax.plot([M["f_at_dpa"]], [1 - M["dTc_over_Tc"]], "s", color="#c0392b", ms=5, label="reported ~10% decrease of $T_c$ at 3.7-3.8 mdpa")
    ax.set_xlabel("removed fraction f")
    ax.set_ylabel(r"$T_c/T_{c0}$ if $T_c \propto T_{BKT}$")
    sec = ax.secondary_xaxis("top", functions=(lambda f: f / M["sites_per_dpa"] * 1e3, lambda d: d * M["sites_per_dpa"] * 1e-3))
    sec.set_xlabel("displacement damage (mdpa) at the implied %.0f sites per dpa" % M["sites_per_dpa"], fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc="lower left")
    tag(ax, "A")
    ax = axs[1]
    for (nu, c, ls, ex) in (("1e+11", "#6baed6", "-", 11), ("1e+13", "#08519c", "--", 13), ("1e+15", "#08306b", ":", 15)):
        a = M["anneal"].get(nu) or M["anneal"].get(nu.replace("+", ""))
        if a:
            ax.plot(a["T_C"], np.array(a["recovery"]) * 100, ls, color=c, lw=1.4,
                    label=r"$\nu_a$ = $10^{%d}$ s$^{-1}$: $E_0$ = %.2f eV, $\sigma_E$ = %.2f eV" % (ex, a["E0_eV"], a["sigma_eV"]))
    ax.plot([150, 400], [25, 60], "o", color="#c0392b", ms=6, label="reported recovery (the two fitted points)")
    ax.set_xlabel(r"annealing temperature ($^\circ$C), 12 h")
    ax.set_ylabel(r"recovered fraction of the $T_c$ decrease (%)")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=5.8, loc="upper left")
    tag(ax, "B")
    save(fig, "Fig8_mapping_annealing.png")


FIGS = {"fig1": fig1, "fig2": fig2, "fig3": fig3, "fig4": fig4, "fig5": fig5, "fig6": fig6, "fig7": fig7, "fig8": fig8}

if __name__ == "__main__":
    for f in sys.argv[1:]:
        FIGS[f]()
