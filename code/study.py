# -*- coding: utf-8 -*-
"""Numerical programme for "Damage-repair kinetics and phase stiffness in a site-diluted XY model".

Studies (each writes results/<study>_raw.json; every run is seeded from SeedSequence([MASTER, study code, ...]), so any
single run can be repeated alone):

  steady        steady-state order U, damaged fraction, local disorder and kill multiplier versus rho, h0 = 0.4 and the
                original repair (theta -> 0), beta = 0, 5, 20, L = 32 and 64, 10 seeds, 2000 + 4000 steps; plus the
                original short protocol (200 + 200 steps, 3 seeds) for comparison
  converge      U(t) from an intact and from a 95%-damaged start, to show how long steady state takes
  branches      two-start tests at fixed rho (the original repair, h0 = 0.4), beta = 5 and 20, L = 32, 48, 64
  ramps         quasi-static rho ramps down and up at seven dwell times, L = 48 (loop area against sweep rate)
  fieldfree     h0 = 0, two starts, symmetric repairs ("keep", "neighbor") and the original "reset0": U, damaged fraction,
                C(r) and the equilibrium-formula stiffness, versus rho
  stiffness     the same with a twist +-phi0: the twist-response stiffness, versus rho, L = 48
  stiff_L       stiffness against lattice size L = 24-96 in the ordered regime
  timescale     U against (rho, lam0, relax) for the three repairs: is rho enough? (h0 = 0.4, L = 32)
  timescale_s   stiffness against lam0 and relax at fixed rho (h0 = 0)
  sensitivity   one parameter varied at a time at rho = 2.5 and 5 (T, g0, h0, beta, L up to 128)
  tauphi        the relaxation time of the phase field in a frozen damage landscape
  fss           quenched site dilution, helicity modulus against f, L = 16-128, T = 0.35 (and 0.25)
  cluster       the same with compact damage clusters of radius 0-3
  validation    control simulations of the twist response (clean lattice, imposed windings, static dilution against an independent
                equilibrium code); validation_long repeats the diluted-lattice comparison with 4x longer runs and 16 maps
  tcurve        the BKT line T_BKT(f): helicity modulus against temperature at fixed dilution (maps damage to Tc)

    python study.py <study> [<study> ...] [--quick]
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xy_engine import Lattice, Params, Quenched  # noqa: E402

RES = os.path.abspath(os.path.join(HERE, "..", "results"))
os.makedirs(RES, exist_ok=True)
MASTER = 20261005
QUICK = "--quick" in sys.argv
NPROC = max(1, (os.cpu_count() or 2) - 1)
CODES = {"steady": 1, "converge": 2, "branches": 3, "ramps": 4, "fieldfree": 5, "stiffness": 6, "stiff_L": 7,
         "timescale": 8, "timescale_s": 9, "sensitivity": 10, "tauphi": 11, "fss": 12, "cluster": 13, "tcurve": 14, "closure0": 15, "stiffness_dmg": 16, "stiff_L2": 17, "closure0b": 18, "stiff_L3": 19, "validation": 20, "validation_long": 21}


def phi0_of(L):
    """Twist amplitude: well below pi/L, beyond which a periodic lattice changes its winding number (phase slip)."""
    return round(0.3 * np.pi / L, 5)


def seed_of(*code):
    return int(np.random.SeedSequence([MASTER] + [int(round(c)) for c in code]).generate_state(1)[0])


def tau_int(x, c=6.0):
    """Sokal automatic-window integrated autocorrelation time, in samples."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if n < 20 or np.var(x) == 0:
        return 0.5
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    ac = ac / ac[0]
    tau = 0.5
    for W in range(1, n):
        tau += ac[W]
        if W >= c * tau:
            return float(tau)
    return float(tau)


def summ(x, every=1):
    x = np.asarray(x, float)
    return {"mean": float(np.nanmean(x)), "sd": float(np.nanstd(x)), "tau": tau_int(x) * every}


def rho_code(r):
    return int(round(r * 100))


# ======================================================================= generic steady-state run
def nss(pd, seed, start="intact", eq=2000, meas=4000, every=5, corr=False, stiff=False, series=False):
    p = Params(**pd)
    lat = Lattice(p, seed)
    if start == "damaged":
        lat.set_damaged(0.95, True)
    for _ in range(eq):
        lat.step()
    N = p.L * p.L
    U, Ui, q, D, mu, muref, cx, sx, cy, sy, wm, wsd = ([] for _ in range(12))
    Call = Cint = None
    nc = 0
    kills = reps = 0
    sum_mu = sum_q = sum_mr = 0.0
    for t in range(meas):
        k, r = lat.step()
        kills += k
        reps += r
        sum_mu += lat.last_mu
        sum_q += lat.last_q
        sum_mr += lat.last_mu_ref
        if t % every == 0:
            U.append(lat.order())
            Ui.append(lat.order_intact())
            q.append(lat.damaged_fraction())
            d, m, mr = lat.disorder_stats(ref_beta=20.0)
            D.append(d)
            mu.append(m)
            muref.append(mr)
            if stiff:
                a, b, c, d2 = lat.bond_terms()
                cx.append(a), sx.append(b), cy.append(c), sy.append(d2)
                w1, w2 = lat.winding_x()
                wm.append(w1), wsd.append(w2)
            if corr:
                ca, ci = lat.correlation()
                Call = ca if Call is None else Call + ca
                Cint = ci if Cint is None else Cint + ci
                nc += 1
    out = {"U": summ(U, every), "Ui": summ(Ui, every), "q": summ(q, every), "D": summ(D, every), "mu": summ(mu, every), "mu_ref20": summ(muref, every),
           "J_kill": kills / (meas * N), "J_rep": reps / (meas * N),
           "mu_kill": sum_mu / meas, "q_pre": sum_q / meas, "mu_ref_kill": sum_mr / meas}
    if stiff:
        out["stiff"] = {"c_x": float(np.mean(cx)), "c_y": float(np.mean(cy)), "s_x": float(np.mean(sx)), "s_y": float(np.mean(sy)),
                        "var_s_x": float(np.var(sx)), "var_s_y": float(np.var(sy)), "tau_s": tau_int(sx) * every,
                        "n": len(sx), "w_mean": float(np.mean(wm)), "w_rowstd": float(np.mean(wsd)),
                        "w_end": float(wm[-1])}
    if corr:
        out["C_all"] = (Call / nc).round(5).tolist()
        out["C_int"] = (Cint / nc).round(5).tolist()
    if series:
        out["series"] = {"U": [round(u, 4) for u in U], "q": [round(u, 4) for u in q], "every": every}
    return out


def base(**kw):
    d = dict(L=32, T=0.35, h0=0.4, g0=0.12, lam0=0.02, beta=20.0, hill_h=6.0, hill_K=0.5, relax=2, repair="reset0", twist=0.0)
    d.update(kw)
    return d


def with_rho(pd, rho):
    d = dict(pd)
    d["R"] = rho * d["lam0"]
    return d


def task_nss(args):
    study, tag, pd, seed_code, kw = args
    rec = nss(pd, seed_of(CODES[study], *seed_code), **kw)
    return {"tag": tag, "params": pd, "seed": seed_code, "kw": {k: v for k, v in kw.items() if k != "series"}, **rec}


def run_tasks(study, tasks, worker=task_nss, cost=None):
    cost = cost or (lambda t: -(t[2]["L"] ** 2) * (t[2]["relax"] if "relax" in t[2] else 1))
    tasks = sorted(tasks, key=cost)
    out = []
    t0 = time.time()
    print("%s: %d tasks" % (study, len(tasks)), flush=True)
    with Pool(NPROC) as pool:
        for i, r in enumerate(pool.imap_unordered(worker, tasks, chunksize=1)):
            out.append(r)
            if (i + 1) % 100 == 0:
                print("  %s %d/%d  (%.0f s)" % (study, i + 1, len(tasks), time.time() - t0), flush=True)
    return out, time.time() - t0


def save(study, runs, meta, secs):
    with open(os.path.join(RES, "%s_raw.json" % study), "w") as f:
        json.dump({"study": study, "meta": meta, "runs": runs}, f)
    print("%s done in %.0f s (%d runs)" % (study, secs, len(runs)), flush=True)


def S(n):          # seeds
    return 3 if QUICK else n


RHOS_MAIN = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0, 40.0]
RHOS_FF = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 40.0]


# ======================================================================= studies
def study_steady():
    tasks = []
    for L in ([32] if QUICK else [32, 64]):
        for beta in (0.0, 5.0, 20.0):
            for rho in RHOS_MAIN:
                for s in range(S(10 if L == 32 else 6)):
                    tasks.append(("steady", "long", with_rho(base(L=L, beta=beta), rho), (L, int(beta), rho_code(rho), s),
                                  dict(eq=2000, meas=4000, every=5)))
    for beta in (0.0, 5.0, 20.0):                       # the original protocol, for comparison
        for rho in (1.0, 2.5, 5.0, 10.0, 40.0):
            for s in range(3):
                tasks.append(("steady", "short", with_rho(base(L=32, beta=beta), rho), (32, int(beta), rho_code(rho), 100 + s),
                              dict(eq=200, meas=200, every=2)))
    runs, t = run_tasks("steady", tasks)
    save("steady", runs, {"rhos": RHOS_MAIN}, t)


def study_converge():
    tasks = []
    for rho in (2.5, 5.0):
        for start in ("intact", "damaged"):
            for s in range(S(8)):
                tasks.append(("converge", "c", with_rho(base(L=32, beta=20.0), rho), (rho_code(rho), start == "damaged", s),
                              dict(start=start, eq=0, meas=8000, every=50, series=True)))
    runs, t = run_tasks("converge", tasks)
    save("converge", runs, {}, t)


def study_branches():
    tasks = []
    rhos = [2.0, 2.5, 3.0, 3.25, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
    for L, ns, betas, rr in ((32, 10, (5.0, 20.0), rhos), (48, 8, (5.0, 20.0), rhos), (64, 6, (20.0,), [3.0, 4.0, 5.0, 6.0, 8.0, 10.0])):
        if QUICK and L > 32:
            continue
        for beta in betas:
            for rho in rr:
                for start in ("intact", "damaged"):
                    for s in range(S(ns)):
                        tasks.append(("branches", "b", with_rho(base(L=L, beta=beta), rho),
                                      (L, int(beta), rho_code(rho), start == "damaged", s),
                                      dict(start=start, eq=3000, meas=3000, every=5)))
    runs, t = run_tasks("branches", tasks)
    save("branches", runs, {"rhos": rhos}, t)


def task_ramp(args):
    pd, seed_code, dwell = args
    p = Params(**with_rho(pd, 16.0))
    lat = Lattice(p, seed_of(CODES["ramps"], *seed_code))
    for _ in range(200):
        lat.step()
    grid = np.linspace(16.0, 1.0, 31)
    nav = max(5, dwell // 5)

    def sweep(rhos):
        out = []
        for r in rhos:
            lat.p = Params(**with_rho(pd, float(r)))
            u = []
            for t in range(dwell):
                lat.step()
                if t >= dwell - nav:
                    u.append(lat.order())
            out.append(float(np.mean(u)))
        return out
    down = sweep(grid)
    up = sweep(grid[::-1])
    return {"tag": "ramp", "params": pd, "seed": seed_code, "dwell": dwell, "rho_grid": grid.tolist(), "down": down, "up": up[::-1]}


def study_ramps():
    dwells = [25, 50, 100] if QUICK else [25, 50, 100, 200, 400, 800, 1600]
    tasks = [(base(L=48 if not QUICK else 32, beta=20.0), (d, s), d) for d in dwells for s in range(S(8))]
    runs, t = run_tasks("ramps", tasks, worker=task_ramp, cost=lambda t: -t[2])
    save("ramps", runs, {"dwells": dwells}, t)


def study_fieldfree():
    tasks = []
    combos = [("keep", 20.0), ("neighbor", 20.0), ("reset0", 20.0), ("keep", 0.0)]
    plan = [(32, 10, combos, RHOS_FF), (48, 10, combos[:2], [3.0, 5.0, 7.0, 10.0, 15.0, 20.0])]
    for L, ns, cmb, rr in plan:
        if QUICK and L > 32:
            continue
        for mode, beta in cmb:
            for rho in rr:
                for start in ("ordered", "damaged"):
                    for s in range(S(ns)):
                        tasks.append(("fieldfree", "ff", with_rho(base(L=L, h0=0.0, beta=beta, repair=mode), rho),
                                      (L, int(beta), ["keep", "neighbor", "reset0"].index(mode), rho_code(rho), start == "damaged", s),
                                      dict(start="intact" if start == "ordered" else "damaged", eq=3000, meas=3000, every=10,
                                           corr=True, stiff=True)))
    runs, t = run_tasks("fieldfree", tasks)
    save("fieldfree", runs, {"rhos": RHOS_FF}, t)


def task_twist(args):
    """Twist response: +-phi0 runs from the same seed; returns the mean currents."""
    study, tag, pd, seed_code, kw = args
    res = {}
    PHI0 = phi0_of(pd["L"])
    for sign in (+1, -1):
        q = dict(pd)
        q["twist"] = sign * PHI0
        res[sign] = nss(q, seed_of(CODES[study], *seed_code), **kw)
    return {"tag": tag, "params": pd, "seed": seed_code, "kw": kw, "phi0": PHI0,
            "plus": res[+1], "minus": res[-1]}


def study_stiffness():
    tasks = []
    for mode, beta in (("keep", 20.0), ("neighbor", 20.0), ("keep", 0.0)):
        for rho in RHOS_FF:
            for s in range(S(10)):
                tasks.append(("stiffness", "tw", with_rho(base(L=32, h0=0.0, beta=beta, repair=mode), rho),
                              (["keep", "neighbor"].index(mode), int(beta), rho_code(rho), False, s),
                              dict(start="intact", eq=3000, meas=6000, every=10, stiff=True)))
    runs, t = run_tasks("stiffness", tasks, worker=task_twist)
    save("stiffness", runs, {"rhos": RHOS_FF, "phi0": phi0_of(32)}, t)


def study_stiff_L2():
    """Stiffness against lattice size across the threshold: does it stay finite (a jump) or vanish as L grows?"""
    tasks = []
    plan = [("keep", [8.0, 9.0], [(32, 10), (48, 10), (64, 8), (96, 6)]), ("neighbor", [6.0, 7.0], [(32, 10), (48, 10), (64, 8)])]
    for mode, rhos, sizes in plan:
        for rho in rhos:
            for L, ns in sizes:
                if QUICK and L > 32:
                    continue
                for s in range(S(ns)):
                    tasks.append(("stiff_L2", "tw", with_rho(base(L=L, h0=0.0, beta=20.0, repair=mode), rho),
                                  (["keep", "neighbor"].index(mode), L, rho_code(rho), s), dict(start="intact", eq=3000, meas=4000, every=10, stiff=True)))
    runs, t = run_tasks("stiff_L2", tasks, worker=task_twist)
    save("stiff_L2", runs, {}, t)


def study_stiff_L3():
    """A finer scan across the two thresholds bracketed by stiff_L2: rho = 7.25, 7.5, 7.75 (keep) and 5.25, 5.5, 5.75 (neighbor)."""
    tasks = []
    plan = [("keep", [7.25, 7.5, 7.75], [(32, 10), (48, 10), (64, 8), (96, 6)]), ("neighbor", [5.25, 5.5, 5.75], [(32, 10), (48, 10), (64, 8)])]
    for mode, rhos, sizes in plan:
        for rho in rhos:
            for L, ns in sizes:
                if QUICK and L > 32:
                    continue
                for s in range(S(ns)):
                    tasks.append(("stiff_L2", "tw", with_rho(base(L=L, h0=0.0, beta=20.0, repair=mode), rho),
                                  (["keep", "neighbor"].index(mode), L, rho_code(rho), 200 + s), dict(start="intact", eq=3000, meas=4000, every=10, stiff=True)))
    runs, t = run_tasks("stiff_L3", tasks, worker=task_twist)
    save("stiff_L3", runs, {}, t)


def study_stiffness_dmg():
    """Twist response of states reached from a 95%-damaged start (random phases): vortices and windings can survive in
    them, so the winding number is recorded and the response is corrected for it (analysis)."""
    tasks = []
    for mode, rr in (("keep", [8.0, 10.0, 12.0, 15.0, 20.0]), ("neighbor", [7.0, 8.0, 10.0, 12.0, 15.0, 20.0])):
        for rho in rr:
            for s in range(S(10)):
                tasks.append(("stiffness", "tw", with_rho(base(L=32, h0=0.0, beta=20.0, repair=mode), rho),
                              (["keep", "neighbor"].index(mode), 20, rho_code(rho), True, 100 + s),
                              dict(start="damaged", eq=3000, meas=6000, every=10, stiff=True)))
    runs, t = run_tasks("stiffness_dmg", tasks, worker=task_twist)
    save("stiffness_dmg", runs, {}, t)


def study_stiff_L():
    tasks = []
    for L, ns in ((24, 10), (32, 10), (48, 10), (64, 8), (96, 6)):
        if QUICK and L > 32:
            continue
        for rho in (5.0, 7.0, 10.0, 20.0):
            for s in range(S(ns)):
                tasks.append(("stiff_L", "tw", with_rho(base(L=L, h0=0.0, beta=20.0, repair="keep"), rho),
                              (L, rho_code(rho), s), dict(start="intact", eq=3000, meas=4000, every=10, stiff=True)))
    runs, t = run_tasks("stiff_L", tasks, worker=task_twist)
    save("stiff_L", runs, {}, t)


LAMS = [0.005, 0.01, 0.02, 0.04, 0.08]
RELAXES = [1, 2, 4, 8]


def study_timescale():
    tasks = []
    for mode in ("reset0", "keep", "neighbor"):
        for rho in (2.5, 5.0):
            for lam in LAMS:
                for rel in RELAXES:
                    for s in range(S(8)):
                        tasks.append(("timescale", "ts", with_rho(base(L=32, beta=20.0, repair=mode, lam0=lam, relax=rel), rho),
                                      (["reset0", "keep", "neighbor"].index(mode), rho_code(rho), int(lam * 10000), rel, s),
                                      dict(eq=1500, meas=3000, every=5)))
    runs, t = run_tasks("timescale", tasks)
    save("timescale", runs, {"lams": LAMS, "relaxes": RELAXES}, t)


def study_timescale_s():
    tasks = []
    for rho in (10.0,):
        for lam in (0.005, 0.02, 0.08):
            for rel in (1, 2, 8):
                for s in range(S(10)):
                    tasks.append(("timescale_s", "tw", with_rho(base(L=32, h0=0.0, beta=20.0, repair="keep", lam0=lam, relax=rel), rho),
                                  (rho_code(rho), int(lam * 10000), rel, s), dict(eq=3000, meas=6000, every=10, stiff=True)))
    runs, t = run_tasks("timescale_s", tasks, worker=task_twist)
    save("timescale_s", runs, {}, t)


def study_sensitivity():
    tasks = []
    variants = [("T", [0.25, 0.35, 0.5, 0.7]), ("g0", [0.0, 0.03, 0.12, 0.3, 0.5]), ("h0", [0.0, 0.2, 0.4, 0.8]),
                ("beta", [0.0, 5.0, 10.0, 20.0, 30.0]), ("L", [24, 32, 48, 64, 96])]
    for name, vals in variants:
        for v in vals:
            for rho in (2.5, 5.0):
                ns = 10 if name != "L" or v <= 64 else 6
                for s in range(S(ns)):
                    pd = base(**{name: v})
                    tasks.append(("sensitivity", name, with_rho(pd, rho), (["T", "g0", "h0", "beta", "L"].index(name), int(v * 100), rho_code(rho), s),
                                  dict(eq=2000, meas=4000, every=5)))
    runs, t = run_tasks("sensitivity", tasks)
    save("sensitivity", runs, {"variants": variants}, t)


def task_tauphi(args):
    pd, seed_code, sweeps = args
    p = Params(**pd)
    lat = Lattice(p, seed_of(CODES["tauphi"], *seed_code))
    for _ in range(2000):
        lat.step()
    U = []
    for _ in range(sweeps):
        lat.metropolis()                  # frozen damage landscape: phases only
        U.append(lat.order())
    return {"tag": "tauphi", "params": pd, "seed": seed_code, "tau_sweeps": tau_int(U), "U": float(np.mean(U))}


def study_tauphi():
    tasks = []
    for h0 in (0.4, 0.0):
        for mode in ("reset0", "keep"):
            if h0 == 0.0 and mode == "reset0":
                continue
            for rho in (1.0, 2.5, 5.0, 10.0):
                for s in range(S(10)):
                    tasks.append((with_rho(base(L=32, beta=20.0, h0=h0, repair=mode), rho), (int(h0 * 10), ["reset0", "keep"].index(mode), rho_code(rho), s), 20000))
    runs, t = run_tasks("tauphi", tasks, worker=task_tauphi, cost=lambda t: 0)
    save("tauphi", runs, {}, t)


FSS_F = [0.10, 0.15] + [round(0.20 + 0.01 * i, 2) for i in range(0, 15)] + [0.38, 0.42, 0.50]


def task_fss(args):
    L, T, f, seed_code, eq, meas, clusters = args
    q = Quenched(L, T, f, seed_of(CODES["cluster" if clusters else "fss"], *seed_code), clusters=clusters)
    for _ in range(eq):
        q.sweep()
    Y, M = [], []
    for t in range(meas):
        q.sweep()
        if t % 4 == 0:
            Y.append(q.helicity())
            M.append(q.magnetization())
    M = np.array(M)
    return {"L": L, "T": T, "f": f, "seed": seed_code, "clusters": clusters, "Y": float(np.mean(Y)), "tau_Y": tau_int(Y) * 4,
            "M2": float(np.mean(M ** 2)), "M4": float(np.mean(M ** 4)), "M": float(np.mean(M))}


def study_fss():
    tasks = []
    sizes = [16, 24, 32] if QUICK else [16, 24, 32, 48, 64, 96, 128]
    for L in sizes:
        ns = {96: 8, 128: 6}.get(L, 10 if L >= 48 else 12)
        fs = FSS_F if L <= 64 else ([f for f in FSS_F if 0.24 <= f <= 0.32] if L == 96 else [0.26, 0.28, 0.30])
        for f in fs:
            for s in range(S(ns)):
                tasks.append((L, 0.35, f, (L, int(round(f * 1000)), 35, s), 2000, 4000, 0))
    for L in ([16, 24] if QUICK else [16, 24, 32, 48]):          # T = 0.25: the transition is at larger f
        for f in [round(0.34 + 0.02 * i, 2) for i in range(0, 11)]:
            for s in range(S(8)):
                tasks.append((L, 0.25, f, (L, int(round(f * 1000)), 25, s), 2000, 4000, 0))
    runs, t = run_tasks("fss", tasks, worker=task_fss, cost=lambda t: -t[0] ** 2)
    save("fss", runs, {"f": FSS_F}, t)


def study_cluster():
    tasks = []
    for radius in (0, 1, 2, 3):
        for L in ([24] if QUICK else [32, 48]):
            for f in [round(0.20 + 0.03 * i, 2) for i in range(0, 13)]:
                for s in range(S(10)):
                    tasks.append((L, 0.35, f, (radius, L, int(round(f * 1000)), s), 2000, 3000, radius))
    runs, t = run_tasks("cluster", tasks, worker=task_fss, cost=lambda t: -t[0] ** 2)
    save("cluster", runs, {}, t)


def study_closure0():
    """The kill multiplier a beta = 20 feedback would apply to lattices in which damage is independent of phase."""
    tasks = []
    for h0, mode, rr in ((0.4, "reset0", RHOS_MAIN), (0.0, "keep", RHOS_FF)):
        for beta in (0.0, 20.0):
            for rho in rr:
                for s in range(S(10)):
                    tasks.append(("closure0", "c0", with_rho(base(L=32, h0=h0, beta=beta, repair=mode), rho),
                                  (int(h0 * 10), int(beta), rho_code(rho), s), dict(eq=2000, meas=3000, every=5)))
    runs, t = run_tasks("closure0", tasks)
    save("closure0", runs, {}, t)


def study_closure0b():
    """The beta = 0 reference of study closure0 extended to small rho, so that it reaches the damaged fractions (up to 0.95) of the beta = 20 states."""
    tasks = []
    for rho in (0.05, 0.1, 0.2, 0.35, 0.5, 0.7):
        for s in range(S(10)):
            tasks.append(("closure0", "c0", with_rho(base(L=32, h0=0.4, beta=0.0, repair="reset0"), rho), (4, 0, rho_code(rho), s), dict(eq=2000, meas=3000, every=5)))
    runs, t = run_tasks("closure0b", tasks)
    save("closure0b", runs, {}, t)


def run_static(pd, seed, g=None, w0=0, eq=1000, meas=3000, every=5):
    """A lattice without damage or repair (lam0 = R = 0), optionally with a prescribed coupling map g and an imposed winding w0 along x:
    the bond terms needed for the twist response and the equilibrium formula, and the winding number."""
    p = Params(**pd)
    lat = Lattice(p, seed)
    if g is not None:
        lat.g = g.copy()
    if w0:
        lat.th = lat.th + 2 * np.pi * w0 * np.arange(p.L)[None, :] / p.L
    for _ in range(eq):
        lat.step()
    cx, sx, cy, sy, wm = [], [], [], [], []
    for t in range(meas):
        lat.step()
        if t % every == 0:
            a, b, c, d = lat.bond_terms()
            cx.append(a), sx.append(b), cy.append(c), sy.append(d)
            wm.append(lat.winding_x()[0])
    return {"c_x": float(np.mean(cx)), "c_y": float(np.mean(cy)), "s_x": float(np.mean(sx)), "s_y": float(np.mean(sy)),
            "var_s_x": float(np.var(sx)), "var_s_y": float(np.var(sy)), "w_mean": float(np.mean(wm)), "n": len(sx)}


def task_study(kind, extra):
    return "validation_long" if extra.get("long") else "validation"


def task_valid(args):
    kind, tag, pd, seed_code, phi0, extra = args
    sd = seed_of(CODES[task_study(kind, extra)], *seed_code)
    out = {"kind": kind, "tag": tag, "params": pd, "seed": seed_code, "phi0": phi0}
    L = pd["L"]
    if kind == "static_eq":                      # independent equilibrium code: Metropolis on the quenched lattice, helicity formula
        q = Quenched(L, pd["T"], extra["f"], sd)
        for _ in range(extra.get("eq", 1000)):
            q.sweep()
        hs = []
        for t in range(extra.get("meas", 3000)):
            q.sweep()
            if t % 5 == 0:
                hs.append(q.helicity())
        out["ups_eq_static"] = float(np.mean(hs))
        out["gmap_seed"] = sd
        return out
    g = None
    if kind == "twist_diluted":                  # same coupling map as the static_eq run of the same seed
        q = Quenched(L, pd["T"], extra["f"], seed_of(CODES[task_study(kind, extra)], *extra["static_code"]))
        g = q.g
    for sign in (+1, -1):
        w0 = extra.get("w0_plus", 0) if sign > 0 else extra.get("w0_minus", 0)
        q2 = dict(pd)
        q2["twist"] = sign * phi0
        out["plus" if sign > 0 else "minus"] = run_static(q2, sd, g=g, w0=w0, eq=extra.get("eq", 1000), meas=extra.get("meas", 3000))
    return out


def study_validation():
    """Control simulations for the twist response (clean lattice, imposed windings, static dilution against an independent equilibrium code)."""
    T = 0.35
    pdc = lambda L: base(L=L, T=T, h0=0.0, lam0=0.0, R=0.0, beta=0.0, repair="keep")
    tasks = []
    for L in (16, 32, 64):
        for c in (0.1, 0.3, 0.6, 1.0, 1.5):
            for s in range(S(8)):
                tasks.append(("twist_clean", "c%g" % c, pdc(L), (1, L, int(c * 10), s), round(c * np.pi / L, 5), {}))
    for f in (0.1, 0.2):
        for s in range(S(8)):
            code = (2, 32, int(f * 100), s)
            tasks.append(("static_eq", "f%g" % f, pdc(32), code, 0.0, {"f": f}))
            tasks.append(("twist_diluted", "f%g" % f, pdc(32), code, phi0_of(32), {"f": f, "static_code": code}))
    for s in range(S(8)):
        for tag, wp, wm in (("w00", 0, 0), ("w11", 1, 1), ("w10", 1, 0)):
            tasks.append(("twist_wound", tag, pdc(32), (3, 32, wp * 10 + wm, s), phi0_of(32), {"w0_plus": wp, "w0_minus": wm}))
    runs, t = run_tasks("validation", tasks, worker=task_valid, cost=lambda t: -(t[2]["L"] ** 2))
    save("validation", runs, {"T": T}, t)



def study_validation_long():
    """Longer runs (4x) of the diluted-lattice comparison, 16 maps, to test whether the small difference between the twist response and the
    static equilibrium formula is a sampling effect."""
    pdc = base(L=32, T=0.35, h0=0.0, lam0=0.0, R=0.0, beta=0.0, repair="keep")
    tasks = []
    for f in (0.1, 0.2):
        for s in range(S(16)):
            code = (2, 32, int(f * 100), s)
            ex = {"f": f, "eq": 3000, "meas": 12000, "long": True}
            tasks.append(("static_eq", "f%g" % f, pdc, code, 0.0, ex))
            tasks.append(("twist_diluted", "f%g" % f, pdc, code, phi0_of(32), dict(ex, static_code=code)))
    runs, t = run_tasks("validation_long", tasks, worker=task_valid)
    save("validation_long", runs, {"T": 0.35}, t)


def study_tcurve():
    """The BKT line T_BKT(f): helicity modulus against temperature at fixed quenched dilution f, L = 16-48."""
    tasks = []
    Ts = [0.40, 0.50, 0.60, 0.70, 0.80, 0.88]
    fs = [0.0, 0.05, 0.10, 0.15, 0.20]
    for L in ([16, 24] if QUICK else [16, 24, 32, 48]):
        for f in fs:
            for T in Ts:
                for s in range(S(8)):
                    tasks.append((L, T, f, (L, int(round(f * 1000)), int(round(T * 100)), s), 1500, 3000, 0))
    runs, t = run_tasks("tcurve", tasks, worker=task_fss, cost=lambda t: -t[0] ** 2)
    save("tcurve", runs, {"T": Ts, "f": fs}, t)


if __name__ == "__main__":
    what = [a for a in sys.argv[1:] if not a.startswith("--")]
    fns = {k: globals()["study_" + k] for k in CODES}
    for w in what:
        print("== %s (%d workers%s)" % (w, NPROC, ", quick" if QUICK else ""), flush=True)
        fns[w]()
