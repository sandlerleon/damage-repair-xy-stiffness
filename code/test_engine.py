# -*- coding: utf-8 -*-
"""Checks on the engine before any result is computed with it.

  1. Same seed, same parameters: xy_engine reproduces the MyUncle engine (the code used for the superconducting preset)
     bit for bit, with the default options (repair = "reset0", twist = 0).
  2. beta = 0: the damaged fraction relaxes to the discrete two-state fixed point 1 - R/(lam + R - lam R), which tends
     to the continuous-time value 1/(1 + rho) as lam, R -> 0.
  3. R = 0 never repairs, lam0 = 0 never damages.
  4. Repair options: "keep" leaves the phase of a repaired site unchanged; "reset0" sets it to 0; "neighbor" sets it to
     the phase of the neighbouring phasors.
  5. The phase-stiffness estimators: on a clean lattice (no damage) the twist response equals the equilibrium
     fluctuation formula within statistical error, and the energy and the disorder of a twisted ordered lattice are
     gauge-consistent (an ordered lattice with the twisted ground state has D = 0 and zero current).
  6. At h0 = 0 with a symmetry-preserving repair the magnetisation phase is uniformly distributed (no preferred phase).

    python test_engine.py          (MyUncle is looked for at ../../_repos/MyUncle or on PYTHONPATH)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xy_engine import Lattice, Params  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, ("  -- " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


# 1. bit-for-bit agreement with MyUncle
for cand in (os.path.join(HERE, "..", "..", "_repos", "MyUncle"), os.path.join(HERE, "..", "..", "myuncle-repo")):
    if os.path.exists(os.path.join(cand, "MyUncle_core.py")):
        sys.path.insert(0, os.path.abspath(cand))
        break
try:
    from MyUncle_core import MC, Config
    have = True
except ImportError:
    have = False
if have:
    for beta, rho, h0 in ((20.0, 5.0, 0.4), (0.0, 2.5, 0.4), (20.0, 2.5, 0.0)):
        cfg = Config(L=24, seed=7, T=0.35, h0=h0, g0=0.12, lam0=0.02, beta=beta, hill_h=6.0, hill_K=0.5)
        cfg = cfg.copy(R=rho * cfg.lam_eff)
        ref = MC(cfg)
        p = Params(L=24, T=0.35, h0=h0, g0=0.12, lam0=0.02, beta=beta).with_rho(rho)
        new = Lattice(p, seed=7)
        for _ in range(60):
            a = ref.step()
            b = new.step()
        same = np.array_equal(ref.th, new.th) and np.array_equal(ref.g, new.g) and a == b
        check("reproduces MyUncle exactly after 60 steps (beta=%g, rho=%g, h0=%g)" % (beta, rho, h0), same)
else:
    print("  [SKIP] MyUncle_core not found; bit-for-bit comparison not run")

# 2. beta = 0 two-state chain
p = Params(L=48, beta=0.0).with_rho(3.0)
lat = Lattice(p, seed=11)
dq = []
for t in range(1500):
    lat.step()
    if t >= 500:
        dq.append(lat.damaged_fraction())
lam = p.lam_eff
target = 1.0 - p.R / (lam + p.R - lam * p.R)
se = np.std(dq) / np.sqrt(len(dq) / 50.0)
check("beta = 0: damaged fraction = 1 - R/(lam + R - lam R) of the discrete two-state chain",
      abs(np.mean(dq) - target) < max(4 * se, 0.003), "%.4f vs %.4f (continuous time %.4f)" % (np.mean(dq), target, 1 / (1 + p.rho)))

# 3. trivial rates
p = Params(L=16).with_rho(0.0)
lat = Lattice(p, seed=3)
lat.set_damaged(0.5)
n0 = int((~lat.intact).sum())
for _ in range(20):
    assert lat.step()[1] == 0
check("R = 0: no site is ever repaired", int((~lat.intact).sum()) >= n0)
lat = Lattice(Params(L=16, lam0=0.0, R=0.5), seed=3)
check("lam0 = 0: no site is ever damaged", sum(lat.step()[0] for _ in range(20)) == 0)

# 4. repair options act on the phase as documented
for mode in ("reset0", "keep", "neighbor"):
    p = Params(L=16, repair=mode, lam0=0.0, R=1.0, beta=0.0)
    lat = Lattice(p, seed=5)
    lat.set_damaged(0.3, random_angles=True)
    before = lat.th.copy()
    rep_mask = ~lat.intact
    nb = None
    if mode == "neighbor":
        z = np.exp(1j * lat.th)
        g = lat.g
        # g after repair: the repaired sites get weight 1; neighbours that are themselves repaired also count 1
        gg = np.ones_like(g)
        nb = np.angle(sum(np.roll(gg * z, s, a) for s, a in ((1, 0), (-1, 0), (1, 1), (-1, 1))))
    lat.site_update()
    if mode == "reset0":
        ok = np.all(lat.th[rep_mask] == 0.0)
    elif mode == "keep":
        ok = np.array_equal(lat.th, before)
    else:
        d = np.angle(np.exp(1j * (lat.th[rep_mask] - nb[rep_mask])))
        ok = np.max(np.abs(d)) < 1e-9
    check("repair = %s acts on the phase as documented" % mode, bool(ok))

# 5. phase-stiffness estimators
q = Params(L=16, T=0.5, h0=0.0, lam0=0.0, R=0.0)        # clean lattice, no damage ever
phi0 = 0.15
ups_eq, ups_tw = [], []
for s in range(8):
    res = {}
    for sign in (0, +1, -1):
        pp = q.replace(twist=sign * phi0)
        lat = Lattice(pp, seed=100 + s)
        for _ in range(400):
            lat.metropolis()
        cs, ss = [], []
        for _ in range(6000):
            lat.metropolis()
            c, sx, _, _ = lat.bond_terms()
            cs.append(c)
            ss.append(sx)
        res[sign] = (np.mean(cs), np.var(ss), np.mean(ss))
    N = 16 * 16
    ups_eq.append(res[0][0] - N / q.T * res[0][1])
    ups_tw.append(-(res[+1][2] - res[-1][2]) / (2 * phi0))
d = np.mean(ups_eq) - np.mean(ups_tw)
err = np.hypot(np.std(ups_eq, ddof=1), np.std(ups_tw, ddof=1)) / np.sqrt(8)
check("clean lattice: twist-response stiffness = equilibrium fluctuation formula", abs(d) < 4 * err + 0.01,
      "%.3f vs %.3f +- %.3f" % (np.mean(ups_tw), np.mean(ups_eq), err))
# twisted ordered lattice: energy-minimal state, D = 0, zero current
p = Params(L=16, T=0.5, h0=0.0, lam0=0.0, R=0.0, twist=0.2)
lat = Lattice(p, seed=1)
lat.th = 0.2 * np.tile(np.arange(16), (16, 1)).astype(float)       # winding: theta = phi x (only consistent mod 2 pi on a ring)
D = lat.local_disorder()
check("twisted ordered state has zero local disorder in the twist gauge (interior sites)", np.max(D[:, 1:-1]) < 1e-12)
_, sx, _, sy = lat.bond_terms()
check("... and zero current away from the seam", abs(sy) < 1e-12)

# 6. symmetry at h0 = 0 with a symmetry-preserving repair
for mode in ("keep", "neighbor"):
    p = Params(L=24, h0=0.0, repair=mode).with_rho(5.0)
    ph = []
    for s in range(40):
        lat = Lattice(p, seed=900 + s)
        lat.th = lat.rng.uniform(-np.pi, np.pi, lat.th.shape)       # random global phase for the start
        for _ in range(150):
            lat.step()
        ph.append(np.angle(np.mean(np.exp(1j * lat.th))))
    r = abs(np.mean(np.exp(1j * np.array(ph))))
    check("h0 = 0, repair = %s: global phase has no preferred value (|<e^{i psi}>| = %.2f over 40 runs)" % (mode, r), r < 0.35)
lat = Lattice(Params(L=24, h0=0.0, repair="reset0").with_rho(5.0), seed=9)
ph = []
for s in range(40):
    lat = Lattice(Params(L=24, h0=0.0, repair="reset0").with_rho(5.0), seed=900 + s)
    lat.th = lat.rng.uniform(-np.pi, np.pi, lat.th.shape)
    for _ in range(150):
        lat.step()
    ph.append(np.angle(np.mean(np.exp(1j * lat.th))))
r0 = abs(np.mean(np.exp(1j * np.array(ph))))
check("h0 = 0, repair = reset0: the repair itself selects the phase 0 (|<e^{i psi}>| = %.2f)" % r0, r0 > 0.8)

print("\n%s" % ("ALL CHECKS PASSED" if not FAILS else "%d FAILED: %s" % (len(FAILS), ", ".join(FAILS))))
sys.exit(1 if FAILS else 0)
