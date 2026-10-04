# -*- coding: utf-8 -*-
"""Kinetic Monte Carlo engine: site-diluted XY lattice with state-dependent damage and repair.

The update rules are those of the MyUncle framework (github.com/sandlerleon/MyUncle, doi:10.5281/zenodo.21223569)
used for the superconducting preset: with the default options this engine reproduces MyUncle trajectories bit for bit
for the same seed (test_engine.py). Three options are added, each switched off by default:

  repair   what a repaired site does to its phase
           "reset0"    theta -> 0 (MyUncle). The reference phase 0 is common to all sites, so this repair is itself a
                       symmetry-breaking field acting at rate R, even when h0 = 0.
           "keep"      the coupling is restored and the phase is left as it is; the dynamics relaxes it.
           "neighbor"  the coupling is restored and the phase is set to that of the (coupling-weighted) neighbouring
                       phasors, i.e. the repaired site rejoins the coherent region around it.
           "keep" and "neighbor" preserve the global U(1) symmetry at h0 = 0.
  twist    a uniform phase twist phi per bond in the x direction (vector potential), so that the stationary response of
           the bond current to the twist measures the phase stiffness without assuming detailed balance.
  observables for U, the damaged fraction, the local disorder and kill multiplier, the bond currents, the terms of the
           helicity modulus and the phase correlation function C(r).

One step = `relax` checkerboard Metropolis sweeps of the phases at fixed integrity, then one synchronous site update:
  kill   an intact site is damaged (g -> g0) with probability lam_i = lam0 [1 + beta D_i^h/(K^h + D_i^h)],
         D_i = 1 - |(1/4) sum_{j in nn(i)} exp(i theta_j)|  (neighbour phases measured in the gauge of the twist)
  repair a damaged site is repaired (g -> 1) with probability R, then its phase is treated as per `repair`
rho = R / lam0 is the repair-to-damage ratio.

Energy: H = -J0 sum_<ij> g_i g_j cos(theta_j - theta_i - phi_ij) - h0 sum_i g_i cos(theta_i), k_B = 1.
"""
from dataclasses import asdict, dataclass

import numpy as np


@dataclass
class Params:
    L: int = 32
    T: float = 0.35
    h0: float = 0.4
    g0: float = 0.12
    lam0: float = 0.02
    beta: float = 20.0
    hill_h: float = 6.0
    hill_K: float = 0.5
    R: float = 0.10
    relax: int = 2
    repair: str = "reset0"
    twist: float = 0.0
    theta_init_sigma: float = 0.10

    @property
    def lam_eff(self):
        return self.lam0

    @property
    def rho(self):
        return self.R / self.lam0

    def with_rho(self, rho):
        d = asdict(self)
        d["R"] = rho * self.lam0
        return Params(**d)

    def replace(self, **kw):
        d = asdict(self)
        d.update(kw)
        return Params(**d)


class Lattice:
    def __init__(self, p: Params, seed: int):
        self.p = p
        self.rng = np.random.default_rng(seed)
        L = p.L
        self.parity = np.add.outer(np.arange(L), np.arange(L)) % 2
        self.th = self.rng.normal(0.0, p.theta_init_sigma, (L, L))
        self.g = np.full((L, L), 1.0)
        self.intact = np.ones((L, L), bool)
        self.last_q, self.last_mu, self.last_mu_ref = 0.0, 1.0, 1.0

    def set_damaged(self, fraction=0.95, random_angles=True):
        L = self.p.L
        v = self.rng.random((L, L)) < fraction
        self.g[v] = self.p.g0
        self.intact[v] = False
        if random_angles:
            self.th = self.rng.uniform(-np.pi, np.pi, (L, L))
        return self

    # ------------------------------------------------------------ dynamics
    def metropolis(self):
        p, th, g = self.p, self.th, self.g
        L, phi = p.L, p.twist
        for col in (0, 1):
            m = self.parity == col
            prop = th + self.rng.uniform(-np.pi, np.pi, (L, L))
            e0 = np.zeros((L, L))
            e1 = np.zeros((L, L))
            for sh, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
                J = g * np.roll(g, sh, ax)
                tj = np.roll(th, sh, ax)
                off = sh * phi if (ax == 1 and phi != 0.0) else 0.0
                e0 -= J * np.cos(th - tj - off)
                e1 -= J * np.cos(prop - tj - off)
            e0 -= p.h0 * g * np.cos(th)
            e1 -= p.h0 * g * np.cos(prop)
            acc = (self.rng.random((L, L)) < np.exp(-(e1 - e0) / p.T)) & m
            th = np.where(acc, prop, th)
        self.th = th

    def _aligned_phasors(self):
        """Phasors of the four neighbours, in the gauge in which the twisted ground state is uniform."""
        z = np.exp(1j * self.th)
        right, left = np.roll(z, -1, 1), np.roll(z, 1, 1)
        if self.p.twist != 0.0:
            right = right * np.exp(-1j * self.p.twist)
            left = left * np.exp(1j * self.p.twist)
        return right, left, np.roll(z, -1, 0), np.roll(z, 1, 0)

    def local_disorder(self):
        if self.p.twist == 0.0:
            z = np.exp(1j * self.th)
            nb = (np.roll(z, 1, 0) + np.roll(z, -1, 0) + np.roll(z, 1, 1) + np.roll(z, -1, 1)) / 4.0
        else:
            a, b, c, d = self._aligned_phasors()
            nb = (a + b + c + d) / 4.0
        return 1.0 - np.abs(nb)

    def kill_multiplier(self, D):
        p = self.p
        return 1.0 + p.beta * D ** p.hill_h / (p.hill_K ** p.hill_h + D ** p.hill_h)

    def site_update(self):
        p = self.p
        D = self.local_disorder()
        mu_arr = self.kill_multiplier(D)
        lam = p.lam_eff * mu_arr
        # what the kill step actually sees: damaged fraction and mean kill multiplier of the intact sites
        self.last_q = 1.0 - float(self.intact.mean())
        self.last_mu = float(mu_arr[self.intact].mean()) if self.intact.any() else float("nan")
        if self.intact.any():           # the multiplier a feedback of reference strength 20 would apply to this configuration
            x = D[self.intact]
            self.last_mu_ref = float(np.mean(1.0 + 20.0 * x ** p.hill_h / (p.hill_K ** p.hill_h + x ** p.hill_h)))
        else:
            self.last_mu_ref = float("nan")
        kill = self.intact & (self.rng.random(D.shape) < lam)
        self.intact &= ~kill
        self.g[kill] = p.g0
        rep = (~self.intact) & (self.rng.random(D.shape) < p.R)
        self.intact |= rep
        self.g[rep] = 1.0
        if p.repair == "reset0":
            self.th[rep] = 0.0
        elif p.repair == "neighbor" and rep.any():
            a, b, c, d = self._aligned_phasors()
            ga, gb = np.roll(self.g, -1, 1), np.roll(self.g, 1, 1)
            gc, gd = np.roll(self.g, -1, 0), np.roll(self.g, 1, 0)
            zn = ga * a + gb * b + gc * c + gd * d
            ok = rep & (np.abs(zn) > 1e-12)
            self.th[ok] = np.angle(zn[ok])
        elif p.repair not in ("keep", "reset0", "neighbor"):
            raise ValueError(p.repair)
        return int(kill.sum()), int(rep.sum())

    def step(self):
        for _ in range(self.p.relax):
            self.metropolis()
        return self.site_update()

    # ------------------------------------------------------------ observables
    def order(self):
        """U = |<exp(i theta)>| over all sites."""
        return float(abs(np.mean(np.exp(1j * self.th))))

    def order_intact(self):
        m = self.intact
        return float(abs(np.mean(np.exp(1j * self.th[m])))) if m.any() else 0.0

    def damaged_fraction(self):
        return float(1.0 - self.intact.mean())

    def disorder_stats(self, ref_beta=None):
        """Mean local disorder and mean kill multiplier over intact sites; with ref_beta, also the mean multiplier that
        a feedback of strength ref_beta would apply to this configuration (used to test the closure at beta = 0)."""
        D = self.local_disorder()
        m = self.intact
        if not m.any():
            return (float("nan"),) * (2 if ref_beta is None else 3)
        out = (float(D[m].mean()), float(self.kill_multiplier(D[m]).mean()))
        if ref_beta is not None:
            p = self.p
            x = D[m]
            out += (float(np.mean(1.0 + ref_beta * x ** p.hill_h / (p.hill_K ** p.hill_h + x ** p.hill_h))),)
        return out

    def bond_terms(self):
        """Per-site means of J cos(d - phi_x) and J sin(d - phi_x) for the x and y bonds (d = theta_j - theta_i)."""
        g, th, phi = self.g, self.th, self.p.twist
        out = []
        for ax, off in ((1, phi), (0, 0.0)):
            J = g * np.roll(g, -1, ax)
            d = np.roll(th, -1, ax) - th - off
            out += [float(np.mean(J * np.cos(d))), float(np.mean(J * np.sin(d)))]
        return tuple(out)          # (c_x, s_x, c_y, s_y)

    def winding_x(self):
        """Winding number of the phases along x, per row: sum of wrapped bond differences / 2 pi. Returns (mean over rows,
        standard deviation over rows); integer and equal in all rows unless vortices are present."""
        d = np.roll(self.th, -1, 1) - self.th
        d = (d + np.pi) % (2 * np.pi) - np.pi
        w = d.sum(axis=1) / (2 * np.pi)
        return float(w.mean()), float(w.std())

    def correlation(self, rmax=None):
        """C(r) = <cos(theta_i - theta_{i+r})> over all sites, and over pairs of intact sites, averaged over x and y."""
        L = self.p.L
        rmax = rmax or L // 2
        th = self.th
        w = self.intact.astype(float)
        call, cint = np.zeros(rmax), np.zeros(rmax)
        for r in range(1, rmax + 1):
            a = c = 0.0
            wsum = 0.0
            for ax in (0, 1):
                cs = np.cos(th - np.roll(th, -r, ax))
                ww = w * np.roll(w, -r, ax)
                a += cs.mean()
                c += (cs * ww).sum()
                wsum += ww.sum()
            call[r - 1] = a / 2.0
            cint[r - 1] = c / wsum if wsum > 0 else np.nan
        return call, cint


# ---------------------------------------------------------------- quenched dilution (field-free finite-size analysis)
class Quenched:
    """Field-free XY model on a quenched site-diluted lattice (removed sites carry g = 0)."""

    def __init__(self, L, T, f, seed, clusters=0, rng_extra=None):
        self.L, self.T = L, T
        self.rng = np.random.default_rng(seed)
        self.parity = np.add.outer(np.arange(L), np.arange(L)) % 2
        self.g = np.ones((L, L))
        if clusters and clusters > 0:
            self.g = self._clustered(L, f, clusters)
        else:
            self.g[self.rng.random((L, L)) < f] = 0.0
        self.th = self.rng.normal(0.0, 0.1, (L, L))

    def _clustered(self, L, f, radius):
        """Damage in compact discs of the given radius (in lattice units), placed at random centres until the
        removed fraction reaches f (the last disc is truncated so that the fraction is exact)."""
        target = int(round(f * L * L))
        g = np.ones((L, L))
        ii, jj = np.meshgrid(np.arange(-radius, radius + 1), np.arange(-radius, radius + 1), indexing="ij")
        disc = [(a, b) for a, b in zip(ii.ravel(), jj.ravel()) if a * a + b * b <= radius * radius + 1e-9]
        removed = 0
        while removed < target:
            ci, cj = self.rng.integers(0, L, 2)
            cells = [((ci + a) % L, (cj + b) % L) for a, b in disc]
            self.rng.shuffle(cells)
            for (x, y) in cells:
                if g[x, y] == 1.0:
                    g[x, y] = 0.0
                    removed += 1
                    if removed >= target:
                        break
        return g

    def sweep(self):
        L, th, g = self.L, self.th, self.g
        for col in (0, 1):
            m = self.parity == col
            prop = th + self.rng.uniform(-np.pi, np.pi, (L, L))
            e0 = np.zeros((L, L))
            e1 = np.zeros((L, L))
            for sh, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
                J = g * np.roll(g, sh, ax)
                tj = np.roll(th, sh, ax)
                e0 -= J * np.cos(th - tj)
                e1 -= J * np.cos(prop - tj)
            acc = (self.rng.random((L, L)) < np.exp(-(e1 - e0) / self.T)) & m
            th = np.where(acc, prop, th)
        self.th = th

    def magnetization(self):
        n = max(self.g.sum(), 1.0)
        return float(abs(np.sum(self.g * np.exp(1j * self.th))) / n)

    def helicity(self):
        """Helicity modulus, averaged over the two lattice directions."""
        L, th, g, T = self.L, self.th, self.g, self.T
        ups = []
        for ax in (0, 1):
            J = g * np.roll(g, -1, ax)
            d = th - np.roll(th, -1, ax)
            ups.append(np.sum(J * np.cos(d)) / (L * L) - np.sum(J * np.sin(d)) ** 2 / (L * L * T))
        return 0.5 * (ups[0] + ups[1])
