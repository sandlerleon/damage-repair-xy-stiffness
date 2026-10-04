# -*- coding: utf-8 -*-
"""Mean-field reduction of the damage-repair kinetics, derived from the implemented transition probabilities.

Variables (all per site):
  q   damaged fraction (a site is "damaged" when g = g0, "intact" when g = 1); p = 1 - q is the intact fraction
  D_i local phase disorder of site i, D_i = 1 - |(1/4) sum_{j in nn(i)} exp(i theta_j)|
  mu(D) = 1 + beta D^h / (K^h + D^h): the kill multiplier; an intact site is damaged with probability lam0 mu(D_i)
  R   repair probability of a damaged site per step;  rho = R / lam0

One step of the code kills first and then repairs, so a site damaged in a step can be repaired in the same step:

    q' = (1 - R) [ q + (1 - q) lam0 mubar ],        mubar = < mu(D_i) >  over intact sites.

Fixed points: q = (1 - R) [q + (1 - q) lam0 mubar]  <=>  rho = (1 - q) mubar / (q + (1 - q) lam0 mubar).
For lam0, R << 1 this is the rate equation, in the time tau = lam0 t,

    dq/dtau = mubar (1 - q) - rho q,        steady state  rho = mubar (1 - q) / q.

q is the DAMAGED fraction. For beta = 0 (mubar = 1) the steady state is q* = 1/(1 + rho) (discrete time: the same with
lam0 corrections) and the intact fraction is rho/(1 + rho); the order U is not a damaged fraction, and a plot of
1/(1 + rho) against U is not a comparison of like with like.

The only approximation is the closure mubar = mu(q): the local disorder of an intact site is replaced by the damaged
fraction (the annealed mean-field identification D = q, motivated by the slaving of the phase order to the intact
fraction, U ~ 1 - q). It is tested against the lattice in analyze.py: the closure can be replaced by the measured
mubar(q) of the stationary lattice, which gives the corresponding reduced steady-state relation.
"""
import numpy as np


def mu(x, beta, h=6.0, K=0.5):
    x = np.asarray(x, float)
    return 1.0 + beta * x ** h / (K ** h + x ** h)


def rho_ode(q, beta, h=6.0, K=0.5):
    """Steady-state relation of the rate equation with the closure mubar = mu(q)."""
    return mu(q, beta, h, K) * (1.0 - q) / q


def rho_map(q, beta, lam0=0.02, h=6.0, K=0.5, mubar=None):
    """Steady-state relation of the discrete-time map (kill then repair, same step)."""
    m = mu(q, beta, h, K) if mubar is None else mubar
    return (1.0 - q) * m / (q + (1.0 - q) * lam0 * m)


def window(beta, h=6.0, K=0.5, lam0=None, n=40000):
    """Bistable window (rho_minus, rho_plus) of the closure model, or None if the steady state is unique.

    The relation rho(q) is S-shaped when it has an interior local minimum followed by a local maximum in q."""
    q = np.linspace(1e-5, 1 - 1e-5, n)
    r = rho_ode(q, beta, h, K) if lam0 is None else rho_map(q, beta, lam0, h, K)
    d = np.sign(np.diff(r))
    idx = np.where(d[:-1] != d[1:])[0]
    mins = [(q[i + 1], r[i + 1]) for i in idx if d[i] < 0 < d[i + 1]]
    maxs = [(q[i + 1], r[i + 1]) for i in idx if d[i] > 0 > d[i + 1]]
    if mins and maxs:
        return float(mins[0][1]), float(maxs[0][1]), float(mins[0][0]), float(maxs[0][0])
    return None


def beta_c(h=6.0, K=0.5, lo=0.0, hi=400.0):
    if window(hi, h, K) is None:
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if window(mid, h, K) is not None:
            hi = mid
        else:
            lo = mid
    return hi


def is_monotonic_decreasing(q, rho, tol=0.0):
    """Is the measured steady-state relation rho(q) strictly decreasing (a single steady state for each rho)?"""
    o = np.argsort(q)
    r = np.asarray(rho)[o]
    return bool(np.all(np.diff(r) < tol)), float(np.max(np.diff(r)))
