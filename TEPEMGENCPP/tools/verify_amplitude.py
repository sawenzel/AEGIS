#!/usr/bin/env python3
"""First-principles check of diffcross.f / DiffCross.h (O2-6340, checkbox 1).

Computes the 5-fold differential cross section for e+e- pair production in
the external-field (two-photon) approach *directly from the Feynman rules*:
explicit Dirac spinors, the two lowest-order diagrams, numerical 2-D
integration over the transverse photon momentum. No part of the code's
analytic reduction (the N1..N18 terms, the Iz/Id/Iv integrals) is reused.

The claim being verified (see ../doc/07-derivation.md):

    NT = Integral d^2k  Sum_spins | ubar(p-) Gamma(k) v(p+) |^2 / (q1^2 q2^2)^2

with  Gamma(k) = w1s (p-s - q1s + m) w2s / [(p- - q1)^2 - m^2]
              + w2s (q1s - p+s + m) w1s / [(q1 - p+)^2 - m^2]

(`Xs` = slash), w1/w2 = (1,0,0,+-beta) the nuclei's coordinate velocities,
and q1 the photon from nucleus 1, whose longitudinal part is fixed by
q1.w1 = 0 and (P - q1).w2 = 0 to

    q1_l = qb (gamma*beta, gamma),   qb = (P.w2) / (2 gamma beta).

The code's dsigma (kbarn/MeV^4/(Zalpha)^4) should then equal

    dsigma = NT * (4/beta^2) * (2pi)^-4 * (1/4) * (1.9733^2/10)

which is exactly the normalisation chain at the end of SUBROUTINE Diffcross.
The reference value comes from ./driver (the validated C++ port in quad).

Usage:  python3 verify_amplitude.py [--points N]
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Dirac algebra (Dirac representation, metric +---)
# ---------------------------------------------------------------------------
I2 = np.eye(2)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
G0 = np.block([[I2, 0 * I2], [0 * I2, -I2]]).astype(complex)
GAM = [G0] + [np.block([[0 * I2, s], [-s, 0 * I2]]).astype(complex)
              for s in (SX, SY, SZ)]


def slash(a):
    """a_mu gamma^mu for a four-vector a = (a0, ax, ay, az); a may carry
    leading broadcast axes (grid batching)."""
    a = np.asarray(a, dtype=complex)
    return (a[..., 0, None, None] * GAM[0]
            - a[..., 1, None, None] * GAM[1]
            - a[..., 2, None, None] * GAM[2]
            - a[..., 3, None, None] * GAM[3])


def mdot(a, b):
    """Minkowski product, +---."""
    return (a[..., 0] * b[..., 0] - a[..., 1] * b[..., 1]
            - a[..., 2] * b[..., 2] - a[..., 3] * b[..., 3])


def spinor_u(p, m):
    """u(p, s) for s = 0, 1; columns normalised to ubar u = 2m."""
    E = p[0]
    sp = p[1] * SX + p[2] * SY + p[3] * SZ
    out = []
    for s in range(2):
        chi = np.zeros(2, dtype=complex)
        chi[s] = 1.0
        top = np.sqrt(E + m) * chi
        bot = sp @ chi / np.sqrt(E + m)
        out.append(np.concatenate([top, bot]))
    return out


def spinor_v(p, m):
    """v(p, s) for s = 0, 1; vbar v = -2m, sum_s v vbar = pslash - m."""
    E = p[0]
    sp = p[1] * SX + p[2] * SY + p[3] * SZ
    out = []
    for s in range(2):
        eta = np.zeros(2, dtype=complex)
        eta[s] = 1.0
        top = sp @ eta / np.sqrt(E + m)
        bot = np.sqrt(E + m) * eta
        out.append(np.concatenate([top, bot]))
    return out


# ---------------------------------------------------------------------------
# The integrand
# ---------------------------------------------------------------------------
def make_integrand(energy, mass, ppt, yp, pmt, ym, dphi):
    """Returns f(kx, ky) = Sum_spins |T|^2 (vectorised over a grid) plus the
    kinematic record used by the derivation cross-check."""
    m = mass
    gamma = energy / (2.0 * 0.938)
    beta = np.sqrt((1 - 1 / gamma) * (1 + 1 / gamma))

    mtp = np.hypot(ppt, m)          # transverse masses
    mtm = np.hypot(pmt, m)
    pp = np.array([mtp * np.cosh(yp), ppt, 0.0, mtp * np.sinh(yp)])
    pm = np.array([mtm * np.cosh(ym), pmt * np.cos(dphi),
                   pmt * np.sin(dphi), mtm * np.sinh(ym)])
    P = pp + pm
    w1 = np.array([1.0, 0.0, 0.0, beta])
    w2 = np.array([1.0, 0.0, 0.0, -beta])

    # Longitudinal photon momentum fixed by the two delta functions:
    # q.w1 = 0 and q.w2 = P.w2  =>  q_l = qb*(gamma*beta, 0, 0, gamma).
    qb = mdot(P, w2) / (2 * gamma * beta)
    qb2 = mdot(P, w1) / (2 * gamma * beta)          # nucleus-2 analogue

    # Stable longitudinal parts of the four denominators, in the same
    # hyperbolic closed forms diffcross.f uses. Evaluating q2^2 = (P-q)^2 or
    # (p- - q)^2 - m^2 componentwise instead loses ~5 digits at forward
    # rapidity (P^0 and q^0 are both ~m_t cosh y_max while the difference is
    # ~qb2*gamma), which is visible at the 1e-6 level in this check -- the
    # very instability the Fortran's split into longitudinal scalar products
    # avoids by construction.
    wy = np.arccosh(gamma)
    wl = 1.0 / gamma
    pmlxql = qb * mtm * np.sinh(wy - ym)
    pplxql = qb * mtp * np.sinh(wy - yp)
    md = m * m + qb * qb - mtm * mtm + 2 * pmlxql   # -> (p- - q)^2 - m^2
    mx = m * m + qb * qb - mtp * mtp + 2 * pplxql   # -> (q - p+)^2 - m^2
    shift_d = -np.array([pm[1], pm[2]])             # the code's kd
    shift_x = -np.array([pp[1], pp[2]])             # the code's kx
    shift_1 = -np.array([P[1], P[2]])               # the code's k1

    ubar_s = [u.conj() @ G0 for u in spinor_u(pm, m)]
    v_s = spinor_v(pp, m)
    w1s, w2s = slash(w1), slash(w2)

    def f(kx, ky):
        kx = np.asarray(kx, dtype=float)
        ky = np.asarray(ky, dtype=float)
        q = np.empty(kx.shape + (4,))
        q[..., 0] = qb * gamma * beta
        q[..., 1] = kx
        q[..., 2] = ky
        q[..., 3] = qb * gamma
        ksq = kx * kx + ky * ky
        q1sq = -(ksq + qb * qb)                     # q1^2, exactly
        q2sq = -((kx + shift_1[0]) ** 2 + (ky + shift_1[1]) ** 2 + qb2 * qb2)
        Dd = -((kx + shift_d[0]) ** 2 + (ky + shift_d[1]) ** 2 + md)
        Dx = -((kx + shift_x[0]) ** 2 + (ky + shift_x[1]) ** 2 + mx)

        num_d = slash(pm - q) + m * np.eye(4)
        num_x = slash(q - pp) + m * np.eye(4)
        Gam = (w1s @ num_d @ w2s) / Dd[..., None, None] \
            + (w2s @ num_x @ w1s) / Dx[..., None, None]

        tot = np.zeros(kx.shape)
        for ub in ubar_s:
            for v in v_s:
                amp = np.einsum('a,...ab,b->...', ub, Gam, v)
                tot += np.abs(amp) ** 2
        return tot / (q1sq * q2sq) ** 2

    # Spike geometry for the integrator: photon-1 pole at k = 0 with width
    # qb = sqrt(m0); photon-2 pole at k = P_perp with width sqrt(m1), and
    # m1 = (P.w1 / (2 gamma beta))^2 -- nucleus 2's own qb (see derivation.md).
    qb2 = mdot(P, w1) / (2 * gamma * beta)
    centers = (np.array([0.0, 0.0]), np.array([P[1], P[2]]))
    widths = (abs(qb), abs(qb2))

    return f, {'gamma': gamma, 'beta': beta, 'qb': qb,
               'centers': centers, 'widths': widths}


def integrate(f, centers, widths, kmax=1e4, nrad=500, nang=384):
    """Integral of f over R^2.

    The integrand has one narrow peak per squared photon propagator: at k = 0
    (width sqrt(m0) = qb) and at k = P_perp (width sqrt(m1), photon 2's own
    qb). A single origin-centred polar grid cannot resolve the off-centre
    spike, so the plane is split into the two Voronoi half-planes of the
    spike centres and each half is integrated in log-radial polar coordinates
    about its own spike.
    """
    x, wgt = np.polynomial.legendre.leggauss(nrad)
    total = 0.0
    for a, b, w in ((centers[0], centers[1], widths[0]),
                    (centers[1], centers[0], widths[1])):
        d = b - a
        dist = np.hypot(*d)
        ang0 = np.arctan2(d[1], d[0])
        rmin = max(w * 1e-6, 1e-12)
        phi = ang0 + np.linspace(0.0, 2 * np.pi, nang, endpoint=False)
        dphi_w = 2 * np.pi / nang
        # Half-plane {u . dhat < dist/2}: radial limit dist/(2 cos(phi-ang0))
        # when that cosine is positive, else the far cutoff.
        cosr = np.cos(phi - ang0)
        L = np.where(cosr > 1e-12, dist / (2 * np.maximum(cosr, 1e-12)), kmax)
        L = np.minimum(L, kmax)
        # Log-radial Gauss-Legendre per angle: r = exp(t), d^2k = r^2 dt dphi.
        t0 = np.log(rmin)
        t1 = np.log(L)                                   # shape (nang,)
        t = 0.5 * (t1[None, :] - t0) * x[:, None] + 0.5 * (t1[None, :] + t0)
        wt = wgt[:, None] * 0.5 * (t1[None, :] - t0)
        R = np.exp(t)
        KX = a[0] + R * np.cos(phi)[None, :]
        KY = a[1] + R * np.sin(phi)[None, :]
        vals = f(KX, KY)
        total += np.sum(wt * vals * R * R) * dphi_w
    return total


def reference_dsigma(energy, mass, ppt, yp, pmt, ym, dphi):
    out = subprocess.run(
        [str(HERE / 'derivation_driver'), 'point', str(energy), str(mass), str(ppt),
         str(yp), str(pmt), str(ym), str(dphi)],
        capture_output=True, text=True, check=True)
    return float(out.stdout.split()[0])


def check_point(energy, mass, ppt, yp, pmt, ym, dphi, verbose=True):
    f, kin = make_integrand(energy, mass, ppt, yp, pmt, ym, dphi)
    nt = integrate(f, kin['centers'], kin['widths'])
    nt2 = integrate(f, kin['centers'], kin['widths'],
                    nrad=700, nang=576)             # convergence control
    beta = kin['beta']
    norm = (4 / beta**2) / (2 * np.pi)**4 / 4 * (1.9733**2 / 10)
    ds_fp = nt2 * norm
    ds_ref = reference_dsigma(energy, mass, ppt, yp, pmt, ym, dphi)
    quad_err = abs(nt2 - nt) / abs(nt2)
    rel = ds_fp / ds_ref - 1
    if verbose:
        print(f"  point (ppt={ppt}, yp={yp}, pmt={pmt}, ym={ym}, dphi={dphi})")
        print(f"    first-principles dsigma = {ds_fp:.12e}")
        print(f"    diffcross (quad) dsigma = {ds_ref:.12e}")
        print(f"    ratio-1 = {rel:.2e}   (quadrature self-error {quad_err:.1e})")
    return rel, quad_err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--energy', type=float, default=5360.0)
    ap.add_argument('--mass', type=float, default=0.511)
    args = ap.parse_args()

    # Benign points (away from the dphi->pi / ym->yp ill-conditioned locus,
    # where the *reference* is trustworthy at quad and double quadrature of
    # the first-principles integrand is comfortable), spanning the production
    # phase space: pt 1..1000 MeV, |y| < 7.
    points = [
        (5.0, 1.0, 3.0, -0.5, 2.0),
        (1.0, 0.0, 1.0, 1.5, 1.0),
        (20.0, 3.0, 7.0, -2.0, 2.5),
        (2.0, -4.0, 0.7, -1.0, 0.7),
        (100.0, 0.5, 40.0, 2.0, 1.2),
        (1.5, 5.5, 1.0, 4.0, 2.9),
    ]
    print(f"energy={args.energy} GeV, mass={args.mass} MeV")
    worst = 0.0
    for ppt, yp, pmt, ym, dphi in points:
        rel, qerr = check_point(args.energy, args.mass, ppt, yp, pmt, ym, dphi)
        worst = max(worst, abs(rel))
    print(f"\nworst |ratio-1| over {len(points)} points: {worst:.2e}")
    ok = worst < 1e-6
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
