#!/usr/bin/env python3
"""Identify and verify the ten Iz/Id/Iv closed forms of diffcross.f.

Each function in diffcross.f (and the C++ port) is a closed-form expression
for a 2-D transverse-momentum integral. This script establishes *which*
integral each one is, by numerically integrating candidate integrands

    1 / [ (k^2+u)^a ((k+x)^2+v)^b ((k+y)^2+w)^c ((k+z)^2+w2)^d ]

over R^2 at generic O(1) arguments and comparing against the closed form
evaluated by ./derivation_driver (quad precision). Agreement at the quadrature tolerance
for exactly one exponent tuple identifies the definition; that tuple is what
../doc/07-derivation.md then states.

Usage: python3 verify_integrals.py
"""
import itertools
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# Generic, well-separated O(1) arguments: nothing symmetric, nothing small.
X = (0.7, 0.3)
Y = (-0.4, 0.9)
Z = (0.5, -0.6)
U, V, W, W2 = 1.3, 0.9, 1.7, 2.1


def driver(cmd, *args):
    out = subprocess.run([str(HERE / 'derivation_driver'), cmd] + [str(a) for a in args],
                         capture_output=True, text=True, check=True)
    return float(out.stdout.split()[0])


def integrate(f, nrad=600, nang=512):
    """Integral over R^2, polar coordinates, r = t/(1-t) map to [0, inf)."""
    t, wgt = np.polynomial.legendre.leggauss(nrad)
    t = 0.5 * (t + 1)
    wgt = 0.5 * wgt
    r = t / (1 - t)
    jac = 1 / (1 - t) ** 2
    phi = np.linspace(0, 2 * np.pi, nang, endpoint=False)
    R, PHI = np.meshgrid(r, phi, indexing='ij')
    vals = f(R * np.cos(PHI), R * np.sin(PHI))
    return np.einsum('i,ij->', wgt * r * jac, vals) * (2 * np.pi / nang)


def denom(kx, ky, shift, mass, power):
    return ((kx + shift[0]) ** 2 + (ky + shift[1]) ** 2 + mass) ** power


def candidate(shifts_masses, powers):
    def f(kx, ky):
        out = 1.0
        for (s, m), p in zip(shifts_masses, powers):
            out = out / denom(kx, ky, s, m, p)
        return out
    return f


def identify(name, ref, sm, max_power=3):
    """Try all exponent tuples up to max_power, report matches."""
    n = len(sm)
    hits = []
    for powers in itertools.product(range(1, max_power + 1), repeat=n):
        if sum(powers) > n + 3:      # keep the search cheap and plausible
            continue
        val = integrate(candidate(sm, powers))
        rel = abs(val / ref - 1)
        if rel < 1e-6:
            hits.append((powers, rel))
    for powers, rel in hits:
        print(f"  {name} = Int d^2k prod 1/(...)^{list(powers)}   "
              f"(rel {rel:.1e})")
    if len(hits) != 1:
        print(f"  {name}: {len(hits)} matches -- NOT uniquely identified")
    return len(hits) == 1


def main():
    ok = True
    # x shifts the *second* denominator: Iz0(x,u,v) pairs (k^2+u), ((k+x)^2+v).
    sm2 = [((0.0, 0.0), U), (X, V)]
    sm3 = [((0.0, 0.0), U), (X, V), (Y, W)]
    sm4 = [((0.0, 0.0), U), (X, V), (Y, W), (Z, W2)]

    print("two denominators (Iz):")
    ok &= identify('Iz0', driver('iz0', *X, U, V), sm2)
    ok &= identify('Iz1', driver('iz1', *X, U, V), sm2)
    ok &= identify('Iz2', driver('iz2', *X, U, V), sm2)

    print("three denominators (Id):")
    args3 = (*X, *Y, U, V, W)
    ok &= identify('Id0', driver('id0', *args3), sm3)
    ok &= identify('Id1', driver('id1', *args3), sm3)
    ok &= identify('Id2', driver('id2', *args3), sm3)
    ok &= identify('Id3', driver('id3', *args3), sm3)

    print("four denominators (Iv):")
    args4 = (*X, *Y, *Z, U, V, W, W2)
    ok &= identify('Iv0', driver('iv0', *args4), sm4)
    ok &= identify('Iv1', driver('iv1', *args4), sm4)
    ok &= identify('Iv2', driver('iv2', *args4), sm4)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
