#!/usr/bin/env python3
"""Exact term-level verification of the 18 N-terms of diffcross.f (O2-6340).

verify_amplitude.py established numerically (1e-12) that the code's NT equals
the spin-summed squared external-field amplitude. This tool upgrades that to
an EXACT statement, in rational arithmetic throughout:

  1. Kinematics is chosen exactly rational (pythagorean transverse masses,
     rational exponentials for every rapidity, a pythagorean azimuth), so all
     longitudinal invariants, shift vectors and propagator masses are exact.
  2. The spin sum is an exact Dirac trace (explicit 4x4 gamma matrices), a
     polynomial in the transverse loop momentum k with rational coefficients.
  3. A reduction engine rewrites any convergent  Int d^2k N(k)/prod D_i^{n_i}
     (D_i = (k+s_i)^2 + m_i) into terminal scalar integrals, using only
     pointwise partial fractions (k is affine in denominator differences),
     reflection parity, angular averaging, and box splitting (1 in span{D_i}
     for four denominators in 2-D). All moves are exact.
  4. The trace integrand and the code's own sum -- coefficients parsed from
     the Fortran source, each integral call mapped onto the physical
     denominators -- are reduced to the same terminal basis. The two
     decompositions differ (terminal triangles/bubbles are not independent
     in 2-D), so the difference is shown to be an EXACT rational combination
     of null relations: IBP identities Int div(v/prod D^n) d^2k = 0 for
     polynomial vector fields v, plus multiply-divide relations obtained by
     reducing D_d * 1/(D_d prod D^n) two ways (these carry the 2-D Gram
     redundancy -- the same one whose degeneration at dphi -> pi causes the
     cancellation catastrophe, doc/07-derivation.md section 8).
     Membership is decided by Gaussian elimination over two independent
     61/62-bit prime fields and then CONFIRMED exactly over Q on the
     selected subsystem.
  5. The whole check runs at two independent kinematic points
     (Schwartz-Zippel), so agreement is not a coincidence of one point.

A successful run therefore proves: at exact rational kinematics, the sum
N1+...+N18 of diffcross.f is identically the integral of the spin-summed
squared amplitude of doc/07-derivation.md section 3 -- not merely equal to
quadrature accuracy.

Needs sympy and numpy. Runtime ~20-60 min. Usage: python3 verify_terms.py
"""

# ======================= reduction engine =======================

import sympy as sp
from sympy import Rational as R

kx, ky = sp.symbols('kx ky')


def dkey(dens):
    return ('int', frozenset(dens.items()))


class Reducer:
    def __init__(self):
        self.out = {}

    def add(self, key, coef):
        if coef == 0:
            return
        self.out[key] = self.out.get(key, R(0)) + coef
        if self.out[key] == 0:
            del self.out[key]

    # ------------------------------------------------------------------
    def reduce(self, N, dens):
        """N: sympy expr in kx,ky (rational coefficients);
        dens: dict {(sx,sy,m): power>0}."""
        N = sp.expand(N)
        if N == 0:
            return
        dens = {d: p for d, p in dens.items() if p > 0}
        nd = len(dens)
        if nd == 0:
            self.add(('div', 'poly'), R(1))  # unnormalised marker
            return
        poly = sp.Poly(N, kx, ky)
        if poly.total_degree() == 0:
            self._constant(N, dens)
            return
        shifts = list(dens)
        if nd >= 3:
            self._eliminate_full(N, dens)
        elif nd == 2:
            self._eliminate_two(N, dens)
        else:
            self._eliminate_one(N, dens)

    # ------------------------------------------------------------------
    def _Dexpr(self, d):
        sx, sy, m = d
        return (kx + sx) ** 2 + (ky + sy) ** 2 + m

    def _eliminate_full(self, N, dens):
        """>=3 denominators. Need two linearly independent shift
        differences; if all shifts are collinear, fall back to the
        two-shift treatment on the collinear axis (not needed at generic
        kinematics; assert instead)."""
        ds = list(dens)
        pivot = sorted(ds)[0]
        others = [d for d in ds if d != pivot]
        # find two others with independent differences
        pair = None
        for i in range(len(others)):
            for j in range(i + 1, len(others)):
                t1 = (others[i][0] - pivot[0], others[i][1] - pivot[1])
                t2 = (others[j][0] - pivot[0], others[j][1] - pivot[1])
                if t1[0] * t2[1] - t1[1] * t2[0] != 0:
                    pair = (others[i], others[j])
                    break
            if pair:
                break
        assert pair, "collinear >=3-shift configuration not supported"
        d1, d2 = pair
        # D_a - D_p = 2 k.(s_a-s_p) + (s_a^2+m_a - s_p^2 - m_p)
        A, B_ = sp.symbols('DA DB')  # stand-ins for D_{d1}, D_{d2}
        P_ = sp.Symbol('DP')         # D_pivot
        sxp, syp, mp_ = pivot

        def cdiff(d):
            sx, sy, m = d
            return (sx * sx + sy * sy + m) - (sxp * sxp + syp * syp + mp_)

        t1 = (d1[0] - sxp, d1[1] - syp)
        t2 = (d2[0] - sxp, d2[1] - syp)
        det = 2 * (t1[0] * t2[1] - t1[1] * t2[0])
        # solve 2 k.t1 = (DA - DP) - c1 ; 2 k.t2 = (DB - DP) - c2
        L1 = (A - P_) - cdiff(d1)
        L2 = (B_ - P_) - cdiff(d2)
        kx_sub = (t2[1] * L1 - t1[1] * L2) / det
        ky_sub = (t1[0] * L2 - t2[0] * L1) / det
        Nsub = sp.expand(N.subs({kx: kx_sub, ky: ky_sub}, simultaneous=True))
        # now a polynomial in DP, DA, DB
        p = sp.Poly(Nsub, P_, A, B_)
        for monom, coef in zip(p.monoms(), p.coeffs()):
            ep, ea, eb = monom
            newd = dict(dens)
            num = sp.Integer(1)
            for d, e in ((pivot, ep), (d1, ea), (d2, eb)):
                cancel = min(e, newd[d])
                newd[d] -= cancel
                if newd[d] == 0:
                    del newd[d]
                for _ in range(e - cancel):
                    num *= self._Dexpr(d)
            self.reduce(coef * num, newd)

    def _eliminate_two(self, N, dens):
        (da, db) = sorted(dens)
        sxa, sya, ma = da
        sxb, syb, mb = db
        tx, ty = sxb - sxa, syb - sya
        t2 = tx * tx + ty * ty
        assert t2 != 0
        # centred coordinates u = k + s_a ; axis through centres -s_a, -s_b
        # is the line through u=0 with direction t. Reflection about it:
        # u -> M u with M = (1/t2) [[tx^2-ty^2, 2 tx ty],[2 tx ty, ty^2-tx^2]]
        ux = kx + sxa
        uy = ky + sya
        rux = ((tx * tx - ty * ty) * ux + 2 * tx * ty * uy) / t2
        ruy = (2 * tx * ty * ux + (ty * ty - tx * tx) * uy) / t2
        Nrefl = N.subs({kx: rux - sxa, ky: ruy - sya}, simultaneous=True)
        Nodd = sp.expand((N - Nrefl) / 2)      # integrates to zero
        Neven = sp.expand((N + Nrefl) / 2)
        if Neven == 0:
            return
        # Neven is polynomial in (L = u.t, S = u^2):
        # substitute via symbols and linear relations
        A, B_ = sp.symbols('DA DB')
        L = sp.Symbol('Lut')
        S = sp.Symbol('Suu')
        # express Neven in terms of L, S: use ux = (L tx + W ty)/t2-ish --
        # simpler: substitute ux, uy by generic and match. Robust approach:
        # write Neven as poly in (ux, uy) via k-sub, then rewrite monomials.
        Nu = sp.expand(Neven.subs({kx: sp.Symbol('ux') - sxa,
                                   ky: sp.Symbol('uy') - sya},
                                  simultaneous=True))
        uxs, uys = sp.Symbol('ux'), sp.Symbol('uy')
        # invariant coordinates: L = ux tx + uy ty ; W = -ux ty + uy tx
        # (W is odd under the reflection; Neven contains only even powers
        #  of W). W^2 = t2 * S - L^2.
        Ls, Ws = sp.symbols('Ls Ws')
        ux_sub = (Ls * tx - Ws * ty) / t2
        uy_sub = (Ls * ty + Ws * tx) / t2
        Nlw = sp.expand(Nu.subs({uxs: ux_sub, uys: uy_sub},
                                simultaneous=True))
        pw = sp.Poly(Nlw, Ws)
        acc = sp.Integer(0)
        for (ew,), coef in zip(pw.monoms(), pw.coeffs()):
            assert ew % 2 == 0, "odd W power survived even projection"
            acc += coef * (t2 * S - Ls ** 2) ** (ew // 2)
        # L = u.t = (D_b - D_a - t2 + m_a - m_b)/2 ; S = u^2 = D_a - m_a
        L_sub = (B_ - A - t2 + ma - mb) / 2
        S_sub = A - ma
        Nd = sp.expand(acc.subs({Ls: L_sub, S: S_sub}, simultaneous=True))
        p = sp.Poly(Nd, A, B_)
        for monom, coef in zip(p.monoms(), p.coeffs()):
            ea, eb = monom
            newd = dict(dens)
            num = sp.Integer(1)
            for d, e in ((da, ea), (db, eb)):
                cancel = min(e, newd[d])
                newd[d] -= cancel
                if newd[d] == 0:
                    del newd[d]
                for _ in range(e - cancel):
                    num *= self._Dexpr(d)
            self.reduce(coef * num, newd)

    def _eliminate_one(self, N, dens):
        (d,) = dens
        sx, sy, m = d
        p_ = dens[d]
        Nu = sp.expand(N.subs({kx: sp.Symbol('ux') - sx,
                               ky: sp.Symbol('uy') - sy}, simultaneous=True))
        uxs, uys = sp.Symbol('ux'), sp.Symbol('uy')
        poly = sp.Poly(Nu, uxs, uys)
        A = sp.Symbol('DA')
        acc = sp.Integer(0)
        from sympy import binomial, factorial2
        for (a, b), coef in zip(poly.monoms(), poly.coeffs()):
            if a % 2 or b % 2:
                continue                       # odd monomial: integrates to 0
            # angular average of ux^a uy^b over the circle:
            # <ux^a uy^b> = u^(a+b) (a-1)!!(b-1)!!/(a+b)!!
            avg = (factorial2(a - 1) * factorial2(b - 1)
                   / factorial2(a + b)) if (a or b) else 1
            acc += coef * avg * (A - m) ** ((a + b) // 2)
        pd = sp.Poly(acc, A)
        for (e,), coef in zip(pd.monoms(), pd.coeffs()):
            power = p_ - e
            if power <= 0:
                self.add(('div', 'power'), coef)
            elif power == 1:
                self.add(('div', f'logtad:{m}'), coef)
            else:
                # Int d^2k / (k^2+m)^power = pi / ((power-1) m^(power-1))
                self.add(('pi',), coef / ((power - 1) * m ** (power - 1)))

    # ------------------------------------------------------------------
    def _constant(self, c, dens):
        nd = len(dens)
        if nd == 4:
            # 1 = sum alpha_i D_i : basis {k^2, kx, ky, 1}
            ds = list(dens)
            M = sp.Matrix([[1, 2 * d[0], 2 * d[1],
                            d[0] ** 2 + d[1] ** 2 + d[2]] for d in ds]).T
            rhs = sp.Matrix([0, 0, 0, 1])
            alpha = M.solve(rhs)
            for d, a in zip(ds, alpha):
                if a == 0:
                    continue
                newd = dict(dens)
                newd[d] -= 1
                if newd[d] == 0:
                    del newd[d]
                self.reduce(c * a, newd)
            return
        if nd == 1:
            (d,) = dens
            m = d[2]
            p_ = dens[d]
            if p_ == 1:
                self.add(('div', f'logtad:{m}'), c)
            else:
                self.add(('pi',), c / ((p_ - 1) * m ** (p_ - 1)))
            return
        # bubble or triangle: terminal
        self.add(dkey(dens), c)


def reduce_integrand(N, dens):
    r = Reducer()
    r.reduce(N, dens)
    return r.out


KX, KY = kx, ky

# ================== null-relation generators ====================

import itertools
import sympy as sp
from sympy import Rational as R


VFIELDS = []
for deg0 in [ (sp.Integer(1), sp.Integer(0)), (sp.Integer(0), sp.Integer(1)),
              (KX, KY) ]:
    VFIELDS.append(deg0)
for mono in (KX, KY, KX*KX, KX*KY, KY*KY):
    VFIELDS.append((mono, sp.Integer(0)))
    VFIELDS.append((sp.Integer(0), mono))
VFIELDS.append((KX*KX*KX + KX*KY*KY, KX*KX*KY + KY*KY*KY))  # k^2 * k

def ibp_relations(seeds):
    rels = []
    for dens in seeds:
        sn = sum(dens.values())
        for (vx, vy) in VFIELDS:
            deg = max(sp.Poly(vx, KX, KY).total_degree() if vx else 0,
                      sp.Poly(vy, KX, KY).total_degree() if vy else 0)
            if 2*sn - deg < 3:      # surface term must vanish
                continue
            divv = sp.diff(vx, KX) + sp.diff(vy, KY)
            r = Reducer()
            if divv:
                r.reduce(sp.expand(divv), dict(dens))
            for d, n in dens.items():
                sx, sy, _ = d
                num = -n * 2 * (vx * (KX + sx) + vy * (KY + sy))
                raised = dict(dens)
                raised[d] += 1
                r.reduce(sp.expand(num), raised)
            if r.out:
                rels.append(r.out)
    return rels



def Dexpr(d):
    sx, sy, m = d
    return (KX + sx)**2 + (KY + sy)**2 + m

def md_relations(seeds, phys):
    rels = []
    for dens in seeds:
        base = Reducer(); base.reduce(sp.Integer(1), dict(dens))
        for d in phys:
            raised = dict(dens)
            raised[d] = raised.get(d, 0) + 1
            r = Reducer(); r.reduce(sp.expand(Dexpr(d)), raised)
            rel = dict(r.out)
            for k, c in base.out.items():
                rel[k] = rel.get(k, 0) - c
                if rel[k] == 0:
                    del rel[k]
            if rel:
                rels.append(rel)
    return rels


# ============ kinematics, trace, Fortran term parsing ===========

import re, sys
import sympy as sp
from sympy import Rational as R


from pathlib import Path
FORTRAN = str(Path(__file__).resolve().parents[2] / 'TEPEMGEN' / 'diffcross.f')

# ---------- exact rational kinematics ----------
def kinematics(m, pmvt, mtm, ppvt, mtp, chp, shp, chm, shm, gam, gb, cph, sph):
    beta = gb/gam; wl = 1/gam
    pm = sp.Matrix([mtm*chm, pmvt*cph, pmvt*sph, mtm*shm])
    pp = sp.Matrix([mtp*chp, ppvt, 0, mtp*shp])
    w1 = sp.Matrix([1,0,0,beta]); w2 = sp.Matrix([1,0,0,-beta])
    mdot = lambda a,b: a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3]
    S = {}
    S['m'] = m
    S['pmlxpml'] = mtm**2; S['pplxppl'] = mtp**2
    S['pmlxppl'] = mtm*mtp*(chp*chm - shp*shm)
    S['w2xpml'] = wl*mtm*(chm*gam + shm*gb)
    S['w2xppl'] = wl*mtp*(chp*gam + shp*gb)
    S['w1xpml'] = wl*mtm*(chm*gam - shm*gb)
    S['w1xppl'] = wl*mtp*(chp*gam - shp*gb)
    S['w1xw1'] = S['w2xw2'] = 1/gam**2
    S['w1xw2'] = 2 - 1/gam**2
    qb = gam*(S['w2xppl']+S['w2xpml'])/(2*gb*gam)
    S['qb'] = qb
    S['w1xql'] = R(0); S['w2xql'] = S['w2xpml']+S['w2xppl']
    S['qlxql'] = -qb**2
    S['pmlxql'] = qb*mtm*(gb*chm - gam*shm)
    S['pplxql'] = qb*mtp*(gb*chp - gam*shp)
    S['pmtxpmt'] = -pmvt**2
    S['pmtxppt'] = -pmvt*ppvt*cph
    S['pptxppt'] = -ppvt**2
    ppt = (ppvt, R(0)); pmt = (pmvt*cph, pmvt*sph)
    V = {}
    V['k1'] = (-ppt[0]-pmt[0], -pmt[1])
    V['kd'] = (-pmt[0], -pmt[1])
    V['kx'] = (-ppt[0], R(0))
    for n in ('k1','kd','kx'):
        V['m'+n] = (-V[n][0], -V[n][1])
    V['mk1d'] = (V['k1'][0]-V['kd'][0], V['k1'][1]-V['kd'][1])
    V['mk1x'] = (V['k1'][0]-V['kx'][0], V['k1'][1]-V['kx'][1])
    V['mkd1'] = (V['kd'][0]-V['k1'][0], V['kd'][1]-V['k1'][1])
    V['mkx1'] = (V['kx'][0]-V['k1'][0], V['kx'][1]-V['k1'][1])
    S['m0'] = -S['qlxql']
    S['m1'] = -S['qlxql'] - S['pplxppl'] - S['pmlxpml'] + 2*(S['pplxql']+S['pmlxql']-S['pmlxppl'])
    S['md'] = m**2 - (S['qlxql'] + S['pmlxpml'] - 2*S['pmlxql'])
    S['mx'] = m**2 - (S['qlxql'] + S['pplxppl'] - 2*S['pplxql'])
    S['r1'] = S['m1'] - S['m0'] + V['k1'][0]**2 + V['k1'][1]**2
    S['rd'] = S['md'] - S['m0'] + V['kd'][0]**2 + V['kd'][1]**2
    S['rx'] = S['mx'] - S['m0'] + V['kx'][0]**2 + V['kx'][1]**2
    # physical denominators (shift, mass) in the engine's convention
    D0 = (R(0), R(0), S['m0'])
    D1 = (V['k1'][0], V['k1'][1], S['m1'])
    Dd = (V['kd'][0], V['kd'][1], S['md'])
    Dx = (V['kx'][0], V['kx'][1], S['mx'])
    masses = [S['m0'], S['m1'], S['md'], S['mx']]
    assert len(set(masses)) == 4, "mass collision breaks call mapping"
    kin = dict(pm=pm, pp=pp, w1=w1, w2=w2, qb=qb, gam=gam, gb=gb, m=m,
               D0=D0, D1=D1, Dd=Dd, Dx=Dx, S=S, V=V)
    return kin

# ---------- exact trace ----------
def traces(kin):
    m = kin['m']
    I2 = sp.eye(2); Z2 = sp.zeros(2,2)
    sx = sp.Matrix([[0,1],[1,0]]); sy = sp.Matrix([[0,-sp.I],[sp.I,0]])
    sz = sp.Matrix([[1,0],[0,-1]])
    blk = lambda a,b,c,d: sp.Matrix(sp.BlockMatrix([[a,b],[c,d]]))
    G = [blk(I2,Z2,Z2,-I2), blk(Z2,sx,-sx,Z2), blk(Z2,sy,-sy,Z2), blk(Z2,sz,-sz,Z2)]
    slash = lambda a: a[0]*G[0]-a[1]*G[1]-a[2]*G[2]-a[3]*G[3]
    E4 = sp.eye(4)
    q = sp.Matrix([kin['qb']*kin['gb'], KX, KY, kin['qb']*kin['gam']])
    A = kin['pm'] - q; B = q - kin['pp']
    Gd = slash(kin['w1'])*(slash(A)+m*E4)*slash(kin['w2'])
    Gx = slash(kin['w2'])*(slash(B)+m*E4)*slash(kin['w1'])
    Gdb = slash(kin['w2'])*(slash(A)+m*E4)*slash(kin['w1'])
    Gxb = slash(kin['w1'])*(slash(B)+m*E4)*slash(kin['w2'])
    Pm = slash(kin['pm'])+m*E4; Pp = slash(kin['pp'])-m*E4
    Ndd = sp.expand(sp.trace(Pm*Gd*Pp*Gdb))
    Ndx = sp.expand(sp.trace(Pm*Gd*Pp*Gxb))
    Nxx = sp.expand(sp.trace(Pm*Gx*Pp*Gxb))
    return Ndd, Ndx, Nxx

# ---------- parse the 18 terms from the Fortran ----------
def parse_terms():
    src = open(FORTRAN).readlines()
    # join continuations within the Diffcross routine
    body = []
    for ln in src:
        if ln[:1].upper() == 'C':
            continue
        code = ln.rstrip('\n')
        if len(code) > 6 and code[5] not in ' 0' and body:
            body[-1] += code[6:]
        else:
            body.append(code[6:] if len(code) > 6 else code)
    text = [l.strip() for l in body]
    terms = []
    pat = re.compile(r'^N(\d+)\s*=\s*(Iz2|Id1|Id2|Id3|Iv1|Iv2)\s*\(([^)]*)\)\s*\*\s*\((.*)\)\s*$')
    for l in text:
        l = ' '.join(l.split())
        mm = pat.match(l)
        if mm:
            terms.append((int(mm.group(1)), mm.group(2),
                          [a.strip() for a in mm.group(3).split(',')],
                          mm.group(4)))
    assert len(terms) == 18, f"parsed {len(terms)} terms"
    return terms

POWERS = {'Iz2': (2,2), 'Id1': (2,1,1), 'Id2': (2,2,1), 'Id3': (2,2,2),
          'Iv1': (2,1,1,1), 'Iv2': (2,2,1,1)}

def call_descriptor(func, args, kin):
    """Map a diffcross call to {physical (s,m): power} in original k coords."""
    p = POWERS[func]
    nv = len(p) - 1                     # number of vector args
    vecs = [kin['V'][a] for a in args[:nv]]
    mass_names = args[nv:]
    masses = [kin['S'][a] for a in mass_names]
    # call-local denominators: ((0,0),masses[0],p[0]), (vecs[i],masses[i+1],p[i+1])
    local = [((R(0), R(0)), masses[0], p[0])] + \
            [(vecs[i], masses[i+1], p[i+1]) for i in range(nv)]
    phys = [kin['D0'], kin['D1'], kin['Dd'], kin['Dx']]
    # find s0 with s_phys = s_local + s0 for a mass-consistent assignment
    for cand in phys:
        for (sl, ml, _) in local:
            if ml == cand[2]:
                s0 = (cand[0]-sl[0], cand[1]-sl[1])
                desc = {}
                ok = True
                for (sl2, ml2, pw) in local:
                    target = (sl2[0]+s0[0], sl2[1]+s0[1])
                    match = [d for d in phys if d[2] == ml2 and
                             d[0] == target[0] and d[1] == target[1]]
                    if not match:
                        ok = False; break
                    desc[match[0]] = pw
                if ok and len(desc) == len(local):
                    return desc
    raise RuntimeError(f"no descriptor for {func}({args})")

def run(kinargs, label):
    kin = kinematics(*kinargs)
    Ndd, Ndx, Nxx = traces(kin)
    D0, D1, Dd, Dx = kin['D0'], kin['D1'], kin['Dd'], kin['Dx']
    r = Reducer()
    r.reduce(Ndd,   {Dd:2, D0:2, D1:2})
    r.reduce(2*Ndx, {Dd:1, Dx:1, D0:2, D1:2})
    r.reduce(Nxx,   {Dx:2, D0:2, D1:2})
    trace_state = r.out

    S = {k: v for k, v in kin['S'].items()}
    code = Reducer()
    for idx, func, args, coefstr in parse_terms():
        coef = sp.sympify(coefstr.lower().replace('d0',''), locals=S)
        coef = sp.nsimplify(coef, rational=True)
        assert coef.is_rational, f"N{idx} coefficient not rational: {coef}"
        desc = call_descriptor(func, args, kin)
        code.reduce(coef, desc)
    code_state = code.out

    keys = set(trace_state) | set(code_state)
    diffs = {k: trace_state.get(k, R(0)) - code_state.get(k, R(0)) for k in keys}
    diffs = {k: v for k, v in diffs.items() if v != 0}
    print(f"[{label}] trace-side terms: {len(trace_state)}, code-side terms: {len(code_state)}, "
          f"differing coefficients: {len(diffs)}")
    return trace_state, code_state, diffs, kin



# ==================== exact membership test =====================
def solve_membership(residual, rels):
    keys = sorted({k for st in rels for k in st} | set(residual), key=repr)
    ki = {k: i for i, k in enumerate(keys)}
    nk, nr = len(keys), len(rels)
    print(f"system: {nk} keys x {nr} relations", flush=True)
    cols = [[R(0)] * nk for _ in range(nr)]
    for j, st in enumerate(rels):
        for k, c in st.items():
            cols[j][ki[k]] = c
    bvec = [R(0)] * nk
    for k, c in residual.items():
        bvec[ki[k]] = c

    def modsolve(p):
        def toF(q):
            return q.p % p * pow(q.q % p, p - 2, p) % p
        A = [[toF(cols[j][i]) for j in range(nr)] + [toF(bvec[i])]
             for i in range(nk)]
        piv_cols, piv_rows = [], []
        row = 0
        for col in range(nr):
            sel = next((r_ for r_ in range(row, nk) if A[r_][col]), None)
            if sel is None:
                continue
            A[row], A[sel] = A[sel], A[row]
            inv = pow(A[row][col], p - 2, p)
            A[row] = [x * inv % p for x in A[row]]
            for r_ in range(nk):
                if r_ != row and A[r_][col]:
                    f = A[r_][col]
                    A[r_] = [(x - f * y) % p for x, y in zip(A[r_], A[row])]
            piv_cols.append(col)
            row += 1
            if row == nk:
                break
        # consistent iff no row with all-zero coefficients but nonzero rhs
        for r_ in range(row, nk):
            if A[r_][nr] % p:
                return None
        lam = {}
        for i, col in enumerate(piv_cols):
            lam[col] = A[i][nr]
        return lam

    for p in (2**61 - 1, (1 << 62) - 57):
        lam = modsolve(p)
        print(f"  membership mod p={p}: {'YES' if lam is not None else 'NO'}",
              flush=True)
        if lam is None:
            return False
    # exact confirmation over Q on the selected columns: membership holds
    # iff appending b to the selected columns does not raise the rank.
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    sel = sorted(lam)
    Msel = [[cols[j][i] for j in sel] for i in range(nk)]
    Maug = [row + [bvec[i]] for i, row in enumerate(Msel)]
    rk = DomainMatrix([[QQ(c.p, c.q) for c in row] for row in Msel],
                      (nk, len(sel)), QQ).rank()
    rk_aug = DomainMatrix([[QQ(c.p, c.q) for c in row] for row in Maug],
                          (nk, len(sel) + 1), QQ).rank()
    ok = rk == rk_aug
    print(f"  exact rank test over Q ({len(sel)} selected relations): "
          f"rank {rk} vs augmented {rk_aug} -> "
          f"{'verified' if ok else 'NOT verified'}", flush=True)
    return ok



# ============================ driver ============================
def verify_point(kinargs, label):
    import itertools
    ts, cs, diffs, kin = run(kinargs, label)
    if not diffs:
        print(f"[{label}] decompositions identical -- EXACT")
        return True
    phys = [kin['D0'], kin['D1'], kin['Dd'], kin['Dx']]
    seeds = []
    for size in (2, 3):
        for subset in itertools.combinations(phys, size):
            for powers in itertools.product((1, 2, 3), repeat=size):
                if sum(powers) > 7:
                    continue
                seeds.append(dict(zip(subset, powers)))
    rels = ibp_relations(seeds) + md_relations(seeds, phys)
    print(f"[{label}] {len(rels)} null relations", flush=True)
    ok = solve_membership(diffs, rels)
    print(f"[{label}] residual is an exact null combination: {ok}", flush=True)
    return ok


def main():
    A = (R(3), R(4), R(5), R(9,4), R(15,4), R(5,4), R(3,4),
         R(5,3), R(-4,3), R(41,9), R(40,9), R(3,5), R(4,5))
    B = (R(5), R(12), R(13), R(15,4), R(25,4), R(5,3), R(4,3),
         R(5,4), R(-3,4), R(25,7), R(24,7), R(5,13), R(12,13))
    ok = all([verify_point(A, "point A"), verify_point(B, "point B")])
    print("PASS: N1..N18 is exactly the reduced squared amplitude"
          if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
