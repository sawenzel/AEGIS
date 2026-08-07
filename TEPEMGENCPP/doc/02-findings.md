# Findings

Established by reading `../TEPEMGEN`, the tracking issue O2-6340 and its
attached figure, and by building and running the unmodified generator.

## 1. The cancellation

`Initdiffcross` sets

```fortran
w1xw1 = 1/gamma**2      w2xw2 = 1/gamma**2      w1xw2 = 2 - 1/gamma**2
```

with `gamma = cm_energy/(2*0.938)`. For PbPb that is gamma ~ 2860, so:

| sqrt(s_NN) [GeV] | gamma | y_beam = arcosh(gamma) | (w1xw1*w2xw2)/(w1xw2)^2 |
|---|---|---|---|
| 5020 | 2676 | 8.585 | 4.9e-15 |
| 5360 | 2857 | 8.651 | 3.8e-15 |
| 5500 | 2932 | 8.677 | 3.4e-15 |

The `N1..N18` coefficient polynomials in `Diffcross` mix both orders freely —
`N4` is literally `(-2*w1xw1*w2xw2 + 8*w1xw2**2)`, i.e. `-3.5e-14 + 32`. Terms
spanning fifteen orders of magnitude are summed into `NT` against a double
epsilon of 2.2e-16, leaving on the order of one significant digit before any
kinematic amplification.

The original authors knew. `diffcross.f` line 27:

> To use this program at LHC energies, please make sure that "double precision"
> variables should better be real*16

That instruction was never followed.

## 2. Where it bites — measured, and it is not where we first guessed

### The quad oracle

`gfortran -freal-8-real-16` promotes every `DOUBLE PRECISION` to `REAL*16`, so
the **unmodified** `diffcross.f` can be compiled a second time and evaluated in
quad. That gives a reference curve before any C++ exists. `tools/` builds both
(`*_d`, `*_q`); they cannot share an executable — same symbol names, different
ABI — so sampling and evaluation are separate programs joined by a file.

Sanity check: both report `gamma=2857.14`, `wy=arcosh(gamma)=8.65072`.

### A hypothesis that was wrong

An earlier reading of this code predicted the corrupted window would bracket
`ym = -arcosh(gamma)`, the beam rapidity, because `Diffcross` evaluates
`cosh(ym +/- wy)`. **Measurement refutes that.** Fixing `yp` and scanning `ym`
at `dphi = pi` puts the worst point at `ym ~ -yp`, tracking `yp` and not `wy`:

| yp | -8.65 | -6.0 | -3.0 | 0.0 | 2.0 | 5.0 |
|---|---|---|---|---|---|---|
| worst at ym | +8.63 | +5.83 | +2.93 | +0.05 | -2.13 | -5.05 |

Those particular points also have enormous *relative* error (up to 1e54) for an
uninteresting reason: `yp + ym = 0` with equal pt and `dphi = pi` is the pair
produced at rest, where the cross section has a zero. Relative error is the
wrong metric next to a zero. It is recorded here so nobody rediscovers it and
mistakes it for a catastrophe.

### What actually governs the conditioning

Over 100 000 points sampled from the real generator at production settings,
binning the double-vs-quad relative error:

| \|ym - yp\| | n | mean rel. err | max |
|---|---|---|---|
| 0.00-0.25 | 11686 | **5.8e-2** | 4.2e+2 |
| 0.50-0.75 | 11457 | 2.0e-3 | 4.6 |
| 1.75-2.00 | 6153 | 6.2e-4 | 2.3e-1 |
| > 3.00 | 7326 | **1.8e-4** | 3.1e-2 |

| \|dphi - pi\| | n | mean rel. err | max |
|---|---|---|---|
| < 0.001 | 1736 | **1.8e-1** | 1.8e+2 |
| < 0.01 | 10581 | 5.5e-2 | 4.2e+2 |
| < 0.05 | 15591 | 2.4e-3 | 1.4 |
| < 0.5 | 16891 | **8.8e-5** | 1.7e-2 |

Conditioning degrades monotonically as `ym -> yp` and `dphi -> pi`: a factor
~300 across the rapidity-difference range and ~2000 across the azimuth range.

The mechanism for the azimuthal half is concrete and has nothing to do with
`gamma`. With `ppt = (ppvt, 0)` and `pmt = (pmvt*cos dphi, pmvt*sin dphi)`, the
cross product `axy = x1*y2 - x2*y1` carries a factor `sin(dphi)`, which
vanishes at `dphi = pi`. The `Id0..Id3` denominators
`B = 4u*axy^2 + A` then lose their stabilising term and collapse onto `A`
alone, and `Id3` divides by `B^3`.

### Why this is the worst possible place for it

**The sampler concentrates events exactly on the ill-conditioned locus.**

- `DsdYmY` is a narrow Gaussian in `YmY = Yp - Ye` peaked at 0
  (`parYmY(1)*exp(-8*Y**2)` for `|Y| < 0.18`).
- `DsdPhi` peaks at `Phi = pi`; the only azimuthal cut in `ee_event`
  (`Phi < 0.03 or Phi > 2*pi - 0.03`) removes the *collinear* configuration and
  leaves the back-to-back one untouched.

So 11.7% of sampled events land in `|ym-yp| < 0.25` and 1.7% within
`|dphi-pi| < 0.001`, where the mean relative error is 18%. This is not a rare
corner the generator wanders into; it is where the generator aims.

## 3. Why it ends in abort(), and why the existing guard cannot help

`Diffcross` ends with

```fortran
IF((setzero.EQ.1).OR.(dsigma.LT.0)) THEN
   dsigma=0D0
   badcount=badcount+1
```

The garbage is positive about as often as negative, so this catches only half
of it. The visible noise in the O2-6340 figure is only the surviving positive
half — a log axis cannot show the zeroed points.

A positive outlier then propagates:

```fortran
Xsect2 = Xsect2 + Wtm2         ->  Xsecttot = Xsect2/Nevnt
Dsect2 = Dsect2 + Wtm2*Wtm2    ->  Dsecttot = sqrt(Dsect2)/Nevnt
```

A single weight overestimated by 1e5 contributes ~1e10 to `Dsect2`. The
convergence test in `TGenEpEmv1::CalcXSection`, `err/xSect < eps`, then cannot
be satisfied; the loop runs to `fMaxXSTest = 1e7` and calls `abort()`.

This is the mechanism by which a corrupted point could produce the reported
"convergence issues and initialization failures" — a variance blowup, not a
quadrature failure.

**Measured, and it does not yet reproduce a failure.** Over 100 000 sampled
points (seed 12345, production settings):

| quantity | value |
|---|---|
| points where the guard fired (`dsigma` forced to 0) | 19 (0.019%) |
| points with rel. err > 1e-9 | 96.2% |
| points with rel. err > 1e-3 | 6.6% |
| points with rel. err > 1e-1 | 249 (0.25%) |
| **bias on sum(dsigma)** — the cross section | **-0.072%** |
| **bias on sum(dsigma^2)** — what `Dsect2` accumulates | **-0.80%** |
| largest single absolute error | 40.0, i.e. 0.04% of the whole sum |
| points with `double > 3x quad` and `double > 1` | 4 |

So: the loss of significance is pervasive (96% of events lose seven digits or
more), but for this sample its *aggregate* effect is a 0.07% shift in the cross
section — comfortably inside the 1% tolerance the generator runs with. The 249
worst-relative-error points contribute 0.086% of the integral, because they sit
where `dsigma` is small.

Two things follow, and they should not be conflated:

1. **The guard is inadequate on its own terms.** It fires on 0.019% of points
   while 0.25% are more than 10% wrong — it catches roughly one in thirteen of
   the badly wrong points, because it only ever sees the negative half.
2. **A catastrophic initialisation failure has not yet been reproduced here.**
   The worst outlier found so far is `double = 8.01` against `quad = 0.045` at
   `yp=-6.82, ym=-6.73, dphi=3.14181` — 180x too large, but small against the
   1084 maximum in the same sample, so not enough to poison `Dsect2`. A larger
   hunt is running; until it lands, the abort chain in this section is a
   *mechanism consistent with the report*, not something observed.

Consequence for the port either way: escalate to higher precision on a
*cancellation monitor*

```
NT = sum(Ni) ;  escalate when |NT| / max|Ni| < ~1e-13
```

never on the sign of `dsigma`. The threshold is to be calibrated from the
Phase 2 study, not guessed.

## 4. A second, unreported numerical bug

`ee_init` integrates the *envelope* functions (not `Diffcross`) with CERNLIB
`Dtrint`. Those envelopes are piecewise with derivative kinks:

- `DsdXpX` at `X = 0.6`
- `DsdYmY` at `|Y| = 0.18` and `|Y| = 4.00`

Adaptive subdivision across a C0 kink converges very slowly, hits `KMX=35`,
calls `MTLPRT('D105.2','TOO HIGH ACCURACY REQUIRED')` and returns 0. `ee_init`
then prints an error and **continues**:

```fortran
IF (Xsec1*Xsec2.le.0.or.Ysec1*Ysec2.le.0.) then
   write(*,*) ' Error: insufficient accuracy of XY-cross sections'
endif                       ! prints and carries on
```

`XYsect` becomes 0, every weight becomes 0. This is *not* what O2-6340 reports
and higher precision does not address it; splitting the integration domain at
the known breakpoints does. Fix it anyway.

There is also a likely structural simplification, to be verified before being
relied on: the two triangles per `Dtrint` call union to the square
`[Xmin,Xmax]^2`, and under `u = Xp+Xe, v = Xp-Xe` the inner `v` integral of
`A*exp(-B|v|) + C*exp(-D|v|)` is closed-form. If that holds, the entire
adaptive-triangle machinery reduces to a 1-D adaptive integral that can be
split at the kinks — removing `dtrint.f` and the `MTLPRT` dependency outright.

## 5. Defects in the integral functions

- `Iz0` computes `arg`, guards `IF(arg.lt.0)`, then **recomputes**
  `LOG((tepxx+u+v+s)**2/(4*u*v))` instead of using `arg`. Harmless today but
  the guard and the evaluated expression can drift apart. The guard also
  addresses the wrong hazard: it catches `u*v < 0`, while the unguarded
  blow-ups are `arg -> 0` (log of zero) and `s -> 0` (division by zero).
  Note `s` is provably real for `u,v > 0` by AM-GM since `tepxx >= 0`; that
  invariant should be asserted in the port rather than rediscovered.
- `setzero` is set deep inside `Iz0` but read only at the very end of
  `Diffcross`, so once a point goes bad the code evaluates an entire tree of
  garbage before discarding it. Early return removes both the waste and the
  global mutable state.
- `Id0..Id3` divide by `B = 4u*axy^2 + A` with `axy = x1*y2 - x2*y1`, and `Id3`
  by `B^3`. `axy -> 0` is the near-collinear configuration — exactly where
  `DsdPhi` peaks. `ee_event` already carries a band-aid for this: a hard 30 mrad
  cut, `if (Phi.lt.0.03.or.Phi.gt.2.*pi-0.03) go to 30`. Whether that cut can be
  narrowed once the arithmetic is sound is a physics improvement worth testing.

## 6. Bugs found while building

- **`TGenEpEmv1::SetMinMaxXSTest(int,int)` is declared in the header and
  defined nowhere in the repository.** The shared library links because nothing
  internal calls it; it fails only for an external caller, at link time. It
  matters here because it is the only knob for bounding the cross-section
  estimation loop, so callers are stuck with the hard-coded `1000 / 1e7` and
  therefore with the `abort()` path above. (The undefined copy constructor and
  `operator=` beside it are deliberate — the old ROOT idiom for suppressing
  copying.)
- **Upstream AEGIS does not build with a current gfortran.** MICROCERN's
  `sortzv.F` fails with `Type mismatch in argument 'a': passed REAL(4) to
  INTEGER(4)`; gfortran >= 10 made such mismatches errors. Needs
  `-fallow-argument-mismatch`.
- **`GeneratorParam` does not configure out of the box.** Its
  `find_package(Pythia6)` ignores the `PYTHIA6_ROOT` environment variable under
  CMake policy CMP0144, so `PYTHIA6_LIBRARY` must be passed explicitly.

## 7. Baseline measurement

Unmodified Fortran, at the settings O2DPG actually uses (`y in [-7,7]`,
pt 1..1000 MeV, PbPb 5.36 TeV, `gRandom` seed 12345):

```
X-section = 3.581881e+04 with 3.581845e+02 error after 25689 trials
INIT ok=1   xsec = 35818.8 b   releps = 0.01      (0.73 s)
```

O2DPG expects `QEDXSecExpected['PbPb'] = 35237.5 b`. The baseline is +1.65%
high, about 1.6 sigma on its own quoted error. **This is the number the port
must reproduce or consciously supersede.**

This seed initialised without trouble, which is consistent with the
intermittency in section 2 — a single successful run proves nothing. Measuring
the *failure rate* over many seeds is the first Phase 0 task. `badcount` is
reachable from a C++ harness without modifying the Fortran, by declaring the
common block directly:

```cpp
extern "C" { extern struct { int badcount; } badpar_; }
```

## 8. Cross-repository consequence

`o2dpg_sim_workflow.py` hard-codes
`QEDXSecExpected = {'PbPb': 35237.5, 'OO': 3.17289, 'NeNe': 7.74633}`, the QED
task greps `xSectionQED` out of `qedgenparam.ini` to check against it, and
`--qed-x-section-ratio` feeds the digitiser's QED pile-up rate.

If correcting the arithmetic moves the cross section — and if the old code was
losing digits, it should — that check fails and the amount of QED background in
digitisation changes. This is a coordinated change needing agreement from
whoever owns the QED normalisation, not a follow-up commit.
