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

**Measured over 100 000 sampled points** (seed 12345, production settings):

| quantity | value |
|---|---|
| points where the guard fired (`dsigma` forced to 0) | 19 (0.019%) |
| points with rel. err > 1e-9 | 96.2% |
| points with rel. err > 1e-3 | 6.6% |
| points with rel. err > 1e-1 | 249 (0.25%) |
| bias on sum(dsigma) | -0.072% |
| bias on sum(dsigma^2) | -0.80% |

At that scale the loss of significance is pervasive but its *aggregate* effect
is small, and the guard is already shown to be unsound on its own terms: it
fires on 0.019% of points while 0.25% are more than 10% wrong, catching about
one badly-wrong point in thirteen, because it only ever sees the negative half.

**CONFIRMED at 12 million events** (4 seeds x 3M). The catastrophe is real,
it is intermittent, and it is exactly the mechanism above:

| seed | xsec bias | **variance bias** | worst outlier (double vs quad) |
|---|---|---|---|
| 777 | +0.016% | -0.026% | 2.98e6 vs ~0.036 |
| 90210 | +2.41% | **+841%** | 75508.7 vs 0.0461 |
| 5150 | +2.37% | **+598%** | 61583.3 vs 0.0361 |

A single point returning 75508.7 where the truth is 0.046 -- **1.6 million
times too large, and 21x the largest legitimate cross section anywhere in that
3M-event sample** -- inflates `sum(dsigma^2)` by 841%. That is the variance
estimate `Dsecttot` is built from, `err/xSect < eps` then cannot be satisfied,
and `CalcXSection` runs to `fMaxXSTest = 1e7` and calls `abort()`.

Every one of these outliers is **positive**, so the Fortran's sign guard never
sees it. Whether a given run draws one is what "sometimes" means: two of four
seeds did.

**The catastrophic points share a signature**, and it is a corner rather than a
generic point: both leptons at *minimum* pt (~1 MeV, the `ptMin` cut), both
near the *rapidity edge* (|y| 6.3-6.9, against the |y| < 7 cut), and
back-to-back in azimuth. They are kept as a regression fixture in
`tools/worst_points.h`. Rounding their coordinates to ten digits does not
remove the catastrophe, so this is a small ill-conditioned region, not a
knife-edge.

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

## 9. Traps met while porting

Recorded because each was silent, or nearly so, and each would have produced a
wrong cross section that still compiled and still looked plausible.

### A tab in the first six columns is a NEW statement

`diffcross.f` line 778 begins with five spaces, a TAB, then `Id3= -( ...`.
Under fixed-form rules with the GNU tab extension, a tab in columns 1-6 means
the next character is column 7 — a new statement — unless it is a digit 1-9,
which marks a continuation. A lexer that only tests "column 6 is non-blank"
reads the tab as a continuation and silently glues

```fortran
      E11= tepxx
```

onto the 40-line `Id3` assignment that follows. Here it happened to produce
code that would not compile; it could as easily have produced code that did.
`tools/f2cpp.py` handles the tab rule and additionally rejects any RHS
containing a stray `=`, which is the fingerprint of this class of merge.

### A decimal literal is parsed at double, then widened

`T(0.938)` builds a **double** literal and only then converts to `T`. gfortran
parses `0.938D0` at the full width of `T`. Since 0.938 is not a dyadic
rational, the two differ by ~1e-17 relative — invisible in the double build,
and enough to dominate the quad comparison. Fixing it improved quad agreement
by a factor of 26 000 (mean 7.4e-17 -> 2.8e-21).

Note the criterion is **binary representability, not digit count**. An earlier
guard rejected literals with more than 17 significant digits and let `0.938`
straight through. Every decimal literal from the Fortran is now emitted as an
exact rational, `T(938)/T(1000)`, which is correct for every `T` and needs no
compiler flags.

### The Fortran's own pi is only good to 25 digits

`PARAMETER (PI=3.141592653589793238462643D0)` — true pi continues
`...38327950288`. Harmless in double, but in `REAL*16` it caps the Fortran at
~1e-25. The port uses `acos(-1)`, evaluated at `T`'s own precision, and is
therefore *more* accurate than the oracle it is validated against. Confirmed
not to be the cause of the residual disagreement by rebuilding with the
truncated constant (`-DTEPEMGEN_FORTRAN_PI`) and observing no material change.

### The cancellation monitor must be formed before normalisation

`|NT| / max|Ni|` has to be taken at the moment of summation. `Diffcross` then
rescales `NT` by `4/beta^2`, by `(2pi)^-4`, by `1/4` and into kbarn — about
2.5e-4 in total — without touching the `Ni`. Taking the ratio afterwards
divides a rescaled numerator by an unscaled denominator and understates
survival by that factor, which made every sampled point look catastrophically
ill-conditioned and put the whole distribution in one histogram bin.

## 10. Phase 2: what precision is actually needed, and what it costs

Measured with `tools/calibrate`, `tools/bench` and `tools/strategy` on 100k
points sampled from the real generator at production settings.

### Cost of each arithmetic type

| type | eps | us/call | vs double |
|---|---|---|---|
| `double` | 2.22e-16 | 11.2 | 1x |
| `long double` (x87) | 1.08e-19 | 27.5 | **2.5x** |
| `__float128` | 1.93e-34 | 1080 | **97x** |

`__float128` is software-emulated. Quad throughout turns the 25 689-trial
initialisation from 0.3 s into 28 s, and the hard-coded `fMaxXSTest = 1e7`
worst case into three hours.

### Effect on the quantities the generator uses

| strategy | xsec bias | variance bias | escalated | us/call |
|---|---|---|---|---|
| double only | **-0.073%** | **-0.796%** | 0% | 11.1 |
| long double only | +0.00001% | +0.00006% | 0% | 26.3 |
| long double -> quad @1e-13 | +0.00000% | +0.00001% | 4.8% | 75.6 |
| long double -> quad @1e-14 | -0.00001% | -0.00009% | 1.1% | 37.6 |
| double -> quad @1e-7 | -0.00000% | -0.00000% | 94.0% | 980.7 |

**The naive reading of the issue — "use `__float128`" — is both too much and
too little.** Too much because quad throughout is 97x. Too little because
escalating *from double* does not work: the monitor never exceeds 1e-4 and
about half the integral comes from points below 1e-12, so holding the error
under 1e-3 needs a 1e-7 threshold, which escalates 94% of points and recovers
almost none of the speed.

Simply moving the working type to `long double` improves the cross-section bias
by a factor of 7000 for 2.5x the cost. That is the single highest-value change
in this whole exercise, and it needs no adaptive machinery at all.

### The monitor has a blind spot

`survival = |NT| / max|Ni|` is **necessary but not sufficient**, for two
reasons found by measurement rather than by design:

1. **It is self-referential.** Each type computes survival from its own terms,
   so a type whose terms are already corrupted reports an over-optimistic
   value. At `ptp=ptm=1 MeV, yp=0, ym=0.05, dphi=pi`:

   | | value | survival |
   |---|---|---|
   | double | 8.0e7 | 3.5e-11 |
   | long double | 4.4e4 | 1.9e-14 |
   | quad (truth) | 5.5e3 | 2.4e-15 |

2. **It cannot see cancellation inside `Iz/Id/Iv`.** The `Id0..Id3`
   denominators `B = 4u*axy^2 + A` collapse as `dphi -> pi`, and the final-sum
   monitor is blind to that. Note the point above is the *peak* of the scan,
   not a zero — this is not the relative-error-near-a-zero artefact of
   section 2.

Consequently the default threshold (1e-13) is chosen **empirically**: it is the
largest value that keeps a deliberately worst-case scan (`dphi = pi`, 400
points across the full rapidity range) within 1e-4 of the quad reference. At
1e-14 that scan still leaves one point 8x wrong. The `eps/survival` heuristic
would have suggested 1e-14 was ample; it is not.

### Production cost, and how to pay for it

`TGenQEDBg` generates `Poisson(fPairsInt)` pairs per event, and `fPairsInt` is
~4200 at PbPb luminosity with a 20 us integration window. With
`NEventsQED ~ 28000` that is ~1e8 `Diffcross` calls per production — so the
per-call cost is not academic:

| | per call | ~1e8 calls |
|---|---|---|
| double (today) | 11 us | ~22 min |
| long double | 27 us | ~52 min |
| long double -> quad @1e-13 | 76 us | ~2.5 h |

The way to pay for this is the redundancy already noted in the plan: `Iv2` ->
`Id2` -> `Id0`/`Iz0`/`Iz1` recomputes the same `Iz` values on the same
arguments many times per point. Computing each `{Iz0,Iz1,Iz2}` triple once and
passing it down is a mechanical change with no effect on the result, and it is
the obvious place to recover the factor the wider type costs. **Not yet done**
— it must come after the Phase 1 equivalence gate, never before.

## 11. Phase 3 spike: the Dtrint reduction holds, and Dtrint is the loose one

`include/Envelopes.h` + `include/Quadrature.h` + `tools/test_norm`.

The two triangles per `Dtrint` call do union to the square `[lo,hi]^2`, and
under `u = Xp+Xe`, `v = Xp-Xe` (Jacobian 1/2) the domain becomes a diamond of
half-width `L(u) = min(u-2lo, 2hi-u)`. The inner `v` integral is closed-form
for both envelopes -- exponentials for `DsdXmX`, error functions plus an
exponential tail for the piecewise `DsdYmY`. So

```
I = 1/2 * integral_{2lo}^{2hi} f(u) * G(L(u)) du
```

A 2-D adaptive triangle integration collapses to one 1-D integral, whose panel
edges can be placed on the known kinks so nothing has to converge across one.
**`dtrint.f` (229 lines of goto-driven CERNLIB) and the `MTLPRT` dependency
can both go.**

Numbers at production settings (`Xmin=0, Xmax=3`; `Ymin=-7, Ymax=7`):

| | reduction | Dtrint | rel diff | Dtrint's own budget |
|---|---|---|---|---|
| XsecX | 2.012153900694043 | 2.011971250182245 | 9.1e-5 | 7.5e-5 |
| XsecY | 224.9720738805846 | 224.9571343935239 | 6.6e-5 | 5.0e-5 |

The disagreement is **Dtrint's**, not the reduction's. `ee_init` passes
`Eps = 0.00005` and `Dtrint` stops at `|SUM0-SUM| <= EPS*(1+|SUM|)`, so ~1e-4
is all it was asked for. Brute-force Simpson on a fine uniform grid — no
adaptivity to go wrong — confirms the reduction for XsecX:

```
n=2000   2.0126131
n=8000   2.0122073
n=20000  2.0121590      -> reduction 2.0121539, not Dtrint's 2.0119713
```

The reduction's own quadrature converges to 2.4e-12 absolute in 12 panels.

So the existing normalisation is accurate to ~1e-4, comfortably inside the 1%
the generator runs with. Replacing `Dtrint` is therefore not a correction of a
wrong number; it removes a **failure mode** — the silent `return 0` on
non-convergence that leaves `ee_init` printing an error and carrying on with
`XYsect = 0`, zeroing every event weight.

A methodological note worth keeping: a nested adaptive 2-D integration was
tried as the independent check and had to be discarded. The inner integral
carries its own error, so the outer integrand is effectively noisy at that
level; with the outer tolerance tighter than the inner it never converges, and
the X case returned 2.238 against 2.012 from two independent methods. Either
loosen the outer tolerance well below the inner, or use a non-adaptive method
as the tie-breaker.


## 12. The monitor needs a second, direct trigger

Section 10 recorded that the cancellation monitor is necessary but not
sufficient. The 12M-event hunt produced a real counter-example rather than a
constructed one, and it forced a design change.

At the seed-90210 catastrophic point the three tiers report:

| | value | survival |
|---|---|---|
| double | 75508.7 | 1.21e-08 |
| long double | **-1.33** | 2.11e-13 |
| quad (truth) | 0.0461 | 7.30e-15 |

`long double` reports survival 2.11e-13 -- comfortably above the 1e-13
threshold that was the default, so **no escalation happened** -- while the true
survival is 7.3e-15 and its answer is wrong by a factor 29 and has the wrong
sign. The monitor is computed from terms that are themselves already
corrupted, so it is most over-optimistic exactly where it matters.

The fix is a second trigger that **measures** the instability instead of
predicting it: evaluate the same point in a *narrower* type and escalate if the
two disagree by more than 1e-6 relative. Two types that disagree cannot both be
right, and unlike the monitor this cannot be fooled by corrupted terms. On the
seed-777 point that is the only thing that catches it -- its `long double`
survival is 9.79e-10, three decades above any sane threshold, while its answer
is -5065 against a true 0.036.

Both triggers are kept. `tools/test_worst` asserts the adaptive path returns
the quad answer at all four known catastrophic points and reports which trigger
fired; two are caught by the monitor, two only by the disagreement check.

Cost: escalation rises from 4.7% to 47.7% of points. That is a deliberate trade
-- CPU is not the binding constraint here, and the alternative is a generator
that intermittently aborts.

## 13. Cross-check against the literature

Two things were asked for here and they are not the same. Reproducing the
issue's figure with the ported code is a *numerical* check and it is done
(`doc/o2-6340-old-vs-new.png`, generated by `tools/make_plot.py`, with the
original Fortran in `REAL*16` as ground truth so the reference does not depend
on the code being validated). Checking the *physics* against the cited paper is
a different and much larger job, and what follows is deliberately less than
that.

### What the citation actually is

The issue cites [hep-ph/0112211](https://arxiv.org/abs/hep-ph/0112211), which
is a 125-page review — Baur, Hencken, Trautmann, Sadovsky and Kharlov,
*Phys. Rept.* **364** (2002) 359 — not the source of the formula in the code.
`diffcross.f` cites Alscher, Hencken, Trautmann and Baur,
*Phys. Rev.* **A55** (1997) 396, and the generator as a whole cites the ALICE
note Hencken, Kharlov, Sadovsky, ALICE-INT-2002-27. The review does not
tabulate a PbPb e+e- total that can be lifted directly.

### What was checked instead

`tools/total_xsec` integrates the differential cross section over the full
phase space, following the recipe in the header of `diffcross.f`:

```
sigma = (Z*alpha)^4 * 2pi * Int dsigma * p+ dp+ * p- dp- * dy+ dy- ddphi
```

| configuration | sigma |
|---|---|
| pt > 0.02 MeV, \|y\| < 12 | 158 +- 10 kb |
| pt > 0.002 MeV | 175 +- 15 kb |
| pt < 20 GeV | **196 +- 18 kb** |
| \|y\| < 16 | 470 +- 290 kb (62% error, not usable) |

against

| reference | value |
|---|---|
| Racah Born asymptotic, point-like, `L = ln(gamma^2)` | 223 kb |
| the same formula with `L = ln(4 gamma^2)` | 289 kb |
| literature, coherent e+e- in UPC PbPb at the LHC | ~200 kb |

The estimate rises monotonically as each truncation is relaxed and settles near
the literature value. The earlier deficit was phase space, not physics.

### What this does and does not establish

It establishes that the overall normalisation, the power of `Z*alpha` and the
MeV/barn unit chain are right. Those are the errors a port or a 25-year-old
code plausibly has, and they would show as factors of a thousand, not thirty
percent.

It does **not** establish that the 5-fold differential cross section correctly
encodes the amplitudes of Alscher et al. Confirming that means re-deriving it,
which is the first checkbox on O2-6340 and remains open. Two further caveats
worth stating plainly:

- The Racah figure is convention-dependent: the same formula gives 223 kb for
  `L = ln(gamma^2)` and 289 kb for `L = ln(4 gamma^2)`. It is quoted here with
  the convention attached rather than as a bare number. It is also a Born,
  point-like result, so the physical value should sit *below* it — Coulomb
  corrections and the nuclear form factor both reduce it.
- The integration is a log-uniform Monte Carlo over an integrand spanning ~20
  decades, which is heavy-tailed. For heavy tails the naive standard error is
  unreliable *and biased low*, because the rare large contributions are exactly
  the ones most often missed. `total_xsec` therefore reports the largest single
  sample as a share of the sum and flags the estimate as a lower bound when one
  sample carries a noticeable part of it. The `|y| < 16` row above, with its
  62% error, is what that failure looks like; it is shown rather than dropped.
