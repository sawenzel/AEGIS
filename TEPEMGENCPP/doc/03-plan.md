# Plan

Each phase ends in a working, buildable tree with a stated acceptance gate. The
Fortran in `../TEPEMGEN` stays in the build as the numerical oracle until
Phase 4; it is not deleted early.

## Design decisions

### Template on the scalar type

```cpp
template <typename T> struct Vec2 { T x, y; };
template <typename T> struct PhysParams { T gamma, beta, m, w1xw1, w1xw2, w2xw2, wl, wy; };
template <typename T> T diffCross(const PhysParams<T>&, T ptPlus, T yPlus,
                                  T ptMinus, T yMinus, T dphi, bool& ok);
```

Explicitly instantiated for `double` and `__float128`. Boost.Multiprecision and
interval types then work by inclusion, without Boost becoming a build
dependency.

Two disciplines make this real rather than decorative:

1. **Call math functions through ADL, never qualified.** `using std::sqrt;
   sqrt(x);` — `std::sqrt(x)` on a Boost multiprecision type silently converts
   to `double` and discards the entire point of the exercise.
2. `__float128` is not a class type, so it needs a small shim mapping to
   `sqrtq`/`logq` from `<quadmath.h>`. Budget for it explicitly.

Restrict the core to a documented operation set — `+ - * /`, `sqrt`, `log`,
`cosh`, `sinh`, `abs`, comparison — so interval arithmetic stays viable.

### Remove the COMMON blocks

| Fortran | C++ |
|---|---|
| `COMMON /PHYSPARAM/` | `PhysParams<T>` by `const&` |
| `COMMON /SETZPARAM/ setzero` | `bool& ok`, early return |
| `COMMON /badpar/ badcount` | counter on the generator object |
| `COMMON /eevent/`, `/eepars/` | members of the sampler class |

This is what buys thread-safety and reentrancy: today two `TGenQEDBg`
instances in one process would silently corrupt each other.

### Mixed precision, not blanket quad

`__float128` is software-emulated; expect 10-50x slowdown. `Diffcross` runs
once per event and `CalcXSection` burns up to 1e7 events, so blanket quad would
make initialisation unusable. Evaluate in `double`, escalate per point on the
cancellation monitor of `02-findings.md` section 3.

### Transcribe the polynomials mechanically

`N1..N18` are machine-generated symbolic output — `N12` alone is 77 lines of
terms. Hand-transcription is not an option: one sign error produces a plausible
wrong cross section that survives every smoke test. Translate by script, and
verify by term count per `N_i` plus random-point evaluation against the Fortran
*before* any refactoring.

The `Iz*/Id*/Iv*` call sites contain deliberate argument permutations that are
trivially fat-fingered and must be transcribed mechanically too:

```fortran
Iz2:  D = (Iz1(x,v,u)/pi)*...     ! (v,u) swapped
Id2:  ... - G1*Iz1(mx,v,u) ...    ! negated x, swapped u/v
Id3:  ... - D2*Iz1(mdyx,w,v) ...  ! negated difference, swapped
```

### Modernity budget

**In:** `constexpr` tables replacing `DATA` statements, `std::array`, structs
over raw `T[2]`, free functions in a namespace, ALICE clang-format and naming,
status returns over global flags.

**Out:** expression templates, policy-based design, virtual dispatch in the
math, a class hierarchy for the integrals, custom allocators, Boost as a hard
dependency. Inheritance stays only where ROOT's `TGenerator` mandates it.

## Phases

### Phase 0 — oracle and harness — **done**

- Unmodified Fortran built and running; baseline in `02-findings.md` section 7.
- **Quad oracle without a port**: `gfortran -freal-8-real-16` recompiles the
  same `diffcross.f` in `REAL*16`. `tools/build_tools.sh` builds both.
- `tools/dump_points` samples from the real sampler (`ee_init_`/`ee_event_`)
  and records the double cross section; `tools/compare_precision` re-evaluates
  in quad and reports the biases on `sum(dsigma)` and `sum(dsigma^2)`.
- `tools/scan_rapidity_{d,q}` reproduces the O2-6340 figure.
- Localisation **refuted and replaced**: it is `|ym-yp|` and `|dphi-pi|` that
  govern conditioning, not `arcosh(gamma)`. See `02-findings.md` section 2.

Still open from this phase, carried forward: spike the `Dtrint` analytic
reduction (needed by Phase 3, not by Phase 1).

### Phase 1 — `diffcross.f` to templated C++, `double` only

Mechanical transcription; COMMON blocks eliminated.

*Gate:* C++ `double` matches golden at ~1e-12 relative. Bit-identical is not
achievable across expression reordering and libm differences — the tolerance is
agreed up front and justified, not tuned until it passes.

### Phase 2 — `__float128`, cancellation monitor, precision study

Quad instantiation and the quadmath shim; escalation implemented and its
threshold calibrated; double-vs-quad evaluated over the golden points to
quantify which kinematic regions lose how many digits.

*Gate:* the O2-6340 figure regenerated from the new code as a committed
regression test.

### Phase 3 — `epemgen.f` sampler, quadrature replaced

Port the samplers **preserving the exact order of `eernd` consumption**. With a
seeded `gRandom` this yields event-by-event identical output against the
Fortran — a far stronger gate than comparing histograms, easy to lose and hard
to recover. Design for it deliberately.

Replace `Dgauss` with a self-contained templated adaptive Gauss-Kronrod
(G7-K15, ~60 lines). Do **not** use GSL or `ROOT::Math::Integrator`: both are
`double`-only, which breaks the precision templating. This is the case where
hand-rolling is the smaller dependency.

*Gates, in order:* (i) event-identical under fixed seed with double math;
(ii) total cross section against the 35818.8 b baseline and O2DPG's 35237.5 b;
(iii) KS/chi2 on the y, pt and dphi distributions, old vs new.

### Phase 4 — remove the Fortran

Delete `*.f`, `TEpEmGen.*`, `TEcommon.h`, the `cfortran.h` dependency and the
`WIN32` name-mangling block. Drop `Fortran` from `project()` and `MICROCERN`
from the link line. Check whether the `pythia6` link and `PYTHIA6_ROOT`
`link_directories` in `TEPEMGEN/CMakeLists.txt` are still needed — they look
like copy-paste, and if so this target stops depending on PYTHIA6 entirely.

Net effect: **TEPEMGEN no longer needs a Fortran compiler.**

`TGenEpEmv1` and `TGenQEDBg` keep their names and signatures, so `QEDLoader.C`
and O2DPG need no change.

## Open decisions

1. **Bit-compatibility policy.** The Phase 3 event-identity gate holds only
   while the C++ computes in `double` with the original expression ordering.
   Once escalation to `__float128` engages, output legitimately differs — that
   is the point. The acceptance criterion must switch from "identical" to
   "statistically equivalent and different in a defensible direction" at a
   defined moment, agreed in advance rather than argued about later.
2. **The cross-section constants.** See `02-findings.md` section 8. Needs a
   named owner before Phase 3 lands.
