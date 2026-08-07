# Scope

## What TEPEMGEN is

2,483 lines, of which 1,679 are Fortran. Small enough to port in full rather
than incrementally wrap.

| File | Lines | Role | Fate |
|---|---|---|---|
| `diffcross.f` | 934 | Exact 5-fold differential cross section: `Initdiffcross`, `Diffcross`, and ten pure integral functions `Iz0-2`, `Id0-3`, `Iv0-2` | port |
| `epemgen.f` | 516 | Sampling layer: envelope parametrisations, rejection / inverse-CDF sampling, weight = exact/envelope | port |
| `dtrint.f` | 229 | CERNLIB D105 adaptive triangle quadrature; used *only* by `ee_init`, 4 calls | delete, see plan Phase 3 |
| `TEpEmGen.{h,cxx}` | 177 | f77 bridge: name mangling, `extern "C"`, `eernd` -> `gRandom` | delete |
| `TEcommon.h` | 32 | `cfortran.h` COMMON-block access | delete |
| `TGenEpEmv1`, `TGenQEDBg` | 537 | ROOT `TGenerator` subclasses | keep, retarget |

External symbols to replace: `Dgauss` (CERNLIB D103, 6 calls, `ee_init` only)
and `MTLPRT` (CERNLIB error printer). Both come from `MICROCERN`, which AEGIS
builds itself (`add_subdirectory(MICROCERN)`) rather than taking from CVMFS —
so there is no external CERNLIB dependency to satisfy, and none to keep.

`eernd` is already C++ (`gRandom->Rndm()`), so the RNG needs no work.

## The compatibility surface

O2 reaches this code through
`$O2_ROOT/share/Generators/external/QEDLoader.C`, driven from
`o2dpg_sim_workflow.py` via `-g extgen`. **Nothing outside AEGIS references the
Fortran or `TEpEmGen`.** The contract is therefore exactly:

- class `TGenEpEmv1` — its `Set*` configuration methods, `Init()`, `GenerateEvent()`
- class `TGenQEDBg` — additionally `SetLumiIntTime()`
- library name `libTEPEMGEN.so`

Those names and signatures do not change. That is what makes the port a
self-contained AEGIS change.

## In scope

- Replace the Fortran numerics with a templated C++ core evaluated in `double`
  by default and `__float128` where precision demands it.
- Keep the arithmetic type a template parameter so Boost.Multiprecision or an
  interval type can be substituted for round-off and stability analysis
  (requested on the tracking issue).
- Remove the COMMON blocks, which are what currently makes the generator
  non-reentrant and thread-hostile.
- Validate against the existing Fortran, which stays in the tree as the oracle.

## Out of scope for this branch

- Any change to O2 or O2DPG. In particular the hard-coded
  `QEDXSecExpected = {'PbPb': 35237.5, 'OO': 3.17289, 'NeNe': 7.74633}` in
  `o2dpg_sim_workflow.py` is *not* touched here — see `02-findings.md` for why
  changing it is a coordinated, cross-repository decision rather than a
  follow-up commit.
- Re-deriving the physics. The published cross section
  ([Phys. Rev. A55 (1997) 396](https://arxiv.org/abs/hep-ph/0112211)) is taken
  as given. Analytically restructuring the terms to avoid the cancellation
  symbolically is noted as a possible future improvement, not attempted here.
- Performance work beyond removing gross redundancy in the integral call tree.

## Build environment

x86_64 EL9 with CVMFS. `scripts/env.sh` pins `CC`, `CXX` **and** `FC` to the
CVMFS `GCC-Toolchain`.

This pinning is load-bearing, not tidiness. On a stock EL9 the system has no
`gfortran` at all, while `/usr/bin` precedes the CVMFS toolchain on `PATH`. The
result is that `gfortran` resolves to the CVMFS GCC 14.2 (nothing else provides
it) but `g++` and `gcc` resolve to the system GCC 11.5 — a silently mixed-ABI
build across the very language boundary this project depends on, with the
Fortran oracle and the C++ port compiled by different compilers. Pinning all
three also means `libquadmath` comes from the same toolchain as the code that
uses it.
