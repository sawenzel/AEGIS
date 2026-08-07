# TEPEMGENCPP — a C++ port of the TEPEMGEN QED-background generator

Tracking issue: [O2-6340](https://its.cern.ch/jira/browse/O2-6340).

`../TEPEMGEN` generates QED-background e+e- pairs for the ALICE MC chain. Its
numerical core is FORTRAN 77 and loses **all** significance to catastrophic
cancellation in part of the kinematic range, which intermittently makes the
generator fail to initialise. This directory holds the replacement: a templated
C++ core that can be evaluated in `double`, `__float128`, or any user-supplied
arithmetic type.

**Status: in progress.** The Fortran in `../TEPEMGEN` remains the numerical
oracle and stays in the build until the port is validated against it. Nothing
outside AEGIS is touched yet; O2/O2DPG integration is deliberately later work.

## Documents

| | |
|---|---|
| [doc/01-scope.md](doc/01-scope.md) | What is being replaced, what must not change, and why |
| [doc/02-findings.md](doc/02-findings.md) | The numerical diagnosis, and bugs found while reading |
| [doc/03-plan.md](doc/03-plan.md) | Phased plan with the acceptance gate for each phase |

## Building

Everything is built on an x86_64 EL9 machine with CVMFS. The scripts pin the
whole toolchain; see [doc/01-scope.md](doc/01-scope.md#build-environment) for
why that pinning is not optional.

```bash
scripts/env.sh bash -c 'echo $CXX'    # sanity: should be the CVMFS GCC 14.2
scripts/build.sh                      # configure + build AEGIS incl. TEPEMGEN
```

`scripts/build.sh` passes two flags that upstream AEGIS needs but does not set
itself — without them the tree does not configure or compile on a current
system at all. Both are explained in `doc/02-findings.md`.
