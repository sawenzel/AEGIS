# Integration test: a real o2-sim job

Run 2026-08-07 on alibicompute01, against `O2sim/v20260807-1`
(O2 `daily-20260807-0700-1`, ROOT 6.36.10, AEGIS `v1.5.9-32`).

## Result

The C++ backend runs a complete `o2-sim` job — generator, transport, output —
and agrees with the Fortran on the cross section to two parts in 100 000.

| backend | `xSectionQED` | trials | `Init()` attempts | wall | outputs |
|---|---|---|---|---|---|
| fortran | 34949.1 b | 25347 | 1 | 9.4 s | Hits + Kine |
| **cpp** | **34948.3 b** | 25340 | 1 | 20.6 s | Hits + Kine |

Difference: **-0.0023%**. The +11 s is entirely `Init()`'s cross-section
estimate — 25 340 calls at 532 us is 13.5 s — and *none* of it is event
generation, because `TGenEpEmv1` produces one pair per event and the test
generated 20. See `04-future-work.md` (c): the per-event cost is paid and
discarded, so a weightless mode would remove even this.

Neither backend failed to converge with this seed. That is expected — the
failure is intermittent (two seeds in four at 3M events) and 25 000 trials is
far short of that. This test shows the port *works* in situ; it is not a test
of the bug.

## How it was injected

`QEDLoader.C` does `gSystem->Load("libTEPEMGEN")`, so the library is chosen by
the dynamic loader:

```bash
export LD_LIBRARY_PATH=<build>/TEPEMGEN:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=<src>/TEPEMGEN:<src>/TEPEMGENCPP/include:$ROOT_INCLUDE_PATH
export TEPEMGEN_BACKEND=cpp
```

The library must be built inside the O2sim environment, so it links the ROOT
that `o2-sim` will load. Building against `AEGIS::v1.5.9-9` gives ROOT 6.36.04
where O2sim has 6.36.10.

The toolchain is not a free choice either: the *shipped* `libTEPEMGEN.so` and
O2's own libraries both link GCC-Toolchain 14.2's libstdc++ (`GLIBCXX_3.4.32`),
even though a bare `g++` on the machine resolves to the system 11.5 by PATH
order. `scripts/env_o2sim.sh` pins `CC`/`CXX`/`FC` accordingly.

## The trap: `ROOT_INCLUDE_PATH` is not optional

**Setting only `LD_LIBRARY_PATH` produces silent memory corruption, not a
link error.**

`TEpEmGen` gained two members (`fBackend`, `fSampler`), which changes the ROOT
class layout. `QEDepem.C` is interpreted by Cling, which parses the *header* it
finds on `ROOT_INCLUDE_PATH` — the CVMFS AEGIS one, with the old layout — while
calls dispatch into the *new* library. Member offsets then disagree, and
`SetPtRange` writes somewhere harmless-looking.

What that looked like in practice:

- with the C++ backend: `init()` refused with `ptMin >= ptMax`, because the pt
  range arrived as 0 and `log10(0) = -inf` on both ends;
- with the Fortran backend: `ee_init` was handed a pt range of zero and the
  rejection sampler ran for **20 minutes** without producing an event.

That second one is worth dwelling on, because it nearly became a false
positive. A generator hanging for twenty minutes in `Init()` is exactly what
O2-6340 describes, and it would have been easy to report it as the bug
reproducing in production. It was not. It was this injection error. With
`ROOT_INCLUDE_PATH` set correctly the same run converges in 9.4 s.

This is only a hazard for *patched-library* workflows. A normal aliBuild
deployment installs headers and library together, so they cannot disagree. But
patching a library into a CVMFS stack is common ALICE practice, and this class
of failure is silent.

If the port is ever to be dropped in without rebuilding consumers, the two
members must move out of the class — into file-static state — so the ROOT
layout is unchanged. That is a real design constraint, not a preference.

## A second bug this exposed

The abort path printed its reason with `printf`, and nothing appeared in the
log. `stdout` is *fully* buffered when redirected to a file, which every batch
job does, so the message died in the buffer when `abort()` fired. The job
aborted with no explanation at all.

Now `fprintf(stderr, ...)` plus an explicit `fflush`. Any diagnostic that
precedes an `abort()` or a crash has to be unbuffered or it is not a
diagnostic.
