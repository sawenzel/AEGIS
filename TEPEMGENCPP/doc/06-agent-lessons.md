# Lessons for the tooling

Feedback from the session that produced this branch (2026-08-07/08), recorded
here because it is about how the work was done rather than about the code, and
because the improvements belong in the `o2dpg-agent` plugin rather than in
AEGIS.

## 1. A credential failure was worked around instead of fixed

**What happened.** An `o2-sim` job aborted with

```
Alien Token Check failed - Please get an alien token before running with
https CCDB endpoint, or alice-ccdb.cern.ch!
```

The response was to switch `--field ccdb` to `--field 5`, removing the CCDB
dependency, and carry on.

**Why that was wrong.** The GRID certificate was on the very host the job ran
on (`~/.globus/usercert_2025.pem`), `alien-token-init` was on `PATH` inside the
loaded environment (`xjalienfs/1.7.0-15/bin/alien-token-init`), and the token
in `/tmp/tokencert_0.pem` was simply three weeks stale. GRID access is a
first-class capability of the agent plugin — `check_setup` and `selftest`
verify it, and it had been configured and verified earlier in the same session.

Worse than the missed shortcut: dodging CCDB silently narrowed the test. The
`--field 5` run was no longer the production configuration, and the change in
what was being validated went unremarked.

**Improvement for the plugin.** When a tool reports a token, certificate or
credential failure, the agent should recognise it as a capability the
installation already has and offer to refresh it — `alien-token-init` for
GRID, `sso_cookie.py` for CERN SSO — rather than reach for a flag that avoids
the dependency. Where refreshing needs a passphrase or 2FA the command must be
handed to the user, but that is still the right first move. Falling back to
avoiding the dependency is acceptable only after the user declines, and the
narrowed scope must be stated.

## 2. Injecting a local library into a CVMFS stack needs the environment adjusted

**The authority is aliBuild + alidist modulefiles.** That is what sets the
ALICE environment up correctly. A locally built library dropped into that stack
is the exception, and three things have to be made to match — only one of which
is obvious:

| | what | if you get it wrong |
|---|---|---|
| library | `LD_LIBRARY_PATH` | the CVMFS library is used and nothing changes |
| **headers** | **`ROOT_INCLUDE_PATH`** | **silent memory corruption** |
| toolchain | `CC`/`CXX`/`FC` and the ROOT version | ABI mismatch, or subtle breakage |

The middle row is the trap. Cling parses headers at runtime for interpreted
macros such as `QEDepem.C`. It found the CVMFS headers, with the old class
layout, while calls dispatched into the new library — so `SetPtRange` wrote at
the wrong offset and the pt range arrived as zero. That produced a generator
that sat in `Init()` for twenty minutes, which is *exactly* the symptom
O2-6340 describes and was nearly reported as the bug reproducing in production.
It was not the bug. Full account in `05-o2dpg-test.md`.

The toolchain row is not a free choice either: shipped ALICE libraries link
GCC-Toolchain's libstdc++ (`GLIBCXX_3.4.32` for 14.2) even though a bare `g++`
on an EL9 host resolves to the system 11.5 by `PATH` order. And the library
must be built against the ROOT the consumer will load — `alienv setenv
O2sim/<tag>`, not a package that pulls an older ROOT. `scripts/env_o2sim.sh`
in this directory does both.

**Improvement for the plugin.** A helper that, given a locally built package
and a target environment, sets `LD_LIBRARY_PATH` and `ROOT_INCLUDE_PATH`
together, pins the toolchain from the shipped libraries (`readelf`/`ldd` on a
reference `.so` rather than a guess), and then *verifies* the local library is
the one loaded before any physics result is trusted. A deliberately invalid
configuration value that only the new code rejects makes a good probe — that is
how the override was confirmed here.

More generally: in a patched-stack scenario, an unexplained hang or a
nonsensical parameter value should be treated as an environment mismatch first
and a real bug second. The reverse order costs a lot and can manufacture a
false positive that looks like exactly the bug you were hunting.

## 3. A consequence for this port specifically

A class carrying a ROOT dictionary (`ClassDef`) cannot gain members and remain
drop-in — that changes the layout and every consumer must be rebuilt. If
`libTEPEMGEN` is ever to be swapped into an existing stack without rebuilding
O2, the `fBackend` and `fSampler` members added to `TEpEmGen` must move to
file-static state. In a normal aliBuild deployment headers and library ship
together and cannot disagree, so this is a constraint on patching workflows
rather than on the release.
