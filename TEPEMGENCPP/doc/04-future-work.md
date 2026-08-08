# Future work

Recorded 2026-08-07. Everything here is beyond the port itself, which is done
and validated. Items (a)-(d) form one coherent piece of work and are listed in
the order they should be settled, because the answer to (a)/(b) decides whether
(d) is necessary or merely nice.

**Decision 2026-08-08 (operator):** the weights are not needed. The QED events
are treated as noise for the detector response — they never enter a physics
analysis through the AODs — so (a)/(b) are resolved as "no", (c) the explicit
weightless mode is the recommended route, and (d) remains attractive even in a
weighted future because a sampler prebuilt from the exact cross section would
make the question moot. With that, the remaining items — a Gram-stable
restructuring of the terms (see 07-derivation.md §10) and the model-level
correction terms (Coulomb corrections, nuclear form factor) — are explicitly
**beyond the intent of O2-6340**: research nice-to-haves that would round the
project off for the internal note, not blockers for production use.

## The finding that motivates all of it

In `TGenEpEmv1::GeneratePair` the per-event weight is filled by
`TEpEmGen::GenerateEvent` and **never read again**. The `TParticle`
constructors take no weight. So at generation time the exact differential cross
section is evaluated at full cost and the result is discarded.

Consequences, and which applies is a physics decision, not a coding one:

- The generated kinematics follow the **envelope parametrisation**
  (`DsdYpY x DsdYmY x DsdXpX x DsdXmX x DsdPhi`), not the exact cross section.
  The exact calculation only fixes the normalisation, through
  `Init()` -> `CalcXSection()` -> `fPairsInt`.
- The doc comment above `GeneratePair` ("the event weight assigned to each
  track") is inherited AliRoot prose. AliRoot's generator machinery could carry
  track weights; this standalone O2 class cannot.
- The envelope is a fit, not an identity. The `DsdPhi` DATA statement records
  `N(Wt>20) = 573` per 1e7 events, so the weight distribution has a real tail
  and the envelope is not a faithful stand-in for the true shape.

## (a) Does the AOD data model carry an MC weight?

Check whether O2's MC tables (`McParticles`, `McCollisions`, and the generator
header) have a per-particle or per-event weight field that could carry `Wtm2`
through to analysis. If there is one, applying the weight is a plumbing job. If
there is not, adding one is a data-model change with a much wider blast radius,
and (d) becomes the more attractive route.

## (b) Apply the weight, if the data model allows it

Only worth doing if (a) says the field exists **and** the physics owner
confirms the events are supposed to follow the exact cross section rather than
the envelope. Applying weights changes every downstream QED-background
distribution, so it is not a silent fix.

Note the interaction with precision: if weights are applied, they must be
*accurate*, and the whole adaptive-precision machinery becomes load-bearing at
runtime. If they are not applied, precision only has to be right during
`Init()`.

## (c) Offer an explicit weightless mode

Independently of (a) and (b), the generator should be able to skip the cross
section evaluation once `Init()` has converged. Today that cost is paid on
every event and thrown away.

Measured, per call:

| implementation | us/call | vs Fortran double |
|---|---|---|
| Fortran double (production today) | 10.0 | 1x |
| C++ double | 11.0 | 1.1x |
| C++ long double | 25.9 | 2.6x |
| **C++ adaptive (shipped)** | **531.8** | **53x** |
| Fortran REAL*16 | 978.9 | 98x |
| C++ `__float128` | 1035.8 | 103x |

`TGenQEDBg` draws `Poisson(~4200)` pairs per event and a production makes ~1e8
calls, so 53x turns roughly 22 minutes into ~19 hours. A weightless mode makes
that cost vanish entirely and confines the precise arithmetic to `Init()` --
25689 calls at 532 us is **14 seconds, once**.

This should be an explicit, named mode rather than an accident of the weight
being unused, so that the choice is visible in the configuration.

## (d) Precompute the sampler once; serve it from CCDB

The current design samples from a crude envelope and corrects with a weight.
The modern equivalent is to build a good sampler **once**, offline, from the
exact cross section, and then sample unweighted at runtime. ROOT's `TFoam` does
exactly this -- adaptive cellular subdivision of an arbitrary N-dimensional
density, an exploration phase, then cheap generation of unweighted events, and
it serialises to a ROOT file.

What makes this far more tractable than "a 5-D table" sounds:

- **Z drops out completely.** It enters in exactly one line, as a `(Z*alpha)^4`
  prefactor, and `Diffcross` does not even take it as an argument. One table
  therefore serves PbPb, XeXe, OO, NeNe and anything else. Verified: the cross
  section scales as Z^4 to 0.0000%.
- Only the beam energy changes the shape, and production uses a handful of
  energies. Measured sensitivity: +2.35% from 5.02 to 5.36 TeV.
- Production cuts are fixed (`y` in [-7,7], pt 1..1000 MeV), so the domain is
  fixed too.
- The expensive quad arithmetic is paid once, offline: ~1e5-1e6 evaluations,
  ten to twenty minutes.

What needs care:

- The density spans ~20 decades and has a sharp ridge at `dphi -> pi`. Build
  the foam in transformed coordinates -- `log pt`, and something like
  `log|dphi - pi|` -- so the structure is smooth where cells are laid down.
  Otherwise the adaptation spends its budget resolving a coordinate artefact.
- The `dphi`-rapidity correlation is not factorisable, which rules out plain
  VEGAS. Cellular/foam handles it; that is the reason to prefer it.
- **A stored table must be versioned and bound to (energy, cuts)**, hashed into
  the object name or held in CCDB with those as metadata. A silently-wrong
  table is a worse failure mode than a slow generator, because nothing crashes
  and nothing warns.

A lighter variant, if a smaller change is wanted: keep the envelope sampling
and tabulate only the **ratio** `dsigma/envelope`. That ratio is near-flat by
construction -- it is what the envelope fit is for -- so a coarse grid
interpolates it well and the interpolation error is second order rather than
spanning decades.

## Writing it up

An **ALICE internal note** is clearly warranted and is the natural successor to
Hencken, Kharlov and Sadovsky, ALICE-INT-2002-27, which documents this
generator. The port, the failure mode and the validation all belong somewhere
citable.

An external paper is a different question and is not ready. The transferable
content is numerical rather than physical:

- A 1997 generator loses essentially all significance at LHC energies, and --
  the interesting part -- the importance sampler concentrates events precisely
  on the ill-conditioned locus. `DsdYmY` peaks at `YmY = 0` and `DsdPhi` at
  `Phi = pi`, which is exactly where the propagator denominators collapse.
- The obvious fix is wrong in both directions. Quad everywhere costs ~100x;
  escalating from double needs a threshold that fires on 94% of points. Moving
  the working type to `long double` improves the bias by 7000x for 2.5x.
- A cancellation monitor computed from already-corrupted terms is
  over-optimistic exactly where it matters, and needs a direct cross-check
  against a narrower type rather than a threshold on itself.

Before that is publishable, three things are missing: the physics derivation
check (the first open checkbox on O2-6340 — since done, see
`07-derivation.md`), the resolution of the discarded
weight above -- which changes what the precision work is *for* -- and ideally
(d), so the paper describes a generator rather than a bug fix. Computer Physics
Communications takes "new version of an existing program" papers when the code
is released; a CHEP proceeding is the lower-effort route.

## Who else uses this code

Searched 2026-08-07. Essentially only ALICE, plus inherited copies:

| repository | what |
|---|---|
| `AliceO2Group/AliceO2`, `O2DPG` | current production use |
| `alisw/AliRoot` | the legacy `AliGenEpEmv1`/`AliGenQEDBg` ancestors |
| `brettviren/ORKA-ILCRoot` | ILCRoot is an AliRoot-derived framework; carries `epemgen.f` and `diffcross.f` verbatim |
| personal forks (`jissong/AliRoot`, `yuanshiming/O2`, ...) | copies, not independent users |

So the bug has been inherited by AliRoot forks and by at least one non-ALICE
framework, but there is no evidence of an active third-party user. That
narrows the audience for an external paper and strengthens the case for an
internal note first.
