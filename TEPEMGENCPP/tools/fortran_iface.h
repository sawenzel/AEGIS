// Direct access to the unmodified Fortran in ../../TEPEMGEN, for use as the
// numerical oracle. Nothing here modifies the Fortran; the COMMON blocks are
// simply declared with their real layout so a C++ harness can read the
// diagnostics the Fortran keeps to itself (notably `badcount`, which is
// otherwise invisible from outside).
//
// PRECISION. Compile this translation unit with -DTEP_QUAD and link against
// diffcross.f compiled with `gfortran -freal-8-real-16` to get the SAME source
// evaluated in REAL*16. That gives a quad-precision reference without porting
// anything, which is what makes the Phase 0 precision study possible before
// Phase 1 exists. The two variants cannot live in one executable -- identical
// symbol names, different ABI -- so they are built as separate binaries.
#pragma once

#ifdef TEP_QUAD
using tepreal = __float128;
#else
using tepreal = double;
#endif

extern "C" {

// SUBROUTINE Initdiffcross(energy, mass)   energy in GeV, mass in MeV
void initdiffcross_(tepreal* energy, tepreal* mass);

// SUBROUTINE Diffcross(ppvt, yp, pmvt, ym, dphi, dsigma)
//   ppvt/pmvt  positron/electron transverse momentum [MeV]
//   yp/ym      positron/electron rapidity
//   dphi       azimuthal angle between them [rad]
//   dsigma     out: kbarn/MeV^4/(Z*alpha)^4
void diffcross_(tepreal* ppvt, tepreal* yp, tepreal* pmvt, tepreal* ym,
                tepreal* dphi, tepreal* dsigma);

// COMMON /badpar/ badcount
// Incremented whenever Diffcross discards a point (setzero, or dsigma < 0).
// Counts only the NEGATIVE half of the corrupted points -- positive garbage
// passes the guard untouched, which is the whole problem. See doc/02-findings.
extern struct { int badcount; } badpar_;

// COMMON /SETZPARAM/ setzero        set by Iz0 on an unphysical branch
extern struct { int setzero; } setzparam_;

// COMMON /PHYSPARAM/ gamma,beta,m,w1xw1,w1xw2,w2xw2,wl,wy
extern struct {
  tepreal gamma, beta, m, w1xw1, w1xw2, w2xw2, wl, wy;
} physparam_;

}  // extern "C"
