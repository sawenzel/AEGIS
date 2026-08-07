// Precision policy for the differential cross section.
//
// The naive reading of O2-6340 is "use __float128". Measurement says that is
// both too much and too little:
//
//   * too much -- __float128 is software-emulated and costs 97x double here
//     (1080 us/call against 11). TGenEpEmv1::Init() needed 25689 trials at
//     production settings, so quad throughout turns a 0.3 s initialisation
//     into 28 s, and the hard-coded fMaxXSTest = 1e7 worst case into 3 hours.
//
//   * too little -- selective escalation from double does not help, because
//     the cancellation is not a rare corner. The monitor never exceeds 1e-4
//     and roughly half the integral comes from points below 1e-12, so holding
//     the error under 1e-3 from a double base needs a 1e-7 threshold, which
//     escalates 94% of points and recovers almost none of the speed.
//
// The measured answer is a middle tier. `long double` (x87, 64-bit mantissa,
// eps 1.08e-19) buys three decades of headroom for 2.5x the cost, which is
// enough for the bulk; __float128 then handles only the ~1% tail. Calibrated
// on 100k points sampled from the real generator (tools/calibrate):
//
//   working type   mean rel err at survival 1e-12    cost
//   double                    7.6e-4                 11 us
//   long double               7.6e-7                 27 us
//
// See doc/02-findings.md and doc/03-plan.md for the full tables.
#ifndef TEPEMGENCPP_ADAPTIVEDIFFCROSS_H
#define TEPEMGENCPP_ADAPTIVEDIFFCROSS_H

#include "DiffCross.h"

#include <cstdint>

namespace o2::aegis::tepemgen
{

/// Evaluates in `Work`, re-evaluating in `Wide` when the cancellation monitor
/// says the `Work` result has no significant digits left.
///
/// LIMITATION, measured and not designed around: the monitor is
/// self-referential. Each type computes survival from its own terms, so a type
/// whose terms are already corrupted reports a survival that is too
/// optimistic. At ptp=ptm=1 MeV, yp=0, ym=0.05, dphi=pi the three tiers report
/// survival 3.5e-11 (double), 1.9e-14 (long double) and 2.4e-15 (quad, the
/// truth), while their values are 8.0e7, 4.4e4 and 5.5e3. The monitor also
/// cannot see cancellation *inside* Iz/Id/Iv -- notably the Id0..Id3
/// denominators B = 4u*axy^2 + A, which collapse as dphi -> pi because axy
/// carries a factor sin(dphi). That is exactly the locus the sampler favours.
///
/// The default threshold of 1e-13 is therefore chosen empirically rather than
/// from an error model: it is the largest value that keeps a deliberately
/// worst-case scan (dphi = pi, 400 points across the full rapidity range)
/// within 1e-4 of the quad reference. At 1e-14 that scan still leaves one
/// point 8x wrong. Do not raise it on the strength of the eps/survival
/// heuristic alone.
///
/// This replaces the Fortran's `IF(dsigma.LT.0) dsigma=0; badcount=badcount+1`.
/// That test is unsound: cancellation produces garbage of either sign, so it
/// catches only the negative half -- measured, it fires on 0.019% of points
/// while 0.25% are more than 10% wrong. Worse, the half it misses flows into
/// the event weight and, squared, into the variance estimate whose ratio to
/// the mean drives the convergence test in TGenEpEmv1::CalcXSection.
template <typename Work = long double, typename Wide = __float128>
class AdaptiveDiffCross
{
 public:
  /// `energy` per nucleon pair in GeV, `mass` the lepton mass in MeV.
  ///
  /// Both are taken as exact rationals rather than decimal literals: a literal
  /// is parsed at double and rounded before the wider type sees it. See the
  /// note in DiffCross.h.
  AdaptiveDiffCross(long numEnergy, long denEnergy, long numMass, long denMass,
                    Work threshold = Work(1e-13))
    : mWork(initDiffCross<Work>(Work(numEnergy) / Work(denEnergy),
                                Work(numMass) / Work(denMass))),
      mWide(initDiffCross<Wide>(Wide(numEnergy) / Wide(denEnergy),
                                Wide(numMass) / Wide(denMass))),
      mThreshold(threshold)
  {
  }

  /// Cross section at one kinematic point, in kbarn/MeV^4/(Z*alpha)^4.
  ///
  /// Unlike the Fortran this never silently returns 0 for a point it could not
  /// evaluate: a point that is genuinely unphysical (`ok == false`, the branch
  /// the Fortran flagged via COMMON /SETZPARAM/) returns 0 and is counted, and
  /// everything else is computed to a precision the monitor vouches for.
  double operator()(double ppvt, double yp, double pmvt, double ym, double dphi)
  {
    ++mCalls;
    const auto r = diffCross<Work>(mWork, static_cast<Work>(ppvt),
                                   static_cast<Work>(yp),
                                   static_cast<Work>(pmvt),
                                   static_cast<Work>(ym),
                                   static_cast<Work>(dphi));
    if (!r.ok) {
      ++mUnphysical;
      return 0.0;
    }
    if (r.survival >= mThreshold) {
      return static_cast<double>(r.dsigma);
    }
    ++mEscalated;
    const auto w = diffCross<Wide>(mWide, static_cast<Wide>(ppvt),
                                   static_cast<Wide>(yp),
                                   static_cast<Wide>(pmvt),
                                   static_cast<Wide>(ym),
                                   static_cast<Wide>(dphi));
    if (!w.ok) {
      ++mUnphysical;
      return 0.0;
    }
    return static_cast<double>(w.dsigma);
  }

  std::int64_t calls() const { return mCalls; }
  std::int64_t escalated() const { return mEscalated; }
  std::int64_t unphysical() const { return mUnphysical; }
  Work threshold() const { return mThreshold; }
  void setThreshold(Work t) { mThreshold = t; }

 private:
  PhysParams<Work> mWork;
  PhysParams<Wide> mWide;
  Work mThreshold;
  std::int64_t mCalls = 0;
  std::int64_t mEscalated = 0;
  std::int64_t mUnphysical = 0;
};

}  // namespace o2::aegis::tepemgen

#endif  // TEPEMGENCPP_ADAPTIVEDIFFCROSS_H
