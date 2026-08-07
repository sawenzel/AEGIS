// The catastrophic points, found by scanning 12 million sampled events
// (4 seeds x 3M) at production settings and keeping every point whose DOUBLE
// cross section exceeds 10000 -- more than twice the largest LEGITIMATE value
// anywhere in those samples (~4400).
//
// These are the events behind O2-6340. Each is positive, so the Fortran guard
// `IF((setzero.EQ.1).OR.(dsigma.LT.0))` never sees it. Squared into Dsect2 it
// inflates the variance estimate by hundreds of percent, the convergence test
// `err/xSect < eps` in TGenEpEmv1::CalcXSection can then never be satisfied,
// the loop runs to fMaxXSTest = 1e7 and calls abort(). Measured over 3M-event
// samples: variance bias +841% (seed 90210) and +598% (seed 5150), against
// -0.03% for a seed that happened not to draw one. That is the reported
// "sometimes".
//
// They share a signature: both leptons at MINIMUM pt (~1 MeV, the ptMin cut),
// both near the rapidity EDGE (|y| 6.3-6.9, the |y|<7 cut), and back-to-back
// in azimuth. A corner of the sampled phase space, not a generic point.
#ifndef TEPEMGENCPP_WORST_POINTS_H
#define TEPEMGENCPP_WORST_POINTS_H

namespace o2::aegis::tepemgen::testdata
{

/// Only the kinematics are stored. The expected values are computed at
/// runtime from the quad instantiation rather than hard-coded: the
/// coordinates here are rounded to 10 digits, so a hard-coded "truth" would
/// be the truth of a slightly different point and would rot silently.
///
/// That rounding is itself informative -- the catastrophe survives it, so this
/// is a small ill-conditioned REGION, not a knife-edge point. At the rounded
/// seed-777 coordinates the double path returns 3.2e6 against a true 0.036,
/// which is worse than at the exact ones.
struct WorstPoint {
  const char* seed;
  double ppvt, yp, pmvt, ym, dphi;
};

/// energy for all of these
inline constexpr double kWorstPointsEnergyGeV = 5360.0;

inline constexpr WorstPoint kWorstPoints[] = {
    {"777", 1.074808674, 6.735805349, 1.124250238, 6.808340873, 3.141516602503257},
    {"90210", 1.030897896, 6.888477533, 1.087895309, 6.676640229, 3.141652469577499},
    {"5150a", 1.003855027, 6.790163181, 1.054785945, 6.8270169, 3.142209904711982},
    {"5150b", 1.167106695, 6.338589285, 1.114665541, 6.88131587, 3.141841463999655},
};

}  // namespace o2::aegis::tepemgen::testdata

#endif  // TEPEMGENCPP_WORST_POINTS_H
