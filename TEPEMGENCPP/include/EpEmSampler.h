// The sampling layer of ../TEPEMGEN/epemgen.f: ee_init and ee_event.
//
// Structure of the original, preserved here because it is what makes the
// weights meaningful: events are drawn from cheap analytic ENVELOPES
// (DsdYpY x DsdYmY x DsdXpX x DsdXmX x DsdPhi, see Envelopes.h), and each is
// then reweighted by the EXACT differential cross section. The envelope only
// has to be a reasonable shape; the exact cross section carries the physics.
//
// Two things are deliberately faithful to the Fortran:
//
//   * The order in which random numbers are consumed. The acceptance gate for
//     this port is statistical, but preserving the stream order means that
//     while the Fortran is still in the tree a divergence can be localised to
//     a single event instead of inferred from a histogram. It costs nothing.
//
//   * The rejection loops, including the `go to` targets, are transcribed as
//     loops with the same accept/reject conditions rather than restructured.
//     Changing which variables are resampled on rejection changes the
//     distribution, subtly and silently.
//
// What is NOT faithful, on purpose:
//
//   * COMMON /eepars/ and /eevent/ become members, so two samplers can coexist.
//   * Dgauss and Dtrint are gone -- see Quadrature.h and Envelopes.h.
//   * ee_init printed an error and CARRIED ON when the normalisation failed,
//     leaving XYsect = 0 and every weight zero. init() returns false instead.
#ifndef TEPEMGENCPP_EPEMSAMPLER_H
#define TEPEMGENCPP_EPEMSAMPLER_H

#include "AdaptiveDiffCross.h"
#include "Envelopes.h"
#include "Quadrature.h"

#include <cmath>
#include <cstdint>
#include <functional>
#include <string>

namespace o2::aegis::tepemgen
{

/// DsdPhi parameters, from the DATA statement in epemgen.f.
/// Single precision, like the other DATA statements -- see f32() in Envelopes.h.
inline constexpr double kParPhi[7] = {f32(5.4825), f32(226.00), f32(11.602),
                                      f32(68.173), f32(1.9509), f32(0.11864e+04),
                                      f32(109.61)};

inline double dsdPhi(double phi)
{
  const double a = std::fabs(phi - 3.14159265358979323846);
  return kParPhi[0] + kParPhi[1] * std::exp(-kParPhi[2] * a) +
         kParPhi[3] * std::exp(-kParPhi[4] * a) +
         kParPhi[5] * std::exp(-kParPhi[6] * a);
}

class EpEmSampler
{
 public:
  struct Config {
    double yMin = -7.0, yMax = 7.0;      ///< rapidity range
    double ptMinMeV = 1.0, ptMaxMeV = 1000.0;
    double cmEnergyGeV = 5360.0;         ///< per nucleon pair
    double z = 82.0;                     ///< atomic number, symmetric systems
  };

  struct Event {
    double yElectron, yPositron;   ///< rapidities
    double xElectron, xPositron;   ///< log10(pt/MeV)
    double phi;                    ///< azimuth between them, radians
    double weight;                 ///< Wtm2
  };

  /// `rng` must return a uniform deviate in the OPEN interval (0,1) -- the
  /// Fortran's eernd rejects 0 and 1 explicitly, and log(0) appears below.
  EpEmSampler(const Config& cfg, std::function<double()> rng)
    : mCfg(cfg), mRng(std::move(rng)),
      // Exact rationals: a decimal literal would be rounded at double before
      // the wider working type sees it. See DiffCross.h.
      mCross(static_cast<long>(cfg.cmEnergyGeV), 1, 5109991, 10000000)
  {
  }

  /// Returns false instead of continuing with a zero normalisation.
  bool init()
  {
    mError.clear();
    mXmin = std::log10(mCfg.ptMinMeV);
    mXmax = std::log10(mCfg.ptMaxMeV);
    if (!(mCfg.yMin < mCfg.yMax)) { mError = "yMin >= yMax"; return false; }
    if (!(mXmin < mXmax)) { mError = "ptMin >= ptMax"; return false; }

    // --- normalisation -----------------------------------------------------
    // The Fortran integrates the envelopes over two triangles with Dtrint.
    // Envelopes.h reduces that to one 1-D integral; see doc/02-findings.md
    // section 11.
    const double xsecX = squareIntegral(
        [](double u) { return dsdXpX(u); }, [](double L) { return intXmX(L); },
        mXmin, mXmax, kXpXKink);
    const double xsecY = squareIntegral(
        [](double u) { return dsdYpY(u); }, [](double L) { return intYmY(L); },
        mCfg.yMin, mCfg.yMax, 0.0);
    if (!(xsecX > 0) || !(xsecY > 0)) {
      mError = "envelope normalisation integral is not positive";
      return false;
    }
    // 2371.5239*(Z/137.035)^4 -- the Fortran's normalisation factor.
    const double za = mCfg.z / f32(137.035);
    mXYsect = f32(2371.5239) * za * za * za * za * xsecX * xsecY;

    // --- rapidity ----------------------------------------------------------
    mYpYmin = 2 * mCfg.yMin;
    mYpYmax = 2 * mCfg.yMax;
    double yp0 = 0.0;
    if (mYpYmin * mYpYmax > 0) yp0 = (mYpYmin > 0) ? mYpYmin : mYpYmax;
    mWYpYmax = dsdYpY(yp0);
    mYmYmax = mCfg.yMax - mCfg.yMin;
    mSgmY1 = 1.0 / std::sqrt(2 * kParYmY[1]);
    mSgmY2 = 1.0 / std::sqrt(2 * kParYmY[3]);

    // The Fortran's three Dgauss calls. Their limits already sit ON the kinks
    // of DsdYmY (0.18 and 4.00), so unlike Dtrint they never had to converge
    // across one -- that is why only Dtrint failed.
    auto ymy = [](double y) { return dsdYmY(y); };
    const double g1 = integrate<double>(ymy, 0.0, kYmYKink1, {}, 1e-12).value;
    const double g2 = integrate<double>(ymy, kYmYKink1, kYmYKink2, {}, 1e-12).value;
    const double ge = integrate<double>(ymy, kYmYKink2, mYmYmax, {}, 1e-12).value;
    const double sy = g1 + g2 + ge;
    if (!(sy > 0)) { mError = "DsdYmY normalisation is not positive"; return false; }
    mGaus1 = g1 / sy;
    mGaus2 = (g1 + g2) / sy;

    // --- azimuth -----------------------------------------------------------
    const double pi = 3.14159265358979323846;
    const double e0 = kParPhi[0] * pi;
    const double e1 = kParPhi[1] * (1 - std::exp(-kParPhi[2] * pi)) / kParPhi[2];
    const double e2 = kParPhi[3] * (1 - std::exp(-kParPhi[4] * pi)) / kParPhi[4];
    const double e3 = kParPhi[5] * (1 - std::exp(-kParPhi[6] * pi)) / kParPhi[6];
    const double sp = e0 + e1 + e2 + e3;
    mExp1 = (e1 + e2 + e3) / sp;
    mExp2 = (e2 + e3) / sp;
    mExp3 = e3 / sp;

    // --- transverse momentum ----------------------------------------------
    mXmXmax = mXmax - mXmin;
    mXpXmin = 2 * mXmin;
    mXpXmax = 2 * mXmax;
    const double x1 = kParXmX[0] * (1 - std::exp(-kParXmX[1] * std::fabs(mXmXmax))) / kParXmX[1];
    const double x2 = kParXmX[2] * (1 - std::exp(-kParXmX[3] * std::fabs(mXmXmax))) / kParXmX[3];
    mExmx1 = x1 / (x1 + x2);

    mIcase = 2;
    mXmed = kXpXKink;
    if (mXpXmax < mXmed) { mIcase = 1; mXmed = mXpXmax; }
    if (mXpXmin > mXmed) { mIcase = 3; mXmed = mXpXmin; }
    if (mIcase == 2) {
      auto xpx = [](double x) { return dsdXpX(x); };
      const double gs = integrate<double>(xpx, mXpXmin, mXmed, {}, 1e-12).value;
      const double ex = integrate<double>(xpx, mXmed, mXpXmax, {}, 1e-12).value;
      const double s = gs + ex;
      if (!(s > 0)) { mError = "DsdXpX normalisation is not positive"; return false; }
      mGaussFrac = gs / s;
    }
    if (mIcase == 1 || mIcase == 2) mSgm1 = 1.0 / std::sqrt(2 * kParXpX[1]);
    if (mIcase == 3 || mIcase == 2)
      mAnorX = 1 - std::exp(-kParXpX[4] * (mXpXmax - mXmed));

    mNevent = 0; mSum = 0; mSum2 = 0;
    return true;
  }

  /// One event. Mirrors ee_event, including which variables are resampled on
  /// each rejection.
  Event next()
  {
    const double pi = 3.14159265358979323846;
    double ypy = 0, ymy = 0, ye = 0, yp = 0;

    // Rapidities. Note the Fortran's `go to 10` from the acceptance test
    // resamples BOTH YpY and YmY, not YmY alone.
    for (;;) {
      do {
        ypy = mYpYmin + (mYpYmax - mYpYmin) * mRng();
      } while (dsdYpY(ypy) < mRng() * mWYpYmax);

      const double r1 = mRng();
      if (r1 < mGaus1) {
        do { ymy = mSgmY1 * normal(); } while (std::fabs(ymy) > kYmYKink1);
      } else if (r1 < mGaus2) {
        do {
          ymy = mSgmY2 * normal();
        } while (std::fabs(ymy) < kYmYKink1 || std::fabs(ymy) > kYmYKink2);
      } else {
        do { ymy = -std::log(mRng()) / kParYmY[5] + kYmYKink2; }
        while (ymy > mYmYmax);
        if (mRng() < 0.5) ymy = -ymy;
      }

      ye = 0.5 * (ypy + ymy);
      yp = 0.5 * (ypy - ymy);
      if (ye >= mCfg.yMin && ye <= mCfg.yMax && yp >= mCfg.yMin && yp <= mCfg.yMax)
        break;
    }

    // Azimuth. The cut removes the COLLINEAR configuration (phi near 0 or
    // 2pi); back-to-back, which is both the peak of DsdPhi and the
    // worst-conditioned region of the cross section, is left untouched.
    double phi = 0;
    for (;;) {
      const double r1 = mRng();
      if (r1 < mExp3) {
        do { phi = -std::log(mRng()) / kParPhi[6]; } while (phi > pi);
      } else if (r1 < mExp2) {
        do { phi = -std::log(mRng()) / kParPhi[4]; } while (phi > pi);
      } else if (r1 < mExp1) {
        do { phi = -std::log(mRng()) / kParPhi[2]; } while (phi > pi);
      } else {
        phi = pi * mRng();
      }
      if (mRng() > 0.5) phi = -phi;
      phi += pi;
      if (phi >= f32(0.03) && phi <= 2 * pi - f32(0.03)) break;
    }

    // Transverse momenta.
    double xe = 0, xp = 0;
    for (;;) {
      int jcase = mIcase;
      if (mIcase == 2) jcase = (mRng() < mGaussFrac) ? 1 : 3;

      double xpx = 0;
      if (jcase == 1) {
        do { xpx = mSgm1 * normal() - kParXpX[2]; }
        while (xpx < mXpXmin || xpx > mXmed);
      } else {
        xpx = -std::log(1.0 - mRng() * mAnorX) / kParXpX[4] + mXmed;
      }

      double xmx = 0;
      if (mRng() < mExmx1) {
        do { xmx = -std::log(mRng()) / kParXmX[1]; }
        while (xmx > std::fabs(mXmXmax));
      } else {
        do { xmx = -std::log(mRng()) / kParXmX[3]; }
        while (xmx > std::fabs(mXmXmax));
      }
      if (mRng() > 0.5) xmx = -xmx;

      xe = 0.5 * (xpx + xmx);
      xp = 0.5 * (xpx - xmx);
      if (xe >= mXmin && xe <= mXmax && xp >= mXmin && xp <= mXmax) {
        mLastXpX = xpx;
        mLastXmX = xmx;
        break;
      }
    }

    const double pte = std::pow(10.0, xe);
    const double ptp = std::pow(10.0, xp);

    // Envelope weight, then the exact cross section. 2.2706950 is the
    // Fortran's magic envelope normalisation.
    double wt = dsdYpY(ypy) * dsdYmY(ymy) * dsdXpX(mLastXpX) *
                dsdXmX(mLastXmX) * dsdPhi(phi);
    wt = f32(2.2706950) / wt;
    const double dsigma = mCross(ptp, yp, pte, ye, phi);
    const double wtm2 = mXYsect * dsigma * (pte * ptp) * (pte * ptp) * wt;

    ++mNevent;
    mSum += wtm2;
    mSum2 += static_cast<long double>(wtm2) * wtm2;
    return {ye, yp, xe, xp, phi, wtm2};
  }

  /// Cross section accumulated so far, in kbarn (as the Fortran's Xsecttot).
  double xSection() const { return mNevent ? static_cast<double>(mSum / mNevent) : 0.0; }
  double xSectionError() const
  {
    return mNevent ? static_cast<double>(std::sqrt((double)mSum2) / mNevent) : 0.0;
  }
  std::int64_t events() const { return mNevent; }
  const std::string& error() const { return mError; }
  const AdaptiveDiffCross<>& cross() const { return mCross; }
  AdaptiveDiffCross<>& cross() { return mCross; }

 private:
  /// Box-Muller, exactly as SUBROUTINE rnormlEP: two deviates in, cos branch
  /// out, second value discarded. Caching it would change the stream.
  double normal()
  {
    const double u1 = mRng();
    const double u2 = mRng();
    return std::cos(2 * 3.14159265358979323846 * u1) * std::sqrt(-2 * std::log(u2));
  }

  /// integral over [lo,hi]^2 of f(a+b)*g(a-b), via the diamond reduction.
  template <typename F, typename G>
  static double squareIntegral(F f, G innerClosedForm, double lo, double hi,
                               double fKink)
  {
    const double u0 = 2 * lo, u1 = 2 * hi;
    auto integrand = [&](double u) {
      return f(u) * innerClosedForm(std::min(u - u0, u1 - u));
    };
    // The diamond half-width kinks at u = lo+hi; f may kink too.
    const auto r = (fKink > u0 && fKink < u1)
                       ? integrate<double>(integrand, u0, u1, {fKink, lo + hi}, 1e-12)
                       : integrate<double>(integrand, u0, u1, {lo + hi}, 1e-12);
    return 0.5 * r.value;
  }

  Config mCfg;
  std::function<double()> mRng;
  AdaptiveDiffCross<> mCross;
  std::string mError;

  double mXmin = 0, mXmax = 0, mXYsect = 0;
  double mYpYmin = 0, mYpYmax = 0, mWYpYmax = 0, mYmYmax = 0;
  double mSgmY1 = 0, mSgmY2 = 0, mGaus1 = 0, mGaus2 = 0;
  double mExp1 = 0, mExp2 = 0, mExp3 = 0;
  double mXmXmax = 0, mXpXmin = 0, mXpXmax = 0, mExmx1 = 0;
  double mXmed = 0, mSgm1 = 0, mAnorX = 0, mGaussFrac = 0;
  int mIcase = 2;
  double mLastXpX = 0, mLastXmX = 0;

  std::int64_t mNevent = 0;
  long double mSum = 0, mSum2 = 0;
};

}  // namespace o2::aegis::tepemgen

#endif  // TEPEMGENCPP_EPEMSAMPLER_H
