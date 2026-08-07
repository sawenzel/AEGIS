// Physics sanity check: integrate the differential cross section over the FULL
// phase space and compare with the literature value for coherent e+e- pair
// production in ultraperipheral PbPb at the LHC, which is of order 200 kb.
//
// This is NOT a derivation check. It cannot tell you whether the 5-fold
// differential cross section in diffcross.f correctly encodes the amplitudes
// of Alscher, Hencken, Trautmann and Baur, Phys. Rev. A55 (1997) 396. What it
// CAN do is catch the failure modes that a port or a 25-year-old code is
// actually likely to have: a wrong overall normalisation, a wrong power of
// Z*alpha, a units slip between MeV and GeV or barn and kbarn. Those would
// show up as orders of magnitude, not percent.
//
// The recipe is the one stated in the header of diffcross.f:
//
//   "to get the total cross section, you have to integrate over
//    Integral dsigma dyp dym ddphi 2 pi ppt dppt pmt dpmt"
//
// with dsigma in kbarn/MeV^4/(Z*alpha)^4. So
//
//   sigma = (Z*alpha)^4 * 2pi * Int dsigma * p+ dp+ * p- dp- * dy+ dy- ddphi
//
// Transverse momenta are sampled log-uniformly, since the integrand spans many
// decades and is concentrated near pt ~ m_e; Int p dp f = Int p^2 f d(ln p).
#include "AdaptiveDiffCross.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>

using namespace o2::aegis::tepemgen;

int main(int argc, char** argv)
{
  const long n = (argc > 1) ? std::atol(argv[1]) : 200000;
  const double energy = (argc > 2) ? std::atof(argv[2]) : 5360.0;
  const double ptLo = (argc > 3) ? std::atof(argv[3]) : 0.02;    // MeV
  const double ptHi = (argc > 4) ? std::atof(argv[4]) : 5000.0;  // MeV
  const double yMax = (argc > 5) ? std::atof(argv[5]) : 12.0;
  const int seed = (argc > 6) ? std::atoi(argv[6]) : 20260807;

  const double Z = 82.0;
  const double za = Z / 137.035;
  const double za4 = za * za * za * za;
  const double pi = 3.14159265358979323846;

  AdaptiveDiffCross<> cross(energy);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 1.0);

  const double lLo = std::log(ptLo), lHi = std::log(ptHi);
  const double volume = (lHi - lLo) * (lHi - lLo)   // d(ln p+) d(ln p-)
                        * (2 * yMax) * (2 * yMax)   // dy+ dy-
                        * (2 * pi);                 // ddphi

  long double sum = 0, sum2 = 0;
  // A log-uniform Monte Carlo over an integrand spanning ~20 decades is
  // heavy-tailed, and for heavy tails the naive standard error is itself
  // unreliable AND biased low -- the rare large contributions are the ones
  // most often missed. Track the single largest contribution as a fraction of
  // the total: if one sample carries a noticeable share, the quoted error is
  // meaningless and the estimate is a lower bound, not a measurement.
  double maxf = 0;
  for (long i = 0; i < n; ++i) {
    const double ptp = std::exp(lLo + (lHi - lLo) * u(rng));
    const double ptm = std::exp(lLo + (lHi - lLo) * u(rng));
    const double yp = -yMax + 2 * yMax * u(rng);
    const double ym = -yMax + 2 * yMax * u(rng);
    const double dphi = 2 * pi * u(rng);

    const double ds = cross(ptp, yp, ptm, ym, dphi);
    // p dp = p^2 d(ln p) for each leg
    const double f = ds * ptp * ptp * ptm * ptm;
    sum += f;
    sum2 += static_cast<long double>(f) * f;
    if (f > maxf) maxf = f;
  }

  const double mean = static_cast<double>(sum / n);
  const double var = static_cast<double>(sum2 / n) - mean * mean;
  const double err = std::sqrt(var > 0 ? var / n : 0.0);

  const double sigma_kb = za4 * 2 * pi * volume * mean;
  const double sigma_kb_err = za4 * 2 * pi * volume * err;

  std::printf("samples            : %ld\n", n);
  std::printf("sqrt(s_NN)         : %g GeV     Z = %g\n", energy, Z);
  std::printf("pt range           : %g .. %g MeV\n", ptLo, ptHi);
  std::printf("rapidity range     : +-%g  (beam rapidity is %.2f)\n", yMax,
              std::acosh(energy / (2 * 0.938)));
  std::printf("escalated to quad  : %.1f%%\n",
              100.0 * cross.escalated() / (cross.calls() ? cross.calls() : 1));
  std::printf("\nsigma(total)       : %.4g +- %.2g kb   (%.1f%% stat)\n",
              sigma_kb, sigma_kb_err,
              sigma_kb != 0 ? 100 * sigma_kb_err / sigma_kb : 0.0);
  const double topShare = (sum > 0) ? static_cast<double>(maxf / sum) : 0.0;
  std::printf("largest single sample: %.2f%% of the whole sum%s\n",
              100 * topShare,
              topShare > 0.01 ? "   <-- TAIL-DOMINATED, treat as a lower bound"
                              : "");
  std::printf("\nreference points\n");
  std::printf("  Racah (Born, point-like, L=ln gamma^2) : ~223 kb\n");
  std::printf("  literature, UPC PbPb at LHC            : ~200 kb\n");
  std::printf("  ratio to Racah                         : %.3f\n", sigma_kb / 223.0);
  return 0;
}
