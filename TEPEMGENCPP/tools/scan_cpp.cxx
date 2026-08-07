// Regression test for O2-6340: reproduce the attached figure from the C++ port.
//
// Scans the differential cross section against electron rapidity at fixed
// remaining kinematics -- the same thing the ticket's plot shows -- evaluating
// with each precision strategy, and checks that the chosen one tracks the quad
// reference where plain double goes to noise.
//
// Emits the curves for plotting AND a pass/fail verdict, so it works as a test
// rather than only as a picture. The figure in the ticket is a scan at
// dphi = pi, which is the ill-conditioned locus (see doc/02-findings.md
// section 2): axy carries a factor sin(dphi) and the Id0..Id3 denominators
// lose their stabilising term there.
//
// Columns: ym  double  longdouble  adaptive  quad  survival_double
#include "AdaptiveDiffCross.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace o2::aegis::tepemgen;

int main(int argc, char** argv)
{
  if (argc < 8) {
    std::fprintf(stderr,
                 "usage: %s E_GeV ptp_MeV ptm_MeV yp dphi ymin ymax [nsteps] [tol]\n",
                 argv[0]);
    return 2;
  }
  const long energy = std::atol(argv[1]);
  const double ptp = std::atof(argv[2]), ptm = std::atof(argv[3]);
  const double yp = std::atof(argv[4]), dphi = std::atof(argv[5]);
  const double ymin = std::atof(argv[6]), ymax = std::atof(argv[7]);
  const int n = (argc > 8) ? std::atoi(argv[8]) : 400;
  const double tol = (argc > 9) ? std::atof(argv[9]) : 1e-4;
  const long double thr = (argc > 10) ? std::atof(argv[10]) : 1e-14L;

  const auto pd = initDiffCross<double>(double(energy), 5109991.0 / 10000000.0);
  const auto pl = initDiffCross<long double>((long double)energy,
                                             5109991.0L / 10000000.0L);
  const auto pq = initDiffCross<__float128>(
      __float128(energy), __float128(5109991) / __float128(10000000));
  AdaptiveDiffCross<long double, __float128> adaptive(energy, 1, 5109991,
                                                      10000000, thr);

  std::printf("# ym double longdouble adaptive quad survival_double\n");
  std::printf("# E=%ld ptp=%g ptm=%g yp=%g dphi=%.15g\n", energy, ptp, ptm, yp,
              dphi);

  double worstDouble = 0, worstAdaptive = 0;
  int nOverDouble = 0, nOverAdaptive = 0, nCompared = 0;

  for (int i = 0; i <= n; ++i) {
    const double ym = ymin + (ymax - ymin) * i / n;
    const auto rd = diffCross<double>(pd, ptp, yp, ptm, ym, dphi);
    const auto rl = diffCross<long double>(pl, (long double)ptp, (long double)yp,
                                           (long double)ptm, (long double)ym,
                                           (long double)dphi);
    const double va = adaptive(ptp, yp, ptm, ym, dphi);
    const auto rq = diffCross<__float128>(pq, __float128(ptp), __float128(yp),
                                          __float128(ptm), __float128(ym),
                                          __float128(dphi));
    const double vq = static_cast<double>(rq.dsigma);

    std::printf("%.10g %.17g %.17g %.17g %.17g %.6e\n", ym, rd.dsigma,
                static_cast<double>(rl.dsigma), va, vq,
                static_cast<double>(rd.survival));

    if (vq > 0 && rq.ok) {
      ++nCompared;
      const double ed = std::fabs(rd.dsigma - vq) / vq;
      const double ea = std::fabs(va - vq) / vq;
      if (ed > worstDouble) worstDouble = ed;
      if (ea > worstAdaptive) worstAdaptive = ea;
      if (ed > tol) ++nOverDouble;
      if (ea > tol) ++nOverAdaptive;
    }
  }

  std::fprintf(stderr, "compared           : %d\n", nCompared);
  std::fprintf(stderr, "double   worst rel : %.3e   (%d over %.0e)\n",
               worstDouble, nOverDouble, tol);
  std::fprintf(stderr, "adaptive worst rel : %.3e   (%d over %.0e)\n",
               worstAdaptive, nOverAdaptive, tol);
  std::fprintf(stderr, "escalated          : %lld / %lld\n",
               (long long)adaptive.escalated(), (long long)adaptive.calls());

  // The point of the exercise: the adaptive curve must stay on the reference
  // across the whole scan. No claim is made about double -- it is expected to
  // fail here, and that failure is the bug being fixed.
  const bool pass = nOverAdaptive == 0;
  std::fprintf(stderr, "\n%s\n", pass ? "PASS" : "FAIL");
  return pass ? 0 : 1;
}
