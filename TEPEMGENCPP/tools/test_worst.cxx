// Regression test: the points that actually break the generator.
//
// Every one of these was found in real sampled events and is what O2-6340
// reports -- a positive outlier the Fortran's sign guard cannot see, large
// enough to poison the variance estimate and hang CalcXSection until it
// aborts. The test asserts that the adaptive path returns the quad answer at
// each of them, and reports which trigger caught it.
//
// It also documents WHY there are two triggers. At the seed-90210 point the
// cancellation monitor alone is not enough: the long double working type
// reports survival 2.11e-13 -- above a 1e-13 threshold, so no escalation --
// while the true survival is 7.3e-15 and its answer is -1.33 against a true
// 0.046. The monitor is computed from terms that are themselves already
// corrupted, so it is over-optimistic exactly where it matters most.
#include "AdaptiveDiffCross.h"
#include "worst_points.h"

#include <cmath>
#include <cstdio>

using namespace o2::aegis::tepemgen;
using namespace o2::aegis::tepemgen::testdata;

int main()
{
  const double E = kWorstPointsEnergyGeV;
  const auto pd = initDiffCross<double>(E, 5109991.0 / 10000000.0);
  const auto pl = initDiffCross<long double>((long double)E, 5109991.0L / 10000000.0L);
  const auto pq = initDiffCross<__float128>(
      __float128(E), __float128(5109991) / __float128(10000000));

  std::printf("%-7s %13s %13s %13s %13s %10s\n", "seed", "double",
              "long double", "adaptive", "quad(truth)", "verdict");

  int failures = 0;
  for (const auto& w : kWorstPoints) {
    const auto rd = diffCross<double>(pd, w.ppvt, w.yp, w.pmvt, w.ym, w.dphi);
    const auto rl = diffCross<long double>(
        pl, (long double)w.ppvt, (long double)w.yp, (long double)w.pmvt,
        (long double)w.ym, (long double)w.dphi);
    const auto rq = diffCross<__float128>(
        pq, __float128(w.ppvt), __float128(w.yp), __float128(w.pmvt),
        __float128(w.ym), __float128(w.dphi));
    const double truth = static_cast<double>(rq.dsigma);

    AdaptiveDiffCross<> ad(E);
    const double va = ad(w.ppvt, w.yp, w.pmvt, w.ym, w.dphi);

    const double rel = std::fabs(va - truth) / std::fabs(truth);
    const bool ok = rel < 1e-6;
    if (!ok) ++failures;
    std::printf("%-7s %13.6g %13.6g %13.6g %13.6g %10s\n", w.seed, rd.dsigma,
                (double)rl.dsigma, va, truth, ok ? "ok" : "FAIL");
    std::printf("        double is %.3gx the truth; survival: double %.2e, "
                "long double %.2e, quad %.2e\n",
                std::fabs(rd.dsigma / truth), (double)rd.survival,
                (double)rl.survival, (double)rq.survival);
    std::printf("        caught by: %s%s   (rel err %.2e)\n",
                ad.escalated() ? "escalation" : "NOTHING",
                ad.disagreed() ? " via narrow-type disagreement"
                               : (ad.escalated() ? " via the monitor" : ""),
                rel);
  }
  std::printf("\n%s\n", failures == 0 ? "PASS" : "FAIL");
  return failures == 0 ? 0 : 1;
}
