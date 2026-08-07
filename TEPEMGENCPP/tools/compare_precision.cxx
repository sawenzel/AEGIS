// Read points produced by dump_points (which carry the double cross section
// computed by the double build), re-evaluate each in REAL*16, and report what
// the loss of significance does to the quantities the generator actually uses.
//
// Quad-only: built with -DTEP_QUAD against diffcross.f compiled
// -freal-8-real-16. Streams, so it runs over tens of millions of points
// without producing an intermediate file.
//
// The headline numbers are the biases on sum(dsigma) and sum(dsigma^2) --
// the cross section and the variance estimate that TGenEpEmv1::CalcXSection
// accumulates and whose ratio drives its convergence test.
#include "fortran_iface.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

struct Outlier {
  double ppvt, yp, pmvt, ym, dphi, vd, vq, ratio;
};

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s E_GeV pointsfile [ntop]\n", argv[0]);
    return 2;
  }
  const int ntop = (argc > 3) ? std::atoi(argv[3]) : 10;

  tepreal energy = static_cast<tepreal>(std::atof(argv[1]));
  tepreal mass = static_cast<tepreal>(0.5109991);
  initdiffcross_(&energy, &mass);

  FILE* f = std::fopen(argv[2], "r");
  if (!f) { std::perror("fopen"); return 1; }

  long n = 0, guard = 0, e1 = 0, e3 = 0, e6 = 0, e9 = 0;
  long double sd = 0, sq = 0, sd2 = 0, sq2 = 0, sabs = 0;
  double maxq = 0;
  std::vector<Outlier> top;   // ranked by |vd - vq|

  char line[1024];
  while (std::fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    double a, b, c, d, e, w, vd;
    int bad;
    if (std::sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %d", &a, &b, &c, &d, &e,
                    &w, &vd, &bad) < 8)
      continue;
    tepreal ta = a, tb = b, tc = c, td = d, te = e, out = 0;
    diffcross_(&ta, &tb, &tc, &td, &te, &out);
    const double vq = static_cast<double>(out);

    ++n;
    guard += bad;
    sd += vd; sq += vq;
    sd2 += static_cast<long double>(vd) * vd;
    sq2 += static_cast<long double>(vq) * vq;
    maxq = std::max(maxq, vq);

    const double adiff = std::fabs(vd - vq);
    sabs += adiff;
    if (vq > 0) {
      const double r = adiff / vq;
      if (r > 1e-9) ++e9;
      if (r > 1e-6) ++e6;
      if (r > 1e-3) ++e3;
      if (r > 1e-1) ++e1;
    }

    // Rank by ABSOLUTE error: a huge relative error next to a zero of the
    // cross section is arithmetically real but physically irrelevant, and
    // ranking by it just surfaces the pair-at-rest configuration every time.
    if (static_cast<int>(top.size()) < ntop || adiff > top.back().ratio) {
      top.push_back({a, b, c, d, e, vd, vq, adiff});
      std::sort(top.begin(), top.end(),
                [](const Outlier& x, const Outlier& y) { return x.ratio > y.ratio; });
      if (static_cast<int>(top.size()) > ntop) top.pop_back();
    }
  }
  std::fclose(f);
  if (!n) { std::fprintf(stderr, "no points read\n"); return 1; }

  std::printf("points                       : %ld\n", n);
  std::printf("guard fired (dsigma -> 0)    : %ld  (%.4f%%)\n", guard,
              100.0 * guard / n);
  std::printf("rel.err > 1e-9               : %.3f%%\n", 100.0 * e9 / n);
  std::printf("rel.err > 1e-6               : %.3f%%\n", 100.0 * e6 / n);
  std::printf("rel.err > 1e-3               : %.3f%%\n", 100.0 * e3 / n);
  std::printf("rel.err > 1e-1               : %.4f%%\n", 100.0 * e1 / n);
  std::printf("\n");
  std::printf("sum(dsigma)   double = %.12Lg\n", sd);
  std::printf("sum(dsigma)   quad   = %.12Lg\n", sq);
  std::printf("  cross-section bias = %+.5f %%\n",
              static_cast<double>(100.0L * (sd - sq) / sq));
  std::printf("sum(dsigma^2) double = %.10Lg\n", sd2);
  std::printf("sum(dsigma^2) quad   = %.10Lg\n", sq2);
  std::printf("  variance bias      = %+.5f %%\n",
              static_cast<double>(100.0L * (sd2 - sq2) / sq2));
  std::printf("sum|err| / integral  = %.5f %%\n",
              static_cast<double>(100.0L * sabs / sq));
  std::printf("max quad dsigma      = %.6g\n", maxq);
  std::printf("\ntop %d by ABSOLUTE error (share of max value in sample):\n",
              ntop);
  std::printf("%12s %8s %8s %10s %14s %14s %10s\n", "abs_err", "yp", "ym",
              "dphi-pi", "double", "quad", "err/maxq");
  for (const auto& o : top)
    std::printf("%12.6g %8.3f %8.3f %10.2e %14.6g %14.6g %10.3f\n", o.ratio,
                o.yp, o.ym, o.dphi - 3.14159265358979, o.vd, o.vq,
                o.ratio / maxq);
  return 0;
}
