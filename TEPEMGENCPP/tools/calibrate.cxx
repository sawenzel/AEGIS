// Calibrate the escalation policy: which working type, and when to escalate?
//
// Pure C++ -- no Fortran linked. The port is validated against the Fortran in
// Phase 1, so the __float128 instantiation is the reference here, and one
// binary can hold every instantiation (which the two Fortran builds could not,
// sharing symbol names with different ABI).
//
// For each candidate working type this bins its error against the quad
// reference by ITS OWN survival value -- the quantity actually available at
// runtime when the escalation decision has to be made. Binning by the quad
// survival would be cheating.
#include "DiffCross.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

using namespace o2::aegis::tepemgen;

namespace
{

struct Point { double a, b, c, d, e; };

constexpr int kBins = 20;

struct Stats {
  long n[kBins] = {};
  long double sum[kBins] = {};
  double mx[kBins] = {};
  long over3[kBins] = {}, over6[kBins] = {};
  long double contrib[kBins] = {};
  long double total = 0;
  long all = 0;
};

int binOf(double surv)
{
  if (!(surv > 0)) return kBins - 1;
  const double l = -std::log10(surv);
  return l < 0 ? 0 : (l >= kBins ? kBins - 1 : static_cast<int>(l));
}

template <typename Fast>
Stats measure(const std::vector<Point>& pts, double energy,
              const std::vector<double>& ref)
{
  const auto p = initDiffCross<Fast>(static_cast<Fast>(energy),
                                     Fast(5109991) / Fast(10000000));
  Stats s;
  for (size_t i = 0; i < pts.size(); ++i) {
    const double vq = ref[i];
    if (!(vq > 0)) continue;
    const auto r = diffCross<Fast>(p, static_cast<Fast>(pts[i].a),
                                   static_cast<Fast>(pts[i].b),
                                   static_cast<Fast>(pts[i].c),
                                   static_cast<Fast>(pts[i].d),
                                   static_cast<Fast>(pts[i].e));
    const double rel = std::fabs(static_cast<double>(r.dsigma) - vq) / vq;
    const int b = binOf(static_cast<double>(r.survival));
    ++s.n[b]; ++s.all;
    s.sum[b] += rel;
    if (rel > s.mx[b]) s.mx[b] = rel;
    if (rel > 1e-3) ++s.over3[b];
    if (rel > 1e-6) ++s.over6[b];
    s.contrib[b] += vq;
    s.total += vq;
  }
  return s;
}

void report(const char* name, double eps, const Stats& s)
{
  std::printf("\n=== working type: %s (eps = %.2e) ===\n", name, eps);
  std::printf("%16s %8s %11s %11s %8s %8s %10s\n", "survival", "n", "mean_rel",
              "max_rel", ">1e-6", ">1e-3", "%integral");
  for (int i = 0; i < kBins; ++i) {
    if (!s.n[i]) continue;
    char lab[32];
    if (i == kBins - 1) std::snprintf(lab, sizeof lab, "< 1e-%d", kBins - 1);
    else std::snprintf(lab, sizeof lab, "1e-%d..1e-%d", i + 1, i);
    std::printf("%16s %8ld %11.3Le %11.3e %7.2f%% %7.2f%% %9.3f%%\n", lab,
                s.n[i], s.sum[i] / s.n[i], s.mx[i], 100.0 * s.over6[i] / s.n[i],
                100.0 * s.over3[i] / s.n[i],
                static_cast<double>(100.0L * s.contrib[i] / s.total));
  }
  std::printf("%12s %12s %16s\n", "threshold", "escalate%", "worst_kept_rel");
  for (int t = 6; t <= 16; ++t) {
    long esc = 0;
    double worstKept = 0;
    for (int i = 0; i < kBins; ++i) {
      if (!s.n[i]) continue;
      if (i + 1 > t) esc += s.n[i];
      else if (s.mx[i] > worstKept) worstKept = s.mx[i];
    }
    char lab[16];
    std::snprintf(lab, sizeof lab, "1e-%d", t);
    std::printf("%12s %11.3f%% %16.3e\n", lab, 100.0 * esc / s.all, worstKept);
  }
}

}  // namespace

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s E_GeV pointsfile [n]\n", argv[0]);
    return 2;
  }
  const double energy = std::atof(argv[1]);
  const long limit = (argc > 3) ? std::atol(argv[3]) : 100000;

  std::vector<Point> pts;
  FILE* f = std::fopen(argv[2], "r");
  if (!f) { std::perror("fopen"); return 1; }
  char line[1024];
  while (std::fgets(line, sizeof line, f) && (long)pts.size() < limit) {
    if (line[0] == '#') continue;
    Point p;
    if (std::sscanf(line, "%lf %lf %lf %lf %lf", &p.a, &p.b, &p.c, &p.d, &p.e) == 5)
      pts.push_back(p);
  }
  std::fclose(f);

  // Reference, computed once.
  const auto pq = initDiffCross<__float128>(
      __float128(energy), __float128(5109991) / __float128(10000000));
  std::vector<double> ref(pts.size());
  for (size_t i = 0; i < pts.size(); ++i)
    ref[i] = static_cast<double>(
        diffCross<__float128>(pq, __float128(pts[i].a), __float128(pts[i].b),
                              __float128(pts[i].c), __float128(pts[i].d),
                              __float128(pts[i].e)).dsigma);

  report("double", std::numeric_limits<double>::epsilon(),
         measure<double>(pts, energy, ref));
  report("long double", std::numeric_limits<long double>::epsilon(),
         measure<long double>(pts, energy, ref));
  return 0;
}
