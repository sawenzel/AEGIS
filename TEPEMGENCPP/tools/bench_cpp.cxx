// Time every C++ variant on the same points, one line each, machine-readable
// for tools/make_cpu_plot.py. Companion to bench_fortran.cxx, which times the
// Fortran in its two precisions.
//
// "adaptive" is the configuration that actually ships: long double working
// type, escalating to __float128 on either the cancellation monitor or
// disagreement with a narrower type. Its escalation rate is reported too,
// since the cost is meaningless without it.
#include "AdaptiveDiffCross.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace o2::aegis::tepemgen;

namespace
{
struct P { double a, b, c, d, e; };

template <typename F>
double timeIt(const std::vector<P>& pts, F&& f)
{
  volatile double sink = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (const auto& p : pts) sink = sink + f(p);
  return std::chrono::duration<double, std::micro>(
             std::chrono::steady_clock::now() - t0).count() / pts.size();
}
}  // namespace

int main(int argc, char** argv)
{
  const long limit = (argc > 2) ? std::atol(argv[2]) : 100000;
  std::vector<P> pts;
  FILE* f = std::fopen(argv[1], "r");
  if (!f) { std::perror("fopen"); return 1; }
  char line[1024];
  while (std::fgets(line, sizeof line, f) && (long)pts.size() < limit) {
    if (line[0] == '#') continue;
    P p;
    if (std::sscanf(line, "%lf %lf %lf %lf %lf", &p.a, &p.b, &p.c, &p.d, &p.e) == 5)
      pts.push_back(p);
  }
  std::fclose(f);
  if (pts.empty()) { std::fprintf(stderr, "no points\n"); return 1; }

  const double E = 5360.0;
  const auto pd = initDiffCross<double>(E, 5109991.0 / 10000000.0);
  const auto pl = initDiffCross<long double>((long double)E, 5109991.0L / 10000000.0L);
  const auto pq = initDiffCross<__float128>(
      __float128(E), __float128(5109991) / __float128(10000000));

  std::printf("cpp_double %.4f\n", timeIt(pts, [&](const P& p) {
    return (double)diffCross<double>(pd, p.a, p.b, p.c, p.d, p.e).dsigma;
  }));
  std::printf("cpp_longdouble %.4f\n", timeIt(pts, [&](const P& p) {
    return (double)diffCross<long double>(pl, (long double)p.a, (long double)p.b,
                                          (long double)p.c, (long double)p.d,
                                          (long double)p.e).dsigma;
  }));
  std::printf("cpp_float128 %.4f\n", timeIt(pts, [&](const P& p) {
    return (double)diffCross<__float128>(pq, __float128(p.a), __float128(p.b),
                                         __float128(p.c), __float128(p.d),
                                         __float128(p.e)).dsigma;
  }));

  AdaptiveDiffCross<> ad(E);
  const double us = timeIt(pts, [&](const P& p) {
    return ad(p.a, p.b, p.c, p.d, p.e);
  });
  std::printf("cpp_adaptive %.4f\n", us);
  std::fprintf(stderr, "adaptive escalated %.1f%% of %lld calls\n",
               100.0 * ad.escalated() / (ad.calls() ? ad.calls() : 1),
               (long long)ad.calls());
  return 0;
}
