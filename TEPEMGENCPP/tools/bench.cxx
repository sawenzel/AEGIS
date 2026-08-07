// How much does __float128 actually cost here?
//
// The calibration says selective escalation saves little -- at production
// settings ~94% of sampled points need the wider type to hold the error under
// 1e-3. So the question stops being "when do we escalate" and becomes "can we
// simply always use quad", which is a question about wall-clock, not accuracy.
//
// Reports per-call cost for each instantiation and what that implies for
// TGenEpEmv1::Init(), which calls Diffcross once per trial and needed ~25.7k
// trials to converge at production settings.
#include "DiffCross.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace o2::aegis::tepemgen;

namespace
{
struct Point { double a, b, c, d, e; };

std::vector<Point> load(const char* path, long limit)
{
  std::vector<Point> v;
  FILE* f = std::fopen(path, "r");
  if (!f) { std::perror("fopen"); std::exit(1); }
  char line[1024];
  while (std::fgets(line, sizeof line, f) && (long)v.size() < limit) {
    if (line[0] == '#') continue;
    Point p;
    if (std::sscanf(line, "%lf %lf %lf %lf %lf", &p.a, &p.b, &p.c, &p.d, &p.e) == 5)
      v.push_back(p);
  }
  std::fclose(f);
  return v;
}
}  // namespace

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s E_GeV pointsfile [n]\n", argv[0]);
    return 2;
  }
  const double energy = std::atof(argv[1]);
  const long n = (argc > 3) ? std::atol(argv[3]) : 200000;
  const auto pts = load(argv[2], n);
  if (pts.empty()) { std::fprintf(stderr, "no points\n"); return 1; }

  const auto pd = initDiffCross<double>(energy, 5109991.0 / 10000000.0);
  const auto pq = initDiffCross<__float128>(
      __float128(energy), __float128(5109991) / __float128(10000000));
  const auto pl = initDiffCross<long double>(
      static_cast<long double>(energy), 5109991.0L / 10000000.0L);

  // Accumulate so the optimiser cannot discard the calls.
  volatile double sinkD = 0;
  volatile double sinkQ = 0;
  volatile double sinkL = 0;

  auto t0 = std::chrono::steady_clock::now();
  for (const auto& p : pts)
    sinkD = sinkD + diffCross<double>(pd, p.a, p.b, p.c, p.d, p.e).dsigma;
  auto t1 = std::chrono::steady_clock::now();
  for (const auto& p : pts)
    sinkQ = sinkQ + static_cast<double>(
        diffCross<__float128>(pq, p.a, p.b, p.c, p.d, p.e).dsigma);
  auto t2 = std::chrono::steady_clock::now();
  for (const auto& p : pts)
    sinkL = sinkL + static_cast<double>(
        diffCross<long double>(pl, p.a, p.b, p.c, p.d, p.e).dsigma);
  auto t3 = std::chrono::steady_clock::now();

  const double us_d = std::chrono::duration<double, std::micro>(t1 - t0).count() / pts.size();
  const double us_q = std::chrono::duration<double, std::micro>(t2 - t1).count() / pts.size();
  const double us_l = std::chrono::duration<double, std::micro>(t3 - t2).count() / pts.size();

  std::printf("points            : %zu\n", pts.size());
  std::printf("double            : %8.3f us/call\n", us_d);
  std::printf("long double       : %8.3f us/call   (%.1fx slower)\n", us_l, us_l / us_d);
  std::printf("__float128        : %8.3f us/call   (%.1fx slower)\n", us_q, us_q / us_d);
  std::printf("\nImplied cost of TGenEpEmv1::Init() at 25689 trials:\n");
  std::printf("  double          : %8.3f s\n", us_d * 25689 / 1e6);
  std::printf("  long double     : %8.3f s\n", us_l * 25689 / 1e6);
  std::printf("  __float128      : %8.3f s\n", us_q * 25689 / 1e6);
  std::printf("Worst case, the hard-coded fMaxXSTest = 1e7 trials:\n");
  std::printf("  double          : %8.1f s\n", us_d * 1e7 / 1e6);
  std::printf("  long double     : %8.1f s\n", us_l * 1e7 / 1e6);
  std::printf("  __float128      : %8.1f s\n", us_q * 1e7 / 1e6);
  return 0;
}
