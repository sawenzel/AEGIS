// Compare precision strategies on the quantities the generator actually uses:
// the cross section (sum of dsigma) and the variance estimate (sum of
// dsigma^2), whose ratio drives the convergence test in
// TGenEpEmv1::CalcXSection and hence the abort path in O2-6340.
//
// Reference is __float128 throughout. Reports bias, escalation rate and
// wall-clock for each candidate, so the trade is visible rather than asserted.
#include "AdaptiveDiffCross.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace o2::aegis::tepemgen;

namespace
{
struct Point { double a, b, c, d, e; };

struct Totals {
  long double sum = 0, sum2 = 0;
  double seconds = 0;
  long long escalated = 0;
};

void row(const char* name, const Totals& t, const Totals& ref, long n)
{
  const double biasSum = static_cast<double>(100.0L * (t.sum - ref.sum) / ref.sum);
  const double biasVar = static_cast<double>(100.0L * (t.sum2 - ref.sum2) / ref.sum2);
  std::printf("%-26s %+11.5f%% %+11.5f%% %9.2f%% %10.2f %8.1f\n", name, biasSum,
              biasVar, 100.0 * t.escalated / n, t.seconds * 1e6 / n,
              t.seconds / ref.seconds * 100.0);
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
  const long n = pts.size();
  if (!n) { std::fprintf(stderr, "no points\n"); return 1; }

  auto run = [&](auto&& eval) {
    Totals t;
    const auto t0 = std::chrono::steady_clock::now();
    for (const auto& p : pts) {
      const double v = eval(p.a, p.b, p.c, p.d, p.e);
      t.sum += v;
      t.sum2 += static_cast<long double>(v) * v;
    }
    t.seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    return t;
  };

  // Reference: quad everywhere.
  AdaptiveDiffCross<__float128, __float128> quad(energy, __float128(2));  // never escalates
  Totals ref = run([&](double a, double b, double c, double d, double e) {
    return quad(a, b, c, d, e);
  });
  ref.escalated = 0;

  std::printf("points: %ld   reference: __float128 throughout\n\n", n);
  std::printf("%-26s %12s %12s %10s %10s %8s\n", "strategy", "xsec bias",
              "var bias", "escalated", "us/call", "%ref");

  {
    AdaptiveDiffCross<double, double> d(energy, double(-1));  // never escalates
    Totals t = run([&](double a, double b, double c, double dd, double e) {
      return d(a, b, c, dd, e);
    });
    row("double only", t, ref, n);
  }
  {
    AdaptiveDiffCross<long double, long double> l(energy, (long double)-1);
    Totals t = run([&](double a, double b, double c, double dd, double e) {
      return l(a, b, c, dd, e);
    });
    row("long double only", t, ref, n);
  }
  for (int t10 = 12; t10 <= 15; ++t10) {
    AdaptiveDiffCross<long double, __float128> ad(
        energy, std::pow((long double)10, -t10));
    Totals t = run([&](double a, double b, double c, double dd, double e) {
      return ad(a, b, c, dd, e);
    });
    t.escalated = ad.escalated();
    char lab[64];
    std::snprintf(lab, sizeof lab, "long double -> quad @1e-%d", t10);
    row(lab, t, ref, n);
  }
  {
    AdaptiveDiffCross<double, __float128> ad(energy, 1e-7);
    Totals t = run([&](double a, double b, double c, double dd, double e) {
      return ad(a, b, c, dd, e);
    });
    t.escalated = ad.escalated();
    row("double -> quad @1e-7", t, ref, n);
  }
  return 0;
}
