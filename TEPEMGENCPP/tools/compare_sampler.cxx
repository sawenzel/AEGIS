// Phase 3 gate: the C++ sampler against epemgen.f.
//
// Two comparisons, in increasing order of what they prove:
//
//  1. STREAM IDENTITY. Both are driven from the same seeded gRandom, run one
//     after the other. If the C++ consumes random numbers in the same order as
//     the Fortran, the kinematics must agree to round-off event by event. This
//     is not the acceptance criterion -- the weights legitimately differ, since
//     the C++ evaluates the cross section in long double with escalation -- but
//     while the Fortran is still in the tree it localises any divergence to a
//     single event instead of leaving a histogram discrepancy to be explained.
//
//  2. STATISTICAL EQUIVALENCE, the actual gate: total cross section, and
//     two-sample Kolmogorov-Smirnov on the rapidity, log10(pt) and azimuth
//     distributions.
#include "EpEmSampler.h"

#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

extern "C" {
void ee_init_(double* ymin, double* ymax, double* ptmin, double* ptmax,
              double* cm_energy, double* z);
void ee_event_(double* ymin, double* ymax, double* ptmin, double* ptmax,
               double* ye, double* yp, double* xe, double* xp, double* phi,
               double* wtm2);
extern struct {
  double Xsect2, Dsect2, Xsecttot, Dsecttot;
  int Nevnt;
} eevent_;
}

using namespace o2::aegis::tepemgen;

namespace
{

struct Sample {
  std::vector<double> ye, yp, xe, xp, phi;
  long double sumW = 0;
  long n = 0;
};

/// Two-sample Kolmogorov-Smirnov. Returns D and the asymptotic p-value.
void ks(std::vector<double> a, std::vector<double> b, double& d, double& p)
{
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());
  size_t i = 0, j = 0;
  d = 0;
  while (i < a.size() && j < b.size()) {
    const double x = std::min(a[i], b[j]);
    while (i < a.size() && a[i] <= x) ++i;
    while (j < b.size() && b[j] <= x) ++j;
    d = std::max(d, std::fabs(double(i) / a.size() - double(j) / b.size()));
  }
  // D == 0 means the samples are identical, which is the *best* outcome here
  // once the streams match. The asymptotic series below is an alternating sum
  // that does not converge at lambda = 0, so it must be special-cased -- left
  // alone it reports p = 0, i.e. "DIFFERENT", for two identical samples.
  if (d == 0.0) { p = 1.0; return; }
  const double ne = std::sqrt(double(a.size()) * b.size() /
                              (double(a.size()) + b.size()));
  const double lam = (ne + 0.12 + 0.11 / ne) * d;
  p = 0;
  for (int k = 1; k <= 100; ++k)
    p += 2 * ((k % 2) ? 1 : -1) * std::exp(-2.0 * k * k * lam * lam);
  p = std::min(1.0, std::max(0.0, p));
}

void report(const char* name, const std::vector<double>& a,
            const std::vector<double>& b)
{
  double d, p;
  ks(a, b, d, p);
  std::printf("  %-10s KS D=%.5f  p=%.4f   %s\n", name, d, p,
              p > 0.01 ? "consistent" : "DIFFERENT");
}

}  // namespace

int main(int argc, char** argv)
{
  const long n = (argc > 1) ? std::atol(argv[1]) : 200000;
  const int seed = (argc > 2) ? std::atoi(argv[2]) : 12345;

  EpEmSampler::Config cfg;   // production settings by default
  double ymin = cfg.yMin, ymax = cfg.yMax;
  double ptmin = cfg.ptMinMeV, ptmax = cfg.ptMaxMeV;
  double energy = cfg.cmEnergyGeV, z = cfg.z;

  // --- Fortran -------------------------------------------------------------
  gRandom = new TRandom3(seed);
  ee_init_(&ymin, &ymax, &ptmin, &ptmax, &energy, &z);
  Sample F;
  for (long i = 0; i < n; ++i) {
    double ye, yp, xe, xp, phi, w;
    ee_event_(&ymin, &ymax, &ptmin, &ptmax, &ye, &yp, &xe, &xp, &phi, &w);
    F.ye.push_back(ye); F.yp.push_back(yp);
    F.xe.push_back(xe); F.xp.push_back(xp); F.phi.push_back(phi);
    F.sumW += w; ++F.n;
  }
  const double fortranXsec = eevent_.Xsecttot;

  // --- C++, same seed, same stream ----------------------------------------
  delete gRandom;
  gRandom = new TRandom3(seed);
  EpEmSampler s(cfg, [] {
    double r;
    do { r = gRandom->Rndm(); } while (r <= 0 || r >= 1);   // as eernd
    return r;
  });
  if (!s.init()) {
    std::fprintf(stderr, "C++ init failed: %s\n", s.error().c_str());
    return 1;
  }
  Sample C;
  for (long i = 0; i < n; ++i) {
    const auto e = s.next();
    C.ye.push_back(e.yElectron); C.yp.push_back(e.yPositron);
    C.xe.push_back(e.xElectron); C.xp.push_back(e.xPositron);
    C.phi.push_back(e.phi);
    C.sumW += e.weight; ++C.n;
  }

  // --- 1. stream identity --------------------------------------------------
  long firstDiff = -1;
  double worstK = 0;
  for (long i = 0; i < n; ++i) {
    const double dk = std::max(std::max(std::fabs(F.ye[i] - C.ye[i]),
                                        std::fabs(F.yp[i] - C.yp[i])),
                               std::max(std::fabs(F.xe[i] - C.xe[i]),
                                        std::fabs(F.phi[i] - C.phi[i])));
    if (dk > worstK) worstK = dk;
    if (dk > 1e-9 && firstDiff < 0) firstDiff = i;
  }
  std::printf("events                 : %ld  (seed %d)\n", n, seed);
  std::printf("\n1. stream identity (kinematics only; weights may differ)\n");
  std::printf("  worst |delta|        : %.3e\n", worstK);
  if (firstDiff >= 0) {
    std::printf("  first divergence at  : event %ld\n", firstDiff);
    std::printf("    fortran ye=%.10g yp=%.10g xe=%.10g phi=%.10g\n",
                F.ye[firstDiff], F.yp[firstDiff], F.xe[firstDiff],
                F.phi[firstDiff]);
    std::printf("    c++     ye=%.10g yp=%.10g xe=%.10g phi=%.10g\n",
                C.ye[firstDiff], C.yp[firstDiff], C.xe[firstDiff],
                C.phi[firstDiff]);
  } else {
    std::printf("  -> streams identical\n");
  }

  // --- 2. statistical equivalence -----------------------------------------
  std::printf("\n2. statistical equivalence (the acceptance gate)\n");
  const double cppXsec = s.xSection();
  std::printf("  x-section fortran    : %.6g kb  (%.6g b)\n", fortranXsec,
              fortranXsec * 1000);
  std::printf("  x-section c++        : %.6g kb  (%.6g b)\n", cppXsec,
              cppXsec * 1000);
  std::printf("  relative difference  : %+.4f%%\n",
              100.0 * (cppXsec - fortranXsec) / fortranXsec);
  std::printf("  escalated to quad    : %.2f%%  unphysical %lld\n",
              100.0 * s.cross().escalated() / std::max<long long>(1, s.cross().calls()),
              (long long)s.cross().unphysical());
  report("y(e-)", F.ye, C.ye);
  report("y(e+)", F.yp, C.yp);
  report("log10 pt", F.xe, C.xe);
  report("phi", F.phi, C.phi);
  return 0;
}
