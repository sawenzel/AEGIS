// Reference evaluator for the O2-6340 derivation check.
//
// Compiles the TEPEMGENCPP port headers (the validated C++ transcription of
// diffcross.f) and exposes, on the command line,
//   * the full 5-fold differential cross section at one kinematic point, and
//   * each individual Iz/Id/Iv integral function,
// evaluated in long double, which on aarch64 is IEEE binary128 (quad). That
// sidesteps the cancellation problem entirely, so this binary serves as the
// numerical oracle both for the first-principles amplitude check
// (verify_amplitude.py) and for the integral-identity check
// (verify_integrals.py).
//
// Build (from tools/):  g++ -O2 -std=c++17 -I ../include derivation_driver.cxx -o derivation_driver
//
// Usage:
//   driver point <energyGeV> <massMeV> <ppt> <yp> <pmt> <ym> <dphi>
//   driver iz0|iz1|iz2  x1 x2  u v
//   driver id0|id1|id2|id3  x1 x2 y1 y2  u v w
//   driver iv0|iv1|iv2  x1 x2 y1 y2 z1 z2  u v w w2
//   driver kin <energyGeV> <massMeV> <ppt> <yp> <pmt> <ym> <dphi>
//       -> prints the internal kinematic quantities (qb, m0, m1, md, mx, k1,
//          kd, kx and the longitudinal scalar products) for cross-checking
//          the derivation's formulas.

#include "DiffCross.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

using T = long double;
using o2::aegis::tepemgen::Vec2;

static T arg(char** argv, int i) { return strtold(argv[i], nullptr); }
static void out(T v) { printf("%.33Lg\n", v); }

int main(int argc, char** argv)
{
  using namespace o2::aegis::tepemgen;
  if (argc < 2) { fprintf(stderr, "usage: see header\n"); return 2; }
  const char* cmd = argv[1];
  bool ok = true;

  if (!strcmp(cmd, "point")) {
    auto phys = initDiffCross<T>(arg(argv, 2), arg(argv, 3));
    auto r = diffCross<T>(phys, arg(argv, 4), arg(argv, 5), arg(argv, 6),
                          arg(argv, 7), arg(argv, 8));
    printf("%.33Lg %d %.6Lg\n", r.dsigma, r.ok ? 1 : 0, r.survival);
    return 0;
  }
  if (!strcmp(cmd, "kin")) {
    // Reproduce the kinematic setup of diffCrossTerms so the derivation's
    // closed formulas can be checked against what the code actually uses.
    auto p = initDiffCross<T>(arg(argv, 2), arg(argv, 3));
    T ppvt = arg(argv, 4), yp = arg(argv, 5), pmvt = arg(argv, 6),
      ym = arg(argv, 7), dphi = arg(argv, 8);
    T m = p.m;
    T ppvl = sqrtl(ppvt * ppvt + m * m);
    T pmvl = sqrtl(pmvt * pmvt + m * m);
    T pmlxpml = pmvl * pmvl, pplxppl = ppvl * ppvl;
    T pmlxppl = pmvl * ppvl * coshl(yp - ym);
    T w2xpml = p.wl * pmvl * coshl(ym + p.wy);
    T w2xppl = p.wl * ppvl * coshl(yp + p.wy);
    T qb = p.gamma * (w2xppl + w2xpml) / sinhl(2 * p.wy);
    T qlxql = -qb * qb;
    T pmlxql = qb * pmvl * sinhl(p.wy - ym);
    T pplxql = qb * ppvl * sinhl(p.wy - yp);
    T m0 = -qlxql;
    T m1 = -qlxql - pplxppl - pmlxpml + 2 * (pplxql + pmlxql - pmlxppl);
    T md = m * m - (qlxql + pmlxpml - 2 * pmlxql);
    T mx = m * m - (qlxql + pplxppl - 2 * pplxql);
    printf("gamma %.21Lg\nbeta %.21Lg\nwy %.21Lg\nqb %.21Lg\n", p.gamma,
           p.beta, p.wy, qb);
    printf("m0 %.21Lg\nm1 %.21Lg\nmd %.21Lg\nmx %.21Lg\n", m0, m1, md, mx);
    printf("k1 %.21Lg %.21Lg\n", -ppvt - pmvt * cosl(dphi), -pmvt * sinl(dphi));
    printf("kd %.21Lg %.21Lg\n", -pmvt * cosl(dphi), -pmvt * sinl(dphi));
    printf("kx %.21Lg %.21Lg\n", -ppvt, 0.0L);
    return 0;
  }

  auto v2 = [&](int i) { return Vec2<T>{arg(argv, i), arg(argv, i + 1)}; };

  if (!strcmp(cmd, "iz0")) { out(Iz0<T>(v2(2), arg(argv, 4), arg(argv, 5), ok)); }
  else if (!strcmp(cmd, "iz1")) { out(Iz1<T>(v2(2), arg(argv, 4), arg(argv, 5), ok)); }
  else if (!strcmp(cmd, "iz2")) { out(Iz2<T>(v2(2), arg(argv, 4), arg(argv, 5), ok)); }
  else if (!strcmp(cmd, "id0")) { out(Id0<T>(v2(2), v2(4), arg(argv, 6), arg(argv, 7), arg(argv, 8), ok)); }
  else if (!strcmp(cmd, "id1")) { out(Id1<T>(v2(2), v2(4), arg(argv, 6), arg(argv, 7), arg(argv, 8), ok)); }
  else if (!strcmp(cmd, "id2")) { out(Id2<T>(v2(2), v2(4), arg(argv, 6), arg(argv, 7), arg(argv, 8), ok)); }
  else if (!strcmp(cmd, "id3")) { out(Id3<T>(v2(2), v2(4), arg(argv, 6), arg(argv, 7), arg(argv, 8), ok)); }
  else if (!strcmp(cmd, "iv0")) { out(Iv0<T>(v2(2), v2(4), v2(6), arg(argv, 8), arg(argv, 9), arg(argv, 10), arg(argv, 11), ok)); }
  else if (!strcmp(cmd, "iv1")) { out(Iv1<T>(v2(2), v2(4), v2(6), arg(argv, 8), arg(argv, 9), arg(argv, 10), arg(argv, 11), ok)); }
  else if (!strcmp(cmd, "iv2")) { out(Iv2<T>(v2(2), v2(4), v2(6), arg(argv, 8), arg(argv, 9), arg(argv, 10), arg(argv, 11), ok)); }
  else { fprintf(stderr, "unknown command %s\n", cmd); return 2; }
  if (!ok) fprintf(stderr, "setzero branch taken\n");
  return 0;
}
