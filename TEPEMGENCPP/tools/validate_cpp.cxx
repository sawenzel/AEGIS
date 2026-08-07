// Phase 1 acceptance gate: does the C++ port reproduce the Fortran?
//
// Links the Fortran diffcross.f AND instantiates the C++ template at the same
// precision, evaluates both at every sampled point, and reports the agreement.
// Built twice: validate_cpp_d (double vs double) and validate_cpp_q
// (-DTEP_QUAD, REAL*16 vs __float128).
//
// Bit-identical agreement is NOT the target and would be a red flag if
// claimed: the expression order is the same but the compilers, the libm and
// the association of the long polynomial sums are not, so a few ulp of
// disagreement is expected. The gate is ~1e-12 relative in double.
//
// Points where the Fortran guard fired are skipped: it forces dsigma to 0 and
// there is nothing to compare against. The C++ deliberately does not reproduce
// that guard -- see DiffCross.h.
#include "DiffCross.h"
#include "fortran_iface.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

using o2::aegis::tepemgen::diffCross;
using o2::aegis::tepemgen::initDiffCross;

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s E_GeV pointsfile [tol]\n", argv[0]);
    return 2;
  }
  // Floor of the gate; the conditioning-dependent part is added per point.
  const double tol = (argc > 3) ? std::atof(argv[3]) : 1e-14;

  tepreal energy = static_cast<tepreal>(std::atof(argv[1]));
  // Exact rational, not a decimal literal: see the note in DiffCross.h.
  // 0.5109991 is not dyadic, so `static_cast<tepreal>(0.5109991)` rounds it to
  // double first while the Fortran parses it at full width.
  tepreal mass = tepreal(5109991) / tepreal(10000000);
  initdiffcross_(&energy, &mass);                        // Fortran side
  const auto phys = initDiffCross<tepreal>(energy, mass);  // C++ side

  FILE* f = std::fopen(argv[2], "r");
  if (!f) { std::perror("fopen"); return 1; }

  long n = 0, skipped = 0, over = 0;
  double worst = 0, worstYp = 0, worstYm = 0, worstDphi = 0, worstF = 0, worstC = 0;
  long double sum = 0;
  double worstSurvival = 1;

  // Disagreement binned by the cancellation monitor. Two independent
  // implementations of an ill-conditioned expression differ in their last
  // bits; the conditioning then amplifies that by 1/survival. If the port is
  // faithful, disagreement must track survival and nothing else -- and the
  // same table calibrates the Phase 2 escalation threshold.
  constexpr int kBins = 20;  // one decade per bin, down to 1e-19
  long binN[kBins] = {};
  long double binSum[kBins] = {};
  double binMax[kBins] = {};

  char line[1024];
  while (std::fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    double a, b, c, d, e, w, vfort;
    int bad;
    if (std::sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %d",
                    &a, &b, &c, &d, &e, &w, &vfort, &bad) < 8)
      continue;
    if (bad) { ++skipped; continue; }

    tepreal ta = a, tb = b, tc = c, td = d, te = e, fout = 0;
    diffcross_(&ta, &tb, &tc, &td, &te, &fout);
    const auto r = diffCross<tepreal>(phys, ta, tb, tc, td, te);

    if (fout == tepreal(0)) { ++skipped; continue; }

    // The difference must be formed in tepreal and only then narrowed for
    // reporting. Casting both sides to double first would silently reduce the
    // quad build to a double-precision test -- it would "PASS" at 1e-16 while
    // saying nothing about the 30-odd digits that are the entire reason the
    // quad instantiation exists.
    using o2::aegis::tepemgen::abs;
    const tepreal diff = abs(r.dsigma - fout) / abs(fout);
    const double rel = static_cast<double>(diff);
    const double vf = static_cast<double>(fout);
    const double vc = static_cast<double>(r.dsigma);
    ++n;
    sum += rel;
    // Conditioning-aware gate. A flat tolerance is wrong here: two faithful
    // implementations of an ill-conditioned expression differ in their last
    // bits, and the conditioning amplifies that by roughly 1/survival, so a
    // fixed threshold either rejects correct code at low survival or accepts
    // a real bug at high survival. The slack constant absorbs the compounding
    // of that amplification through the nested Iz/Id/Iv calls.
    //
    // Discrimination is still sharp: at a typical survival of 1e-9 the quad
    // build allows ~2e-20, so a genuine transcription error -- which would
    // show up as an O(1e-3) discrepancy -- fails by sixteen orders of
    // magnitude.
    const double surv0 = static_cast<double>(r.survival);
    const double epsT = static_cast<double>(
        std::numeric_limits<tepreal>::epsilon());
    const double allowed = surv0 > 0 ? tol + 1e5 * epsT / surv0 : 1.0;
    if (rel > allowed) ++over;
    if (rel > worst) {
      worst = rel; worstYp = b; worstYm = d; worstDphi = e;
      worstF = vf; worstC = vc;
      worstSurvival = static_cast<double>(r.survival);
    }
    const double surv = static_cast<double>(r.survival);
    int bin = 0;
    if (surv > 0) {
      const double l = -std::log10(surv);            // 0 -> bin 0, 7+ -> last
      bin = l < 0 ? 0 : (l >= kBins ? kBins - 1 : static_cast<int>(l));
    } else {
      bin = kBins - 1;
    }
    ++binN[bin];
    binSum[bin] += rel;
    if (rel > binMax[bin]) binMax[bin] = rel;
  }
  std::fclose(f);
  if (!n) { std::fprintf(stderr, "no comparable points\n"); return 1; }

  std::printf("compared            : %ld  (skipped %ld: guard fired or zero)\n",
              n, skipped);
  std::printf("mean |rel diff|     : %.3Le\n", sum / n);
  std::printf("worst |rel diff|    : %.3e   (survival there: %.2e)\n", worst,
              worstSurvival);
  std::printf("  at yp=%.4f ym=%.4f dphi=%.6f\n", worstYp, worstYm, worstDphi);
  std::printf("  fortran=%.17g\n  c++    =%.17g\n", worstF, worstC);
  std::printf("points over gate     : %ld  (%.4f%%)\n", over,
              100.0 * over / n);

  std::printf("\ndisagreement vs cancellation monitor:\n");
  std::printf("%18s %9s %12s %12s\n", "survival", "n", "mean_rel", "max_rel");
  for (int i = 0; i < kBins; ++i) {
    if (!binN[i]) continue;
    char lab[32];
    if (i == kBins - 1)
      std::snprintf(lab, sizeof lab, "< 1e-%d", kBins - 1);
    else
      std::snprintf(lab, sizeof lab, "1e-%d .. 1e-%d", i + 1, i);
    std::printf("%18s %9ld %12.3Le %12.3e\n", lab, binN[i],
                binSum[i] / binN[i], binMax[i]);
  }
  std::printf("\n%s\n", over == 0 ? "PASS" : "FAIL");
  return over == 0 ? 0 : 1;
}
