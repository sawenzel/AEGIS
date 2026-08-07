// The sampling envelopes of ../TEPEMGEN/epemgen.f, and the normalisation
// integrals ee_init computes from them.
//
// ee_init integrates these -- NOT the exact cross section -- with CERNLIB
// Dtrint over two triangles. That is where the second, unreported bug lives
// (doc/02-findings.md section 4): the envelopes have derivative kinks, adaptive
// subdivision across a kink crawls, Dtrint hits KMX=35, calls MTLPRT and
// returns 0, and ee_init prints an error and carries on with XYsect = 0.
//
// This header removes the problem rather than tightening a tolerance. The two
// triangles per call union to the square [lo,hi]^2, and the integrand is a
// product f(Xp+Xe)*g(Xp-Xe). Under u = Xp+Xe, v = Xp-Xe (Jacobian 1/2) the
// square becomes a diamond whose half-width is L(u) = min(u-2lo, 2hi-u), so
//
//     I = 1/2 * integral_{2lo}^{2hi} f(u) * G(L(u)) du,
//     G(L) = integral_{-L}^{L} g(v) dv
//
// and G is available in closed form for both envelopes. A 2-D adaptive
// triangle integration collapses to one 1-D integral, which can then be split
// at the known kinks so every panel is smooth. dtrint.f (229 lines of
// goto-driven CERNLIB) and the MTLPRT dependency both go away.
#ifndef TEPEMGENCPP_ENVELOPES_H
#define TEPEMGENCPP_ENVELOPES_H

#include <cmath>

namespace o2::aegis::tepemgen
{

// --- parameters, verbatim from the DATA statements in epemgen.f -------------
//
// SINGLE PRECISION, and that is not a typo. The DATA statements carry no D0
// suffix -- `data parYpY / -4.8584, 0.11403E-01, ... /` -- so in Fortran these
// are REAL*4 literals widened to REAL*8 in the implicit real*8 context. The
// values the Fortran actually uses therefore differ from their decimal
// appearance in the 8th digit: -4.8584 is really -4.8583998680114746.
//
// f32() reproduces that. Physically it is irrelevant -- these are fit
// parameters good to four digits, and the envelope is only a proposal
// distribution whose imperfection the exact-cross-section weight corrects
// exactly. But without it the C++ and the Fortran diverge by ~1e-8 per event,
// which is far too large to dismiss as round-off and costs an afternoon to
// explain. Measured: it is the whole of the 7.3e-8 stream discrepancy.
//
// The commented-out first row of each DATA statement is the pre-1998 scaling;
// the active row is used here, as in the Fortran.

/// Round a decimal through REAL*4, as an unsuffixed Fortran literal does.
constexpr double f32(double x) { return static_cast<double>(static_cast<float>(x)); }
inline constexpr double kParXpX[5] = {f32(8.3668), f32(1.2004), f32(0.47225),
                                      f32(9.6951), f32(2.3814)};
inline constexpr double kParXmX[4] = {f32(9.5038), f32(39.040), f32(1.9492),
                                      f32(3.5660)};
inline constexpr double kParYpY[6] = {f32(-4.8584), f32(0.11403e-01), f32(4.0649),
                                      f32(0.38920e-04), f32(4.0283), f32(0.19238e-01)};
inline constexpr double kParYmY[6] = {f32(1.80), f32(8.), f32(1.7438),
                                      f32(0.17867), f32(24.188), f32(1.3097)};

/// Kink in DsdXpX; also the branch point in ee_event's Icase logic.
inline constexpr double kXpXKink = f32(0.6);
/// Kinks in DsdYmY.
inline constexpr double kYmYKink1 = f32(0.18);
inline constexpr double kYmYKink2 = f32(4.00);

// --- envelopes --------------------------------------------------------------

inline double dsdXpX(double x)
{
  return x < kXpXKink
             ? kParXpX[0] * std::exp(-kParXpX[1] * (x + kParXpX[2]) * (x + kParXpX[2]))
             : kParXpX[3] * std::exp(-kParXpX[4] * x);
}

inline double dsdXmX(double x)
{
  const double a = std::fabs(x);
  return kParXmX[0] * std::exp(-kParXmX[1] * a) +
         kParXmX[2] * std::exp(-kParXmX[3] * a);
}

inline double dsdYpY(double y)
{
  const double y2 = y * y;
  return kParYpY[0] * std::exp(-kParYpY[1] * y2) +
         kParYpY[2] * std::exp(-kParYpY[3] * y2 * y2) +
         kParYpY[4] * std::exp(-kParYpY[5] * y2);
}

inline double dsdYmY(double y)
{
  const double a = std::fabs(y);
  if (a < kYmYKink1) return kParYmY[0] * std::exp(-kParYmY[1] * a * a);
  if (a < kYmYKink2) return kParYmY[2] * std::exp(-kParYmY[3] * a * a);
  return kParYmY[4] * std::exp(-kParYmY[5] * a);
}

// --- closed-form inner integrals -------------------------------------------

/// integral_{-L}^{L} dsdXmX(v) dv, exactly. Sum of two symmetric exponentials.
inline double intXmX(double L)
{
  if (L <= 0) return 0.0;
  return 2.0 * (kParXmX[0] * (1.0 - std::exp(-kParXmX[1] * L)) / kParXmX[1] +
                kParXmX[2] * (1.0 - std::exp(-kParXmX[3] * L)) / kParXmX[3]);
}

/// integral_{-L}^{L} dsdYmY(v) dv, exactly.
///
/// Piecewise: two Gaussian pieces (hence erf) and an exponential tail. The
/// kinks at 0.18 and 4.00 are handled by construction here -- they are what
/// defeats the adaptive quadrature in the Fortran.
inline double intYmY(double L)
{
  if (L <= 0) return 0.0;
  const double rtpi = 1.7724538509055160273;  // sqrt(pi)
  // integral_0^t of c*exp(-k*x^2) = c*sqrt(pi)/(2*sqrt(k))*erf(sqrt(k)*t)
  auto gauss = [rtpi](double c, double k, double t) {
    return c * rtpi / (2.0 * std::sqrt(k)) * std::erf(std::sqrt(k) * t);
  };
  double half = 0.0;
  const double t1 = std::min(L, kYmYKink1);
  half += gauss(kParYmY[0], kParYmY[1], t1);
  if (L > kYmYKink1) {
    const double t2 = std::min(L, kYmYKink2);
    half += gauss(kParYmY[2], kParYmY[3], t2) -
            gauss(kParYmY[2], kParYmY[3], kYmYKink1);
  }
  if (L > kYmYKink2) {
    half += kParYmY[4] / kParYmY[5] *
            (std::exp(-kParYmY[5] * kYmYKink2) - std::exp(-kParYmY[5] * L));
  }
  return 2.0 * half;
}

}  // namespace o2::aegis::tepemgen

#endif  // TEPEMGENCPP_ENVELOPES_H
