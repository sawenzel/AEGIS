// A small adaptive Gauss-Kronrod (G7-K15) integrator, templated on the scalar.
//
// Replaces CERNLIB Dgauss (D103, 6 calls in ee_init) and, together with the
// analytic reduction in Envelopes.h, Dtrint (D105, 229 lines of goto-driven
// CERNLIB plus the MTLPRT dependency).
//
// Why not GSL or ROOT::Math::Integrator: both are double-only. That would break
// the precision templating which is the whole point of this port, and pull in
// MathMore besides. Sixty lines of Gauss-Kronrod is the smaller dependency.
//
// Unlike Dtrint this reports failure instead of returning 0 and printing --
// ee_init currently carries on with XYsect = 0 after such a "failure", which
// silently zeroes every event weight.
#ifndef TEPEMGENCPP_QUADRATURE_H
#define TEPEMGENCPP_QUADRATURE_H

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <vector>

namespace o2::aegis::tepemgen
{

/// Gauss-Kronrod 7-15 nodes and weights on [-1,1], symmetric half only.
namespace detail
{
inline constexpr double kGK15x[8] = {
    0.000000000000000000, 0.207784955007898468, 0.405845151377397167,
    0.586087235467691130, 0.741531185599394440, 0.864864423359769073,
    0.949107912342758525, 0.991455371120812639};
inline constexpr double kGK15w[8] = {
    0.209482141084727828, 0.204432940075298892, 0.190350578064785410,
    0.169004726639267903, 0.140653259715525919, 0.104790010322250184,
    0.063092092629978553, 0.022935322010529225};
/// Gauss-7 weights, on the odd-indexed Kronrod nodes.
inline constexpr double kG7w[4] = {0.417959183673469388, 0.381830050505118945,
                                   0.279705391489276668, 0.129484966168869693};
}  // namespace detail

template <typename T>
struct QuadResult {
  T value{};
  T error{};
  bool converged = false;
  int intervals = 0;
};

/// One G7-K15 panel: returns the Kronrod estimate and |K - G| as the error.
template <typename T, typename F>
void gaussKronrod15(F&& f, T a, T b, T& value, T& error)
{
  const T c = (a + b) / T(2);
  const T h = (b - a) / T(2);
  T k = T(0), g = T(0);
  for (int i = 0; i < 8; ++i) {
    const T x = T(detail::kGK15x[i]);
    const T fs = (i == 0) ? f(c) : f(c - h * x) + f(c + h * x);
    k += T(detail::kGK15w[i]) * fs;
    if (i % 2 == 0) g += T(detail::kG7w[i / 2]) * fs;
  }
  value = k * h;
  const T d = (k - g) * h;
  error = d < T(0) ? -d : d;
}

/// Adaptive integration over [a,b], bisecting the worst panel each round.
///
/// `breakpoints` are inserted as initial panel edges. Passing the known kinks
/// of the integrand there is not an optimisation but the fix: adaptive
/// subdivision across a C0 kink converges so slowly that Dtrint gives up.
template <typename T, typename F>
QuadResult<T> integrate(F&& f, T a, T b,
                        std::initializer_list<T> breakpoints = {},
                        T relTol = T(1e-10), int maxIntervals = 2000)
{
  struct Panel { T a, b, v, e; };
  std::vector<Panel> panels;

  std::vector<T> edges{a};
  for (T p : breakpoints)
    if (p > a && p < b) edges.push_back(p);
  edges.push_back(b);
  for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
    if (!(edges[i + 1] > edges[i])) continue;
    Panel p{edges[i], edges[i + 1], T(0), T(0)};
    gaussKronrod15(f, p.a, p.b, p.v, p.e);
    panels.push_back(p);
  }

  QuadResult<T> r;
  while (static_cast<int>(panels.size()) < maxIntervals) {
    T total{}, err{};
    std::size_t worst = 0;
    for (std::size_t i = 0; i < panels.size(); ++i) {
      total += panels[i].v;
      err += panels[i].e;
      if (panels[i].e > panels[worst].e) worst = i;
    }
    const T mag = total < T(0) ? -total : total;
    if (err <= relTol * mag || err == T(0)) {
      r.value = total; r.error = err; r.converged = true;
      r.intervals = static_cast<int>(panels.size());
      return r;
    }
    Panel p = panels[worst];
    const T m = (p.a + p.b) / T(2);
    Panel l{p.a, m, T(0), T(0)}, rr{m, p.b, T(0), T(0)};
    gaussKronrod15(f, l.a, l.b, l.v, l.e);
    gaussKronrod15(f, rr.a, rr.b, rr.v, rr.e);
    panels[worst] = l;
    panels.push_back(rr);
  }

  for (const auto& p : panels) { r.value += p.v; r.error += p.e; }
  r.converged = false;
  r.intervals = static_cast<int>(panels.size());
  return r;
}

}  // namespace o2::aegis::tepemgen

#endif  // TEPEMGENCPP_QUADRATURE_H
