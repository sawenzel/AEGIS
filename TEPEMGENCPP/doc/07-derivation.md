# Derivation of the TEPEMGEN differential cross section (O2-6340, checkbox 1)

2026-08-08. This document derives, from the QED Feynman rules in the
external-field approximation, the five-fold differential cross section that
`diffcross.f` (AEGIS/TEPEMGEN) and its C++ port `DiffCross.h`
(this directory) implement, and verifies the
implementation against the derivation numerically. It addresses the first
open checkbox of [O2-6340] — "review of the underlying physics" — for the
part that can be settled independently: that the code computes the
lowest-order two-photon amplitude of Alscher, Hencken, Trautmann and Baur,
Phys. Rev. **A55** (1997) 396, with the stated kinematic conventions and
normalisation.

**Result.** The code's `dsigma` equals, to 1e-12 relative at generic
phase-space points (quadrature-limited),

```
dsigma = (4/beta^2) (2pi)^-4 (1/4) (hbar c)^2 *
         Integral d^2k  Sum_spins | ubar(p-) Gamma(k) v(p+) |^2 / (q1^2 q2^2)^2
```

with the two-diagram vertex function Gamma defined in §3 — i.e. the 18-term
sum `N1..N18` *is* the spin-summed square of the external-field Born
amplitude, term-for-term normalisation included. Verified by
`tools/verify_amplitude.py` against `tools/derivation_driver` (the C++ port in quad);
the check reuses **nothing** of the code's algebra — explicit Dirac spinors
and numerical 2-D integration only.

## 1. Setup and notation

Metric (+,−,−,−), ħ = c = 1, lepton mass m (MeV). Two point-like nuclei of
charge Ze move on straight lines along ±z with velocity β and Lorentz factor
γ = E_NN/(2 m_p) (the code's `gamma = energy/(2*0.938)`), at transverse
impact parameters b₁, b₂. Their *coordinate* velocities are

```
w1 = (1, 0, 0, +beta),   w2 = (1, 0, 0, -beta)
```

these are exactly the code's `w1, w2`:

```
w1·w1 = w2·w2 = 1 − beta² = 1/gamma²          (code: w1xw1, w2xw2)
w1·w2 = 1 + beta² = 2 − 1/gamma²              (code: w1xw2)
```

Writing longitudinal 2-vectors in hyperbolic form, `w1 = wl (cosh wy, sinh wy)`
with `wl = 1/gamma`, `wy = arcosh gamma` (the beam rapidity): the code's `wl`,
`wy`.

The produced leptons have transverse momenta **p**₊⊥ = (p₊t, 0),
**p**₋⊥ = p₋t (cos Δφ, sin Δφ), rapidities y±, transverse masses
m⊥± = √(p±t² + m²). Every four-vector product in the code is split into a
*longitudinal* (t,z)-Minkowski part (suffix `l`) and a *transverse* part
(suffix `t`, carrying the metric's minus sign: `pmtxpmt = −p₋t²`). The
longitudinal products are evaluated in hyperbolic closed form
(`pmlxppl = m⊥₊ m⊥₋ cosh(y₊−y₋)` etc.); §7 explains why that is essential
and not a stylistic choice.

## 2. The fields of the nuclei

The current of nucleus j is J^μ(x) = Ze w_j^μ δ²(x⊥ − b_j) δ(z − β_j t).
With the convention Ã(q) = ∫d⁴x e^{iq·x} A(x) and □A = J (Lorenz gauge,
Heaviside–Lorentz, e² = 4πα):

```
J̃_j^μ(q) = 2π Ze w_j^μ δ(q·w_j) e^{−i q⊥·b_j}
Ã_j^μ(q) = − (1/q²) J̃_j^μ(q)
```

The delta function is the covariant statement that a boosted Coulomb field
transfers no energy in its own rest frame; it is what makes the longitudinal
photon kinematics *fixed* rather than integrated (§4).

## 3. The second-order amplitude

To lowest order the pair e⁻(p₋,s₋) e⁺(p₊,s₊) is produced by absorbing one
photon from each nucleus. The second-order S-matrix element, with
P = p₊ + p₋ and q the momentum supplied by nucleus 1, is

```
S_fi = i e² ∫ d⁴q/(2π)⁴  ū(p₋) { A̸₁(q)   S_F(p₋−q)  A̸₂(P−q)
                                + A̸₂(P−q) S_F(q−p₊)  A̸₁(q) } v(p₊)
```

(the two orderings of the two vertices; S_F(p) = (p̸+m)/(p²−m²)). Inserting
the fields and defining the vertex function

```
Gamma(q) = w̸₁ (p̸₋ − q̸ + m) w̸₂ / [(p₋−q)² − m²]
         + w̸₂ (q̸ − p̸₊ + m) w̸₁ / [(q−p₊)² − m²]
```

gives, after the two delta functions collapse the longitudinal q-integrals
(Jacobian: ∂(q·w₁, q·w₂)/∂(q⁰,q_z) = 2β),

```
S_fi = i e²(Ze)²/(2β) · e^{−iP⊥·b₂} ∫ d²k/(2π)²  e^{−ik·b}  T(k),
T(k) = ū(p₋) Gamma(q(k)) v(p₊) / (q² (P−q)²),      b = b₁ − b₂,
```

where k = q⊥ and the longitudinal part of q is fixed (§4).

## 4. The longitudinal photon kinematics — the code's `qb`, `q_l`

The constraints are q·w₁ = 0 and (P−q)·w₂ = 0, i.e. q·w₂ = P·w₂. The unique
longitudinal 2-vector orthogonal to w₁ is ∝ (sinh wy, cosh wy), so

```
q_l = qb (sinh wy, cosh wy) = qb (γβ, γ),
q_l·w₂ = qb wl sinh 2wy = 2qb γβ   ⇒   qb = P·w₂ / (2γβ)
```

which is the code's `qb = gamma*(w2xppl+w2xpml)/sinh(2*wy)` exactly:
P·w₂ = w2xppl + w2xpml (w₂ *is* the coordinate velocity) and
sinh 2wy = 2 sinh wy cosh wy = 2γ²β, so γ/sinh 2wy = 1/(2γβ). It is
spacelike: q_l·q_l = −qb². The derived products

```
q_l·p₋_l = qb m⊥₋ sinh(wy − y₋)     (code: pmlxql)
q_l·p₊_l = qb m⊥₊ sinh(wy − y₊)     (code: pplxql)
q_l·w₂  = w2xpml + w2xppl           (code: w2xql),   q_l·w₁ = 0
```

match `Diffcross` line for line. By the mirrored argument for nucleus 2,
q₂ = P − q has fixed longitudinal part with

```
qb₂ = P·w₁ / (2γβ),        q₂_l·q₂_l = −qb₂².
```

## 5. The four denominators — the code's `m0, m1, md, mx, k1, kd, kx`

With q = (q_l, k), all four denominators are negatives of *positive-definite*
quadratics in k:

```
q²          = q_l² − k²             = −(k² + m0),          m0 = qb²
(P−q)²                              = −((k+k1)² + m1),     k1 = −P⊥,  m1 = qb₂²
(p₋−q)² − m² = q² − 2 q·p₋          = −((k+kd)² + md),     kd = −p₋⊥
(q−p₊)² − m² = q² − 2 q·p₊          = −((k+kx)² + mx),     kx = −p₊⊥
md = m² + qb² − m⊥₋² + 2 q_l·p₋_l,   mx = m² + qb² − m⊥₊² + 2 q_l·p₊_l
```

These are exactly the code's `m0`, `m1` (expanded form), `md`, `mx` and shift
vectors. Two useful facts the code does not state:

- `m1 = qb₂²` — photon 2's own `qb` squared. (Proof: q₂_l ⊥ w₂ and
  q₂_l·w₁ = P·w₁, mirroring §4.) Verified numerically at 1e-14.
- All four "masses" are positive on the physical domain (md > 0 because
  y₋ < wy always at the cuts in use), so the 2-D integrals have no poles and
  need no iε — the amplitude at fixed pair kinematics is real up to a global
  phase.

## 6. Impact-parameter integration: why single 2-D integrals suffice

The generator's cross section is the impact-parameter integral of the pair
production probability, σ = ∫d²b P(b). By Parseval,

```
∫d²b | ∫d²k/(2π)² e^{−ik·b} T(k) |²  =  (1/(2π)²) ∫d²k |T(k)|²
```

— the b-integral diagonalises the double transverse integral. This is why
`Diffcross` contains only *single* ∫d²k integrals of products of the four
denominators (§5), with total exponent up to (2,2,1,1): the two photon
propagators appear squared (|amplitude|²), the lepton propagators once per
diagram factor (direct², crossed², interference).

## 7. Assembling the cross section; the exact normalisation chain

With relativistic normalisation, dΠ = d³p₊ d³p₋ / ((2π)⁶ 4E₊E₋) and
d³p/(2E) = ½ p_t dp_t dφ dy, integrating the overall azimuth (2π):

```
dσ/(dp₊t dp₋t dy₊ dy₋ dΔφ) / (p₊t p₋t)
  = |e²(Ze)²/(2β)|² (2π)^−2 (2π)^−6 (1/4) ∫d²k Σ_spins |T(k)|²
  = (Zα)⁴ · (4/β²) · (2π)^−4 · (1/4) · ∫d²k Σ_spins |T(k)|²
```

using e⁴(Ze)⁴ = (4π)⁴(Zα)⁴ and (4π)⁴/(2π)⁴ = 16. This is *exactly* the
code's final factor chain

```fortran
NT=NT*4D0/beta**2          ! (4pi)^4 / (2pi)^4 / 4beta^2  ("w/u" comment)
NT=NT/(2*pi)**6*(2*pi)**2  ! two d^3p/(2pi)^3 and the Parseval (2pi)^-2
NT=NT/4D0                  ! 1/2E+ 1/2E-  ->  (1/2)(1/2) in the y,pt measure
NT=NT*(1.9733D0)**2/10D0   ! MeV^-2 -> kbarn: (hbar c)^2 = 0.38939 kb MeV^2
```

so the identification is `NT = ∫d²k Σ_spins |T(k)|²`, with the p_t factors
and the global 2π left in the integration measure as the header comment says.
The `(Zα)⁴` stays outside, as documented (`dsigma` is per (Zα)⁴ — this is
also why Z never appears in `Diffcross` and one table serves every collision
system).

## 8. The reduction to Iz/Id/Iv, and what each integral is

Σ_spins |ū Γ v|² is a Dirac trace: a polynomial in scalar products of
{w₁, w₂, p₊, p₋, q}, divided by the denominators of §5. Splitting each
scalar product into its fixed longitudinal part plus k-dependent transverse
part and reducing the k-polynomial numerators against the denominators
(2-D Passarino–Veltman-style) expresses NT as coefficient polynomials in the
*longitudinal* products times ten scalar 2-D integrals. Those are the
`N1..N18` terms. The ten scalar integrals were identified numerically
(`tools/verify_integrals.py`, agreement 1.7e-14, exponent tuple unique in
each search space):

```
Iz0(x;u,v)        = ∫d²k  (k²+u)^−1 ((k+x)²+v)^−1
Iz1               =        (k²+u)^−2 ((k+x)²+v)^−1        = −∂u Iz0
Iz2               =        (k²+u)^−2 ((k+x)²+v)^−2        = ∂u∂v Iz0
Id0(x,y;u,v,w)    = ∫d²k  (k²+u)^−1 ((k+x)²+v)^−1 ((k+y)²+w)^−1
Id1               =        (k²+u)^−2 (…)^−1 (…)^−1        = −∂u Id0
Id2               =        (k²+u)^−2 ((k+x)²+v)^−2 (…)^−1 = ∂u∂v Id0
Id3               =        all three squared              = −∂u∂v∂w Id0
Iv0(x,y,z;u,v,w,w2)= ∫d²k (k²+u)^−1 (…)^−1 (…)^−1 (…)^−1
Iv1               =        (k²+u)^−2, rest single         = −∂u Iv0
Iv2               =        (k²+u)^−2 ((k+x)²+v)^−2, rest single
```

The family is systematic: Iz/Id/Iv = 2/3/4 denominators, index = number of
squared ones. The base integral has the classic Feynman-parameter closed
form: with t² ≡ x²,

```
Iz0 = ∫₀¹ dt / [t(1−t)x² + (1−t)u + t v] · π
    = (π/s) log[ (x²+u+v+s)² / (4uv) ],    s = √((x²+u+v)² − 4uv)
```

(the discriminant identity (x²+v−u)² + 4x²u = (x²+u+v)² − 4uv gives the
code's `s`). Everything else follows by differentiating under the integral
sign and by the partial-fraction recursions the code implements (`Id0` as a
combination of three `Iz0`, `Iv0` as a combination of four `Id0`, …).

**The Gram-determinant origin of the Δφ → π instability.** Every `Id`/`Iv`
reduction divides by `B = 4u·axy² + A`, where `axy = x₁y₂ − x₂y₁` is the
cross product of two shift vectors — a 2-D Gram determinant. The shifts are
built from p₊⊥ and p₋⊥, so at Δφ = π (and at Δφ = 0) they are collinear and
every axy vanishes; `Id3` divides by B³. The *integrals themselves* stay
finite — only the reduction representation degenerates. This is the textbook
small-Gram-determinant instability of tensor reduction, and it is the
analytic reason the port's empirical finding (02-findings.md §2) shows
conditioning degrading as Δφ → π: the catastrophic cancellation lives in the
reduced representation, not in the physics. A future analytic reformulation
should target exactly these B-denominators (e.g. expansions around the
collinear limit, or a rotated reduction basis), which is the concrete content
of the "restructure the terms symbolically" idea noted in 01-scope.md.

## 9. Numerical verification

Two scripts in `tools/`, both using `derivation_driver.cxx` (the C++ port compiled at
`long double` = IEEE quad on aarch64) as the reference:

- `verify_integrals.py` — identifies each Iz/Id/Iv closed form against its
  defining integral (table in §8), generic O(1) arguments, agreement
  **1.7e-14**, unique in the exponent search space.
- `verify_amplitude.py` — the full first-principles check: explicit Dirac
  spinors (Dirac representation, ūu = 2m), the Γ of §3, stable denominators
  of §5, two-patch log-polar quadrature (one patch per squared photon
  propagator: peaks of width qb at k = 0 and qb₂ at k = P⊥). Agreement with
  the port's quad `dsigma` at six production-phase-space points:
  **1e-12 … 1e-16** (5 points), 4e-9 (one point, quadrature-limited).

Two numerical remarks worth keeping:

- The check itself initially disagreed at 1.6e-6 at forward rapidity because
  it computed (P−q)² componentwise: P⁰ ≈ q⁰ ≈ 224 MeV while the difference
  is ≈ 0.013 MeV. The fix is the hyperbolic closed forms of §5 — the same
  device the 1997 Fortran uses. At LHC γ the *inputs* to the term polynomials
  are only obtainable stably this way; the instability O2-6340 is about lives
  further downstream, in the N1..N18 summation (§8).
- The two-photon-pole structure (widths qb, qb₂ ~ 1e-4 MeV at PbPb) is also
  why any *numerical* cross-check of this integral needs the pole positions
  fed to the quadrature; a naive 2-D grid misses the k = P⊥ spike entirely.

## 10. What is established, and what deliberately is not

Established:

- `Diffcross` computes the lowest-order (one photon per nucleus, "Born")
  external-field amplitude for point-like nuclei, with exact longitudinal
  kinematics and exact normalisation — §§2-7 constitute the derivation and
  the numerical identity of §9 confirms the 18-term sum against it.
- The impact-parameter-integrated form (§6): what is computed is
  σ = ∫d²b P₁(b) for single-pair production; the generator's Poisson
  machinery for multiple pairs sits on top, outside `Diffcross`.
- The scalar-integral basis and its closed forms (§8).

Not established here, by choice:

- **Term-by-term symbolic identity of N1..N18** with the reduced trace. The
  numeric identity of the sum at machine precision across phase space makes
  a term-level error that cancels pointwise implausible; a symbolic
  re-reduction (sympy gamma algebra + the §8 recursions) is the natural next
  step if wanted, and would also produce the restructured, Gram-stable form.
- **Physics beyond the code's scope**: Coulomb corrections (higher orders in
  Zα — known to reduce the total by ~10-20% at LHC), nuclear form factors
  (point-like here; the pt cuts make this defensible for the QED-background
  use), and the equivalent-single-photon vs. two-photon interference terms
  that vanish for straight-line trajectories. These are properties of the
  Alscher et al. *model*, shared by the Fortran and the port alike; changing
  them is a physics decision on O2-6340, not a porting question.

[O2-6340]: https://its.cern.ch/jira/browse/O2-6340
