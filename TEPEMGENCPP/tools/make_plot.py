#!/usr/bin/env python3
"""Reproduce the O2-6340 figure, comparing the old and new implementations.

The plot in the issue showed the differential cross section against electron
rapidity computed twice, in double and in REAL*16, with the double curve
collapsing into noise over part of the range. This is the same comparison with
the ported code added, and with the residual error shown underneath so the
claim is visible rather than asserted.

The reference curve is the ORIGINAL Fortran recompiled with
`gfortran -freal-8-real-16`, not the C++ quad instantiation. That keeps ground
truth independent of the thing being validated.

Usage:  make_plot.py <fort_d.dat> <fort_q.dat> <cpp.dat> <out.png>
  fort_d/fort_q : ym  dsigma  badcount_delta          (tools/scan_rapidity_*)
  cpp           : ym  double  longdouble  adaptive  quad  survival  (tools/scan_cpp)
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path, cols):
    out = [[] for _ in cols]
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        for k, c in enumerate(cols):
            out[k].append(float(f[c]))
    return out


def main():
    fd, fq, cpp, png = sys.argv[1:5]
    y, old = load(fd, [0, 1])
    _, ref = load(fq, [0, 1])
    _, new = load(cpp, [0, 3])

    # A log axis cannot show the sign, and the broken curve goes negative about
    # as often as it goes positive -- which is the whole reason the Fortran's
    # `IF(dsigma.LT.0)` guard catches only half the damage. Plot magnitudes and
    # mark the negative points, rather than silently dropping them the way the
    # original figure had to.
    def mag(v):
        return [abs(x) if x != 0 else float("nan") for x in v]

    negx = [yy for yy, v in zip(y, old) if v < 0]
    negy = [abs(v) for v in old if v < 0]

    fig, (ax, axe) = plt.subplots(
        2, 1, figsize=(9, 7.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.6], "hspace": 0.08})

    ax.plot(y, mag(old), lw=0.9, color="#1f77b4",
            label="old: Fortran, double  (in production today)")
    ax.plot(y, mag(new), lw=2.6, color="#2ca02c", alpha=0.75,
            label="new: C++ port, adaptive precision")
    ax.plot(y, mag(ref), lw=1.1, color="#d62728", ls="--",
            label="reference: same Fortran in REAL*16")
    if negx:
        ax.plot(negx, negy, "v", ms=4.5, color="#1f77b4", mfc="none", lw=0.6,
                label=f"old result was NEGATIVE ({len(negx)} points)")

    ax.set_yscale("log")
    ax.set_ylabel(r"|d$\sigma$|   [kbarn/MeV$^4$/(Z$\alpha$)$^4$]")
    ax.set_title("QED pair generator: differential cross section vs electron rapidity\n"
                 r"PbPb 5.36 TeV,  $p_T^\pm$ = 1 MeV,  $y_+$ = 0,  "
                 r"$\Delta\varphi = \pi$   (the ill-conditioned locus)",
                 fontsize=11)
    ax.legend(loc="lower center", fontsize=9, framealpha=0.95)
    ax.grid(alpha=0.25, which="both", lw=0.4)

    def relerr(v):
        return [abs(a - r) / abs(r) if r != 0 else float("nan")
                for a, r in zip(v, ref)]

    axe.plot(y, relerr(old), lw=0.9, color="#1f77b4")
    axe.plot(y, relerr(new), lw=1.6, color="#2ca02c")
    axe.axhline(1.0, color="0.4", lw=0.8, ls=":")
    axe.text(-9.8, 1.4, "100% error", fontsize=8, color="0.35")
    axe.set_yscale("log")
    axe.set_xlabel("electron rapidity  $y_-$")
    axe.set_ylabel("relative error\nvs reference")
    axe.grid(alpha=0.25, which="both", lw=0.4)

    worst_old = max(x for x in relerr(old) if x == x)
    worst_new = max(x for x in relerr(new) if x == x)
    axe.text(0.99, 0.90,
             f"worst error   old: {worst_old:.1e}    new: {worst_new:.1e}",
             transform=axe.transAxes, ha="right", va="top", fontsize=9)

    fig.savefig(png, dpi=140, bbox_inches="tight")
    print(f"wrote {png}")
    print(f"worst relative error -- old {worst_old:.3e}, new {worst_new:.3e}")
    print(f"points where the old result was negative: {len(negx)}")


if __name__ == "__main__":
    main()
