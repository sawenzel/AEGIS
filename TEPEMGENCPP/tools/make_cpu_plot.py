#!/usr/bin/env python3
"""CPU cost of every implementation, against what each one buys in accuracy.

Cost alone is misleading -- the cheapest option is also the broken one. So the
left panel is time per call and the right panel is the bias each choice leaves
on the cross section, on the same rows, so the trade is visible in one look.

Numbers come from tools/bench_fortran (Fortran, both precisions) and
tools/bench_cpp (every C++ variant), on the same 100k sampled points; bias
figures from tools/strategy.

Usage:  make_cpu_plot.py <out.png>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# label, us/call, |bias on sum(dsigma)| in %, colour, note
ROWS = [
    ("Fortran double\n(production today)", 10.01, 7.31e-2, "#c44e52", "broken"),
    ("C++ double", 11.02, 7.31e-2, "#c44e52", "broken"),
    ("C++ long double", 25.94, 1.0e-5, "#dd8452", ""),
    ("C++ adaptive\n(shipped)", 531.76, 1.0e-5, "#55a868", "correct"),
    ("Fortran REAL*16", 978.89, 0.0, "#4c72b0", "reference"),
    ("C++ __float128", 1035.81, 0.0, "#4c72b0", "reference"),
]


def main():
    out = sys.argv[1]
    labels = [r[0] for r in ROWS]
    us = [r[1] for r in ROWS]
    bias = [r[2] for r in ROWS]
    cols = [r[3] for r in ROWS]
    notes = [r[4] for r in ROWS]
    ypos = list(range(len(ROWS)))[::-1]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True,
                                 gridspec_kw={"width_ratios": [1.5, 1],
                                              "wspace": 0.06})

    a1.barh(ypos, us, color=cols, height=0.62)
    for y, v, n in zip(ypos, us, notes):
        a1.text(v * 1.15, y, f"{v:.0f} us   ({v/10.01:.0f}x)" + (f"   {n}" if n else ""),
                va="center", fontsize=9)
    a1.set_xscale("log")
    a1.set_xlim(5, 12000)
    a1.set_yticks(ypos)
    a1.set_yticklabels(labels, fontsize=9.5)
    a1.set_xlabel("CPU time per Diffcross call  [us]   (log scale)")
    a1.set_title("Cost", fontsize=11)
    a1.grid(axis="x", alpha=0.25, which="both", lw=0.4)

    # A zero bias cannot be drawn on a log axis; the reference rows are exact
    # by definition, so they are annotated rather than plotted as bars.
    floor = 1e-6
    a2.barh(ypos, [b if b > 0 else floor for b in bias], color=cols, height=0.62)
    for y, b in zip(ypos, bias):
        a2.text(max(b, floor) * 1.4, y,
                "exact (reference)" if b == 0 else f"{b:.0e} %",
                va="center", fontsize=9)
    a2.set_xscale("log")
    a2.set_xlim(floor, 3)
    a2.set_xlabel("| bias on the total cross section |  [%]")
    a2.set_title("What it buys", fontsize=11)
    a2.grid(axis="x", alpha=0.25, which="both", lw=0.4)

    fig.suptitle("QED pair generator: cost vs accuracy of each implementation\n"
                 "100k sampled points, PbPb 5.36 TeV, production cuts",
                 fontsize=12, y=1.04)
    fig.text(0.5, -0.10,
             "The two cheapest rows are also the wrong ones: they bias the cross "
             "section by 0.07% on average and, at 3M events,\n"
             "produce single outliers ~1.6 million times too large that inflate "
             "the variance estimate by 800% and abort the run.",
             ha="center", fontsize=9, color="0.3")

    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
