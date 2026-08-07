// Sample kinematic points from the REAL generator and write them out, so the
// cross section at those points can afterwards be evaluated in double and in
// quad by eval_points_{d,q} and the two compared.
//
// Why not evaluate here: this binary links libTEPEMGEN, whose diffcross.f is
// the double build. The quad build has the same symbol names and a different
// ABI, so it cannot be linked into the same executable. Sampling and
// evaluation are therefore separate programs communicating through a file --
// which is also what makes the sample reusable as frozen golden data.
//
// Calls ee_init_/ee_event_ directly rather than going through TGenEpEmv1:
// the class exposes the per-event kinematics only via a protected method, and
// TEpEmGen::CalcXSection would abort() on us mid-sample.
//
// Output columns:  ppvt yp pmvt ym dphi weight
#include "TRandom3.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

extern "C" {
void ee_init_(double* ymin, double* ymax, double* ptmin, double* ptmax,
              double* cm_energy, double* z);
void ee_event_(double* ymin, double* ymax, double* ptmin, double* ptmax,
               double* ye, double* yp, double* xe, double* xp, double* phi,
               double* wtm2);
// This binary links the double build of diffcross.f, so it can record the
// double cross section alongside each point. compare_precision then needs only
// one pass and one intermediate file instead of two 200 MB ones.
void diffcross_(double* ppvt, double* yp, double* pmvt, double* ym,
                double* dphi, double* dsigma);
extern struct { int badcount; } badpar_;
}

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: %s nevents seed outfile\n", argv[0]);
    return 2;
  }
  const long n = std::atol(argv[1]);
  const int seed = std::atoi(argv[2]);
  const char* out = argv[3];

  gRandom = new TRandom3(seed);

  // The settings O2DPG actually uses: QEDGenParam.yMin=-7 yMax=7
  // ptMin=0.001 ptMax=1. (GeV -> MeV here), PbPb 5.36 TeV, Z=82.
  double ymin = -7., ymax = 7., ptmin = 1., ptmax = 1000.;
  double energy = 5360., z = 82.;
  ee_init_(&ymin, &ymax, &ptmin, &ptmax, &energy, &z);

  FILE* f = std::fopen(out, "w");
  if (!f) { std::perror("fopen"); return 1; }
  std::fprintf(f, "# ppvt yp pmvt ym dphi weight dsigma_double bad\n");
  std::fprintf(f, "# nevents=%ld seed=%d ymin=%g ymax=%g ptmin=%g ptmax=%g E=%g Z=%g\n",
               n, seed, ymin, ymax, ptmin, ptmax, energy, z);

  const int bad0 = badpar_.badcount;
  for (long i = 0; i < n; ++i) {
    double ye = 0, yp = 0, xe = 0, xp = 0, phi = 0, w = 0;
    ee_event_(&ymin, &ymax, &ptmin, &ptmax, &ye, &yp, &xe, &xp, &phi, &w);
    // ee_event returns log10(pt/MeV); Diffcross is called by it as
    // Diffcross(Ptp, Yp, Pte, Ye, Phi, ...) -- keep that argument order.
    double ptp = std::pow(10., xp);
    double pte = std::pow(10., xe);
    double ds = 0;
    const int b0 = badpar_.badcount;
    diffcross_(&ptp, &yp, &pte, &ye, &phi, &ds);
    std::fprintf(f, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %d\n", ptp, yp,
                 pte, ye, phi, w, ds, badpar_.badcount - b0);
  }
  std::fclose(f);
  std::fprintf(stderr, "wrote %ld points to %s; badcount over sample = %d\n", n,
               out, badpar_.badcount - bad0);
  return 0;
}
