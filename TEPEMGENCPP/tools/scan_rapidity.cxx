// Scan the differential cross section against electron rapidity at fixed
// remaining kinematics -- the figure attached to O2-6340.
//
// Built twice from this one source: `scan_rapidity_d` (double) and
// `scan_rapidity_q` (-DTEP_QUAD, against diffcross.f compiled with
// -freal-8-real-16). Diffing the two outputs quantifies the significance lost
// to cancellation, and locating the corrupted window as a function of beam
// energy tests the arcosh(gamma) prediction in doc/02-findings.md section 2.
//
// Output columns:  ym  dsigma  badcount_delta
#include "fortran_iface.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  if (argc != 9) {
    std::fprintf(stderr,
                 "usage: %s E_GeV ptp_MeV ptm_MeV yp dphi ymin ymax nsteps\n",
                 argv[0]);
    return 2;
  }
  const double E = std::atof(argv[1]);
  const double ptp = std::atof(argv[2]);
  const double ptm = std::atof(argv[3]);
  const double yp = std::atof(argv[4]);
  const double dphi = std::atof(argv[5]);
  const double ymin = std::atof(argv[6]);
  const double ymax = std::atof(argv[7]);
  const int n = std::atoi(argv[8]);

  tepreal energy = static_cast<tepreal>(E);
  tepreal mass = static_cast<tepreal>(0.5109991);  // electron mass, as ee_init
  initdiffcross_(&energy, &mass);

  // gamma and wy=arcosh(gamma) come straight from the Fortran, so the
  // predicted danger zone is reported by the code under test rather than
  // recomputed here.
  std::fprintf(stderr, "# E=%g GeV  gamma=%.6g  wy=arcosh(gamma)=%.6g\n", E,
               static_cast<double>(physparam_.gamma),
               static_cast<double>(physparam_.wy));
  std::printf("# ym dsigma badcount_delta\n");
  std::printf("# E=%g gamma=%.10g wy=%.10g ptp=%g ptm=%g yp=%g dphi=%g\n", E,
              static_cast<double>(physparam_.gamma),
              static_cast<double>(physparam_.wy), ptp, ptm, yp, dphi);

  for (int i = 0; i <= n; ++i) {
    const double ym = ymin + (ymax - ymin) * i / n;
    tepreal a = static_cast<tepreal>(ptp);
    tepreal b = static_cast<tepreal>(yp);
    tepreal c = static_cast<tepreal>(ptm);
    tepreal d = static_cast<tepreal>(ym);
    tepreal e = static_cast<tepreal>(dphi);
    tepreal out = 0;
    const int before = badpar_.badcount;
    diffcross_(&a, &b, &c, &d, &e, &out);
    std::printf("%.10g %.17g %d\n", ym, static_cast<double>(out),
                badpar_.badcount - before);
  }
  return 0;
}
