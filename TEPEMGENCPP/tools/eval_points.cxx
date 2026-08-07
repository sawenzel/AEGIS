// Evaluate Diffcross at a list of kinematic points produced by dump_points.
// Built twice: eval_points_d (double) and eval_points_q (-DTEP_QUAD, linked
// against diffcross.f compiled -freal-8-real-16). Diffing the two outputs is
// the precision measurement.
//
// Output columns:  ppvt yp pmvt ym dphi dsigma bad
#include "fortran_iface.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: %s E_GeV pointsfile\n", argv[0]);
    return 2;
  }
  tepreal energy = static_cast<tepreal>(std::atof(argv[1]));
  tepreal mass = static_cast<tepreal>(0.5109991);
  initdiffcross_(&energy, &mass);

  FILE* f = std::fopen(argv[2], "r");
  if (!f) { std::perror("fopen"); return 1; }

  std::printf("# ppvt yp pmvt ym dphi dsigma bad\n");
  char line[1024];
  while (std::fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    double a, b, c, d, e, w;
    if (std::sscanf(line, "%lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &w) < 5)
      continue;
    tepreal ta = a, tb = b, tc = c, td = d, te = e, out = 0;
    const int before = badpar_.badcount;
    diffcross_(&ta, &tb, &tc, &td, &te, &out);
    std::printf("%.17g %.17g %.17g %.17g %.17g %.17g %d\n", a, b, c, d, e,
                static_cast<double>(out), badpar_.badcount - before);
  }
  std::fclose(f);
  return 0;
}
