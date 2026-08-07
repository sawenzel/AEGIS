// Time the FORTRAN diffcross.f, in whichever precision it was compiled with.
//
// Built twice, against diffcross_d.o and diffcross_q.o, so the REAL*16 cost is
// measured rather than assumed from the C++ __float128 number. They are not
// the same thing: gfortran's REAL*16 and g++'s __float128 both end up in
// libquadmath, but the surrounding code, the argument passing and the
// intrinsics are compiled differently.
//
// Prints one line, machine-readable, for tools/make_cpu_plot.py.
#include "fortran_iface.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s label pointsfile [n]\n", argv[0]);
    return 2;
  }
  const char* label = argv[1];
  const long limit = (argc > 3) ? std::atol(argv[3]) : 100000;

  tepreal energy = 5360;
  tepreal mass = tepreal(5109991) / tepreal(10000000);
  initdiffcross_(&energy, &mass);

  struct P { double a, b, c, d, e; };
  std::vector<P> pts;
  FILE* f = std::fopen(argv[2], "r");
  if (!f) { std::perror("fopen"); return 1; }
  char line[1024];
  while (std::fgets(line, sizeof line, f) && (long)pts.size() < limit) {
    if (line[0] == '#') continue;
    P p;
    if (std::sscanf(line, "%lf %lf %lf %lf %lf", &p.a, &p.b, &p.c, &p.d, &p.e) == 5)
      pts.push_back(p);
  }
  std::fclose(f);
  if (pts.empty()) { std::fprintf(stderr, "no points\n"); return 1; }

  volatile double sink = 0;
  const auto t0 = std::chrono::steady_clock::now();
  for (const auto& p : pts) {
    tepreal a = p.a, b = p.b, c = p.c, d = p.d, e = p.e, out = 0;
    diffcross_(&a, &b, &c, &d, &e, &out);
    sink = sink + static_cast<double>(out);
  }
  const double us =
      std::chrono::duration<double, std::micro>(
          std::chrono::steady_clock::now() - t0).count() / pts.size();

  std::printf("%s %.4f\n", label, us);
  return 0;
}
