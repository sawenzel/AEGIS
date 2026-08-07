#!/bin/bash
# Build the Phase 0 oracle tools.
#
# diffcross.f is compiled TWICE from the same untouched source:
#   *_d  ordinary DOUBLE PRECISION
#   *_q  -freal-8-real-16, i.e. every DOUBLE PRECISION becomes REAL*16
# The two share symbol names but not ABI, so they go into separate binaries.
# This yields a quad reference before any C++ port exists -- see
# tools/fortran_iface.h.
#
# dump_points is different: it links the real libTEPEMGEN (double) because it
# needs the sampler ee_event_ and its eernd -> gRandom bridge.
#
# Usage:  tools/build_tools.sh
set -eo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TEP=$(cd "$HERE/../../TEPEMGEN" && pwd)
LIB=${LIB_DIR:-$(cd "$HERE/../../.." && pwd)/build/TEPEMGEN}
OUT=${OUT_DIR:-$(cd "$HERE/../../.." && pwd)/build/tepemgencpp-tools}
mkdir -p "$OUT"

"$HERE/../scripts/env.sh" bash -c "
set -eo pipefail
cd '$OUT'
FFLAGS='-O2 -fallow-argument-mismatch -std=legacy'

\$FC \$FFLAGS                    -c '$TEP/diffcross.f' -o diffcross_d.o
\$FC \$FFLAGS -freal-8-real-16   -c '$TEP/diffcross.f' -o diffcross_q.o

# Standalone tools: their own copy of diffcross, no ROOT, no libTEPEMGEN.
for prog in scan_rapidity eval_points; do
  \$CXX -std=c++17 -O2 -I'$HERE'            -c '$HERE'/\$prog.cxx -o \${prog}_d.o
  \$CXX -std=c++17 -O2 -I'$HERE' -DTEP_QUAD -c '$HERE'/\$prog.cxx -o \${prog}_q.o
  \$CXX \${prog}_d.o diffcross_d.o -o \${prog}_d -lgfortran -lquadmath
  \$CXX \${prog}_q.o diffcross_q.o -o \${prog}_q -lgfortran -lquadmath
done

# Quad only: reads the double result out of the points file and compares.
\$CXX -std=c++17 -O2 -I'$HERE' -DTEP_QUAD -c '$HERE'/compare_precision.cxx -o compare_precision.o
\$CXX compare_precision.o diffcross_q.o -o compare_precision -lgfortran -lquadmath

# Phase 1 gate: C++ port against the Fortran, at matching precision.
\$CXX -std=c++17 -O2 -I'$HERE' -I'$HERE/../include'            -c '$HERE'/validate_cpp.cxx -o validate_cpp_d.o
\$CXX -std=c++17 -O2 -I'$HERE' -I'$HERE/../include' -DTEP_QUAD -c '$HERE'/validate_cpp.cxx -o validate_cpp_q.o
\$CXX validate_cpp_d.o diffcross_d.o -o validate_cpp_d -lgfortran -lquadmath
\$CXX validate_cpp_q.o diffcross_q.o -o validate_cpp_q -lgfortran -lquadmath

# Calibration: pure C++, both instantiations in one binary (no Fortran, so no
# symbol clash). The quad instantiation is the reference once Phase 1 passed.
\$CXX -std=c++17 -O2 -I'$HERE' -I'$HERE/../include' '$HERE'/calibrate.cxx -o calibrate -lquadmath

\$CXX -std=c++17 -O2 -I'$HERE' -I'$HERE/../include' '$HERE'/bench.cxx -o bench -lquadmath

\$CXX -std=c++17 -O2 -I'$HERE' -I'$HERE/../include' '$HERE'/strategy.cxx -o strategy -lquadmath

\$CXX -std=c++17 -O2 -I'$HERE' -I'$HERE/../include' '$HERE'/scan_cpp.cxx -o scan_cpp -lquadmath

\$CXX -std=c++17 -O2 -I'$HERE' -I'$HERE/../include' '$HERE'/test_norm.cxx -o test_norm -lquadmath

# Sampler: needs the real library (ee_event_ plus the eernd -> gRandom bridge).
\$CXX -std=c++17 -O2 '$HERE'/dump_points.cxx -o dump_points \
     \$(root-config --cflags --libs) -L'$LIB' -lTEPEMGEN -Wl,-rpath,'$LIB' -lgfortran
"

echo "tools in: $OUT"
ls -1 "$OUT" | grep -vE '\.o$' || true
