#!/bin/bash
# Configure and build AEGIS (including TEPEMGEN) for the C++ port work.
#
# Two flags below are needed by UPSTREAM AEGIS, which does not set them itself
# and consequently does not build on a current system:
#
#   -fallow-argument-mismatch
#       MICROCERN/sortzv.F is 1970s CERNLIB and passes REAL(4) where INTEGER(4)
#       is declared. gfortran >= 10 turned that from a warning into an error.
#
#   -DPYTHIA6_LIBRARY=...
#       GeneratorParam's find_package(Pythia6) ignores the PYTHIA6_ROOT
#       environment variable under CMake policy CMP0144, so the library has to
#       be named explicitly or configuration fails outright.
#
# Usage:  scripts/build.sh [extra cmake args...]
set -eo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC=$(cd "$HERE/../.." && pwd)            # the AEGIS source tree
BUILD=${BUILD_DIR:-$SRC/../build}

"$HERE/env.sh" bash -c "
  set -eo pipefail
  cmake -S '$SRC' -B '$BUILD' -Wno-dev \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE:-RelWithDebInfo} \
        -DCMAKE_Fortran_FLAGS=-fallow-argument-mismatch \
        -DPYTHIA6_LIBRARY=\$PYTHIA6_ROOT/lib/libpythia6.so $* &&
  cmake --build '$BUILD' -j ${JOBS:-8}
"

echo "built: $BUILD/TEPEMGEN/libTEPEMGEN.so"
