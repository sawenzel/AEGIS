#!/bin/bash
# Same toolchain pinning as env.sh, but inside the O2sim environment, so the
# library is built against the ROOT that o2-sim will actually load (6.36.10,
# not the 6.36.04 that AEGIS::v1.5.9-9 pulls in).
#
# GCC-Toolchain 14.2 is the right compiler here and not a guess: the SHIPPED
# libTEPEMGEN.so and O2 itself both link that toolchain libstdc++
# (GLIBCXX_3.4.32). A system-g++ build would be the odd one out.
source /cvmfs/alice.cern.ch/etc/login.sh >/dev/null 2>&1
set -eo pipefail
GCCTC=/cvmfs/alice.cern.ch/el9-x86_64/Packages/GCC-Toolchain/v14.2.0-alice2-1
export CC="$GCCTC/bin/gcc" CXX="$GCCTC/bin/g++" FC="$GCCTC/bin/gfortran"
exec alienv setenv O2sim/v20260807-1 -c "$@"
