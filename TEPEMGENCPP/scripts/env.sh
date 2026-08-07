#!/bin/bash
# Build environment for the TEPEMGEN C++ port. x86_64 EL9 with CVMFS.
#
# Pins CC/CXX/FC explicitly to the CVMFS GCC-Toolchain. This is load-bearing,
# not tidiness: a stock EL9 has no gfortran at all, and /usr/bin precedes the
# CVMFS toolchain on PATH. So gfortran resolves to the CVMFS GCC 14.2 (nothing
# else provides it) while g++/gcc resolve to the system GCC 11.5 -- a silently
# mixed build across the one language boundary this project depends on, with
# the Fortran oracle and the C++ port compiled by different compilers.
# Pinning all three also keeps libquadmath and its callers in one toolchain.
#
# NOTE: no `set -u` around the CVMFS login script -- it references unset vars
# and dies under it, which manifests as this script producing no output at all.
#
# Usage:  scripts/env.sh <command...>
#   e.g.  scripts/env.sh bash -c 'echo $CXX'

source /cvmfs/alice.cern.ch/etc/login.sh >/dev/null 2>&1
set -eo pipefail

GCCTC=${GCCTC:-/cvmfs/alice.cern.ch/el9-x86_64/Packages/GCC-Toolchain/v14.2.0-alice2-1}
export CC="$GCCTC/bin/gcc" CXX="$GCCTC/bin/g++" FC="$GCCTC/bin/gfortran"

# AEGIS::<version> pulls in ROOT and pythia6, which is all the tree needs;
# MICROCERN is built from source inside AEGIS itself.
exec alienv setenv "${AEGIS_PACKAGE:-AEGIS::v1.5.9-9}" -c "$@"
