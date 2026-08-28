#!/bin/bash
set -euo pipefail

ROOT=/home/maorm/RICH
MADVORO_DIR="${ROOT}/source/3D/tessellation/voronoi"
BUILD_DIR="${ROOT}/build/gnuReleaseMPI/dependencies/madvoro"
OUT="${MADVORO_DIR}/build_periodic/example_periodic_boundaries_mpi"

DEFINES=$(grep CXX_DEFINES "${BUILD_DIR}/CMakeFiles/madvoro.dir/flags.make" | sed 's/CXX_DEFINES = //')

mpicxx -std=c++17 -O2 ${DEFINES} \
  -I"${MADVORO_DIR}" \
  -I"${MADVORO_DIR}/range/finders" \
  -I"${MADVORO_DIR}/examples" \
  -I"${ROOT}/source/3D/tessellation" \
  -I"${ROOT}/source/utils" \
  -I"${ROOT}/source/opt/vcl" \
  -I/software/x86_64/5.14.0/boost/1.78.0/include \
  -fopenmp \
  "${MADVORO_DIR}/examples/periodic_boundaries/main.cpp" \
  "${BUILD_DIR}/libmadvoro.a" \
  -o "${OUT}"

echo "Built ${OUT}"
