#!/bin/bash
set -e

if [ -z "$MADVORO_PREFIX" ] || [ -z "$BOOST_DIR" ]; then
    echo "Usage: MADVORO_PREFIX=... BOOST_DIR=... [VTK_DIR=...] [HDF5_DIR=...] bash build.sh"
    echo ""
    echo "  MADVORO_PREFIX - MadVoro installation prefix (containing lib/ and include/)"
    echo "  BOOST_DIR      - Boost installation prefix"
    echo "  VTK_DIR        - (optional) VTK installation prefix, enables ToVTK output"
    echo "  HDF5_DIR       - (optional) HDF5 installation prefix (with C++ support), enables ToHDF5 output"
    exit 1
fi

EXTRA_FLAGS=""
EXTRA_LIBS=""

if [ -n "$HDF5_DIR" ]; then
    EXTRA_FLAGS="$EXTRA_FLAGS -DMADVORO_WITH_HDF5"
    EXTRA_LIBS="$EXTRA_LIBS -L${HDF5_DIR}/lib -lhdf5_cpp -lhdf5_hl_cpp -lhdf5 -lhdf5_hl"
fi

if [ -n "$VTK_DIR" ]; then
    EXTRA_FLAGS="$EXTRA_FLAGS -DMADVORO_WITH_VTK"
    VTK_LIBDIR="${VTK_DIR}/lib64"
    if [ ! -d "$VTK_LIBDIR" ]; then
        VTK_LIBDIR="${VTK_DIR}/lib"
    fi
    EXTRA_LIBS="$EXTRA_LIBS -L${VTK_LIBDIR}"
fi

mpicxx main.cpp \
    -I "${MADVORO_PREFIX}/include" \
    -I "${BOOST_DIR}/include" \
    $EXTRA_FLAGS \
    -L "${MADVORO_PREFIX}/lib" -lmadvoro \
    $EXTRA_LIBS \
    -g -O2 \
    -lpthread -ldl -lz \
    -o test
