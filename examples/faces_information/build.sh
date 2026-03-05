#!/bin/bash
set -e

if [ -z "$MADVORO_PREFIX" ] || [ -z "$BOOST_DIR" ]; then
    echo "Usage: MADVORO_PREFIX=... BOOST_DIR=... bash build.sh"
    echo ""
    echo "  MADVORO_PREFIX - MadVoro installation prefix (containing lib/ and include/)"
    echo "  BOOST_DIR      - Boost installation prefix"
    exit 1
fi

g++ main.cpp \
    -I "${MADVORO_PREFIX}/include" \
    -I "${BOOST_DIR}/include" \
    -L "${MADVORO_PREFIX}/lib" -lmadvoro \
    -g -O2 \
    -lpthread -ldl -lz \
    -o test
