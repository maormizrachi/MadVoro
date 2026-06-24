#!/usr/bin/env bash
# install_deps.sh — Clone MadVoro's external dependencies for standalone builds.
#
# When MadVoro is used inside RICH, these packages are already present as
# submodules.  Standalone users should run this script once before building.
#
# Usage:
#   ./install_deps.sh [DEST_DIR]
#
# DEST_DIR defaults to ./deps

set -euo pipefail

DEST="${1:-deps}"
mkdir -p "$DEST"

clone_if_missing() {
    local name="$1" url="$2" branch="${3:-}"
    local target="$DEST/$name"
    if [ -d "$target" ]; then
        echo "[skip] $name already exists at $target"
    else
        echo "[clone] $name -> $target"
        if [ -n "$branch" ]; then
            git clone --depth 1 -b "$branch" "$url" "$target"
        else
            git clone --depth 1 "$url" "$target"
        fi
    fi
}

echo "=== Installing MadVoro dependencies into $DEST ==="
echo ""

# MPI utilities (serialization, exchange, collectives)
clone_if_missing "mpi_utils" "https://github.com/maormizrachi/mpi_utils.git"

# Spatial data structures (OctTree, KDTree, RangeTree, BoundingBox, Sphere)
clone_if_missing "spatial_ds" "https://github.com/maormizrachi/spatial_ds.git"

# Mesh decomposition (Hilbert ordering, load balancing, points manager,
# environment agents, indexing kernels)
clone_if_missing "MeshDecomposer3D" "https://github.com/maormizrachi/MeshDecomposer3D.git"

echo ""
echo "=== Done ==="
echo ""
echo "Pass the following to CMake:"
echo "  -DMADVORO_DEPS_DIR=$(cd "$DEST" && pwd)"
echo ""
echo "Or set the include paths manually:"
echo "  -I$DEST/mpi_utils -I$DEST/spatial_ds -I$DEST/MeshDecomposer3D"
