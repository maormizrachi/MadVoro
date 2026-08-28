#!/bin/bash
#SBATCH --job-name=madvoro-periodic
#SBATCH --partition=bigrun
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --output=/home/maorm/RICH/source/3D/tessellation/voronoi/examples/periodic_boundaries/periodic_boundaries-%j.out
#SBATCH --error=/home/maorm/RICH/source/3D/tessellation/voronoi/examples/periodic_boundaries/periodic_boundaries-%j.err

set -euo pipefail

ROOT=/home/maorm/RICH
MADVORO_DIR="${ROOT}/source/3D/tessellation/voronoi"
"${MADVORO_DIR}/examples/periodic_boundaries/build_mpi_test.sh"

TEST_BIN="${MADVORO_DIR}/build_periodic/example_periodic_boundaries_mpi"

echo "Running periodic boundary tests with ${SLURM_NTASKS} MPI ranks"
export OMP_NUM_THREADS=1
mpirun -np "${SLURM_NTASKS}" "${TEST_BIN}"
