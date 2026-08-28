#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <Voronoi3D.hpp>
#include "../Vector3D.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#include "range/SentPointsContainer.hpp"
#endif

static constexpr double REL_TOL = 1e-6;

static bool ApproxEqual(double a, double b, double relTol = REL_TOL)
{
    double scale = std::max(std::abs(a), std::abs(b));
    if(scale < 1.0)
    {
        scale = 1.0;
    }
    return std::abs(a - b) <= relTol * scale;
}

#ifdef MADVORO_WITH_MPI
static bool TestSentPointsContainer()
{
    SentPointsContainer container;
    std::vector<size_t> first = container.addPointsAsSent(1, std::vector<size_t>({42}), ZERO_IMAGE_CODE);
    std::vector<size_t> second = container.addPointsAsSent(1, std::vector<size_t>({42}), ZERO_IMAGE_CODE);
    Vector3D txMinus(-1.0, 0.0, 0.0);
    Vector3D tyMinus(0.0, -1.0, 0.0);
    std::vector<size_t> third = container.addPointsAsSent(1, std::vector<size_t>({42}), ImageCode(txMinus));
    std::vector<size_t> fourth = container.addPointsAsSent(1, std::vector<size_t>({42}), ImageCode(tyMinus));

    if(first.size() != 1 || first[0] != 42)
    {
        std::cerr << "SentPointsContainer: expected first insert to return {42}" << std::endl;
        return false;
    }
    if(!second.empty())
    {
        std::cerr << "SentPointsContainer: expected duplicate zero-code insert to return empty" << std::endl;
        return false;
    }
    if(third.size() != 1 || third[0] != 42)
    {
        std::cerr << "SentPointsContainer: expected x-minus image insert to return {42}" << std::endl;
        return false;
    }
    if(fourth.size() != 1 || fourth[0] != 42)
    {
        std::cerr << "SentPointsContainer: expected y-minus image insert to return {42}" << std::endl;
        return false;
    }
    return true;
}
#endif

static bool CellHasPhysicalBoundary(const MadVoro::Voronoi3D<Vector3D> &voronoi, size_t cellIndex)
{
    const MadVoro::face_vec &faces = voronoi.GetCellFaces(cellIndex);
    for(size_t faceIdx : faces)
    {
        if(voronoi.BoundaryFace(faceIdx))
        {
            return true;
        }
    }
    return false;
}

static bool TestOnePointPeriodicSerial()
{
    Vector3D ll(0.0, 0.0, 0.0);
    Vector3D ur(1.0, 1.0, 1.0);
    MadVoro::Voronoi3D<Vector3D> voronoi(ll, ur);
    voronoi.SetPeriodicBoundaries(true, true, true);

    std::vector<Vector3D> points = {Vector3D(0.5, 0.5, 0.5)};
    voronoi.Build(points);

    double boxVolume = 1.0;
    double cellVolume = voronoi.GetVolume(0);
    if(!ApproxEqual(cellVolume, boxVolume))
    {
        std::cerr << "One-point periodic: expected cell volume " << boxVolume << ", got " << cellVolume << std::endl;
        return false;
    }
    if(CellHasPhysicalBoundary(voronoi, 0))
    {
        std::cerr << "One-point periodic: cell has physical boundary faces" << std::endl;
        return false;
    }
    return true;
}

static bool TestLattice2x2x2Serial()
{
    Vector3D ll(0.0, 0.0, 0.0);
    Vector3D ur(1.0, 1.0, 1.0);
    MadVoro::Voronoi3D<Vector3D> voronoi(ll, ur);
    voronoi.SetPeriodicBoundaries(true, true, true);

    std::vector<Vector3D> points;
    for(int ix = 0; ix < 2; ++ix)
    {
        for(int iy = 0; iy < 2; ++iy)
        {
            for(int iz = 0; iz < 2; ++iz)
            {
                points.emplace_back(0.25 + 0.5 * ix, 0.25 + 0.5 * iy, 0.25 + 0.5 * iz);
            }
        }
    }
    voronoi.Build(points);

    double expectedCellVolume = 1.0 / 8.0;
    for(size_t i = 0; i < points.size(); ++i)
    {
        if(!ApproxEqual(voronoi.GetVolume(i), expectedCellVolume))
        {
            std::cerr << "2x2x2 lattice: cell " << i << " volume " << voronoi.GetVolume(i)
                      << " != " << expectedCellVolume << std::endl;
            return false;
        }
    }
    return true;
}

#ifdef MADVORO_WITH_MPI
static bool TestOwnerWrapping(const MadVoro::Voronoi3D<Vector3D> &voronoi)
{
    Vector3D wrapped = voronoi.WrapPeriodicPoint(Vector3D(-0.03, 0.5, 0.5));
    if(!ApproxEqual(wrapped.x, 0.97))
    {
        std::cerr << "WrapPeriodicPoint: expected x=0.97, got " << wrapped.x << std::endl;
        return false;
    }
    int ownerWrapped = voronoi.GetOwner(wrapped);
    int ownerShifted = voronoi.GetOwner(Vector3D(4.97, 0.5, 0.5));
    if(ownerWrapped != ownerShifted)
    {
        std::cerr << "GetOwner: owner mismatch for equivalent periodic coordinates" << std::endl;
        return false;
    }
    return true;
}

static bool TestPeriodicParallel(MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    Vector3D ll(0.0, 0.0, 0.0);
    Vector3D ur(1.0, 1.0, 1.0);
    MadVoro::Voronoi3D<Vector3D> voronoi(ll, ur);
    voronoi.SetPeriodicBoundaries(true, true, true);

    std::vector<Vector3D> localPoints;
    for(int ix = 0; ix < 2; ++ix)
    {
        for(int iy = 0; iy < 2; ++iy)
        {
            for(int iz = 0; iz < 2; ++iz)
            {
                size_t globalIdx = static_cast<size_t>(ix * 4 + iy * 2 + iz);
                if(globalIdx % static_cast<size_t>(size) == static_cast<size_t>(rank))
                {
                    localPoints.emplace_back(0.25 + 0.5 * ix, 0.25 + 0.5 * iy, 0.25 + 0.5 * iz);
                }
            }
        }
    }

    voronoi.BuildParallel(localPoints);

    double localVolume = 0.0;
    for(size_t i = 0; i < voronoi.GetPointNo(); ++i)
    {
        localVolume += voronoi.GetVolume(i);
    }
    double totalVolume = 0.0;
    MPI_Allreduce(&localVolume, &totalVolume, 1, MPI_DOUBLE, MPI_SUM, comm);

    bool pass = ApproxEqual(totalVolume, 1.0);
    if(!pass && rank == 0)
    {
        std::cerr << "Periodic parallel: total volume " << totalVolume << " != 1.0 (size=" << size << ")" << std::endl;
    }

    bool ownerPass = TestOwnerWrapping(voronoi);
    if(!ownerPass && rank == 0)
    {
        std::cerr << "Periodic parallel: owner wrapping test failed" << std::endl;
    }

    int localOk = (pass && ownerPass) ? 1 : 0;
    int globalOk = 0;
    MPI_Allreduce(&localOk, &globalOk, 1, MPI_INT, MPI_MIN, comm);
    return globalOk == 1;
}
#endif

int main(int argc, char *argv[])
{
    int failed = 0;

#ifdef MADVORO_WITH_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(!TestSentPointsContainer())
    {
        if(rank == 0)
        {
            std::cerr << "FAILED: SentPointsContainer image-code deduplication" << std::endl;
        }
        ++failed;
    }
    else if(rank == 0)
    {
        std::cout << "PASSED: SentPointsContainer image-code deduplication" << std::endl;
    }
#endif

#ifndef MADVORO_WITH_MPI
    if(!TestOnePointPeriodicSerial())
    {
        std::cerr << "FAILED: one-point fully periodic serial test" << std::endl;
        ++failed;
    }
    else
    {
        std::cout << "PASSED: one-point fully periodic serial test" << std::endl;
    }

    if(!TestLattice2x2x2Serial())
    {
        std::cerr << "FAILED: 2x2x2 periodic lattice serial test" << std::endl;
        ++failed;
    }
    else
    {
        std::cout << "PASSED: 2x2x2 periodic lattice serial test" << std::endl;
    }
#endif

#ifdef MADVORO_WITH_MPI
    if(!TestPeriodicParallel(MPI_COMM_WORLD))
    {
        if(rank == 0)
        {
            std::cerr << "FAILED: periodic parallel test" << std::endl;
        }
        ++failed;
    }
    else if(rank == 0)
    {
        std::cout << "PASSED: periodic parallel test" << std::endl;
    }
    MPI_Finalize();
#endif

    if(failed > 0)
    {
        return EXIT_FAILURE;
    }
    std::cout << "All periodic boundary tests passed." << std::endl;
    return EXIT_SUCCESS;
}
