#ifndef FOX_MPI_UTILS_HPP
#define FOX_MPI_UTILS_HPP

#include <vector>
#include <mpi.h>
#include <madvoro/Voronoi3D.hpp>

using namespace MadVoro;

std::pair<std::vector<Vector3D>, std::vector<double>> SpreadPointsToProcessors(const std::vector<Vector3D> &points, const std::vector<double> &isInside)
{
    struct PointData
    {
        double x, y, z;
        double isInside;
    };

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t numPoints = points.size();
    
    // broadcast points number to all
    MPI_Bcast(&numPoints, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    size_t pointsPerRank = numPoints / size;
    size_t remainingPoints = numPoints % size;
    size_t myPointsNum = pointsPerRank + ((rank == size - 1)? remainingPoints : 0);
    std::vector<PointData> recvBuff(myPointsNum);
    
    if(rank == 0)
    {
        std::vector<int> offsets(size, 0), counts(size, 0);
        for(int i = 0; i < size; i ++)
        {
            counts[i] = pointsPerRank + ((i == size - 1)? remainingPoints : 0); 
            counts[i] *= sizeof(PointData);
            if(i > 0)
            {
                offsets[i] = offsets[i - 1] + counts[i - 1];
            }
        }

        std::vector<PointData> sendBuff;
        for(size_t i = 0; i < numPoints; i ++)
        {
            sendBuff.push_back({points[i].x, points[i].y, points[i].z, isInside[i]});
        }

        MPI_Scatterv(sendBuff.data(), counts.data(), offsets.data(), MPI_BYTE, recvBuff.data(), myPointsNum * sizeof(PointData), MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Scatterv(NULL, NULL, NULL, MPI_BYTE, recvBuff.data(), myPointsNum * sizeof(PointData), MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    
    std::vector<Vector3D> myPoints;
    std::vector<double> myIsInside;

    for(size_t i = 0; i < myPointsNum; i ++)
    {
        myPoints.push_back(Vector3D(recvBuff[i].x, recvBuff[i].y, recvBuff[i].z));
        myIsInside.push_back(recvBuff[i].isInside);
    }
    return {myPoints, myIsInside};
}

std::pair<std::vector<Vector3D>, std::vector<double>> GetPointsAfterBuildExchange(const Voronoi3D &voronoi, const std::vector<Vector3D> &originalPoints, const std::vector<double> &isInside)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t N = voronoi.GetPointNo(); // my points num after exchange
    std::vector<Vector3D> myNewPoints = voronoi.getMeshPoints();
    myNewPoints.resize(N);
    std::vector<double> myNewIsInside;

    // first points in the list of points are always points that were mine, and stayed at mine.
    size_t pointsSendToSelf = voronoi.GetSelfIndex().size();
    for(size_t i = 0; i < pointsSendToSelf; i++)
    {
        size_t pointIdxBeforeChange = voronoi.GetSelfIndex()[i];
        myNewIsInside.push_back(isInside[pointIdxBeforeChange]);
    }

    std::vector<double> toSend;
    std::vector<int> sendCounts(size, 0), sendDispls(size, 0);

    // next are points from other processors. The processors I communicated with are saved in `voronoi.GetSentProcs()`
    // The points I sent to rank `voronoi.GetSentProcs()[i]`, are saved in `voronoi.GetSentPoints()[i]`.

    const std::vector<int> &communicatedRanks = voronoi.GetSentProcs();

    // the send buffer is eventially sent to 'Alltoallv', so we need to arrange it by ranks
    for(int otherRank = 0; otherRank < size; otherRank++)
    {
        size_t i = std::distance(communicatedRanks.begin(), std::find(communicatedRanks.begin(), communicatedRanks.end(), otherRank));
        if(i == communicatedRanks.size())
        {
            continue; // nothing sent to rank
        }

        // the points I received from the rank are in
        const std::vector<size_t> &pointsSentToRank = voronoi.GetSentPoints()[i]; 
        size_t countForRank = pointsSentToRank.size();
        sendCounts[otherRank] = countForRank;
        for(size_t j = 0; j < countForRank; j++)
        {
            size_t oldIndexInMine = pointsSentToRank[j];
            toSend.push_back(isInside[oldIndexInMine]);
        }
    }
    
    // synchronize send and recv counts
    std::vector<int> recvCounts(size);
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> recvDispls(size, 0);
    size_t recvTotal = recvCounts[0];
    
    // calculate displacements
    for(int i = 1; i < size; i++)
    {
        sendDispls[i] = sendDispls[i - 1] + sendCounts[i - 1];
        recvDispls[i] = recvDispls[i - 1] + recvCounts[i - 1];
        recvTotal += recvCounts[i];
    }
    
    // prepare joint recv buffer
    std::vector<double> toRecv(recvTotal);
    MPI_Alltoallv(toSend.data(), sendCounts.data(), sendDispls.data(), MPI_DOUBLE, toRecv.data(), recvCounts.data(), recvDispls.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    
    // the list `voronoi.GetSentProcs()` is not only the list of processors I communicated with, but also the ranks whom I received from.
    for(size_t i = 0; i < communicatedRanks.size(); i++)
    {
        int rankCommunicated = communicatedRanks[i];
        size_t countForRank = recvCounts[rankCommunicated];
        int displacement = recvDispls[rankCommunicated];
        for(size_t j = 0; j < countForRank; j++)
        {
            assert(displacement + j < toRecv.size());
            myNewIsInside.push_back(toRecv[displacement + j]);
        }
    }
    
    return {myNewPoints, myNewIsInside};
}

#endif // FOX_MPI_UTILS_HPP