// Where a rank's build is written: the mesh's own arrays (by reference), and what only the build gives
#ifndef MADVORO_BUILDER_DISTRIBUTED_BUILD_HPP
#define MADVORO_BUILDER_DISTRIBUTED_BUILD_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "types.hpp"

namespace MadVoro
{

//! \brief The arrays the build fills belong to the mesh: they must outlive the build call, and the builder keeps no reference to them after it
template <typename PointT>
struct DistributedBuild
{
    // the cells' points, 4 unused points, then the ghosts
    std::vector<PointT> &points;
    std::vector<face_vec> &cellFaces;
    std::vector<std::pair<std::size_t, std::size_t>> &faceNeighbors;
    std::vector<point_vec> &faceVertices;
    std::vector<double> &faceArea;
    std::vector<PointT> &faceCM;
    std::vector<PointT> &vertices;
    std::vector<double> &volume;
    // the cells' centroids, and of a mirror ghost the mirror of its cell's
    std::vector<PointT> &CM;
    #ifdef MADVORO_WITH_MPI
        // written by a build that searches other ranks: the ranks it shares points with, the rank's points each has (indices among the rank's points), and where each one's points are among the ghosts
        std::vector<int> &duplicatedProcs;
        std::vector<std::vector<std::size_t>> &duplicatedPoints;
        std::vector<std::vector<std::size_t>> &ghostIndices;
    #endif // MADVORO_WITH_MPI

    // a ghost that copies one of the rank's points (not built, or a periodic image) -> that point's index among the rank's points
    AllPointsMap ghostRankIndex = AllPointsMap();
    // the index of the first of the ghosts the build started with
    std::size_t carriedGhostsIndex = 0;
    // the index of each point other ranks sent, in the order they came
    std::vector<std::size_t> receivedIndex = std::vector<std::size_t>();
    std::size_t rounds = 0;
};

} // namespace MadVoro

#endif // MADVORO_BUILDER_DISTRIBUTED_BUILD_HPP
