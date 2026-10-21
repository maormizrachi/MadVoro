/** \file MergeSubdomains.hpp
   \brief Merges the tessellations of a rank's thread subdomains into one
   \author Maor Mizrachi
*/
#ifndef MADVORO_SUBDOMAINS_MERGE_SUBDOMAINS_HPP
#define MADVORO_SUBDOMAINS_MERGE_SUBDOMAINS_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../types.hpp"

namespace MadVoro
{

template <typename PointT>
class Voronoi3D;

/*! \brief The tessellations of a Voronoi3D's thread subdomains, kept between builds
    to reuse their memory. A copy is empty: the copied tessellation makes its own
    subdomains on its next build. */
template <typename PointT>
struct SubdomainMeshes
{
    SubdomainMeshes() = default;

    SubdomainMeshes(const SubdomainMeshes &) {}

    SubdomainMeshes &operator=(const SubdomainMeshes &)
    {
        this->meshes.clear();
        return *this;
    }

    std::vector<std::unique_ptr<Voronoi3D<PointT>>> meshes;
};

//! \brief A tessellation merged from the tessellations of thread subdomains.
template <typename PointT>
struct SubdomainMerge
{
    //! The rank's points, the four bounding points, the points other ranks sent, then the mirror ghosts of every subdomain
    std::vector<PointT> points;
    //! The faces of each of the rank's points
    std::vector<face_vec> cellFaces;
    //! The two points of each face
    std::vector<std::pair<std::size_t, std::size_t>> faceNeighbors;
    //! The vertices of each face, in `vertices`
    std::vector<point_vec> faceVertices;
    std::vector<double> faceArea;
    std::vector<PointT> faceCM;
    //! Face vertices (the circumcenters of every subdomain's tetrahedra)
    std::vector<PointT> vertices;
    //! The volume of each of the rank's points
    std::vector<double> volume;
    //! The centroid of each point in `points` (of a mirror ghost: the mirror of its cell's centroid)
    std::vector<PointT> CM;
};

/*! \brief Merges the tessellations of thread subdomains.
    Each subdomain is a partial build of the rank's points followed by the points other
    ranks sent (the all-points indices of these follow the rank's points): its own
    points are some of the rank's points, and its ghosts are other points of the rank
    (the points of other subdomains), points other ranks sent, or mirror ghosts of its
    own points. A face between cells of two subdomains is built by both; the copy of the
    first subdomain that built it is kept, and it is listed by both cells. Each cell's
    volume and centroid are those of the subdomain that owns it.
    \param subdomains The subdomains, in order
    \param points The rank's points (`pointsNum` of them; together, the subdomains own
           each of them once), four places for the bounding points, then the points
           other ranks sent */
template <typename PointT>
SubdomainMerge<PointT> MergeSubdomains(const std::vector<const Voronoi3D<PointT> *> &subdomains, const std::vector<PointT> &points,
                                       std::size_t pointsNum)
{
    constexpr std::size_t INVALID = std::numeric_limits<std::size_t>::max();
    constexpr std::size_t BOUNDING_POINTS = 4;

    SubdomainMerge<PointT> merged;
    merged.points = points;
    if(not subdomains.empty())
    {
        const Voronoi3D<PointT> &first = *subdomains.front();
        const std::vector<PointT> &firstPoints = first.getMeshPoints();
        std::copy(firstPoints.begin() + first.GetPointNo(), firstPoints.begin() + first.GetPointNo() + BOUNDING_POINTS,
                  merged.points.begin() + pointsNum);
    }
    merged.CM.resize(merged.points.size());
    merged.cellFaces.resize(pointsNum);
    merged.volume.resize(pointsNum);
    // a point's index among all the points of a subdomain, in the merged tessellation
    const auto fromAllPoints = [pointsNum](std::size_t index) { return (index < pointsNum) ? index : index + BOUNDING_POINTS; };

    // Faces between cells of two subdomains, by their (smaller, larger) cell pair
    std::unordered_map<std::uint64_t, std::size_t> sharedFaces;
    const auto sharedFaceKey = [pointsNum](std::size_t a, std::size_t b)
    {
        return static_cast<std::uint64_t>(std::min(a, b)) * pointsNum + std::max(a, b);
    };

    for(const Voronoi3D<PointT> *subdomainPtr : subdomains)
    {
        const Voronoi3D<PointT> &subdomain = *subdomainPtr;
        const std::size_t ownNum = subdomain.GetPointNo();
        const std::vector<PointT> &subPoints = subdomain.getMeshPoints();
        const AllPointsMap &indicesInPoints = subdomain.GetIndicesInAllPoints();

        // The subdomain's points in the merged tessellation
        std::vector<std::size_t> toMerged(subPoints.size(), INVALID);
        for(std::size_t j = 0; j < ownNum; ++j)
        {
            toMerged[j] = fromAllPoints(indicesInPoints.at(j));
        }
        for(std::size_t j = ownNum + BOUNDING_POINTS; j < subPoints.size(); ++j)
        {
            const auto it = indicesInPoints.find(j);
            if(it != indicesInPoints.cend())
            {
                // a point of another subdomain, or a point another rank sent
                toMerged[j] = fromAllPoints(it->second);
            }
            else
            {
                // a mirror ghost of one of the subdomain's points
                toMerged[j] = merged.points.size();
                merged.points.push_back(subPoints[j]);
                merged.CM.push_back(subdomain.GetCellCM(j));
            }
        }

        // Faces, with their vertices moved past the vertices of the previous subdomains
        const std::size_t vertexOffset = merged.vertices.size();
        const std::vector<std::pair<std::size_t, std::size_t>> &faceNeighbors = subdomain.GetAllFaceNeighbors();
        std::vector<std::size_t> faceToMerged(faceNeighbors.size());
        for(std::size_t f = 0; f < faceNeighbors.size(); ++f)
        {
            const std::size_t a = toMerged[faceNeighbors[f].first];
            const std::size_t b = toMerged[faceNeighbors[f].second];
            const bool shared = (a < pointsNum) and (b < pointsNum) and
                                ((faceNeighbors[f].first >= ownNum) or (faceNeighbors[f].second >= ownNum));
            if(shared)
            {
                const auto [it, inserted] = sharedFaces.emplace(sharedFaceKey(a, b), merged.faceNeighbors.size());
                if(not inserted)
                {
                    // already built by a previous subdomain
                    faceToMerged[f] = it->second;
                    continue;
                }
            }
            faceToMerged[f] = merged.faceNeighbors.size();
            merged.faceNeighbors.emplace_back(a, b);
            point_vec faceVertices = subdomain.GetPointsInFace(f);
            for(std::size_t &vertex : faceVertices)
            {
                vertex += vertexOffset;
            }
            merged.faceVertices.push_back(std::move(faceVertices));
            merged.faceArea.push_back(subdomain.GetArea(f));
            merged.faceCM.push_back(subdomain.FaceCM(f));
        }
        const std::vector<PointT> &vertices = subdomain.GetFacePoints();
        merged.vertices.insert(merged.vertices.end(), vertices.begin(), vertices.end());

        // Cells
        for(std::size_t j = 0; j < ownNum; ++j)
        {
            const std::size_t cell = toMerged[j];
            face_vec &cellFaces = merged.cellFaces[cell];
            for(const std::size_t &face : subdomain.GetCellFaces(j))
            {
                cellFaces.push_back(faceToMerged[face]);
            }
            merged.volume[cell] = subdomain.GetVolume(j);
            merged.CM[cell] = subdomain.GetCellCM(j);
        }
    }

    // A face between two subdomains that only one of them built (e.g. a face so small
    // that the other dropped it) is a face of both cells
    for(const std::pair<const std::uint64_t, std::size_t> &sharedFace : sharedFaces)
    {
        const std::size_t face = sharedFace.second;
        const std::array<std::size_t, 2> cells = {merged.faceNeighbors[face].first, merged.faceNeighbors[face].second};
        for(const std::size_t cell : cells)
        {
            face_vec &cellFaces = merged.cellFaces[cell];
            if(std::find(cellFaces.begin(), cellFaces.end(), face) == cellFaces.end())
            {
                cellFaces.push_back(face);
            }
        }
    }
    return merged;
}

} // namespace MadVoro

#endif // MADVORO_SUBDOMAINS_MERGE_SUBDOMAINS_HPP
