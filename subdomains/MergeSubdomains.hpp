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
#include "../exception/MadVoroException.hpp"
#include "../geometry/HostGeometryExecutor.hpp"
#include "../range/RangeQueryData.h"

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
    //! The built points, the four bounding points, the points other ranks sent, the ghosts
    //! that are copies of the rank's points (rankPointGhosts), then the mirror ghosts of
    //! every subdomain
    std::vector<PointT> points;
    //! The all-points index of each ghost that is a copy of one of the rank's points: a
    //! point not built (in a partial build), or a periodic image of a point
    std::vector<std::size_t> rankPointGhosts;
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
    points are some of the built points, and its ghosts are other points of the rank
    (the points of other subdomains, or points not built), periodic images of the rank's
    points, points other ranks sent, or mirror ghosts of its own points. Ghosts that are
    copies of the rank's points (not built, or images) are placed once each, ordered by
    point and image, after the points other ranks sent. A face between cells of two subdomains is built by both; the copy of the
    first subdomain that built it is kept, and it is listed by both cells. Each cell's
    volume and centroid are those of the subdomain that owns it.
    The subdomains are merged in parallel (HostParallelFor), with the result of merging
    them one after another: each finds its points and the faces it keeps on its own,
    prefix sums place every subdomain's faces, vertices and mirror ghosts after those of
    the subdomains before it, and each then writes its part.
    \param subdomains The subdomains, in order
    \param points The built points (`pointsNum` of them; together, the subdomains own
           each of them once), four places for the bounding points, then the points
           other ranks sent
    \param mergedIndexOf Each point's index in `points` by its all-points index (the
           rank's `rankPointsNum` points, then the points other ranks sent), INVALID for
           a point of the rank that is not built
    \param storage Arrays whose memory the result reuses (e.g. the last merge's); their
           contents are discarded */
template <typename PointT>
SubdomainMerge<PointT> MergeSubdomains(const std::vector<const Voronoi3D<PointT> *> &subdomains, const std::vector<PointT> &points,
                                       std::size_t pointsNum, const std::vector<std::size_t> &mergedIndexOf, std::size_t rankPointsNum,
                                       SubdomainMerge<PointT> storage = SubdomainMerge<PointT>())
{
    constexpr std::size_t INVALID = std::numeric_limits<std::size_t>::max();
    constexpr std::size_t BOUNDING_POINTS = 4;
    const std::size_t count = subdomains.size();

    SubdomainMerge<PointT> merged = std::move(storage);
    merged.cellFaces.clear();
    merged.faceNeighbors.clear();
    merged.faceVertices.clear();
    merged.faceArea.clear();
    merged.faceCM.clear();
    merged.vertices.clear();
    merged.volume.clear();
    merged.CM.clear();
    merged.rankPointGhosts.clear();
    merged.points.assign(points.begin(), points.end());
    if(not subdomains.empty())
    {
        const Voronoi3D<PointT> &first = *subdomains.front();
        const std::vector<PointT> &firstPoints = first.getMeshPoints();
        std::copy(firstPoints.begin() + first.GetPointNo(), firstPoints.begin() + first.GetPointNo() + BOUNDING_POINTS,
                  merged.points.begin() + pointsNum);
    }
    const std::size_t rankGhostsBegin = merged.points.size();
    merged.cellFaces.resize(pointsNum);
    merged.volume.resize(pointsNum);
    // a ghost that is a copy of one of the rank's points, by its all-points index and image
    const auto rankGhostKey = [](std::size_t index, int imageCode)
    {
        return static_cast<std::uint64_t>(index) * NUM_IMAGE_CODES + static_cast<std::uint64_t>(imageCode);
    };
    // a face between cells of two subdomains, by its (smaller, larger) cell pair
    const auto sharedFaceKey = [pointsNum](std::size_t a, std::size_t b)
    {
        return static_cast<std::uint64_t>(std::min(a, b)) * pointsNum + std::max(a, b);
    };

    // What each subdomain finds on its own
    struct Part
    {
        std::vector<std::size_t> toMerged;   // its points in the merged tessellation (mirrors: INVALID until placed)
        std::vector<std::size_t> mirrors;    // its mirror ghosts, in order
        std::vector<std::pair<std::uint64_t, std::size_t>> rankGhosts; // its copies of the rank's points: (key, point)
        std::unordered_map<std::uint64_t, std::size_t> sharedFaces; // its first face of each pair of cells of two subdomains
        std::vector<std::size_t> keptIndex;  // each face's index among the faces it keeps (INVALID: a copy)
        std::vector<std::pair<std::size_t, std::size_t>> copyOf; // a copy's (subdomain, face) that is kept
        std::size_t keptCount = 0;
        std::vector<std::pair<std::size_t, std::size_t>> missing; // (cell, face) only the other cell's subdomain built
        bool unknownPoint = false;
    };
    std::vector<Part> parts(count);
    std::vector<std::size_t> cellSubdomain(pointsNum, INVALID);

    // 1. Each subdomain's points, and its faces between cells of two subdomains
    Geometry::HostParallelFor(count, [&](std::size_t s)
    {
        const Voronoi3D<PointT> &subdomain = *subdomains[s];
        Part &part = parts[s];
        const std::size_t ownNum = subdomain.GetPointNo();
        const std::size_t subPointsNum = subdomain.getMeshPoints().size();
        const AllPointsMap &indicesInPoints = subdomain.GetIndicesInAllPoints();
        part.toMerged.assign(subPointsNum, INVALID);
        part.mirrors.clear();
        part.rankGhosts.clear();
        for(std::size_t j = 0; j < ownNum; ++j)
        {
            const auto it = indicesInPoints.find(j);
            if(it == indicesInPoints.cend() or it->second >= rankPointsNum or mergedIndexOf[it->second] >= pointsNum)
            {
                part.unknownPoint = true;
                continue;
            }
            part.toMerged[j] = mergedIndexOf[it->second];
            cellSubdomain[part.toMerged[j]] = s;
        }
        for(std::size_t j = ownNum + BOUNDING_POINTS; j < subPointsNum; ++j)
        {
            const auto it = indicesInPoints.find(j);
            if(it == indicesInPoints.cend())
            {
                part.mirrors.push_back(j); // a mirror ghost of one of the subdomain's points
                continue;
            }
            const std::size_t index = it->second;
            const PointT &translation = subdomain.GetPeriodicImageTranslation(j);
            const bool image = (ImageCode(translation) != ZERO_IMAGE_CODE);
            if(index >= mergedIndexOf.size())
                part.unknownPoint = true;
            else if(index >= rankPointsNum or (not image and mergedIndexOf[index] != INVALID))
                part.toMerged[j] = mergedIndexOf[index]; // a point another rank sent, or a point of another subdomain
            else
                part.rankGhosts.emplace_back(rankGhostKey(index, ImageCode(translation)), j); // a point not built, or an image
        }
        const std::vector<std::pair<std::size_t, std::size_t>> &faceNeighbors = subdomain.GetAllFaceNeighbors();
        for(std::size_t f = 0; f < faceNeighbors.size(); ++f)
        {
            const std::size_t a = part.toMerged[faceNeighbors[f].first];
            const std::size_t b = part.toMerged[faceNeighbors[f].second];
            if(a < pointsNum and b < pointsNum and ((faceNeighbors[f].first >= ownNum) or (faceNeighbors[f].second >= ownNum)))
                part.sharedFaces.emplace(sharedFaceKey(a, b), f);
        }
    });
    for(const Part &part : parts)
    {
        if(part.unknownPoint)
            throw MadVoro::Exception::MadVoroException("MergeSubdomains: a subdomain's point is not among the rank's points");
    }

    // The ghosts that are copies of the rank's points, once each, in the order of their keys
    std::vector<std::uint64_t> rankGhostKeys;
    for(const Part &part : parts)
    {
        for(const auto &[key, j] : part.rankGhosts)
            rankGhostKeys.push_back(key);
    }
    std::sort(rankGhostKeys.begin(), rankGhostKeys.end());
    rankGhostKeys.erase(std::unique(rankGhostKeys.begin(), rankGhostKeys.end()), rankGhostKeys.end());
    merged.points.resize(rankGhostsBegin + rankGhostKeys.size());
    merged.rankPointGhosts.resize(rankGhostKeys.size());
    for(std::size_t g = 0; g < rankGhostKeys.size(); ++g)
        merged.rankPointGhosts[g] = static_cast<std::size_t>(rankGhostKeys[g] / NUM_IMAGE_CODES);
    Geometry::HostParallelFor(count, [&](std::size_t s)
    {
        Part &part = parts[s];
        const std::vector<PointT> &subPoints = subdomains[s]->getMeshPoints();
        for(const auto &[key, j] : part.rankGhosts)
        {
            const std::size_t index = rankGhostsBegin + (std::lower_bound(rankGhostKeys.begin(), rankGhostKeys.end(), key) - rankGhostKeys.begin());
            part.toMerged[j] = index;
            merged.points[index] = subPoints[j]; // (the same point in every subdomain that has it)
        }
    });
    const std::size_t mirrorsBegin = merged.points.size();

    // 2. The faces each subdomain keeps: a face between cells of two subdomains is kept by
    //    the first subdomain (in order) that built it, as its first face of that pair
    Geometry::HostParallelFor(count, [&](std::size_t s)
    {
        const Voronoi3D<PointT> &subdomain = *subdomains[s];
        Part &part = parts[s];
        const std::size_t ownNum = subdomain.GetPointNo();
        const std::vector<std::pair<std::size_t, std::size_t>> &faceNeighbors = subdomain.GetAllFaceNeighbors();
        part.keptIndex.assign(faceNeighbors.size(), INVALID);
        part.copyOf.assign(faceNeighbors.size(), {INVALID, INVALID});
        part.keptCount = 0;
        for(std::size_t f = 0; f < faceNeighbors.size(); ++f)
        {
            const std::size_t a = part.toMerged[faceNeighbors[f].first];
            const std::size_t b = part.toMerged[faceNeighbors[f].second];
            if(a < pointsNum and b < pointsNum and ((faceNeighbors[f].first >= ownNum) or (faceNeighbors[f].second >= ownNum)))
            {
                // only the subdomains of its two cells can build it
                const std::uint64_t key = sharedFaceKey(a, b);
                std::size_t owner = s;
                for(const std::size_t cell : {a, b})
                {
                    const std::size_t other = cellSubdomain[cell];
                    if(other < owner and parts[other].sharedFaces.count(key) != 0)
                        owner = other;
                }
                const std::size_t first = parts[owner].sharedFaces.at(key);
                if(owner != s or first != f)
                {
                    part.copyOf[f] = {owner, first};
                    continue;
                }
            }
            part.keptIndex[f] = part.keptCount++;
        }
    });

    // 3. Where each subdomain's faces, vertices and mirror ghosts go
    std::vector<std::size_t> faceOffset(count + 1, 0), vertexOffset(count + 1, 0), mirrorOffset(count + 1, 0);
    for(std::size_t s = 0; s < count; ++s)
    {
        faceOffset[s + 1] = faceOffset[s] + parts[s].keptCount;
        vertexOffset[s + 1] = vertexOffset[s] + subdomains[s]->GetFacePoints().size();
        mirrorOffset[s + 1] = mirrorOffset[s] + parts[s].mirrors.size();
    }
    merged.points.resize(mirrorsBegin + mirrorOffset[count]);
    merged.CM.resize(merged.points.size());
    merged.faceNeighbors.resize(faceOffset[count]);
    merged.faceVertices.resize(faceOffset[count]);
    merged.faceArea.resize(faceOffset[count]);
    merged.faceCM.resize(faceOffset[count]);
    merged.vertices.resize(vertexOffset[count]);

    // 4. Each subdomain writes its part
    Geometry::HostParallelFor(count, [&](std::size_t s)
    {
        const Voronoi3D<PointT> &subdomain = *subdomains[s];
        Part &part = parts[s];
        const std::size_t ownNum = subdomain.GetPointNo();
        const std::vector<PointT> &subPoints = subdomain.getMeshPoints();
        for(std::size_t m = 0; m < part.mirrors.size(); ++m)
        {
            const std::size_t j = part.mirrors[m];
            const std::size_t index = mirrorsBegin + mirrorOffset[s] + m;
            part.toMerged[j] = index;
            merged.points[index] = subPoints[j];
            merged.CM[index] = subdomain.GetCellCM(j);
        }

        const std::vector<std::pair<std::size_t, std::size_t>> &faceNeighbors = subdomain.GetAllFaceNeighbors();
        std::vector<std::size_t> faceToMerged(faceNeighbors.size());
        for(std::size_t f = 0; f < faceNeighbors.size(); ++f)
        {
            if(part.keptIndex[f] == INVALID)
            {
                const auto [owner, first] = part.copyOf[f];
                faceToMerged[f] = faceOffset[owner] + parts[owner].keptIndex[first];
                continue;
            }
            const std::size_t face = faceOffset[s] + part.keptIndex[f];
            faceToMerged[f] = face;
            const std::size_t a = part.toMerged[faceNeighbors[f].first];
            const std::size_t b = part.toMerged[faceNeighbors[f].second];
            merged.faceNeighbors[face] = {a, b};
            point_vec faceVertices = subdomain.GetPointsInFace(f);
            for(std::size_t &vertex : faceVertices)
            {
                vertex += vertexOffset[s];
            }
            merged.faceVertices[face] = std::move(faceVertices);
            merged.faceArea[face] = subdomain.GetArea(f);
            merged.faceCM[face] = subdomain.FaceCM(f);
            // a face between two subdomains that the other did not build (e.g. a face so
            // small that it dropped it) is a face of both cells
            if(a < pointsNum and b < pointsNum and ((faceNeighbors[f].first >= ownNum) or (faceNeighbors[f].second >= ownNum)))
            {
                const std::uint64_t key = sharedFaceKey(a, b);
                for(const std::size_t cell : {a, b})
                {
                    const std::size_t other = cellSubdomain[cell];
                    if(other != s and other < count and parts[other].sharedFaces.count(key) == 0)
                        part.missing.emplace_back(cell, face);
                }
            }
        }
        const std::vector<PointT> &vertices = subdomain.GetFacePoints();
        std::copy(vertices.begin(), vertices.end(), merged.vertices.begin() + vertexOffset[s]);

        // Cells
        for(std::size_t j = 0; j < ownNum; ++j)
        {
            const std::size_t cell = part.toMerged[j];
            face_vec &cellFaces = merged.cellFaces[cell];
            for(const std::size_t &face : subdomain.GetCellFaces(j))
            {
                cellFaces.push_back(faceToMerged[face]);
            }
            merged.volume[cell] = subdomain.GetVolume(j);
            merged.CM[cell] = subdomain.GetCellCM(j);
        }
    });

    // 5. The faces only one of their two cells' subdomains built, in a fixed order
    std::vector<std::pair<std::size_t, std::size_t>> missing;
    for(const Part &part : parts)
    {
        missing.insert(missing.end(), part.missing.begin(), part.missing.end());
    }
    std::sort(missing.begin(), missing.end());
    for(const auto &[cell, face] : missing)
    {
        face_vec &cellFaces = merged.cellFaces[cell];
        if(std::find(cellFaces.begin(), cellFaces.end(), face) == cellFaces.end())
        {
            cellFaces.push_back(face);
        }
    }
    return merged;
}

} // namespace MadVoro

#endif // MADVORO_SUBDOMAINS_MERGE_SUBDOMAINS_HPP
