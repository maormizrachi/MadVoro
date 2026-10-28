// Merges the tessellations of a rank's thread subdomains into one
#ifndef MADVORO_SUBDOMAINS_MERGE_SUBDOMAINS_HPP
#define MADVORO_SUBDOMAINS_MERGE_SUBDOMAINS_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../types.hpp"
#include "../exception/MadVoroException.hpp"
#include "../geometry/FlatVoronoi.hpp"
#include "../geometry/HostGeometryExecutor.hpp"
#include "../range/RangeQueryData.h"
#include "../builder/BuildDomain.hpp"
#include "../builder/DistributedBuild.hpp"

namespace MadVoro
{

// A thread subdomain's tessellation
template <typename PointT>
struct SubdomainTessellation
{
    const std::vector<PointT> *points; // its cells' points, then its ghosts
    std::size_t cellsNum;
    const AllPointsMap *rankIndex; // a point -> its index among the rank's points and the points other ranks sent (none: a mirror)
    const Geometry::FlatVoronoi *voronoi;
    const std::vector<PointT> *CM; // of a mirror ghost: the mirror of its cell's centroid
};

/** \brief A face built by two subdomains is kept from the first. `points`: the built points, 4 unused points and the points other ranks sent; mergedIndexOf: the merged index of the rank's points and the points other ranks sent (INVALID for rank points not built). The ghosts copying rank points follow `points`, then the mirrors. */
template <typename PointT>
void MergeSubdomains(const std::vector<SubdomainTessellation<PointT>> &subdomains, const BuildDomain<PointT> &domain, const std::vector<PointT> &points, std::size_t pointsNum, const std::vector<std::size_t> &mergedIndexOf, std::size_t rankPointsNum, DistributedBuild<PointT> &merged)
{
    constexpr std::size_t INVALID = std::numeric_limits<std::size_t>::max();
    const std::size_t count = subdomains.size();

    merged.cellFaces.clear();
    merged.faceNeighbors.clear();
    merged.faceVertices.clear();
    merged.faceArea.clear();
    merged.faceCM.clear();
    merged.vertices.clear();
    merged.volume.clear();
    merged.CM.clear();
    merged.ghostRankIndex.clear();
    merged.points.assign(points.begin(), points.end());
    const std::size_t rankGhostsBegin = merged.points.size();
    merged.cellFaces.resize(pointsNum);
    merged.volume.resize(pointsNum);
    const auto rankGhostKey = [](std::size_t index, int imageCode)
    {
        return static_cast<std::uint64_t>(index) * NUM_IMAGE_CODES + static_cast<std::uint64_t>(imageCode);
    };
    const auto sharedFaceKey = [pointsNum](std::size_t a, std::size_t b)
    {
        return static_cast<std::uint64_t>(std::min(a, b)) * pointsNum + std::max(a, b);
    };

    struct Part
    {
        std::vector<std::size_t> toMerged;
        std::vector<std::size_t> mirrors;
        std::vector<std::pair<std::uint64_t, std::size_t>> rankGhosts;
        std::unordered_map<std::uint64_t, std::size_t> sharedFaces; // first face of each cell pair across subdomains
        std::vector<std::size_t> keptIndex; // INVALID: a copy
        std::vector<std::pair<std::size_t, std::size_t>> copyOf; // (subdomain, face) of the kept copy
        std::size_t keptCount = 0;
        std::vector<std::pair<std::size_t, std::size_t>> missing; // (cell, face) only the other cell's subdomain built
        bool unknownPoint = false;
    };
    std::vector<Part> parts(count);
    std::vector<std::size_t> cellSubdomain(pointsNum, INVALID);

    // 1. Each subdomain's points and its faces shared with other subdomains
    Geometry::HostParallelFor(count, [&](std::size_t s)
    {
        const SubdomainTessellation<PointT> &subdomain = subdomains[s];
        const Geometry::FlatVoronoi &voronoi = *subdomain.voronoi;
        Part &part = parts[s];
        const std::size_t ownNum = subdomain.cellsNum;
        const std::size_t subPointsNum = subdomain.points->size();
        const AllPointsMap &rankIndex = *subdomain.rankIndex;
        part.toMerged.assign(subPointsNum, INVALID);
        part.mirrors.clear();
        part.rankGhosts.clear();
        for(std::size_t j = 0; j < ownNum; ++j)
        {
            const auto it = rankIndex.find(j);
            if(it == rankIndex.cend() or it->second >= rankPointsNum or mergedIndexOf[it->second] >= pointsNum)
            {
                part.unknownPoint = true;
                continue;
            }
            part.toMerged[j] = mergedIndexOf[it->second];
            cellSubdomain[part.toMerged[j]] = s;
        }
        for(std::size_t j = ownNum; j < subPointsNum; ++j)
        {
            const auto it = rankIndex.find(j);
            if(it == rankIndex.cend())
            {
                part.mirrors.push_back(j);
                continue;
            }
            const std::size_t index = it->second;
            const int imageCode = ImageCode(domain.PeriodicTranslation((*subdomain.points)[j]));
            const bool image = (imageCode != ZERO_IMAGE_CODE);
            if(index >= mergedIndexOf.size())
            {
                part.unknownPoint = true;
            }
            else if(index >= rankPointsNum or (not image and mergedIndexOf[index] != INVALID))
            {
                part.toMerged[j] = mergedIndexOf[index];
            }
            else
            {
                part.rankGhosts.emplace_back(rankGhostKey(index, imageCode), j); // a point not built, or an image
            }
        }
        for(std::size_t f = 0; f < voronoi.FaceCount(); ++f)
        {
            const std::size_t first = voronoi.faceNeighbor0[f], second = voronoi.faceNeighbor1[f];
            const std::size_t a = part.toMerged[first];
            const std::size_t b = part.toMerged[second];
            if(a < pointsNum and b < pointsNum and ((first >= ownNum) or (second >= ownNum)))
            {
                part.sharedFaces.emplace(sharedFaceKey(a, b), f);
            }
        }
    });
    for(const Part &part : parts)
    {
        if(part.unknownPoint)
        {
            throw MadVoro::Exception::MadVoroException("MergeSubdomains: a subdomain's point is not among the rank's points");
        }
    }

    // Ghosts copying rank points, once each, ordered by key
    std::vector<std::uint64_t> rankGhostKeys;
    for(const Part &part : parts)
    {
        for(const auto &[key, j] : part.rankGhosts)
        {
            rankGhostKeys.push_back(key);
        }
    }
    std::sort(rankGhostKeys.begin(), rankGhostKeys.end());
    rankGhostKeys.erase(std::unique(rankGhostKeys.begin(), rankGhostKeys.end()), rankGhostKeys.end());
    merged.points.resize(rankGhostsBegin + rankGhostKeys.size());
    for(std::size_t g = 0; g < rankGhostKeys.size(); ++g)
    {
        merged.ghostRankIndex.emplace_hint(merged.ghostRankIndex.end(), rankGhostsBegin + g, static_cast<std::size_t>(rankGhostKeys[g] / NUM_IMAGE_CODES));
    }
    Geometry::HostParallelFor(count, [&](std::size_t s)
    {
        Part &part = parts[s];
        const std::vector<PointT> &subPoints = *subdomains[s].points;
        for(const auto &[key, j] : part.rankGhosts)
        {
            const std::size_t index = rankGhostsBegin + (std::lower_bound(rankGhostKeys.begin(), rankGhostKeys.end(), key) - rankGhostKeys.begin());
            part.toMerged[j] = index;
            merged.points[index] = subPoints[j];
        }
    });
    const std::size_t mirrorsBegin = merged.points.size();

    // 2. A face shared by two subdomains is kept by the first subdomain that built it
    Geometry::HostParallelFor(count, [&](std::size_t s)
    {
        const Geometry::FlatVoronoi &voronoi = *subdomains[s].voronoi;
        Part &part = parts[s];
        const std::size_t ownNum = subdomains[s].cellsNum;
        part.keptIndex.assign(voronoi.FaceCount(), INVALID);
        part.copyOf.assign(voronoi.FaceCount(), {INVALID, INVALID});
        part.keptCount = 0;
        for(std::size_t f = 0; f < voronoi.FaceCount(); ++f)
        {
            const std::size_t first = voronoi.faceNeighbor0[f], second = voronoi.faceNeighbor1[f];
            const std::size_t a = part.toMerged[first];
            const std::size_t b = part.toMerged[second];
            if(a < pointsNum and b < pointsNum and ((first >= ownNum) or (second >= ownNum)))
            {
                // only the subdomains of its two cells can build it
                const std::uint64_t key = sharedFaceKey(a, b);
                std::size_t owner = s;
                for(const std::size_t cell : {a, b})
                {
                    const std::size_t other = cellSubdomain[cell];
                    if(other < owner and parts[other].sharedFaces.count(key) != 0)
                    {
                        owner = other;
                    }
                }
                const std::size_t kept = parts[owner].sharedFaces.at(key);
                if(owner != s or kept != f)
                {
                    part.copyOf[f] = {owner, kept};
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
        vertexOffset[s + 1] = vertexOffset[s] + subdomains[s].voronoi->vertices.size();
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
        const SubdomainTessellation<PointT> &subdomain = subdomains[s];
        const Geometry::FlatVoronoi &voronoi = *subdomain.voronoi;
        Part &part = parts[s];
        const std::size_t ownNum = subdomain.cellsNum;
        const std::vector<PointT> &subPoints = *subdomain.points;
        for(std::size_t m = 0; m < part.mirrors.size(); ++m)
        {
            const std::size_t j = part.mirrors[m];
            const std::size_t index = mirrorsBegin + mirrorOffset[s] + m;
            part.toMerged[j] = index;
            merged.points[index] = subPoints[j];
            merged.CM[index] = (*subdomain.CM)[j];
        }

        std::vector<std::size_t> faceToMerged(voronoi.FaceCount());
        for(std::size_t f = 0; f < voronoi.FaceCount(); ++f)
        {
            if(part.keptIndex[f] == INVALID)
            {
                const auto [owner, kept] = part.copyOf[f];
                faceToMerged[f] = faceOffset[owner] + parts[owner].keptIndex[kept];
                continue;
            }
            const std::size_t face = faceOffset[s] + part.keptIndex[f];
            faceToMerged[f] = face;
            const std::size_t first = voronoi.faceNeighbor0[f], second = voronoi.faceNeighbor1[f];
            const std::size_t a = part.toMerged[first];
            const std::size_t b = part.toMerged[second];
            merged.faceNeighbors[face] = {a, b};
            point_vec &faceVertices = merged.faceVertices[face];
            faceVertices.clear();
            for(std::size_t k = voronoi.faceVertexOffsets[f]; k < voronoi.faceVertexOffsets[f + 1]; ++k)
            {
                faceVertices.push_back(voronoi.faceVertices[k] + vertexOffset[s]);
            }
            merged.faceArea[face] = voronoi.faceArea[f];
            merged.faceCM[face] = PointT(voronoi.faceCM[f].x, voronoi.faceCM[f].y, voronoi.faceCM[f].z);
            // a shared face the other subdomain did not build (e.g. dropped as tiny) is a face of both cells
            if(a < pointsNum and b < pointsNum and ((first >= ownNum) or (second >= ownNum)))
            {
                const std::uint64_t key = sharedFaceKey(a, b);
                for(const std::size_t cell : {a, b})
                {
                    const std::size_t other = cellSubdomain[cell];
                    if(other != s and other < count and parts[other].sharedFaces.count(key) == 0)
                    {
                        part.missing.emplace_back(cell, face);
                    }
                }
            }
        }
        for(std::size_t v = 0; v < voronoi.vertices.size(); ++v)
        {
            merged.vertices[vertexOffset[s] + v] = PointT(voronoi.vertices[v].x, voronoi.vertices[v].y, voronoi.vertices[v].z);
        }

        for(std::size_t j = 0; j < ownNum; ++j)
        {
            const std::size_t cell = part.toMerged[j];
            face_vec &cellFaces = merged.cellFaces[cell];
            cellFaces.clear();
            for(std::size_t k = voronoi.cellFaceOffsets[j]; k < voronoi.cellFaceOffsets[j + 1]; ++k)
            {
                cellFaces.push_back(faceToMerged[voronoi.cellFaces[k]]);
            }
            merged.volume[cell] = voronoi.cellVolume[j];
            merged.CM[cell] = PointT(voronoi.cellCM[j].x, voronoi.cellCM[j].y, voronoi.cellCM[j].z);
        }
    });

    // 5. Faces only one of their two cells' subdomains built, in a fixed order
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
}

} // namespace MadVoro

#endif // MADVORO_SUBDOMAINS_MERGE_SUBDOMAINS_HPP
