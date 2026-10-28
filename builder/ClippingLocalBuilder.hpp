// Builds each cell on its own by clipping a cube with the domain's walls and the bisectors of its neighbors (Ray et al. 2018), in Kokkos' default execution space or on the host
#ifndef MADVORO_BUILDER_CLIPPING_LOCAL_BUILDER_HPP
#define MADVORO_BUILDER_CLIPPING_LOCAL_BUILDER_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <cstdio>
#include <unistd.h>

#include "builder/LocalBuilder.hpp"
#include "geometry/ClipExecutor.hpp"
#include "geometry/FlatVoronoi.hpp"
#include "exception/MadVoroException.hpp"

namespace MadVoro
{

template <typename PointT>
class ClippingLocalBuilder : public LocalBuilder<PointT>
{
public:
    // `onHost`: on the host, as in a thread of thread subdomains
    explicit ClippingLocalBuilder(bool onHost) : onHost(onHost)
    {
    }

    void Build(const std::vector<PointT> &points, const BuildDomain<PointT> &domain, const std::vector<std::size_t> &order) override
    {
        this->cellsNum = points.size();
        this->points.resize(points.size());
        for(std::size_t i = 0; i < points.size(); ++i)
        {
            this->points[i] = Geometry::ToFlatPoint(points[i]);
        }
        this->beyondWall.assign(points.size(), 0);
        this->order = order;
        this->SetWalls(domain);
        // every cell starts as the cube of the domain and the cells, widened by their width on each side
        Geometry::FlatPoint low = Geometry::ToFlatPoint(domain.ll), high = Geometry::ToFlatPoint(domain.ur);
        for(const Geometry::FlatPoint &p : this->points)
        {
            low = Geometry::FlatPoint{std::min(low.x, p.x), std::min(low.y, p.y), std::min(low.z, p.z)};
            high = Geometry::FlatPoint{std::max(high.x, p.x), std::max(high.y, p.y), std::max(high.z, p.z)};
        }
        double width = std::max({high.x - low.x, high.y - low.y, high.z - low.z});
        if(not (width > 0))
        {
            width = 1;
        }
        this->mirrorOf.assign(this->cellsNum * this->walls.count, Geometry::FLAT_INVALID_INDEX);
        const Geometry::FlatPoint margin{width, width, width};
        this->WithExecutor([&](auto &executor) { executor.Reset(this->points, this->cellsNum, low - margin, high + margin, this->walls); });
        this->firstNew = 0;
        this->changed = true;
    }

    void Add(const std::vector<PointT> &ghosts) override
    {
        if(ghosts.empty())
        {
            return;
        }
        for(const PointT &ghost : ghosts)
        {
            this->points.push_back(Geometry::ToFlatPoint(ghost));
            this->beyondWall.push_back(this->IsBeyondWall(this->points.back()) ? 1 : 0);
        }
        this->changed = true;
    }

    void AddMirrors(const std::vector<PointT> &mirrors, const std::vector<MirrorSource> &sources) override
    {
        const std::size_t first = this->points.size();
        this->Add(mirrors);
        for(std::size_t i = 0; i < sources.size(); ++i)
        {
            if(sources[i].cell == MirrorSource::NOT_A_MIRROR or sources[i].cell >= this->cellsNum or sources[i].face >= this->wallOfFace.size() or this->wallOfFace[sources[i].face] < 0)
            {
                continue;
            }
            Geometry::flat_index_t &mirror = this->mirrorOf[sources[i].cell * this->walls.count + this->wallOfFace[sources[i].face]];
            if(mirror == Geometry::FLAT_INVALID_INDEX)
            {
                mirror = first + i;
            }
        }
    }

    void CellVertices(std::size_t cell, std::vector<PointT> &vertices) override
    {
        this->Update();
        if(this->fetchedOf[cell] == Geometry::FLAT_INVALID_INDEX)
        {
            this->FetchVertices(std::vector<std::size_t>{cell});
        }
        const flat_index_t k = this->fetchedOf[cell];
        vertices.clear();
        for(flat_index_t v = this->fetchedOffsets[k]; v < this->fetchedOffsets[k + 1]; ++v)
        {
            const Geometry::FlatPoint &vertex = this->fetchedPositions[v];
            vertices.push_back(PointT(vertex.x, vertex.y, vertex.z));
        }
    }

    double CellRadius(std::size_t cell, const PointT &) override
    {
        this->Update();
        return std::sqrt(this->radius2[cell]);
    }

    void FetchVertices(const std::vector<std::size_t> &cells) override
    {
        this->Update();
        for(const flat_index_t &cell : this->fetchedCells)
        {
            this->fetchedOf[cell] = Geometry::FLAT_INVALID_INDEX;
        }
        this->fetchedCells.assign(cells.begin(), cells.end());
        this->fetchedOffsets.assign(1, 0);
        for(std::size_t k = 0; k < cells.size(); ++k)
        {
            this->fetchedOf[cells[k]] = k;
            this->fetchedOffsets.push_back(this->fetchedOffsets.back() + this->vertexCounts[cells[k]]);
        }
        this->WithExecutor([&](auto &executor) { executor.Vertices(this->fetchedCells, this->fetchedOffsets, this->fetchedPositions); });
    }

    void Finish(Geometry::RawFacesView &raw) override
    {
        this->Update();
        raw = Geometry::RawFacesView();
        raw.order = this->order.data();
        if(this->cellsNum == 0)
        {
            return;
        }
        auto mark = std::chrono::high_resolution_clock::now();
        const Geometry::ClipFinishReport report = this->WithExecutor([&](auto &executor) { return executor.Finish(this->order, this->mirrorOf, raw); });
        this->Lap(mark, this->times[3]);
        if(report.error != Geometry::ClipFinishError::None)
        {
            this->ThrowFinishError(report);
        }
        // MADVORO_CLIP_TIMING=1 prints where the build's time went
        static const bool timing = std::getenv("MADVORO_CLIP_TIMING") != nullptr;
        if(timing)
        {
            std::cout << "Clipping: cells " << this->cellsNum << " points " << this->points.size() << " updates " << this->updates << " grids " << this->times[0] << " s, clipping " << this->times[1] << " s (" << this->resolvedWithExpansions << " cells again in expansions, " << this->resolved << " on the host, " << this->resolvedForVertices << " for a vertex), radii " << this->times[2] << " s, raw faces " << this->times[3] << " s" << std::endl;
        }
        this->times.fill(0);
        this->updates = 0;
        this->resolved = 0;
        this->resolvedForVertices = 0;
        this->resolvedWithExpansions = 0;
    }

    void ShrinkToFit(void) override
    {
        this->points.shrink_to_fit();
        this->beyondWall.shrink_to_fit();
        this->vertexCounts.shrink_to_fit();
        this->radius2.shrink_to_fit();
        this->fetchedOf.shrink_to_fit();
        this->fetchedCells.shrink_to_fit();
        this->fetchedOffsets.shrink_to_fit();
        this->fetchedPositions.shrink_to_fit();
    }

    void ReleaseMemory(void) override
    {
        this->points = std::vector<Geometry::FlatPoint>();
        this->beyondWall = std::vector<std::uint8_t>();
        this->mirrorOf = std::vector<Geometry::flat_index_t>();
        this->order = std::vector<std::size_t>();
        this->vertexCounts = std::vector<std::uint8_t>();
        this->radius2 = std::vector<double>();
        this->fetchedOf = std::vector<flat_index_t>();
        this->fetchedCells = std::vector<flat_index_t>();
        this->fetchedOffsets = std::vector<flat_index_t>();
        this->fetchedPositions = std::vector<Geometry::FlatPoint>();
        this->WithExecutor([](auto &executor) { executor.Release(); });
        this->cellsNum = 0;
    }

private:
    using flat_index_t = Geometry::flat_index_t;

    // the domain's faces that are not periodic, with their normals out of it
    void SetWalls(const BuildDomain<PointT> &domain)
    {
        this->walls = Geometry::ClipWalls();
        Geometry::FlatPoint inside{0.0, 0.0, 0.0};
        std::size_t vertexCount = 0;
        for(const Face3D<PointT> &face : domain.box)
        {
            for(const PointT &vertex : face.vertices)
            {
                inside = inside + Geometry::ToFlatPoint(vertex);
                ++vertexCount;
            }
        }
        if(vertexCount > 0)
        {
            inside = (1.0 / static_cast<double>(vertexCount)) * inside;
        }
        this->wallOfFace.assign(domain.box.size(), -1);
        for(std::size_t f = 0; f < domain.box.size(); ++f)
        {
            const Face3D<PointT> &face = domain.box[f];
            const Geometry::FlatPoint origin = Geometry::ToFlatPoint(face.vertices[0]);
            Geometry::FlatPoint normal = Geometry::Cross(Geometry::ToFlatPoint(face.vertices[1]) - origin, Geometry::ToFlatPoint(face.vertices[2]) - origin);
            normal = (1.0 / std::sqrt(Geometry::Dot(normal, normal))) * normal;
            const double components[3] = {std::abs(normal.x), std::abs(normal.y), std::abs(normal.z)};
            const int axis = static_cast<int>(std::max_element(components, components + 3) - components);
            if(domain.periodic[axis])
            {
                continue;
            }
            if(this->walls.count == Geometry::CLIP_MAX_WALLS)
            {
                MadVoro::Exception::MadVoroException eo("ClippingLocalBuilder: the domain has more walls than the capacity");
                eo.addEntry("capacity", Geometry::CLIP_MAX_WALLS);
                throw eo;
            }
            if(Geometry::Dot(inside - origin, normal) > 0)
            {
                normal = -1.0 * normal;
            }
            this->wallOfFace[f] = this->walls.count;
            this->walls.point[this->walls.count] = origin;
            this->walls.normal[this->walls.count] = normal;
            ++this->walls.count;
        }
    }

    // a ghost beyond a wall is a mirror: inside, the walls cut as much as the mirrors
    bool IsBeyondWall(const Geometry::FlatPoint &point) const
    {
        for(int wall = 0; wall < this->walls.count; ++wall)
        {
            if(Geometry::Dot(point - this->walls.point[wall], this->walls.normal[wall]) > 0)
            {
                return true;
            }
        }
        return false;
    }

    template <typename F>
    decltype(auto) WithExecutor(F &&f)
    {
        #if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
            if(not this->onHost)
            {
                return f(this->deviceExecutor);
            }
        #endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST
        return f(this->hostExecutor);
    }

    // clips the cells by the points added since the last update
    void Update(void)
    {
        if(not this->changed)
        {
            return;
        }
        auto mark = std::chrono::high_resolution_clock::now();
        double clippingTime = 0;
        const Geometry::ClipGrid all = Geometry::MakeClipGrid(this->points, 0, this->points.size(), this->beyondWall);
        const bool someAdded = (this->firstNew > 0) and (this->firstNew < this->points.size());
        const Geometry::ClipGrid added = someAdded ? Geometry::MakeClipGrid(this->points, this->firstNew, this->points.size(), this->beyondWall) : Geometry::ClipGrid();
        this->Lap(mark, this->times[0]);
        this->WithExecutor([&](auto &executor)
        {
            executor.Update(all, someAdded ? &added : nullptr, this->firstNew);
            this->resolved += executor.Resolved();
            this->resolvedForVertices += executor.ResolvedForVertices();
            this->resolvedWithExpansions += executor.ResolvedWithExpansions();
            clippingTime = this->times[1];
            this->Lap(mark, this->times[1]);
            clippingTime = this->times[1] - clippingTime;
            executor.Radii(this->vertexCounts, this->radius2);
            this->Lap(mark, this->times[2]);
        });
        this->fetchedOf.assign(this->cellsNum, Geometry::FLAT_INVALID_INDEX);
        this->fetchedCells.clear();
        // MADVORO_CLIP_TIMING=2 also prints each update
        static const bool eachUpdate = std::getenv("MADVORO_CLIP_TIMING") != nullptr and std::string(std::getenv("MADVORO_CLIP_TIMING")) == "2";
        if(eachUpdate)
        {
            std::printf("Clipping update %zu: pid %d cells %zu added %zu clipping %.4f s\n", this->updates, static_cast<int>(getpid()), this->cellsNum, this->points.size() - this->firstNew, clippingTime);
        }
        ++this->updates;
        this->firstNew = this->points.size();
        this->changed = false;
        // a cell that clipped to no vertices failed
        for(flat_index_t c = 0; c < this->cellsNum; ++c)
        {
            if(this->vertexCounts[c] == 0)
            {
                MadVoro::Exception::MadVoroException eo("ClippingLocalBuilder: a cell has no vertices");
                eo.addEntry("cell", c);
                throw eo;
            }
        }
    }

    static void Lap(std::chrono::high_resolution_clock::time_point &mark, double &total)
    {
        const auto now = std::chrono::high_resolution_clock::now();
        total += std::chrono::duration<double>(now - mark).count();
        mark = now;
    }

    void ThrowFinishError(const Geometry::ClipFinishReport &report) const
    {
        std::string message;
        switch(report.error)
        {
            case Geometry::ClipFinishError::CellOverflow:
                message = "a cell has more planes or vertices than the capacity";
                break;
            case Geometry::ClipFinishError::CellUnresolved:
                message = "a cell has a vertex that was not placed or a removal that was not decided";
                break;
            case Geometry::ClipFinishError::ReachesCube:
                message = "a cell reaches the bounding cube (its ghosts are missing)";
                break;
            case Geometry::ClipFinishError::BrokenLoop:
                message = "the vertices around a face do not form one loop";
                break;
            case Geometry::ClipFinishError::NoMirror:
                message = "a face on a wall has no mirror among the ghosts";
                break;
            default:
                message = "building the tessellation failed";
                break;
        }
        MadVoro::Exception::MadVoroException eo("ClippingLocalBuilder: " + message);
        eo.addEntry("cell", report.cell);
        eo.addEntry("point", PointT(this->points[report.cell].x, this->points[report.cell].y, this->points[report.cell].z));
        throw eo;
    }

    bool onHost;
    std::size_t cellsNum = 0;
    std::vector<Geometry::FlatPoint> points; // the cells' points, then the ghosts
    std::vector<std::uint8_t> beyondWall; // per point: a mirror, which clips nothing
    Geometry::ClipWalls walls;
    std::vector<int> wallOfFace; // per face of the domain's box: its wall (-1 for a periodic face)
    std::vector<Geometry::flat_index_t> mirrorOf; // [cell * walls + wall]: the cell's mirror ghost across the wall
    std::vector<std::size_t> order; // the cells' Hilbert order
    std::size_t firstNew = 0; // the points from it on were added after the last update
    bool changed = false;
    std::vector<std::uint8_t> vertexCounts; // each cell's vertex count, from the last update
    std::vector<double> radius2; // the squared distance from each cell's point to its farthest vertex, from the last update
    std::vector<flat_index_t> fetchedOf; // per cell: its index among the fetched cells (FLAT_INVALID_INDEX if not fetched since the last update)
    std::vector<flat_index_t> fetchedCells, fetchedOffsets; // the cells whose vertices were fetched since the last update, and where each one's are
    std::vector<Geometry::FlatPoint> fetchedPositions;
    // grids, clipping, radii, tessellation (seconds)
    std::array<double, 4> times{};
    std::size_t updates = 0;
    std::size_t resolved = 0, resolvedForVertices = 0, resolvedWithExpansions = 0;
    Geometry::HostClipExecutor hostExecutor;
    #if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
        Geometry::KokkosClipExecutor deviceExecutor;
    #endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST
};

} // namespace MadVoro

#endif // MADVORO_BUILDER_CLIPPING_LOCAL_BUILDER_HPP
