// Builds cells from a Delaunay tetrahedralization: its circumcenters are the vertices
#ifndef MADVORO_BUILDER_DELAUNAY_LOCAL_BUILDER_HPP
#define MADVORO_BUILDER_DELAUNAY_LOCAL_BUILDER_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include <boost/container/flat_set.hpp>

#include "types.hpp"
#include "builder/LocalBuilder.hpp"
#include "delaunay/Delaunay3D.hpp"
#include "geometry/FlatMesh.hpp"
#include "geometry/FlatVoronoi.hpp"
#include "geometry/GeometryExecutor.hpp"
#include "geometry/HiPrecisionCircumcenter.hpp"
#include "geometry/kernels/Circumcenter.hpp"
#include "utils/container_utils.hpp"
#include "exception/MadVoroException.hpp"

namespace MadVoro
{

template <typename PointT>
class DelaunayLocalBuilder : public LocalBuilder<PointT>
{
public:
    // `onHost`: geometry on the host executor, as in a thread of thread subdomains
    explicit DelaunayLocalBuilder(bool onHost) : onHost(onHost)
    {
    }

    void Build(const std::vector<PointT> &points, const BuildDomain<PointT> &domain, const std::vector<std::size_t> &order) override
    {
        this->cellsNum = points.size();
        this->pointTetras.clear();
        this->R.clear();
        this->centers.clear();
        if(not points.empty())
        {
            this->R.reserve(points.size() * 11);
            this->centers.reserve(points.size() * 11);
        }
        this->del.Clean();
        this->flatMesh = Geometry::FlatMesh();
        this->order = order;
        // a builder with no cells still needs the bounding tetrahedron
        std::pair<PointT, PointT> boundingBox = std::make_pair(domain.ll, domain.ur);
        for(const PointT &point : points)
        {
            boundingBox.first.x = std::min(boundingBox.first.x, point.x);
            boundingBox.second.x = std::max(boundingBox.second.x, point.x);
            boundingBox.first.y = std::min(boundingBox.first.y, point.y);
            boundingBox.second.y = std::max(boundingBox.second.y, point.y);
            boundingBox.first.z = std::min(boundingBox.first.z, point.z);
            boundingBox.second.z = std::max(boundingBox.second.z, point.z);
        }
        this->del.Build(points, boundingBox.second, boundingBox.first, this->order);
        this->R.resize(this->del.tetras_.size());
        ContainerOps::conditional_shrink(this->R);
        this->centers.resize(this->R.size());
        ContainerOps::conditional_shrink(this->centers);
        this->SetPointTetras();
        this->changed = false;
    }

    void Add(const std::vector<PointT> &ghosts) override
    {
        if(ghosts.empty())
        {
            return;
        }
        this->del.BuildExtra(ghosts);
        this->changed = true;
    }

    void CellVertices(std::size_t cell, std::vector<PointT> &vertices) override
    {
        this->Update();
        const tetra_vec &tetras = this->pointTetras[cell];
        vertices.resize(tetras.size());
        for(std::size_t k = 0; k < tetras.size(); ++k)
        {
            const std::size_t tetra = tetras[k];
            if(this->R[tetra] < 0)
            {
                this->ComputeCircumcenter(tetra);
            }
            vertices[k] = this->centers[tetra];
        }
    }

    double CellRadius(std::size_t cell, const PointT &point) override
    {
        this->Update();
        double radius2 = 0;
        for(const std::size_t &tetra : this->pointTetras[cell])
        {
            if(this->R[tetra] < 0)
            {
                this->ComputeCircumcenter(tetra);
            }
            const PointT d = this->centers[tetra] - point;
            radius2 = std::max(radius2, ScalarProd(d, d));
        }
        return std::sqrt(radius2);
    }

    void Finish(Geometry::RawFacesView &raw) override
    {
        this->Update();
        Geometry::FlatVoronoi &voronoi = this->vertexStore;
        voronoi = Geometry::FlatVoronoi();
        raw = Geometry::RawFacesView();
        raw.order = this->order.data();
        if(this->cellsNum == 0)
        {
            return;
        }
        // PointTetras' order decides the face numbering
        this->flatMesh = Geometry::MakeFlatMesh(this->del, this->pointTetras);
        const std::size_t tetrasNum = this->del.tetras_.size();
        voronoi.ResizeCellsAndTetras(this->cellsNum, tetrasNum);
        const Geometry::FlatVoronoiView view = voronoi.View();

        // a positive R marks a valid cached center
        const bool reuseCenters = this->del.newTetras_.empty() and this->del.changed_tetras_.empty() and this->R.size() >= tetrasNum;
        std::vector<std::uint8_t> computeCenter(tetrasNum, 0);
        Geometry::HostParallelFor(tetrasNum, [&](std::size_t i)
        {
            if(not this->HasCell(this->del.tetras_[i]))
            {
                return;
            }
            if(reuseCenters and this->R[i] > 0)
            {
                view.vertices[i] = Geometry::ToFlatPoint(this->centers[i]);
            }
            else
            {
                computeCenter[i] = 1;
            }
        });
        this->WithExecutor([&](auto &executor)
        {
            executor.SetMesh(this->flatMesh);
            executor.Circumcenters(computeCenter, voronoi);
        });
        Geometry::HostParallelFor(tetrasNum, [&](std::size_t i)
        {
            if(computeCenter[i] != 0)
            {
                if(view.tetraStatus[i] == Geometry::GeometryStatus::Ok)
                {
                    this->centers[i] = PointT(view.vertices[i].x, view.vertices[i].y, view.vertices[i].z);
                }
                return;
            }
            #ifndef NDEBUG
                if(this->HasCell(this->del.tetras_[i]))
                {
                    const Geometry::FlatMeshView mesh = this->flatMesh.View();
                    Geometry::FlatPoint check{0.0, 0.0, 0.0};
                    double radius = 0;
                    const Tetrahedron &tetra = this->del.tetras_[i];
                    if(Geometry::Circumcenter(mesh.points[tetra.points[0]], mesh.points[tetra.points[1]], mesh.points[tetra.points[2]], mesh.points[tetra.points[3]], check, radius) == Geometry::GeometryStatus::Ok)
                    {
                        assert(check.x == view.vertices[i].x and check.y == view.vertices[i].y and check.z == view.vertices[i].z);
                    }
                }
            #endif // NDEBUG
        });
        for(std::size_t i = 0; i < tetrasNum; ++i)
        {
            if(view.tetraStatus[i] == Geometry::GeometryStatus::NeedsHost)
            {
                this->HiPrecisionCircumcenter(i);
                view.vertices[i] = Geometry::ToFlatPoint(this->centers[i]);
                view.tetraStatus[i] = Geometry::GeometryStatus::Ok;
            }
        }
        this->WithExecutor([&voronoi](auto &executor) { executor.SetVertices(voronoi); });

        const Geometry::FaceBuildReport report = this->WithExecutor([&](auto &executor) { return executor.BuildRawFaces(this->order, raw); });
        if(report.result.error != Geometry::FaceBuildError::None)
        {
            this->ThrowFaceError(report);
        }
        // the bounding points sit between the cells and the ghosts
        raw.gapBegin = this->cellsNum;
        raw.gapSize = BOUNDING_POINTS;
    }

    void ShrinkToFit(void) override
    {
        this->pointTetras.shrink_to_fit();
        this->R.shrink_to_fit();
        this->centers.shrink_to_fit();
        this->del.ShrinkToFit();
    }

    void ReleaseMemory(void) override
    {
        ContainerOps::release_container_memory(this->pointTetras);
        ContainerOps::release_container_memory(this->R);
        ContainerOps::release_container_memory(this->centers);
        this->flatMesh = Geometry::FlatMesh();
        this->vertexStore = Geometry::FlatVoronoi();
        this->WithExecutor([](auto &executor) { executor.Release(); });
        this->del.ReleaseMemory();
    }

    //! \brief The mesh the last Finish built the cells from
    const Geometry::FlatMesh &GetFlatMesh(void) const
    {
        return this->flatMesh;
    }

private:
    static constexpr std::size_t BOUNDING_POINTS = 4;

    template <typename F>
    decltype(auto) WithExecutor(F &&f)
    {
        #if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
            if(this->onHost)
            {
                return f(this->hostExecutor);
            }
        #endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST
        return f(this->executor);
    }

    bool HasCell(const Tetrahedron &tetra) const
    {
        for(std::size_t i = 0; i < 4; ++i)
        {
            if(tetra.points[i] < this->cellsNum)
            {
                return true;
            }
        }
        return false;
    }

    // the tetrahedra the added points changed, before the next use
    void Update(void)
    {
        if(not this->changed)
        {
            return;
        }
        this->R.resize(this->del.tetras_.size(), RADIUS_UNINITIALIZED);
        ContainerOps::conditional_shrink(this->R);
        this->centers.resize(this->R.size());
        ContainerOps::conditional_shrink(this->centers);
        this->SetPointTetras();
        this->changed = false;
    }

    // updates each cell's tetrahedra after the Delaunay changed
    void SetPointTetras(void)
    {
        std::vector<Tetrahedron> &tetras = this->del.tetras_;
        this->pointTetras.resize(this->cellsNum);
        ContainerOps::conditional_shrink(this->pointTetras);
        // empty tetrahedra are not relevant
        for(const std::size_t &empty : this->del.empty_tetras_)
        {
            Tetrahedron &tetra = tetras[empty];
            for(std::size_t i = 0; i < 4; ++i)
            {
                if(tetra.points[i] < this->cellsNum)
                {
                    tetra_vec &tetrasOfPoint = this->pointTetras[tetra.points[i]];
                    auto it = std::find(tetrasOfPoint.begin(), tetrasOfPoint.end(), empty);
                    if(it != tetrasOfPoint.end())
                    {
                        *it = tetrasOfPoint.back();
                        tetrasOfPoint.pop_back();
                    }
                }
                tetra.points[i] = std::numeric_limits<std::size_t>::max();
                tetra.neighbors[i] = std::numeric_limits<std::size_t>::max();
            }
        }
        for(const std::pair<std::size_t, Tetrahedron> &changedTetra : this->del.changed_tetras_)
        {
            const std::size_t &tetraIndex = changedTetra.first;
            const Tetrahedron &tetra = changedTetra.second;
            for(std::size_t i = 0; i < 4; ++i)
            {
                if(tetra.points[i] < this->cellsNum)
                {
                    tetra_vec &tetrasOfPoint = this->pointTetras[tetra.points[i]];
                    auto it = std::find(tetrasOfPoint.begin(), tetrasOfPoint.end(), tetraIndex);
                    if(it != tetrasOfPoint.end())
                    {
                        *it = tetrasOfPoint.back();
                        tetrasOfPoint.pop_back();
                    }
                }
            }
        }
        for(const std::size_t &i : this->del.newTetras_)
        {
            Tetrahedron &tetra = tetras[i];
            this->R[i] = RADIUS_UNINITIALIZED;
            if(not tetra.newTetra)
            {
                continue;
            }
            tetra.newTetra = false;
            for(int j = 0; j < 4; ++j)
            {
                if(tetra.points[j] < this->cellsNum)
                {
                    this->pointTetras[tetra.points[j]].push_back(i);
                }
            }
        }
        this->del.changed_tetras_.clear();
        this->del.newTetras_.clear();
    }

    // a tetrahedron's circumcenter and circumradius, on first use
    void ComputeCircumcenter(std::size_t index)
    {
        const Tetrahedron &tetra = this->del.tetras_[index];
        Geometry::FlatPoint center{0.0, 0.0, 0.0};
        double radius = 0;
        const Geometry::GeometryStatus status = Geometry::Circumcenter(Geometry::ToFlatPoint(this->del.points_[tetra.points[0]]), Geometry::ToFlatPoint(this->del.points_[tetra.points[1]]), Geometry::ToFlatPoint(this->del.points_[tetra.points[2]]), Geometry::ToFlatPoint(this->del.points_[tetra.points[3]]), center, radius);
        this->centers[index] = PointT(center.x, center.y, center.z);
        this->R[index] = (status == Geometry::GeometryStatus::Ok) ? radius : this->HiPrecisionCircumcenter(index);
        if(std::isnan(this->R[index]) or not std::isfinite(this->R[index]))
        {
            MadVoro::Exception::MadVoroException eo("DelaunayLocalBuilder: the circumradius of a tetrahedron is invalid");
            eo.addEntry("Radius", this->R[index]);
            eo.addEntry("Tetra Index", index);
            eo.addEntry("Tetra Points Indices", std::vector({tetra.points[0], tetra.points[1], tetra.points[2], tetra.points[3]}));
            eo.addEntry("Tetra Points", std::vector({this->del.points_[tetra.points[0]], this->del.points_[tetra.points[1]], this->del.points_[tetra.points[2]], this->del.points_[tetra.points[3]]}));
            eo.addEntry("Cells", this->cellsNum);
            throw eo;
        }
    }

    // the circumcenter of a tetrahedron the double precision kernel cannot certify; returns the circumradius
    double HiPrecisionCircumcenter(std::size_t index)
    {
        const Tetrahedron &tetra = this->del.tetras_[index];
        std::array<Geometry::FlatPoint, 4> points;
        for(int k = 0; k < 4; ++k)
        {
            points[k] = Geometry::ToFlatPoint(this->del.points_[tetra.points[k]]);
        }
        Geometry::FlatPoint center{0.0, 0.0, 0.0};
        const double radius = Geometry::HiPrecisionCircumcenter(points, center);
        this->centers[index] = PointT(center.x, center.y, center.z);
        return radius;
    }

    void ThrowFaceError(const Geometry::FaceBuildReport &report) const
    {
        const std::size_t point = report.point;
        const std::size_t other = report.result.other;
        std::string message = "DelaunayLocalBuilder: face construction failed";
        switch(report.result.error)
        {
            case Geometry::FaceBuildError::BadNextLoopTetra:
                message = "Bad NextLoopTetra";
                break;
            case Geometry::FaceBuildError::FanTooLong:
                message = "DelaunayLocalBuilder: the tetrahedron fan around a face visited more tetrahedra than the point has";
                break;
            case Geometry::FaceBuildError::BigTetraNeighbor:
                message = "Neighboring big tet point";
                break;
            case Geometry::FaceBuildError::DegenerateFace:
                message = "MakeRightHandFace: every edge of the face is degenerate";
                break;
            case Geometry::FaceBuildError::NoSecondEdge:
                message = "MakeRightHandFace: could not find a second non degenerate edge of the face";
                break;
            default:
                break;
        }
        MadVoro::Exception::MadVoroException eo(message);
        if(report.result.error == Geometry::FaceBuildError::FanTooLong)
        {
            eo.addEntry("fan size", report.result.fanSize);
            eo.addEntry("tetras of point", this->pointTetras[point].size());
        }
        eo.addEntry("point", point);
        eo.addEntry("point_other", other);
        eo.addEntry("Cells", this->cellsNum);
        eo.addEntry("Point", this->del.points_[point]);
        if(other < this->del.points_.size())
        {
            eo.addEntry("Point other", this->del.points_[other]);
        }
        throw eo;
    }

    bool onHost;
    std::size_t cellsNum = 0;
    // the cells' points, the 4 bounding points, then the ghosts
    Delaunay3D<PointT> del;
    std::vector<std::size_t> order; // the cells' Hilbert order
    std::vector<tetra_vec> pointTetras; // each cell's tetrahedra
    std::vector<double> R; // each tetrahedron's circumradius (negative: not computed)
    std::vector<PointT> centers; // each tetrahedron's circumcenter
    bool changed = false; // points were added since the tetrahedra lists were updated
    Geometry::FlatMesh flatMesh;
    Geometry::FlatVoronoi vertexStore; // the last Finish's vertices (the circumcenters), which its raw faces refer to
    Geometry::GeometryExecutor executor;
    #if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
        Geometry::HostGeometryExecutor hostExecutor;
    #endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST
};

} // namespace MadVoro

#endif // MADVORO_BUILDER_DELAUNAY_LOCAL_BUILDER_HPP
