// Builds the cells of given points, knowing nothing of ranks or ghost searches
#ifndef MADVORO_BUILDER_LOCAL_BUILDER_HPP
#define MADVORO_BUILDER_LOCAL_BUILDER_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include "builder/BuildDomain.hpp"
#include "geometry/RawFaces.hpp"

namespace MadVoro
{

//! \brief What a mirror ghost mirrors: one of the cells, across one of the domain's box faces (NOT_A_MIRROR for a ghost that is not a mirror)
struct MirrorSource
{
    static constexpr std::size_t NOT_A_MIRROR = static_cast<std::size_t>(-1);
    std::size_t cell = NOT_A_MIRROR;
    std::size_t face = 0;
};

template <typename PointT>
class LocalBuilder
{
public:
    virtual ~LocalBuilder() = default;

    //! \brief Starts a build of the cells of `points` (indices [0, points.size())) in `domain`, in `order` (the cells' build order: a Hilbert order)
    virtual void Build(const std::vector<PointT> &points, const BuildDomain<PointT> &domain, const std::vector<std::size_t> &order) = 0;

    //! \brief Adds ghosts: points that are not built, indexed after all the points given so far
    virtual void Add(const std::vector<PointT> &ghosts) = 0;

    //! \brief Adds ghosts that mirror cells across the box's faces (sources[i] for mirrors[i])
    virtual void AddMirrors(const std::vector<PointT> &mirrors, const std::vector<MirrorSource> &)
    {
        this->Add(mirrors);
    }

    //! \brief The current vertices of a cell
    virtual void CellVertices(std::size_t cell, std::vector<PointT> &vertices) = 0;

    //! \brief The distance from a cell's point to its farthest current vertex
    virtual double CellRadius(std::size_t cell, const PointT &point) = 0;

    //! \brief Brings the current vertices of these cells at once, for the CellVertices calls that follow
    virtual void FetchVertices(const std::vector<std::size_t> &)
    {
    }

    //! \brief The built cells as raw faces: the vertices, and each face once with its cell, its other point and its vertex loop (the shared stage finishes them; they stay valid until the next Build)
    virtual void Finish(Geometry::RawFacesView &raw) = 0;

    virtual void ShrinkToFit(void) = 0;

    virtual void ReleaseMemory(void) = 0;
};

} // namespace MadVoro

#endif // MADVORO_BUILDER_LOCAL_BUILDER_HPP
