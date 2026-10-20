#ifndef MADVORO_GEOMETRY_KERNELS_CELL_MOMENTS_HPP
#define MADVORO_GEOMETRY_KERNELS_CELL_MOMENTS_HPP

#include "../FlatMesh.hpp"
#include "../FlatVoronoi.hpp"
#include "../../portability/Math.hpp"
#include "TetraVolume.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace detail
{
    // M_PI, under a name that MeshDecomposer3D's PI macro does not replace
    constexpr double CELL_WIDTH_PI = 3.14159265358979323846;

    // Adds the tetrahedron's volume, and its centroid weighted by that volume.
    MADVORO_HOST_DEVICE_TEMPLATE
    template<typename VolumeFn>
    MADVORO_INLINE_FUNCTION GeometryStatus AddTetraMoments(const FlatPoint &p0, const FlatPoint &p1, const FlatPoint &p2,
                                                           const FlatPoint &p3, const VolumeFn &tetraVolume,
                                                           double &volume, FlatPoint &weightedCM)
    {
        double vol = 0;
        if(tetraVolume(p0, p1, p2, p3, vol) != GeometryStatus::Ok)
            return GeometryStatus::NeedsHost;
        double x = 0, y = 0, z = 0;
        x += p0.x; y += p0.y; z += p0.z;
        x += p1.x; y += p1.y; z += p1.z;
        x += p2.x; y += p2.y; z += p2.z;
        x += p3.x; y += p3.y; z += p3.z;
        const FlatPoint tetraCM{x * 0.25, y * 0.25, z * 0.25};
        volume += vol;
        weightedCM = weightedCM + FlatPoint{tetraCM.x * vol, tetraCM.y * vol, tetraCM.z * vol};
        return GeometryStatus::Ok;
    }

    // Sums the moments of the cone from `apex` over every face of the cell.
    MADVORO_HOST_DEVICE_TEMPLATE
    template<typename VolumeFn>
    MADVORO_INLINE_FUNCTION GeometryStatus ConeMoments(flat_index_t cell, const FlatPoint &apex,
                                                       const ConstFlatVoronoiView &voronoi, const VolumeFn &tetraVolume,
                                                       double &volume, FlatPoint &weightedCM)
    {
        volume = 0;
        weightedCM = FlatPoint{0.0, 0.0, 0.0};
        const flat_index_t faceCount = voronoi.CellFaceCount(cell);
        for(flat_index_t k = 0; k < faceCount; ++k)
        {
            const flat_index_t face = voronoi.CellFace(cell, k);
            const flat_index_t vertexCount = voronoi.FaceVertexCount(face);
            const FlatPoint &v0 = voronoi.vertices[voronoi.FaceVertex(face, 0)];
            for(flat_index_t j = 0; j + 2 < vertexCount; ++j)
            {
                if(AddTetraMoments(v0, voronoi.vertices[voronoi.FaceVertex(face, j + 1)],
                                   voronoi.vertices[voronoi.FaceVertex(face, j + 2)], apex, tetraVolume,
                                   volume, weightedCM) != GeometryStatus::Ok)
                    return GeometryStatus::NeedsHost;
            }
        }
        return GeometryStatus::Ok;
    }
} // namespace detail

/*! \brief Volume and centroid of a cell, from cones over its faces with the
    generator as apex. When the centroid lands more than 0.4 of the cell's width
    from the generator (a high aspect ratio cell), the cones are recomputed with
    that centroid as apex. Faces are visited in the cell's face order, so the sums
    are formed in the same order as a per-face accumulation over all faces.
    \param tetraVolume Callable (p0, p1, p2, p3, double &volume) -> GeometryStatus;
           TetraVolume in kernels, an exact orient3d on the host.
    \return NeedsHost when tetraVolume could not certify a volume. */
MADVORO_HOST_DEVICE_TEMPLATE
template<typename VolumeFn>
MADVORO_INLINE_FUNCTION GeometryStatus CellMoments(flat_index_t cell, const FlatPoint &generator,
                                                   const ConstFlatVoronoiView &voronoi, const VolumeFn &tetraVolume,
                                                   double &volume, FlatPoint &cm)
{
    FlatPoint weightedCM{0.0, 0.0, 0.0};
    if(detail::ConeMoments(cell, generator, voronoi, tetraVolume, volume, weightedCM) != GeometryStatus::Ok)
        return GeometryStatus::NeedsHost;
    const double inverse = 1.0 / volume;
    cm = FlatPoint{weightedCM.x * inverse, weightedCM.y * inverse, weightedCM.z * inverse};

    const double width = Portable::Pow(3 * volume * 0.25 / detail::CELL_WIDTH_PI, 0.3333333333);
    const FlatPoint offset = cm - generator;
    if(Portable::Sqrt(offset.x * offset.x + offset.y * offset.y + offset.z * offset.z) > 0.4 * width)
    {
        const FlatPoint apex = cm;
        if(detail::ConeMoments(cell, apex, voronoi, tetraVolume, volume, weightedCM) != GeometryStatus::Ok)
            return GeometryStatus::NeedsHost;
        const double inverse2 = 1.0 / volume;
        cm = FlatPoint{weightedCM.x * inverse2, weightedCM.y * inverse2, weightedCM.z * inverse2};
    }
    return GeometryStatus::Ok;
}

//! \brief TetraVolume as a callable, for CellMoments.
struct CertifiedTetraVolume
{
    MADVORO_INLINE_FUNCTION GeometryStatus operator()(const FlatPoint &p0, const FlatPoint &p1, const FlatPoint &p2,
                                                      const FlatPoint &p3, double &volume) const
    {
        return TetraVolume(p0, p1, p2, p3, volume);
    }
};

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_CELL_MOMENTS_HPP
