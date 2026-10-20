#ifndef MADVORO_GEOMETRY_KERNELS_FACE_BUILD_HPP
#define MADVORO_GEOMETRY_KERNELS_FACE_BUILD_HPP

#include <cstdint>

#include "../FlatMesh.hpp"
#include "../FlatVoronoi.hpp"
#include "../../portability/Math.hpp"

namespace MadVoro
{
namespace Geometry
{

//! \brief Why building the faces of a point failed. Each value except ScratchOverflow is an error the host reports.
enum class FaceBuildError : std::uint8_t
{
    None = 0,
    BadNextLoopTetra,   // a tetrahedron of the ring around an edge has no next tetrahedron
    FanTooLong,         // the ring around an edge visited more tetrahedra than the point has
    BigTetraNeighbor,   // a face separates the point from a bounding-tetrahedron point
    DegenerateFace,     // every edge of a face is degenerate
    NoSecondEdge,       // a face has no second non degenerate edge
    ScratchOverflow     // the scratch arrays are too small for this point (only with fixed-size scratch)
};

struct FaceBuildResult
{
    FaceBuildError error = FaceBuildError::None;
    flat_index_t other = FLAT_INVALID_INDEX;  //!< The neighbor whose face failed
    flat_index_t fanSize = 0;                 //!< Size of the ring when it failed
};

/*! \brief Scratch arrays for building the faces of one point. `capacity` bounds a
    ring of tetrahedra around one edge (and so the vertices of one face), and
    `acceptedCapacity` bounds the neighbors of the point. On the host they can be
    sized from the mesh so they never overflow. */
struct FaceScratch
{
    flat_index_t *fan = nullptr;
    flat_index_t *face = nullptr;
    flat_index_t *backup = nullptr;
    flat_index_t *flip = nullptr;
    FlatPoint *points = nullptr;
    double *diffs = nullptr;
    double *areas = nullptr;
    flat_index_t capacity = 0;
    flat_index_t *accepted = nullptr;
    flat_index_t acceptedCapacity = 0;
};

/*! \brief Scratch arrays of fixed size, for kernels that cannot allocate. A point
    whose rings or neighbors do not fit reports ScratchOverflow. */
template<int Capacity, int AcceptedCapacity>
struct FixedFaceScratch
{
    flat_index_t fan[Capacity];
    flat_index_t face[Capacity];
    flat_index_t backup[Capacity];
    flat_index_t flip[Capacity];
    FlatPoint points[Capacity];
    double diffs[Capacity];
    double areas[Capacity];
    flat_index_t accepted[AcceptedCapacity];

    MADVORO_INLINE_FUNCTION FaceScratch View()
    {
        FaceScratch view;
        view.fan = this->fan;
        view.face = this->face;
        view.backup = this->backup;
        view.flip = this->flip;
        view.points = this->points;
        view.diffs = this->diffs;
        view.areas = this->areas;
        view.capacity = Capacity;
        view.accepted = this->accepted;
        view.acceptedCapacity = AcceptedCapacity;
        return view;
    }
};

namespace face_detail
{
    MADVORO_INLINE_FUNCTION double Norm(const FlatPoint &v)
    {
        return Portable::Sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    MADVORO_INLINE_FUNCTION FlatPoint Scale(const FlatPoint &v, double s)
    {
        return FlatPoint{v.x * s, v.y * s, v.z * s};
    }

    //! \brief The next tetrahedron of the ring around edge (N0, N1), or false when there is none.
    MADVORO_INLINE_FUNCTION bool NextLoopTetra(const FlatTetra &tetra, flat_index_t lastTetra, flat_index_t N0,
                                               flat_index_t N1, flat_index_t &next)
    {
        int i = 0;
        for(; i < 4; i++)
        {
            const flat_index_t point = tetra.points[i];
            if(point != N0 && point != N1 && tetra.neighbors[i] != lastTetra)
                break;
        }
        if(i >= 4)
            return false;
        next = tetra.neighbors[i];
        return true;
    }

    /*! \brief Keeps the ring vertices that are not duplicates of the previous one.
        \return The largest squared edge length, at least R */
    MADVORO_INLINE_FUNCTION double CleanDuplicates(const flat_index_t *fan, flat_index_t N, const FlatPoint *centers,
                                                   flat_index_t *res, flat_index_t &resSize, double R,
                                                   double *diffs, FlatPoint *vtemp)
    {
        resSize = 0;
        for(flat_index_t i = 0; i < N; ++i)
            vtemp[i] = centers[fan[i]];
        for(flat_index_t i = N - 1; i > 0; --i)
            vtemp[i] = vtemp[i] - vtemp[i - 1];
        vtemp[0] = vtemp[0] - centers[fan[N - 1]];
        for(flat_index_t i = 0; i < N; ++i)
        {
            diffs[i] = Dot(vtemp[i], vtemp[i]);
            R = Portable::Max(R, diffs[i]);
        }
        for(flat_index_t i = 0; i < N; ++i)
            if(diffs[i] > R * 1e-16)
                res[resSize++] = fan[i];
        return R;
    }

    //! \brief Area and centroid of a face, from a fan of triangles around its first vertex.
    MADVORO_INLINE_FUNCTION void CalcFaceAreaCM(const flat_index_t *face, flat_index_t N, const FlatPoint *centers,
                                                FlatPoint *points, double *Atemp, double &Area, FlatPoint &CM)
    {
        for(flat_index_t i = 0; i < N; i++)
            points[i] = centers[face[i]];
        const flat_index_t Nloop = N - 2;
        Area = 0;
        for(flat_index_t i = 0; i < Nloop; i++)
        {
            const FlatPoint temp4{points[i + 1].x - points[0].x, points[i + 1].y - points[0].y, points[i + 1].z - points[0].z};
            const FlatPoint temp5{points[i + 2].x - points[0].x, points[i + 2].y - points[0].y, points[i + 2].z - points[0].z};
            const FlatPoint temp3 = Cross(temp4, temp5);
            Atemp[i] = 0.3333333333333333 * 0.5 * Portable::Sqrt(Dot(temp3, temp3));
        }
        double x = 0, y = 0, z = 0;
        for(flat_index_t i = 0; i < Nloop; i++)
        {
            const double A = Atemp[i];
            x += A * points[0].x;
            y += A * points[0].y;
            z += A * points[0].z;
            x += A * points[i + 1].x;
            y += A * points[i + 1].y;
            z += A * points[i + 1].z;
            x += A * points[i + 2].x;
            y += A * points[i + 2].y;
            z += A * points[i + 2].z;
            Area += 3.0 * A;
        }
        CM = Scale(FlatPoint{x, y, z}, 1.0 / (Area + Portable::DOUBLE_MIN * 100));
    }

    //! \brief Orders the face vertices right handed with respect to `point`.
    MADVORO_INLINE_FUNCTION FaceBuildError MakeRightHandFace(flat_index_t *face, flat_index_t N, const FlatPoint &point,
                                                             const FlatPoint *centers, flat_index_t *temp,
                                                             double areascale)
    {
        flat_index_t counter = 0;
        FlatPoint V1 = centers[face[counter + 1]] - centers[face[counter]];
        const double AScale = 1e-14 * areascale;
        // These loops skip degenerate edges, bounded so that a face whose every edge
        // is below the area scale does not index past the face.
        while(Dot(V1, V1) < AScale)
        {
            ++counter;
            if(counter >= N)
                return FaceBuildError::DegenerateFace;
            V1 = centers[face[(counter + 1) % N]] - centers[face[counter % N]];
        }
        FlatPoint V2 = centers[face[(counter + 2) % N]] - centers[face[(counter + 1) % N]];
        while(Dot(V2, V2) < AScale)
        {
            ++counter;
            if(counter >= 2 * N)
                return FaceBuildError::NoSecondEdge;
            V2 = centers[face[(counter + 2) % N]] - centers[face[(counter + 1) % N]];
        }
        // Do we need to flip handedness?
        if(Dot(Cross(V1, V2), point - centers[face[0]]) > 0)
        {
            for(flat_index_t j = 0; j < N; ++j)
                temp[j] = face[j];
            for(flat_index_t i = 0; i < N; ++i)
                face[i] = temp[N - i - 1];
        }
        return FaceBuildError::None;
    }

    MADVORO_INLINE_FUNCTION void Erase(flat_index_t *list, flat_index_t &size, flat_index_t i)
    {
        for(flat_index_t j = i; j + 1 < size; ++j)
            list[j] = list[j + 1];
        --size;
    }

    MADVORO_INLINE_FUNCTION FlatPoint CornerNormal(const flat_index_t *face, flat_index_t i, flat_index_t N,
                                                   const FlatPoint *centers)
    {
        return Cross(centers[face[i]] - centers[face[(N + i - 1) % N]],
                     centers[face[(i + 1) % N]] - centers[face[(N + i - 1) % N]]);
    }

    /*! \brief Removes face vertices that are collinear with their neighbors or bend
        against the face normal.
        \return False when fewer than three vertices would remain (the face is dropped) */
    MADVORO_INLINE_FUNCTION bool CleanSameLine(flat_index_t *face, flat_index_t &size, const FlatPoint *centers,
                                               double *area_vec_temp, flat_index_t *old)
    {
        const flat_index_t N = size;
        const double small_fraction = 1e-14;
        // Find correct normal
        for(flat_index_t i = 0; i < N; ++i)
        {
            area_vec_temp[i] = Norm(CornerNormal(face, i, N, centers));
            old[i] = face[i];
        }

        double max_value = area_vec_temp[0];
        double second_max_value = max_value;
        flat_index_t max_index = 0, second_max_index = 0;
        for(flat_index_t i = 1; i < N; ++i)
        {
            if(area_vec_temp[i] > max_value)
            {
                second_max_value = max_value;
                max_value = area_vec_temp[i];
                second_max_index = max_index;
                max_index = i;
            }
            else
            {
                if(area_vec_temp[i] > second_max_value)
                {
                    second_max_value = area_vec_temp[i];
                    second_max_index = i;
                }
            }
        }

        const double area_scale = area_vec_temp[max_index];
        FlatPoint good_normal = CornerNormal(face, max_index, N, centers);
        good_normal = Scale(good_normal, 1.0 / Norm(good_normal));

        // Nindeces is only refreshed after an erase that does not end the loop, as before
        flat_index_t Nindeces = size;
        for(flat_index_t i = 0; i < Nindeces; ++i)
        {
            FlatPoint normal_temp = CornerNormal(face, i, Nindeces, centers);
            const double area = Norm(normal_temp);
            normal_temp = Scale(normal_temp, 1.0 / (100 * Portable::DOUBLE_MIN + area));
            if((area < area_scale * small_fraction) || (Dot(normal_temp, good_normal) < 0.9999))
            {
                Erase(face, size, i);
                if(i == size - 1)
                    break;
                --i;
                Nindeces = size;
            }
        }

        if(Nindeces < 3)
        {
            for(flat_index_t i = 0; i < N; ++i)
                face[i] = old[i];
            size = N;
            max_index = second_max_index;
            good_normal = CornerNormal(face, max_index, N, centers);
            good_normal = Scale(good_normal, 1.0 / Norm(good_normal));

            for(flat_index_t i = 0; i < Nindeces; i++)
            {
                FlatPoint normal_temp = CornerNormal(face, i, Nindeces, centers);
                const double area = Norm(normal_temp);
                normal_temp = Scale(normal_temp, 1.0 / area);
                if((area < area_scale * small_fraction) || Dot(normal_temp, good_normal) < 0.9999)
                {
                    Erase(face, size, i);
                    if(i == size - 1)
                        break;
                    --i;
                    Nindeces = size;
                }
            }
        }

        return Nindeces >= 3;
    }
} // namespace face_detail

//! \brief Counts the faces and face vertices a point creates (first pass).
struct FaceCountSink
{
    flat_index_t faces = 0;
    flat_index_t vertices = 0;

    MADVORO_INLINE_FUNCTION void Add(flat_index_t, flat_index_t, const flat_index_t *, flat_index_t vertexCount,
                                     double, const FlatPoint &)
    {
        ++this->faces;
        this->vertices += vertexCount;
    }
};

/*! \brief Writes the faces a point creates, starting at the face and vertex offsets
    the first pass assigned to it (second pass). */
struct FaceFillSink
{
    FlatVoronoiView voronoi;
    flat_index_t face = 0;
    flat_index_t vertex = 0;

    MADVORO_INLINE_FUNCTION void Add(flat_index_t point, flat_index_t other, const flat_index_t *faceVertices,
                                     flat_index_t vertexCount, double area, const FlatPoint &cm)
    {
        this->voronoi.faceNeighbor0[this->face] = point;
        this->voronoi.faceNeighbor1[this->face] = other;
        this->voronoi.faceArea[this->face] = area;
        this->voronoi.faceCM[this->face] = cm;
        this->voronoi.faceVertexOffsets[this->face] = this->vertex;
        for(flat_index_t k = 0; k < vertexCount; ++k)
            this->voronoi.faceVertices[this->vertex + k] = faceVertices[k];
        this->vertex += vertexCount;
        ++this->face;
    }
};

/*! \brief Builds the faces between `point` and each neighbor with a larger index,
    so every face is built by exactly one of its two points. For each such edge it
    walks the ring of tetrahedra around it, then drops duplicate vertices, faces
    that are too small and collinear vertices, and orders the face right handed
    with respect to `point`. A rejected face is attempted again from the next
    tetrahedron that shares the edge, since its ring then starts elsewhere.
    Faces are handed to `sink` in the order they are accepted.
    \param centers Circumcenters, indexed by tetrahedron
    \param looseAreaCut Per point: 1 for a point outside the domain that is not a
           periodic image, whose faces use the looser small-area cut */
MADVORO_HOST_DEVICE_TEMPLATE
template<typename Sink>
MADVORO_INLINE_FUNCTION FaceBuildResult BuildPointFaces(flat_index_t point, const FlatMeshView &mesh,
                                                        const FlatPoint *centers, const std::uint8_t *looseAreaCut,
                                                        const FaceScratch &s, Sink &sink)
{
    FaceBuildResult result;
    flat_index_t acceptedCount = 0;
    const flat_index_t ntet = mesh.PointTetraCount(point);
    for(flat_index_t j = 0; j < ntet; ++j)
    {
        const flat_index_t tetcheck = mesh.PointTetra(point, j);
        for(int k = 0; k < 4; ++k)
        {
            const flat_index_t other = mesh.tetras[tetcheck].points[k];
            if(other == point || other < point)
                continue;
            // Did we already build this face?
            bool built = false;
            for(flat_index_t a = 0; a < acceptedCount; ++a)
                built = built || (s.accepted[a] == other);
            if(built)
                continue;

            // The ring of tetrahedra around the edge. The walk only ends by returning
            // to its first tetrahedron, so it is bounded by the point's tetrahedra.
            flat_index_t fanSize = 0;
            s.fan[fanSize++] = tetcheck;
            flat_index_t next = 0;
            if(!face_detail::NextLoopTetra(mesh.tetras[tetcheck], tetcheck, point, other, next))
            {
                result.error = FaceBuildError::BadNextLoopTetra;
                result.other = other;
                return result;
            }
            flat_index_t cur = next;
            flat_index_t last = tetcheck;
            while(next != tetcheck)
            {
                if(fanSize > ntet)
                {
                    result.error = FaceBuildError::FanTooLong;
                    result.other = other;
                    result.fanSize = fanSize;
                    return result;
                }
                if(fanSize >= s.capacity)
                {
                    result.error = FaceBuildError::ScratchOverflow;
                    result.other = other;
                    return result;
                }
                s.fan[fanSize++] = cur;
                if(!face_detail::NextLoopTetra(mesh.tetras[cur], last, point, other, next))
                {
                    result.error = FaceBuildError::BadNextLoopTetra;
                    result.other = other;
                    return result;
                }
                last = cur;
                cur = next;
            }
            // Is face too small?
            if(fanSize < 3)
                continue;
            const FlatPoint edge = mesh.points[point] - mesh.points[other];
            flat_index_t faceSize = 0;
            const double Asize = face_detail::CleanDuplicates(s.fan, fanSize, centers, s.face, faceSize,
                                                              Dot(edge, edge), s.diffs, s.points);
            if(faceSize < 3)
                continue;
            double area = 0;
            FlatPoint cm{0.0, 0.0, 0.0};
            face_detail::CalcFaceAreaCM(s.face, faceSize, centers, s.points, s.diffs, area, cm);
            if(area < (Asize * (looseAreaCut[other] ? 1e-14 : 1e-15)))
                continue;
            if(other >= mesh.ownedCount && other < (mesh.ownedCount + 4))
            {
                result.error = FaceBuildError::BigTetraNeighbor;
                result.other = other;
                return result;
            }
            // Make faces right handed
            const FaceBuildError handError =
                face_detail::MakeRightHandFace(s.face, faceSize, mesh.points[point], centers, s.flip, area);
            if(handError != FaceBuildError::None)
            {
                result.error = handError;
                result.other = other;
                return result;
            }
            if(!face_detail::CleanSameLine(s.face, faceSize, centers, s.areas, s.backup))
                continue;
            sink.Add(point, other, s.face, faceSize, area, cm);
            if(acceptedCount >= s.acceptedCapacity)
            {
                result.error = FaceBuildError::ScratchOverflow;
                result.other = other;
                return result;
            }
            s.accepted[acceptedCount++] = other;
        }
    }
    return result;
}

/*! \brief Moves a face centroid onto the plane that bisects the two points of the
    face (this prevents large face velocities for close by points). */
MADVORO_INLINE_FUNCTION FlatPoint ProjectFaceCM(const FlatPoint &cm, const FlatPoint &p0, const FlatPoint &p1)
{
    FlatPoint mid = p0 + p1;
    mid = FlatPoint{mid.x * 0.5, mid.y * 0.5, mid.z * 0.5};
    const FlatPoint norm = p1 - p0;
    const double s = Dot(cm - mid, norm);
    const double d = Dot(norm, norm);
    const FlatPoint shift{norm.x * s / d, norm.y * s / d, norm.z * s / d};
    return cm - shift;
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_FACE_BUILD_HPP
