// Clips each Voronoi cell from a cube by the bisectors of its neighbors, as in Ray et al. 2018: a cell is the triangulation of its planes (each vertex is where 3 planes meet)
#ifndef MADVORO_GEOMETRY_KERNELS_CLIP_CELL_HPP
#define MADVORO_GEOMETRY_KERNELS_CLIP_CELL_HPP

#include <array>
#include <cstdint>

#include "../FlatMesh.hpp"
#include "../FlatVoronoi.hpp"
#include "../HiPrecisionCircumcenter.hpp"
#include "Circumcenter.hpp"
#include "Expansions.hpp"
#include "../../portability/Math.hpp"
#include "../../utils/Predicates3D.hpp"

namespace MadVoro
{
namespace Geometry
{

constexpr int CLIP_MAX_PLANES = 96;
constexpr int CLIP_MAX_VERTICES = 192;
// a plane from it on is a face of the bounding cube, CLIP_CUBE_PLANE + (2 * axis + (upper face ? 1 : 0)): the bisector of the cell's point and its mirror across the face
constexpr flat_index_t CLIP_CUBE_PLANE = FLAT_INVALID_INDEX - 6;
// the domain's rigid walls: CLIP_WALL_PLANE + wall, the bisector of the cell's point and its mirror across the wall
constexpr int CLIP_MAX_WALLS = 16;
constexpr flat_index_t CLIP_WALL_PLANE = CLIP_CUBE_PLANE - CLIP_MAX_WALLS;
constexpr double CLIP_INFINITE = 1e300;
// what decides what floating point cannot: nothing (the cell NeedsHost), exact arithmetic (on the host), or expansions in fixed arrays (in kernels; a cell they do not fit, or with a vertex floating point cannot place, NeedsHost)
constexpr int CLIP_FLOATING = 0;
constexpr int CLIP_EXACT = 1;
constexpr int CLIP_EXPANSIONS = 2;

enum class ClipStatus : std::uint8_t
{
    Ok = 0,
    Overflow, // more planes or vertices than the capacity
    NeedsHost, // floating point could not decide which vertices a point removes (the host clips the cell exactly)
    NeedsHostVertex // floating point could not place a vertex (the host clips the cell exactly)
};

//! \brief A vertex: the cell's planes (slots) that meet there, counterclockwise seen from outside the cell
struct ClipTriangle
{
    std::uint8_t plane[3];
};

//! \brief Planes every cell is clipped by: a point on each and its unit normal, pointing out of the domain
struct ClipWalls
{
    int count = 0;
    FlatPoint point[CLIP_MAX_WALLS];
    FlatPoint normal[CLIP_MAX_WALLS];
};

//! \brief The cells being clipped (fixed capacity per cell), cells [firstCell, firstCell + cellCount) of `points`
struct ClipCellsView
{
    const FlatPoint *points = nullptr; // the cells' points, then the ghosts
    flat_index_t firstCell = 0;
    flat_index_t cellCount = 0;
    FlatPoint cubeLow{0.0, 0.0, 0.0}, cubeHigh{0.0, 0.0, 0.0};
    ClipWalls walls;
    flat_index_t *planes = nullptr; // [cell * CLIP_MAX_PLANES + slot]: the point of a bisector, a wall or a cube face
    std::uint8_t *planeCount = nullptr;
    ClipTriangle *triangles = nullptr; // [cell * CLIP_MAX_VERTICES + vertex]
    FlatPoint *vertices = nullptr; // [cell * CLIP_MAX_VERTICES + vertex], relative to the cell's point
    std::uint8_t *triangleCount = nullptr;
    double *radius2 = nullptr; // the squared distance from the cell's point to its farthest vertex
    double *known = nullptr; // the cell was clipped by every point (given before the last update) this close
    ClipStatus *status = nullptr;
};

//! \brief Points sorted into a uniform grid
struct ClipGridView
{
    const flat_index_t *cellStart = nullptr; // [grid cell .. grid cell + 1): the grid cell's points in `points`
    const flat_index_t *points = nullptr;
    FlatPoint low{0.0, 0.0, 0.0};
    double size = 1;
    int dims[3] = {0, 0, 0};
};

namespace clip_detail
{
    constexpr double ROUND_OFF = 1.1102230246251565e-16;
    constexpr double ORIENT_BOUND = (7.0 + 56.0 * ROUND_OFF) * ROUND_OFF;
    constexpr double INSPHERE_BOUND = (16.0 + 224.0 * ROUND_OFF) * ROUND_OFF;

    MADVORO_INLINE_FUNCTION double Component(const FlatPoint &p, int axis)
    {
        return (axis == 0) ? p.x : ((axis == 1) ? p.y : p.z);
    }

    MADVORO_INLINE_FUNCTION double Min(double a, double b)
    {
        return (b < a) ? b : a;
    }

    MADVORO_INLINE_FUNCTION int Clamp(int value, int low, int high)
    {
        return (value < low) ? low : ((value > high) ? high : value);
    }

    // the sign of Shewchuk's orient3d where the floating point value is certain, else 0
    MADVORO_INLINE_FUNCTION int Orient3dSign(const FlatPoint &a, const FlatPoint &b, const FlatPoint &c, const FlatPoint &d)
    {
        const double adx = a.x - d.x, bdx = b.x - d.x, cdx = c.x - d.x;
        const double ady = a.y - d.y, bdy = b.y - d.y, cdy = c.y - d.y;
        const double adz = a.z - d.z, bdz = b.z - d.z, cdz = c.z - d.z;
        const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy, cdxady = cdx * ady, adxcdy = adx * cdy, adxbdy = adx * bdy, bdxady = bdx * ady;
        const double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
        const double permanent = (Portable::Abs(bdxcdy) + Portable::Abs(cdxbdy)) * Portable::Abs(adz) + (Portable::Abs(cdxady) + Portable::Abs(adxcdy)) * Portable::Abs(bdz) + (Portable::Abs(adxbdy) + Portable::Abs(bdxady)) * Portable::Abs(cdz);
        const double bound = ORIENT_BOUND * permanent;
        return (det > bound) ? 1 : ((-det > bound) ? -1 : 0);
    }

    // the sign of Shewchuk's insphere where the floating point value is certain, else 0
    MADVORO_INLINE_FUNCTION int InsphereSign(const FlatPoint &a, const FlatPoint &b, const FlatPoint &c, const FlatPoint &d, const FlatPoint &e)
    {
        const double aex = a.x - e.x, bex = b.x - e.x, cex = c.x - e.x, dex = d.x - e.x;
        const double aey = a.y - e.y, bey = b.y - e.y, cey = c.y - e.y, dey = d.y - e.y;
        const double aez = a.z - e.z, bez = b.z - e.z, cez = c.z - e.z, dez = d.z - e.z;
        const double aexbey = aex * bey, bexaey = bex * aey, bexcey = bex * cey, cexbey = cex * bey;
        const double cexdey = cex * dey, dexcey = dex * cey, dexaey = dex * aey, aexdey = aex * dey;
        const double aexcey = aex * cey, cexaey = cex * aey, bexdey = bex * dey, dexbey = dex * bey;
        const double ab = aexbey - bexaey, bc = bexcey - cexbey, cd = cexdey - dexcey, da = dexaey - aexdey, ac = aexcey - cexaey, bd = bexdey - dexbey;
        const double abc = aez * bc - bez * ac + cez * ab;
        const double bcd = bez * cd - cez * bd + dez * bc;
        const double cda = cez * da + dez * ac + aez * cd;
        const double dab = dez * ab + aez * bd + bez * da;
        const double alift = aex * aex + aey * aey + aez * aez;
        const double blift = bex * bex + bey * bey + bez * bez;
        const double clift = cex * cex + cey * cey + cez * cez;
        const double dlift = dex * dex + dey * dey + dez * dez;
        const double det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
        const double permanent = ((Portable::Abs(cexdey) + Portable::Abs(dexcey)) * Portable::Abs(bez) + (Portable::Abs(dexbey) + Portable::Abs(bexdey)) * Portable::Abs(cez) + (Portable::Abs(bexcey) + Portable::Abs(cexbey)) * Portable::Abs(dez)) * alift
            + ((Portable::Abs(dexaey) + Portable::Abs(aexdey)) * Portable::Abs(cez) + (Portable::Abs(aexcey) + Portable::Abs(cexaey)) * Portable::Abs(dez) + (Portable::Abs(cexdey) + Portable::Abs(dexcey)) * Portable::Abs(aez)) * blift
            + ((Portable::Abs(aexbey) + Portable::Abs(bexaey)) * Portable::Abs(dez) + (Portable::Abs(bexdey) + Portable::Abs(dexbey)) * Portable::Abs(aez) + (Portable::Abs(dexaey) + Portable::Abs(aexdey)) * Portable::Abs(bez)) * clift
            + ((Portable::Abs(bexcey) + Portable::Abs(cexbey)) * Portable::Abs(aez) + (Portable::Abs(cexaey) + Portable::Abs(aexcey)) * Portable::Abs(bez) + (Portable::Abs(aexbey) + Portable::Abs(bexaey)) * Portable::Abs(cez)) * dlift;
        const double bound = INSPHERE_BOUND * permanent;
        return (det > bound) ? 1 : ((-det > bound) ? -1 : 0);
    }

    MADVORO_INLINE_FUNCTION int Sign(double value)
    {
        return (value > 0) ? 1 : ((value < 0) ? -1 : 0);
    }

    // a double-double: hi + lo, |lo| at most half an ulp of hi
    struct DoubleDouble
    {
        double hi;
        double lo;
    };

    MADVORO_INLINE_FUNCTION DoubleDouble TwoSum(double a, double b)
    {
        const double sum = a + b;
        const double bVirtual = sum - a;
        return DoubleDouble{sum, (a - (sum - bVirtual)) + (b - bVirtual)};
    }

    MADVORO_INLINE_FUNCTION DoubleDouble FastTwoSum(double a, double b)
    {
        const double sum = a + b;
        return DoubleDouble{sum, b - (sum - a)};
    }

    MADVORO_INLINE_FUNCTION DoubleDouble Exact(double a)
    {
        return DoubleDouble{a, 0.0};
    }

    // a - b exactly
    MADVORO_INLINE_FUNCTION DoubleDouble ExactDifference(double a, double b)
    {
        return TwoSum(a, -b);
    }

    MADVORO_INLINE_FUNCTION DoubleDouble operator+(const DoubleDouble &a, const DoubleDouble &b)
    {
        DoubleDouble s = TwoSum(a.hi, b.hi);
        const DoubleDouble t = TwoSum(a.lo, b.lo);
        s = FastTwoSum(s.hi, s.lo + t.hi);
        return FastTwoSum(s.hi, s.lo + t.lo);
    }

    MADVORO_INLINE_FUNCTION DoubleDouble operator-(const DoubleDouble &a)
    {
        return DoubleDouble{-a.hi, -a.lo};
    }

    MADVORO_INLINE_FUNCTION DoubleDouble operator-(const DoubleDouble &a, const DoubleDouble &b)
    {
        return a + (-b);
    }

    MADVORO_INLINE_FUNCTION DoubleDouble operator*(const DoubleDouble &a, const DoubleDouble &b)
    {
        const double product = a.hi * b.hi;
        const double error = Portable::Fma(a.hi, b.hi, -product);
        return FastTwoSum(product, error + (a.hi * b.lo + a.lo * b.hi));
    }

    // where the bisectors of p with the 3 points meet (relative to p), in double-double from their exact offsets; false if the planes do not meet in a point that double-double can certify
    MADVORO_INLINE_FUNCTION bool BisectorsMeetDD(const FlatPoint &p, const FlatPoint *points, FlatPoint &vertex)
    {
        DoubleDouble rows[3][3];
        DoubleDouble values[3];
        double scale = 1;
        for(int i = 0; i < 3; ++i)
        {
            rows[i][0] = ExactDifference(points[i].x, p.x);
            rows[i][1] = ExactDifference(points[i].y, p.y);
            rows[i][2] = ExactDifference(points[i].z, p.z);
            values[i] = rows[i][0] * rows[i][0] + rows[i][1] * rows[i][1] + rows[i][2] * rows[i][2];
            values[i] = DoubleDouble{0.5 * values[i].hi, 0.5 * values[i].lo};
            scale *= Portable::Sqrt(values[i].hi * 2);
        }
        const auto cross = [](const DoubleDouble *a, const DoubleDouble *b, DoubleDouble *out)
        {
            out[0] = a[1] * b[2] - a[2] * b[1];
            out[1] = a[2] * b[0] - a[0] * b[2];
            out[2] = a[0] * b[1] - a[1] * b[0];
        };
        DoubleDouble c12[3], c20[3], c01[3];
        cross(rows[1], rows[2], c12);
        cross(rows[2], rows[0], c20);
        cross(rows[0], rows[1], c01);
        const DoubleDouble det = rows[0][0] * c12[0] + rows[0][1] * c12[1] + rows[0][2] * c12[2];
        if(not (Portable::Abs(det.hi) > 1e-28 * scale))
        {
            return false;
        }
        double coordinates[3];
        for(int k = 0; k < 3; ++k)
        {
            const DoubleDouble numerator = values[0] * c12[k] + values[1] * c20[k] + values[2] * c01[k];
            coordinates[k] = numerator.hi / det.hi;
        }
        vertex = FlatPoint{coordinates[0], coordinates[1], coordinates[2]};
        return true;
    }

    // the relative error of a double-double operation is a small multiple of ROUND_OFF^2; the bounds allow 256 of them on the permanent
    constexpr double DOUBLE_DOUBLE_BOUND = 256.0 * ROUND_OFF * ROUND_OFF;

    // the sign of orient3d in double-double where it is certain, else 0
    MADVORO_INLINE_FUNCTION int Orient3dSignDD(const FlatPoint &a, const FlatPoint &b, const FlatPoint &c, const FlatPoint &d)
    {
        const DoubleDouble adx = ExactDifference(a.x, d.x), bdx = ExactDifference(b.x, d.x), cdx = ExactDifference(c.x, d.x);
        const DoubleDouble ady = ExactDifference(a.y, d.y), bdy = ExactDifference(b.y, d.y), cdy = ExactDifference(c.y, d.y);
        const DoubleDouble adz = ExactDifference(a.z, d.z), bdz = ExactDifference(b.z, d.z), cdz = ExactDifference(c.z, d.z);
        const DoubleDouble det = adz * (bdx * cdy - cdx * bdy) + bdz * (cdx * ady - adx * cdy) + cdz * (adx * bdy - bdx * ady);
        const double permanent = (Portable::Abs(bdx.hi * cdy.hi) + Portable::Abs(cdx.hi * bdy.hi)) * Portable::Abs(adz.hi) + (Portable::Abs(cdx.hi * ady.hi) + Portable::Abs(adx.hi * cdy.hi)) * Portable::Abs(bdz.hi) + (Portable::Abs(adx.hi * bdy.hi) + Portable::Abs(bdx.hi * ady.hi)) * Portable::Abs(cdz.hi);
        const double bound = DOUBLE_DOUBLE_BOUND * permanent * 1.01;
        const double value = det.hi + det.lo;
        return (value > bound) ? 1 : ((-value > bound) ? -1 : 0);
    }

    // the sign of insphere in double-double where it is certain, else 0
    MADVORO_INLINE_FUNCTION int InsphereSignDD(const FlatPoint &a, const FlatPoint &b, const FlatPoint &c, const FlatPoint &d, const FlatPoint &e)
    {
        const DoubleDouble aex = ExactDifference(a.x, e.x), bex = ExactDifference(b.x, e.x), cex = ExactDifference(c.x, e.x), dex = ExactDifference(d.x, e.x);
        const DoubleDouble aey = ExactDifference(a.y, e.y), bey = ExactDifference(b.y, e.y), cey = ExactDifference(c.y, e.y), dey = ExactDifference(d.y, e.y);
        const DoubleDouble aez = ExactDifference(a.z, e.z), bez = ExactDifference(b.z, e.z), cez = ExactDifference(c.z, e.z), dez = ExactDifference(d.z, e.z);
        const DoubleDouble ab = aex * bey - bex * aey, bc = bex * cey - cex * bey, cd = cex * dey - dex * cey, da = dex * aey - aex * dey, ac = aex * cey - cex * aey, bd = bex * dey - dex * bey;
        const DoubleDouble abc = aez * bc - bez * ac + cez * ab;
        const DoubleDouble bcd = bez * cd - cez * bd + dez * bc;
        const DoubleDouble cda = cez * da + dez * ac + aez * cd;
        const DoubleDouble dab = dez * ab + aez * bd + bez * da;
        const DoubleDouble alift = aex * aex + aey * aey + aez * aez;
        const DoubleDouble blift = bex * bex + bey * bey + bez * bez;
        const DoubleDouble clift = cex * cex + cey * cey + cez * cez;
        const DoubleDouble dlift = dex * dex + dey * dey + dez * dez;
        const DoubleDouble det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
        const double ax = Portable::Abs(aex.hi), ay = Portable::Abs(aey.hi), az = Portable::Abs(aez.hi), bx = Portable::Abs(bex.hi), by = Portable::Abs(bey.hi), bz = Portable::Abs(bez.hi);
        const double cx = Portable::Abs(cex.hi), cy = Portable::Abs(cey.hi), cz = Portable::Abs(cez.hi), dx = Portable::Abs(dex.hi), dy = Portable::Abs(dey.hi), dz = Portable::Abs(dez.hi);
        const double permanent = ((cx * dy + dx * cy) * bz + (dx * by + bx * dy) * cz + (bx * cy + cx * by) * dz) * (ax * ax + ay * ay + az * az)
            + ((dx * ay + ax * dy) * cz + (ax * cy + cx * ay) * dz + (cx * dy + dx * cy) * az) * (bx * bx + by * by + bz * bz)
            + ((ax * by + bx * ay) * dz + (bx * dy + dx * by) * az + (dx * ay + ax * dy) * bz) * (cx * cx + cy * cy + cz * cz)
            + ((bx * cy + cx * by) * az + (cx * ay + ax * cy) * bz + (ax * by + bx * ay) * cz) * (dx * dx + dy * dy + dz * dz);
        const double bound = DOUBLE_DOUBLE_BOUND * permanent * 1.01;
        const double value = det.hi + det.lo;
        return (value > bound) ? 1 : ((-value > bound) ? -1 : 0);
    }

    /** \brief Whether point q is strictly closer than p to the vertex where the bisectors of p with a, b and c meet: q strictly inside the sphere through p, a, b and c. 0 or 1, or -1 when `Arithmetic` cannot decide */
    template<int Arithmetic>
    MADVORO_INLINE_FUNCTION int CloserToVertex(const FlatPoint &a, const FlatPoint &b, const FlatPoint &c, const FlatPoint &p, const FlatPoint &q)
    {
        int orientation = Orient3dSign(a, b, c, p);
        int inside = InsphereSign(a, b, c, p, q);
        // where double precision cannot decide, double-double usually can (ties that only round-off breaks)
        if(orientation == 0)
        {
            orientation = Orient3dSignDD(a, b, c, p);
        }
        if(inside == 0)
        {
            inside = InsphereSignDD(a, b, c, p, q);
        }
        if constexpr(Arithmetic == CLIP_EXPANSIONS)
        {
            if(orientation == 0)
            {
                orientation = ExpansionOrient3dSign(a, b, c, p);
            }
            if(inside == 0)
            {
                inside = ExpansionInsphereSign(a, b, c, p, q);
                if(inside == EXPANSION_OVERFLOW)
                {
                    return -1;
                }
            }
            // (a flat tetrahedron has no sphere: such a vertex is not removed)
            return (orientation * inside > 0) ? 1 : 0;
        }
        else if constexpr(Arithmetic == CLIP_EXACT)
        {
            if(orientation == 0)
            {
                orientation = Sign(MadVoro::orient3d(std::array<FlatPoint, 4>{a, b, c, p}));
            }
            if(inside == 0)
            {
                inside = Sign(MadVoro::insphere(std::array<FlatPoint, 5>{a, b, c, p, q}));
            }
            // (a flat tetrahedron has no sphere: such a vertex is not removed)
            return (orientation * inside > 0) ? 1 : 0;
        }
        else
        {
            if(orientation == 0)
            {
                return -1;
            }
            if(inside == 0)
            {
                // on the sphere to round-off
                return -1;
            }
            return (orientation * inside > 0) ? 1 : 0;
        }
    }

    // the point whose bisector with the cell's point the plane is (for a wall or a cube face: the point's mirror across it)
    MADVORO_INLINE_FUNCTION FlatPoint PointOf(const ClipCellsView &cells, flat_index_t cell, flat_index_t plane)
    {
        if(plane < CLIP_WALL_PLANE)
        {
            return cells.points[plane];
        }
        const FlatPoint &p = cells.points[cell];
        if(plane < CLIP_CUBE_PLANE)
        {
            const int wall = static_cast<int>(plane - CLIP_WALL_PLANE);
            const FlatPoint &normal = cells.walls.normal[wall];
            const double twice = 2 * Dot(p - cells.walls.point[wall], normal);
            return FlatPoint{p.x - twice * normal.x, p.y - twice * normal.y, p.z - twice * normal.z};
        }
        const int face = static_cast<int>(plane - CLIP_CUBE_PLANE);
        const int axis = face / 2;
        const double offset = ((face % 2) == 1) ? 2 * (Component(cells.cubeHigh, axis) - Component(p, axis)) : -2 * (Component(p, axis) - Component(cells.cubeLow, axis));
        return FlatPoint{p.x + ((axis == 0) ? offset : 0.0), p.y + ((axis == 1) ? offset : 0.0), p.z + ((axis == 2) ? offset : 0.0)};
    }

    // where the planes meet with a wall or cube face: on the first point's bisectors with the other points, on the walls and on the cube faces (false if they hardly meet)
    MADVORO_INLINE_FUNCTION bool PlanesVertex(const ClipCellsView &cells, const flat_index_t order[4], FlatPoint &center)
    {
        const FlatPoint &first = cells.points[order[0]];
        FlatPoint rows[3];
        double values[3];
        for(int i = 0; i < 3; ++i)
        {
            const flat_index_t plane = order[i + 1];
            if(plane < CLIP_WALL_PLANE)
            {
                rows[i] = cells.points[plane] - first;
                values[i] = 0.5 * Dot(rows[i], rows[i]);
            }
            else if(plane < CLIP_CUBE_PLANE)
            {
                const int wall = static_cast<int>(plane - CLIP_WALL_PLANE);
                rows[i] = cells.walls.normal[wall];
                values[i] = Dot(cells.walls.normal[wall], cells.walls.point[wall] - first);
            }
            else
            {
                const int face = static_cast<int>(plane - CLIP_CUBE_PLANE);
                const int axis = face / 2;
                rows[i] = FlatPoint{(axis == 0) ? 1.0 : 0.0, (axis == 1) ? 1.0 : 0.0, (axis == 2) ? 1.0 : 0.0};
                values[i] = Component(((face % 2) == 1) ? cells.cubeHigh : cells.cubeLow, axis) - Component(first, axis);
            }
        }
        const FlatPoint c12 = Cross(rows[1], rows[2]), c20 = Cross(rows[2], rows[0]), c01 = Cross(rows[0], rows[1]);
        const double det = Dot(rows[0], c12);
        const double scale = Portable::Sqrt(Dot(rows[0], rows[0]) * Dot(rows[1], rows[1]) * Dot(rows[2], rows[2]));
        if(not (Portable::Abs(det) > 1e-8 * scale))
        {
            return false;
        }
        const double inverse = 1.0 / det;
        center = FlatPoint{first.x + (values[0] * c12.x + values[1] * c20.x + values[2] * c01.x) * inverse, first.y + (values[0] * c12.y + values[1] * c20.y + values[2] * c01.y) * inverse, first.z + (values[0] * c12.z + values[1] * c20.z + values[2] * c01.z) * inverse};
        return true;
    }

    /** \brief The vertex of the 4 points (indices in increasing order: a mirror is of the first point) where the bisectors of the first with the others meet, the same in each cell that has it. False if floating point cannot place it and `Exact` is false */
    template<bool Exact>
    MADVORO_INLINE_FUNCTION bool CanonicalVertex(const ClipCellsView &cells, const flat_index_t order[4], FlatPoint &center)
    {
        // with walls or cube faces, from their planes themselves
        if(order[3] >= CLIP_WALL_PLANE and PlanesVertex(cells, order, center))
        {
            return true;
        }
        FlatPoint points[4];
        for(int i = 0; i < 4; ++i)
        {
            points[i] = PointOf(cells, order[0], order[i]);
        }
        double radius = 0;
        if(Circumcenter(points[0], points[1], points[2], points[3], center, radius) == GeometryStatus::Ok)
        {
            return true;
        }
        if constexpr(Exact)
        {
            HiPrecisionCircumcenter(std::array<FlatPoint, 4>{points[0], points[1], points[2], points[3]}, center);
            return true;
        }
        else
        {
            return false;
        }
    }

    /** \brief Where the bisectors of the cell's point with the points of planes a, b and c meet (relative to the cell's point): on the planes the cell's predicates decide by (for a wall or a cube face, the bisector with the point's own mirror), so that each position agrees with the decisions even where the planes hardly meet. False if floating point cannot place it and `Arithmetic` is not exact */
    template<int Arithmetic>
    MADVORO_INLINE_FUNCTION bool VertexOf(const ClipCellsView &cells, flat_index_t cell, flat_index_t a, flat_index_t b, flat_index_t c, FlatPoint &vertex)
    {
        const FlatPoint &p = cells.points[cell];
        const flat_index_t planes[3] = {a, b, c};
        FlatPoint rows[3];
        double values[3];
        for(int i = 0; i < 3; ++i)
        {
            rows[i] = PointOf(cells, cell, planes[i]) - p;
            values[i] = 0.5 * Dot(rows[i], rows[i]);
        }
        const FlatPoint c12 = Cross(rows[1], rows[2]), c20 = Cross(rows[2], rows[0]), c01 = Cross(rows[0], rows[1]);
        const double det = Dot(rows[0], c12);
        const double scale = Portable::Sqrt(Dot(rows[0], rows[0]) * Dot(rows[1], rows[1]) * Dot(rows[2], rows[2]));
        if(Portable::Abs(det) > 1e-8 * scale)
        {
            const double inverse = 1.0 / det;
            vertex = FlatPoint{(values[0] * c12.x + values[1] * c20.x + values[2] * c01.x) * inverse, (values[0] * c12.y + values[1] * c20.y + values[2] * c01.y) * inverse, (values[0] * c12.z + values[1] * c20.z + values[2] * c01.z) * inverse};
        }
        else
        {
            if constexpr(Arithmetic == CLIP_EXACT)
            {
                FlatPoint center{0.0, 0.0, 0.0};
                HiPrecisionCircumcenter(std::array<FlatPoint, 4>{p, PointOf(cells, cell, a), PointOf(cells, cell, b), PointOf(cells, cell, c)}, center);
                vertex = center - p;
            }
            else if constexpr(Arithmetic == CLIP_EXPANSIONS)
            {
                const FlatPoint points[3] = {PointOf(cells, cell, a), PointOf(cells, cell, b), PointOf(cells, cell, c)};
                if(not BisectorsMeetDD(p, points, vertex))
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        return Portable::Abs(vertex.x) < CLIP_INFINITE and Portable::Abs(vertex.y) < CLIP_INFINITE and Portable::Abs(vertex.z) < CLIP_INFINITE;
    }

    MADVORO_INLINE_FUNCTION void UpdateRadius(const ClipCellsView &cells, flat_index_t cell)
    {
        const flat_index_t local = cell - cells.firstCell;
        const FlatPoint *vertices = cells.vertices + local * CLIP_MAX_VERTICES;
        double radius2 = 0;
        for(int t = 0; t < cells.triangleCount[local]; ++t)
        {
            radius2 = Portable::Max(radius2, Dot(vertices[t], vertices[t]));
        }
        cells.radius2[local] = radius2;
    }

    // frees the plane slots no vertex uses; false if none could be freed
    MADVORO_INLINE_FUNCTION bool CompactPlanes(const ClipCellsView &cells, flat_index_t cell)
    {
        const flat_index_t local = cell - cells.firstCell;
        flat_index_t *planes = cells.planes + local * CLIP_MAX_PLANES;
        ClipTriangle *triangles = cells.triangles + local * CLIP_MAX_VERTICES;
        std::uint8_t used[CLIP_MAX_PLANES];
        for(int k = 0; k < cells.planeCount[local]; ++k)
        {
            used[k] = 0;
        }
        for(int t = 0; t < cells.triangleCount[local]; ++t)
        {
            for(int k = 0; k < 3; ++k)
            {
                used[triangles[t].plane[k]] = 1;
            }
        }
        std::uint8_t remap[CLIP_MAX_PLANES];
        int count = 0;
        for(int k = 0; k < cells.planeCount[local]; ++k)
        {
            if(used[k] != 0)
            {
                remap[k] = static_cast<std::uint8_t>(count);
                planes[count] = planes[k];
                ++count;
            }
        }
        for(int t = 0; t < cells.triangleCount[local]; ++t)
        {
            for(int k = 0; k < 3; ++k)
            {
                triangles[t].plane[k] = remap[triangles[t].plane[k]];
            }
        }
        const bool freed = count < cells.planeCount[local];
        cells.planeCount[local] = static_cast<std::uint8_t>(count);
        return freed;
    }
} // namespace clip_detail

//! \brief The cube, with its 8 vertices
MADVORO_INLINE_FUNCTION void InitClipCell(const ClipCellsView &cells, flat_index_t cell)
{
    const flat_index_t local = cell - cells.firstCell;
    flat_index_t *planes = cells.planes + local * CLIP_MAX_PLANES;
    ClipTriangle *triangles = cells.triangles + local * CLIP_MAX_VERTICES;
    FlatPoint *vertices = cells.vertices + local * CLIP_MAX_VERTICES;
    for(int face = 0; face < 6; ++face)
    {
        planes[face] = CLIP_CUBE_PLANE + static_cast<flat_index_t>(face);
    }
    cells.planeCount[local] = 6;
    const FlatPoint &p = cells.points[cell];
    int t = 0;
    for(int sx = 0; sx < 2; ++sx)
    {
        for(int sy = 0; sy < 2; ++sy)
        {
            for(int sz = 0; sz < 2; ++sz)
            {
                const std::uint8_t a = static_cast<std::uint8_t>(sx), b = static_cast<std::uint8_t>(2 + sy), c = static_cast<std::uint8_t>(4 + sz);
                // counterclockwise from outside: the normals' triple product is positive
                const bool positive = ((sx + sy + sz) % 2) == 1;
                triangles[t].plane[0] = a;
                triangles[t].plane[1] = positive ? b : c;
                triangles[t].plane[2] = positive ? c : b;
                vertices[t] = FlatPoint{((sx == 1) ? cells.cubeHigh.x : cells.cubeLow.x) - p.x, ((sy == 1) ? cells.cubeHigh.y : cells.cubeLow.y) - p.y, ((sz == 1) ? cells.cubeHigh.z : cells.cubeLow.z) - p.z};
                ++t;
            }
        }
    }
    cells.triangleCount[local] = 8;
    cells.known[local] = 0;
    cells.status[local] = ClipStatus::Ok;
    clip_detail::UpdateRadius(cells, cell);
}

/** \brief Clips a cell by the bisector of its point and point `q`: removes the vertices q is strictly closer to than the cell's point. Decided on the points themselves, in `Arithmetic` where floating point cannot (a cell it cannot decide NeedsHost) */
template<int Arithmetic>
MADVORO_INLINE_FUNCTION void ClipCellBy(const ClipCellsView &cells, flat_index_t cell, flat_index_t q)
{
    const FlatPoint &p = cells.points[cell];
    const FlatPoint qPoint = clip_detail::PointOf(cells, cell, q);
    const FlatPoint offset = qPoint - p;
    if(not (Dot(offset, offset) > 0))
    {
        // a copy of the cell's point
        return;
    }
    const flat_index_t local = cell - cells.firstCell;
    flat_index_t *planes = cells.planes + local * CLIP_MAX_PLANES;
    ClipTriangle *triangles = cells.triangles + local * CLIP_MAX_VERTICES;
    FlatPoint *vertices = cells.vertices + local * CLIP_MAX_VERTICES;
    const int count = cells.triangleCount[local];
    // a point the cell already has a plane of removes nothing
    for(int k = 0; k < cells.planeCount[local]; ++k)
    {
        if(planes[k] == q)
        {
            return;
        }
    }
    // a vertex's position decides outside a margin wider than its round-off (to 1e-6 of its distance from the point, scaled by q's)
    const double qDistance = Portable::Sqrt(Dot(offset, offset));
    std::uint8_t removed[CLIP_MAX_VERTICES];
    int removedCount = 0;
    int undecided = 0;
    for(int t = 0; t < count; ++t)
    {
        const FlatPoint fromQ = vertices[t] - offset;
        const double farther = Dot(fromQ, fromQ) - Dot(vertices[t], vertices[t]);
        const double margin = 1e-5 * qDistance * Portable::Sqrt(Dot(vertices[t], vertices[t]));
        removed[t] = (farther < -margin) ? 1 : ((farther > margin) ? 0 : 2);
        removedCount += (removed[t] == 1) ? 1 : 0;
        undecided += (removed[t] == 2) ? 1 : 0;
    }
    if(removedCount == 0 and undecided == 0)
    {
        return;
    }
    for(int t = 0; t < count; ++t)
    {
        if(removed[t] != 2)
        {
            continue;
        }
        const FlatPoint a = clip_detail::PointOf(cells, cell, planes[triangles[t].plane[0]]);
        const FlatPoint b = clip_detail::PointOf(cells, cell, planes[triangles[t].plane[1]]);
        const FlatPoint c = clip_detail::PointOf(cells, cell, planes[triangles[t].plane[2]]);
        const int closer = clip_detail::CloserToVertex<Arithmetic>(a, b, c, p, qPoint);
        if(closer < 0)
        {
            cells.status[local] = ClipStatus::NeedsHost;
            return;
        }
        removed[t] = static_cast<std::uint8_t>(closer);
        removedCount += removed[t];
    }
    if(removedCount == 0)
    {
        return;
    }
    if(cells.planeCount[local] >= CLIP_MAX_PLANES and not clip_detail::CompactPlanes(cells, cell))
    {
        cells.status[local] = ClipStatus::Overflow;
        return;
    }
    // the edges of the removed vertices whose other vertex stays: the new face's boundary
    std::uint8_t edgeFrom[3 * CLIP_MAX_VERTICES], edgeTo[3 * CLIP_MAX_VERTICES];
    int edgeCount = 0;
    for(int t = 0; t < count; ++t)
    {
        if(removed[t] == 0)
        {
            continue;
        }
        for(int k = 0; k < 3; ++k)
        {
            edgeFrom[edgeCount] = triangles[t].plane[k];
            edgeTo[edgeCount] = triangles[t].plane[(k + 1) % 3];
            ++edgeCount;
        }
    }
    std::uint8_t boundaryFrom[3 * CLIP_MAX_VERTICES], boundaryTo[3 * CLIP_MAX_VERTICES];
    int boundaryCount = 0;
    for(int e = 0; e < edgeCount; ++e)
    {
        bool inner = false;
        for(int f = 0; f < edgeCount and not inner; ++f)
        {
            inner = (edgeFrom[f] == edgeTo[e]) and (edgeTo[f] == edgeFrom[e]);
        }
        if(not inner)
        {
            boundaryFrom[boundaryCount] = edgeFrom[e];
            boundaryTo[boundaryCount] = edgeTo[e];
            ++boundaryCount;
        }
    }
    if(count - removedCount + boundaryCount > CLIP_MAX_VERTICES)
    {
        cells.status[local] = ClipStatus::Overflow;
        return;
    }
    // the new vertices, where the boundary's edges meet the new plane
    FlatPoint created[CLIP_MAX_VERTICES];
    for(int e = 0; e < boundaryCount; ++e)
    {
        if(not clip_detail::VertexOf<Arithmetic>(cells, cell, planes[boundaryFrom[e]], planes[boundaryTo[e]], q, created[e]))
        {
            cells.status[local] = ClipStatus::NeedsHostVertex;
            return;
        }
    }
    const std::uint8_t slot = cells.planeCount[local];
    planes[slot] = q;
    cells.planeCount[local] = static_cast<std::uint8_t>(slot + 1);
    int kept = 0;
    for(int t = 0; t < count; ++t)
    {
        if(removed[t] == 0)
        {
            triangles[kept] = triangles[t];
            vertices[kept] = vertices[t];
            ++kept;
        }
    }
    for(int e = 0; e < boundaryCount; ++e)
    {
        triangles[kept].plane[0] = boundaryFrom[e];
        triangles[kept].plane[1] = boundaryTo[e];
        triangles[kept].plane[2] = slot;
        vertices[kept] = created[e];
        ++kept;
    }
    cells.triangleCount[local] = static_cast<std::uint8_t>(kept);
    clip_detail::UpdateRadius(cells, cell);
}

//! \brief The cube clipped by the walls
template<int Arithmetic>
MADVORO_INLINE_FUNCTION void StartClipCell(const ClipCellsView &cells, flat_index_t cell)
{
    InitClipCell(cells, cell);
    for(int wall = 0; wall < cells.walls.count; ++wall)
    {
        ClipCellBy<Arithmetic>(cells, cell, CLIP_WALL_PLANE + static_cast<flat_index_t>(wall));
    }
}

namespace clip_detail
{
    // visits the grid cells in shells around the point's until the shells cover the cell's 2R (or the grid); returns the distance they cover. `NewOnly`: every point of the grid is new
    template<bool NewOnly, int Arithmetic>
    MADVORO_INLINE_FUNCTION double ClipByGrid(const ClipCellsView &cells, flat_index_t cell, const ClipGridView &grid, flat_index_t firstNew, double known)
    {
        const flat_index_t local = cell - cells.firstCell;
        const FlatPoint p = cells.points[cell];
        int g[3];
        for(int a = 0; a < 3; ++a)
        {
            g[a] = Clamp(static_cast<int>((Component(p, a) - Component(grid.low, a)) / grid.size), 0, grid.dims[a] - 1);
        }
        const double known2 = known * known;
        for(int s = 0;; ++s)
        {
            const int from0 = Clamp(g[0] - s, 0, grid.dims[0] - 1), to0 = Clamp(g[0] + s, 0, grid.dims[0] - 1);
            const int from1 = Clamp(g[1] - s, 0, grid.dims[1] - 1), to1 = Clamp(g[1] + s, 0, grid.dims[1] - 1);
            for(int i = from0; i <= to0; ++i)
            {
                for(int j = from1; j <= to1; ++j)
                {
                    const bool side = (i == g[0] - s) or (i == g[0] + s) or (j == g[1] - s) or (j == g[1] + s);
                    // the shell's cells in this column: all of them on its sides, else its two ends
                    for(int end = 0; end < (side ? (2 * s + 1) : 2); ++end)
                    {
                        const int k = side ? (g[2] - s + end) : ((end == 0) ? (g[2] - s) : (g[2] + s));
                        if(k < 0 or k >= grid.dims[2] or (s == 0 and end > 0))
                        {
                            continue;
                        }
                        const flat_index_t gridCell = (static_cast<flat_index_t>(i) * grid.dims[1] + j) * grid.dims[2] + k;
                        for(flat_index_t m = grid.cellStart[gridCell]; m < grid.cellStart[gridCell + 1]; ++m)
                        {
                            const flat_index_t q = grid.points[m];
                            if(q == cell or cells.status[local] != ClipStatus::Ok)
                            {
                                continue;
                            }
                            const FlatPoint offset = cells.points[q] - p;
                            const double distance2 = Dot(offset, offset);
                            // a point beyond 2R cannot cut the cell; one within `known` (given before) already did
                            if(distance2 < 4 * cells.radius2[local] and (NewOnly or q >= firstNew or distance2 > known2))
                            {
                                ClipCellBy<Arithmetic>(cells, cell, q);
                            }
                        }
                    }
                }
            }
            // the distance from the point to the region the shells cover (points lie only inside the grid)
            double covered = CLIP_INFINITE;
            for(int a = 0; a < 3; ++a)
            {
                if(g[a] - s > 0)
                {
                    covered = Min(covered, Component(p, a) - (Component(grid.low, a) + (g[a] - s) * grid.size));
                }
                if(g[a] + s < grid.dims[a] - 1)
                {
                    covered = Min(covered, Component(grid.low, a) + (g[a] + s + 1) * grid.size - Component(p, a));
                }
            }
            if(covered >= CLIP_INFINITE or covered * covered >= 4 * cells.radius2[local] or cells.status[local] != ClipStatus::Ok)
            {
                return covered;
            }
        }
    }
} // namespace clip_detail

/** \brief Clips a cell by the points it may not have been clipped by: those from `firstNew` on (also in `added`, if not null), and, unless the cell is secure (2R within its known distance), those of `all` beyond its known distance */
template<int Arithmetic>
MADVORO_INLINE_FUNCTION void UpdateClipCell(const ClipCellsView &cells, flat_index_t cell, const ClipGridView &all, const ClipGridView *added, flat_index_t firstNew)
{
    const flat_index_t local = cell - cells.firstCell;
    if(cells.status[local] != ClipStatus::Ok)
    {
        return;
    }
    const double known = cells.known[local];
    if(known >= CLIP_INFINITE or 4 * cells.radius2[local] <= known * known)
    {
        if(added != nullptr and added->dims[0] > 0)
        {
            clip_detail::ClipByGrid<true, Arithmetic>(cells, cell, *added, firstNew, known);
        }
        return;
    }
    const double covered = clip_detail::ClipByGrid<false, Arithmetic>(cells, cell, all, firstNew, known);
    cells.known[local] = Portable::Max(known, covered);
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_CLIP_CELL_HPP
