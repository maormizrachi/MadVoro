#ifndef MADVORO_DELAUNAY2D_HPP
#define MADVORO_DELAUNAY2D_HPP

#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include "../utils/Predicates3D.hpp"
#include "../exception/MadVoroException.hpp"

namespace MadVoro
{
    namespace Plane2D
    {
        constexpr double epsilon = 1.1102230246251565e-16; // 2^-53
        constexpr double ccwErrBoundA = (3.0 + 16.0 * epsilon) * epsilon;
        constexpr double iccErrBoundA = (10.0 + 96.0 * epsilon) * epsilon;

        /*! \brief Exact orientation of a, b, c in the plane: > 0 counterclockwise, < 0 clockwise, 0 collinear. Shewchuk's
          floating filter, then orient3d of the points lifted to z = 0 with (c.x, c.y, 1) as the fourth point, which as a
          polynomial is -orient2d, so the 3D exact arithmetic decides the uncertain cases.
         */
        template <typename PointT>
        double Orient(double ax, double ay, double bx, double by, double cx, double cy)
        {
            const double left = (ax - cx) * (by - cy), right = (ay - cy) * (bx - cx);
            const double det = left - right;
            const double errbound = ccwErrBoundA * (std::abs(left) + std::abs(right));
            if(det > errbound or -det > errbound)
            {
                return det;
            }
            const std::array<PointT, 4> points = {PointT(ax, ay, 0), PointT(bx, by, 0), PointT(cx, cy, 0), PointT(cx, cy, 1)};
            return -orient3d(points);
        }

        /*! \brief Exact in-circle test: > 0 if d is inside the circle through the counterclockwise a, b, c, < 0 outside,
          0 on it. Shewchuk's floating filter, then insphere of a, b, c, (c.x, c.y, 1), d lifted to z = 0: expanding that
          determinant along the z column leaves the in-circle determinant, with the sign reversed.
         */
        template <typename PointT>
        double InCircle(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy)
        {
            const double adx = ax - dx, ady = ay - dy, bdx = bx - dx, bdy = by - dy, cdx = cx - dx, cdy = cy - dy;
            const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy, cdxady = cdx * ady, adxcdy = adx * cdy, adxbdy = adx * bdy, bdxady = bdx * ady;
            const double alift = adx * adx + ady * ady, blift = bdx * bdx + bdy * bdy, clift = cdx * cdx + cdy * cdy;
            const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
            const double permanent = (std::abs(bdxcdy) + std::abs(cdxbdy)) * alift + (std::abs(cdxady) + std::abs(adxcdy)) * blift
                + (std::abs(adxbdy) + std::abs(bdxady)) * clift;
            const double errbound = iccErrBoundA * permanent;
            if(det > errbound or -det > errbound)
            {
                return det;
            }
            const std::array<PointT, 5> points = {PointT(ax, ay, 0), PointT(bx, by, 0), PointT(cx, cy, 0), PointT(cx, cy, 1), PointT(dx, dy, 0)};
            return -insphere(points);
        }
    } // namespace Plane2D

    /*! \brief Incremental Delaunay triangulation of points in a plane (Lawson flips, exact predicates), for the wedge's
      plane build (Voronoi3D::SetPlaneMesh). It starts from a large triangle around a box; its three corners are vertices
      0, 1, 2 (IsSuper) and every inserted point must lie inside it. Each vertex carries a label (the caller's index).
      Triangles are only split or flipped, never removed, so a triangle index stays valid; a triangle created or changed
      since the ghost search last queried it has checkBig set, as Delaunay3D does for tetrahedra.
     */
    template <typename PointT>
    class Delaunay2D
    {
    public:
        static constexpr size_t none = std::numeric_limits<size_t>::max();

        struct Triangle
        {
            std::array<size_t, 3> vertices;  //!< counterclockwise
            std::array<size_t, 3> neighbors; //!< neighbors[i] shares the edge opposite vertices[i]; none outside
            bool checkBig = true;
        };

        //! Clears the triangulation and starts from a triangle around [xMin, xMax] x [yMin, yMax]
        void Reset(double xMin, double xMax, double yMin, double yMax)
        {
            x_.clear();
            y_.clear();
            labels_.clear();
            triangles_.clear();
            vertexTriangle_.clear();
            const double cx = 0.5 * (xMin + xMax), cy = 0.5 * (yMin + yMax);
            const double size = std::max(xMax - xMin, yMax - yMin);
            AddVertex(cx - 20 * size, cy - 10 * size, none);
            AddVertex(cx + 20 * size, cy - 10 * size, none);
            AddVertex(cx, cy + 20 * size, none);
            triangles_.push_back(Triangle{{0, 1, 2}, {none, none, none}, true});
            vertexTriangle_ = {0, 0, 0};
            last_ = 0;
        }

        //! Inserts a point and returns its vertex; throws if it repeats a vertex or lies outside the starting triangle
        size_t Insert(double x, double y, size_t label)
        {
            const size_t vertex = AddVertex(x, y, label);
            std::array<double, 3> orientations;
            const size_t t = Locate(x, y, orientations);
            size_t zeros = 0, zeroEdge = 0;
            for(size_t i = 0; i < 3; ++i)
            {
                if(orientations[i] == 0)
                {
                    ++zeros;
                    zeroEdge = i;
                }
            }
            if(zeros >= 2)
            {
                const size_t other = triangles_[t].vertices[3 - zeroEdge - FirstZero(orientations)];
                MadVoro::Exception::MadVoroException eo("Delaunay2D::Insert: the point repeats a vertex");
                eo.addEntry("x", x);
                eo.addEntry("y", y);
                eo.addEntry("Label", label);
                eo.addEntry("Label of the vertex", labels_[other]);
                throw eo;
            }
            std::vector<size_t> &stack = flipStack_;
            stack.clear();
            if(zeros == 0)
            {
                SplitTriangle(t, vertex, stack);
            }
            else
            {
                SplitEdge(t, zeroEdge, vertex, stack);
            }
            while(not stack.empty())
            {
                const size_t tri = stack.back();
                stack.pop_back();
                FlipIfIllegal(tri, vertex, stack);
            }
            return vertex;
        }

        size_t VertexCount(void) const { return x_.size(); }

        bool IsSuper(size_t vertex) const { return vertex < 3; }

        double X(size_t vertex) const { return x_[vertex]; }

        double Y(size_t vertex) const { return y_[vertex]; }

        size_t Label(size_t vertex) const { return labels_[vertex]; }

        const std::vector<Triangle> &Triangles(void) const { return triangles_; }

        std::vector<Triangle> &Triangles(void) { return triangles_; }

        //! The triangles around a vertex, counterclockwise (a corner of the starting triangle has an open star)
        void Star(size_t vertex, std::vector<size_t> &res) const
        {
            res.clear();
            const size_t first = vertexTriangle_[vertex];
            size_t t = first;
            do
            {
                res.push_back(t);
                const Triangle &tri = triangles_[t];
                const size_t k = IndexOf(tri, vertex);
                t = tri.neighbors[(k + 1) % 3];
            } while(t != first and t != none and res.size() <= triangles_.size());
        }

        //! Circumcentre and circumradius of a triangle (relative to its first vertex, for accuracy)
        void Circumcircle(size_t triangle, double &cx, double &cy, double &radius) const
        {
            const std::array<size_t, 3> &v = triangles_[triangle].vertices;
            const double bx = x_[v[1]] - x_[v[0]], by = y_[v[1]] - y_[v[0]], ccx = x_[v[2]] - x_[v[0]], ccy = y_[v[2]] - y_[v[0]];
            const double b2 = bx * bx + by * by, c2 = ccx * ccx + ccy * ccy;
            const double d = 2 * (bx * ccy - by * ccx);
            const double ux = (ccy * b2 - by * c2) / d, uy = (bx * c2 - ccx * b2) / d;
            cx = x_[v[0]] + ux;
            cy = y_[v[0]] + uy;
            radius = std::sqrt(ux * ux + uy * uy);
        }

        //! Exact orientation of vertices a, b, c
        double Orient(size_t a, size_t b, size_t c) const
        {
            return Plane2D::Orient<PointT>(x_[a], y_[a], x_[b], y_[b], x_[c], y_[c]);
        }

        //! Exact in-circle test of vertex d against the circle of triangle t
        double InCircle(size_t t, size_t d) const
        {
            const std::array<size_t, 3> &v = triangles_[t].vertices;
            return Plane2D::InCircle<PointT>(x_[v[0]], y_[v[0]], x_[v[1]], y_[v[1]], x_[v[2]], y_[v[2]], x_[d], y_[d]);
        }

    private:
        std::vector<double> x_, y_;
        std::vector<size_t> labels_;
        std::vector<Triangle> triangles_;
        std::vector<size_t> vertexTriangle_; // a triangle of each vertex
        std::vector<size_t> flipStack_;
        size_t last_ = 0; // where the next point location starts

        size_t AddVertex(double x, double y, size_t label)
        {
            x_.push_back(x);
            y_.push_back(y);
            labels_.push_back(label);
            vertexTriangle_.push_back(none);
            return x_.size() - 1;
        }

        static size_t IndexOf(const Triangle &tri, size_t vertex)
        {
            return tri.vertices[0] == vertex ? 0 : (tri.vertices[1] == vertex ? 1 : 2);
        }

        static size_t FirstZero(const std::array<double, 3> &orientations)
        {
            return orientations[0] == 0 ? 0 : (orientations[1] == 0 ? 1 : 2);
        }

        /*! Visibility walk from the last triangle: crosses an edge the point is strictly behind. Returns the triangle that
          contains the point (closed) and the orientation of the point against each edge (the edge opposite vertex i).
         */
        size_t Locate(double x, double y, std::array<double, 3> &orientations) const
        {
            size_t t = last_;
            size_t start = 0;
            const size_t limit = 4 * triangles_.size() + 16;
            for(size_t steps = 0; steps < limit; ++steps)
            {
                const Triangle &tri = triangles_[t];
                size_t next = none;
                for(size_t s = 0; s < 3; ++s)
                {
                    const size_t i = (start + s) % 3;
                    const size_t a = tri.vertices[(i + 1) % 3], b = tri.vertices[(i + 2) % 3];
                    orientations[i] = Plane2D::Orient<PointT>(x_[a], y_[a], x_[b], y_[b], x, y);
                    if(orientations[i] < 0)
                    {
                        next = tri.neighbors[i];
                        if(next == none)
                        {
                            MadVoro::Exception::MadVoroException eo("Delaunay2D::Insert: the point is outside the starting triangle");
                            eo.addEntry("x", x);
                            eo.addEntry("y", y);
                            throw eo;
                        }
                        break;
                    }
                }
                if(next == none)
                {
                    return t;
                }
                t = next;
                start = (start + 1) % 3; // rotating the first edge tried keeps the walk from cycling
            }
            MadVoro::Exception::MadVoroException eo("Delaunay2D::Insert: point location did not converge");
            eo.addEntry("x", x);
            eo.addEntry("y", y);
            throw eo;
        }

        size_t NewTriangle(void)
        {
            triangles_.push_back(Triangle{{none, none, none}, {none, none, none}, true});
            return triangles_.size() - 1;
        }

        //! The neighbor pointer of `outside` that pointed to `from` now points to `to`
        void Relink(size_t outside, size_t from, size_t to)
        {
            if(outside == none)
            {
                return;
            }
            std::array<size_t, 3> &neighbors = triangles_[outside].neighbors;
            for(size_t i = 0; i < 3; ++i)
            {
                if(neighbors[i] == from)
                {
                    neighbors[i] = to;
                    return;
                }
            }
        }

        void Set(size_t t, size_t v0, size_t v1, size_t v2, size_t n0, size_t n1, size_t n2)
        {
            Triangle &tri = triangles_[t];
            tri.vertices = {v0, v1, v2};
            tri.neighbors = {n0, n1, n2};
            tri.checkBig = true;
            vertexTriangle_[v0] = t;
            vertexTriangle_[v1] = t;
            vertexTriangle_[v2] = t;
            last_ = t;
        }

        // The point strictly inside t = (v0, v1, v2): three triangles, each with p in place of one vertex
        void SplitTriangle(size_t t, size_t p, std::vector<size_t> &stack)
        {
            const Triangle old = triangles_[t];
            const size_t v0 = old.vertices[0], v1 = old.vertices[1], v2 = old.vertices[2];
            const size_t n0 = old.neighbors[0], n1 = old.neighbors[1], n2 = old.neighbors[2];
            const size_t t0 = t, t1 = NewTriangle(), t2 = NewTriangle();
            Set(t0, p, v1, v2, n0, t1, t2);
            Set(t1, v0, p, v2, t0, n1, t2);
            Set(t2, v0, v1, p, t0, t1, n2);
            Relink(n1, t, t1);
            Relink(n2, t, t2);
            stack.push_back(t0);
            stack.push_back(t1);
            stack.push_back(t2);
        }

        // The point on the edge of t opposite vertex i, shared with u: four triangles
        void SplitEdge(size_t t, size_t i, size_t p, std::vector<size_t> &stack)
        {
            const Triangle oldT = triangles_[t];
            const size_t a = oldT.vertices[i], b = oldT.vertices[(i + 1) % 3], c = oldT.vertices[(i + 2) % 3];
            const size_t nb = oldT.neighbors[(i + 1) % 3], nc = oldT.neighbors[(i + 2) % 3]; // across c-a and a-b
            const size_t u = oldT.neighbors[i];
            if(u == none)
            {
                MadVoro::Exception::MadVoroException eo("Delaunay2D::Insert: the point is on the edge of the starting triangle");
                eo.addEntry("x", x_[p]);
                eo.addEntry("y", y_[p]);
                throw eo;
            }
            const Triangle oldU = triangles_[u];
            size_t j = 0;
            while(oldU.neighbors[j] != t)
            {
                ++j;
            }
            // u = (w, c, b) counterclockwise, w opposite the shared edge
            const size_t w = oldU.vertices[j];
            const size_t mb = oldU.neighbors[IndexOf(oldU, b)], mc = oldU.neighbors[IndexOf(oldU, c)]; // across w-c and b-w
            const size_t ta = t, tb = NewTriangle(), ua = u, ub = NewTriangle();
            Set(ta, a, b, p, ub, tb, nc);
            Set(tb, a, p, c, ua, nb, ta);
            Set(ua, w, c, p, tb, ub, mb);
            Set(ub, w, p, b, ta, mc, ua);
            Relink(nb, t, tb);
            Relink(mc, u, ub);
            stack.push_back(ta);
            stack.push_back(tb);
            stack.push_back(ua);
            stack.push_back(ub);
        }

        // Lawson flip of the edge of t opposite p if the vertex across it is inside t's circle
        void FlipIfIllegal(size_t t, size_t p, std::vector<size_t> &stack)
        {
            const Triangle oldT = triangles_[t];
            const size_t k = IndexOf(oldT, p);
            const size_t u = oldT.neighbors[k];
            if(u == none)
            {
                return;
            }
            const Triangle oldU = triangles_[u];
            size_t j = 0;
            while(oldU.neighbors[j] != t)
            {
                ++j;
            }
            const size_t w = oldU.vertices[j];
            if(not (InCircle(t, w) > 0))
            {
                return;
            }
            // t = (p, a, b), u = (w, b, a); after the flip (p, a, w) and (p, w, b)
            const size_t a = oldT.vertices[(k + 1) % 3], b = oldT.vertices[(k + 2) % 3];
            const size_t ntb = oldT.neighbors[(k + 2) % 3], nta = oldT.neighbors[(k + 1) % 3]; // across p-a and b-p
            const size_t nub = oldU.neighbors[IndexOf(oldU, b)], nua = oldU.neighbors[IndexOf(oldU, a)]; // across a-w and w-b
            Set(t, p, a, w, nub, u, ntb);
            Set(u, p, w, b, nua, nta, t);
            Relink(nub, u, t);
            Relink(nta, t, u);
            stack.push_back(t);
            stack.push_back(u);
        }
    };
} // namespace MadVoro

#endif // MADVORO_DELAUNAY2D_HPP
