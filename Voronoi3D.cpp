#include <optional>
#include "Voronoi3D.hpp"
#include <spatial_ds/utils/BoundingBox.hpp>
#include "3D/elementary/Mat33.hpp"
#include "utils/Predicates3D.hpp"
#include "misc/universal_error.hpp"
#include "misc/utils.hpp"
#include "misc/io3D.hpp"
#include "geometry/Intersections.hpp"
#include "misc/int2str.hpp"

#ifdef RICH_MPI
#include "mpi/mpi_commands.hpp"

#include "3D/range/finders/BruteForce.hpp"
#include "3D/range/finders/RangeTree.hpp"
#include "3D/range/finders/OctTree.hpp"
#include "3D/range/finders/KDTree.hpp"
#include "3D/range/finders/GroupRangeTree.hpp"
#include <MeshDecomposer3D/environment/hilbert/HilbertTreeEnvAgent.hpp>

#include <MeshDecomposer3D/kernels/Rectangle.hpp>
#include <MeshDecomposer3D/kernels/SameRectangle.hpp>

#endif // RICH_MPI

// #define VORONOI_DEBUG

bool PointInPoly(Tessellation3D const &tess, Vector3D const &point, std::size_t index)
{
    face_vec const &faces = tess.GetCellFaces(index);
    vector<Vector3D> const &points = tess.GetFacePoints();
    std::size_t N = faces.size();
    std::array<Vector3D, 4> vec;
    for (std::size_t i = 0; i < N; ++i)
    {
        double R = fastsqrt(tess.GetArea(faces[i]));
        size_t N1 = 0;
        size_t N2 = 0;
        Vector3D V1, V2;
        size_t counter = 0;
        point_vec const &InFace = tess.GetPointsInFace(faces[i]);
        size_t NinFace = InFace.size();
        N1 = 1;
        V1 = points[InFace[(counter + 1) % NinFace]] - points[InFace[0]];
        while (fastabs(V1) < 0.01 * R)
        {
            ++counter;
            assert(counter < NinFace);
            V1 = points[InFace[(counter + 1) % NinFace]] - points[InFace[0]];
            ++N1;
        }
        V2 = points[InFace[(counter + 2) % NinFace]] - points[InFace[N1]];
        N2 = (counter + 2) % NinFace;
        while (fastabs(V2) < 0.01 * R || fastabs(CrossProduct(V1, V2)) < 0.0001 * tess.GetArea(faces[i]))
        {
            ++counter;
            if (counter > 2 * NinFace)
                break;
            V2 = points[InFace[(counter + 2) % NinFace]] - points[InFace[N1]];
            N2 = (counter + 2) % NinFace;
        }
        if (counter > 2 * NinFace)
        {
            std::cout << "Weird face in PointInPoly, cell " << index << " face " << faces[i] << " i " << i << " face area " << tess.GetArea(faces[i]) << std::endl;
            for (size_t j = 0; j < NinFace; ++j)
                std::cout << "Point j " << points[InFace[j]].x << "," << points[InFace[j]].y << "," << points[InFace[j]].z << std::endl;
            Vector3D normal = tess.GetFaceNeighbors(faces[i]).second == index ? tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).second) - tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).first) : tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).first) - tess.GetMeshPoint(tess.GetFaceNeighbors(faces[i]).second);
            if (ScalarProd(normal, point - points[InFace[0]]) < 0)
                return false;
        }
        else
        {
            vec[0] = points[InFace[0]];
            vec[1] = points[InFace.at(N1)];
            vec[2] = points[InFace.at(N2)];
            vec[3] = tess.GetMeshPoint(index);
            double s1 = orient3d(vec);
            vec[3] = point;
            double s2 = orient3d(vec);
            if (s1 * s2 < -0)
                return false;
        }
    }
    return true;
}

bool PointInPoly(std::vector<Face> const& faces, Vector3D const &point)
{
    std::size_t const N = faces.size();
    std::array<Vector3D, 4> vec;
    vec[3] = point;
    for (std::size_t i = 0; i < N; ++i)
    {
        vec[0] = faces[i].vertices[0];
        vec[1] = faces[i].vertices[1];
        vec[2] = faces[i].vertices[2];
        double const s = orient3d(vec);
        if (s > 0)
            return false;
    }
    return true;
}

namespace
{
#ifdef RICH_MPI
    void GetPastDuplicate(size_t point, vector<size_t> &res, vector<vector<size_t>> const &sorted_to_duplicate,
                                                vector<size_t> const &procs)
    {
        res.clear();
        for (size_t i = 0; i < procs.size(); ++i)
        {
            if (std::binary_search(sorted_to_duplicate[i].begin(), sorted_to_duplicate[i].end(), point))
                res.push_back(procs[i]);
        }
    }
#endif
    boost::multiprecision::cpp_dec_float_50 Calc33Det(std::array<boost::multiprecision::cpp_dec_float_50, 9> const &points)
    {
        return points[0] * (points[4] * points[8] - points[5] * points[7]) + points[1] * (points[5] * points[6] - points[3] * points[8]) + points[2] * (points[3] * points[7] - points[4] * points[6]);
    }
}

namespace
{
    bool ShouldCalcTetraRadius(Tetrahedron const &T, size_t Norg)
    {
        for (size_t i = 0; i < 4; ++i)
            if (T.points[i] < Norg)
                return true;
        return false;
    }

    void FirstCheckList(std::stack<std::size_t> &check_stack, vector<unsigned char> &future_check, size_t Norg,
                                            Delaunay3D const &del, vector<tetra_vec> const &PointsInTetra)
    {
        //        check_stack.empty();
        future_check.resize(Norg, 0);
        size_t Ntetra = del.tetras_.size();
        vector<unsigned char> tetra_check(Ntetra, 0);

        for (size_t i = 0; i < Ntetra; ++i)
        {
            Tetrahedron const &tetra = del.tetras_[i];
            for (size_t j = 0; j < 4; ++j)
            {
                if (tetra.points[j] >= Norg)
                {
                    for (size_t k = 0; k < 4; ++k)
                    {
                        size_t tetcheck = tetra.points[k];
                        if (tetra.points[k] < Norg)
                        {
                            size_t ntet = PointsInTetra[tetcheck].size();
                            for (size_t z = 0; z < ntet; ++z)
                                tetra_check[PointsInTetra[tetcheck][z]] = 1;
                        }
                    }
                    break;
                }
            }
        }
        for (size_t i = 0; i < Ntetra; ++i)
        {
            if (tetra_check[i] == 1)
            {
                Tetrahedron const &tetra = del.tetras_[i];
                for (size_t j = 0; j < 4; ++j)
                {
                    if (tetra.points[j] < Norg)
                        future_check[tetra.points[j]] = 1;
                }
            }
        }
        for (size_t i = 0; i < Norg; ++i)
            if (future_check[i] == 1)
                check_stack.push(i);
    }

    vector<Face> BuildBox(Vector3D const &ll, Vector3D const &ur)
    {
        double dx = ur.x - ll.x;
        double dy = ur.y - ll.y;
        double dz = ur.z - ll.z;
        vector<Face> res(6);
        vector<Vector3D> points;
        points.push_back(ll);
        points.push_back(ll + Vector3D(dx, 0, 0));
        points.push_back(ll + Vector3D(dx, dy, 0));
        points.push_back(ll + Vector3D(0, dy, 0));
        points.push_back(ll + Vector3D(0, 0, dz));
        points.push_back(ll + Vector3D(dx, 0, dz));
        points.push_back(ll + Vector3D(dx, dy, dz));
        points.push_back(ll + Vector3D(0, dy, dz));
        points.push_back(ur);
        res[0].vertices.push_back(points[0]);
        res[0].vertices.push_back(points[1]);
        res[0].vertices.push_back(points[2]);
        res[0].vertices.push_back(points[3]);
        res[1].vertices.push_back(points[0]);
        res[1].vertices.push_back(points[4]);
        res[1].vertices.push_back(points[5]);
        res[1].vertices.push_back(points[1]);
        res[2].vertices.push_back(points[3]);
        res[2].vertices.push_back(points[7]);
        res[2].vertices.push_back(points[4]);
        res[2].vertices.push_back(points[0]);
        res[3].vertices.push_back(points[2]);
        res[3].vertices.push_back(points[6]);
        res[3].vertices.push_back(points[7]);
        res[3].vertices.push_back(points[3]);
        res[4].vertices.push_back(points[1]);
        res[4].vertices.push_back(points[5]);
        res[4].vertices.push_back(points[6]);
        res[4].vertices.push_back(points[2]);
        res[5].vertices.push_back(points[5]);
        res[5].vertices.push_back(points[4]);
        res[5].vertices.push_back(points[7]);
        res[5].vertices.push_back(points[6]);
        return res;
    }

#ifdef RICH_MPI
    vector<Vector3D> GetBoxNormals(Vector3D const &ll, Vector3D const &ur, vector<Face> const& box_faces_)
    {
        const vector<Face> faces = box_faces_.empty() ? BuildBox(ll, ur) : box_faces_;
        vector<Vector3D> res(faces.size());
        size_t N = res.size();
        for (size_t i = 0; i < N; ++i)
        {
            CrossProduct(faces[i].vertices[2] - faces[i].vertices[0], faces[i].vertices[1] - faces[i].vertices[0], res[i]);
            res[i] *= 1.0 /abs(res[i]);
        }
        return res;
    }

    size_t BoxIndex(vector<Vector3D> const &fnormals, const Vector3D &normal)
    {
        double max_angle = ScalarProd(fnormals[0], normal);
        size_t loc = 0;
        size_t N = fnormals.size();
        for (size_t i = 1; i < N; i++)
        {
            double temp = ScalarProd(fnormals[i], normal);
            if (temp > max_angle)
            {
                max_angle = temp;
                loc = i;
            }
        }
        return loc;
    }
#endif

    double CleanDuplicates(std::array<size_t, 128> const &indeces, const vector<Vector3D> &points,
                                                 boost::container::small_vector<size_t, 8> &res, double R,
                                                 std::array<double, 128> &diffs,
                                                 std::array<Vector3D, 128> &vtemp, const size_t N)
    {
        res.clear();
        for (size_t i = 0; i < N; ++i)
            vtemp[i] = points[indeces[i]];
        for (size_t i = N - 1; i > 0; --i)
        {
            vtemp[i].x -= vtemp[i - 1].x;
            vtemp[i].y -= vtemp[i - 1].y;
            vtemp[i].z -= vtemp[i - 1].z;
        }
        vtemp[0] -= points[indeces[N - 1]];
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(max \
                                                     : R)
#endif
        for (size_t i = 0; i < N; ++i)
        {
            diffs[i] = ScalarProd(vtemp[i], vtemp[i]);
            R = std::max(R, diffs[i]);
        }
        for (size_t i = 0; i < N; ++i)
            if (diffs[i] > R * 1e-16)
                res.push_back(indeces[i]);
        return R;
    }

    bool CleanSameLine(boost::container::small_vector<size_t, 8> &indeces, vector<Vector3D> const& face_points, std::array<double, 128> &area_vec_temp)
    {
        point_vec old;
        size_t const N = indeces.size();
        double const small_fraction = 1e-14;
        // double const medium_fraction = 3e-1;
        // Find correct normal
        Vector3D good_normal;
        for(size_t i = 0; i < N; ++i)
        {
            area_vec_temp[i] = fastabs(CrossProduct(face_points[indeces[i]] - face_points[indeces[(N + i - 1) % N]], face_points[indeces[(i + 1) % N]] 
            - face_points[indeces[(N + i - 1) % N]]));
            old.push_back(indeces[i]);
        }

        double max_value = area_vec_temp[0];
        double second_max_value = max_value;
        size_t max_index = 0, second_max_index = 0;
        for(size_t i = 1; i < N; ++i)
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

        double const area_scale = area_vec_temp[max_index];
        good_normal = CrossProduct(face_points[indeces[max_index]] - face_points[indeces[(N + max_index - 1) % N]], face_points[indeces[(max_index + 1) % N]] - face_points[indeces[(N + max_index - 1) % N]]);
        good_normal *= 1.0 / fastabs(good_normal);

        size_t Nindeces = indeces.size();
        for(size_t i = 0; i < Nindeces; ++i)
        {
            Vector3D normal_temp = CrossProduct(face_points[indeces[i]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]], face_points[indeces[(i + 1) % Nindeces]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]]);
            double const area = fastabs(normal_temp);
            normal_temp *= 1.0 / (100 * std::numeric_limits<double>::min() + area);
            if((area < area_scale * small_fraction) || (ScalarProd(normal_temp, good_normal) < 0.9999))
            {
                indeces.erase(indeces.begin() + i);
                if(i == indeces.size() - 1)
                    break;
                --i;
                Nindeces = indeces.size();
            }
        }

        if(Nindeces < 3)
        {
            indeces = old;
            max_index = second_max_index;
            good_normal = CrossProduct(face_points[indeces[max_index]] - face_points[indeces[(N + max_index - 1) % N]], face_points[indeces[(max_index + 1) % N]] - face_points[indeces[(N + max_index - 1) % N]]);
            good_normal *= 1.0 / fastabs(good_normal);

            for(size_t i = 0; i < Nindeces; i++)
            {
                Vector3D normal_temp = CrossProduct(face_points[indeces[i]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]], face_points[indeces[(i + 1) % Nindeces]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]]);
                double const area = fastabs(normal_temp);
                normal_temp *= 1.0 / area;
                if((area < area_scale * small_fraction) || ScalarProd(normal_temp, good_normal) < 0.9999)
                {
                    indeces.erase(indeces.begin() + i);
                    if(i == indeces.size() - 1)
                        break;
                    --i;
                    Nindeces = indeces.size();
                }
            }
        }

        if(Nindeces < 3)
        {
            Nindeces = N;
            indeces = old;
            UniversalError eo("Bad CleanSameLine");
            eo.addEntry("N", N);
            eo.addEntry("good normal x", good_normal.x);
            eo.addEntry("good normal y", good_normal.y);
            eo.addEntry("good normal z", good_normal.z);
            for(size_t i = 0; i < N; ++i)
            {
                eo.addEntry("index", old[i]);
                eo.addEntry("area_vec_temp", area_vec_temp[i]);
                eo.addEntry("point " + std::to_string(i) + " x", face_points[indeces[i]].x);
                eo.addEntry("point " + std::to_string(i) + " y", face_points[indeces[i]].y);
                eo.addEntry("point " + std::to_string(i) + " z", face_points[indeces[i]].z);
                Vector3D normal_temp = CrossProduct(face_points[indeces[i]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]], face_points[indeces[(i + 1) % Nindeces]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]]);
                normal_temp *= (1.0 / abs(normal_temp));
                eo.addEntry("normal " + std::to_string(i) + " x", normal_temp.x);
                eo.addEntry("normal " + std::to_string(i) + " y", normal_temp.y);
                eo.addEntry("normal " + std::to_string(i) + " z", normal_temp.z);
                eo.addEntry("dot", ScalarProd(good_normal, normal_temp));
            }
            throw eo;
        }

        return true;
    }


    void MakeRightHandFace(boost::container::small_vector<size_t, 8> &indeces, Vector3D const &point, vector<Vector3D> const &face_points,
                                                 std::array<size_t, 128> &temp, double areascale)
    {
        Vector3D V1, V2;
        size_t counter = 0;
        const size_t N = indeces.size();
        V1 = face_points[indeces[counter + 1]];
        V1 -= face_points[indeces[counter]];
        double AScale = 1e-14 * areascale;
        while (ScalarProd(V1, V1) < AScale)
        {
            ++counter;
            assert(counter < N);
            V1 = face_points[indeces[(counter + 1) % N]];
            V1 -= face_points[indeces[counter]];
        }
        V2 = face_points[indeces[(counter + 2) % N]];
        V2 -= face_points[indeces[(counter + 1) % N]];
        while (ScalarProd(V2, V2) < AScale)
        {
            ++counter;
            assert(counter < 2 * N);
            V2 = face_points[indeces[(counter + 2) % N]];
            V2 -= face_points[indeces[(counter + 1) % N]];
        }
        // Do we need to flip handness?
        if (ScalarProd(CrossProduct(V1, V2), point - face_points[indeces[0]]) > 0)
        {
            const size_t Ninner = indeces.size();
#ifdef __INTEL_COMPILER
#pragma omp simd early_exit
#endif
            for (size_t j = 0; j < Ninner; ++j)
                temp[j] = indeces[j];
#ifdef __INTEL_COMPILER
#pragma omp simd early_exit
#endif
            for (size_t i = 0; i < N; ++i)
                indeces[i] = temp[(N - i - 1)];
        }
    }

    size_t NextLoopTetra(Tetrahedron const &cur_tetra, size_t last_tetra, size_t N0, size_t N1)
    {
        size_t i = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (; i < 4; i++)
        {
            size_t point = cur_tetra.points[i];
            if (point != N0 && point != N1 && cur_tetra.neighbors[i] != last_tetra)
                break;
        }
        if(i >= 4)
            throw UniversalError("Bad NextLoopTetra");
        return cur_tetra.neighbors[i];
    }

    void CalcFaceAreaCM(boost::container::small_vector<size_t, 8> const &indeces, std::vector<Vector3D> const &allpoints,
                                            std::array<Vector3D, 128> &points, double &Area, Vector3D &CM,
                                            std::array<double, 128> &Atemp)
    {
        //CM.Set(0.0, 0.0, 0.0);
        size_t Nloop = indeces.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (size_t i = 0; i < Nloop; i++)
            points[i] = allpoints[indeces[i]];
        Nloop -= 2;
        Area = 0;
        //Vector3D temp3, temp4, temp5;
        for (size_t i = 0; i < Nloop; i++)
        {
            //temp4.Set(points[i + 1].x - points[0].x, points[i + 1].y - points[0].y, points[i + 1].z - points[0].z);
            Vector3D temp4(points[i + 1].x - points[0].x, points[i + 1].y - points[0].y, points[i + 1].z - points[0].z);
            //temp5.Set(points[i + 2].x - points[0].x, points[i + 2].y - points[0].y, points[i + 2].z - points[0].z);
            Vector3D temp5(points[i + 2].x - points[0].x, points[i + 2].y - points[0].y, points[i + 2].z - points[0].z);
            Vector3D temp3;
            CrossProduct(temp4, temp5, temp3);
            Atemp[i] = 0.3333333333333333 * 0.5 * fastsqrt(ScalarProd(temp3, temp3));
        }
        double x = 0, y = 0, z = 0;
#ifdef __INTEL_COMPILER
#pragma vector aligned
        //#pragma omp simd reduction(+:x, y, z, Area)
#endif
        for (size_t i = 0; i < Nloop; i++)
        {
            double A = Atemp[i];
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
        CM.Set(x, y, z);
        CM *= (1.0 / (Area + std::numeric_limits<double>::min() * 100)); //prevent overflow
    }

    bool PointInDomain(Vector3D const &ll, Vector3D const &ur, Vector3D const &point)
    {
        if (point.x > ll.x && point.x < ur.x && point.y > ll.y && point.y < ur.y && point.z > ll.z && point.z < ur.z)
            return true;
        else
            return false;
    }

    Vector3D MirrorPoint(Face const &face, Vector3D const &point)
    {
        Vector3D normal = CrossProduct(face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]);
        normal = normal / abs(normal);
        return point - (2 * ScalarProd(point - face.vertices[0], normal)) * normal;
    }
}

size_t Voronoi3D::SetPointTetras(void)
{
    std::vector<Tetrahedron> &tetras = this->del_.tetras_;
    std::vector<std::pair<size_t, Tetrahedron>> &changed_tetras = this->del_.changed_tetras_;
    boost::container::flat_set<size_t> const &empty_tetras = this->del_.empty_tetras_;
    std::vector<size_t> &newTetras = this->del_.newTetras_;
    PointTetras_.resize(this->Norg_);
    PointTetras_.shrink_to_fit();
    // static vector<tetra_vec> tmpPointTetras;

    #ifdef USE_VCL_VECTORIZATION
        Vec4uq _Norg(this->Norg_);
    #endif // USE_VCL_VECTORIZATION

    size_t Ntetra = tetras.size();
    size_t bigtet(0);
    bool has_good, has_big;
    // change empty tetras to be not relevant
    for (boost::container::flat_set<size_t>::const_iterator it = empty_tetras.begin(); it != empty_tetras.end(); ++it)
    {
        Tetrahedron &tetra = tetras[*it];
        for (size_t i = 0; i < 4; ++i)
        {
            if(tetra.points[i] < this->Norg_)
            {
                tetra_vec &points_tetras = PointTetras_[tetra.points[i]];
                auto it2 = std::find(points_tetras.begin(), points_tetras.end(), *it);
                if(it2 != points_tetras.end()) {
                    *it2 = points_tetras.back();
                    points_tetras.pop_back();
                }
            }
            tetras[*it].points[i] = std::numeric_limits<std::size_t>::max();
            tetras[*it].neighbors[i] = std::numeric_limits<std::size_t>::max();
        }
    }

    for (const std::pair<size_t, Tetrahedron> &tetraInfo : changed_tetras)
    {
        const size_t &tetraIndex = tetraInfo.first;
        const Tetrahedron &tetra = tetraInfo.second;

        for (size_t i = 0; i < 4; ++i)
        {
            if(tetra.points[i] < this->Norg_)
            {
                tetra_vec &points_tetras = PointTetras_[tetra.points[i]];
                auto it2 = std::find(points_tetras.begin(), points_tetras.end(), tetraIndex);
                if(it2 != points_tetras.end())
                {
                    *it2 = points_tetras.back();
                    points_tetras.pop_back();
                }
            }
        }
    }

    size_t counter = 0;

    // std::sort(newTetras.begin(), newTetras.end());
    // auto it = std::unique(newTetras.begin(), newTetras.end());
    // newTetras.resize(std::distance(newTetras.begin(), it));
    // size_t redundentNum = Ntetra - newTetras.size();
    // #ifdef RICH_MPI
    //     int rank;
    //     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //     MPI_Allreduce(MPI_IN_PLACE, &redundentNum, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    //     if(rank == 0)
    //     {
    //         std::cout << "Redundent number max: " << redundentNum << std::endl;
    //     }
    // #endif // RICH_MPI
    
    // for (size_t i = 0; i < Ntetra; ++i)
    for(size_t i : newTetras)
    {
        Tetrahedron &tet = tetras[i];
        this->R_[i] = RADIUS_UNINITIALIZED;

        if(not tet.newTetra)
        {
            continue;
        }
        assert(tet.newTetra);
        tet.newTetra = false;

        has_good = false;
        has_big = false;
        counter++;

        #ifdef USE_VCL_VECTORIZATION
            Vec4uq _points(tet.points[0], tet.points[1], tet.points[2], tet.points[3]);
            Vec4qb cmp = (_points < _Norg);
        #endif // USE_VCL_VECTORIZATION
        
        for(int j = 0; j < 4; ++j)
        {
        #ifdef USE_VCL_VECTORIZATION
            if(cmp[j])
        #else // USE_VCL_VECTORIZATION
            if(tet.points[j] < Norg)
        #endif // USE_VCL_VECTORIZATION
            {
                has_good = true;
                PointTetras_[tet.points[j]].push_back(i);
            }
            else
            {
                has_big = true;
            }
        }
        if(has_big and has_good)
        {
            bigtet = i;
        }
    }

    changed_tetras.clear();
    newTetras.clear();

    // std::cout << "Tetra counter: " << counter << " / " << Ntetra << std::endl;
    return bigtet;
}

Voronoi3D::Voronoi3D(std::vector<Face> const& box_faces) : Voronoi3D()
{
    this->box_faces_ = box_faces;
    size_t const Nfaces = box_faces.size();
    if(Nfaces < 4)
        throw UniversalError("Zero face vector in Voronoi3D constructor");
    ll_ = box_faces[0].vertices[0];
    ur_ = ll_;
    for(size_t i = 0; i < Nfaces; ++i)
    {
        size_t const Nvertices = box_faces[i].vertices.size();
        for(size_t j = 0; j < Nvertices; ++j)
        {
            ll_.x = std::min(ll_.x, box_faces[i].vertices[j].x);
            ll_.y = std::min(ll_.y, box_faces[i].vertices[j].y);
            ll_.z = std::min(ll_.z, box_faces[i].vertices[j].z);
            ur_.x = std::max(ur_.x, box_faces[i].vertices[j].x);
            ur_.y = std::max(ur_.y, box_faces[i].vertices[j].y);
            ur_.z = std::max(ur_.z, box_faces[i].vertices[j].z);
        }
    }
}

Voronoi3D::Voronoi3D(Vector3D const &ll, Vector3D const &ur) : ll_(ll), ur_(ur), Norg_(0), bigtet_(0), set_temp_(std::set<int>()), stack_temp_(std::stack<int>()),
                                                              del_(Delaunay3D()), PointTetras_(vector<tetra_vec>()), R_(vector<double>()), tetra_centers_(vector<Vector3D>()),
                                                              FacesInCell_(vector<face_vec>()),
                                                              PointsInFace_(vector<point_vec>()),
                                                              FaceNeighbors_(vector<std::pair<std::size_t, std::size_t>>()),
                                                              CM_(vector<Vector3D>()), Face_CM_(vector<Vector3D>()),
                                                              volume_(vector<double>()), area_(vector<double>()), 
                                                              #ifdef RICH_MPI
                                                                sentprocs_(vector<int>()), sentpoints_(vector<vector<std::size_t>>()),  duplicatedprocs_(vector<int>()), 
                                                                duplicated_points_(vector<vector<std::size_t>>()), Nghost_(vector<vector<std::size_t>>()), self_index_(vector<std::size_t>()), 
                                                              #endif // RICH_MPI
                                                              temp_points_(std::array<Vector3D, 4>()), temp_points2_(std::array<Vector3D, 5>()), box_faces_(std::vector<Face>()),
                                                              #ifdef RICH_MPI
                                                                pointsManager(nullptr),
                                                                allMyPoints(), 
                                                              #endif // RICH_MPI
                                                              indicesInAllMyPoints()
{
    this->box_faces_ = BuildBox(this->ll_, this->ur_);
    #ifdef RICH_MPI
        // initialize points manager
        this->pointsManager = std::shared_ptr<HilbertPointsManager<Vector3D, VoronoiPayload>>(new HilbertPointsManager<Vector3D, VoronoiPayload>(this->ll_, this->ur_));
    #endif // RICH_MPI
}

Voronoi3D::Voronoi3D() : Voronoi3D(Vector3D(), Vector3D())
{}

void Voronoi3D::CalcRigidCM(std::size_t face_index)
{
    Vector3D normal = normalize(del_.points_[FaceNeighbors_[face_index].first] - del_.points_[FaceNeighbors_[face_index].second]);
    std::size_t real, other;
    if (FaceNeighbors_[face_index].first >= Norg_)
    {
        real = FaceNeighbors_[face_index].second;
        other = FaceNeighbors_[face_index].first;
    }
    else
    {
        real = FaceNeighbors_[face_index].first;
        other = FaceNeighbors_[face_index].second;
    }
    CM_[other] = CM_[real] - 2 * normal * ScalarProd(normal, CM_[real] - tetra_centers_[PointsInFace_[face_index][0]]);
}

vector<Vector3D> Voronoi3D::CreateBoundaryPoints(vector<std::pair<std::size_t, std::size_t>> const &to_duplicate,
                                                 vector<vector<size_t>> &past_duplicate)
{
    size_t Ncheck = to_duplicate.size();
    vector<std::pair<std::size_t, std::size_t>> to_add;
    to_add.reserve(Ncheck);
    vector<Face> faces = box_faces_.empty() ? BuildBox(ll_, ur_) : box_faces_;
    vector<Vector3D> res;
    bool first_time = past_duplicate.empty();
    if (first_time)
        past_duplicate.resize(faces.size());
    for (std::size_t i = 0; i < Ncheck; ++i)
    {
        if (first_time || !std::binary_search(past_duplicate[to_duplicate[i].first].begin(),
                                                                                    past_duplicate[to_duplicate[i].first].end(), to_duplicate[i].second))
        {
            res.push_back(MirrorPoint(faces[to_duplicate[i].first], del_.points_[to_duplicate[i].second]));
            to_add.push_back(to_duplicate[i]);
        }
    }
    for (size_t i = 0; i < to_add.size(); ++i)
        past_duplicate[to_add[i].first].push_back(to_add[i].second);
    for (size_t i = 0; i < past_duplicate.size(); ++i)
        std::sort(past_duplicate[i].begin(), past_duplicate[i].end());
    return res;
}

#ifdef RICH_MPI
    vector<vector<std::size_t>> const &Voronoi3D::GetGhostIndeces(void) const
    {
        return Nghost_;
    }
#endif // RICH_MPI

/**
 * gets a point index, and returns the maximal radius of the tetrahedra containing that point.
 * @param index the index of the point (within the points list)
*/
double Voronoi3D::GetMaxRadius(const size_t &index) const
{
    std::size_t N = PointTetras_[index].size();
    double res = 0;
    #ifdef __INTEL_COMPILER
    #pragma ivdep
    #endif
    for(std::size_t i = 0; i < N; ++i)
    {
        res = std::max(res, GetRadius(PointTetras_[index][i]));
    }
    return res;
}

/**
 * gets a point index, and returns the minimal radius of the tetrahedra containing that point.
 * @param index the index of the point (within the points list)
*/
double Voronoi3D::GetMinRadius(const size_t &index) const
{
    std::size_t N = PointTetras_[index].size();
    double res = std::numeric_limits<double>::max();
    #ifdef __INTEL_COMPILER
    #pragma ivdep
    #endif
    for(std::size_t i = 0; i < N; ++i)
    {
        res = std::min(res, GetRadius(PointTetras_[index][i]));
    }
    return res;
}

/**
 * if the initial box does not exist, builds its faces according to the leftmost and rightmost points.
 * If it does, does not build the faces again.
 * @return the normals to the faces
*/
void Voronoi3D::InitialBoxBuild(std::vector<Face> &box, std::vector<Vector3D> &normals)
{
    box = box_faces_.empty() ? BuildBox(this->ll_, this->ur_) : this->box_faces_;
    size_t Nfaces = box.size();
    normals.resize(Nfaces);

    // calculates the normals for each one of the box's faces
    for (size_t i = 0; i < Nfaces; ++i)
    {
        normals[i] = CrossProduct(box[i].vertices[1] - box[i].vertices[0], box[i].vertices[2] - box[i].vertices[0]);
        normals[i] *= (1.0 / fastsqrt(ScalarProd(normals[i], normals[i])));
    }
}

/**
 * \author Maor Mizrachi
 * \brief Initializes internal data structures, for the voronoi build
*/
void Voronoi3D::BuildInitialize(size_t num_points)
{
    ++buildGeneration_;
    // assert(num_points > 0);
    // Clear data
    PointTetras_.clear();
    PointTetras_.shrink_to_fit();
    R_.clear();
    if(num_points > 0) R_.reserve(num_points * 11);
    tetra_centers_.clear();
    if(num_points > 0) tetra_centers_.reserve(num_points * 11);
    // Voronoi Data
    del_.Clean();
    FacesInCell_.clear();
    FacesInCell_.shrink_to_fit();
    PointsInFace_.clear();
    PointsInFace_.shrink_to_fit();
    FaceNeighbors_.clear();
    FaceNeighbors_.shrink_to_fit();
    CM_.clear();
    CM_.shrink_to_fit();
    Face_CM_.clear();
    Face_CM_.shrink_to_fit();
    volume_.clear();
    volume_.shrink_to_fit();
    area_.clear();
    area_.shrink_to_fit();
    Norg_ = num_points;
    #ifdef RICH_MPI
        duplicatedprocs_.clear();
        duplicatedprocs_.shrink_to_fit();
        duplicated_points_.clear();
        duplicated_points_.shrink_to_fit();
        Nghost_.clear();
        Nghost_.shrink_to_fit();
    #endif // RICH_MPI
}

#ifdef RICH_MPI
    #ifdef VORONOI_DEBUG
    namespace
    {
        template<typename T>
        void reportDuplications(const std::vector<T> &vector)
        {
            for(size_t i = 0; i < vector.size(); i++)
            {
                for(size_t j = 0; j < vector.size(); j++)
                {
                    if(i == j) continue;
                    if(vector[i] == vector[j])
                    {
                        std::cout << "duplication found in indices " << i << " and " << j << ": " << vector[i] << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, 2050);
                    }
                }
            }
        }
    }
    #endif // VORONOI_DEBUG

/**
 * \author Maor Mizrachi
 * \brief Checks if a certain point is under my responsibility
*/
bool Voronoi3D::PointInMyDomain(const Vector3D &point) const
{
    size_t containingCell = this->GetContainingCell(point);
    return this->IsPointInCell(point, containingCell);
}

inline int Voronoi3D::GetOwner(const Vector3D &point) const
{
    return this->pointsManager->getEnvironmentAgent()->getOwner(point);
}

inline void Voronoi3D::SetImbalanceTolerance(double tolerance)
{
    this->pointsManager->setImbalanceTolerance(tolerance);
}

std::tuple<std::vector<Vector3D>, std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>>
    Voronoi3D::InitialGhostPointsExchange(const MPI_Comm &comm) const
{
    int size;
    MPI_Comm_size(comm, &size);

    std::vector<std::vector<Vector3D>> sentPoints(size);
    std::vector<std::vector<size_t>> sentPointsIndices(size);

    
    for(int _rank = 0; _rank < size; _rank++)
    {
        size_t rankIndex = std::distance(this->real_duplicated_proc.begin(), std::find(this->real_duplicated_proc.begin(), this->real_duplicated_proc.end(), _rank));
        if(rankIndex != this->real_duplicated_proc.size())
        {
            // found
            std::vector<size_t> &sentIndices = sentPointsIndices[_rank];
            std::vector<Vector3D> &pointsToSend = sentPoints[_rank];

            // rank _rank is duplicated
            for(size_t pointIdx : this->real_duplicated_points[rankIndex])
            {
                if(pointIdx < this->Norg_)
                {
                    sentIndices.push_back(pointIdx);
                    pointsToSend.push_back(this->allMyPoints[pointIdx]);
                }
            }
        }
    }

    bool hasAnythingToSend = false;
    for(int _rank = 0; _rank < size && !hasAnythingToSend; _rank++)
        hasAnythingToSend = !sentPoints[_rank].empty();

    int globalHasData = 0;
    int localHasData = hasAnythingToSend ? 1 : 0;
    MPI_Allreduce(&localHasData, &globalHasData, 1, MPI_INT, MPI_MAX, comm);

    if(globalHasData == 0)
    {
        return std::tuple(std::vector<Vector3D>{},
                          sentPointsIndices,
                          std::vector<std::vector<size_t>>(size));
    }

    std::vector<std::vector<Vector3D>> recvPoints = MPI_Exchange_all_to_all_sparse(sentPoints, comm);

    std::vector<Vector3D> ghostPoints;
    std::vector<std::vector<size_t>> recvPointsIndices(size);
    size_t totalSize = 0;

    for(int _rank = 0; _rank < size; _rank++)
    {
        std::vector<Vector3D> &receivedFromRank = recvPoints[_rank];
        std::vector<size_t> &receivedIndices = recvPointsIndices[_rank];
        size_t receiving = receivedFromRank.size();
        for(size_t i = 0; i < receiving; i++)
        {
            receivedIndices.push_back(totalSize);
            totalSize++;
        }
        ghostPoints.insert(ghostPoints.end(), receivedFromRank.begin(), receivedFromRank.end());
    }

    return std::tuple(ghostPoints, sentPointsIndices, recvPointsIndices);
}

void Voronoi3D::FilterRealGhostPoints()
{
    this->real_duplicated_proc.clear();
    this->real_duplicated_points.clear();

    // std::vector<bool> isNecessaryRecvPoint(this->del_.points_.size(), false);
    // for(size_t pointIdx = 0; pointIdx < this->Norg_; pointIdx++)
    // {
    //     for(const size_t &neighborIdx : this->GetNeighbors(pointIdx))
    //     {
    //         isNecessaryRecvPoint[neighborIdx] = true;
    //     }
    // }

    // auto ifRecvCopyLambda = [&isNecessaryRecvPoint](const size_t &ghostPointIdx){return isNecessaryRecvPoint[ghostPointIdx];};

    for(size_t i = 0; i < this->duplicatedprocs_.size(); i++)
    {
        int _rank = this->duplicatedprocs_[i];
        this->real_duplicated_proc.push_back(_rank);
        std::vector<size_t> newSend;
        // check for any original sent point, if it has neighbors that belong to rank `_rank`. If yes, the point is necessary to be sent
        for(const size_t &pointIdxInBuild : this->duplicated_points_[i])
        {
            if(pointIdxInBuild >= this->Norg_)
            {
                // point was not participating in the last built
                continue;
            }
            bool foundNeighbor = false;
            for(const size_t &neighborIdx : this->GetNeighbors(pointIdxInBuild))
            {
                if(std::find(this->Nghost_[i].cbegin(), this->Nghost_[i].cend(), neighborIdx) != this->Nghost_[i].cend())
                {
                    // found a neighbor of `pointIdx` which is a ghost point of mine
                    foundNeighbor = true;
                    break;
                }
            }
            if(foundNeighbor)
            {
                size_t pointIdxInAll = this->indicesInAllMyPoints[pointIdxInBuild];
                newSend.push_back(pointIdxInAll);
            }
        }
        this->real_duplicated_points.emplace_back(newSend);
    }
}

/**
 * \author Maor Mizrachi
 * \brief Updates the duplicated points array
*/
void Voronoi3D::UpdateDuplicatedPoints(const std::vector<int> &sentProc, const std::vector<std::vector<size_t>> &sentPoints)
{
    for(size_t i = 0; i < sentProc.size(); i++)
    {
      int _rank = sentProc[i];
      size_t rankIdx = std::find(this->duplicatedprocs_.begin(), this->duplicatedprocs_.end(), _rank) - this->duplicatedprocs_.begin();
      if(rankIdx == this->duplicatedprocs_.size())
      {
        // TODO: necessary? If `rankIdx` didn't appear in `this->duplicatedprocs_`, we will delete it in the next part
        // new rank in this->duplicatedprocs_, initialize it
        this->duplicatedprocs_.push_back(_rank);
        this->duplicated_points_.emplace_back(std::vector<size_t>());
        this->Nghost_.emplace_back(std::vector<size_t>());
      }
      for(const size_t &pointIdx : sentPoints[i])
      {
        this->duplicated_points_[rankIdx].push_back(pointIdx);
      }
    }
}

/**
 * \author Maor Mizrachi
 * \brief Ensures that the duplicated and ghost arrays contain only the points from/to ranks which are intersecting (sent iff received)
*/
void Voronoi3D::EnsureSymmetry(const std::vector<int> &sentProc, const std::vector<std::vector<int>> &recvProcLists)
{
    for(size_t i = 0; i < this->duplicatedprocs_.size(); i++)
    {
        int _rank =  this->duplicatedprocs_[i];
        bool notAppearingInSent = (std::find(sentProc.begin(), sentProc.end(), _rank) == sentProc.end());
        bool notAppearingInAllRecv = std::all_of(recvProcLists.cbegin(), recvProcLists.cend(), [_rank](const std::vector<int> &recvProcList){return std::find(recvProcList.cbegin(), recvProcList.cend(), _rank) == recvProcList.cend();});
        
        if(notAppearingInSent or notAppearingInAllRecv)
        {
            // not in the intersection, remove the rank
            this->duplicatedprocs_.erase(this->duplicatedprocs_.begin() + i);
            this->duplicated_points_.erase(this->duplicated_points_.begin() + i);
            this->Nghost_.erase(this->Nghost_.begin() + i);
            i--;
        }
    }
}

// void Voronoi3D::InitialExchange(const std::vector<Vector3D> &points, std::vector<int> &sentProc, std::vector<std::vector<size_t>> &sentPoints, const MPI_Comm &comm)
// {
//     const EnvironmentAgent *envAgent = this->pointsManager->getEnvironmentAgent();
//     bool supportsFurthestClosestRanks;
//     std::function<HilbertCurveEnvironmentAgent::DistancesVector(const Vector3D&)> getFurthestClosestRanks;

//     // check if has 'smartAgent' (an agent that can caluclate distances of ranks as well)            
//     const DistributedOctEnvironmentAgent *distribuedOctEnvAgent = dynamic_cast<const DistributedOctEnvironmentAgent*>(envAgent);
//     if(distribuedOctEnvAgent != nullptr)
//     {
//         supportsFurthestClosestRanks = true;
//         getFurthestClosestRanks = [distribuedOctEnvAgent](const Vector3D &point){return distribuedOctEnvAgent->getClosestFurthestPointsByRanks(point);};
//     }
//     const HilbertTreeEnvironmentAgent *hilbertTreeEnvAgent = dynamic_cast<const HilbertTreeEnvironmentAgent*>(envAgent);
//     if(hilbertTreeEnvAgent != nullptr)
//     {
//         supportsFurthestClosestRanks = true;
//         getFurthestClosestRanks = [hilbertTreeEnvAgent](const Vector3D &point){return hilbertTreeEnvAgent->getClosestFurthestPointsByRanks(point);};
//     }

//     if(not supportsFurthestClosestRanks)
//     {
//         return;
//     }
    
//     int rank, size;
//     MPI_Comm_rank(comm, &rank);
//     MPI_Comm_size(comm, &size);

//     size_t counter = 0;

//     for(size_t pointIdx = 0; pointIdx < points.size(); pointIdx++)
//     {
//         bool isBorderPoint = false;
//         for(const size_t &tetraIdx : this->PointTetras_[pointIdx])
//         {
//             const Tetrahedron &tet = this->del_.tetras_[tetraIdx];
//             isBorderPoint = (tet.points[0] >= this->Norg_) or (tet.points[1] >= this->Norg_) or (tet.points[2] >= this->Norg_) or (tet.points[3] >= this->Norg_);
//             if(isBorderPoint)
//             {
//                 break;
//             }
//         }
//         if(!isBorderPoint)
//         {
//             continue;
//         }
//         int closestRank = std::numeric_limits<int>::max();
//         double closestDistance = std::numeric_limits<double>::max();
//         auto distances = getFurthestClosestRanks(points[pointIdx]);
//         for(int _rank = 0; _rank < size; _rank++)
//         {
//             if(_rank == rank)
//             {
//                 continue;
//             }
//             if(distances[_rank].first < closestDistance)
//             {
//                 closestDistance = distances[_rank].first;
//                 closestRank = _rank;
//             }
//         }
//         size_t rankIdx = std::distance(sentProc.begin(), std::find(sentProc.begin(), sentProc.end(), closestRank));
//         if(rankIdx == sentProc.size())
//         {
//             sentProc.push_back(closestRank);
//             sentPoints.emplace_back(std::vector<size_t>());
//         }
//         sentPoints[rankIdx].push_back(pointIdx);
//         counter++;
//     }

//     std::vector<size_t> sendLengths(size, 0);
//     std::vector<std::vector<Vector3D>> toSend;
//     toSend.resize(sentProc.size());
    
//     std::vector<MPI_Request> requests;
//     requests.reserve(4 * sentProc.size()); // heuristic

//     std::vector<size_t> recvLengths(size, 0);

//     for(size_t i = 0; i < sentProc.size(); i++)
//     {
//         int _rank = sentProc[i];
//         sendLengths[_rank] = sentPoints[i].size();
//         toSend[i].reserve(sentPoints[i].size());
//         for(size_t &pointIdx : sentPoints[i])
//         {
//             toSend[i].emplace_back(points[pointIdx]);
//         }
//     }

//     MPI_Alltoall(&sendLengths[0], sizeof(size_t), MPI_BYTE, &recvLengths[0], sizeof(size_t), MPI_BYTE, comm);

//     size_t totalLength = 0; 
//     for(int _rank = 0; _rank < size; _rank++)
//     {
//         if(recvLengths[_rank] > 0)
//         {
//             totalLength += recvLengths[_rank];
//             size_t rankIdx = std::distance(sentProc.begin(), std::find(sentProc.begin(), sentProc.end(), _rank));
//             if(rankIdx != sentProc.size())
//             {
//                 // rank has already been found
//                 continue;
//             }
//             sentProc.push_back(_rank);
//             sentPoints.emplace_back(std::vector<size_t>());
//         }
//     }

//     std::vector<PointData> almostExtraPoints;
//     almostExtraPoints.resize(totalLength);
//     size_t insertedSoFar = 0; 
//     for(const int &_rank : sentProc)
//     {
//         if(recvLengths[_rank] > 0)
//         {
//             // std::cout << "rank " << rank << " is receiving " << recvLengths[_rank] << " from rank " << _rank << ", insertedSoFar is " << insertedSoFar << "(total length: " << totalLength << ")" << std::endl;
//             requests.push_back(MPI_REQUEST_NULL);
//             MPI_Irecv(&almostExtraPoints[insertedSoFar], sizeof(PointData) * recvLengths[_rank], MPI_BYTE, _rank, INITIAL_SENDRECV_TAG, comm, &requests[requests.size() - 1]);
//             size_t dupRankIdx = std::distance(this->duplicatedprocs_.begin(), std::find(this->duplicatedprocs_.begin(), this->duplicatedprocs_.end(), _rank));
//             if(dupRankIdx == this->duplicatedprocs_.size())
//             {
//                 // new rank in this->duplicatedprocs_, initialize it
//                 this->duplicatedprocs_.push_back(_rank);
//                 this->duplicated_points_.emplace_back(std::vector<size_t>());
//                 this->Nghost_.emplace_back(std::vector<size_t>());
//             }
//             for(size_t i = 0; i < recvLengths[_rank]; i++)
//             {
//                 // batchInfo.pointsFromRanks[_rank][i] holds an index of point, but this point will be added to my delaunay, so
//                 // its index there will be this->del_.points_.size() + batchInfo.pointsFromRanks[_rank][i]
//                 this->Nghost_[dupRankIdx].push_back(this->del_.points_.size() + insertedSoFar + i);
//             }
//         }
//     }

//     for(size_t i = 0; i < toSend.size(); i++)
//     {
//         int _rank = sentProc[i];
//         // std::cout << "rank " << rank << " is sending " << toSend[i].size() << " to rank " << _rank << std::endl;
//         requests.push_back(MPI_REQUEST_NULL);
//         MPI_Isend(&toSend[i][0], sizeof(_3DPoint) * toSend[i].size(), MPI_BYTE, _rank, INITIAL_SENDRECV_TAG, comm, &requests[requests.size() - 1]);
//     }

//     if(!requests.empty())
//     {
//         MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
//     }   

//     std::vector<Vector3D> extraPoints;
//     for(const _3DPoint &_point : almostExtraPoints)
//     {
//         extraPoints.emplace_back(Vector3D(_point.x, _point.y, _point.z));
//     }

//     this->del_.BuildExtra(extraPoints);

//     this->R_.resize(this->del_.tetras_.size());
//     std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
//     this->tetra_centers_.resize(this->R_.size());
//     this->bigtet_ = SetPointTetras();
// }


/**
 * \author Maor Mizrachi
 * \brief Sets the ghost points arrays (duplicatedprocs_, duplicated_points_, Nghost_)
*/
void Voronoi3D::SetGhostArray(const std::vector<int> &recvProc, const std::vector<std::vector<size_t>> &recvPoints)
{
    for(size_t i = 0; i < recvProc.size(); i++)
    {
        int _rank = recvProc[i];
        const std::vector<size_t> &receivedFromRank = recvPoints[i];
        size_t rankIdx = std::find(this->duplicatedprocs_.begin(), this->duplicatedprocs_.end(), _rank) - this->duplicatedprocs_.begin();
        if(rankIdx == this->duplicatedprocs_.size())
        {
            // new rank in this->duplicatedprocs_, initialize it
            this->duplicatedprocs_.push_back(_rank);
            this->duplicated_points_.emplace_back(std::vector<size_t>());
            this->Nghost_.emplace_back(std::vector<size_t>());
        }
        for(const size_t &RelativePointIdx : receivedFromRank)
        {
            // batchInfo.pointsFromRanks[_rank][i] holds an index of point, but this point will be added to my delaunay, so
            // its index there will be this->del_.points_.size() + batchInfo.pointsFromRanks[_rank][i]
            this->Nghost_[rankIdx].push_back(this->del_.points_.size() + RelativePointIdx);
        }
    }
}

/**
 * \author Maor Mizrachi
 * \brief Makes load rebalancing if needed, if needed, and initializing the environment agent (the object which is responsible for dividing the space to ranks)
*/
std::vector<Vector3D> Voronoi3D::PrepareToBuildParallel(const std::vector<Vector3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing, bool suppressExchange)
{
    if(this->radiuses.size() < allPoints.size())
    {
        this->radiuses.resize(allPoints.size(), RADIUS_UNINITIALIZED);
    }
    if(this->all_CM.size() < allPoints.size())
    {
        this->all_CM.resize(allPoints.size());
    }
    
    std::chrono::high_resolution_clock::time_point start, end;
    
    assert(this->pointsManager != nullptr);

    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int canDoRebalance = ((not suppressRebalancing) and (indicesToBuild.size() == allPoints.size()))? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &canDoRebalance, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    bool allowRebalance = (canDoRebalance == 1);

    int canDoExchange = (not suppressExchange) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &canDoExchange, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    bool allowExchange = (canDoExchange == 1);

    if(suppressExchange)
        this->allPointsWeights.resize(allPoints.size(), 1.0);

    std::vector<VoronoiPayload> payloads(allPoints.size());
    for(size_t i = 0; i < allPoints.size(); i++)
        payloads[i] = VoronoiPayload(this->radiuses[i], this->all_CM[i]);

    PointsExchangeResult<Vector3D, VoronoiPayload> exchangeResult = this->pointsManager->update(allPoints, allowExchange? allWeights : this->allPointsWeights, payloads, indicesToBuild, allowRebalance, allowExchange);

    this->allMyPoints = std::move(exchangeResult.newPoints);
    this->radiuses.resize(exchangeResult.newPayloads.size());
    this->all_CM.resize(exchangeResult.newPayloads.size());
    for(size_t i = 0; i < exchangeResult.newPayloads.size(); i++)
    {
        this->radiuses[i] = exchangeResult.newPayloads[i].radius;
        this->all_CM[i] = exchangeResult.newPayloads[i].CM;
    }
    this->sentprocs_ = std::move(exchangeResult.sentProcessors);
    this->sentpoints_ = std::move(exchangeResult.sentIndicesToProcessors);
    this->self_index_ = std::move(exchangeResult.indicesToSelf);
    this->allPointsWeights = std::move(exchangeResult.newWeights);

    assert(this->allMyPoints.size() == this->allPointsWeights.size());

    std::vector<Vector3D> new_points;
    this->indicesInAllMyPoints = Tessellation3D::AllPointsMap();
    size_t numOfSelfPoints = this->self_index_.size();

    // this loop determines the points list, and the matching indices of each point (in the points list) to the long points list
    size_t allPointsSize = this->allMyPoints.size();
    for(size_t pointIdx = 0; pointIdx < allPointsSize; pointIdx++)
    {
        if(exchangeResult.participatingIndices.at(pointIdx))
        {
            this->indicesInAllMyPoints[new_points.size()] = pointIdx;
            new_points.push_back(this->allMyPoints[pointIdx]);
        }
    }
    
    this->BuildInitialize(new_points.size());

    return new_points;
}

/**
 * Shuffling or adding more points. The input for this function is a list of points, and a masks list. The mask list is a list of indices, that says
 * for each point in the points list, what's its matching points in the old points list (the current tesselation). If the point is new, the mask should be higher
 * then the current number of points in the tesselation.
*/
void Voronoi3D::PreparePoints(const std::vector<Vector3D> &points, const std::vector<size_t> &mask)
{
    if(points.size() != mask.size())
    {
        UniversalError eo("In Voronoi3D::PreparePoints, mask size is not equal to the points size");
        eo.addEntry("Mask size", mask.size());
        eo.addEntry("Points size", points.size());
        throw eo;
    }
    size_t originalPointsNum = this->allMyPoints.size();
    size_t newPointsNum = points.size();

    std::vector<IndexedVector3D> oldPoints;
    for(size_t i = 0; i < newPointsNum; i++)
    {
        size_t matchingPointIdx = mask[i];
        if(matchingPointIdx < originalPointsNum)
        {
            // this point has a matching old point
            oldPoints.emplace_back(points[i], matchingPointIdx);
        }
    }

    std::vector<double> newRadiuses(newPointsNum, RADIUS_UNINITIALIZED);
    if(!oldPoints.empty())
    {
        OctTree<IndexedVector3D> oldPointsTree(this->ll_, this->ur_, oldPoints);
        newRadiuses = std::vector<double>(newPointsNum);
        for(size_t i = 0; i < newPointsNum; i++)
        {
            size_t matchingPointIdx = mask[i];
            double radius;

            if(matchingPointIdx >= originalPointsNum)
            {
                // the point is a new, but we take its initial radius to be the same as the closest point's radius
                if(newPointsNum > 1)
                {
                    size_t closestPointIdx = oldPointsTree.closestPoint(points[i]).getIndex();
                    radius = this->radiuses.at(closestPointIdx);
                }
                else
                    radius = abs(this->ll_ - this->ur_);
            }
            else
            {
                // old point, copy the radius
                radius = this->radiuses[matchingPointIdx];
            }
            newRadiuses[i] = radius;
        }
    }
    // copy the radiuses into the radiuses array
    this->radiuses = std::move(newRadiuses);
}

void Voronoi3D::UpdatePointsTree(const std::vector<Vector3D> &activePoints)
{    
    auto start = std::chrono::high_resolution_clock::now();
    
    Vector3D width = this->ur_ - this->ll_;
    this->allMyPointsTree = std::make_shared<OctTree<IndexedVector3D>>(IndexedVector3D(this->ll_ - width * 0.00001, std::numeric_limits<size_t>::max()),
                                                                    IndexedVector3D(this->ur_ + width * 0.00001, std::numeric_limits<size_t>::max()));
    size_t allPointsNum = this->allMyPoints.size();
    for(size_t pointIdx = 0; pointIdx < allPointsNum; pointIdx++)
    {
        const Vector3D &point = activePoints[pointIdx];
        this->allMyPointsTree->insert(IndexedVector3D(point.x, point.y, point.z, pointIdx));
    }

    if(this->allMyPoints.size() == activePoints.size())
    {
        // not a real parital build
        this->myPointsTree = this->allMyPointsTree;
    }
    else
    {
        this->myPointsTree = std::make_shared<OctTree<IndexedVector3D>>(IndexedVector3D(this->ll_ - width * 0.00001, std::numeric_limits<size_t>::max()),
                                                                        IndexedVector3D(this->ur_ + width * 0.00001, std::numeric_limits<size_t>::max()));
        for(size_t pointIdx = 0; pointIdx < this->Norg_; pointIdx++)
        {
            const Vector3D &point = this->allMyPoints[this->indicesInAllMyPoints[pointIdx]];
            this->myPointsTree->insert(IndexedVector3D(point.x, point.y, point.z, pointIdx));
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
    {
        std::cout << "Time for tree: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }
}
 
bool Voronoi3D::DidRebalance(void) const
{
    return this->pointsManager->didRebalance();
}

void Voronoi3D::BuildPartiallyParallel(const std::vector<Vector3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing, bool suppressExchange)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(const size_t &idx : indicesToBuild)
    {
        if(idx >= allPoints.size())
        {
            UniversalError eo("BuildParatiallyParallel: Illegal point was given");
            eo.addEntry("Index", idx);
            eo.addEntry("Size", allPoints.size());
            throw eo;
        }
    }

    std::chrono::high_resolution_clock::time_point start, end;
    START_TIMER_PREEMPTIVE("Prepare to build");

    start = std::chrono::high_resolution_clock::now();

    std::vector<Vector3D> activePoints = this->PrepareToBuildParallel(allPoints, allWeights, indicesToBuild, suppressRebalancing, suppressExchange);
    end = std::chrono::high_resolution_clock::now();
    if(rank == 0)
    {
        std::cout << "Time for preparing (hilbert tree + load balancing): " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    start = std::chrono::high_resolution_clock::now();

    std::vector<size_t> order;

    // build delaunay
    if(not activePoints.empty())
    {
        std::pair<Vector3D, Vector3D> bounding_box = std::make_pair(activePoints[0], activePoints[0]);
        for(const Vector3D &point : activePoints)
        {
            bounding_box.first.x = std::min(bounding_box.first.x, point.x);
            bounding_box.second.x = std::max(bounding_box.second.x, point.x);
            bounding_box.first.y = std::min(bounding_box.first.y, point.y);
            bounding_box.second.y = std::max(bounding_box.second.y, point.y);
            bounding_box.first.z = std::min(bounding_box.first.z, point.z);
            bounding_box.second.z = std::max(bounding_box.second.z, point.z);
        }

        // Ensure bounding box covers at least the domain box to prevent
        // degenerate big tets for collinear/coplanar point sets
        bounding_box.first.x = std::min(bounding_box.first.x, this->ll_.x);
        bounding_box.first.y = std::min(bounding_box.first.y, this->ll_.y);
        bounding_box.first.z = std::min(bounding_box.first.z, this->ll_.z);
        bounding_box.second.x = std::max(bounding_box.second.x, this->ur_.x);
        bounding_box.second.y = std::max(bounding_box.second.y, this->ur_.y);
        bounding_box.second.z = std::max(bounding_box.second.z, this->ur_.z);

        // performs internal tesselation:
        // std::cout << "checking duplications..." << std::endl;
        // reportDuplications(new_points);
        order = HilbertOrder3D(activePoints);
        
        // initial build for the points
        this->del_.Build(activePoints, bounding_box.second, bounding_box.first, order);
    }

    // updates the radiuses array of the tetrahedra, as well as the lists for each point what tetras it belongs to
    this->R_.resize(this->del_.tetras_.size());
    this->R_.shrink_to_fit();
    // std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
    this->tetra_centers_.resize(this->R_.size());
    this->tetra_centers_.shrink_to_fit();
    this->bigtet_ = SetPointTetras();

    // MPI_Barrier(MPI_COMM_WORLD);

    end = std::chrono::high_resolution_clock::now();
    if(rank == 0)
    {
        std::cout << "Time for initial build: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    if(this->radiuses.size() < this->Norg_)
    {
        UniversalError eo("Voronoi3D:BuildPartiallyParallel: wrong size of radiuses array");
        eo.addEntry("Rank", rank);
        eo.addEntry("radiuses.size()", this->radiuses.size());
        eo.addEntry("Norg_", this->Norg_);
        throw eo;
    }

    START_TIMER_PREEMPTIVE("Tree construction");

    this->UpdatePointsTree(activePoints);
    
    START_TIMER_PREEMPTIVE("Radiuses calculation");
    start = std::chrono::high_resolution_clock::now();
    this->UpdateRadiuses(activePoints);    
    end = std::chrono::high_resolution_clock::now();

    if(rank == 0)
    {
        std::cout << "Time for radiuses: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    START_TIMER_PREEMPTIVE("Range agent update");

    start = std::chrono::high_resolution_clock::now();
    this->UpdateRangeFinder();
    end = std::chrono::high_resolution_clock::now();

    if(rank == 0)
    {
        std::cout << "Time for range agent: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    START_TIMER_PREEMPTIVE("Bringing ghosts");

    this->BringGhostPointsToBuild(MPI_COMM_WORLD);

    START_TIMER_PREEMPTIVE("Building voronoi");

    start = std::chrono::high_resolution_clock::now();
    
    CM_.resize(del_.points_.size());
    CM_.shrink_to_fit();
    volume_.resize(Norg_);
    volume_.shrink_to_fit();

    if(not activePoints.empty())
    {
        // Create Voronoi
        BuildVoronoi(order);
    }


    // todo: why?
    // std::vector<double>().swap(this->R_);
    // std::vector<tetra_vec>().swap(this->PointTetras_);

    this->UpdateCMs();

    MPI_exchange_data(*this, this->volume_, true);

    end = std::chrono::high_resolution_clock::now();

    if(rank == 0)
    {
        std::cout << "Time for build Voronoi from Delaunay: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    // save the list of the real ghost points
    // this->FilterRealGhostPoints();
}

boost::container::flat_map<size_t, std::pair<rank_t, size_t>> GetGhostInfo(const std::vector<rank_t> &duplicatedProcs, const std::vector<std::vector<size_t>> &duplicatedPoints, const std::vector<std::vector<size_t>> &Nghost, const boost::container::flat_map<size_t, std::pair<rank_t, size_t>> &whereNow)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<std::vector<std::pair<rank_t, size_t>>> toSend(size);

    for(size_t i = 0; i < duplicatedProcs.size(); i++)
    {
        rank_t _rank = duplicatedProcs[i];
        const std::vector<size_t> &dupPoints = duplicatedPoints[i]; // points I sent to `_rank`
        for(size_t j = 0; j < dupPoints.size(); j++)
        {
            size_t pointIdx = dupPoints[j];
            auto it = whereNow.find(pointIdx);
            assert(it != whereNow.end());
            toSend[_rank].push_back(it->second);
        }
    }

    std::vector<std::vector<std::pair<rank_t, size_t>>> received = MPI_Exchange_all_to_all_sparse(toSend, MPI_COMM_WORLD);
    boost::container::flat_map<size_t, std::pair<rank_t, size_t>> ghostsInfo;
    for(size_t i = 0; i < duplicatedProcs.size(); i++)
    {
        rank_t _rank = duplicatedProcs[i];
        const std::vector<size_t> &NghostOfRank = Nghost[i];
        const std::vector<std::pair<rank_t, size_t>> &receivedFromRank = received[_rank];
        assert(NghostOfRank.size() == receivedFromRank.size());
        for(size_t j = 0; j < NghostOfRank.size(); j++)
        {
            size_t oldIndex = NghostOfRank[j];
            ghostsInfo.insert({oldIndex, receivedFromRank[j]});
        }
    }
    return ghostsInfo;
}

boost::container::flat_map<size_t, std::pair<rank_t, size_t>> GetRemoteIndices(const std::vector<size_t> &selfIndex, const std::vector<rank_t> &sentProcs, const std::vector<std::vector<size_t>> &sentPoints)
{    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<size_t> sentPointsNum(size, 0), recvPointsNum(size);
    for(size_t i = 0; i < sentProcs.size(); i++)
    {
        rank_t _rank = sentProcs[i];
        sentPointsNum[_rank] = sentPoints[i].size();
    }
    MPI_Alltoall(sentPointsNum.data(), 1, MPI_UNSIGNED_LONG_LONG, recvPointsNum.data(), 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    // Offsets must follow processesSend order (matching dataExchange output ordering)
    std::vector<size_t> newPointsOffsets(size, 0);
    size_t currentOffset = selfIndex.size();
    for(size_t i = 0; i < sentProcs.size(); i++)
    {
        rank_t _rank = sentProcs[i];
        newPointsOffsets[_rank] = currentOffset;
        currentOffset += recvPointsNum[_rank];
    }

    std::vector<size_t> remoteOffsets(size);
    MPI_Alltoall(newPointsOffsets.data(), 1, MPI_UNSIGNED_LONG_LONG, remoteOffsets.data(), 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

    boost::container::flat_map<size_t, std::pair<rank_t, size_t>> whereNow;

    // first my points
    for(size_t i = 0; i < selfIndex.size(); i++)
    {
        size_t previousIdx = selfIndex[i];
        size_t newIdx = i;
        whereNow[previousIdx] = std::make_pair(rank, newIdx);
    }

    // now other
    for(size_t i = 0; i < sentProcs.size(); i++)
    {
        rank_t _rank = sentProcs[i];
        const std::vector<size_t> &sentPointsOfRank = sentPoints[i];
        for(size_t j = 0; j < sentPointsOfRank.size(); j++)
        {
            size_t previousIdx = sentPointsOfRank[j];
            size_t newIdx = remoteOffsets[_rank] + j;
            whereNow[previousIdx] = std::make_pair(_rank, newIdx);
        }
    }
    return whereNow;
}

void Voronoi3D::MockMesh(void)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    boost::container::flat_map<size_t, std::pair<rank_t, size_t>> whereNow = GetRemoteIndices(this->self_index_, this->sentprocs_, this->sentpoints_); // previous point -> holder (rank + index)
    assert(whereNow.size() == this->Norg_);
    boost::container::flat_map<size_t, std::pair<rank_t, size_t>> ghostsInfo = GetGhostInfo(this->duplicatedprocs_, this->duplicated_points_, this->Nghost_, whereNow); // previous ghost -> holder (rank + index)

    std::vector<Vector3D> delaunayPoints = this->del_.points_;

    std::vector<Vector3D> new_points;
    this->indicesInAllMyPoints = Tessellation3D::AllPointsMap();
    size_t numOfSelfPoints = this->self_index_.size();

    // this loop determines the points list, and the matching indices of each point (in the points list) to the long points list
    size_t allPointsSize = this->allMyPoints.size();
    for(size_t pointIdx = 0; pointIdx < allPointsSize; pointIdx++)
    {
        this->indicesInAllMyPoints[new_points.size()] = pointIdx;
        new_points.push_back(this->allMyPoints[pointIdx]);
    }
    
    // now find what ghosts should be sent
    using GhostAsk = std::pair<rank_t, size_t>; // receiver rank, point index on sender
    std::vector<std::vector<GhostAsk>> askToSend(size); // [sender][(receiver, point)...]
    std::vector<std::vector<Vector3D>> mirrorsToSend(size);
    std::vector<Vector3D> allMirrors;

    auto addAskToSend = [&](rank_t sender, rank_t receiver, size_t pointIndex)
    {
        if(sender == receiver)
        {
            return;
        }
        askToSend[sender].emplace_back(receiver, pointIndex);
    };

    size_t previousN = this->Norg_;
    for(size_t i = 0; i < previousN; i++)
    {
        // looking for previous point `i`
        const auto &[newOwner, newIndex] = whereNow.at(i);
        const std::vector<size_t> &neighbors = this->GetNeighbors(i);
        for(size_t neighbor : neighbors)
        {
            // looking for previous neighbor `neighbor`
            if(neighbor < previousN)
            {
                // original point and neighbor are both local
                const auto &[neighborNewOwner, neighborNewIndex] = whereNow.at(neighbor);
                addAskToSend(neighborNewOwner, newOwner, neighborNewIndex);
            }
            else
            {
                assert(neighbor >= previousN + 4); // 4 for the big tetra
                // neighbor is either a ghost or a mirror
                auto it = ghostsInfo.find(neighbor);
                if(it != ghostsInfo.end())
                {
                    // neighbor is a ghost
                    const auto &[neighborNewOwner, neighborNewIndex] = it->second;
                    addAskToSend(neighborNewOwner, newOwner, neighborNewIndex);
                }
                else
                {
                    // mirror
                    mirrorsToSend[newOwner].push_back(delaunayPoints[neighbor]);
                    // allMirrors.push_back(delaunayPoints[neighbor]);
                }
            }
        }
    }

    // mirrorsToSend = std::vector<std::vector<Vector3D>>(size, allMirrors); // todo: remove
    
    this->BuildInitialize(new_points.size());
    std::vector<size_t> order;
    if(not new_points.empty())
    {
        order = HilbertOrder3D(new_points);
        std::pair<Vector3D, Vector3D> bounding_box = std::make_pair(new_points[0], new_points[0]);
        for(const Vector3D &point : new_points)
        {
            bounding_box.first.x = std::min(bounding_box.first.x, point.x);
            bounding_box.second.x = std::max(bounding_box.second.x, point.x);
            bounding_box.first.y = std::min(bounding_box.first.y, point.y);
            bounding_box.second.y = std::max(bounding_box.second.y, point.y);
            bounding_box.first.z = std::min(bounding_box.first.z, point.z);
            bounding_box.second.z = std::max(bounding_box.second.z, point.z);
        }
        
        bounding_box.first.x = std::min(bounding_box.first.x, this->ll_.x);
        bounding_box.first.y = std::min(bounding_box.first.y, this->ll_.y);
        bounding_box.first.z = std::min(bounding_box.first.z, this->ll_.z);
        bounding_box.second.x = std::max(bounding_box.second.x, this->ur_.x);
        bounding_box.second.y = std::max(bounding_box.second.y, this->ur_.y);
        bounding_box.second.z = std::max(bounding_box.second.z, this->ur_.z);

        this->del_.Build(new_points, bounding_box.second, bounding_box.first, order);
    }
    // updates the radiuses array of the tetrahedra, as well as the lists for each point what tetras it belongs to
    this->R_.resize(this->del_.tetras_.size());
    this->R_.shrink_to_fit();
    // std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
    this->tetra_centers_.resize(this->R_.size());
    this->tetra_centers_.shrink_to_fit();
    this->bigtet_ = SetPointTetras();

#ifdef TIMING
    auto _prof = [&](const char *label, double localTime)
    {
        struct { double val; int rank; } localIn{localTime, rank}, globalMax;
        MPI_Allreduce(&localIn, &globalMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        double globalSum = 0;
        MPI_Allreduce(&localTime, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(rank == 0)
            std::cout << label << ": max=" << globalMax.val << "s (rank " << globalMax.rank
                      << "), avg=" << globalSum / size << "s" << std::endl;
    };
    auto _t0 = std::chrono::high_resolution_clock::now();
    auto _t1 = _t0;
#endif

    this->UpdatePointsTree(new_points);

#ifdef TIMING
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    this->UpdateRadiuses(new_points);
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("UpdateRadiuses", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    this->UpdateRangeFinder();
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("UpdateRangeFinder", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    size_t countMirrors = 0;
    auto _mirrorsResult = MPI_Exchange_sparse_by_rank(mirrorsToSend, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 1);
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("MirrorsExchange", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    for(const auto &[_senderRank, incomingMirrors] : _mirrorsResult)
    {
        if(this->Norg_ > 0)
        {
            this->del_.BuildExtra(incomingMirrors);
            countMirrors += incomingMirrors.size();
        }
    }
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("BuildExtra(mirrors)", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    for(std::vector<GhostAsk> &asksForSender : askToSend)
    {
        if(asksForSender.size() > 1)
        {
            std::sort(asksForSender.begin(), asksForSender.end());
            asksForSender.erase(std::unique(asksForSender.begin(), asksForSender.end()), asksForSender.end());
        }
    }
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("BuildAskToSend", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    std::vector<std::pair<rank_t, std::vector<GhostAsk>>> whatIshouldSend =
        MPI_Exchange_sparse_by_rank(askToSend, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 2);
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("AskToSendExchange", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    std::vector<std::vector<Vector3D>> whatIShouldSendToRanks(size);

    std::vector<rank_t> newDuplicatedProcs;
    std::vector<std::vector<size_t>> newDuplicatedPoints;
    boost::container::flat_map<rank_t, boost::container::flat_set<size_t>> sentToProcessors;
    boost::container::flat_map<rank_t, size_t> procsToIndices;

    for(const auto &[_requestingRank, sendInfo] : whatIshouldSend)
    {
        for(const auto &[_rank, pointIdx] : sendInfo)
        {
            if(_rank == rank)
            {
                continue;
            }
            auto it = procsToIndices.find(_rank);
            size_t idx = newDuplicatedProcs.size();
            if(it == procsToIndices.end())
            {
                newDuplicatedProcs.push_back(_rank);
                newDuplicatedPoints.emplace_back();
                procsToIndices[_rank] = idx;
            }
            else
            {
                idx = it->second;
            }
            std::vector<size_t> &newDuplicatedPointsOfRank = newDuplicatedPoints[idx];

            boost::container::flat_set<size_t> &sentToRank = sentToProcessors[_rank];
            if(sentToRank.insert(pointIdx).second)
            {
                whatIShouldSendToRanks[_rank].push_back(new_points[pointIdx]);
                newDuplicatedPointsOfRank.push_back(pointIdx);
            }
        }
    }
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("BuildSendToRanks", std::chrono::duration<double>(_t1 - _t0).count());
#endif

    // todo: first build delaunay
    size_t currentIndex = this->del_.points_.size();
    std::vector<Vector3D> buildExtra;

    std::vector<std::vector<size_t>> newNghost(newDuplicatedProcs.size());

#ifdef TIMING
    _t0 = std::chrono::high_resolution_clock::now();
#endif
    std::vector<std::pair<rank_t, std::vector<Vector3D>>> receivedByRanks =
        MPI_Exchange_sparse_by_rank(whatIShouldSendToRanks, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 3);
#ifdef TIMING
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("FinalPointsExchange", std::chrono::duration<double>(_t1 - _t0).count());
#endif
    for(const auto &[_rank, receivedFromRank] : receivedByRanks)
    {
        if(receivedFromRank.empty())
        {
            continue;
        }
        assert(_rank != rank); // not me
        size_t idx = std::distance(newDuplicatedProcs.begin(), std::find(newDuplicatedProcs.begin(), newDuplicatedProcs.end(), _rank));
        if(idx == newDuplicatedProcs.size())
        {
            UniversalError eo("Voronoi3D::MockMesh: received from an unexpected rank");
            eo.addEntry("My rank", rank);
            eo.addEntry("Received from rank", _rank);
            throw eo;
        }
        assert(idx != newDuplicatedProcs.size());
        std::vector<size_t> &NghostOfRank = newNghost[idx];

        for(const Vector3D &point : receivedFromRank)
        {
            buildExtra.push_back(point);
            NghostOfRank.push_back(currentIndex);
            currentIndex++;
        }
    }

    this->del_.BuildExtra(buildExtra);
    this->R_.resize(this->del_.tetras_.size());
    this->R_.shrink_to_fit();
    std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
    this->tetra_centers_.resize(this->R_.size());
    this->tetra_centers_.shrink_to_fit();
    this->bigtet_ = this->SetPointTetras();

    this->duplicatedprocs_ = std::move(newDuplicatedProcs);
    this->duplicated_points_ = std::move(newDuplicatedPoints);
    this->Nghost_ = std::move(newNghost);

    CM_.resize(del_.points_.size());
    CM_.shrink_to_fit();
    volume_.resize(Norg_);
    volume_.shrink_to_fit();

    if(not new_points.empty())
    {
        // Create Voronoi
        this->BuildVoronoi(order);
    }

    this->UpdateCMs();
}

void Voronoi3D::SetLoadBalancer(std::shared_ptr<LoadBalancer<Vector3D>> loadBalancer)
{
    assert(this->pointsManager != nullptr);
    this->pointsManager->setLoadBalancer(loadBalancer);

    std::vector<Vector3D> previousPoints = this->allMyPoints;
    std::vector<size_t> allIndices(previousPoints.size());
    std::iota(allIndices.begin(), allIndices.end(), 0);

    std::vector<VoronoiPayload> payloads(previousPoints.size());
    for(size_t i = 0; i < previousPoints.size(); i++)
        payloads[i] = VoronoiPayload(this->radiuses[i], this->all_CM[i]);

    PointsExchangeResult<Vector3D, VoronoiPayload> exchangeResult = this->pointsManager->update(previousPoints, this->allPointsWeights, payloads, allIndices, false);

    this->allMyPoints = std::move(exchangeResult.newPoints);
    this->radiuses.resize(exchangeResult.newPayloads.size());
    this->all_CM.resize(exchangeResult.newPayloads.size());
    for(size_t i = 0; i < exchangeResult.newPayloads.size(); i++)
    {
        this->radiuses[i] = exchangeResult.newPayloads[i].radius;
        this->all_CM[i] = exchangeResult.newPayloads[i].CM;
    }
    this->sentprocs_ = std::move(exchangeResult.sentProcessors);
    this->sentpoints_ = std::move(exchangeResult.sentIndicesToProcessors);
    this->self_index_ = std::move(exchangeResult.indicesToSelf);
    this->allPointsWeights = std::move(exchangeResult.newWeights);

    this->MockMesh();
}

void Voronoi3D::PresetLoadBalancer(std::shared_ptr<LoadBalancer<Vector3D>> loadBalancer)
{
    assert(this->pointsManager != nullptr);
    this->pointsManager->setLoadBalancer(loadBalancer);
}

void Voronoi3D::Rebalance(const std::vector<double> &weights)
{
    if(this->pointsManager == nullptr)
    {
        UniversalError eo("Voronoi3D::Rebalance: pointsManager is nullptr");
        throw eo;
    }

    if(weights.size() != this->allMyPoints.size())
    {
        UniversalError eo("Voronoi3D::Rebalance: weights.size() != allMyPoints.size()");
        eo.addEntry("weights.size()", weights.size());
        eo.addEntry("allMyPoints.size()", this->allMyPoints.size());
        throw eo;
    }
    
    std::vector<Vector3D> previousPoints = this->allMyPoints;
    std::vector<size_t> allIndices(previousPoints.size());
    std::iota(allIndices.begin(), allIndices.end(), 0);
    std::vector<VoronoiPayload> payloads(previousPoints.size());
    for(size_t i = 0; i < previousPoints.size(); i++)
        payloads[i] = VoronoiPayload(this->radiuses[i], this->all_CM[i]);

    PointsExchangeResult<Vector3D, VoronoiPayload> exchangeResult = this->pointsManager->update(previousPoints, weights, payloads, allIndices);

    this->allMyPoints = std::move(exchangeResult.newPoints);
    this->radiuses.resize(exchangeResult.newPayloads.size());
    this->all_CM.resize(exchangeResult.newPayloads.size());
    for(size_t i = 0; i < exchangeResult.newPayloads.size(); i++)
    {
        this->radiuses[i] = exchangeResult.newPayloads[i].radius;
        this->all_CM[i] = exchangeResult.newPayloads[i].CM;
    }
    this->sentprocs_ = std::move(exchangeResult.sentProcessors);
    this->sentpoints_ = std::move(exchangeResult.sentIndicesToProcessors);
    this->self_index_ = std::move(exchangeResult.indicesToSelf);
    this->allPointsWeights = std::move(exchangeResult.newWeights);

    this->MockMesh();
}

#endif // RICH_MPI

/**
 * \author Maor Mizrachi
 * \brief Gets a point, its radius, a box and the normals to the box's faces, and returns the faces indices that the sphere (around `point`, in the given `radius`) intersects
*/
std::vector<size_t> CheckToMirror(const Sphere<Vector3D> &sphere, const std::vector<Face> &box, const std::vector<Vector3D> &normals)
{
    std::vector<size_t> facesItCuts;
    // std::cout << "point = " << point << ", radius = " << radius << std::endl;
    for(size_t i = 0; i < box.size(); i++)
    {
        // check for intersecting the sphere with radius `radius` around `point`, with the `i`th face of `box`
        if(FaceSphereIntersections(box[i], sphere, normals[i]))
        {
            // intersects! mirror the point
            facesItCuts.push_back(i);
        }
    }
    return facesItCuts;
}

void Voronoi3D::UpdateCMs(void)
{
    // first, calculate CM for active local points
    this->CalcAllCM(); // Now this->CM_ calculates correct CM for all active points, and maybe for more
    size_t FaceNeighborsSize = FaceNeighbors_.size();
    for(std::size_t i = 0; i < FaceNeighborsSize; ++i)
    {
        if(this->BoundaryFace(i))
        {
            this->CalcRigidCM(i);
        }
    }

    this->SyncPartialBuildData(this->CM_, this->all_CM);
}

void Voronoi3D::UpdateRadiuses(const std::vector<Vector3D> &points)
{
    // use an oct tree to fast calculate the distance to closest point
    OctTree<Vector3D> myOctTree(this->ll_, this->ur_, this->allMyPoints.begin(), this->allMyPoints.end());
    size_t const N = this->indicesInAllMyPoints.size();
    for(const std::pair<size_t, size_t> &indices : this->indicesInAllMyPoints)
    {
        const size_t &pointIndexOnBuild = indices.first;
        const Vector3D &point = points[pointIndexOnBuild];
        size_t pointIndexAmongAll = indices.second; 
        if(this->radiuses[pointIndexAmongAll] <= 0)
        {
            // point does not have a radius from a previous timestep. Initialize a radius
            if(N > 1)
                this->radiuses[pointIndexAmongAll] = this->allMyPointsTree->closestPointDistance(point, false); // todo second closest
            else
                this->radiuses[pointIndexAmongAll] = abs(this->ll_ - this->ur_);
        }
    }
}

void Voronoi3D::UpdateRangeFinder()
{
    this->rangeFinder = std::make_shared<OctTreeFinder>(this->allMyPointsTree.get(), this->allMyPoints);
    // if(this->rangeFinder.get() == nullptr)
    // {
    //     //BruteForceFinder rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    //     //RangeTreeFinder rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    //     this->rangeFinder = std::make_shared<RangeFinder>(this->allMyPoints.begin(), this->allMyPoints.end(), this->ll_, this->ur_);
    //     //KDTreeFinder rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_, this->ll_, this->ur_);
    //     //GroupRangeTreeFinder<256> rangeFinder(this->del_.points_.begin(), this->del_.points_.begin() + this->Norg_);
    // }
    // else
    // {
    //     // there are three steps:
    //     /*
    //     1. test which points from the  last timestep that might have changed.
    //     These are the active points from last time step (the points in `indicesInAllMyPoints`),
    //     and the REAL ghost points.
    //     */
    //      // todo: this part should not happen here (since this function is running after the exchange)
    //     std::vector<size_t> pointsToRemove = this->indicesInAllMyPoints;
    //     for(const std::vector<size_t> &ghostPoints : this->real_duplicated_points)
    //     {
    //         pointsToRemove.insert(pointsToRemove.end(), ghostPoints.cbegin(), ghostPoints.cend());
    //     }
    //     /*
    //     2. There are points that I owned in the last timestep, but I don't own anymore in this one. Remove them. 
    //     */
        
    //     /*
    //     3. There are points that I did not own in the last timestep, but I do own now. Add them.
    //     */

    //     std::vector<Vector3D> pointsToAdd;
    //     // todo
    //     this->rangeFinder->replacePoints();
    // }
}

/**
 * \author Maor Mizrachi
 * \brief Creates a batch for a cycle (iteration) in the ghost points bringing loop
*/
std::pair<std::vector<SmallRangeQueryData>, std::vector<BigRangeQueryData>> Voronoi3D::CreateBatches(boost::container::flat_set<size_t> &smallPoints, boost::container::flat_set<size_t> &largePoints, const boost::container::flat_map<size_t, size_t> &firstLargeIteration, std::vector<double> &currentRadiuses, size_t iterations)
{
    std::vector<SmallRangeQueryData> smallQueries;
    std::vector<BigRangeQueryData> bigQueries;
    boost::container::flat_set<size_t> tetraToCancel;

    if(iterations == 1)
    {        
        // at first iteration, run an initial query, all the points are small
        for(const size_t &pointIdx : smallPoints)
        {
            const Vector3D &point = this->del_.points_[pointIdx];
            smallQueries.emplace_back();
            SmallRangeQueryData &query = smallQueries.back();
            query.pointIdx = pointIdx;
            query.center = {point.x, point.y, point.z};
            query.radius = currentRadiuses[pointIdx];
            query.maxPointsToGet = RANGE_MAX_POINTS_TO_GET + 1;

            if(currentRadiuses[pointIdx] <= 0)
            {
                UniversalError eo("Radius for a certain point is <= 0 (in 'Voronoi3D::CreateBatches')");
                eo.addEntry("Point Index", pointIdx);
                eo.addEntry("Its radius", currentRadiuses[pointIdx]);
                throw eo;
            }
        }
    }
    else
    {
        // treat large points
        for(const size_t &pointIdx : largePoints)
        {
            const Vector3D &point = this->del_.points_[pointIdx];

            // std::cout << "firstLargeIteration.at(pointIdx) is " << firstLargeIteration.at(pointIdx) << std::endl;
            bool askOnlyClose = (iterations == firstLargeIteration.at(pointIdx)); // on the first iteration as large, ask only the close ranks

            for(const size_t &tetraIdx : this->PointTetras_[pointIdx])
            {
                if(not this->del_.tetras_[tetraIdx].checkBig)
                {
                    continue; // tetra does not need to be checked
                }
                const Vector3D &center = this->tetra_centers_[tetraIdx];
                double radius = this->GetRadius(tetraIdx);
                // from each big tetrahedron, ask each one of the intersecting ranks to give us the closest point it has to our point
                // for a large point queries, if the iteration number is 2, we ask only the near ranks to give their closest point.
                // From the 3rd iteration, we ask all the intersecting ranks to give their closest point.
                if(not askOnlyClose)
                {
                    // do not cancel the tetra if we ask only close
                    tetraToCancel.insert(tetraIdx);
                }
                bigQueries.emplace_back();
                BigRangeQueryData &query = bigQueries.back();
                query.pointIdx = pointIdx;
                query.center = {center.x, center.y, center.z};
                query.radius = radius;
                query.originalPoint = {point.x, point.y, point.z};
                query.askOnlyClose = askOnlyClose;
                // add the tetra to the list of tetrahedra to clear (mark as 'not new')
            }
        }

        // treat small points
        for(const size_t &pointIdx : smallPoints)
        {
            // submit one query whichf is a union of the others
            const Vector3D &point = this->del_.points_[pointIdx];
            double radius = currentRadiuses[pointIdx];
           
            if(currentRadiuses[pointIdx] <= 0)
            {
                UniversalError eo("Radius for a certain point is <= 0 (in 'Voronoi3D::CreateBatches')");
                eo.addEntry("Point Index", pointIdx);
                eo.addEntry("Point", point);
                eo.addEntry("Current Radius", currentRadiuses[pointIdx]);
                throw eo;
            }

            smallQueries.emplace_back();
            SmallRangeQueryData &query = smallQueries.back();
            query.pointIdx = pointIdx;
            query.center = {point.x, point.y, point.z};
            query.radius = radius;
            query.maxPointsToGet = RANGE_MAX_POINTS_TO_GET + 1;
        }
    }

    for(const size_t &tetraIdx : tetraToCancel)
    {
        this->del_.tetras_[tetraIdx].checkBig = false;
    }
    return {smallQueries, bigQueries};
}

/**
 * \author Maor Mizrachi
 * \brief Gets a list of query, and tests for creating mirror points. In the end of this procedure, `mirroredPoints` contains pairs of <faceIdx, pointIdx>, of points that should be mirrored, in relative to which faces
*/
template<typename QueryDataType>
std::vector<std::pair<size_t, size_t>> MirrorPoints(const std::vector<QueryDataType> &queries, const std::vector<Face> &box, const std::vector<Vector3D> &normals)
{   
    static_assert(std::is_convertible<QueryDataType*, RangeQueryData*>::value, "MirrorPoints: QueryDataType must inherit 'RangeQueryData'");

    std::vector<std::pair<size_t, size_t>> mirroredPoints;
    for(const QueryDataType &query : queries)
    {
        // check for mirroring:
        Sphere<Vector3D> sphere(Vector3D(query.center), query.radius);
        size_t pointIdx = query.pointIdx;

        std::vector<size_t> facesItCuts = CheckToMirror(sphere, box, normals);

        for(const size_t &faceIdx : facesItCuts)
        {
            mirroredPoints.push_back(std::make_pair(faceIdx, pointIdx));
        }
    }
    return mirroredPoints;
}

void Voronoi3D::BringSelfGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                        BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                        boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                                        boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints,
                                        RangeFinder::_set<size_t> &selfIgnorePoints)
{
    std::chrono::high_resolution_clock::time_point start1, end1, start2, end2, start3, end3;

    std::vector<Vector3D> newPoints;
    {
        start1 = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<size_t>> selfSmallQueriesAnswers = smallRangeAgent.selfBatchAnswer(smallQueries, selfIgnorePoints);
        size_t i = 0;
        for(const std::vector<size_t> &newSmallQueriesPoints : selfSmallQueriesAnswers)
        {
            const SmallRangeQueryData &query = smallQueries[i];
            for(const size_t &pointIdxInAll : newSmallQueriesPoints)
            {
                this->indicesInAllMyPoints[this->del_.points_.size() + newPoints.size()] = pointIdxInAll;
                newPoints.push_back(this->allMyPoints[pointIdxInAll]);
            }
            numOfResultsForSmallPoints[query.pointIdx] = newSmallQueriesPoints.size();
            i++;
        }
        end1 = std::chrono::high_resolution_clock::now();
    }
    {
        start2 = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<size_t>> selfBigQueriesAnswers = bigRangeAgent.selfBatchAnswer(bigQueries, selfIgnorePoints);
        size_t i = 0;
        for(const std::vector<size_t> &newBigQueriesPoints : selfBigQueriesAnswers)
        {
            const BigRangeQueryData &query = bigQueries[i];
            for(const size_t &pointIdxInAll : newBigQueriesPoints)
            {
                this->indicesInAllMyPoints[this->del_.points_.size() + newPoints.size()] = pointIdxInAll;
                newPoints.push_back(this->allMyPoints[pointIdxInAll]);
            }
            numOfResultsForBigPoints[query.pointIdx] = newBigQueriesPoints.size();
            i++;
        }
        end2 = std::chrono::high_resolution_clock::now();
    }

    start3 = std::chrono::high_resolution_clock::now();
    this->del_.BuildExtra(newPoints);
    end3 = std::chrono::high_resolution_clock::now();
    #ifdef TIMING
        #ifdef RICH_MPI
            rank_t rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if(rank == 0)
        #endif // RICH_MPI
            {
                std::cout << "Time for small: " << std::chrono::duration<double>(end1 - start1).count() << " seconds" <<
                    ", for large: " << std::chrono::duration<double>(end2 - start2).count() << " seconds" <<
                    ", and Delaunay construction time: " << std::chrono::duration<double>(end3 - start3).count() << " seconds" << std::endl;
            }
    #endif // TIMING
}

#ifdef RICH_MPI
    void Voronoi3D::BringRemoteGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                        BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                        boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                                        boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints)
    {
        std::chrono::high_resolution_clock::time_point start1, end1, start2, end2;

        std::vector<Vector3D> newPoints;
        // large points queries
        {
            start1 = std::chrono::high_resolution_clock::now();
            QueryBatchInfo<BigRangeQueryData, Vector3D> bigBatchInfo = bigRangeAgent.runBatch(bigQueries);
            newPoints.reserve(bigBatchInfo.result.size());
            for(const QueryInfo<BigRangeQueryData, Vector3D> &ans : bigBatchInfo.queriesAnswers)
            {
                numOfResultsForBigPoints[ans.data.pointIdx] += ans.finalResults.size();
            }
            newPoints.insert(newPoints.end(), bigBatchInfo.result.begin(), bigBatchInfo.result.end());
            this->SetGhostArray(bigRangeAgent.getRecvProc(), bigRangeAgent.getRecvPoints());
            this->del_.BuildExtra(newPoints);
            end1 = std::chrono::high_resolution_clock::now();
        }
        // small points queries
        {        
            start2 = std::chrono::high_resolution_clock::now();
            QueryBatchInfo<SmallRangeQueryData, Vector3D> smallBatchInfo = smallRangeAgent.runBatch(smallQueries);
            for(const QueryInfo<SmallRangeQueryData, Vector3D> &ans : smallBatchInfo.queriesAnswers)
            {
                numOfResultsForSmallPoints[ans.data.pointIdx] += ans.finalResults.size();
            }
            newPoints.clear();
            newPoints.reserve(smallBatchInfo.result.size());
            newPoints.insert(newPoints.end(), smallBatchInfo.result.begin(), smallBatchInfo.result.end());
            this->SetGhostArray(smallRangeAgent.getRecvProc(), smallRangeAgent.getRecvPoints());
            this->del_.BuildExtra(newPoints);
            end2 = std::chrono::high_resolution_clock::now();
        }    

        #ifdef TIMING
        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0)
        {
            std::cout << "Time for small: " << std::chrono::duration<double>(end2 - start2).count() << " seconds" <<
                ", for large: " << std::chrono::duration<double>(end1 - start1).count() << " seconds" << std::endl;
        }
        #endif // TIMING
    }
#endif // RICH_MPI

/**
 * \author Maor Mizrachi
 * \brief Calculates the points for next iteration, and determines each one's type (small or big)
*/
std::pair<boost::container::flat_set<size_t>, boost::container::flat_set<size_t>>
Voronoi3D::DetermineNextIterationPoints(size_t iterations,
                                            boost::container::flat_map<size_t, size_t> &firstLargeIteration,
                                            std::vector<double> &currentRadiuses,
                                            const boost::container::flat_map<size_t, size_t> &resultOfSmallPoints,
                                            const boost::container::flat_map<size_t, size_t> &resultOfBigPoints)
{
    boost::container::flat_set<size_t> newSmallPoints, newLargePoints;
    
    // small
    for(const std::pair<size_t, size_t> &pointIdxResult : resultOfSmallPoints)
    {
        const size_t &pointIdxInBuild = pointIdxResult.first;
        size_t pointIdxInAllPoints = this->indicesInAllMyPoints[pointIdxInBuild];
        const size_t &resultSize = pointIdxResult.second;

        // small query                
        if(resultSize > RANGE_MAX_POINTS_TO_GET)
        {
            // the result is too big, we should consider this point as large
            newLargePoints.insert(pointIdxInBuild);
            firstLargeIteration[pointIdxInBuild] = iterations + 1; // the first large iteration for `pointIdx` is the next one
            this->radiuses[pointIdxInAllPoints] = LARGE_POINTS_SHRINK_RADIUS_RATIO * currentRadiuses[pointIdxInBuild];
        }
        else
        {
            double maxRadius = this->GetMaxRadius(pointIdxInBuild);

            if(currentRadiuses[pointIdxInBuild] < 2 * maxRadius)
            {
                // point is not yet done!
                currentRadiuses[pointIdxInBuild] *= RADIUSES_GROWING_FACTOR; // increase radius by 'RADIUSES_GROWING_FACTOR'
                newSmallPoints.insert(pointIdxInBuild);
            }
            else
            {
                // point is finished, set a radius for next iteration
                this->radiuses[pointIdxInAllPoints] = RADIUSES_GROWING_FACTOR * (2 * maxRadius);
            }
        }
    }

    // big
    for(const std::pair<size_t, size_t> &pointIdxResult : resultOfBigPoints)
    {
        const size_t &pointIdx = pointIdxResult.first;
        const size_t &resultSize = pointIdxResult.second;

        // query is large, check if it returned non empty. If yes, we are not yet done
        if((iterations == firstLargeIteration.at(pointIdx)) or (resultSize > 0))
        {
            newLargePoints.insert(pointIdx);
        }
    }

    return std::pair(newSmallPoints, newLargePoints);
}

/**
 * \author Maor Mizrachi
 * \brief The algorithm follows arepro paper (https://www.mpa-garching.mpg.de/~volker/arepo/arepo_paper.pdf), section 2.4.
*/
#ifdef RICH_MPI
    void Voronoi3D::BringGhostPointsToBuild(const MPI_Comm &comm)
#else // RICH_MPI
    void Voronoi3D::BringGhostPointsToBuild()
#endif // RICH_MPI
{
    int rank = 0, size = 1;
    #ifdef RICH_MPI
    const bool serialMode = (this->pointsManager == nullptr); // TODO: fix - wrong!
    if(!serialMode)
    {
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
    }
    else
    {
        std::cout << "In serial mode" << std::endl;
    }
    #endif // RICH_MPI

    std::vector<Face> box;
    std::vector<Vector3D> normals;
    this->InitialBoxBuild(box, normals);
    
    boost::container::flat_set<size_t> smallPoints; // indices of 'small' points
    boost::container::flat_set<size_t> largePoints; // indices of 'large' points
    boost::container::flat_map<size_t, size_t> firstLargeIteration;

    size_t totalBigQueries = 0;
    size_t totalSmallQueries = 0;

    std::vector<double> currentRadiuses(this->Norg_, RADIUS_UNINITIALIZED);
    for(const std::pair<size_t, size_t> &indices : this->indicesInAllMyPoints)
    {
        size_t pointIndexInBuild = indices.first;
        size_t pointIndexInAll = indices.second;
        smallPoints.insert(pointIndexInBuild);
        // largePoints.insert(pointIndexInBuild);
        // firstLargeIteration.insert({pointIndexInBuild, 0});

        if(pointIndexInBuild >= this->Norg_)
        {
            std::cout << "Error: given point index " << pointIndexInBuild << ", while built only with " << this->Norg_ << " points" << std::endl;
        }
        currentRadiuses[pointIndexInBuild] = this->radiuses[pointIndexInAll];
    }

    #ifdef RICH_MPI
    std::vector<int> alreadyRecvProcs;
    std::optional<SentPointsContainer> optPointsContainer;
    if (!serialMode)
    {
        auto [ghostPointsFromLastBuild, alreadySentPoints1, alreadyRecvPoints1] = this->InitialGhostPointsExchange(comm);
        std::vector<int> alreadySentProcs;
        std::vector<std::vector<size_t>> alreadySentPoints2, alreadyRecvPoints2;
        for(rank_t _rank = 0; _rank < size; _rank++)
        {
            if(alreadySentPoints1[_rank].size() > 0)
            {
                alreadySentProcs.push_back(_rank);
                alreadySentPoints2.emplace_back(std::move(alreadySentPoints1[_rank]));
            }
            if(alreadyRecvPoints1[_rank].size() > 0)
            {
                alreadyRecvProcs.push_back(_rank);
                alreadyRecvPoints2.emplace_back(std::move(alreadyRecvPoints1[_rank]));
            }
        }
        this->SetGhostArray(alreadyRecvProcs, alreadyRecvPoints2);
        this->del_.BuildExtra(ghostPointsFromLastBuild);
        this->R_.resize(this->del_.tetras_.size(), RADIUS_UNINITIALIZED);
        this->R_.shrink_to_fit();
        this->tetra_centers_.resize(this->R_.size());
        this->tetra_centers_.shrink_to_fit();
        this->bigtet_ = SetPointTetras();
        optPointsContainer.emplace(alreadySentProcs, alreadySentPoints2);
    }
    else
    {
        optPointsContainer.emplace(); // empty container for serial mode
    }
    SentPointsContainer &pointsContainer = *optPointsContainer;

    // In serial mode: envAgent is null (safe -- only used by talk agent for remote queries,
    // which never fire since sendToSelf=false and size=1).
    const std::shared_ptr<EnvironmentAgent<Vector3D>> envAgent = serialMode ? nullptr : this->pointsManager->getEnvironmentAgent();
    const MPI_Comm &agentComm = serialMode ? MPI_COMM_SELF : comm;
    BigRangeAgent bigRangeAgent(this->rangeFinder.get(), envAgent, pointsContainer, agentComm);
    SmallRangeAgent smallRangeAgent(this->rangeFinder.get(), envAgent, pointsContainer, agentComm);
#else // RICH_MPI
    BigRangeAgent bigRangeAgent(this->rangeFinder.get());
    SmallRangeAgent smallRangeAgent(this->rangeFinder.get());
#endif // RICH_MPI


    #ifdef RICH_MPI
        MPI_Request finishedReq;
        int I_finished = 0;
    #endif // RICH_MPI
    int finished;

    // Use a set for O(1) lookup instead of O(n) linear search in the loop
    std::set<std::pair<size_t, size_t>> allMirroredSet;

    size_t iterations = 0;

    bool considerOwnPoints = (size == 1) or (this->Norg_ != this->allMyPoints.size()); // there are points which we ignore in this step, so we have, in the range searching, communicate with us as well
    RangeFinder::_set<size_t> selfIgnorePoints;
    for(const std::pair<size_t, size_t> &indices : this->indicesInAllMyPoints)
    {
        selfIgnorePoints.insert(indices.second);
    }

    auto start = std::chrono::high_resolution_clock::now();
    
    START_TIMER_PREEMPTIVE("Main Loop");
    size_t total_new_points = 0;
    while(true) // loop is not really infinite (has 'break')
    {
        auto start_iter = std::chrono::high_resolution_clock::now();
        
        boost::container::flat_map<size_t, size_t> numOfResultsForSmallPoints;
        boost::container::flat_map<size_t, size_t> numOfResultsForBigPoints;
                
        size_t smallPointsNum = smallPoints.size();
        size_t largePointsNum = largePoints.size();
        #ifdef RICH_MPI
        if (!serialMode)
        {
            MPI_Reduce((rank == 0)? MPI_IN_PLACE : &smallPointsNum, &smallPointsNum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
            MPI_Reduce((rank == 0)? MPI_IN_PLACE : &largePointsNum, &largePointsNum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        }
        #endif // RICH_MPI
        iterations++;
        size_t averageGP = this->del_.points_.size();
        #ifdef RICH_MPI
        if (!serialMode)
            MPI_Reduce((rank == 0)? MPI_IN_PLACE : &averageGP, &averageGP, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
        #endif // RICH_MPI

        averageGP /= size;
        totalBigQueries += largePointsNum;
        totalSmallQueries += smallPointsNum;

        if(rank == 0) std::cout << "iteration " << iterations << " (" << smallPointsNum << " small points, " << largePointsNum << " large points, average ghost points: " << averageGP << ")" << std::endl;

        auto start1 = std::chrono::high_resolution_clock::now();
        auto [smallQueries, bigQueries] = this->CreateBatches(smallPoints, largePoints, firstLargeIteration, currentRadiuses, iterations);
        std::vector<std::pair<size_t, size_t>> mirroredPoints = MirrorPoints(smallQueries, box, normals);
        std::vector<std::pair<size_t, size_t>> moreMirroredPoints = MirrorPoints(bigQueries, box, normals);
        mirroredPoints.insert(mirroredPoints.end(), moreMirroredPoints.begin(), moreMirroredPoints.end());
        auto end1 = std::chrono::high_resolution_clock::now();
        
        #ifdef TIMING
        if(rank == 0)
        {
            std::cout << "Time for creating batches and mirrors: " << std::chrono::duration<double>(end1 - start1).count() << " seconds" << std::endl;
        }
        #endif // TIMING

        #ifdef RICH_MPI
        if (!serialMode)
        {

            I_finished = (smallQueries.empty() and bigQueries.empty())? 1 : 0;
            MPI_Iallreduce(&I_finished, &finished, 1, MPI_INT, MPI_SUM, comm, &finishedReq);
        }
        else
        {
            finished = (smallQueries.empty() and bigQueries.empty())? 1 : 0;
        }
        #else // RICH_MPI
            finished = (smallQueries.empty() and bigQueries.empty())? 1 : 0;
        #endif // RICH_MPI

        if(considerOwnPoints)
        {
            this->BringSelfGhostPoints(bigQueries, smallQueries, bigRangeAgent, smallRangeAgent, numOfResultsForBigPoints, numOfResultsForSmallPoints, selfIgnorePoints);
        }

        #ifdef RICH_MPI
        if (!serialMode)
            this->BringRemoteGhostPoints(bigQueries, smallQueries, bigRangeAgent, smallRangeAgent, numOfResultsForBigPoints, numOfResultsForSmallPoints);
        #endif // RICH_MPI

        auto start2 = std::chrono::high_resolution_clock::now();
        std::vector<Vector3D> newPoints;
        newPoints.reserve(mirroredPoints.size());
        for(const std::pair<size_t, size_t> &pairFacePoint : mirroredPoints)
        {
            // check if we have already mirrored this point with this face - O(log n) lookup
            if(allMirroredSet.find(pairFacePoint) == allMirroredSet.end())
            {
                allMirroredSet.insert(pairFacePoint); // remember we mirrored this point with this face
                newPoints.push_back(MirrorPoint(box[pairFacePoint.first], this->del_.points_[pairFacePoint.second]));
            }
        }
        this->del_.BuildExtra(newPoints);
        auto end2 = std::chrono::high_resolution_clock::now();

        #ifdef TIMING
        if(rank == 0)
        {
            std::cout << "Time for building mirrors: " << std::chrono::duration<double>(end2 - start2).count() << " seconds" << std::endl;
        }
        #endif // TIMING

        auto start3 = std::chrono::high_resolution_clock::now();

        this->R_.resize(this->del_.tetras_.size(), RADIUS_UNINITIALIZED);
        this->R_.shrink_to_fit();
        // std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
        this->tetra_centers_.resize(this->R_.size());
        this->tetra_centers_.shrink_to_fit();
        this->bigtet_ = SetPointTetras();

        auto end3 = std::chrono::high_resolution_clock::now();
        
        #ifdef TIMING
        if(rank == 0)
        {
            std::cout << "Time for calculating tetras: " << std::chrono::duration<double>(end3 - start3).count() << " seconds" << std::endl;
        }
        #endif // TIMING

        size_t new_points = 0;
        #ifdef RICH_MPI        
        if (!serialMode)
        {
            size_t new_points_until_now = std::accumulate(this->Nghost_.cbegin(), this->Nghost_.cend(), 0, [](const size_t &a, const std::vector<size_t> &b){return a + b.size();});
            new_points = new_points_until_now - total_new_points;
            total_new_points = new_points_until_now;
            MPI_Allreduce(MPI_IN_PLACE, &new_points, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        }
        #else // RICH_MPI
            new_points = newPoints.size();
        #endif // RICH_MPI

        auto end_iter = std::chrono::high_resolution_clock::now();

        if(rank == 0) std::cout << "added new points: " << new_points << ", total time: " << std::chrono::duration<double>(end_iter - start_iter).count() << std::endl;

        std::tie(smallPoints, largePoints) = this->DetermineNextIterationPoints(iterations, firstLargeIteration, currentRadiuses, numOfResultsForSmallPoints, numOfResultsForBigPoints);

        // #ifdef RICH_MPI
        //     std::tie(smallPoints, largePoints) = this->DetermineNextIterationPoints(iterations, firstLargeIteration, currentRadiuses, selfSmallQueriesAnswers, selfBigQueriesAnswers, smallBatchInfo.queriesAnswers, bigBatchInfo.queriesAnswers);
        // #else // RICH_MPI
        // #endif // RICH_MPI

        #ifdef RICH_MPI
        if (!serialMode)
            MPI_Wait(&finishedReq, MPI_STATUS_IGNORE);
        #endif // RICH_MPI

        if(finished == size)
        {
            break;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    if(rank == 0)
    {
        std::cout << "Time for bringing ghosts: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    if(rank == 0)
    {
        std::cout << "Total small queries: " << totalSmallQueries << ", total big queries: " << totalBigQueries << std::endl;
    }
    
    START_TIMER_PREEMPTIVE("Organizing sent/recv and ghosts arrays");
    #ifdef RICH_MPI   
    if (!serialMode)
    {     
        const std::vector<std::vector<size_t>> &sentPoints = pointsContainer.getSentData();
        const std::vector<int> &sentProc = pointsContainer.getSentProc();

        // calculate this->duplicated_points_
        this->UpdateDuplicatedPoints(sentProc, sentPoints);
        // remove whomever that does not appear both in my sent vector and receive vector (because if one appears in only one, it means that we either sent it a point, or received one, but has no used of it at all (otherwise it would require a symetric call))
        // this->EnsureSymmetry(sentProc, {alreadyRecvProcs, smallRangeAgent.getRecvProc(), bigRangeAgent.getRecvProc()});    // todo: uncomment
        this->EnsureSymmetry(sentProc, {alreadyRecvProcs, smallRangeAgent.getRecvProc(), bigRangeAgent.getRecvProc()});    
    }
    #endif // RICH_MPI
}

#ifdef RICH_MPI
    vector<vector<std::size_t>> &Voronoi3D::GetGhostIndeces(void)
    {
        return Nghost_;
    }
#endif // RICH_MPI

void Voronoi3D::CalcAllCM(void)
{
    std::array<Vector3D, 4> tetra;
    size_t Nfaces = FaceNeighbors_.size();
    assert(Nfaces == 0 or Nfaces >= 4);
    Vector3D vtemp;
    std::vector<Vector3D> vectemp;
    double vol;
    for (size_t i = 0; i < Nfaces; ++i)
    {
        size_t N0 = FaceNeighbors_[i].first;
        size_t N1 = FaceNeighbors_[i].second;
        size_t Npoints = PointsInFace_[i].size();
        vectemp.resize(Npoints);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (size_t j = 0; j < Npoints; ++j)
            vectemp[j] = tetra_centers_[PointsInFace_[i][j]];
        Npoints -= 2;
        tetra[0] = vectemp[0];

        for (std::size_t j = 0; j < Npoints; ++j)
        {
            tetra[1] = vectemp[j + 1];
            tetra[2] = vectemp[j + 2];
            if (N1 < Norg_)
            {
                tetra[3] = del_.points_[N1];
                vol = std::abs(GetTetraVolume(tetra));
                GetTetraCM(tetra, vtemp);
                volume_[N1] += vol;
                vtemp *= vol;
                CM_[N1] += vtemp;
            }
            tetra[3] = del_.points_[N0];
            vol = std::abs(GetTetraVolume(tetra));
            GetTetraCM(tetra, vtemp);
            volume_[N0] += vol;
            vtemp *= vol;
            CM_[N0] += vtemp;
        }
    }
#ifdef __INTEL_COMPILER
    //#pragma vector aligned
#pragma omp simd
#endif
    for (size_t i = 0; i < Norg_; ++i)
        CM_[i] *= (1.0 / volume_[i]);
    // Recalc points with high aspect ratio
    for (size_t i = 0; i < Norg_; ++i)
    {
        if (fastabs(CM_[i] - del_.points_[i]) > 0.4 * GetWidth(i))
        {
            tetra[3] = CM_[i];
            CM_[i] = Vector3D();
            volume_[i] = 0;
            Nfaces = FacesInCell_[i].size();
            for (size_t k = 0; k < Nfaces; ++k)
            {
                size_t Face = FacesInCell_[i][k];
                size_t Npoints = PointsInFace_[Face].size();
                tetra[0] = tetra_centers_[PointsInFace_[Face][0]];
                for (std::size_t j = 0; j < Npoints - 2; ++j)
                {
                    tetra[1] = tetra_centers_[PointsInFace_[Face][j + 1]];
                    tetra[2] = tetra_centers_[PointsInFace_[Face][j + 2]];
                    double vol2 = std::abs(GetTetraVolume(tetra));
                    volume_[i] += vol2;
                    GetTetraCM(tetra, vtemp);
                    CM_[i] += vol2 * vtemp;
                }
            }
            CM_[i] *= (1.0 / volume_[i]);
        }
    }
}

std::pair<Vector3D, Vector3D> Voronoi3D::GetBoxCoordinates(void) const
{
    return std::pair<Vector3D, Vector3D>(ll_, ur_);
}

void Voronoi3D::BuildNoBox(vector<Vector3D> const &points, vector<vector<Vector3D>> const &ghosts, vector<size_t> toduplicate)
{
    assert(points.size() > 0);
    // Clear data
    PointTetras_.clear();
    PointTetras_.shrink_to_fit();
    R_.clear();
    R_.reserve(points.size());
    tetra_centers_.clear();
    tetra_centers_.reserve(points.size() * 7);
    del_.Clean();
    // Voronoi Data
    FacesInCell_.clear();
    FacesInCell_.shrink_to_fit();
    PointsInFace_.clear();
    PointsInFace_.shrink_to_fit();
    FaceNeighbors_.clear();
    FaceNeighbors_.shrink_to_fit();
    CM_.clear();
    CM_.shrink_to_fit();
    Face_CM_.clear();
    Face_CM_.shrink_to_fit();
    volume_.clear();
    volume_.shrink_to_fit();
    area_.clear();
    area_.shrink_to_fit();
    Norg_ = points.size();
    #ifdef RICH_MPI
        duplicatedprocs_.clear();
        duplicatedprocs_.shrink_to_fit();
        duplicated_points_.clear();
        duplicated_points_.shrink_to_fit();
        Nghost_.clear();
        Nghost_.shrink_to_fit();
    #endif // RICH_MPI

    std::vector<size_t> order = HilbertOrder3D(points);

    del_.Build(points, ur_, ll_, order);

    for (size_t i = 0; i < ghosts.size(); ++i)
    {
        del_.BuildExtra(ghosts[i]);
    }
    vector<std::pair<size_t, size_t>> duplicate(6);
    for (size_t j = 0; j < toduplicate.size(); ++j)
    {
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
        for (size_t i = 0; i < 6; ++i)
            duplicate[i] = std::pair<size_t, size_t>(i, toduplicate[j]);
        vector<vector<size_t>> past_duplicates;
        vector<Vector3D> extra_points = CreateBoundaryPoints(duplicate, past_duplicates);
        del_.BuildExtra(extra_points);
    }

    R_.resize(del_.tetras_.size());
    R_.shrink_to_fit();
    // std::fill(R_.begin(), R_.end(), RADIUS_UNINITIALIZED);
    tetra_centers_.resize(R_.size());
    tetra_centers_.shrink_to_fit();
    bigtet_ = SetPointTetras();

    CM_.resize(Norg_);
    CM_.shrink_to_fit();
    volume_.resize(Norg_, 0);
    volume_.shrink_to_fit();
    // Create Voronoi
    BuildVoronoi(order);

    CalcAllCM();
    CM_.resize(del_.points_.size());
    CM_.shrink_to_fit();
    for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
        if (BoundaryFace(i))
            CalcRigidCM(i);
}

void Voronoi3D::BuildDebug(int rank)
{
    std::vector<size_t> order = read_vecst("order_" + std::to_string(rank) + ".bin");
    std::vector<Vector3D> points = read_vec3d("points0_" + std::to_string(rank) + ".bin");
    Norg_ = points.size();
    std::vector<Vector3D> bb = read_vec3d("bb_" + std::to_string(rank) + ".bin");
    del_.Build(points, bb[1], bb[0], order);
    points = read_vec3d("points1_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);
    points = read_vec3d("points2_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);
    points = read_vec3d("points3_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);
    points = read_vec3d("points4_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);

    bigtet_ = SetPointTetras();

    R_.resize(del_.tetras_.size());
    // std::fill(R_.begin(), R_.end(), RADIUS_UNINITIALIZED);
    tetra_centers_.resize(R_.size());

    CM_.resize(del_.points_.size());
    volume_.resize(Norg_, 0);
    // Create Voronoi
    BuildVoronoi(order);

    std::vector<double>().swap(R_);
    std::vector<tetra_vec>().swap(PointTetras_);
    std::vector<Tetrahedron>().swap(del_.tetras_);

    CalcAllCM();
    for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
        if (BoundaryFace(i))
            CalcRigidCM(i);
}

void Voronoi3D::BuildPartially(const std::vector<Vector3D> &allPoints, const std::vector<size_t> &indicesToBuild)
{
    if(this->radiuses.size() < allPoints.size())
    {
        this->radiuses.resize(allPoints.size(), RADIUS_UNINITIALIZED);
    }
    
    if(this->all_CM.size() < allPoints.size())
    {
        this->all_CM.resize(allPoints.size());
    }
    
    this->allMyPoints = allPoints;
    // this->radiuses.resize(allPoints.size(), RADIUS_UNINITIALIZED);
    std::vector<Vector3D> activePoints = VectorValues(allPoints, indicesToBuild);
    
    size_t pointsCounter = 0;
    this->indicesInAllMyPoints = Tessellation3D::AllPointsMap();
    for(const size_t &pointIdx : indicesToBuild)
    {
        this->indicesInAllMyPoints[pointIdx] = pointsCounter;
        pointsCounter++;
    }

    this->BuildInitialize(activePoints.size());
    // Norg_ = points.size();
    // // Clear data
    // PointTetras_.clear();
    // R_.clear();
    // R_.reserve(this->Norg_ * 11);
    // tetra_centers_.clear();
    // tetra_centers_.reserve(this->Norg_ * 11);
    // // Voronoi Data
    // FacesInCell_.clear();
    // PointsInFace_.clear();
    // FaceNeighbors_.clear();
    // CM_.clear();
    // Face_CM_.clear();
    // volume_.clear();
    // area_.clear();
    // Nghost_.clear();

    std::chrono::high_resolution_clock::time_point start, end;
    start = std::chrono::high_resolution_clock::now();
    
    std::vector<size_t> order = HilbertOrder3D(activePoints);

    // build delaunay
    if(not activePoints.empty())
    {
        std::pair<Vector3D, Vector3D> bounding_box = std::make_pair(activePoints[0], activePoints[0]);
        for(const Vector3D &point : activePoints)
        {
            bounding_box.first.x = std::min(bounding_box.first.x, point.x);
            bounding_box.second.x = std::max(bounding_box.second.x, point.x);
            bounding_box.first.y = std::min(bounding_box.first.y, point.y);
            bounding_box.second.y = std::max(bounding_box.second.y, point.y);
            bounding_box.first.z = std::min(bounding_box.first.z, point.z);
            bounding_box.second.z = std::max(bounding_box.second.z, point.z);
        }

        bounding_box.first.x = std::min(bounding_box.first.x, this->ll_.x);
        bounding_box.first.y = std::min(bounding_box.first.y, this->ll_.y);
        bounding_box.first.z = std::min(bounding_box.first.z, this->ll_.z);
        bounding_box.second.x = std::max(bounding_box.second.x, this->ur_.x);
        bounding_box.second.y = std::max(bounding_box.second.y, this->ur_.y);
        bounding_box.second.z = std::max(bounding_box.second.z, this->ur_.z);
        
        // performs internal tesselation:
        // std::cout << "checking duplications..." << std::endl;
        // reportDuplications(new_points);
        // order = HilbertOrder3D(activePoints);
        
        // initial build for the points
        this->del_.Build(activePoints, bounding_box.second, bounding_box.first, order);    
    }        
    // updates the radiuses array of the tetrahedra, as well as the lists for each point what tetras it belongs to
    this->R_.resize(this->del_.tetras_.size());
    this->R_.shrink_to_fit();
    // std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
    this->tetra_centers_.resize(this->R_.size());
    this->tetra_centers_.shrink_to_fit();
    this->bigtet_ = SetPointTetras();
    
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time for initial build: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    this->allMyPointsTree = std::make_shared<OctTree<IndexedVector3D>>(IndexedVector3D(this->ll_, std::numeric_limits<size_t>::max()),
                                                                    IndexedVector3D(this->ur_, std::numeric_limits<size_t>::max()));
    size_t allPointsNum = this->allMyPoints.size();
    for(size_t pointIdx = 0; pointIdx < allPointsNum; pointIdx++)
    {
        const Vector3D &point = activePoints[pointIdx];
        this->allMyPointsTree->insert(IndexedVector3D(point.x, point.y, point.z, pointIdx));
    }

    if(this->allMyPoints.size() == activePoints.size())
    {
        // not a real parital build
        this->myPointsTree = this->allMyPointsTree;
    }
    else
    {
        this->myPointsTree = std::make_shared<OctTree<IndexedVector3D>>(IndexedVector3D(this->ll_, std::numeric_limits<size_t>::max()),
                                                                        IndexedVector3D(this->ur_, std::numeric_limits<size_t>::max()));
        for(size_t pointIdx = 0; pointIdx < this->Norg_; pointIdx++)
        {
            const Vector3D &point = activePoints[pointIdx];
            this->myPointsTree->insert(IndexedVector3D(point.x, point.y, point.z, pointIdx));
        }
    }

    this->UpdateRadiuses(activePoints);

    this->UpdateRangeFinder();
    
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Time for data structures initialization: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;

    #ifdef RICH_MPI
        this->BringGhostPointsToBuild(MPI_COMM_SELF);
    #else // RICH_MPI
        this->BringGhostPointsToBuild(); 
    #endif // RICH_MPI

    // vector<std::pair<std::size_t, std::size_t>> ghost_index = SerialFirstIntersections();
    // vector<vector<size_t>> past_duplicates;
    // vector<Vector3D> extra_points = CreateBoundaryPoints(ghost_index, past_duplicates);

    // del_.BuildExtra(extra_points);

    // R_.resize(del_.tetras_.size());
    // std::fill(R_.begin(), R_.end(), RADIUS_UNINITIALIZED);
    // tetra_centers_.resize(R_.size());
    // bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);
    // ghost_index = SerialFindIntersections(true);
    // extra_points = CreateBoundaryPoints(ghost_index, past_duplicates);
    // del_.BuildExtra(extra_points);

    // R_.resize(del_.tetras_.size());
    // std::fill(R_.begin(), R_.end(), RADIUS_UNINITIALIZED);
    // tetra_centers_.resize(R_.size());
    // bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);
    // ghost_index = SerialFindIntersections(false);
    // extra_points = CreateBoundaryPoints(ghost_index, past_duplicates);
    // del_.BuildExtra(extra_points);
    // bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

    // std::vector<std::pair<size_t, size_t>>().swap(ghost_index);
    // std::vector<std::vector<size_t>>().swap(past_duplicates);
    // std::vector<Vector3D>().swap(extra_points);

    CM_.resize(del_.points_.size());
    CM_.shrink_to_fit();
    volume_.resize(Norg_);
    volume_.shrink_to_fit();

    // Create Voronoi
    BuildVoronoi(order);

    this->UpdateCMs();
}

void Voronoi3D::BuildVoronoi(std::vector<size_t> const &order)
{
    FacesInCell_.resize(Norg_);
    area_.resize(Norg_ * 10);
    Face_CM_.resize(Norg_ * 10);
    FaceNeighbors_.resize(Norg_ * 10);
    PointsInFace_.resize(Norg_ * 10);

    std::array<size_t, 128> temp, temp3;
    // Build all voronoi points
    std::size_t Ntetra = del_.tetras_.size();
    for (size_t i = 0; i < Ntetra; ++i)
        if (ShouldCalcTetraRadius(del_.tetras_[i], Norg_))
            CalcTetraRadiusCenter(i);
    // Organize the faces and assign them to cells
    std::array<double, 128> diffs, Atempvec;

    size_t FaceCounter = 0;
    boost::container::flat_set<size_t> neigh_set;
    point_vec *temp_points_in_face;
    std::array<Vector3D, 128> clean_vec;
    std::array<double, 128> area_vec_temp;

    //std::vector<Vector3D, boost::alignment::aligned_allocator<Vector3D, 32> > clean_vec;
    for (size_t i = 0; i < Norg_; ++i)
    {
        neigh_set.clear();
        neigh_set.reserve(20);
        size_t point = order[i];
        size_t ntet = PointTetras_[point].size();
        // for each point loop over its tetras
        for (size_t j = 0; j < ntet; ++j)
        {
            const size_t tetcheck = PointTetras_[point][j];
            for (size_t k = 0; k < 4; ++k)
            {
                size_t point_other = del_.tetras_[tetcheck].points[k];
                if (point_other != point && point_other > point)
                {
                    // Did we already build this face?
                    if (neigh_set.find(point_other) == neigh_set.end())
                    {
                        size_t temp_size = 0;
                        // Find all tetras for face
                        temp[0] = tetcheck;
                        ++temp_size;
                        size_t next_check = NextLoopTetra(del_.tetras_[tetcheck], tetcheck, point, point_other);
                        size_t cur_check = next_check;
                        size_t last_check = tetcheck;
                        while (next_check != tetcheck)
                        {
                            Tetrahedron const &tet_check = del_.tetras_[cur_check];
                            temp[temp_size] = cur_check;
                            ++temp_size;
                            next_check = NextLoopTetra(tet_check, last_check, point, point_other);
                            last_check = cur_check;
                            cur_check = next_check;
                        }
                        // Is face too small?
                        if (temp_size < 3)
                            continue;
                        temp_points_in_face = &PointsInFace_[FaceCounter];
                        //temp_points_in_face->reserve(8);
                        double Asize = CleanDuplicates(temp, tetra_centers_, *temp_points_in_face, ScalarProd(del_.points_[point] - del_.points_[point_other], del_.points_[point] - del_.points_[point_other]), diffs, clean_vec, temp_size);
                        if (temp_points_in_face->size() < 3)
                            continue;
                        CalcFaceAreaCM(*temp_points_in_face, tetra_centers_, clean_vec, area_[FaceCounter],
                                                     Face_CM_[FaceCounter], Atempvec);
                        if (area_[FaceCounter] < (Asize * (IsPointOutsideBox(point_other) ? 1e-14 : 1e-15)))
                            continue;
                        if (point_other >= Norg_ && point_other < (Norg_ + 4))
                        {
                            UniversalError eo("Neighboring big tet point");
                            eo.addEntry("point", point);
                            eo.addEntry("point_other", point_other);
                            eo.addEntry("Norg", Norg_);
                            eo.addEntry("Bounding Box", std::make_pair(this->ll_, this->ur_));
                            eo.addEntry("bigtet v0", this->del_.points_[Norg_]);
                            eo.addEntry("bigtet v1", this->del_.points_[Norg_ + 1]);
                            eo.addEntry("bigtet v2", this->del_.points_[Norg_ + 2]);
                            eo.addEntry("bigtet v3", this->del_.points_[Norg_ + 3]);
                            throw eo;
                        }
                        // Make faces right handed
                        MakeRightHandFace(*temp_points_in_face, del_.points_[point], tetra_centers_, temp3, area_[FaceCounter]);
                        try
                        {
                            CleanSameLine(*temp_points_in_face, tetra_centers_, area_vec_temp);
                        }
                        catch(UniversalError &eo)
                        {
                            eo.addEntry("Points0", del_.points_[point]);
                            eo.addEntry("Points1", del_.points_[point_other]);
                            eo.addEntry("diff", abs(del_.points_[point] - del_.points_[point_other]));
                            eo.addEntry("Asize", Asize);
                            int rank = 0;
                            #ifdef RICH_MPI
                               MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                            #endif // RICH_MPI
                            // throw eo;
                            continue;
                        }
                        FaceNeighbors_[FaceCounter].first = point;
                        FaceNeighbors_[FaceCounter].second = point_other;

                        FacesInCell_[point].push_back(FaceCounter);
                        if (point_other < Norg_)
                        {
                            FacesInCell_[point_other].push_back(FaceCounter);
                        }
                        neigh_set.insert(point_other);
                        ++FaceCounter;
                        // realloc memory if needed
                        if (FaceCounter == FaceNeighbors_.size())
                        {
                            area_.resize(static_cast<size_t>(static_cast<double>(area_.size()) * 1.25));
                            Face_CM_.resize(static_cast<size_t>(static_cast<double>(Face_CM_.size()) * 1.25));
                            FaceNeighbors_.resize(static_cast<size_t>(static_cast<double>(FaceNeighbors_.size()) * 1.25));
                            PointsInFace_.resize(static_cast<size_t>(static_cast<double>(PointsInFace_.size()) * 1.25));
                        }
                    }
                }
            }
        }
    }

    // Fix Face CM (this prevents large face velocities for close by points)
    size_t Nfaces = FaceCounter;
    Vector3D mid, norm;
    for (size_t i = 0; i < Nfaces; ++i)
    {
        mid = del_.points_[FaceNeighbors_[i].first];
        mid += del_.points_[FaceNeighbors_[i].second];
        mid *= 0.5;
        norm = del_.points_[FaceNeighbors_[i].second];
        norm -= del_.points_[FaceNeighbors_[i].first];
        Face_CM_[i] -= ScalarProd(Face_CM_[i] - mid, norm) * norm / ScalarProd(norm, norm);
    }

    area_.resize(FaceCounter);
    area_.shrink_to_fit();
    Face_CM_.resize(FaceCounter);
    Face_CM_.shrink_to_fit();
    FaceNeighbors_.resize(FaceCounter);
    FaceNeighbors_.shrink_to_fit();
    PointsInFace_.resize(FaceCounter);
    PointsInFace_.shrink_to_fit();
    for(size_t i = 0; i < Norg_; ++i)
    {
        FacesInCell_[i].shrink_to_fit();
    }
}

inline double Voronoi3D::GetRadius(const size_t &index) const
{ 
    R_[index] = (R_[index] < 0)? CalcTetraRadiusCenter(index) : R_[index];
    if(std::isnan(this->R_[index]) or not std::isfinite(this->R_[index]))
    {
        UniversalError eo("Voronoi3D:GetRadius: Radius is invalid");
        size_t N_points = this->del_.points_.size();
		bool found = false;
        for(size_t i = 0; i < N_points; ++i)
		{
			for(size_t j = 0; j < N_points; ++j)
			{
				if(i != j and this->del_.points_[i] == this->del_.points_[j])
				{
					eo.Append2ErrorMessage(" - Duplicated point found");
					eo.addEntry("Point1", this->del_.points_[i]);
					eo.addEntry("Point2", this->del_.points_[j]);
                    eo.addEntry("Point Index 1", i);
                    eo.addEntry("Point Index 2", j);
                    found = true;
                    break;
				}
			}
            if(found)
            {
                break;
            }
		}
		eo.Append2ErrorMessage(" - Though no duplicated points found");
        eo.addEntry("Radius", this->R_[index]);
        eo.addEntry("Tetra Index", index);
        const Tetrahedron &tet = this->del_.tetras_[index];
        eo.addEntry("Tetra Points Indices", std::vector({tet.points[0], tet.points[1], tet.points[2], tet.points[3]}));
        eo.addEntry("Tetra Points", std::vector({this->del_.points_[tet.points[0]], this->del_.points_[tet.points[1]], this->del_.points_[tet.points[2]], this->del_.points_[tet.points[3]]}));
        eo.addEntry("Norg", this->Norg_);
        throw eo;
    }
    return this->R_[index];
}

void Voronoi3D::FindIntersectionsSingle(vector<Face> const &box, std::size_t point, Sphere<Vector3D> &sphere,
                                                                                vector<size_t> &intersecting_faces, std::vector<double> &Rtemp, std::vector<Vector3D> &vtemp)
{
    intersecting_faces.clear();
    std::size_t N = PointTetras_[point].size();
    Rtemp.resize(N);
    vtemp.resize(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        Rtemp[i] = GetRadius(PointTetras_[point][i]);
        vtemp[i] = tetra_centers_[PointTetras_[point][i]];
    }
    size_t bsize = box.size();
    for (std::size_t j = 0; j < bsize; ++j)
    {
        Vector3D normal = CrossProduct(box[j].vertices[1] - box[j].vertices[0], box[j].vertices[2] - box[j].vertices[0]);
        normal *= (1.0 / fastsqrt(ScalarProd(normal, normal)));
        for (std::size_t i = 0; i < N; ++i)
        {
            sphere.radius = Rtemp[i];
            sphere.center = vtemp[i];
            if (FaceSphereIntersections(box[j], sphere, normal))
            {
                intersecting_faces.push_back(j);
                break;
            }
        }
    }
}

void Voronoi3D::GetPointToCheck(std::size_t point, vector<unsigned char> const &checked, vector<std::size_t> &res)
{
    res.clear();
    std::size_t ntetra = PointTetras_[point].size();
    for (std::size_t i = 0; i < ntetra; ++i)
    {
        size_t tetra = PointTetras_[point][i];
        for (std::size_t j = 0; j < 4; ++j)
            if (del_.tetras_[tetra].points[j] < Norg_ && checked[del_.tetras_[tetra].points[j]] == 0)
                res.push_back(del_.tetras_[tetra].points[j]);
    }
    std::sort(res.begin(), res.end());
    res = unique(res);
}

std::size_t Voronoi3D::GetFirstPointToCheck(void) const
{
    std::size_t i;
    Tetrahedron const &tet = del_.tetras_[bigtet_];
    for (i = 0; i < 4; ++i)
        if (tet.points[i] < Norg_)
            break;
    if (i < 4)
        return tet.points[i];
    else
        throw UniversalError("Can't find first point to start boundary search");
}

vector<std::pair<std::size_t, std::size_t>> Voronoi3D::SerialFirstIntersections(void)
{
    vector<Face> box;
    vector<Vector3D> normals;
    this->InitialBoxBuild(box, normals);
    size_t Nfaces = box.size();

    //    vector<std::size_t> point_neigh;
    vector<std::pair<std::size_t, std::size_t>> res;
    Sphere<Vector3D> sphere;
    vector<unsigned char> will_check(Norg_, 0);
    std::size_t cur_loc;
    std::stack<std::size_t> check_stack;
    FirstCheckList(check_stack, will_check, Norg_, del_, PointTetras_);
    std::vector<double> vdist(Nfaces);
    std::vector<Vector3D> vtemp(Nfaces);
    while (!check_stack.empty())
    {
        cur_loc = check_stack.top();
        check_stack.pop();
        double inv_max = 0;
        size_t max_loc = 0;
        size_t j = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
        for (; j < Nfaces; ++j)
        {
            vtemp[j] = del_.points_[cur_loc];
            vtemp[j] -= box[j].vertices[0];
            double sprod = vtemp[j].x * normals[j].x + vtemp[j].y * normals[j].y + vtemp[j].z * normals[j].z;
            vdist[j] = 1.0 / std::abs(sprod);
        }
        j = 0;
        for (; j < Nfaces; ++j)
        {
            if (vdist[j] > inv_max)
            {
                inv_max = vdist[j];
                max_loc = j;
            }
        }
        res.push_back(std::pair<std::size_t, std::size_t>(max_loc, cur_loc));
    }
    return res;
}

vector<std::pair<std::size_t, std::size_t>> Voronoi3D::SerialFindIntersections(bool first_run)
{
    if (Norg_ < 50)
    {
        vector<std::pair<std::size_t, std::size_t>> res;
        size_t const Nfaces = box_faces_.empty() ? 6 : box_faces_.size();
        res.reserve(Norg_ * Nfaces);
        for (size_t i = 0; i < Norg_; ++i)
            for (size_t j = 0; j < Nfaces; ++j)
                res.push_back(std::pair<std::size_t, std::size_t>(j, i));
        return res;
    }
    std::stack<std::size_t> check_stack;
    vector<Face> box = box_faces_.empty() ? BuildBox(ll_, ur_) : box_faces_;
    vector<std::size_t> point_neigh;
    vector<std::pair<std::size_t, std::size_t>> res;
    Sphere<Vector3D> sphere;
    vector<unsigned char> checked(Norg_, 0), will_check(Norg_, 0);
    std::size_t cur_loc;
    if (first_run)
    {
        FirstCheckList(check_stack, will_check, Norg_, del_, PointTetras_);
        //cur_loc = check_stack.top();
        check_stack.pop();
    }
    else
    {
        cur_loc = GetFirstPointToCheck();
        check_stack.push(cur_loc);
        will_check[cur_loc] = true;
    }
    vector<size_t> intersecting_faces;
    std::vector<double> Rtemp;
    std::vector<Vector3D> vtemp;
    while (!check_stack.empty())
    {
        cur_loc = check_stack.top();
        check_stack.pop();
        checked[cur_loc] = true;
        // Does sphere have any intersections?
        bool added = false;
        FindIntersectionsSingle(box, cur_loc, sphere, intersecting_faces, Rtemp, vtemp);
        if (!intersecting_faces.empty())
        {
            added = true;
            for (std::size_t j = 0; j < intersecting_faces.size(); ++j)
                res.push_back(std::pair<std::size_t, std::size_t>(intersecting_faces[j], cur_loc));
        }
        if (added && !first_run)
        {
            GetPointToCheck(cur_loc, checked, point_neigh);
            std::size_t Nneigh = point_neigh.size();
            for (std::size_t j = 0; j < Nneigh; ++j)
                if (point_neigh[j] < Norg_ && !will_check[point_neigh[j]])
                {
                    check_stack.push(point_neigh[j]);
                    will_check[point_neigh[j]] = true;
                }
        }
    }
    return res;
}

double Voronoi3D::CalcTetraRadiusCenter(const size_t &index) const
{
    Vector3D v2(del_.points_[del_.tetras_[index].points[1]]);
    v2 -= del_.points_[del_.tetras_[index].points[0]];
    Vector3D v3(del_.points_[del_.tetras_[index].points[2]]);
    v3 -= del_.points_[del_.tetras_[index].points[0]];
    Vector3D v4(del_.points_[del_.tetras_[index].points[3]]);
    v4 -= del_.points_[del_.tetras_[index].points[0]];

    Mat33<double> m_a(v2.x, v2.y, v2.z,
                                        v3.x, v3.y, v3.z,
                                        v4.x, v4.y, v4.z);
    double a = m_a.determinant();
    if(std::abs(a) < 100 * std::numeric_limits<double>::min())
        return CalcTetraRadiusCenterHiPrecision(index);
    Mat33<double> m_Dx(ScalarProd(v2, v2), v2.y, v2.z,
                                         ScalarProd(v3, v3), v3.y, v3.z,
                                         ScalarProd(v4, v4), v4.y, v4.z);
    double DDx = m_Dx.determinant();

    Mat33<double> m_Dy(ScalarProd(v2, v2), v2.x, v2.z,
                                         ScalarProd(v3, v3), v3.x, v3.z,
                                         ScalarProd(v4, v4), v4.x, v4.z);
    double DDy = -m_Dy.determinant();

    Mat33<double> m_Dz(ScalarProd(v2, v2), v2.x, v2.y,
                                         ScalarProd(v3, v3), v3.x, v3.y,
                                         ScalarProd(v4, v4), v4.x, v4.y);
    double DDz = m_Dz.determinant();
    Vector3D center = Vector3D(DDx / (2 * a), DDy / (2 * a), DDz / (2 * a)) + del_.points_[del_.tetras_[index].points[0]];
    tetra_centers_[index] = center;
    double Rres = 0.5 * std::sqrt(DDx * DDx + DDy * DDy + DDz * DDz) / std::abs(a);
    // Sanity check
    /*double Rcheck0 = fastabs(del_.points_[del_.tetras_[index].points[0]] - center);
        double Rcheck1 = fastabs(del_.points_[del_.tetras_[index].points[1]] - center);
        double Rcheck2 = fastabs(del_.points_[del_.tetras_[index].points[2]] - center);
        double Rcheck3 = fastabs(del_.points_[del_.tetras_[index].points[3]] - center);*/
    Vector3D v1(del_.points_[del_.tetras_[index].points[0]]);
    double Rcheck0 = fastabs(v1 - center);
    v2 += v1;
    v2 -= center;
    double Rcheck1 = fastabs(v2);
    v3 += v1;
    v3 -= center;
    double Rcheck2 = fastabs(v3);
    v4 += v1;
    v4 -= center;
    double Rcheck3 = fastabs(v4);
    double tol = 1 + 1e-6;
    if (((Rcheck0 + Rcheck1 + Rcheck2 + Rcheck3) * tol < (4 * Rcheck0)) || ((Rcheck0 + Rcheck1 + Rcheck2 + Rcheck3) > (tol * 4 * Rcheck0)))
        return CalcTetraRadiusCenterHiPrecision(index);
    if (Rcheck0 > tol * Rres || Rcheck0 * tol < Rres)
        return CalcTetraRadiusCenterHiPrecision(index);
    double const a_tol = 1e-8;
    if (std::abs(a) < Rres * Rres * Rres * a_tol)
        return CalcTetraRadiusCenterHiPrecision(index);
    return Rres;
}

double Voronoi3D::CalcTetraRadiusCenterHiPrecision(const size_t &index) const
{
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V0;
    V0[0] = del_.points_[del_.tetras_[index].points[0]].x;
    V0[1] = del_.points_[del_.tetras_[index].points[0]].y;
    V0[2] = del_.points_[del_.tetras_[index].points[0]].z;
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V2;
    V2[0] = del_.points_[del_.tetras_[index].points[1]].x;
    V2[1] = del_.points_[del_.tetras_[index].points[1]].y;
    V2[2] = del_.points_[del_.tetras_[index].points[1]].z;
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V3;
    V3[0] = del_.points_[del_.tetras_[index].points[2]].x;
    V3[1] = del_.points_[del_.tetras_[index].points[2]].y;
    V3[2] = del_.points_[del_.tetras_[index].points[2]].z;
    std::array<boost::multiprecision::cpp_dec_float_50, 3> V4;
    V4[0] = del_.points_[del_.tetras_[index].points[3]].x;
    V4[1] = del_.points_[del_.tetras_[index].points[3]].y;
    V4[2] = del_.points_[del_.tetras_[index].points[3]].z;
    V2[0] -= V0[0];
    V2[1] -= V0[1];
    V2[2] -= V0[2];
    V3[0] -= V0[0];
    V3[1] -= V0[1];
    V3[2] -= V0[2];
    V4[0] -= V0[0];
    V4[1] -= V0[1];
    V4[2] -= V0[2];
    std::array<boost::multiprecision::cpp_dec_float_50, 9> mat;
    mat[0] = V2[0];
    mat[1] = V2[1];
    mat[2] = V2[2];
    mat[3] = V3[0];
    mat[4] = V3[1];
    mat[5] = V3[2];
    mat[6] = V4[0];
    mat[7] = V4[1];
    mat[8] = V4[2];
    boost::multiprecision::cpp_dec_float_50 ba = Calc33Det(mat);
    mat[0] = V2[0] * V2[0] + V2[1] * V2[1] + V2[2] * V2[2];
    mat[1] = V2[1];
    mat[2] = V2[2];
    mat[3] = V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2];
    mat[4] = V3[1];
    mat[5] = V3[2];
    mat[6] = V4[0] * V4[0] + V4[1] * V4[1] + V4[2] * V4[2];
    mat[7] = V4[1];
    mat[8] = V4[2];
    boost::multiprecision::cpp_dec_float_50 bDx = Calc33Det(mat);
    mat[0] = V2[0] * V2[0] + V2[1] * V2[1] + V2[2] * V2[2];
    mat[1] = V2[0];
    mat[2] = V2[2];
    mat[3] = V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2];
    mat[4] = V3[0];
    mat[5] = V3[2];
    mat[6] = V4[0] * V4[0] + V4[1] * V4[1] + V4[2] * V4[2];
    mat[7] = V4[0];
    mat[8] = V4[2];
    boost::multiprecision::cpp_dec_float_50 bDy = -Calc33Det(mat);
    mat[0] = V2[0] * V2[0] + V2[1] * V2[1] + V2[2] * V2[2];
    mat[1] = V2[0];
    mat[2] = V2[1];
    mat[3] = V3[0] * V3[0] + V3[1] * V3[1] + V3[2] * V3[2];
    mat[4] = V3[0];
    mat[5] = V3[1];
    mat[6] = V4[0] * V4[0] + V4[1] * V4[1] + V4[2] * V4[2];
    mat[7] = V4[0];
    mat[8] = V4[1];
    boost::multiprecision::cpp_dec_float_50 bDz = Calc33Det(mat);
    boost::multiprecision::cpp_dec_float_50 temp = (bDx / (2 * ba) + V0[0]);
    tetra_centers_[index].x = temp.convert_to<double>();
    temp = (bDy / (2 * ba) + V0[1]);
    tetra_centers_[index].y = temp.convert_to<double>();
    temp = (bDz / (2 * ba) + V0[2]);
    tetra_centers_[index].z = temp.convert_to<double>();
    temp = (boost::multiprecision::sqrt(bDx * bDx + bDy * bDy + bDz * bDz) / ba);
    return 0.5 * temp.convert_to<double>();
}

void Voronoi3D::GetTetraCM(std::array<Vector3D, 4> const &points, Vector3D &CM) const
{
    double x = 0, y = 0, z = 0;
    //CM.Set(0, 0, 0);
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+ \
                                                     : x, y, z)
#endif
    for (std::size_t i = 0; i < 4; i++)
    {
        x += points[i].x;
        y += points[i].y;
        z += points[i].z;
    }
    CM.Set(x, y, z);
    CM *= 0.25;
}

double Voronoi3D::GetTetraVolume(std::array<Vector3D, 4> const &points) const
{
    return std::abs(orient3d(points)) / 6.0;
}

/*
void Voronoi3D::CalcCellCMVolume(std::size_t index)
{
    volume_[index] = 0;
    CM_[index] = Vector3D();
    std::size_t Nfaces = FacesInCell_[index].size();
    std::array<Vector3D, 4> tetra;
    tetra[3] = del_.points_[index];
    Vector3D vtemp;
    for (std::size_t i = 0; i < Nfaces; ++i)
        {
            std::size_t face = FacesInCell_[index][i];
            std::size_t Npoints = PointsInFace_[face].size();
            tetra[0] = tetra_centers_[PointsInFace_[face][0]];
            double fvol = 0;
            for (std::size_t j = 0; j < Npoints - 2; ++j)
	{
	    tetra[1] = tetra_centers_[PointsInFace_[face][j + 1]];
	    tetra[2] = tetra_centers_[PointsInFace_[face][j + 2]];
	    double vol = GetTetraVolume(tetra);
	    fvol += std::abs(vol);
	    GetTetraCM(tetra, vtemp);
	    CM_[index] += std::abs(vol)*vtemp;
	}
            volume_[index] += fvol;
        }
    CM_[index] = CM_[index] / volume_[index];
}
*/

void Voronoi3D::output(std::string const &filename) const
{

    std::ofstream file_handle(filename.c_str(), std::ios::out | std::ios::binary);
    assert(file_handle.is_open());
    binary_write_single_int(static_cast<int>(Norg_), file_handle);

    // Points
    for (std::size_t i = 0; i < Norg_; ++i)
    {
        binary_write_single_double(del_.points_[i].x, file_handle);
        binary_write_single_double(del_.points_[i].y, file_handle);
        binary_write_single_double(del_.points_[i].z, file_handle);
    }

    binary_write_single_int(static_cast<int>(tetra_centers_.size()), file_handle);
    // Face Points
    for (std::size_t i = 0; i < tetra_centers_.size(); ++i)
    {
        binary_write_single_double(tetra_centers_[i].x, file_handle);
        binary_write_single_double(tetra_centers_[i].y, file_handle);
        binary_write_single_double(tetra_centers_[i].z, file_handle);
    }

    // Faces in cell
    for (std::size_t i = 0; i < Norg_; ++i)
    {
        binary_write_single_int(static_cast<int>(FacesInCell_[i].size()), file_handle);
        for (std::size_t j = 0; j < FacesInCell_[i].size(); ++j)
            binary_write_single_int(static_cast<int>(FacesInCell_[i][j]), file_handle);
    }

    // Points in Face
    binary_write_single_int(static_cast<int>(PointsInFace_.size()), file_handle);
    for (std::size_t i = 0; i < PointsInFace_.size(); ++i)
    {
        binary_write_single_int(static_cast<int>(PointsInFace_[i].size()), file_handle);
        for (std::size_t j = 0; j < PointsInFace_[i].size(); ++j)
            binary_write_single_int(static_cast<int>(PointsInFace_[i][j]), file_handle);
    }

    file_handle.close();
}

#ifdef RICH_MPI
void Voronoi3D::output_buildextra(std::string const &filename) const
{
    std::ofstream file_handle(filename.c_str(), std::ios::out | std::ios::binary);
    assert(file_handle.is_open());
    size_t stemp = Norg_;
    binary_write_single_int(static_cast<int>(stemp), file_handle);
    stemp = del_.points_.size();
    binary_write_single_int(static_cast<int>(stemp), file_handle);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Points
    for (std::size_t i = 0; i < stemp; ++i)
    {
        binary_write_single_double(del_.points_[i].x, file_handle);
        binary_write_single_double(del_.points_[i].y, file_handle);
        binary_write_single_double(del_.points_[i].z, file_handle);
    }

    binary_write_single_int(static_cast<int>(duplicatedprocs_.size()), file_handle);
    // Procs
    assert(duplicatedprocs_.size() == Nghost_.size());
    for (size_t i = 0; i < duplicatedprocs_.size(); ++i)
    {
        binary_write_single_int(static_cast<int>(duplicatedprocs_[i]), file_handle);
        binary_write_single_int(static_cast<int>(Nghost_[i].size()), file_handle);
        for (size_t j = 0; j < Nghost_[i].size(); ++j)
            binary_write_single_int(static_cast<int>(Nghost_[i][j]), file_handle);
    }
    file_handle.close();
}
#endif

size_t Voronoi3D::GetContainingCell(const Vector3D &point) const
{
    return this->myPointsTree->closestPoint(point).getIndex();
}

std::size_t Voronoi3D::GetPointNo(void) const
{
    return Norg_;
}

const Vector3D &Voronoi3D::GetMeshPoint(std::size_t index) const
{
    return del_.points_[index];
}

double Voronoi3D::GetArea(std::size_t index) const
{
    return area_[index];
}

Vector3D const &Voronoi3D::GetCellCM(std::size_t index) const
{
    return this->CM_[index];
}

std::size_t Voronoi3D::GetTotalFacesNumber(void) const
{
    return FaceNeighbors_.size();
}

double Voronoi3D::GetWidth(std::size_t index) const
{
    return std::pow(3 * volume_[index] * 0.25 / M_PI, 0.3333333333);
}

double Voronoi3D::GetVolume(std::size_t index) const
{
    return volume_[index];
}

face_vec const &Voronoi3D::GetCellFaces(std::size_t index) const
{
    return FacesInCell_[index];
}

vector<Vector3D> &Voronoi3D::accessMeshPoints(void)
{
    return del_.points_;
}

const vector<Vector3D> &Voronoi3D::getMeshPoints(void) const
{
    return del_.points_;
}

const Tessellation3D::AllPointsMap &Voronoi3D::GetIndicesInAllPoints(void) const
{
    return this->indicesInAllMyPoints;
}

const std::vector<Vector3D> &Voronoi3D::getAllPoints(void) const
{
    return this->allMyPoints;
}

std::vector<Vector3D> &Voronoi3D::getAllPoints(void)
{
    return this->allMyPoints;
}

vector<std::size_t> Voronoi3D::GetNeighbors(std::size_t index) const
{
    const size_t N = FacesInCell_[index].size();
    vector<size_t> res(N);
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
    for (size_t i = 0; i < N; ++i)
    {
        size_t face = FacesInCell_[index][i];
        res[i] = FaceNeighbors_[face].first == index ? FaceNeighbors_[face].second : FaceNeighbors_[face].first;
    }
    return res;
}

void Voronoi3D::GetNeighbors(size_t index, vector<size_t> &res) const
{
    std::size_t N = FacesInCell_[index].size();
    res.resize(N);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
    for (std::size_t i = 0; i < N; ++i)
    {
        std::size_t face = FacesInCell_[index][i];
        res[i] = FaceNeighbors_[face].first == index ? FaceNeighbors_[face].second : FaceNeighbors_[face].first;
    }
}

Tessellation3D *Voronoi3D::clone(void) const
{
    return new Voronoi3D(*this);
}

Voronoi3D::Voronoi3D(Voronoi3D const &other) : ll_(other.ll_), ur_(other.ur_), Norg_(other.Norg_), bigtet_(other.bigtet_),
                                                set_temp_(other.set_temp_), stack_temp_(other.stack_temp_), del_(other.del_), PointTetras_(other.PointTetras_), R_(other.R_),
                                                tetra_centers_(other.tetra_centers_), FacesInCell_(other.FacesInCell_), PointsInFace_(other.PointsInFace_),
                                                FaceNeighbors_(other.FaceNeighbors_), CM_(other.CM_), Face_CM_(other.Face_CM_), volume_(other.volume_), area_(other.area_),
                                                #ifdef RICH_MPI
                                                    sentprocs_(other.sentprocs_), sentpoints_(other.sentpoints_), duplicatedprocs_(other.duplicatedprocs_), duplicated_points_(other.duplicated_points_),
                                                    Nghost_(other.Nghost_), self_index_(other.self_index_),
                                                #endif // RICH_MPI
                                                temp_points_(std::array<Vector3D, 4>()), temp_points2_(std::array<Vector3D, 5>()), box_faces_(other.box_faces_),
                                                #ifdef RICH_MPI
                                                    pointsManager(other.pointsManager->clone()), indexingToSave(other.indexingToSave),
                                                    rangeFinder(other.rangeFinder), radiuses(other.radiuses), allMyPoints(other.allMyPoints), allPointsWeights(other.allPointsWeights),
                                                #endif // RICH_MPI
                                                indicesInAllMyPoints(other.indicesInAllMyPoints)
                                                {}

bool Voronoi3D::NearBoundary(std::size_t index) const
{
    std::size_t N = FacesInCell_[index].size();
    for (std::size_t i = 0; i < N; ++i)
    {
        if (BoundaryFace(FacesInCell_[index][i]))
            return true;
    }
    return false;
}

bool Voronoi3D::IsPointInCell(const Vector3D &point, size_t cellIndex, bool verbose) const
{
    if(cellIndex >= this->Norg_)
    {
        UniversalError eo("Voronoi3D::IsPointInCell: cell index out of range");
        eo.addEntry("point", point);
        eo.addEntry("cellIndex", cellIndex);
        eo.addEntry("Norg", this->Norg_);
        throw eo;
    }
    std::optional<UniversalError> verboseInfo;
    if(verbose)
    {
        verboseInfo.emplace("Voronoi3D::IsPointInCell verbose");
        verboseInfo->addEntry("point", point);
        verboseInfo->addEntry("cellIndex", cellIndex);
        verboseInfo->addEntry("cellCenter", this->del_.points_[cellIndex]);

        Vector3D cellLL(std::numeric_limits<double>::max());
        Vector3D cellUR(std::numeric_limits<double>::lowest());
        const std::vector<Vector3D> &vertices = this->GetFacePoints();
        for(size_t faceIdx : this->FacesInCell_[cellIndex])
        {
            for(size_t vertIdx : this->GetAllPointsInFace()[faceIdx])
            {
                const Vector3D &v = vertices[vertIdx];
                cellLL.x = std::min(cellLL.x, v.x);
                cellLL.y = std::min(cellLL.y, v.y);
                cellLL.z = std::min(cellLL.z, v.z);
                cellUR.x = std::max(cellUR.x, v.x);
                cellUR.y = std::max(cellUR.y, v.y);
                cellUR.z = std::max(cellUR.z, v.z);
            }
        }
        BoundingBox<Vector3D> cellBox(cellLL, cellUR);  
        verboseInfo->addEntry("cellBBox", cellBox);
        verboseInfo->addEntry("is point in bounding box of cell?", cellBox.contains(point)? "yes": "no");
    }

    double width = this->GetWidth(cellIndex);
    for(size_t faceIdx : this->FacesInCell_[cellIndex])
    {
        const Vector3D &p1 = this->del_.points_[this->GetFaceNeighbors(faceIdx).first];
        const Vector3D &p2 = this->del_.points_[this->GetFaceNeighbors(faceIdx).second];
        Vector3D normal = (this->GetFaceNeighbors(faceIdx).second == cellIndex)? p2 - p1 : p1 - p2;
        Vector3D p = (p1 + p2) * 0.5;
        double dot = ScalarProd(normal, point - p);
        if(verboseInfo)
        {
            size_t neighbor = (this->GetFaceNeighbors(faceIdx).first == cellIndex)? this->GetFaceNeighbors(faceIdx).second : this->GetFaceNeighbors(faceIdx).first;
            verboseInfo->addEntry("face " + std::to_string(faceIdx) + " neighbor", neighbor);
            std::string neighborOwner = "???";
            if(neighbor < this->Norg_)
            {
                neighborOwner = "me";
            }
            else
            {
                for(size_t i = 0; i < this->Nghost_.size(); i++)
                {
                    const std::vector<size_t> &ghosts = this->Nghost_[i];
                    if(std::find(ghosts.begin(), ghosts.end(), neighbor) != ghosts.end())
                    {
                        neighborOwner = "rank " + std::to_string(this->duplicatedprocs_[i]);
                        break;
                    }
                } 
            }
            verboseInfo->addEntry("face " + std::to_string(faceIdx) + ", neighbor " + std::to_string(neighbor) + " owner", neighborOwner);
            verboseInfo->addEntry("face " + std::to_string(faceIdx) + " dot", dot);
            verboseInfo->addEntry("face " + std::to_string(faceIdx) + " neighbor point", this->GetMeshPoint(neighbor));
            verboseInfo->addEntry("face " + std::to_string(faceIdx) + " neighbor distance", abs(this->GetMeshPoint(neighbor) - point));
        }
        if(not verboseInfo and dot < -1e-12 * width)
        {
            return false;
        }
    }
    if(verboseInfo)
    {
        throw *verboseInfo;
    }
    return true;
}

bool Voronoi3D::IsPointOutsideBox(size_t index) const
{
    if(box_faces_.empty())
        return !PointInDomain(ll_, ur_, del_.points_[index]);
    else
        return !PointInPoly(box_faces_, del_.points_[index]);
}

bool Voronoi3D::IsPointOutsideBox(const Vector3D &point) const
{
    if(box_faces_.empty())
        return !PointInDomain(ll_, ur_, point);
    else
        return !PointInPoly(box_faces_, point);
}

bool Voronoi3D::BoundaryFace(std::size_t index) const
{
    if (FaceNeighbors_[index].first >= Norg_ || FaceNeighbors_[index].second >= Norg_)
    {
#ifdef RICH_MPI
        if(box_faces_.empty())
        {
            if (PointInDomain(ll_, ur_, del_.points_[std::max(FaceNeighbors_[index].first, FaceNeighbors_[index].second)]))
                return false;
            else
                return true;
        }
        else
            if(PointInPoly(box_faces_, del_.points_[std::max(FaceNeighbors_[index].first, FaceNeighbors_[index].second)]))
                return false;
            else
#endif
            return true;
    }
    else
        return false;
}

#ifdef RICH_MPI
bool Voronoi3D::CheckContinuityOfZone(void) const
{
    std::vector<bool> reached(this->Norg_, false);
    reached[0] = true;
    std::stack<std::size_t> stack;
    stack.push(0);
    while(not stack.empty())
    {
        size_t cell = stack.top();
        stack.pop();
        for(const auto &neighbor : this->GetNeighbors(cell))
        {
            if(neighbor >= this->Norg_)
            {
                continue; // ghost
            }
            if(not reached[neighbor])
            {
                reached[neighbor] = true;
                stack.push(neighbor);
            }
        }
    }
    return std::all_of(reached.cbegin(), reached.cend(), [](const bool &b){return b;});
}
#endif // RICH_MPI

#ifdef RICH_MPI
    vector<vector<std::size_t>> &Voronoi3D::GetDuplicatedPoints(void)
    {
        return duplicated_points_;
    }

    vector<vector<std::size_t>> const &Voronoi3D::GetDuplicatedPoints(void) const
    {
        return duplicated_points_;
    }
#endif // RICH_MPI

std::size_t Voronoi3D::GetTotalPointNumber(void) const
{
    return del_.points_.size();
}

vector<Vector3D> &Voronoi3D::GetAllCM(void)
{
    return CM_;
}

vector<Vector3D> Voronoi3D::GetAllCM(void) const
{
    return CM_;
}

void Voronoi3D::GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point) const
{
    result.clear();
    result.reserve(70);
    vector<std::size_t> neigh = GetNeighbors(point);
    result = neigh;
    std::size_t N = neigh.size();
    std::sort(neigh.begin(), neigh.end());
    vector<std::size_t> temp;
    for (std::size_t i = 0; i < N; ++i)
    {
        if (neigh[i] < Norg_)
        {
            temp = GetNeighbors(neigh[i]);
            result.insert(result.end(), temp.begin(), temp.end());
        }
    }
    std::sort(result.begin(), result.end());
    result = unique(result);
    result = RemoveList(result, neigh);
    RemoveVal(result, point);
}

vector<boost::container::small_vector<size_t, 8>> &Voronoi3D::GetAllPointsInFace(void)
{
    return PointsInFace_;
}

vector<boost::container::small_vector<size_t, 8>> const& Voronoi3D::GetAllPointsInFace(void)const
{
    return PointsInFace_;
}

size_t &Voronoi3D::GetPointNo(void)
{
    return Norg_;
}

size_t Voronoi3D::GetAllPointsNo(void) const
{
    return this->allMyPoints.size();
}

std::vector<std::pair<size_t, size_t>> &Voronoi3D::GetAllFaceNeighbors(void)
{
    return FaceNeighbors_;
}

const std::vector<std::pair<size_t, size_t>> &Voronoi3D::GetAllFaceNeighbors(void) const
{
    return FaceNeighbors_;
}

vector<double> &Voronoi3D::GetAllVolumes(void)
{
    return volume_;
}

vector<double> Voronoi3D::GetAllVolumes(void) const
{
    return volume_;
}

Vector3D Voronoi3D::Normal(std::size_t faceindex) const
{
    return del_.points_[FaceNeighbors_[faceindex].second] - del_.points_[FaceNeighbors_[faceindex].first];
}

bool Voronoi3D::IsGhostPoint(std::size_t index) const
{
    return index >= Norg_;
}

const Vector3D &Voronoi3D::FaceCM(std::size_t index) const
{
    return Face_CM_[index];
}

Vector3D Voronoi3D::CalcFaceVelocity(std::size_t index, Vector3D const &v0, Vector3D const &v1) const
{
    std::size_t p0 = FaceNeighbors_[index].first;
    std::size_t p1 = FaceNeighbors_[index].second;
    Vector3D r0 = GetMeshPoint(p0);
    Vector3D r1 = GetMeshPoint(p1);
    Vector3D r_diff = r1 - r0;
    double abs_r_diff = ScalarProd(r_diff, r_diff);

    Vector3D f = FaceCM(index);
    r1 += r0;
    r1 *= 0.5;
    f -= r1;
    Vector3D delta_w = ScalarProd((v0 - v1), f) * r_diff / abs_r_diff;
#ifdef RICH_DEBUG
    double dw_abs = fastabs(delta_w);
#endif // RICH_DEBUG
    Vector3D w = (v0 + v1) * 0.5;
#ifdef RICH_DEBUG
    //    double w_abs = std::max(fastabs(v0),fastabs(v1));
#endif // RICH_DEBUG
    //if (dw_abs > w_abs)
    //	delta_w *= (1 + (std::atan(dw_abs / w_abs) - 0.25 * M_PI)*2) * (w_abs / dw_abs);
#ifdef RICH_DEBUG
    if (!std::isfinite(dw_abs))
    {
        r0 = GetMeshPoint(p0);
        r1 = GetMeshPoint(p1);
        f = FaceCM(index);
        UniversalError eo("Bad Face velocity");
        eo.addEntry("Face index", index);
        eo.addEntry("Neigh 0", p0);
        eo.addEntry("Neigh 1", p1);
        eo.addEntry("Neigh 0 x", r0.x);
        eo.addEntry("Neigh 0 y", r0.y);
        eo.addEntry("Neigh 0 z", r0.z);
        eo.addEntry("Neigh 0 CMx", CM_[p0].x);
        eo.addEntry("Neigh 0 CMy", CM_[p0].y);
        eo.addEntry("Neigh 0 CMz", CM_[p0].z);
        eo.addEntry("Neigh 1 x", r1.x);
        eo.addEntry("Neigh 1 y", r1.y);
        eo.addEntry("Neigh 1 z", r1.z);
        eo.addEntry("Neigh 1 CMx", CM_[p1].x);
        eo.addEntry("Neigh 1 CMy", CM_[p1].y);
        eo.addEntry("Neigh 1 CMz", CM_[p1].z);
        eo.addEntry("Face CMx", f.x);
        eo.addEntry("Face CMy", f.y);
        eo.addEntry("Face CMz", f.z);
        eo.addEntry("V0x", v0.x);
        eo.addEntry("V0y", v0.y);
        eo.addEntry("V0z", v0.z);
        eo.addEntry("V1x", v1.x);
        eo.addEntry("V1y", v1.y);
        eo.addEntry("V1z", v1.z);
        throw eo;
    }
#endif
    w += delta_w;
    return w;
}

vector<double> &Voronoi3D::GetAllArea(void)
{
    return area_;
}

vector<Vector3D> &Voronoi3D::GetAllFaceCM(void)
{
    return Face_CM_;
}

vector<face_vec> &Voronoi3D::GetAllCellFaces(void)
{
    return FacesInCell_;
}

vector<face_vec> const& Voronoi3D::GetAllCellFaces(void)const
{
    return FacesInCell_;
}

vector<Vector3D> &Voronoi3D::GetFacePoints(void)
{
    return tetra_centers_;
}

vector<Vector3D> const &Voronoi3D::GetFacePoints(void) const
{
    return tetra_centers_;
}

point_vec const &Voronoi3D::GetPointsInFace(std::size_t index) const
{
    return PointsInFace_[index];
}

const std::pair<std::size_t, std::size_t> &Voronoi3D::GetFaceNeighbors(std::size_t face_index) const
{
    return FaceNeighbors_[face_index];
}

#ifdef RICH_MPI
    vector<int> Voronoi3D::GetDuplicatedProcs(void) const
    {
        return duplicatedprocs_;
    }

    vector<int> Voronoi3D::GetSentProcs(void) const
    {
        return sentprocs_;
    }

    vector<vector<std::size_t>> const &Voronoi3D::GetSentPoints(void) const
    {
        return sentpoints_;
    }

    vector<std::size_t> const &Voronoi3D::GetSelfIndex(void) const
    {
        return self_index_;
    }

    vector<int> &Voronoi3D::GetSentProcs(void)
    {
        return sentprocs_;
    }

    vector<vector<std::size_t>> &Voronoi3D::GetSentPoints(void)
    {
        return sentpoints_;
    }

    vector<std::size_t> &Voronoi3D::GetSelfIndex(void)
    {
        return self_index_;
    }
#endif // RICH_MPI

#ifdef RICH_MPI
void Voronoi3D::SetKernel(const std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>> &indexing)
{
    MPI_Barrier(MPI_COMM_WORLD); // everyone should set the kernel
    this->indexingToSave = indexing;
    HilbertPointsManager<Vector3D, VoronoiPayload> *hilbertPointsManager = dynamic_cast<HilbertPointsManager<Vector3D, VoronoiPayload> *>(this->pointsManager.get());
    if(hilbertPointsManager == nullptr)
    {
        // points manager is not a 'HilbertPointsManager', or was not initialized yet
        return;
    }
    hilbertPointsManager->setIndexing(indexing);
}

std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>> Voronoi3D::GetKernel() const
{
    return this->indexingToSave;
}
#endif // RICH_MPI

void Voronoi3D::SetBox(const Vector3D &ll, const Vector3D &ur)
{
    this->ll_ = ll;
    this->ur_ = ur;
    this->box_faces_ = BuildBox(this->ll_, this->ur_);
    #ifdef RICH_MPI
        this->pointsManager = std::make_shared<HilbertPointsManager<Vector3D, VoronoiPayload>>(this->ll_, this->ur_, MPI_COMM_WORLD);
        // this->radiuses.clear();
    #endif // RICH_MPI
}

#ifdef RICH_MPI
void Voronoi3D::SetBox(Vector3D const &ll, Vector3D const &ur, const std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>> &newIndexing)
{
    this->SetBox(ll, ur);
    this->SetKernel(newIndexing);
}

const std::vector<double> &Voronoi3D::GetPointsBuildWeights() const
{
    return this->allPointsWeights;
}

const std::shared_ptr<EnvironmentAgent<Vector3D>> Voronoi3D::GetEnvironmentAgent() const
{
    return this->pointsManager->getEnvironmentAgent();
}

#endif // RICH_MPI
