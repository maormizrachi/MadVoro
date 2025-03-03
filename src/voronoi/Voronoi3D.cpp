#include "Voronoi3D.hpp"
#include <algorithm>
#include <cfloat>
#include <stack>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <set>
#include <array>
#include <tuple>
#include <limits>
#include <numeric>
#include <chrono>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/container/small_vector.hpp>
#include "io/hdf5/WriteVoronoiHDF5.hpp"
#include "io/vtk/WriteVoronoiVTK.hpp"
#include "utils/print/all.hpp"

#ifdef MADVORO_WITH_VCL
  #include <vectorclass.h>
#endif // MADVORO_WITH_VCL

#ifdef MADVORO_WITH_MPI
  #include <mpi.h>
#endif // MADVORO_WITH_MPI

#include "delaunay/Delaunay3D.hpp"
#include "geometry/Intersections.hpp"
#include "utils/Predicates3D.hpp"
#include "elementary/Face3D.hpp"
#include "elementary/Point3D.hpp"
#include "elementary/Mat33.hpp"
#include "utils/Predicates3D.hpp"
#include "utils/utils.hpp"
#include "io/io3D.hpp"
#include "exception/InvalidArgumentException.hpp"

#ifdef MADVORO_WITH_MPI
  #include "mpi/serialize/mpi_commands.hpp"
#endif

// finders
#include "hilbert/HilbertOrder3D.hpp"
#include "range/SmallRangeAgent.hpp"
#include "range/BigRangeAgent.hpp"
#include "range/finders/BruteForce.hpp"
#include "range/finders/RangeTree.hpp"
#include "range/finders/OctTree.hpp"
#include "range/finders/KDTree.hpp"
#include "range/finders/GroupRangeTree.hpp"

#ifdef MADVORO_WITH_MPI
  // env agents
  #include "environment/EnvironmentAgent.h"
  #include "environment/hilbert/DistributedOctEnvAgent.hpp"
  #include "environment/hilbert/HilbertTreeEnvAgent.hpp"
  #include "environment/hilbert/HilbertEnvAgent.hpp"
#endif // MADVORO_WITH_MPI

#ifdef MADVORO_WITH_MPI
  #include "voronoi/pointsManager/HilbertPointsManager.hpp"
  #define INITIAL_SENDRECV_TAG 1105
#endif 

#define LARGE_POINTS_SHRINK_RADIUS_RATIO 0.95
#define RANGE_MAX_POINTS_TO_GET 15 // 15
#define RADIUSES_GROWING_FACTOR 1.1
#define RADIUS_UNINITIALIZED -1

using namespace MadVoro;

/* ========================= IMPLEMENTATION ========================= */
namespace MadVoro
{  
    typedef std::array<std::size_t, 4> b_array_4;
    typedef std::array<std::size_t, 3> b_array_3;
    
    //! \brief Container for neighbouring tetrahedra
    typedef boost::container::small_vector<size_t, 40> tetra_vec;

    class Voronoi3D::Voronoi3DImpl
    {
        using AllPointsMap = boost::container::flat_map<size_t, size_t>;
        using IndexedPointsTree = DataStructure::OctTree<Range::IndexedPoint3D>;
    
    private:
        Point3D ll_, ur_;
        std::size_t Norg_, bigtet_;

        std::set<int> set_temp_;
        std::stack<int> stack_temp_;

        void FindIntersectionsSingle(vector<Face3D> const& box, std::size_t point, Geometry::Sphere<Point3D> &sphere,
                vector<size_t> &intersecting_faces, std::vector<double> &Rtemp, std::vector<Point3D> &vtemp);

        std::size_t GetFirstPointToCheck(void)const;

        void GetPointToCheck(std::size_t point, vector<unsigned char> const& checked, vector<std::size_t> &res);
        
        void CalcRigidCM(std::size_t face_index);

        void GetTetraCM(std::array<Point3D, 4> const& points, Point3D &CM) const;

        double GetTetraVolume(std::array<Point3D, 4> const& points)const;

        double GetRadius(const size_t &index) const;

        void CalcAllCM(void);

        vector<std::pair<std::size_t, std::size_t> > SerialFindIntersections(bool first_run);

        vector<std::pair<std::size_t, std::size_t> > SerialFirstIntersections(void);

        double CalcTetraRadiusCenterHiPrecision(const size_t &index) const;

        double CalcTetraRadiusCenter(const size_t &index) const;

        vector<Point3D> CreateBoundaryPoints(vector<std::pair<std::size_t, std::size_t> > const& to_duplicate,
                vector<vector<size_t> > &past_duplicate);
        void BuildVoronoi(std::vector<size_t> const& order);

        void InitialBoxBuild(std::vector<Face3D> &box, std::vector<Point3D> &normals);

        void BringSelfGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                    BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                    boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                                    boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints,
                                    boost::container::flat_set<size_t> &selfIgnorePoints);
        
        #ifdef MADVORO_WITH_MPI
        void BringGhostPointsToBuild(const MPI_Comm &comm);
        #else
        void BringGhostPointsToBuild();
        #endif // MADVORO_WITH_MPI

        std::pair<std::vector<SmallRangeQueryData>, std::vector<BigRangeQueryData>> CreateBatches(boost::container::flat_set<size_t> &smallPoints, boost::container::flat_set<size_t> &largePoints, const boost::container::flat_map<size_t, size_t> &firstLargeIteration, std::vector<double> &currentRadiuses, size_t iterations);

        std::pair<boost::container::flat_set<size_t>, boost::container::flat_set<size_t>>
        DetermineNextIterationPoints(size_t iterations,
                                        boost::container::flat_map<size_t, size_t> &firstLargeIteration,
                                        std::vector<double> &currentRadiuses,
                                        const boost::container::flat_map<size_t, size_t> &resultOfSmallPoints,
                                        const boost::container::flat_map<size_t, size_t> &resultOfBigPoints
                                    );

        void UpdateRadiuses(const std::vector<Point3D> &points);

        void UpdateCMs(void);
        
        void UpdateRangeFinder(void);

        #ifdef MADVORO_WITH_MPI
        std::vector<Point3D> PrepareToBuildParallel(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing);
        void FilterRealGhostPoints();
        void UpdateDuplicatedPoints(const std::vector<int> &sentProc, const std::vector<std::vector<size_t>> &sentPoints);
        void EnsureSymmetry(const std::vector<int> &sentProc, const std::vector<std::vector<int>> &recvProcLists);
        std::tuple<std::vector<Point3D>, std::vector<int>, std::vector<std::vector<size_t>>, std::vector<int>, std::vector<std::vector<size_t>>> InitialGhostPointsExchange(const MPI_Comm &comm = MPI_COMM_WORLD) const;
        void InitialExchange(const std::vector<Point3D> &points, std::vector<int> &sentProc, std::vector<std::vector<size_t>> &sentPoints, const MPI_Comm &comm = MPI_COMM_WORLD);
        void SetGhostArray(const std::vector<int> &recvProc, const std::vector<std::vector<size_t>> &recvPoints);  
        void BringRemoteGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                            BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                            boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                                            boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints);
        #endif // MADVORO_WITH_MPI

        Delaunay3D del_;
        //vector<vector<std::size_t> > PointTetras_; // The tetras containing each point
        vector<tetra_vec> PointTetras_; // The tetras containing each point
        mutable vector<double> R_; // The radius of the sphere of each tetra
        mutable vector<Point3D> tetra_centers_;
        // Voronoi Data
        //vector<vector<std::size_t> > FacesInCell_;
        vector<face_vec > FacesInCell_;
        std::vector<point_vec > PointsInFace_; // Right hand with regard to first neighbor
        //vector<vector<std::size_t> > PointsInFace_; // Right hand with regard to first neighbor
        vector<std::pair<std::size_t, std::size_t> > FaceNeighbors_;
        vector<Point3D> all_CM;
        vector<Point3D> CM_, Face_CM_; // center of masses
        vector<double> volume_; // volumes of each one of the tetrahedra
        vector<double> area_; // surface area of each one of the tetrahedra
        
        #ifdef MADVORO_WITH_MPI
        vector<int> sentprocs_;
        vector<vector<std::size_t>> sentpoints_; // if rank `i` is inside index `j` in `sentprocs_`, then the points in sentpoints_[j] are the points I sent to rank `i` in the initial points exchange in build
        vector<int> duplicatedprocs_; 
        vector<vector<std::size_t>> duplicated_points_;  // if rank `i` is inside index `j` in `duplicatedprocs_`, then Nghost_[j] includes all the points in `i`'s delaunay, which are actually mine
        vector<int> real_duplicated_proc;
        vector<vector<std::size_t>> real_duplicated_points; // indices of points which are a real ghost points
        vector<vector<std::size_t>> Nghost_; // if rank `i` is inside index `j` in `duplicatedprocs_`, then Nghost_[j] includes all the points in my delaunay, which are belongs, originally, to i
        vector<std::size_t> self_index_; // indexes of the points which are truely mine (inside the points list)
        #endif // MADVORO_WITH_MPI

        Voronoi3DImpl();
        Voronoi3DImpl(Voronoi3DImpl const &other);
        std::array<Point3D, 4> temp_points_;
        std::array<Point3D, 5> temp_points2_;
        std::vector<Face3D> box_faces_;

        std::shared_ptr<IndexedPointsTree> myPointsTree;
        std::shared_ptr<IndexedPointsTree> allMyPointsTree;
        #ifdef MADVORO_WITH_MPI
        std::shared_ptr<PointsManager> pointsManager;
        #endif // MADVORO_WITH_MPI

        std::shared_ptr<Range::RangeFinder> rangeFinder;
        std::vector<Point3D> allMyPoints;
        std::vector<double> allPointsWeights;
        std::vector<double> radiuses;

        AllPointsMap indicesInAllMyPoints; // the indices of the points in `del_.points_`, in the list of all points
        bool verbosity;

    public:
        #ifdef MADVORO_WITH_MPI
        const std::vector<double> &GetPointsBuildWeights() const;
        
        const EnvironmentAgent *GetEnvironmentAgent() const;
        
        std::vector<Point3D> BuildParallel(const std::vector<Point3D> &points, const std::vector<double> &weights, bool suppressRebalancing = false)
        {
            std::vector<size_t> indicesToBuild(points.size());
            std::iota(indicesToBuild.begin(), indicesToBuild.end(), 0);
            return this->BuildPartiallyParallel(points, weights, indicesToBuild, suppressRebalancing);
        }

        inline std::vector<Point3D> BuildParallel(const std::vector<Point3D> &points, bool suppressRebalancing = false)
        {
            return this->BuildParallel(points, std::vector<double>(points.size(), 1.0), suppressRebalancing);
        }
        #endif // MADVORO_WITH_MPI

        #ifdef MADVORO_WITH_MPI
        vector<int>& GetSentProcs(void);

        vector<vector<size_t> >& GetSentPoints(void);

        vector<size_t>& GetSelfIndex(void);
        #endif // MADVORO_WITH_MPI

        vector<Point3D>& GetAllFaceCM(void);

        const vector<Point3D>& GetAllFaceCM(void) const;

        const Point3D &FaceCM(std::size_t index)const;

        Voronoi3DImpl(Point3D const& ll, Point3D const& ur);

        Voronoi3DImpl(std::vector<Face3D> const& box_faces);

        void output(std::string const& filename)const;

        void BuildInitialize(size_t num_points);

        void BuildPartially(const std::vector<Point3D> &allPoints, const std::vector<size_t> &indicesToBuild);

        void Build(const std::vector<Point3D> &points);

    #ifdef MADVORO_WITH_MPI
        /*! \brief Output extra build
        \param filename Output file name
        */
        void output_buildextra(std::string const& filename) const;

        void PreparePoints(const std::vector<Point3D> &points, const std::vector<size_t> &mask);

        std::vector<Point3D> BuildPartiallyParallel(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing = false);

        bool PointInMyDomain(const Point3D &point) const;

        int GetOwner(const Point3D &point) const;
    #endif // MADVORO_WITH_MPI

        void BuildDebug(int rank);

        double GetMaxRadius(const size_t &index) const;

        double GetMinRadius(const size_t &index) const;

        size_t GetContainingCell(const Point3D &point) const;

        std::size_t GetPointNo(void) const;

        const Point3D &GetMeshPoint(std::size_t index) const;

        double GetArea(std::size_t index) const;

        Point3D const& GetCellCM(std::size_t index) const;

        std::size_t GetTotalFacesNumber(void) const;

        double GetWidth(std::size_t index) const;

        double GetVolume(std::size_t index) const;

        face_vec const& GetCellFaces(std::size_t index) const;
        
        vector<Point3D>& accessMeshPoints(void);

        const vector<Point3D>& getMeshPoints(void) const;

        const AllPointsMap &GetIndicesInAllPoints(void) const;

        const std::vector<Point3D> &getAllPoints(void) const;

        std::vector<Point3D> &getAllPoints(void);

        size_t GetAllPointsNo(void) const;

        vector<std::size_t> GetNeighbors(std::size_t index)const;

        Voronoi3DImpl* clone(void) const;

        bool NearBoundary(std::size_t index) const;

        bool BoundaryFace3D(std::size_t index) const;

        #ifdef MADVORO_WITH_MPI
        vector<vector<std::size_t> >& GetDuplicatedPoints(void);

        vector<vector<std::size_t> >const& GetDuplicatedPoints(void)const;

        vector<int> GetDuplicatedProcs(void)const;

        vector<int> GetSentProcs(void)const;

        vector<vector<std::size_t> > const& GetSentPoints(void)const;

        vector<std::size_t> const& GetSelfIndex(void) const;

        #endif // MADVORO_WITH_MPI
        std::size_t GetTotalPointNumber(void)const;

        vector<Point3D> & GetAllCM(void);

        vector<Point3D > GetAllCM(void)const;

        void GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point)const;

        Point3D Normal(std::size_t faceindex)const;

        bool IsGhostPoint(std::size_t index)const;

        Point3D CalcFaceVelocity(std::size_t index, Point3D const& v0, Point3D const& v1)const;

        vector<Point3D>& GetFacePoints(void);

        vector<double>& GetAllArea(void);

        vector<Point3D>const& GetFacePoints(void) const;

        vector<face_vec >& GetAllCellFaces(void);

        vector<face_vec >const& GetAllCellFaces(void) const;

        point_vec const& GetPointsInFace(std::size_t index) const;

        const std::pair<std::size_t, std::size_t> &GetFaceNeighbors(std::size_t face_index) const;

        #ifdef MADVORO_WITH_MPI
        vector<vector<std::size_t> > const& GetGhostIndeces(void) const;

        vector<vector<std::size_t> >& GetGhostIndeces(void);
        #endif // MADVORO_WITH_MPI

        void GetNeighbors(size_t index, vector<size_t> &res) const;

        std::pair<Point3D, Point3D> GetBoxCoordinates(void) const;

        void BuildNoBox(vector<Point3D> const& points, vector<vector<Point3D> > const& ghosts,vector<size_t> toduplicate);

        vector<double>& GetAllVolumes(void);

        vector<double> GetAllVolumes(void)const;

        std::vector<std::pair<size_t, size_t>> &GetAllFaceNeighbors(void);

        const std::vector<std::pair<size_t, size_t>> &GetAllFaceNeighbors(void) const;

        vector<point_vec > & GetAllPointsInFace(void);

        vector<point_vec > const& GetAllPointsInFace(void) const;

        size_t& GetPointNo(void);

        bool IsPointOutsideBox(size_t index) const;

        void SetBox(Point3D const& ll, Point3D const& ur);

        std::vector<Face3D> GetBoxFaces(void) const {return box_faces_;}

        std::vector<Face3D>& ModifyBoxFaces(void) {return box_faces_;}

        template<typename T>
        void SyncPartialBuildData(std::vector<T> &partialBuildData, std::vector<T> &allBuildData) const;
    
        inline void SetVerbosity(bool value){this->verbosity = value;}

        bool PointInPolyTess(Point3D const &point, std::size_t index);
    };
}

  template<typename T>
  inline void MadVoro::Voronoi3D::Voronoi3DImpl::SyncPartialBuildData(std::vector<T> &partialBuildData, std::vector<T> &allBuildData) const
  {
    size_t Norg = this->GetPointNo();
    if(partialBuildData.size() < Norg)
    {
      MadVoro::Exception::MadVoroException eo("Voronoi3D::SyncPartialBuildData: Partial build data has lower size than the number of points");
      eo.addEntry("Partial build data size", partialBuildData.size());
      eo.addEntry("Number of points", Norg);
      throw eo;
    }
    const Voronoi3DImpl::AllPointsMap &indicesInAllMyPoints = this->GetIndicesInAllPoints();
    
    allBuildData.resize(this->GetAllPointsNo());

    // update CMs of active local points in all points CM vector
    for(size_t i = 0; i < Norg; i++)
    {
        size_t pointIdx = indicesInAllMyPoints.at(i);
        allBuildData[pointIdx] = partialBuildData[i];
    }

    // update the CM of local non active points
    size_t sizeOfMeshPoints = this->getMeshPoints().size();
    partialBuildData.resize(sizeOfMeshPoints);
    for(size_t i = Norg; i < sizeOfMeshPoints; i++)
    {
        // check if the point is mine. i.e, appears in `indicesInAllMyPoints`. Just copy the CM from there.
        bool pointIsMine = (indicesInAllMyPoints.find(i) != indicesInAllMyPoints.cend());
        if(pointIsMine)
        {
            size_t pointIdx = indicesInAllMyPoints.at(i);
            partialBuildData[i] = allBuildData[pointIdx];
        }
    }

    #ifdef MADVORO_WITH_MPI
        // update the CM of active and not active, but non local points

        std::vector<std::vector<T>> incoming = MPI::MPI_exchange_data_indexed(this->GetDuplicatedProcs(), allBuildData, this->GetDuplicatedPoints());
        size_t incomingSize = incoming.size();
        const std::vector<std::vector<size_t>> &Nghost = this->GetGhostIndeces();
        assert(this->GetDuplicatedProcs().size() == Nghost.size());
        assert(incomingSize == Nghost.size());
        for (size_t i = 0; i < incomingSize; ++i)
        {
            size_t _size = incoming[i].size();
            assert(_size == Nghost[i].size());
            for (size_t j = 0; j < _size; ++j)
            {
                partialBuildData[Nghost.at(i).at(j)] = incoming[i][j];
            }
        }
    #endif // MADVORO_WITH_MPI
  }


bool MadVoro::Voronoi3D::Voronoi3DImpl::PointInPolyTess(Point3D const &point, std::size_t index)
{
    face_vec const &faces = this->GetCellFaces(index);
    vector<Point3D> const &points = this->GetFacePoints();
    std::size_t N = faces.size();
    std::array<Point3D, 4> vec;
    for (std::size_t i = 0; i < N; ++i)
    {
        double R = Utils::fastsqrt(this->GetArea(faces[i]));
        size_t N1 = 0;
        size_t N2 = 0;
        Point3D V1, V2;
        size_t counter = 0;
        point_vec const &InFace3D = this->GetPointsInFace(faces[i]);
        size_t NinFace3D = InFace3D.size();
        N1 = 1;
        V1 = points[InFace3D[(counter + 1) % NinFace3D]] - points[InFace3D[0]];
        while (fastabs(V1) < 0.01 * R)
        {
            ++counter;
            assert(counter < NinFace3D);
            V1 = points[InFace3D[(counter + 1) % NinFace3D]] - points[InFace3D[0]];
            ++N1;
        }
        V2 = points[InFace3D[(counter + 2) % NinFace3D]] - points[InFace3D[N1]];
        N2 = (counter + 2) % NinFace3D;
        while (fastabs(V2) < 0.01 * R || fastabs(CrossProduct(V1, V2)) < 0.0001 * this->GetArea(faces[i]))
        {
            ++counter;
            if (counter > 2 * NinFace3D)
                break;
            V2 = points[InFace3D[(counter + 2) % NinFace3D]] - points[InFace3D[N1]];
            N2 = (counter + 2) % NinFace3D;
        }
        if (counter > 2 * NinFace3D)
        {
            std::cout << "Weird face in PointInPoly, cell " << index << " face " << faces[i] << " i " << i << " face area " << this->GetArea(faces[i]) << std::endl;
            for (size_t j = 0; j < NinFace3D; ++j)
                std::cout << "Point j " << points[InFace3D[j]].x << "," << points[InFace3D[j]].y << "," << points[InFace3D[j]].z << std::endl;
            Point3D normal = this->GetFaceNeighbors(faces[i]).second == index ? this->GetMeshPoint(this->GetFaceNeighbors(faces[i]).second) - this->GetMeshPoint(this->GetFaceNeighbors(faces[i]).first) : this->GetMeshPoint(this->GetFaceNeighbors(faces[i]).first) - this->GetMeshPoint(this->GetFaceNeighbors(faces[i]).second);
            if (ScalarProd(normal, point - points[InFace3D[0]]) < 0)
                return false;
        }
        else
        {
            vec[0] = points[InFace3D[0]];
            vec[1] = points[InFace3D.at(N1)];
            vec[2] = points[InFace3D.at(N2)];
            vec[3] = this->GetMeshPoint(index);
            double s1 = orient3d(vec);
            vec[3] = point;
            double s2 = orient3d(vec);
            if (s1 * s2 < -0)
                return false;
        }
    }
    return true;
}

namespace MadVoro
{
    bool PointInPoly(std::vector<Face3D> const& faces, Point3D const &point)
    {
        std::size_t const N = faces.size();
        std::array<Point3D, 4> vec;
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
}

namespace
{
#ifdef MADVORO_WITH_MPI
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

    vector<Face3D> BuildBox(Point3D const &ll, Point3D const &ur)
    {
        double dx = ur.x - ll.x;
        double dy = ur.y - ll.y;
        double dz = ur.z - ll.z;
        vector<Face3D> res(6);
        vector<Point3D> points;
        points.push_back(ll);
        points.push_back(ll + Point3D(dx, 0, 0));
        points.push_back(ll + Point3D(dx, dy, 0));
        points.push_back(ll + Point3D(0, dy, 0));
        points.push_back(ll + Point3D(0, 0, dz));
        points.push_back(ll + Point3D(dx, 0, dz));
        points.push_back(ll + Point3D(dx, dy, dz));
        points.push_back(ll + Point3D(0, dy, dz));
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

#ifdef MADVORO_WITH_MPI
    vector<Point3D> GetBoxNormals(Point3D const &ll, Point3D const &ur, vector<Face3D> const& box_faces_)
    {
        const vector<Face3D> faces = box_faces_.empty() ? BuildBox(ll, ur) : box_faces_;
        vector<Point3D> res(faces.size());
        size_t N = res.size();
        for (size_t i = 0; i < N; ++i)
        {
            CrossProduct(faces[i].vertices[2] - faces[i].vertices[0], faces[i].vertices[1] - faces[i].vertices[0], res[i]);
            res[i] *= 1.0 /abs(res[i]);
        }
        return res;
    }

    size_t BoxIndex(vector<Point3D> const &fnormals, const Point3D &normal)
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

    double CleanDuplicates(std::array<size_t, 128> const &indeces, const vector<Point3D> &points,
                                                 boost::container::small_vector<size_t, 8> &res, double R,
                                                 std::array<double, 128> &diffs,
                                                 std::array<Point3D, 128> &vtemp, const size_t N)
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

    size_t SetPointTetras(vector<tetra_vec> &PointTetras, size_t Norg, vector<Tetrahedron> &tetras,
                                                boost::container::flat_set<size_t> const &empty_tetras)
    {
        PointTetras.clear();
        PointTetras.resize(Norg);

        #ifdef MADVORO_WITH_VCL
            Vec4uq _Norg(Norg);
        #endif // MADVORO_WITH_VCL

        size_t Ntetra = tetras.size();
        size_t bigtet(0);
        bool has_good, has_big;
        // change empty tetras to be not relevant
        for (boost::container::flat_set<size_t>::const_iterator it = empty_tetras.begin(); it != empty_tetras.end(); ++it)
        {
#ifdef __INTEL_COMPILER
#pragma omp simd early_exit
#endif
            for (size_t i = 0; i < 4; ++i)
            {
                tetras[*it].points[i] = std::numeric_limits<std::size_t>::max();
                tetras[*it].neighbors[i] = std::numeric_limits<std::size_t>::max();
            }
        }

        for (size_t i = 0; i < Ntetra; ++i)
        {
            const Tetrahedron &tet = tetras[i];

            has_good = false;
            has_big = false;

            #ifdef MADVORO_WITH_VCL
                Vec4uq _points(tet.points[0], tet.points[1], tet.points[2], tet.points[3]);
                Vec4qb cmp = (_points < _Norg);
            #endif // MADVORO_WITH_VCL
            
            for(int j = 0; j < 4; ++j)
            {
                #ifdef MADVORO_WITH_VCL
                    if(cmp[j])
                #else // MADVORO_WITH_VCL
                    if(tet.points[j] < Norg)
                #endif // MADVORO_WITH_VCL
                {
                    has_good = true;
                    PointTetras[tet.points[j]].push_back(i);
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
        return bigtet;
    }

    bool CleanSameLine(boost::container::small_vector<size_t, 8> &indeces, vector<Point3D> const& face_points, std::array<double, 128> &area_vec_temp)
    {
        point_vec old;
        size_t const N = indeces.size();
        double const small_fraction = 1e-14;
        // double const medium_fraction = 3e-1;
        // Find correct normal
        Point3D good_normal;
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
            Point3D normal_temp = CrossProduct(face_points[indeces[i]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]], face_points[indeces[(i + 1) % Nindeces]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]]);
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
                Point3D normal_temp = CrossProduct(face_points[indeces[i]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]], face_points[indeces[(i + 1) % Nindeces]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]]);
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
            MadVoro::Exception::MadVoroException eo("Bad CleanSameLine");
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
                Point3D normal_temp = CrossProduct(face_points[indeces[i]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]], face_points[indeces[(i + 1) % Nindeces]] - face_points[indeces[(Nindeces + i - 1) % Nindeces]]);
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

    void MakeRightHandFace(boost::container::small_vector<size_t, 8> &indeces, Point3D const &point, vector<Point3D> const &face_points,
                                                 std::array<size_t, 128> &temp, double areascale)
    {
        Point3D V1, V2;
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
            throw MadVoro::Exception::MadVoroException("Bad NextLoopTetra");
        return cur_tetra.neighbors[i];
    }

    void CalcFaceAreaCM(boost::container::small_vector<size_t, 8> const &indeces, std::vector<Point3D> const &allpoints,
                                            std::array<Point3D, 128> &points, double &Area, Point3D &CM,
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
        //Point3D temp3, temp4, temp5;
        for (size_t i = 0; i < Nloop; i++)
        {
            //temp4.Set(points[i + 1].x - points[0].x, points[i + 1].y - points[0].y, points[i + 1].z - points[0].z);
            Point3D temp4(points[i + 1].x - points[0].x, points[i + 1].y - points[0].y, points[i + 1].z - points[0].z);
            //temp5.Set(points[i + 2].x - points[0].x, points[i + 2].y - points[0].y, points[i + 2].z - points[0].z);
            Point3D temp5(points[i + 2].x - points[0].x, points[i + 2].y - points[0].y, points[i + 2].z - points[0].z);
            Point3D temp3;
            CrossProduct(temp4, temp5, temp3);
            Atemp[i] = 0.3333333333333333 * 0.5 * Utils::fastsqrt(ScalarProd(temp3, temp3));
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

    bool PointInDomain(Point3D const &ll, Point3D const &ur, Point3D const &point)
    {
        if (point.x > ll.x && point.x < ur.x && point.y > ll.y && point.y < ur.y && point.z > ll.z && point.z < ur.z)
            return true;
        else
            return false;
    }

    Point3D MirrorPoint(Face3D const &face, Point3D const &point)
    {
        Point3D normal = CrossProduct(face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]);
        normal = normal / abs(normal);
        return point - (2 * ScalarProd(point - face.vertices[0], normal)) * normal;
    }
}

MadVoro::Voronoi3D::Voronoi3DImpl::Voronoi3DImpl(std::vector<Face3D> const& box_faces) : Voronoi3DImpl()
{
    this->box_faces_ = box_faces;
    size_t const Nfaces = box_faces.size();
    if(Nfaces < 4)
        throw MadVoro::Exception::InvalidArgumentException("Zero face vector in Voronoi3D constructor");
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

MadVoro::Voronoi3D::Voronoi3DImpl::Voronoi3DImpl(Point3D const &ll, Point3D const &ur)
{
    this->ll_ = ll;
    this->ur_ = ur;
    this->Norg_ = 0;
    this->bigtet_ = 0;
    this->set_temp_ = std::set<int>();
    this->stack_temp_ = std::stack<int>();
    this->del_ = Delaunay3D();
    this->PointTetras_ = vector<tetra_vec>();
    this->R_ = vector<double>();
    this->tetra_centers_ = vector<Point3D>();
    this->FacesInCell_ = vector<face_vec>();
    this->PointsInFace_ = vector<point_vec>();
    this->FaceNeighbors_ = vector<std::pair<std::size_t, std::size_t>>();
    this->CM_ = vector<Point3D>();
    this->Face_CM_ = vector<Point3D>();
    this->volume_ = vector<double>();
    this->area_ = vector<double>();
    #ifdef MADVORO_WITH_MPI
    this->sentprocs_ = vector<int>();
    this->sentpoints_ = std::vector<std::vector<size_t>>();
    this->duplicatedprocs_ = std::vector<int>();
    this->duplicated_points_ = std::vector<std::vector<size_t>>();
    this->Nghost_ = std::vector<std::vector<size_t>>();
    this->self_index_ = std::vector<std::size_t>();
    #endif // MADVORO_WITH_MPI
    this->temp_points_ = std::array<Point3D, 4>();
    this->temp_points2_ = std::array<Point3D, 5>();
    this->box_faces_ = std::vector<Face3D>();
    #ifdef MADVORO_WITH_MPI
    this->pointsManager = std::shared_ptr<PointsManager>();
    this->allMyPoints = std::vector<Point3D>();
    #endif // MADVORO_WITH_MPI
    this->indicesInAllMyPoints = AllPointsMap();
    this->verbosity = false;
}

MadVoro::Voronoi3D::Voronoi3DImpl::Voronoi3DImpl() : Voronoi3DImpl(Point3D(), Point3D())
{}

void MadVoro::Voronoi3D::Voronoi3DImpl::CalcRigidCM(std::size_t face_index)
{
    Point3D normal = normalize(del_.points_[FaceNeighbors_[face_index].first] - del_.points_[FaceNeighbors_[face_index].second]);
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

vector<Point3D> MadVoro::Voronoi3D::Voronoi3DImpl::CreateBoundaryPoints(vector<std::pair<std::size_t, std::size_t>> const &to_duplicate,
                                                 vector<vector<size_t>> &past_duplicate)
{
    size_t Ncheck = to_duplicate.size();
    vector<std::pair<std::size_t, std::size_t>> to_add;
    to_add.reserve(Ncheck);
    vector<Face3D> faces = box_faces_.empty() ? BuildBox(ll_, ur_) : box_faces_;
    vector<Point3D> res;
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

#ifdef MADVORO_WITH_MPI
    vector<vector<std::size_t>> const &MadVoro::Voronoi3D::Voronoi3DImpl::GetGhostIndeces(void) const
    {
        return Nghost_;
    }
#endif // MADVORO_WITH_MPI

/**
 * gets a point index, and returns the maximal radius of the tetrahedra containing that point.
 * @param index the index of the point (within the points list)
*/
double MadVoro::Voronoi3D::Voronoi3DImpl::GetMaxRadius(const size_t &index) const
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
double MadVoro::Voronoi3D::Voronoi3DImpl::GetMinRadius(const size_t &index) const
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
void MadVoro::Voronoi3D::Voronoi3DImpl::InitialBoxBuild(std::vector<Face3D> &box, std::vector<Point3D> &normals)
{
    box = box_faces_.empty() ? BuildBox(this->ll_, this->ur_) : this->box_faces_;
    size_t Nfaces = box.size();
    normals.resize(Nfaces);

    // calculates the normals for each one of the box's faces
    for (size_t i = 0; i < Nfaces; ++i)
    {
        normals[i] = CrossProduct(box[i].vertices[1] - box[i].vertices[0], box[i].vertices[2] - box[i].vertices[0]);
        normals[i] *= (1.0 / Utils::fastsqrt(ScalarProd(normals[i], normals[i])));
    }
}

/**
 * \author Maor Mizrachi
 * \brief Initializes internal data structures, for the voronoi build
*/
void MadVoro::Voronoi3D::Voronoi3DImpl::BuildInitialize(size_t num_points)
{
    // assert(num_points > 0);
    // Clear data
    PointTetras_.clear();
    R_.clear();
    if(num_points > 0) R_.reserve(num_points * 11);
    tetra_centers_.clear();
    if(num_points > 0) tetra_centers_.reserve(num_points * 11);
    // Voronoi Data
    del_.Clean();
    FacesInCell_.clear();
    PointsInFace_.clear();
    FaceNeighbors_.clear();
    CM_.clear();
    Face_CM_.clear();
    volume_.clear();
    area_.clear();
    Norg_ = num_points;
    #ifdef MADVORO_WITH_MPI
        duplicatedprocs_.clear();
        duplicated_points_.clear();
        Nghost_.clear();
    #endif // MADVORO_WITH_MPI
}

#ifdef MADVORO_WITH_MPI
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
bool MadVoro::Voronoi3D::Voronoi3DImpl::PointInMyDomain(const Point3D &point) const
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(this->pointsManager.get() != nullptr);
    assert(this->pointsManager->getEnvironmentAgent() != nullptr);
    return (this->GetOwner(point) == rank);
}

inline int MadVoro::Voronoi3D::Voronoi3DImpl::GetOwner(const Point3D &point) const
{
    return this->pointsManager->getEnvironmentAgent()->getOwner(point);
}

std::tuple<std::vector<Point3D>, std::vector<int>, std::vector<std::vector<size_t>>, std::vector<int>, std::vector<std::vector<size_t>>>
    MadVoro::Voronoi3D::Voronoi3DImpl::InitialGhostPointsExchange(const MPI_Comm &comm) const
{
    int size;
    MPI_Comm_size(comm, &size);

    std::vector<int> toSendSizes(size);
    std::vector<int> sendDisplacements(size);
    std::vector<_3DPoint> pointsToSend;

    std::vector<int> sentProcs;
    std::vector<std::vector<size_t>> sentPoints;

    for(int _rank = 0; _rank < size; _rank++)
    {
        size_t rankIndex = std::distance(this->real_duplicated_proc.begin(), std::find(this->real_duplicated_proc.begin(), this->real_duplicated_proc.end(), _rank));
        if(rankIndex == this->real_duplicated_proc.size())
        {
            // rank _rank is not duplicated
            toSendSizes[_rank] = 0;
        }
        else
        {
            std::vector<size_t> sentIndices;
            // rank _rank is duplicated
            toSendSizes[_rank] = 0;
            for(size_t pointIdx : this->real_duplicated_points[rankIndex])
            {
                if(pointIdx < this->Norg_)
                {
                    sentIndices.push_back(pointIdx);
                    toSendSizes[_rank]++;
                    pointsToSend.push_back(_3DPoint(this->allMyPoints[pointIdx]));
                }
            }
            if(not sentIndices.empty())
            {
                sentProcs.push_back(_rank);
                sentPoints.emplace_back(sentIndices);
            }
        }
        toSendSizes[_rank] *= sizeof(_3DPoint);
        if(_rank == 0)
        {
            sendDisplacements[_rank] = 0;
        }
        else
        {
            sendDisplacements[_rank] = sendDisplacements[_rank - 1] + toSendSizes[_rank - 1];
        }
    }

    std::vector<int> toRecvSizes(size);
    MPI_Alltoall(toSendSizes.data(), 1, MPI_INT, toRecvSizes.data(), 1, MPI_INT, comm);
    std::vector<int> recvDisplacements(size);
    size_t totalSize = 0;
    std::vector<int> recvProcs;
    std::vector<std::vector<size_t>> recvPoints;

    for(int _rank = 0; _rank < size; _rank++)
    {
        if(_rank == 0)
        {
            recvDisplacements[_rank] = 0;
        }
        else
        {
            recvDisplacements[_rank] = recvDisplacements[_rank - 1] + toRecvSizes[_rank - 1];
        }
        size_t receiving = toRecvSizes[_rank] / sizeof(_3DPoint);
        if(receiving > 0)
        {
            std::vector<size_t> receivedIndices; 
            for(size_t i = 0; i < receiving; i++)
            {
                receivedIndices.push_back(totalSize);
                totalSize++;
            }
            recvProcs.push_back(_rank);
            recvPoints.emplace_back(receivedIndices);
        }
    }
    std::vector<_3DPoint> almostGhostPoints(totalSize);
    MPI_Alltoallv(pointsToSend.data(), toSendSizes.data(), sendDisplacements.data(), MPI_BYTE, almostGhostPoints.data(), toRecvSizes.data(), recvDisplacements.data(), MPI_BYTE, comm);
    
    std::vector<Point3D> ghostPoints;
    ghostPoints.reserve(almostGhostPoints.size());
    for(const _3DPoint &point : almostGhostPoints)
    {
        ghostPoints.emplace_back(point.x, point.y, point.z);
    }
    return std::tuple(ghostPoints, sentProcs, sentPoints, recvProcs, recvPoints);
}

void MadVoro::Voronoi3D::Voronoi3DImpl::FilterRealGhostPoints()
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
void MadVoro::Voronoi3D::Voronoi3DImpl::UpdateDuplicatedPoints(const std::vector<int> &sentProc, const std::vector<std::vector<size_t>> &sentPoints)
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
void MadVoro::Voronoi3D::Voronoi3DImpl::EnsureSymmetry(const std::vector<int> &sentProc, const std::vector<std::vector<int>> &recvProcLists)
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

void MadVoro::Voronoi3D::Voronoi3DImpl::InitialExchange(const std::vector<Point3D> &points, std::vector<int> &sentProc, std::vector<std::vector<size_t>> &sentPoints, const MPI_Comm &comm)
{
    const EnvironmentAgent *envAgent = this->pointsManager->getEnvironmentAgent();
    bool supportsFurthestClosestRanks;
    std::function<HilbertCurveEnvironmentAgent::DistancesVector(const _3DPoint&)> getFurthestClosestRanks;

    // check if has 'smartAgent' (an agent that can caluclate distances of ranks as well)            
    const DistributedOctEnvironmentAgent *distribuedOctEnvAgent = dynamic_cast<const DistributedOctEnvironmentAgent*>(envAgent);
    if(distribuedOctEnvAgent != nullptr)
    {
        supportsFurthestClosestRanks = true;
        getFurthestClosestRanks = [distribuedOctEnvAgent](const _3DPoint &point){return distribuedOctEnvAgent->getClosestFurthestPointsByRanks(point);};
    }
    const HilbertTreeEnvironmentAgent *hilbertTreeEnvAgent = dynamic_cast<const HilbertTreeEnvironmentAgent*>(envAgent);
    if(hilbertTreeEnvAgent != nullptr)
    {
        supportsFurthestClosestRanks = true;
        getFurthestClosestRanks = [hilbertTreeEnvAgent](const _3DPoint &point){return hilbertTreeEnvAgent->getClosestFurthestPointsByRanks(point);};
    }

    if(not supportsFurthestClosestRanks)
    {
        return;
    }
    
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    size_t counter = 0;

    for(size_t pointIdx = 0; pointIdx < points.size(); pointIdx++)
    {
        bool isBorderPoint = false;
        for(const size_t &tetraIdx : this->PointTetras_[pointIdx])
        {
            const Tetrahedron &tet = this->del_.tetras_[tetraIdx];
            isBorderPoint = (tet.points[0] >= this->Norg_) or (tet.points[1] >= this->Norg_) or (tet.points[2] >= this->Norg_) or (tet.points[3] >= this->Norg_);
            if(isBorderPoint)
            {
                break;
            }
        }
        if(!isBorderPoint)
        {
            continue;
        }
        int closestRank = std::numeric_limits<int>::max();
        double closestDistance = std::numeric_limits<double>::max();
        auto distances = getFurthestClosestRanks(points[pointIdx]);
        for(int _rank = 0; _rank < size; _rank++)
        {
            if(_rank == rank)
            {
                continue;
            }
            if(distances[_rank].first < closestDistance)
            {
                closestDistance = distances[_rank].first;
                closestRank = _rank;
            }
        }
        size_t rankIdx = std::distance(sentProc.begin(), std::find(sentProc.begin(), sentProc.end(), closestRank));
        if(rankIdx == sentProc.size())
        {
            sentProc.push_back(closestRank);
            sentPoints.emplace_back(std::vector<size_t>());
        }
        sentPoints[rankIdx].push_back(pointIdx);
        counter++;
    }

    std::vector<size_t> sendLengths(size, 0);
    std::vector<std::vector<_3DPoint>> toSend;
    toSend.resize(sentProc.size());
    
    std::vector<MPI_Request> requests;
    requests.reserve(4 * sentProc.size()); // heuristic

    std::vector<size_t> recvLengths(size, 0);

    for(size_t i = 0; i < sentProc.size(); i++)
    {
        int _rank = sentProc[i];
        sendLengths[_rank] = sentPoints[i].size();
        toSend[i].reserve(sentPoints[i].size());
        for(size_t &pointIdx : sentPoints[i])
        {
            toSend[i].emplace_back(_3DPoint(points[pointIdx].x, points[pointIdx].y, points[pointIdx].z));
        }
    }

    MPI_Alltoall(&sendLengths[0], sizeof(size_t), MPI_BYTE, &recvLengths[0], sizeof(size_t), MPI_BYTE, comm);

    size_t totalLength = 0; 
    for(int _rank = 0; _rank < size; _rank++)
    {
        if(recvLengths[_rank] > 0)
        {
            totalLength += recvLengths[_rank];
            size_t rankIdx = std::distance(sentProc.begin(), std::find(sentProc.begin(), sentProc.end(), _rank));
            if(rankIdx != sentProc.size())
            {
                // rank has already been found
                continue;
            }
            sentProc.push_back(_rank);
            sentPoints.emplace_back(std::vector<size_t>());
        }
    }

    std::vector<_3DPoint> almostExtraPoints;
    almostExtraPoints.resize(totalLength);
    size_t insertedSoFar = 0; 
    for(const int &_rank : sentProc)
    {
        if(recvLengths[_rank] > 0)
        {
            // std::cout << "rank " << rank << " is receiving " << recvLengths[_rank] << " from rank " << _rank << ", insertedSoFar is " << insertedSoFar << "(total length: " << totalLength << ")" << std::endl;
            requests.push_back(MPI_REQUEST_NULL);
            MPI_Irecv(&almostExtraPoints[insertedSoFar], sizeof(_3DPoint) * recvLengths[_rank], MPI_BYTE, _rank, INITIAL_SENDRECV_TAG, comm, &requests[requests.size() - 1]);
            size_t dupRankIdx = std::distance(this->duplicatedprocs_.begin(), std::find(this->duplicatedprocs_.begin(), this->duplicatedprocs_.end(), _rank));
            if(dupRankIdx == this->duplicatedprocs_.size())
            {
                // new rank in this->duplicatedprocs_, initialize it
                this->duplicatedprocs_.push_back(_rank);
                this->duplicated_points_.emplace_back(std::vector<size_t>());
                this->Nghost_.emplace_back(std::vector<size_t>());
            }
            for(size_t i = 0; i < recvLengths[_rank]; i++)
            {
                // batchInfo.pointsFromRanks[_rank][i] holds an index of point, but this point will be added to my delaunay, so
                // its index there will be this->del_.points_.size() + batchInfo.pointsFromRanks[_rank][i]
                this->Nghost_[dupRankIdx].push_back(this->del_.points_.size() + insertedSoFar + i);
            }
        }
    }

    for(size_t i = 0; i < toSend.size(); i++)
    {
        int _rank = sentProc[i];
        // std::cout << "rank " << rank << " is sending " << toSend[i].size() << " to rank " << _rank << std::endl;
        requests.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&toSend[i][0], sizeof(_3DPoint) * toSend[i].size(), MPI_BYTE, _rank, INITIAL_SENDRECV_TAG, comm, &requests[requests.size() - 1]);
    }

    if(!requests.empty())
    {
        MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
    }   

    std::vector<Point3D> extraPoints;
    for(const _3DPoint &_point : almostExtraPoints)
    {
        extraPoints.emplace_back(Point3D(_point.x, _point.y, _point.z));
    }

    this->del_.BuildExtra(extraPoints);

    this->R_.resize(this->del_.tetras_.size());
    std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
    this->tetra_centers_.resize(this->R_.size());
    this->bigtet_ = SetPointTetras(this->PointTetras_, this->Norg_, this->del_.tetras_, this->del_.empty_tetras_);
}


/**
 * \author Maor Mizrachi
 * \brief Sets the ghost points arrays (duplicatedprocs_, duplicated_points_, Nghost_)
*/
void MadVoro::Voronoi3D::Voronoi3DImpl::SetGhostArray(const std::vector<int> &recvProc, const std::vector<std::vector<size_t>> &recvPoints)
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
std::vector<Point3D> MadVoro::Voronoi3D::Voronoi3DImpl::PrepareToBuildParallel(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing)
{
    if(this->radiuses.size() < allPoints.size())
    {
        this->radiuses.resize(allPoints.size(), RADIUS_UNINITIALIZED);
    }
    if(this->all_CM.size() < allPoints.size())
    {
        this->all_CM.resize(allPoints.size());
    }
    
    if(this->pointsManager.get() == nullptr)
    {
        // initialize points manager
        this->pointsManager = std::shared_ptr<HilbertPointsManager>(new HilbertPointsManager(this->ll_, this->ur_));
    }

    int canDoRebalance = ((not suppressRebalancing) and (indicesToBuild.size() == allPoints.size()))? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &canDoRebalance, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    bool allowRebalance = (canDoRebalance == 1);

    PointsExchangeResult exchangeResult = this->pointsManager->update(allPoints, allWeights, indicesToBuild, this->radiuses, this->all_CM, allowRebalance); // does rebalancing (if necessary) and exchanging

    this->allMyPoints = std::move(exchangeResult.newPoints);
    this->radiuses = std::move(exchangeResult.newRadiuses);
    this->all_CM = std::move(exchangeResult.newCMs);
    this->sentprocs_ = std::move(exchangeResult.sentProcessors);
    this->sentpoints_ = std::move(exchangeResult.sentIndicesToProcessors);
    this->self_index_ = std::move(exchangeResult.indicesToSelf);
    this->allPointsWeights = std::move(exchangeResult.newWeights);

    assert(this->allMyPoints.size() == this->allPointsWeights.size());

    std::vector<Point3D> new_points;
    this->indicesInAllMyPoints = AllPointsMap();
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
void MadVoro::Voronoi3D::Voronoi3DImpl::PreparePoints(const std::vector<Point3D> &points, const std::vector<size_t> &mask)
{
    if(points.size() != mask.size())
    {
        MadVoro::Exception::InvalidArgumentException eo("In MadVoro::Voronoi3D::PreparePoints, mask size is not equal to the points size");
        eo.addEntry("Mask size", mask.size());
        eo.addEntry("Points size", points.size());
        throw eo;
    }
    size_t originalPointsNum = this->allMyPoints.size();
    size_t newPointsNum = points.size();

    std::vector<Range::IndexedPoint3D> oldPoints;
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
        IndexedPointsTree oldPointsTree(this->ll_, this->ur_, oldPoints);
        newRadiuses = std::vector<double>(newPointsNum);
        for(size_t i = 0; i < newPointsNum; i++)
    {
            size_t matchingPointIdx = mask[i];
            double radius;

            if(matchingPointIdx >= originalPointsNum)
            {
                // the point is a new, but we take its initial radius to be the same as the closest point's radius
                size_t closestPointIdx = oldPointsTree.closestPoint(points[i], false /* Dont include self */).getIndex();
                radius = this->radiuses.at(closestPointIdx);
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

std::vector<Point3D> MadVoro::Voronoi3D::Voronoi3DImpl::BuildPartiallyParallel(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(const size_t &idx : indicesToBuild)
    {
        if(idx >= allPoints.size())
        {
            MadVoro::Exception::InvalidArgumentException eo("BuildParatiallyParallel: Illegal point was given");
            eo.addEntry("Index", idx);
            eo.addEntry("Size", allPoints.size());
            throw eo;
        }
    }

    std::vector<Point3D> activePoints = this->PrepareToBuildParallel(allPoints, allWeights, indicesToBuild, suppressRebalancing);
    
    std::vector<size_t> order;

    // build delaunay
    if(not activePoints.empty())
    {
        std::pair<Point3D, Point3D> bounding_box = std::make_pair(activePoints[0], activePoints[0]);
        for(const Point3D &point : activePoints)
        {
            bounding_box.first.x = std::min(bounding_box.first.x, point.x);
            bounding_box.second.x = std::max(bounding_box.second.x, point.x);
            bounding_box.first.y = std::min(bounding_box.first.y, point.y);
            bounding_box.second.y = std::max(bounding_box.second.y, point.y);
            bounding_box.first.z = std::min(bounding_box.first.z, point.z);
            bounding_box.second.z = std::max(bounding_box.second.z, point.z);
        }

        // performs internal tesselation:
        // std::cout << "checking duplications..." << std::endl;
        // reportDuplications(new_points);
        order = HilbertOrder3D(activePoints);
        
        // initial build for the points
        this->del_.Build(activePoints, bounding_box.second, bounding_box.first, order);
    }

    // updates the radiuses array of the tetrahedra, as well as the lists for each point what tetras it belongs to
    this->R_.resize(this->del_.tetras_.size());
    std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
    this->tetra_centers_.resize(this->R_.size());
    this->bigtet_ = SetPointTetras(this->PointTetras_, this->Norg_, this->del_.tetras_, this->del_.empty_tetras_);

    if(this->radiuses.size() < this->Norg_)
    {
        MadVoro::Exception::InvalidArgumentException eo("Voronoi3D:PrepareToBuildParallel: wrong size of radiuses array");
        eo.addEntry("Rank", rank);
        eo.addEntry("radiuses.size()", this->radiuses.size());
        eo.addEntry("Norg_", this->Norg_);
        throw eo;
    }

    this->allMyPointsTree = std::make_shared<IndexedPointsTree>(Range::IndexedPoint3D(this->ll_, std::numeric_limits<size_t>::max()),
                                                                Range::IndexedPoint3D(this->ur_, std::numeric_limits<size_t>::max()));
    size_t allPointsNum = this->allMyPoints.size();
    for(size_t pointIdx = 0; pointIdx < allPointsNum; pointIdx++)
    {
        const Point3D &point = activePoints[pointIdx];
        this->allMyPointsTree->insert(Range::IndexedPoint3D(point.x, point.y, point.z, pointIdx));
    }

    if(this->allMyPoints.size() == activePoints.size())
    {
        // not a real parital build
        this->myPointsTree = this->allMyPointsTree;
    }
    else
    {
        this->myPointsTree = std::make_shared<IndexedPointsTree>(Range::IndexedPoint3D(this->ll_, std::numeric_limits<size_t>::max()),
                                                                Range::IndexedPoint3D(this->ur_, std::numeric_limits<size_t>::max()));
        for(size_t pointIdx = 0; pointIdx < this->Norg_; pointIdx++)
        {
            const Point3D &point = activePoints[pointIdx];
            this->myPointsTree->insert(Range::IndexedPoint3D(point.x, point.y, point.z, pointIdx));
        }
    }

    this->UpdateRadiuses(activePoints);    

    this->UpdateRangeFinder();

    this->BringGhostPointsToBuild(MPI_COMM_WORLD);

    CM_.resize(del_.points_.size());
    volume_.resize(Norg_);

    if(not activePoints.empty())
    {
        // Create Voronoi
        BuildVoronoi(order);
    }

    // todo: why?
    // std::vector<double>().swap(this->R_);
    // std::vector<tetra_vec>().swap(this->PointTetras_);

    this->UpdateCMs();

    // save the list of the real ghost points
    // this->FilterRealGhostPoints();

    return this->allMyPoints;
}
#endif // MADVORO_WITH_MPI

/**
 * \author Maor Mizrachi
 * \brief Gets a point, its radius, a box and the normals to the box's faces, and returns the faces indices that the sphere (around `point`, in the given `radius`) intersects
*/
std::vector<size_t> CheckToMirror(const MadVoro::Geometry::Sphere<Point3D> &sphere, const std::vector<Face3D> &box, const std::vector<Point3D> &normals)
{
    std::vector<size_t> facesItCuts;
    // std::cout << "point = " << point << ", radius = " << radius << std::endl;
    for(size_t i = 0; i < box.size(); i++)
    {
        // check for intersecting the sphere with radius `radius` around `point`, with the `i`th face of `box`
        if(Face3DSphereIntersections(box[i], sphere, normals[i]))
        {
            // intersects! mirror the point
            facesItCuts.push_back(i);
        }
    }
    return facesItCuts;
}

void MadVoro::Voronoi3D::Voronoi3DImpl::UpdateCMs(void)
{
    // first, calculate CM for active local points
    this->CalcAllCM(); // Now this->CM_ calculates correct CM for all active points, and maybe for more
    size_t Face3DNeighborsSize = FaceNeighbors_.size();
    for(std::size_t i = 0; i < Face3DNeighborsSize; ++i)
    {
        if(this->BoundaryFace3D(i))
        {
            this->CalcRigidCM(i);
        }
    }

    this->SyncPartialBuildData(this->CM_, this->all_CM);
}

void MadVoro::Voronoi3D::Voronoi3DImpl::UpdateRadiuses(const std::vector<Point3D> &points)
{
    // use an oct tree to fast calculate the distance to closest point
    DataStructure::OctTree<Point3D> myOctTree(this->ll_, this->ur_, this->allMyPoints.begin(), this->allMyPoints.end());
    for(const std::pair<size_t, size_t> &indices : this->indicesInAllMyPoints)
    {
        const size_t &pointIndexOnBuild = indices.first;
        const Point3D &point = points[pointIndexOnBuild];
        size_t pointIndexAmongAll = indices.second; 
        if(this->radiuses[pointIndexAmongAll] <= 0)
        {
            // point does not have a radius from a previous timestep. Initialize a radius
            this->radiuses[pointIndexAmongAll] = Utils::fastsqrt(this->allMyPointsTree->closestPointDistance(point, false)); // todo second closest
        }
    }
}

void MadVoro::Voronoi3D::Voronoi3DImpl::UpdateRangeFinder()
{
    this->rangeFinder = std::make_shared<Range::OctTreeFinder>(this->allMyPointsTree.get(), this->allMyPoints);
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

    //     std::vector<Point3D> pointsToAdd;
    //     // todo
    //     this->rangeFinder->replacePoints();
    // }
}

/**
 * \author Maor Mizrachi
 * \brief Creates a batch for a cycle (iteration) in the ghost points bringing loop
*/
std::pair<std::vector<SmallRangeQueryData>, std::vector<BigRangeQueryData>> MadVoro::Voronoi3D::Voronoi3DImpl::CreateBatches(boost::container::flat_set<size_t> &smallPoints, boost::container::flat_set<size_t> &largePoints, const boost::container::flat_map<size_t, size_t> &firstLargeIteration, std::vector<double> &currentRadiuses, size_t iterations)
{
    std::vector<SmallRangeQueryData> smallQueries;
    std::vector<BigRangeQueryData> bigQueries;
    boost::container::flat_set<size_t> tetraToCancel;

    if(iterations == 1)
    {        
        // at first iteration, run an initial query, all the points are small
        for(const size_t &pointIdx : smallPoints)
        {
            const Point3D &point = this->del_.points_[pointIdx];
            smallQueries.emplace_back();
            SmallRangeQueryData &query = smallQueries.back();
            query.pointIdx = pointIdx;
            query.center = {point.x, point.y, point.z};
            query.radius = currentRadiuses[pointIdx];
            query.maxPointsToGet = RANGE_MAX_POINTS_TO_GET + 1;

            if(currentRadiuses[pointIdx] <= 0)
            {
                MadVoro::Exception::MadVoroException eo("Radius for a certain point is <= 0 (in 'MadVoro::Voronoi3D::CreateBatches')");
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
            const Point3D &point = this->del_.points_[pointIdx];

            // std::cout << "firstLargeIteration.at(pointIdx) is " << firstLargeIteration.at(pointIdx) << std::endl;
            bool askOnlyClose = (iterations == firstLargeIteration.at(pointIdx)); // on the first iteration as large, ask only the close ranks

            for(const size_t &tetraIdx : this->PointTetras_[pointIdx])
            {
                if(not this->del_.tetras_[tetraIdx].newTetra)
                {
                    continue; // tetra does not need to be checked
                }
                const Point3D &center = this->tetra_centers_[tetraIdx];
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
            const Point3D &point = this->del_.points_[pointIdx];
            double radius = currentRadiuses[pointIdx];
           
            if(currentRadiuses[pointIdx] <= 0)
            {
                MadVoro::Exception::MadVoroException eo("Radius for a certain point is <= 0 (in 'MadVoro::Voronoi3D::CreateBatches')");
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
        this->del_.tetras_[tetraIdx].newTetra = false;
    }
    return {smallQueries, bigQueries};
}

/**
 * \author Maor Mizrachi
 * \brief Gets a list of query, and tests for creating mirror points. In the end of this procedure, `mirroredPoints` contains pairs of <faceIdx, pointIdx>, of points that should be mirrored, in relative to which faces
*/
template<typename QueryDataType>
std::vector<std::pair<size_t, size_t>> MirrorPoints(const std::vector<QueryDataType> &queries, const std::vector<Face3D> &box, const std::vector<Point3D> &normals)
{   
    static_assert(std::is_convertible<QueryDataType*, RangeQueryData*>::value, "MirrorPoints: QueryDataType must inherit 'RangeQueryData'");

    std::vector<std::pair<size_t, size_t>> mirroredPoints;
    for(const QueryDataType &query : queries)
    {
        // check for mirroring:
        MadVoro::Geometry::Sphere<Point3D> sphere(Point3D(query.center), query.radius);
        size_t pointIdx = query.pointIdx;

        std::vector<size_t> facesItCuts = CheckToMirror(sphere, box, normals);

        for(const size_t &faceIdx : facesItCuts)
        {
            mirroredPoints.push_back(std::make_pair(faceIdx, pointIdx));
        }
    }
    return mirroredPoints;
}

void MadVoro::Voronoi3D::Voronoi3DImpl::BringSelfGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                                            BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                                            boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                                                            boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints,
                                                            boost::container::flat_set<size_t> &selfIgnorePoints)
{
    std::vector<Point3D> newPoints;
    {
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
    }
    {
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
    }
    this->del_.BuildExtra(newPoints);
}

#ifdef MADVORO_WITH_MPI
    void MadVoro::Voronoi3D::Voronoi3DImpl::BringRemoteGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                                                    BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                                                    boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                                                                    boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints)
    {
        std::vector<Point3D> newPoints;
        // large points queries
        {
            QueryAgent::QueryBatchInfo<BigRangeQueryData, _3DPoint> bigBatchInfo = bigRangeAgent.runBatch(bigQueries);
            newPoints.reserve(bigBatchInfo.result.size());
            for(const QueryAgent::QueryInfo<BigRangeQueryData, _3DPoint> &ans : bigBatchInfo.queriesAnswers)
            {
                numOfResultsForBigPoints[ans.data.pointIdx] += ans.finalResults.size();
            }
            for(const _3DPoint &_point : bigBatchInfo.result)
            {
                newPoints.push_back(Point3D(_point.x, _point.y, _point.z));
            }
            this->SetGhostArray(bigRangeAgent.getRecvProc(), bigRangeAgent.getRecvPoints());
            this->del_.BuildExtra(newPoints);
        }
        // small points queries
        {        
            QueryAgent::QueryBatchInfo<SmallRangeQueryData, _3DPoint> smallBatchInfo = smallRangeAgent.runBatch(smallQueries);
            for(const QueryAgent::QueryInfo<SmallRangeQueryData, _3DPoint> &ans : smallBatchInfo.queriesAnswers)
            {
                numOfResultsForSmallPoints[ans.data.pointIdx] += ans.finalResults.size();
            }
            newPoints.clear();
            newPoints.reserve(smallBatchInfo.result.size());
            for(const _3DPoint &_point : smallBatchInfo.result)
            {
                newPoints.push_back(Point3D(_point.x, _point.y, _point.z));
            }
            this->SetGhostArray(smallRangeAgent.getRecvProc(), smallRangeAgent.getRecvPoints());
            this->del_.BuildExtra(newPoints);
        }    
    }
#endif // MADVORO_WITH_MPI

/**
 * \author Maor Mizrachi
 * \brief Calculates the points for next iteration, and determines each one's type (small or big)
*/
std::pair<boost::container::flat_set<size_t>, boost::container::flat_set<size_t>>
MadVoro::Voronoi3D::Voronoi3DImpl::DetermineNextIterationPoints(size_t iterations,
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
#ifdef MADVORO_WITH_MPI
    void MadVoro::Voronoi3D::Voronoi3DImpl::BringGhostPointsToBuild(const MPI_Comm &comm)
#else // MADVORO_WITH_MPI
    void MadVoro::Voronoi3D::Voronoi3DImpl::BringGhostPointsToBuild()
#endif // MADVORO_WITH_MPI
{
    int rank = 0, size = 1;
    #ifdef MADVORO_WITH_MPI
        
        int mpiInitialized;
        MPI_Initialized(&mpiInitialized);

        if(mpiInitialized)
        {
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &size);
        }
    #endif // MADVORO_WITH_MPI

    std::vector<Face3D> box;
    std::vector<Point3D> normals;
    this->InitialBoxBuild(box, normals);
    
    boost::container::flat_set<size_t> smallPoints; // indices of 'small' points
    boost::container::flat_set<size_t> largePoints; // indices of 'large' points
    boost::container::flat_map<size_t, size_t> firstLargeIteration;

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

    #ifdef MADVORO_WITH_MPI
        auto [ghostPointsFromLastBuild, alreadySentProc, alreadySentPoints, alreadyRecvProcs, alreadyRecvPoints] = this->InitialGhostPointsExchange();
        this->SetGhostArray(alreadyRecvProcs, alreadyRecvPoints);
        this->del_.BuildExtra(ghostPointsFromLastBuild);
        this->R_.resize(this->del_.tetras_.size(), RADIUS_UNINITIALIZED);
        std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
        this->tetra_centers_.resize(this->R_.size());
        this->bigtet_ = SetPointTetras(this->PointTetras_, this->Norg_, this->del_.tetras_, this->del_.empty_tetras_);
        
        SentPointsContainer pointsContainer(alreadySentProc, alreadySentPoints);
    #endif // MADVORO_WITH_MPI
    
    #ifdef MADVORO_WITH_MPI
        BigRangeAgent bigRangeAgent(this->rangeFinder.get(), this->pointsManager->getEnvironmentAgent(), pointsContainer, comm);
        SmallRangeAgent smallRangeAgent(this->rangeFinder.get(), this->pointsManager->getEnvironmentAgent(), pointsContainer, comm);
    #else // MADVORO_WITH_MPI
        BigRangeAgent bigRangeAgent(this->rangeFinder.get());
        SmallRangeAgent smallRangeAgent(this->rangeFinder.get());
    #endif // MADVORO_WITH_MPI

    #ifdef MADVORO_WITH_MPI
        MPI_Request finishedReq;
        int I_finished = 0;
    #endif // MADVORO_WITH_MPI
    int finished;

    std::vector<std::pair<size_t, size_t>> allMirrored;

    size_t iterations = 0;

    bool considerOwnPoints = (size == 1) or (this->Norg_ != this->allMyPoints.size()); // there are points which we ignore in this step, so we have, in the range searching, communicate with us as well
    boost::container::flat_set<size_t> selfIgnorePoints;
    for(const std::pair<size_t, size_t> &indices : this->indicesInAllMyPoints)
    {
        selfIgnorePoints.insert(indices.second);
    }

    while(true) // loop is not really infinite (has 'break')
    {
        boost::container::flat_map<size_t, size_t> numOfResultsForSmallPoints;
        boost::container::flat_map<size_t, size_t> numOfResultsForBigPoints;
        
        size_t smallPointsNum = smallPoints.size();
        size_t largePointsNum = largePoints.size();
        #ifdef MADVORO_WITH_MPI
            if(mpiInitialized)
            {
                MPI_Reduce((rank == 0)? MPI_IN_PLACE : &smallPointsNum, &smallPointsNum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
                MPI_Reduce((rank == 0)? MPI_IN_PLACE : &largePointsNum, &largePointsNum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
            }
        #endif // MADVORO_WITH_MPI
        iterations++;
        if(rank == 0 and this->verbosity) std::cout << "iteration " << iterations << " (" << smallPointsNum << " small points, " << largePointsNum << " large points)" << std::endl;
        auto [smallQueries, bigQueries] = this->CreateBatches(smallPoints, largePoints, firstLargeIteration, currentRadiuses, iterations);
        std::vector<std::pair<size_t, size_t>> mirroredPoints = MirrorPoints(smallQueries, box, normals);
        std::vector<std::pair<size_t, size_t>> moreMirroredPoints = MirrorPoints(bigQueries, box, normals);
        mirroredPoints.insert(mirroredPoints.end(), moreMirroredPoints.begin(), moreMirroredPoints.end());

        #ifdef MADVORO_WITH_MPI
            I_finished = (smallQueries.empty() and bigQueries.empty())? 1 : 0;
            if(mpiInitialized)
            {
                MPI_Iallreduce(&I_finished, &finished, 1, MPI_INT, MPI_SUM, comm, &finishedReq);
            }
        #else // MADVORO_WITH_MPI
            finished = (smallQueries.empty() and bigQueries.empty())? 1 : 0;
        #endif // MADVORO_WITH_MPI

        if(considerOwnPoints)
        {
            this->BringSelfGhostPoints(bigQueries, smallQueries, bigRangeAgent, smallRangeAgent, numOfResultsForBigPoints, numOfResultsForSmallPoints, selfIgnorePoints);
        }

        #ifdef MADVORO_WITH_MPI
            this->BringRemoteGhostPoints(bigQueries, smallQueries, bigRangeAgent, smallRangeAgent, numOfResultsForBigPoints, numOfResultsForSmallPoints);
        #endif // MADVORO_WITH_MPI
        
        std::vector<Point3D> newPoints;
        newPoints.reserve(mirroredPoints.size());
        allMirrored.reserve(allMirrored.size() + mirroredPoints.size());
        for(const std::pair<size_t, size_t> &pairFace3DPoint : mirroredPoints)
        {
            // check if we have already mirrored this point with this face
            if(std::find(allMirrored.begin(), allMirrored.end(), pairFace3DPoint) == allMirrored.end())
            {
                allMirrored.push_back(pairFace3DPoint); // remember we mirrored this point with this face
                newPoints.push_back(MirrorPoint(box[pairFace3DPoint.first], this->del_.points_[pairFace3DPoint.second]));
            }
        }
        this->del_.BuildExtra(newPoints);

        this->R_.resize(this->del_.tetras_.size(), RADIUS_UNINITIALIZED);
        std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
        this->tetra_centers_.resize(this->R_.size());
        this->bigtet_ = SetPointTetras(this->PointTetras_, this->Norg_, this->del_.tetras_, this->del_.empty_tetras_);
        
        std::tie(smallPoints, largePoints) = this->DetermineNextIterationPoints(iterations, firstLargeIteration, currentRadiuses, numOfResultsForSmallPoints, numOfResultsForBigPoints);

        // #ifdef MADVORO_WITH_MPI
        //     std::tie(smallPoints, largePoints) = this->DetermineNextIterationPoints(iterations, firstLargeIteration, currentRadiuses, selfSmallQueriesAnswers, selfBigQueriesAnswers, smallBatchInfo.queriesAnswers, bigBatchInfo.queriesAnswers);
        // #else // MADVORO_WITH_MPI
        // #endif // MADVORO_WITH_MPI

        #ifdef MADVORO_WITH_MPI
            if(mpiInitialized)
            {
                MPI_Wait(&finishedReq, MPI_STATUS_IGNORE);
            }
        #endif // MADVORO_WITH_MPI

        if(finished == size)
        {
            break;
        }
    }

    #ifdef MADVORO_WITH_MPI        
        const std::vector<std::vector<size_t>> &sentPoints = pointsContainer.getSentData();
        const std::vector<int> &sentProc = pointsContainer.getSentProc();

        // calculate this->duplicated_points_
        this->UpdateDuplicatedPoints(sentProc, sentPoints);
        // remove whomever that does not appear both in my sent vector and receive vector (because if one appears in only one, it means that we either sent it a point, or received one, but has no used of it at all (otherwise it would require a symetric call))
        // this->EnsureSymmetry(sentProc, {alreadyRecvProcs, smallRangeAgent.getRecvProc(), bigRangeAgent.getRecvProc()});    // todo: uncomment
        this->EnsureSymmetry(sentProc, {alreadyRecvProcs, smallRangeAgent.getRecvProc(), bigRangeAgent.getRecvProc()});    
    #endif // MADVORO_WITH_MPI
}

#ifdef MADVORO_WITH_MPI
    vector<vector<std::size_t>> &MadVoro::Voronoi3D::Voronoi3DImpl::GetGhostIndeces(void)
    {
        return Nghost_;
    }
#endif // MADVORO_WITH_MPI

void MadVoro::Voronoi3D::Voronoi3DImpl::CalcAllCM(void)
{
    std::array<Point3D, 4> tetra;
    size_t Nfaces = FaceNeighbors_.size();
    assert(Nfaces == 0 or Nfaces >= 4);
    Point3D vtemp;
    std::vector<Point3D> vectemp;
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
            CM_[i] = Point3D();
            volume_[i] = 0;
            Nfaces = FacesInCell_[i].size();
            for (size_t k = 0; k < Nfaces; ++k)
            {
                size_t Face3D = FacesInCell_[i][k];
                size_t Npoints = PointsInFace_[Face3D].size();
                tetra[0] = tetra_centers_[PointsInFace_[Face3D][0]];
                for (std::size_t j = 0; j < Npoints - 2; ++j)
                {
                    tetra[1] = tetra_centers_[PointsInFace_[Face3D][j + 1]];
                    tetra[2] = tetra_centers_[PointsInFace_[Face3D][j + 2]];
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

std::pair<Point3D, Point3D> MadVoro::Voronoi3D::Voronoi3DImpl::GetBoxCoordinates(void) const
{
    return std::pair<Point3D, Point3D>(ll_, ur_);
}

void MadVoro::Voronoi3D::Voronoi3DImpl::BuildNoBox(vector<Point3D> const &points, vector<vector<Point3D>> const &ghosts, vector<size_t> toduplicate)
{
    assert(points.size() > 0);
    // Clear data
    PointTetras_.clear();
    R_.clear();
    R_.reserve(points.size());
    tetra_centers_.clear();
    tetra_centers_.reserve(points.size() * 7);
    del_.Clean();
    // Voronoi Data
    FacesInCell_.clear();
    PointsInFace_.clear();
    FaceNeighbors_.clear();
    CM_.clear();
    Face_CM_.clear();
    volume_.clear();
    area_.clear();
    Norg_ = points.size();
    #ifdef MADVORO_WITH_MPI
        duplicatedprocs_.clear();
        duplicated_points_.clear();
        Nghost_.clear();
    #endif // MADVORO_WITH_MPI

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
        vector<Point3D> extra_points = CreateBoundaryPoints(duplicate, past_duplicates);
        del_.BuildExtra(extra_points);
    }

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), RADIUS_UNINITIALIZED);
    tetra_centers_.resize(R_.size());
    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

    CM_.resize(Norg_);
    volume_.resize(Norg_, 0);
    // Create Voronoi
    BuildVoronoi(order);

    CalcAllCM();
    CM_.resize(del_.points_.size());
    for (std::size_t i = 0; i < FaceNeighbors_.size(); ++i)
        if (BoundaryFace3D(i))
            CalcRigidCM(i);
}

void MadVoro::Voronoi3D::Voronoi3DImpl::BuildDebug(int rank)
{
    std::vector<size_t> order = IO::read_vecst("order_" + std::to_string(rank) + ".bin");
    std::vector<Point3D> points = IO::read_vec3d("points0_" + std::to_string(rank) + ".bin");
    Norg_ = points.size();
    std::vector<Point3D> bb = IO::read_vec3d("bb_" + std::to_string(rank) + ".bin");
    del_.Build(points, bb[1], bb[0], order);
    points = IO::read_vec3d("points1_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);
    points = IO::read_vec3d("points2_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);
    points = IO::read_vec3d("points3_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);
    points = IO::read_vec3d("points4_" + std::to_string(rank) + ".bin");
    del_.BuildExtra(points);

    bigtet_ = SetPointTetras(PointTetras_, Norg_, del_.tetras_, del_.empty_tetras_);

    R_.resize(del_.tetras_.size());
    std::fill(R_.begin(), R_.end(), RADIUS_UNINITIALIZED);
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
        if (BoundaryFace3D(i))
            CalcRigidCM(i);
}

void MadVoro::Voronoi3D::Voronoi3DImpl::BuildPartially(const std::vector<Point3D> &allPoints, const std::vector<size_t> &indicesToBuild)
{
    #ifdef MADVORO_WITH_MPI
        int mpiInitialized;
        MPI_Initialized(&mpiInitialized);

        if(not mpiInitialized)
        {
            std::cerr << "MadVoro is compiled with MPI. You MUST use `MPI_Init()` at the beginning of your program." << std::endl;
            std::cerr << "You can use `BuildParallel` for parallel build, or continue with sequential build (this function)" << std::endl;
            throw Exception::MadVoroException("MPI is not initialized");
        }
    #endif // MADVORO_WITH_MPI

    if(this->radiuses.size() < allPoints.size())
    {
        this->radiuses.resize(allPoints.size(), RADIUS_UNINITIALIZED);
    }

    if(this->all_CM.size() < allPoints.size())
    {
        this->all_CM.resize(allPoints.size());
    }

    #ifdef MADVORO_WITH_MPI
        if(this->pointsManager.get() == nullptr)
        {
            // initialize points manager
            this->pointsManager = std::shared_ptr<HilbertPointsManager>(new HilbertPointsManager(this->ll_, this->ur_, MPI_COMM_SELF));
        }

        this->pointsManager->update(allPoints, std::vector<double>(allPoints.size(), 1), indicesToBuild, this->radiuses, this->all_CM); // does rebalancing (if necessary) and exchanging    
    #endif // MADVORO_WITH_MPI
    
    this->allMyPoints = allPoints;
    // this->radiuses.resize(allPoints.size(), RADIUS_UNINITIALIZED);
    std::vector<Point3D> activePoints;
    activePoints.reserve(indicesToBuild.size());
    for(const size_t &pointIdx : indicesToBuild)
    {
        activePoints.push_back(allPoints[pointIdx]);
    }

    size_t pointsCounter = 0;
    for(const size_t &pointIdx : indicesToBuild)
    {
        this->indicesInAllMyPoints[pointIdx] = pointsCounter;
        pointsCounter++;
    }

    this->BuildInitialize(activePoints.size());

    std::vector<size_t> order = HilbertOrder3D(activePoints);

    // build delaunay
    if(not activePoints.empty())
    {
        std::pair<Point3D, Point3D> bounding_box = std::make_pair(activePoints[0], activePoints[0]);
        for(const Point3D &point : activePoints)
        {
            bounding_box.first.x = std::min(bounding_box.first.x, point.x);
            bounding_box.second.x = std::max(bounding_box.second.x, point.x);
            bounding_box.first.y = std::min(bounding_box.first.y, point.y);
            bounding_box.second.y = std::max(bounding_box.second.y, point.y);
            bounding_box.first.z = std::min(bounding_box.first.z, point.z);
            bounding_box.second.z = std::max(bounding_box.second.z, point.z);
        }
        if(activePoints.size() == 1)
        {
            bounding_box.first = this->ll_;
            bounding_box.second = this->ur_;
        }
        // performs internal tesselation:
        order = HilbertOrder3D(activePoints);
        
        // initial build for the points
        this->del_.Build(activePoints, bounding_box.second, bounding_box.first, order);
    }

    // updates the radiuses array of the tetrahedra, as well as the lists for each point what tetras it belongs to
    this->R_.resize(this->del_.tetras_.size());
    std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
    this->tetra_centers_.resize(this->R_.size());
    this->bigtet_ = SetPointTetras(this->PointTetras_, this->Norg_, this->del_.tetras_, this->del_.empty_tetras_);

    this->allMyPointsTree = std::make_shared<IndexedPointsTree>(Range::IndexedPoint3D(this->ll_, std::numeric_limits<size_t>::max()),
                                                                Range::IndexedPoint3D(this->ur_, std::numeric_limits<size_t>::max()));
    size_t allPointsNum = this->allMyPoints.size();
    for(size_t pointIdx = 0; pointIdx < allPointsNum; pointIdx++)
    {
        const Point3D &point = activePoints[pointIdx];
        this->allMyPointsTree->insert(Range::IndexedPoint3D(point.x, point.y, point.z, pointIdx));
    }

    if(this->allMyPoints.size() == activePoints.size())
    {
        // not a real parital build
        this->myPointsTree = this->allMyPointsTree;
    }
    else
    {
        this->myPointsTree = std::make_shared<IndexedPointsTree>(Range::IndexedPoint3D(this->ll_, std::numeric_limits<size_t>::max()),
                                                                Range::IndexedPoint3D(this->ur_, std::numeric_limits<size_t>::max()));
        for(size_t pointIdx = 0; pointIdx < this->Norg_; pointIdx++)
        {
            const Point3D &point = activePoints[pointIdx];
            this->myPointsTree->insert(Range::IndexedPoint3D(point.x, point.y, point.z, pointIdx));
        }
    }

    this->UpdateRadiuses(activePoints);

    this->UpdateRangeFinder();

    #ifdef MADVORO_WITH_MPI
        this->BringGhostPointsToBuild(MPI_COMM_SELF);
    #else // MADVORO_WITH_MPI
        this->BringGhostPointsToBuild(); 
    #endif // MADVORO_WITH_MPI

    CM_.resize(del_.points_.size());
    volume_.resize(Norg_);

    // Create Voronoi
    BuildVoronoi(order);

    this->UpdateCMs();
}

void MadVoro::Voronoi3D::Voronoi3DImpl::Build(const std::vector<Point3D> &points)
{
    std::vector<size_t> indicesToBuild(points.size());
    std::iota(indicesToBuild.begin(), indicesToBuild.end(), 0);
    this->BuildPartially(points, indicesToBuild);
}

void MadVoro::Voronoi3D::Voronoi3DImpl::BuildVoronoi(std::vector<size_t> const &order)
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

    size_t Face3DCounter = 0;
    boost::container::flat_set<size_t> neigh_set;
    point_vec *temp_points_in_face;
    std::array<Point3D, 128> clean_vec;
    std::array<double, 128> area_vec_temp;

    //std::vector<Point3D, boost::alignment::aligned_allocator<Point3D, 32> > clean_vec;
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
                        temp_points_in_face = &PointsInFace_[Face3DCounter];
                        //temp_points_in_face->reserve(8);
                        double Asize = CleanDuplicates(temp, tetra_centers_, *temp_points_in_face, ScalarProd(del_.points_[point] - del_.points_[point_other], del_.points_[point] - del_.points_[point_other]), diffs, clean_vec, temp_size);
                        if (temp_points_in_face->size() < 3)
                            continue;
                        CalcFaceAreaCM(*temp_points_in_face, tetra_centers_, clean_vec, area_[Face3DCounter],
                                                     Face_CM_[Face3DCounter], Atempvec);
                        if (area_[Face3DCounter] < (Asize * (IsPointOutsideBox(point_other) ? 1e-14 : 1e-15)))
                            continue;
                        if (point_other >= Norg_ && point_other < (Norg_ + 4))
                        {
                            MadVoro::Exception::MadVoroException eo("Neighboring big tet point");
                            throw eo;
                        }
                        // Make faces right handed
                        MakeRightHandFace(*temp_points_in_face, del_.points_[point], tetra_centers_, temp3, area_[Face3DCounter]);
                        try
                        {
                            CleanSameLine(*temp_points_in_face, tetra_centers_, area_vec_temp);
                        }
                        catch(Exception::MadVoroException &eo)
                        {
                            eo.addEntry("Points0", del_.points_[point]);
                            eo.addEntry("Points1", del_.points_[point_other]);
                            eo.addEntry("diff", abs(del_.points_[point] - del_.points_[point_other]));
                            eo.addEntry("Asize", Asize);
                            int rank = 0;
                            #ifdef MADVORO_WITH_MPI
                               MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                            #endif // MADVORO_WITH_MPI
                            // throw eo;
                            continue;
                        }
                        FaceNeighbors_[Face3DCounter].first = point;
                        FaceNeighbors_[Face3DCounter].second = point_other;

                        FacesInCell_[point].push_back(Face3DCounter);
                        if (point_other < Norg_)
                        {
                            FacesInCell_[point_other].push_back(Face3DCounter);
                        }
                        neigh_set.insert(point_other);
                        ++Face3DCounter;
                        // realloc memory if needed
                        if (Face3DCounter == FaceNeighbors_.size())
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

    // Fix Face3D CM (this prevents large face velocities for close by points)
    size_t Nfaces = Face3DCounter;
    Point3D mid, norm;
    for (size_t i = 0; i < Nfaces; ++i)
    {
        mid = del_.points_[FaceNeighbors_[i].first];
        mid += del_.points_[FaceNeighbors_[i].second];
        mid *= 0.5;
        norm = del_.points_[FaceNeighbors_[i].second];
        norm -= del_.points_[FaceNeighbors_[i].first];
        Face_CM_[i] -= ScalarProd(Face_CM_[i] - mid, norm) * norm / ScalarProd(norm, norm);
    }

    area_.resize(Face3DCounter);
    area_.shrink_to_fit();
    Face_CM_.resize(Face3DCounter);
    Face_CM_.shrink_to_fit();
    FaceNeighbors_.resize(Face3DCounter);
    FaceNeighbors_.shrink_to_fit();
    PointsInFace_.resize(Face3DCounter);
    PointsInFace_.shrink_to_fit();
    for (size_t i = 0; i < Norg_; ++i)
        FacesInCell_[i].shrink_to_fit();
}

inline double MadVoro::Voronoi3D::Voronoi3DImpl::GetRadius(const size_t &index) const
{ 
    R_[index] = (R_[index] < 0)? CalcTetraRadiusCenter(index) : R_[index];
    if(std::isnan(this->R_[index]) or not std::isfinite(this->R_[index]))
    {
        MadVoro::Exception::MadVoroException eo("Voronoi3D:GetRadius: Radius is invalid");
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
        std::vector<Point3D> tetPoints = std::vector({this->del_.points_[tet.points[0]], this->del_.points_[tet.points[1]], this->del_.points_[tet.points[2]], this->del_.points_[tet.points[3]]});
        eo.addEntry("Tetra Points", tetPoints);
        eo.addEntry("Norg", this->Norg_);
        throw eo;
    }
    return this->R_[index];
}

void MadVoro::Voronoi3D::Voronoi3DImpl::FindIntersectionsSingle(vector<Face3D> const &box, std::size_t point, MadVoro::Geometry::Sphere<Point3D> &sphere,
                                                                                vector<size_t> &intersecting_faces, std::vector<double> &Rtemp, std::vector<Point3D> &vtemp)
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
        Point3D normal = CrossProduct(box[j].vertices[1] - box[j].vertices[0], box[j].vertices[2] - box[j].vertices[0]);
        normal *= (1.0 / Utils::fastsqrt(ScalarProd(normal, normal)));
        for (std::size_t i = 0; i < N; ++i)
        {
            sphere.radius = Rtemp[i];
            sphere.center = vtemp[i];
            if (Face3DSphereIntersections(box[j], sphere, normal))
            {
                intersecting_faces.push_back(j);
                break;
            }
        }
    }
}

void MadVoro::Voronoi3D::Voronoi3DImpl::GetPointToCheck(std::size_t point, vector<unsigned char> const &checked, vector<std::size_t> &res)
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
    res = Utils::unique(res);
}

std::size_t MadVoro::Voronoi3D::Voronoi3DImpl::GetFirstPointToCheck(void) const
{
    std::size_t i;
    Tetrahedron const &tet = del_.tetras_[bigtet_];
    for (i = 0; i < 4; ++i)
        if (tet.points[i] < Norg_)
            break;
    if (i < 4)
        return tet.points[i];
    else
        throw MadVoro::Exception::MadVoroException("Can't find first point to start boundary search");
}

vector<std::pair<std::size_t, std::size_t>> MadVoro::Voronoi3D::Voronoi3DImpl::SerialFirstIntersections(void)
{
    vector<Face3D> box;
    vector<Point3D> normals;
    this->InitialBoxBuild(box, normals);
    size_t Nfaces = box.size();

    //    vector<std::size_t> point_neigh;
    vector<std::pair<std::size_t, std::size_t>> res;
    MadVoro::Geometry::Sphere<Point3D> sphere;
    vector<unsigned char> will_check(Norg_, 0);
    std::size_t cur_loc;
    std::stack<std::size_t> check_stack;
    FirstCheckList(check_stack, will_check, Norg_, del_, PointTetras_);
    std::vector<double> vdist(Nfaces);
    std::vector<Point3D> vtemp(Nfaces);
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

vector<std::pair<std::size_t, std::size_t>> MadVoro::Voronoi3D::Voronoi3DImpl::SerialFindIntersections(bool first_run)
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
    vector<Face3D> box = box_faces_.empty() ? BuildBox(ll_, ur_) : box_faces_;
    vector<std::size_t> point_neigh;
    vector<std::pair<std::size_t, std::size_t>> res;
    MadVoro::Geometry::Sphere<Point3D> sphere;
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
    std::vector<Point3D> vtemp;
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

double MadVoro::Voronoi3D::Voronoi3DImpl::CalcTetraRadiusCenter(const size_t &index) const
{
    Point3D v2(del_.points_[del_.tetras_[index].points[1]]);
    v2 -= del_.points_[del_.tetras_[index].points[0]];
    Point3D v3(del_.points_[del_.tetras_[index].points[2]]);
    v3 -= del_.points_[del_.tetras_[index].points[0]];
    Point3D v4(del_.points_[del_.tetras_[index].points[3]]);
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
    Point3D center = Point3D(DDx / (2 * a), DDy / (2 * a), DDz / (2 * a)) + del_.points_[del_.tetras_[index].points[0]];
    tetra_centers_[index] = center;
    double Rres = 0.5 * std::sqrt(DDx * DDx + DDy * DDy + DDz * DDz) / std::abs(a);
    // Sanity check
    /*double Rcheck0 = fastabs(del_.points_[del_.tetras_[index].points[0]] - center);
        double Rcheck1 = fastabs(del_.points_[del_.tetras_[index].points[1]] - center);
        double Rcheck2 = fastabs(del_.points_[del_.tetras_[index].points[2]] - center);
        double Rcheck3 = fastabs(del_.points_[del_.tetras_[index].points[3]] - center);*/
    Point3D v1(del_.points_[del_.tetras_[index].points[0]]);
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

double MadVoro::Voronoi3D::Voronoi3DImpl::CalcTetraRadiusCenterHiPrecision(const size_t &index) const
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

void MadVoro::Voronoi3D::Voronoi3DImpl::GetTetraCM(std::array<Point3D, 4> const &points, Point3D &CM) const
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

double MadVoro::Voronoi3D::Voronoi3DImpl::GetTetraVolume(std::array<Point3D, 4> const &points) const
{
    return std::abs(orient3d(points)) / 6.0;
}

/*
void MadVoro::Voronoi3D::CalcCellCMVolume(std::size_t index)
{
    volume_[index] = 0;
    CM_[index] = Point3D();
    std::size_t Nfaces = FacesInCell_[index].size();
    std::array<Point3D, 4> tetra;
    tetra[3] = del_.points_[index];
    Point3D vtemp;
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

void MadVoro::Voronoi3D::Voronoi3DImpl::output(std::string const &filename) const
{

    std::ofstream file_handle(filename.c_str(), std::ios::out | std::ios::binary);
    assert(file_handle.is_open());
    IO::binary_write_single_int(static_cast<int>(Norg_), file_handle);

    // Points
    for (std::size_t i = 0; i < Norg_; ++i)
    {
        IO::binary_write_single_double(del_.points_[i].x, file_handle);
        IO::binary_write_single_double(del_.points_[i].y, file_handle);
        IO::binary_write_single_double(del_.points_[i].z, file_handle);
    }

    IO::binary_write_single_int(static_cast<int>(tetra_centers_.size()), file_handle);
    // Face3D Points
    for (std::size_t i = 0; i < tetra_centers_.size(); ++i)
    {
        IO::binary_write_single_double(tetra_centers_[i].x, file_handle);
        IO::binary_write_single_double(tetra_centers_[i].y, file_handle);
        IO::binary_write_single_double(tetra_centers_[i].z, file_handle);
    }

    // Faces in cell
    for (std::size_t i = 0; i < Norg_; ++i)
    {
        IO::binary_write_single_int(static_cast<int>(FacesInCell_[i].size()), file_handle);
        for (std::size_t j = 0; j < FacesInCell_[i].size(); ++j)
            IO::binary_write_single_int(static_cast<int>(FacesInCell_[i][j]), file_handle);
    }

    // Points in Face3D
    IO::binary_write_single_int(static_cast<int>(PointsInFace_.size()), file_handle);
    for (std::size_t i = 0; i < PointsInFace_.size(); ++i)
    {
        IO::binary_write_single_int(static_cast<int>(PointsInFace_[i].size()), file_handle);
        for (std::size_t j = 0; j < PointsInFace_[i].size(); ++j)
            IO::binary_write_single_int(static_cast<int>(PointsInFace_[i][j]), file_handle);
    }

    file_handle.close();
}

#ifdef MADVORO_WITH_MPI
void MadVoro::Voronoi3D::Voronoi3DImpl::output_buildextra(std::string const &filename) const
{
    std::ofstream file_handle(filename.c_str(), std::ios::out | std::ios::binary);
    assert(file_handle.is_open());
    size_t stemp = Norg_;
    IO::binary_write_single_int(static_cast<int>(stemp), file_handle);
    stemp = del_.points_.size();
    IO::binary_write_single_int(static_cast<int>(stemp), file_handle);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Points
    for (std::size_t i = 0; i < stemp; ++i)
    {
        IO::binary_write_single_double(del_.points_[i].x, file_handle);
        IO::binary_write_single_double(del_.points_[i].y, file_handle);
        IO::binary_write_single_double(del_.points_[i].z, file_handle);
    }

    IO::binary_write_single_int(static_cast<int>(duplicatedprocs_.size()), file_handle);
    // Procs
    assert(duplicatedprocs_.size() == Nghost_.size());
    for (size_t i = 0; i < duplicatedprocs_.size(); ++i)
    {
        IO::binary_write_single_int(static_cast<int>(duplicatedprocs_[i]), file_handle);
        IO::binary_write_single_int(static_cast<int>(Nghost_[i].size()), file_handle);
        for (size_t j = 0; j < Nghost_[i].size(); ++j)
            IO::binary_write_single_int(static_cast<int>(Nghost_[i][j]), file_handle);
    }
    file_handle.close();
}
#endif

size_t MadVoro::Voronoi3D::Voronoi3DImpl::GetContainingCell(const Point3D &point) const
{
    return this->myPointsTree->closestPoint(point).getIndex();
}

std::size_t MadVoro::Voronoi3D::Voronoi3DImpl::GetPointNo(void) const
{
    return Norg_;
}

const Point3D &MadVoro::Voronoi3D::Voronoi3DImpl::GetMeshPoint(std::size_t index) const
{
    return del_.points_[index];
}

double MadVoro::Voronoi3D::Voronoi3DImpl::GetArea(std::size_t index) const
{
    return area_[index];
}

Point3D const &MadVoro::Voronoi3D::Voronoi3DImpl::GetCellCM(std::size_t index) const
{
    return this->CM_[index];
}

std::size_t MadVoro::Voronoi3D::Voronoi3DImpl::GetTotalFacesNumber(void) const
{
    return FaceNeighbors_.size();
}

double MadVoro::Voronoi3D::Voronoi3DImpl::GetWidth(std::size_t index) const
{
    return std::pow(3 * volume_[index] * 0.25 / M_PI, 0.3333333333);
}

double MadVoro::Voronoi3D::Voronoi3DImpl::GetVolume(std::size_t index) const
{
    return volume_[index];
}

face_vec const &MadVoro::Voronoi3D::Voronoi3DImpl::GetCellFaces(std::size_t index) const
{
    return FacesInCell_[index];
}

vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::accessMeshPoints(void)
{
    return del_.points_;
}

const vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::getMeshPoints(void) const
{
    return del_.points_;
}

const MadVoro::Voronoi3D::AllPointsMap &MadVoro::Voronoi3D::Voronoi3DImpl::GetIndicesInAllPoints(void) const
{
    return this->indicesInAllMyPoints;
}

const std::vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::getAllPoints(void) const
{
    return this->allMyPoints;
}

std::vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::getAllPoints(void)
{
    return this->allMyPoints;
}

vector<std::size_t> MadVoro::Voronoi3D::Voronoi3DImpl::GetNeighbors(std::size_t index) const
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

void MadVoro::Voronoi3D::Voronoi3DImpl::GetNeighbors(size_t index, vector<size_t> &res) const
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

MadVoro::Voronoi3D::Voronoi3DImpl *MadVoro::Voronoi3D::Voronoi3DImpl::clone(void) const
{
    return new Voronoi3DImpl(*this);
}

MadVoro::Voronoi3D::Voronoi3DImpl::Voronoi3DImpl(Voronoi3D::Voronoi3DImpl const &other) : ll_(other.ll_), ur_(other.ur_), Norg_(other.Norg_), bigtet_(other.bigtet_),
                                                set_temp_(other.set_temp_), stack_temp_(other.stack_temp_), del_(other.del_), PointTetras_(other.PointTetras_), R_(other.R_),
                                                tetra_centers_(other.tetra_centers_), FacesInCell_(other.FacesInCell_), PointsInFace_(other.PointsInFace_),
                                                FaceNeighbors_(other.FaceNeighbors_), CM_(other.CM_), Face_CM_(other.Face_CM_), volume_(other.volume_), area_(other.area_),
                                                #ifdef MADVORO_WITH_MPI
                                                    sentprocs_(other.sentprocs_), sentpoints_(other.sentpoints_), duplicatedprocs_(other.duplicatedprocs_), duplicated_points_(other.duplicated_points_),
                                                    Nghost_(other.Nghost_), self_index_(other.self_index_),
                                                #endif // MADVORO_WITH_MPI
                                                temp_points_(std::array<Point3D, 4>()), temp_points2_(std::array<Point3D, 5>()), box_faces_(other.box_faces_),
                                                #ifdef MADVORO_WITH_MPI
                                                    pointsManager(other.pointsManager), rangeFinder(other.rangeFinder), radiuses(other.radiuses),
                                                    allMyPoints(other.allMyPoints), allPointsWeights(other.allPointsWeights),
                                                #endif // MADVORO_WITH_MPI
                                                indicesInAllMyPoints(other.indicesInAllMyPoints), verbosity(other.verbosity)
                                                {}

bool MadVoro::Voronoi3D::Voronoi3DImpl::NearBoundary(std::size_t index) const
{
    std::size_t N = FacesInCell_[index].size();
    for (std::size_t i = 0; i < N; ++i)
    {
        if (BoundaryFace3D(FacesInCell_[index][i]))
            return true;
    }
    return false;
}

bool MadVoro::Voronoi3D::Voronoi3DImpl::IsPointOutsideBox(size_t index) const
{
    if(box_faces_.empty())
        return !PointInDomain(ll_, ur_, del_.points_[index]);
    else
        return !PointInPoly(box_faces_, del_.points_[index]);
}

bool MadVoro::Voronoi3D::Voronoi3DImpl::BoundaryFace3D(std::size_t index) const
{
    if (FaceNeighbors_[index].first >= Norg_ || FaceNeighbors_[index].second >= Norg_)
    {
#ifdef MADVORO_WITH_MPI
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
    return false;
}

#ifdef MADVORO_WITH_MPI
    vector<vector<std::size_t>> &MadVoro::Voronoi3D::Voronoi3DImpl::GetDuplicatedPoints(void)
    {
        return duplicated_points_;
    }

    vector<vector<std::size_t>> const &MadVoro::Voronoi3D::Voronoi3DImpl::GetDuplicatedPoints(void) const
    {
        return duplicated_points_;
    }
#endif // MADVORO_WITH_MPI

std::size_t MadVoro::Voronoi3D::Voronoi3DImpl::GetTotalPointNumber(void) const
{
    return del_.points_.size();
}

vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllCM(void)
{
    return CM_;
}

vector<Point3D> MadVoro::Voronoi3D::Voronoi3DImpl::GetAllCM(void) const
{
    return CM_;
}

void MadVoro::Voronoi3D::Voronoi3DImpl::GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point) const
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
    result = Utils::unique(result);
    result = Utils::RemoveList(result, neigh);
    Utils::RemoveVal(result, point);
}

vector<boost::container::small_vector<size_t, 8>> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllPointsInFace(void)
{
    return PointsInFace_;
}

vector<boost::container::small_vector<size_t, 8>> const& MadVoro::Voronoi3D::Voronoi3DImpl::GetAllPointsInFace(void)const
{
    return PointsInFace_;
}

size_t &MadVoro::Voronoi3D::Voronoi3DImpl::GetPointNo(void)
{
    return Norg_;
}

size_t MadVoro::Voronoi3D::Voronoi3DImpl::GetAllPointsNo(void) const
{
    return this->allMyPoints.size();
}

std::vector<std::pair<size_t, size_t>> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllFaceNeighbors(void)
{
    return FaceNeighbors_;
}

const std::vector<std::pair<size_t, size_t>> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllFaceNeighbors(void) const
{
    return FaceNeighbors_;
}

vector<double> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllVolumes(void)
{
    return volume_;
}

vector<double> MadVoro::Voronoi3D::Voronoi3DImpl::GetAllVolumes(void) const
{
    return volume_;
}

Point3D MadVoro::Voronoi3D::Voronoi3DImpl::Normal(std::size_t faceindex) const
{
    return del_.points_[FaceNeighbors_[faceindex].second] - del_.points_[FaceNeighbors_[faceindex].first];
}

bool MadVoro::Voronoi3D::Voronoi3DImpl::IsGhostPoint(std::size_t index) const
{
    return index >= Norg_;
}

const Point3D &MadVoro::Voronoi3D::Voronoi3DImpl::FaceCM(std::size_t index) const
{
    return Face_CM_[index];
}

Point3D MadVoro::Voronoi3D::Voronoi3DImpl::CalcFaceVelocity(std::size_t index, Point3D const &v0, Point3D const &v1) const
{
    std::size_t p0 = FaceNeighbors_[index].first;
    std::size_t p1 = FaceNeighbors_[index].second;
    Point3D r0 = GetMeshPoint(p0);
    Point3D r1 = GetMeshPoint(p1);
    Point3D r_diff = r1 - r0;
    double abs_r_diff = ScalarProd(r_diff, r_diff);

    Point3D f = FaceCM(index);
    r1 += r0;
    r1 *= 0.5;
    f -= r1;
    Point3D delta_w = ScalarProd((v0 - v1), f) * r_diff / abs_r_diff;
#ifdef MADVORO_DEBUG
    double dw_abs = fastabs(delta_w);
#endif // MADVORO_DEBUG
    Point3D w = (v0 + v1) * 0.5;
#ifdef MADVORO_DEBUG
    //    double w_abs = std::max(fastabs(v0),fastabs(v1));
#endif // MADVORO_DEBUG
    //if (dw_abs > w_abs)
    //	delta_w *= (1 + (std::atan(dw_abs / w_abs) - 0.25 * M_PI)*2) * (w_abs / dw_abs);
#ifdef MADVORO_DEBUG
    if (!std::isfinite(dw_abs))
    {
        r0 = GetMeshPoint(p0);
        r1 = GetMeshPoint(p1);
        f = FaceCM(index);
        MadVoro::Exception::MadVoroException eo("Bad Face3D velocity");
        eo.addEntry("Face3D index", index);
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
        eo.addEntry("Face3D CMx", f.x);
        eo.addEntry("Face3D CMy", f.y);
        eo.addEntry("Face3D CMz", f.z);
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

vector<double> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllArea(void)
{
    return area_;
}

const vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllFaceCM(void) const
{
    return Face_CM_;
}

vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllFaceCM(void)
{
    return Face_CM_;
}

vector<face_vec> &MadVoro::Voronoi3D::Voronoi3DImpl::GetAllCellFaces(void)
{
    return FacesInCell_;
}

vector<face_vec> const& MadVoro::Voronoi3D::Voronoi3DImpl::GetAllCellFaces(void)const
{
    return FacesInCell_;
}

vector<Point3D> &MadVoro::Voronoi3D::Voronoi3DImpl::GetFacePoints(void)
{
    return tetra_centers_;
}

vector<Point3D> const &MadVoro::Voronoi3D::Voronoi3DImpl::GetFacePoints(void) const
{
    return tetra_centers_;
}

point_vec const &MadVoro::Voronoi3D::Voronoi3DImpl::GetPointsInFace(std::size_t index) const
{
    return PointsInFace_[index];
}

const std::pair<std::size_t, std::size_t> &MadVoro::Voronoi3D::Voronoi3DImpl::GetFaceNeighbors(std::size_t face_index) const
{
    return FaceNeighbors_[face_index];
}

#ifdef MADVORO_WITH_MPI
    vector<int> MadVoro::Voronoi3D::Voronoi3DImpl::GetDuplicatedProcs(void) const
    {
        return duplicatedprocs_;
    }

    vector<int> MadVoro::Voronoi3D::Voronoi3DImpl::GetSentProcs(void) const
    {
        return sentprocs_;
    }

    vector<vector<std::size_t>> const &MadVoro::Voronoi3D::Voronoi3DImpl::GetSentPoints(void) const
    {
        return sentpoints_;
    }

    vector<std::size_t> const &MadVoro::Voronoi3D::Voronoi3DImpl::GetSelfIndex(void) const
    {
        return self_index_;
    }

    vector<int> &MadVoro::Voronoi3D::Voronoi3DImpl::GetSentProcs(void)
    {
        return sentprocs_;
    }

    vector<vector<std::size_t>> &MadVoro::Voronoi3D::Voronoi3DImpl::GetSentPoints(void)
    {
        return sentpoints_;
    }

    vector<std::size_t> &MadVoro::Voronoi3D::Voronoi3DImpl::GetSelfIndex(void)
    {
        return self_index_;
    }
#endif // MADVORO_WITH_MPI

void MadVoro::Voronoi3D::Voronoi3DImpl::SetBox(const Point3D &ll, const Point3D &ur)
{
    this->ll_ = ll;
    this->ur_ = ur;
    #ifdef MADVORO_WITH_MPI
        this->pointsManager = std::shared_ptr<PointsManager>();
        // this->radiuses.clear();
    #endif // MADVORO_WITH_MPI
}

#ifdef MADVORO_WITH_MPI
    const std::vector<double> &MadVoro::Voronoi3D::Voronoi3DImpl::GetPointsBuildWeights() const
    {
        return this->allPointsWeights;
    }

    const EnvironmentAgent *MadVoro::Voronoi3D::Voronoi3DImpl::GetEnvironmentAgent() const
    {
        return this->pointsManager->getEnvironmentAgent();
    }
#endif // MADVORO_WITH_MPI

#ifdef MADVORO_WITH_HDF5

#endif // MADVORO_WITH_HDF5

/* ========================= Interface ========================= */

template<typename T>
std::vector<typename T::impl> translateImpl(const std::vector<T> &values)
{
    std::vector<typename T::impl> result;
    for(const auto &value : values)
    {
        result.push_back(T::impl(value));
    }
    return result;
}

template<typename T>
std::vector<T> translateToImpl(const std::vector<typename T::impl> &values)
{
    std::vector<T> result;
    for(const auto &value : values)
    {
        result.push_back(value.pImpl->get());
    }
    return result;
}


Vector3D pointToVector(const Point3D &v)
{
    return Vector3D(v.x, v.y, v.z);
}

Point3D vectorToPoint(const Vector3D &p)
{
    return Point3D(p.x, p.y, p.z);
}

std::vector<Vector3D> pointsToVectors(const std::vector<Point3D> &vectors)
{
    std::vector<Vector3D> result;
    result.reserve(vectors.size());
    for(const Point3D &v : vectors)
    {
        result.emplace_back(v.x, v.y, v.z);
    }
    return result;
}

std::vector<Point3D> vectorsToPoints(const std::vector<Vector3D> &points)
{
    std::vector<Point3D> result;
    result.reserve(points.size());
    for(const Vector3D &p : points)
    {
        result.emplace_back(p.x, p.y, p.z);
    }
    return result;
}

std::vector<Face3D> facesToFaces3D(const std::vector<Face> &faces)
{
    std::vector<Face3D> result;
    for(const Face &face : faces)
    {
        result.emplace_back();
        Face3D &f = result.back();
        f.neighbors = face.neighbors;
        for(const Vector3D &p : face.vertices)
        {
            f.vertices.push_back(vectorToPoint(p));
        }
    }
    return result;
}

std::vector<Face> faces3DToFaces(const std::vector<Face3D> &faces)
{
    std::vector<Face> result;
    for(const Face3D &face : faces)
    {
        result.emplace_back();
        Face &f = result.back();
        f.neighbors = face.neighbors;
        for(const Point3D &p : face.vertices)
        {
            f.vertices.push_back(pointToVector(p));
        }
    }
    return result;
}

MadVoro::Voronoi3D::Voronoi3D(const Vector3D &ll, const Vector3D &ur): pImpl(new Voronoi3DImpl(vectorToPoint(ll), vectorToPoint(ur)))
{}

MadVoro::Voronoi3D::Voronoi3D(const std::vector<Face> &box_faces): pImpl(new Voronoi3DImpl(facesToFaces3D(box_faces)))
{}

MadVoro::Voronoi3D::~Voronoi3D()
{
    delete this->pImpl;
}

std::vector<Vector3D> MadVoro::Voronoi3D::GetAllFaceCM(void) const
{
    return pointsToVectors(this->pImpl->GetAllFaceCM());
}

void MadVoro::Voronoi3D::BuildPartially(const std::vector<Vector3D> &allPoints, const std::vector<std::size_t> &indicesToBuild)
{
    this->pImpl->BuildPartially(vectorsToPoints(allPoints), indicesToBuild);
}

void MadVoro::Voronoi3D::Build(const std::vector<Vector3D> &points)
{
    this->pImpl->Build(vectorsToPoints(points));
}

#ifdef MADVORO_WITH_MPI
    const std::vector<double> &MadVoro::Voronoi3D::GetPointsBuildWeights() const
    {
        return this->pImpl->GetPointsBuildWeights();
    }
        
    std::vector<Vector3D> MadVoro::Voronoi3D::BuildParallel(const std::vector<Vector3D> &points, const std::vector<double> &weights, bool suppressRebalancing)
    {
        return pointsToVectors(this->pImpl->BuildParallel(vectorsToPoints(points), weights, suppressRebalancing));
    }

    std::vector<Vector3D> MadVoro::Voronoi3D::BuildParallel(const std::vector<Vector3D> &points, bool suppressRebalancing)
    {
        return pointsToVectors(this->pImpl->BuildParallel(vectorsToPoints(points), suppressRebalancing));
    }

    std::vector<Vector3D> MadVoro::Voronoi3D::BuildPartiallyParallel(const std::vector<Vector3D> &allPoints, const std::vector<double> &allWeights, const std::vector<std::size_t> &indicesToBuild, bool suppressRebalancing)
    {
        return pointsToVectors(this->pImpl->BuildPartiallyParallel(vectorsToPoints(allPoints), allWeights, indicesToBuild, suppressRebalancing));
    }

    bool MadVoro::Voronoi3D::PointInMyDomain(const Vector3D &point) const
    {
        return this->pImpl->PointInMyDomain(vectorToPoint(point));
    }

    int MadVoro::Voronoi3D::GetOwner(const Vector3D &point) const
    {
        return this->pImpl->GetOwner(vectorToPoint(point));
    }
#endif // MADVORO_WITH_MPI

std::size_t MadVoro::Voronoi3D::GetContainingCell(const Vector3D &point) const
{
    return this->pImpl->GetContainingCell(vectorToPoint(point));
}

Vector3D MadVoro::Voronoi3D::FaceCM(std::size_t index) const
{
    return pointToVector(this->pImpl->FaceCM(index));
}

std::size_t MadVoro::Voronoi3D::GetPointNo(void) const
{
    return this->pImpl->GetPointNo();
}

Vector3D MadVoro::Voronoi3D::GetMeshPoint(std::size_t index) const
{
    return pointToVector(this->pImpl->GetMeshPoint(index));
}

double MadVoro::Voronoi3D::GetArea(std::size_t faceIndex) const
{
    return this->pImpl->GetArea(faceIndex);
}

Vector3D MadVoro::Voronoi3D::GetCellCM(std::size_t index) const
{
    return pointToVector(this->pImpl->GetCellCM(index));
}

std::size_t MadVoro::Voronoi3D::GetTotalFacesNumber(void) const
{
    return this->pImpl->GetTotalFacesNumber();
}

double MadVoro::Voronoi3D::GetWidth(std::size_t index) const
{
    return this->pImpl->GetWidth(index);
}

double MadVoro::Voronoi3D::GetVolume(std::size_t index) const
{
    return this->pImpl->GetVolume(index);
}

const face_vec &MadVoro::Voronoi3D::GetCellFaces(std::size_t index) const
{
    return this->pImpl->GetCellFaces(index);
}

std::vector<Vector3D> MadVoro::Voronoi3D::getMeshPoints(void) const
{
    return pointsToVectors(this->pImpl->getMeshPoints());
}

const MadVoro::Voronoi3D::AllPointsMap &MadVoro::Voronoi3D::GetIndicesInAllPoints(void) const
{
    return this->pImpl->GetIndicesInAllPoints();
}

std::vector<Vector3D> MadVoro::Voronoi3D::getAllPoints(void) const
{
    return pointsToVectors(this->pImpl->getAllPoints());
}

std::size_t MadVoro::Voronoi3D::GetAllPointsNo(void) const
{
    return this->pImpl->GetAllPointsNo();
}

std::vector<std::size_t> MadVoro::Voronoi3D::GetNeighbors(std::size_t index) const
{
    return this->pImpl->GetNeighbors(index);
};

bool MadVoro::Voronoi3D::NearBoundary(std::size_t index) const
{
    return this->pImpl->NearBoundary(index);
};

bool MadVoro::Voronoi3D::BoundaryFace(std::size_t index) const
{
    return this->pImpl->BoundaryFace3D(index);
};

#ifdef MADVORO_WITH_MPI
    // Communication methods
    const std::vector<std::vector<std::size_t>> &MadVoro::Voronoi3D::GetDuplicatedPoints(void) const
    {
        return this->pImpl->GetDuplicatedPoints();
    }

    std::vector<int> MadVoro::Voronoi3D::GetDuplicatedProcs(void) const
    {
        return this->pImpl->GetDuplicatedProcs();
    }

    std::vector<int> MadVoro::Voronoi3D::GetSentProcs(void) const
    {
        return this->pImpl->GetSentProcs();
    }

    const std::vector<std::vector<std::size_t>> &MadVoro::Voronoi3D::GetSentPoints(void) const
    {
        return this->pImpl->GetSentPoints();
    }

    const std::vector<std::size_t> &MadVoro::Voronoi3D::GetSelfIndex(void) const
    {
        return this->pImpl->GetSelfIndex();
    }

    const std::vector<std::vector<std::size_t>> &MadVoro::Voronoi3D::GetGhostIndeces(void) const
    {
        return this->pImpl->GetGhostIndeces();
    }

    std::vector<std::vector<std::size_t>> &MadVoro::Voronoi3D::GetGhostIndeces(void)
    {
        return this->pImpl->GetGhostIndeces();
    }
#endif // MADVORO_WITH_MPI

std::size_t MadVoro::Voronoi3D::GetTotalPointNumber(void) const
{
    return this->pImpl->GetTotalPointNumber();
};

std::vector<Vector3D> MadVoro::Voronoi3D::GetAllCM(void) const
{
    return pointsToVectors(this->pImpl->GetAllCM());
}

Vector3D MadVoro::Voronoi3D::Normal(std::size_t faceindex) const
{
    return pointToVector(this->pImpl->Normal(faceindex));
}

bool MadVoro::Voronoi3D::IsGhostPoint(std::size_t index) const
{
    return this->pImpl->IsGhostPoint(index);
}

std::vector<Vector3D> MadVoro::Voronoi3D::GetFacePoints(void) const
{
    return pointsToVectors(this->pImpl->GetFacePoints());
}

const std::vector<face_vec> &MadVoro::Voronoi3D::GetAllCellFaces(void) const
{
    return this->pImpl->GetAllCellFaces();
}

const point_vec &MadVoro::Voronoi3D::GetPointsInFace(std::size_t index) const
{
    return this->pImpl->GetPointsInFace(index);
}

const std::pair<std::size_t, std::size_t> &MadVoro::Voronoi3D::GetFaceNeighbors(std::size_t face_index) const
{
    return this->pImpl->GetFaceNeighbors(face_index);
}

void MadVoro::Voronoi3D::GetNeighbors(size_t index, std::vector<std::size_t> &res) const
{
    return this->pImpl->GetNeighbors(index, res);
}

std::pair<Vector3D, Vector3D> MadVoro::Voronoi3D::GetBoxCoordinates(void) const
{
    const std::pair<Point3D, Point3D> res = this->pImpl->GetBoxCoordinates();
    return std::make_pair(pointToVector(res.first), pointToVector(res.second));
}

std::vector<double> MadVoro::Voronoi3D::GetAllVolumes(void) const
{
    return this->pImpl->GetAllVolumes();
}

const std::vector<std::pair<std::size_t, std::size_t>> &MadVoro::Voronoi3D::GetAllFaceNeighbors(void) const
{
    return this->pImpl->GetAllFaceNeighbors();
}

const std::vector<point_vec> &MadVoro::Voronoi3D::GetAllPointsInFace(void) const
{
    return this->pImpl->GetAllPointsInFace();
}

void MadVoro::Voronoi3D::SetBox(const Vector3D &ll, const Vector3D &ur)
{
    this->pImpl->SetBox(vectorToPoint(ll), vectorToPoint(ur));
}

std::vector<Face> MadVoro::Voronoi3D::GetBoxFaces(void) const
{
    return faces3DToFaces(this->pImpl->GetBoxFaces());
}

void MadVoro::Voronoi3D::SetVerbosity(bool value)
{
    this->pImpl->SetVerbosity(value);
}

#ifdef MADVORO_WITH_HDF5
    void MadVoro::Voronoi3D::ToHDF5(const std::string &fileName, const std::vector<std::string> &fieldNames, const std::vector<std::vector<double>> &fieldValues)
    {
        MadVoro::IO::WriteVoronoiHDF5(*this, fileName, fieldValues, fieldNames);
    }
#endif // MADVORO_WITH_HDF5

#ifdef MADVORO_WITH_VTK
    void MadVoro::Voronoi3D::ToVTK(const std::string &fileName, const std::vector<std::string> &fieldNames, const std::vector<std::vector<double>> &fieldValues)
    {
        MadVoro::IO::WriteVoronoiVTK(*this, fileName, fieldValues, fieldNames);
    }
#endif // MADVORO_WITH_VTK