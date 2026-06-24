#ifndef VORONOI3D_FULL_HPP
#define VORONOI3D_FULL_HPP

#include <vector>
#include <memory>
#include <set>
#include <stack>
#include <array>
#include <string>
#include <unordered_set>
#include <numeric>
#include <cassert>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>

#include "delaunay/Delaunay3D.hpp"
#include "elementary/Face3D.hpp"
#include "elementary/Point3D.hpp"
#include "madvoro/types.hpp"
#include "ds/utils/Sphere.hpp"
#include "range/SmallRangeAgent.hpp"
#include "range/BigRangeAgent.hpp"
#include "exception/MadVoroException.hpp"

#ifdef MADVORO_WITH_MPI
  #include "mpi/serialize/mpi_commands.hpp"
#endif

namespace MadVoro
{
  class LoadBalancer;
  class EnvironmentAgent;
  class PointsManager;

  namespace Range { class RangeFinder; struct IndexedPoint3D; }
  namespace DataStructure { template<typename T> class OctTree; }
}

namespace MadVoro
{
  typedef std::array<std::size_t, 4> b_array_4;
  typedef std::array<std::size_t, 3> b_array_3;
  typedef boost::container::small_vector<std::size_t, 40> tetra_vec;

  class Voronoi3DFull
  {
    using IndexedPointsTree = DataStructure::OctTree<Range::IndexedPoint3D>;

  private:
    Point3D ll_, ur_;
    std::size_t Norg_, bigtet_;

    std::set<int> set_temp_;
    std::stack<int> stack_temp_;

    void FindIntersectionsSingle(std::vector<Face3D> const& box, std::size_t point, Geometry::Sphere<Point3D> &sphere,
            std::vector<std::size_t> &intersecting_faces, std::vector<double> &Rtemp, std::vector<Point3D> &vtemp);

    std::size_t GetFirstPointToCheck(void)const;

    void GetPointToCheck(std::size_t point, std::vector<unsigned char> const& checked, std::vector<std::size_t> &res);

    void CalcRigidCM(std::size_t face_index);

    void GetTetraCM(std::array<Point3D, 4> const& points, Point3D &CM) const;

    double GetTetraVolume(std::array<Point3D, 4> const& points)const;

    double GetRadius(const std::size_t &index) const;

    void CalcAllCM(void);

    std::vector<std::pair<std::size_t, std::size_t>> SerialFindIntersections(bool first_run);

    std::vector<std::pair<std::size_t, std::size_t>> SerialFirstIntersections(void);

    double CalcTetraRadiusCenterHiPrecision(const std::size_t &index) const;

    double CalcTetraRadiusCenter(const std::size_t &index) const;

    std::vector<Point3D> CreateBoundaryPoints(std::vector<std::pair<std::size_t, std::size_t>> const& to_duplicate,
            std::vector<std::vector<std::size_t>> &past_duplicate);
    void BuildVoronoi(std::vector<std::size_t> const& order);

    void InitialBoxBuild(std::vector<Face3D> &box, std::vector<Point3D> &normals);

    void BringSelfGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                boost::container::flat_map<std::size_t, std::size_t> &numOfResultsForBigPoints,
                                boost::container::flat_map<std::size_t, std::size_t> &numOfResultsForSmallPoints,
                                std::unordered_set<std::size_t> &selfIgnorePoints);

    #ifdef MADVORO_WITH_MPI
    void BringGhostPointsToBuild(const MPI_Comm &comm);
    #else
    void BringGhostPointsToBuild();
    #endif

    std::pair<std::vector<SmallRangeQueryData>, std::vector<BigRangeQueryData>> CreateBatches(
        boost::container::flat_set<std::size_t> &smallPoints, boost::container::flat_set<std::size_t> &largePoints,
        const boost::container::flat_map<std::size_t, std::size_t> &firstLargeIteration,
        std::vector<double> &currentRadiuses, std::size_t iterations);

    std::pair<boost::container::flat_set<std::size_t>, boost::container::flat_set<std::size_t>>
    DetermineNextIterationPoints(std::size_t iterations,
                                    boost::container::flat_map<std::size_t, std::size_t> &firstLargeIteration,
                                    std::vector<double> &currentRadiuses,
                                    const boost::container::flat_map<std::size_t, std::size_t> &resultOfSmallPoints,
                                    const boost::container::flat_map<std::size_t, std::size_t> &resultOfBigPoints);

    void UpdateRadiuses(const std::vector<Point3D> &points);

    void UpdateCMs(void);

    void UpdateRangeFinder(void);

    std::size_t SetPointTetras(void);

    #ifdef MADVORO_WITH_MPI
    std::vector<Point3D> PrepareToBuildParallel(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<std::size_t> &indicesToBuild, bool suppressRebalancing);
    void FilterRealGhostPoints();
    void UpdateDuplicatedPoints(const std::vector<int> &sentProc, const std::vector<std::vector<std::size_t>> &sentPoints);
    void EnsureSymmetry(const std::vector<int> &sentProc, const std::vector<std::vector<int>> &recvProcLists);
    std::tuple<std::vector<Point3D>, std::vector<int>, std::vector<std::vector<std::size_t>>, std::vector<int>, std::vector<std::vector<std::size_t>>> InitialGhostPointsExchange(const MPI_Comm &comm = MPI_COMM_WORLD) const;
    void InitialExchange(const std::vector<Point3D> &points, std::vector<int> &sentProc, std::vector<std::vector<std::size_t>> &sentPoints, const MPI_Comm &comm = MPI_COMM_WORLD);
    void SetGhostArray(const std::vector<int> &recvProc, const std::vector<std::vector<std::size_t>> &recvPoints);
    void BringRemoteGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                        BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                        boost::container::flat_map<std::size_t, std::size_t> &numOfResultsForBigPoints,
                                        boost::container::flat_map<std::size_t, std::size_t> &numOfResultsForSmallPoints);
    #endif

    Delaunay3D del_;
    std::vector<tetra_vec> PointTetras_;
    mutable std::vector<double> R_;
    mutable std::vector<Point3D> tetra_centers_;
    std::vector<face_vec> FacesInCell_;
    std::vector<point_vec> PointsInFace_;
    std::vector<std::pair<std::size_t, std::size_t>> FaceNeighbors_;
    std::vector<Point3D> all_CM;
    std::vector<Point3D> CM_, Face_CM_;
    std::vector<double> volume_;
    std::vector<double> area_;

    #ifdef MADVORO_WITH_MPI
    std::vector<int> sentprocs_;
    std::vector<std::vector<std::size_t>> sentpoints_;
    std::vector<int> duplicatedprocs_;
    std::vector<std::vector<std::size_t>> duplicated_points_;
    std::vector<int> real_duplicated_proc;
    std::vector<std::vector<std::size_t>> real_duplicated_points;
    std::vector<std::vector<std::size_t>> Nghost_;
    std::vector<std::size_t> self_index_;
    #endif

    Voronoi3DFull();
    Voronoi3DFull(Voronoi3DFull const &other);
    std::array<Point3D, 4> temp_points_;
    std::array<Point3D, 5> temp_points2_;
    std::vector<Face3D> box_faces_;

    std::shared_ptr<IndexedPointsTree> myPointsTree;
    std::shared_ptr<IndexedPointsTree> allMyPointsTree;
    #ifdef MADVORO_WITH_MPI
    std::shared_ptr<PointsManager> pointsManager;
    #endif

    std::shared_ptr<Range::RangeFinder> rangeFinder;
    std::vector<Point3D> allMyPoints;
    std::vector<double> allPointsWeights;
    std::vector<double> radiuses;

    AllPointsMap indicesInAllMyPoints;
    bool verbosity;

  public:
    #ifdef MADVORO_WITH_MPI
    const std::vector<double> &GetPointsBuildWeights() const;

    const EnvironmentAgent *GetEnvironmentAgent() const;

    std::vector<Point3D> BuildParallel(const std::vector<Point3D> &points, const std::vector<double> &weights, bool suppressRebalancing = false)
    {
        std::vector<std::size_t> indicesToBuild(points.size());
        std::iota(indicesToBuild.begin(), indicesToBuild.end(), 0);
        return this->BuildPartiallyParallel(points, weights, indicesToBuild, suppressRebalancing);
    }

    inline std::vector<Point3D> BuildParallel(const std::vector<Point3D> &points, bool suppressRebalancing = false)
    {
        return this->BuildParallel(points, std::vector<double>(points.size(), 1.0), suppressRebalancing);
    }
    #endif

    #ifdef MADVORO_WITH_MPI
    std::vector<int>& GetSentProcs(void);
    std::vector<std::vector<std::size_t>>& GetSentPoints(void);
    std::vector<std::size_t>& GetSelfIndex(void);
    #endif

    std::vector<Point3D>& GetAllFaceCM(void);
    const std::vector<Point3D>& GetAllFaceCM(void) const;
    const Point3D &FaceCM(std::size_t index)const;

    Voronoi3DFull(Point3D const& ll, Point3D const& ur);
    Voronoi3DFull(std::vector<Face3D> const& box_faces);

    void output(std::string const& filename)const;
    void BuildInitialize(std::size_t num_points);
    void BuildPartially(const std::vector<Point3D> &allPoints, const std::vector<std::size_t> &indicesToBuild);
    void Build(const std::vector<Point3D> &points);

    #ifdef MADVORO_WITH_MPI
    void output_buildextra(std::string const& filename) const;
    void PreparePoints(const std::vector<Point3D> &points, const std::vector<std::size_t> &mask);
    std::vector<Point3D> BuildPartiallyParallel(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<std::size_t> &indicesToBuild, bool suppressRebalancing = false);
    void MockMesh(void);
    void SetLoadBalancer(std::shared_ptr<LoadBalancer> loadBalancer);
    void Rebalance(const std::vector<double> &weights);
    void SetImbalanceTolerance(double tolerance);
    bool PointInMyDomain(const Point3D &point) const;
    int GetOwner(const Point3D &point) const;
    #endif

    void BuildDebug(int rank);
    double GetMaxRadius(const std::size_t &index) const;
    double GetMinRadius(const std::size_t &index) const;
    std::size_t GetContainingCell(const Point3D &point) const;
    std::size_t GetPointNo(void) const;
    const Point3D &GetMeshPoint(std::size_t index) const;
    double GetArea(std::size_t index) const;
    Point3D const& GetCellCM(std::size_t index) const;
    std::size_t GetTotalFacesNumber(void) const;
    double GetWidth(std::size_t index) const;
    double GetVolume(std::size_t index) const;
    face_vec const& GetCellFaces(std::size_t index) const;
    std::vector<Point3D>& accessMeshPoints(void);
    const std::vector<Point3D>& getMeshPoints(void) const;
    const AllPointsMap &GetIndicesInAllPoints(void) const;
    const std::vector<Point3D> &getAllPoints(void) const;
    std::vector<Point3D> &getAllPoints(void);
    std::size_t GetAllPointsNo(void) const;
    std::vector<std::size_t> GetNeighbors(std::size_t index)const;
    Voronoi3DFull* clone(void) const;
    bool NearBoundary(std::size_t index) const;
    bool BoundaryFace3D(std::size_t index) const;

    #ifdef MADVORO_WITH_MPI
    std::vector<std::vector<std::size_t>>& GetDuplicatedPoints(void);
    std::vector<std::vector<std::size_t>> const& GetDuplicatedPoints(void)const;
    std::vector<int> GetDuplicatedProcs(void)const;
    std::vector<int> GetSentProcs(void)const;
    std::vector<std::vector<std::size_t>> const& GetSentPoints(void)const;
    std::vector<std::size_t> const& GetSelfIndex(void) const;
    #endif

    std::size_t GetTotalPointNumber(void)const;
    std::vector<Point3D> & GetAllCM(void);
    std::vector<Point3D> GetAllCM(void)const;
    void GetNeighborNeighbors(std::vector<std::size_t> &result, std::size_t point)const;
    Point3D Normal(std::size_t faceindex)const;
    bool IsGhostPoint(std::size_t index)const;
    Point3D CalcFaceVelocity(std::size_t index, Point3D const& v0, Point3D const& v1)const;
    std::vector<Point3D>& GetFacePoints(void);
    std::vector<double>& GetAllArea(void);
    std::vector<Point3D>const& GetFacePoints(void) const;
    std::vector<face_vec>& GetAllCellFaces(void);
    std::vector<face_vec> const& GetAllCellFaces(void) const;
    point_vec const& GetPointsInFace(std::size_t index) const;
    const std::pair<std::size_t, std::size_t> &GetFaceNeighbors(std::size_t face_index) const;

    #ifdef MADVORO_WITH_MPI
    std::vector<std::vector<std::size_t>> const& GetGhostIndeces(void) const;
    std::vector<std::vector<std::size_t>>& GetGhostIndeces(void);
    #endif

    void GetNeighbors(std::size_t index, std::vector<std::size_t> &res) const;
    std::pair<Point3D, Point3D> GetBoxCoordinates(void) const;
    void BuildNoBox(std::vector<Point3D> const& points, std::vector<std::vector<Point3D>> const& ghosts, std::vector<std::size_t> toduplicate);
    std::vector<double>& GetAllVolumes(void);
    std::vector<double> GetAllVolumes(void)const;
    std::vector<std::pair<std::size_t, std::size_t>> &GetAllFaceNeighbors(void);
    const std::vector<std::pair<std::size_t, std::size_t>> &GetAllFaceNeighbors(void) const;
    std::vector<point_vec> & GetAllPointsInFace(void);
    std::vector<point_vec> const& GetAllPointsInFace(void) const;
    std::size_t& GetPointNo(void);
    bool IsPointOutsideBox(std::size_t index) const;
    void SetBox(Point3D const& ll, Point3D const& ur);
    std::vector<Face3D> GetBoxFaces(void) const { return box_faces_; }
    std::vector<Face3D>& ModifyBoxFaces(void) { return box_faces_; }

    template<typename T>
    void SyncPartialBuildData(std::vector<T> &partialBuildData, std::vector<T> &allBuildData) const;

    inline void SetVerbosity(bool value) { this->verbosity = value; }

    bool PointInPolyTess(Point3D const &point, std::size_t index);

    #ifdef MADVORO_WITH_HDF5
    void ToHDF5(const std::string &fileName, const std::vector<std::string> &fieldNames, const std::vector<std::vector<double>> &fieldValues);
    #endif
    #ifdef MADVORO_WITH_VTK
    void ToVTK(const std::string &fileName, const std::vector<std::string> &fieldNames, const std::vector<std::vector<double>> &fieldValues);
    #endif
  };

  template<typename T>
  inline void Voronoi3DFull::SyncPartialBuildData(std::vector<T> &partialBuildData, std::vector<T> &allBuildData) const
  {
    std::size_t Norg = this->GetPointNo();
    if(partialBuildData.size() < Norg)
    {
      MadVoro::Exception::MadVoroException eo("Voronoi3D::SyncPartialBuildData: Partial build data has lower size than the number of points");
      eo.addEntry("Partial build data size", partialBuildData.size());
      eo.addEntry("Number of points", Norg);
      throw eo;
    }
    const AllPointsMap &indicesInAllMyPoints = this->GetIndicesInAllPoints();

    allBuildData.resize(this->GetAllPointsNo());

    for(std::size_t i = 0; i < Norg; i++)
    {
        std::size_t pointIdx = indicesInAllMyPoints.at(i);
        allBuildData[pointIdx] = partialBuildData[i];
    }

    std::size_t sizeOfMeshPoints = this->getMeshPoints().size();
    partialBuildData.resize(sizeOfMeshPoints);
    for(std::size_t i = Norg; i < sizeOfMeshPoints; i++)
    {
        bool pointIsMine = (indicesInAllMyPoints.find(i) != indicesInAllMyPoints.cend());
        if(pointIsMine)
        {
            std::size_t pointIdx = indicesInAllMyPoints.at(i);
            partialBuildData[i] = allBuildData[pointIdx];
        }
    }

    #ifdef MADVORO_WITH_MPI
        std::vector<std::vector<T>> incoming = MPI::MPI_exchange_data_indexed(this->GetDuplicatedProcs(), allBuildData, this->GetDuplicatedPoints());
        std::size_t incomingSize = incoming.size();
        const std::vector<std::vector<std::size_t>> &Nghost = this->GetGhostIndeces();
        assert(this->GetDuplicatedProcs().size() == Nghost.size());
        assert(incomingSize == Nghost.size());
        for (std::size_t i = 0; i < incomingSize; ++i)
        {
            std::size_t _size = incoming[i].size();
            assert(_size == Nghost[i].size());
            for (std::size_t j = 0; j < _size; ++j)
            {
                partialBuildData[Nghost.at(i).at(j)] = incoming[i][j];
            }
        }
    #endif
  }
}

#endif // VORONOI3D_FULL_HPP
