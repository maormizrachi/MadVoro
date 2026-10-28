/** \file Voronoi3D.hpp
   \brief A 3D Voronoi
   \author Elad Steinberg, Maor Mizrachi
*/
#ifndef VORONOI3D_HPP
#define VORONOI3D_HPP 1

#if  _MSC_VER
#define _USE_MATH_DEFINES
#endif // _MSC_VER

#include <algorithm>
#include <numeric>
#include <cfloat>
#include <stack>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <set>
#include <map>
#include <array>
#include <tuple>
#include <limits>
#include <unordered_set>
#include <type_traits>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/container/small_vector.hpp>

#ifdef USE_VCL_VECTORIZATION
  #include <vectorclass.h>
#endif // USE_VCL_VECTORIZATION

#ifdef MADVORO_WITH_MPI
  #include <mpi.h>
#endif // MADVORO_WITH_MPI

#include <optional>
#include <exception>
#include "VoronoiPayload.hpp"
#include "types.hpp"
#include "elementary/Face3D.hpp"
#include "delaunay/Delaunay3D.hpp"
#include "geometry/Intersections.hpp"
#include "geometry/GeometryExecutor.hpp"
#include "geometry/kernels/Circumcenter.hpp"
#include "geometry/kernels/MirrorCM.hpp"
#include "utils/Predicates3D.hpp"
#include <MeshDecomposer3D/hilbert/HilbertOrder3D.hpp>
#include "subdomains/MergeSubdomains.hpp"
#include "subdomains/SubdomainSearch.hpp"
#include "ghosts/GhostBringer.hpp"
#include "builder/BuildDomain.hpp"
#include "builder/ClippingLocalBuilder.hpp"
#include "builder/DelaunayLocalBuilder.hpp"
#include "builder/DistributedBuilder.hpp"
#include "utils/container_utils.hpp"
#include <spatial_ds/utils/BoundingBox.hpp>
#include "exception/MadVoroException.hpp"
#include "exception/InvalidArgumentException.hpp"
#include "elementary/PointOps.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi_utils/mpi_exchange.hpp>
#include <MeshDecomposer3D/kernels/Rectangle.hpp>
#include <MeshDecomposer3D/kernels/SameRectangle.hpp>
#endif

// finders
#include "range/finders/BruteForce.hpp"
#include "range/finders/RangeTreeFinder.hpp"
#include "range/finders/OctTreeFinder.hpp"
#include "range/finders/KDTreeFinder.hpp"
#include "range/finders/GroupRangeTreeFinder.hpp"

#ifdef MADVORO_WITH_MPI
  // env agents
  #include <MeshDecomposer3D/environment/hilbert/DistributedOctEnvAgent.hpp>
  #include <MeshDecomposer3D/environment/hilbert/HilbertTreeEnvAgent.hpp>
#endif // MADVORO_WITH_MPI

#ifdef MADVORO_WITH_MPI
  #include <MeshDecomposer3D/points_manager/HilbertPointsManager.hpp>
  #define INITIAL_SENDRECV_TAG 1105
#endif 

#ifdef TIMING
#define VORONOI_TIMING(name) (auto name = std::chrono::high_resolution_clock::now();)
#define VORONOI_REPORT_TIMING(caption,start,end) ({auto duration = std::chrono::duration<double>(end - start).count(); if(rank == 0) std::cout << "Time for " << caption << ": " << duration << " seconds" << std::endl;})
#else
#define VORONOI_TIMING(name) (void)0;
#define VORONOI_REPORT_TIMING(caption,start,end) (void)0;
#endif // TIMING

#ifndef START_TIMER_PREEMPTIVE
#define START_TIMER_PREEMPTIVE(msg) (void)0
#endif

namespace MadVoro {

using namespace MadVoro::fallback;

template <typename PointT>
inline typename PointT::coord_type WrapPeriodicScalar(typename PointT::coord_type x, typename PointT::coord_type lo, typename PointT::coord_type hi)
{
    typename PointT::coord_type L = hi - lo;
    typename PointT::coord_type r = std::fmod(x - lo, L);
    if(r < typename PointT::coord_type(0))
    {
        r += L;
    }
    typename PointT::coord_type y = lo + r;
    if(not (y < hi))
    {
        y = lo;
    }
    return y;
}

typedef std::array<std::size_t, 4> b_array_4;
typedef std::array<std::size_t, 3> b_array_3;

//! \brief A three dimensional voronoi tessellation
template <typename PointT>
class Voronoi3D
{
public:
  using Face_T = Face3D<PointT>;
  using Point_T = PointT;
  using IVec = IndexedVector<PointT>;

private:
  PointT ll_, ur_;
  std::size_t Norg_;

  std::set<int> set_temp_;
  std::stack<int> stack_temp_;

  //  void CalcCellCMVolume(std::size_t index);

  vector<PointT> CreateBoundaryPoints(const std::vector<PointT> &points, vector<std::pair<std::size_t, std::size_t> > const& to_duplicate, vector<vector<size_t> > &past_duplicate);

  void InitialBoxBuild(std::vector<Face3D<PointT>> &box, std::vector<PointT> &normals);

  BuildDomain<PointT> Domain(void);

  // where a build writes: the mesh's own arrays
  DistributedBuild<PointT> Output(void);

  // what the mesh keeps of a build besides its arrays: which ghosts copy the rank's points, the periodic images' translations, the trees
  void AfterBuild(const DistributedBuild<PointT> &built);

  // the ghosts' centroids and volumes (from the ghost records), after a build
  void SyncGhostValues(void);

  PointT InferPeriodicTranslation(const PointT &p) const;

  void ShiftPeriodicImageCMs(void);
  
  #ifdef MADVORO_WITH_MPI
    std::vector<PointT> PrepareToBuildParallel(const std::vector<PointT> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing, bool suppressExchange);

    void FilterRealGhostPoints();
    void InitialExchange(const std::vector<PointT> &points, std::vector<int> &sentProc, std::vector<std::vector<size_t>> &sentPoints, const MPI_Comm &comm = MPI_COMM_WORLD);

  #endif // MADVORO_WITH_MPI

  // the built points, 4 unused points, then the ghosts
  std::vector<PointT> meshPoints_;
  vector<PointT> tetra_centers_; // the vertices
  // Voronoi Data
  //vector<vector<std::size_t> > FacesInCell_;
  vector<face_vec > FacesInCell_;
  std::vector<point_vec > PointsInFace_; // Right hand with regard to first neighbor
  //vector<vector<std::size_t> > PointsInFace_; // Right hand with regard to first neighbor
  vector<std::pair<std::size_t, std::size_t> > FaceNeighbors_;
  vector<PointT> all_CM;
  vector<PointT> CM_, Face_CM_; // center of masses
  vector<double> volume_; // volumes of each one of the tetrahedra
  vector<double> area_; // surface area of each one of the tetrahedra
  DistributedBuilder<PointT> builder_{[](bool onHost) -> std::unique_ptr<LocalBuilder<PointT>> { return std::make_unique<DelaunayLocalBuilder<PointT>>(onHost); }};
  
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

  Voronoi3D();

public:
  Voronoi3D(Voronoi3D const &other);

private:
  std::array<PointT, 4> temp_points_;
  std::array<PointT, 5> temp_points2_;
  std::vector<Face3D<PointT>> box_faces_;
  std::array<bool, 3> periodic_ = {false, false, false};
  std::vector<PointT> periodicImageTranslation_;
  mutable std::vector<std::size_t> periodicImagePhysical_;

  std::shared_ptr<OctTree<IVec>> myPointsTree;
  std::shared_ptr<OctTree<IVec>> allMyPointsTree;
  // Thread subdomains (0 or 1 builds as one domain)
  size_t threadSubdomainCount_ = 0;
  bool reportedUnusedSubdomains_ = false;

  void ReportUnusedThreadSubdomains(void);

  // builds the prepared cells; knownGhosts and knownMirrors: every ghost is known (other ranks' points, then mirrors)
  void BuildCells(bool otherRanks, const std::vector<PointT> *knownGhosts = nullptr, const std::vector<PointT> *knownMirrors = nullptr);

  std::vector<PointT> PreparePartialBuild(const std::vector<PointT> &allPoints, const std::vector<size_t> &indicesToBuild);
  #ifdef MADVORO_WITH_MPI
    std::shared_ptr<PointsManager<PointT, VoronoiPayload<PointT>>> pointsManager;
    std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>> indexingToSave = std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>>();
  #endif // MADVORO_WITH_MPI

  std::vector<PointT> allMyPoints;
  std::vector<double> allPointsWeights;
  std::vector<double> radiuses;

  AllPointsMap indicesInAllMyPoints; // the indices of the points in `meshPoints_`, in the list of all points

  template<typename T>
  void SyncPartialBuildData(std::vector<T> &partialBuildData, std::vector<T> &allBuildData) const;

  size_t buildGeneration_ = 0;

public:

  #ifdef MADVORO_WITH_MPI
    const std::vector<double> &GetPointsBuildWeights() const;
    
    const std::shared_ptr<EnvironmentAgent<PointT>> GetEnvironmentAgent() const;
    
    void SetKernel(const std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>> &indexing = std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>>());
    
    std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>> GetKernel() const;
    
    void SetBox(PointT const &ll, PointT const &ur, const std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>> &newIndexing);
  #endif // MADVORO_WITH_MPI

  #ifdef MADVORO_WITH_MPI
    vector<int> &GetSentProcs(void);

    vector<vector<size_t>> &GetSentPoints(void);

    vector<size_t> &GetSelfIndex(void);
  #endif // MADVORO_WITH_MPI

  vector<PointT> &GetAllFaceCM(void);

  /*! \brief Calculates the centre of mass of a face
    \param index Face index
    \return Face centre of mass
   */
  const PointT &FaceCM(std::size_t index) const;

  /*! \brief class constructor
    \param ll Lower left
    \param ur Upper right
   */
  Voronoi3D(const PointT &ll, const PointT &ur);

  Voronoi3D(const std::vector<Face3D<PointT>> &box_faces);

  void BuildInitialize(size_t num_points);

  void Build(vector<PointT> const& points)
  {
    std::vector<size_t> indicesToBuild(points.size());
    std::iota(indicesToBuild.begin(), indicesToBuild.end(), 0);
    this->BuildPartially(points, indicesToBuild);
  }

  size_t GetBuildGeneration(void) const { return buildGeneration_; }

  /** \brief Builds from `count` thread subdomains (0 or 1 builds as one domain). */
  void SetThreadSubdomains(size_t count) { threadSubdomainCount_ = count; }

  size_t GetThreadSubdomains(void) const { return threadSubdomainCount_; }

  //! \brief How cells are built: from a Delaunay tetrahedralization, or each on its own by clipping (in Kokkos' default execution space, e.g. a GPU)
  enum class CellBuilder
  {
      Delaunay,
      Clipping
  };

  void SetCellBuilder(CellBuilder cellBuilder);

  CellBuilder GetCellBuilder(void) const { return cellBuilder_; }

private:
  CellBuilder cellBuilder_ = CellBuilder::Delaunay;

public:

  void ReleaseMemory(void);

  /** \brief Returns spare container capacity to the allocator (call after a rebuild). */
  void ShrinkToFit(void);

  void BuildPartially(const std::vector<PointT> &allPoints, const std::vector<size_t> &indicesToBuild);

  void SetPeriodicBoundaries(bool x, bool y, bool z);

  std::array<bool, 3> GetPeriodicBoundaries() const;

  bool HasPeriodicBoundaries() const;

  void WrapPeriodicPoint(PointT &p) const;

  PointT WrapPeriodicPoint(const PointT &p) const;

  bool IsPeriodicImage(size_t delaunayIndex) const;

  const PointT &GetPeriodicImageTranslation(size_t delaunayIndex) const;

  size_t ResolvePeriodicImageIndex(size_t meshIndex) const;

  bool IsPointInCell(const PointT &point, size_t cellIndex, bool verbose = false) const;

  /*! \brief Like IsPointInCell, but returns the neighbor across the first violated face.
    \param point The point to test
    \param cellIndex The cell to test against (must be < Norg)
    \return (true, cellIndex) if inside; (false, neighborIndex) if a face was violated
   */
  std::pair<bool, size_t> FindViolatedFaceNeighbor(const PointT &point, size_t cellIndex) const;

#ifdef MADVORO_WITH_MPI
  void PreparePoints(const std::vector<PointT> &points, const std::vector<size_t> &mask);

  inline void SetPointsManager(std::shared_ptr<PointsManager<PointT, VoronoiPayload<PointT>>> pointsManager){this->pointsManager = pointsManager;};

  inline std::shared_ptr<PointsManager<PointT, VoronoiPayload<PointT>>> GetPointsManager() const{return this->pointsManager;};

  void MockMesh(void);
  
  void SetLoadBalancer(std::shared_ptr<LoadBalancer<PointT>> loadBalancer);

  void PresetLoadBalancer(std::shared_ptr<LoadBalancer<PointT>> loadBalancer);
  
  void Rebalance(const std::vector<double> &weights);

  inline bool ShouldRebalance(const std::vector<double> &weights) const {return this->pointsManager->shouldRebalance(weights);}

  inline bool ShouldRebalance(void) const {return this->pointsManager->shouldRebalance();}
    
  inline std::shared_ptr<LoadBalancer<PointT>> GetLoadBalancer(void) {return this->pointsManager->getLoadBalancer();};

  inline const std::shared_ptr<LoadBalancer<PointT>> GetLoadBalancer(void) const {return this->pointsManager->getLoadBalancer();};

  std::vector<PointT> BuildPartiallyParallel(const std::vector<PointT> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing = false, bool suppressExchange = false);

  std::vector<PointT> BuildPartiallyParallel(const std::vector<PointT> &allPoints, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing = false, bool suppressExchange = false)
  {
    return this->BuildPartiallyParallel(allPoints, std::vector<double>(allPoints.size(), 1.0), indicesToBuild, suppressRebalancing, suppressExchange);
  }

  std::vector<PointT> BuildParallel(const std::vector<PointT> &points, const std::vector<double> &weights, bool suppressRebalancing = false, bool suppressExchange = false)
  {
    std::vector<size_t> indicesToBuild(points.size());
    std::iota(indicesToBuild.begin(), indicesToBuild.end(), 0);
    return this->BuildPartiallyParallel(points, weights, indicesToBuild, suppressRebalancing, suppressExchange);
  }

  std::vector<PointT> BuildParallel(const std::vector<PointT> &points, bool suppressRebalancing = false, bool suppressExchange = false)
  {
    return this->BuildParallel(points, std::vector<double>(points.size(), 1.0), suppressRebalancing, suppressExchange);
  }
  
  bool DidRebalance(void) const;

  bool CheckContinuityOfZone(void) const;
  
  bool PointInMyDomain(const PointT &point) const;

  int GetOwner(const PointT &point) const;

  void SetImbalanceTolerance(double tolerance);
#endif // MADVORO_WITH_MPI

  double GetMaxRadius(const size_t &index) const;

  double GetMinRadius(const size_t &index) const;

  size_t GetContainingCell(const PointT &point) const;

  std::size_t GetPointNo(void) const;

  /*! \brief Get mehs point position
    \param index Index
    \return Position of point
   */
  const PointT &GetMeshPoint(std::size_t index) const;

  /*! \brief Calculate face area
    \param index Face index
    \return Face area
   */
  double GetArea(std::size_t index) const;

  /*! \brief Get cell centre of mass
    \param index Point index
    \return Centre of mass
   */
  PointT const& GetCellCM(std::size_t index) const;

  std::size_t GetTotalFacesNumber(void) const;

  /*! \brief Get cell size
    \param index Point index
    \return Cell width
   */
  double GetWidth(std::size_t index) const;

  /*! \brief Get cell volume
    \param index Point index
    \return Cell volume
   */
  double GetVolume(std::size_t index) const;

  /*! \brief Get cell faces
    \param index Point index
    \return List of bounding faces
   */ 
  face_vec const& GetCellFaces(std::size_t index) const;
  
  vector<PointT>& accessMeshPoints(void);

  const vector<PointT>& getMeshPoints(void) const;

  const AllPointsMap &GetIndicesInAllPoints(void) const;

  const std::vector<PointT> &getAllPoints(void) const;

  std::vector<PointT> &getAllPoints(void);

  size_t GetAllPointsNo(void) const;

  /*! \brief Get neighbours
    \param index Point index
    \return List of indices of neighbouring points
   */
  vector<std::size_t> GetNeighbors(std::size_t index)const;

  /*! \brief Checs if a point is near a boundary
    \param index Point index
    \return True if point is near a boundary
   */
  bool NearBoundary(std::size_t index) const;

  /*! \brief Checks if a face is on the boundary
    \param index Face index
    \return True if the face is on the boundary
   */
  bool BoundaryFace(std::size_t index) const;

  #ifdef MADVORO_WITH_MPI
  vector<vector<std::size_t> >& GetDuplicatedPoints(void);

  vector<vector<std::size_t> >const& GetDuplicatedPoints(void)const;

    /*! \brief Get Duplicated processe
      \return List of duplicated points
    */
    vector<int> GetDuplicatedProcs(void)const;

  /*! \brief Get a list of parallel processes to which points have been sent
    \return List of process numbers
   */
  vector<int> GetSentProcs(void)const;

  /*! \brief List of point sent to parallel processes, partitioned by processor number
    \return List of list of indices
   */
  vector<vector<std::size_t> > const& GetSentPoints(void)const;

  /*! \brief Get indices of all real cells
    \return List of indices of all real cells
   */
  vector<std::size_t> const& GetSelfIndex(void) const;

  #endif // MADVORO_WITH_MPI
  std::size_t GetTotalPointNumber(void)const;

  vector<PointT> & GetAllCM(void);

  vector<PointT > GetAllCM(void)const;

  /*! \brief Get neighbours of neighbours
    \param result Result
    \param point Point index
   */
  void GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point)const;

  /*! \brief Calculate normal vector to face
    \param faceindex Index of face
    \return Vector normal to face
   */
  PointT Normal(std::size_t faceindex)const;

  /*! \brief Checks if a point is a ghost
    \param index Index
    \return bool if the point is a ghost
   */
  bool IsGhostPoint(std::size_t index)const;

  /*! \brief Calculate the face velocity
    \param index Face index
    \param v0 Velocity of one neighbour
    \param v1 Velocity of second neighbour
    \return Face velocity
   */
  PointT CalcFaceVelocity(std::size_t index, PointT const& v0, PointT const& v1)const;

  vector<PointT>& GetFacePoints(void);

  vector<double>& GetAllArea(void);

  vector<PointT>const& GetFacePoints(void) const;

  //! \brief The flat Delaunay mesh of the last build as one domain
  const Geometry::FlatMesh &GetFlatMesh(void) const;

  //! \brief The flat tessellation of the last build as one domain (the ghosts follow the cells, with no unused points between them)
  const Geometry::FlatVoronoi &GetFlatVoronoi(void) const;

  vector<face_vec >& GetAllCellFaces(void);

  vector<face_vec >const& GetAllCellFaces(void) const;

  /*! \brief Get the points in face
    \param index Face index
    \return Indices of points in face
   */
  point_vec const& GetPointsInFace(std::size_t index) const;

  /*! \brief Get the neighbours across a face
    \param face_index Index of face
    \return Indices of neighbour across face
   */
  const std::pair<std::size_t, std::size_t> &GetFaceNeighbors(std::size_t face_index) const;

  #ifdef MADVORO_WITH_MPI
    /*! \brief Get the indices of ghost points
      \return List of list of ghost points
    */
    vector<vector<std::size_t> > const& GetGhostIndeces(void) const;

    /*! \brief Get the indices of ghost points
      \return List of list of ghost points
    */
    vector<vector<std::size_t> >& GetGhostIndeces(void);
  #endif // MADVORO_WITH_MPI

  void GetNeighbors(size_t index, vector<size_t> &res)const;

  /*! \brief Get the positions of opposite corners on the bounding box
    \return Pair of points on opposite corners
   */
  std::pair<PointT, PointT> GetBoxCoordinates(void)const;

  /*! \brief Build tessellation without a box
    \param points Mesh generating points
    \param ghosts Ghost points
    \param toduplicate Indices of points to be duplicated
   */
  void BuildNoBox(vector<PointT> const& points, vector<vector<PointT> > const& ghosts,vector<size_t> toduplicate);

  vector<double>& GetAllVolumes(void);

  vector<double> GetAllVolumes(void)const;

  /*! \brief Get all face neighbours
    \return List of pairs of indices to neighbours
   */
  std::vector<std::pair<size_t, size_t>> &GetAllFaceNeighbors(void);

  /*! \brief Get all face neighbours
    \return List of pairs of indices to neighbours
   */
  const std::vector<std::pair<size_t, size_t>> &GetAllFaceNeighbors(void) const;

  /*! \brief List all points in face
    \return List of all points in face
   */
  vector<point_vec > & GetAllPointsInFace(void);

  vector<point_vec > const& GetAllPointsInFace(void) const;

  /*! \brief Get the number of points
    \return The number of points
   */
  size_t& GetPointNo(void);

  /*! \brief Check whether a point is inside the computational domain
    \param index Index of point
    \return bool True if point is inside
   */
  bool IsPointOutsideBox(size_t index)const;

  bool IsPointOutsideBox(const PointT &point) const;

  /*! \brief Adjust position of the boundary
    \param ll Lower left corner
    \param ur Upper right corner
  */
  void SetBox(PointT const& ll, PointT const& ur);

  const std::vector<Face3D<PointT>> &GetBoxFaces(void) const {return this->box_faces_;}

  std::vector<Face3D<PointT>> GetBoxFaces(void) {return this->box_faces_;}

  std::vector<Face3D<PointT>>& ModifyBoxFaces(void) {return this->box_faces_;}
};

#ifdef MADVORO_WITH_MPI
#include <mpi_utils/mpi_alltoall.hpp>
#endif // MADVORO_WITH_MPI

template <typename PointT>
point_vec_v_t<PointT> VectorValues(std::vector<PointT> const& v, point_vec const& index)
{
    point_vec_v_t<PointT> res(index.size());
    for (size_t i = 0; i < index.size(); ++i)
        res[i] = v[index[i]];
    return res;
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
}

namespace
{
    template <typename PointT>
    vector<Face3D<PointT>> BuildBox(PointT const &ll, PointT const &ur)
    {
        double dx = ur.x - ll.x;
        double dy = ur.y - ll.y;
        double dz = ur.z - ll.z;
        vector<Face3D<PointT>> res(6);
        vector<PointT> points;
        points.push_back(ll);
        points.push_back(ll + PointT(dx, 0, 0));
        points.push_back(ll + PointT(dx, dy, 0));
        points.push_back(ll + PointT(0, dy, 0));
        points.push_back(ll + PointT(0, 0, dz));
        points.push_back(ll + PointT(dx, 0, dz));
        points.push_back(ll + PointT(dx, dy, dz));
        points.push_back(ll + PointT(0, dy, dz));
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
    template <typename PointT>
    vector<PointT> GetBoxNormals(PointT const &ll, PointT const &ur, vector<Face3D<PointT>> const& box_faces_)
    {
        const bool need_build = box_faces_.empty();
        const vector<Face3D<PointT>> built_faces = need_build ? BuildBox(ll, ur) : vector<Face3D<PointT>>();
        const vector<Face3D<PointT>> &faces = need_build ? built_faces : box_faces_;
        vector<PointT> res(faces.size());
        size_t N = res.size();
        for (size_t i = 0; i < N; ++i)
        {
            CrossProduct(faces[i].vertices[2] - faces[i].vertices[0], faces[i].vertices[1] - faces[i].vertices[0], res[i]);
            res[i] *= 1.0 /abs(res[i]);
        }
        return res;
    }

    template <typename PointT>
    size_t BoxIndex(vector<PointT> const &fnormals, const PointT &normal)
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
}

template <typename PointT>
Voronoi3D<PointT>::Voronoi3D(std::vector<Face3D<PointT>> const& box_faces) : Voronoi3D<PointT>()
{
    this->box_faces_ = box_faces;
    size_t const Nfaces = box_faces.size();
    if(Nfaces < 4)
        throw MadVoro::Exception::MadVoroException("Zero face vector in Voronoi3D constructor");
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

template <typename PointT>
Voronoi3D<PointT>::Voronoi3D(PointT const &ll, PointT const &ur) : ll_(ll), ur_(ur), Norg_(0), set_temp_(std::set<int>()), stack_temp_(std::stack<int>()),
                                                              tetra_centers_(vector<PointT>()),
                                                              FacesInCell_(vector<face_vec>()),
                                                              PointsInFace_(vector<point_vec>()),
                                                              FaceNeighbors_(vector<std::pair<std::size_t, std::size_t>>()),
                                                              CM_(vector<PointT>()), Face_CM_(vector<PointT>()),
                                                              volume_(vector<double>()), area_(vector<double>()), 
                                                              #ifdef MADVORO_WITH_MPI
                                                                sentprocs_(vector<int>()), sentpoints_(vector<vector<std::size_t>>()),  duplicatedprocs_(vector<int>()), 
                                                                duplicated_points_(vector<vector<std::size_t>>()), Nghost_(vector<vector<std::size_t>>()), self_index_(vector<std::size_t>()), 
                                                              #endif // MADVORO_WITH_MPI
                                                              temp_points_(std::array<PointT, 4>()), temp_points2_(std::array<PointT, 5>()), box_faces_(std::vector<Face3D<PointT>>()),
                                                              #ifdef MADVORO_WITH_MPI
                                                                pointsManager(nullptr),
                                                                allMyPoints(), 
                                                              #endif // MADVORO_WITH_MPI
                                                              indicesInAllMyPoints()
{
    this->box_faces_ = BuildBox(this->ll_, this->ur_);
    #ifdef MADVORO_WITH_MPI
        // initialize points manager
        this->pointsManager = std::shared_ptr<HilbertPointsManager<PointT, VoronoiPayload<PointT>>>(new HilbertPointsManager<PointT, VoronoiPayload<PointT>>(this->ll_, this->ur_));
    #endif // MADVORO_WITH_MPI
}

template <typename PointT>
Voronoi3D<PointT>::Voronoi3D() : Voronoi3D<PointT>(PointT(), PointT())
{}

template <typename PointT>
void Voronoi3D<PointT>::SetCellBuilder(CellBuilder cellBuilder)
{
    this->cellBuilder_ = cellBuilder;
    if(cellBuilder == CellBuilder::Clipping)
    {
        this->builder_ = DistributedBuilder<PointT>([](bool onHost) -> std::unique_ptr<LocalBuilder<PointT>> { return std::make_unique<ClippingLocalBuilder<PointT>>(onHost); });
    }
    else
    {
        this->builder_ = DistributedBuilder<PointT>([](bool onHost) -> std::unique_ptr<LocalBuilder<PointT>> { return std::make_unique<DelaunayLocalBuilder<PointT>>(onHost); });
    }
}

template <typename PointT>
const Geometry::FlatMesh &Voronoi3D<PointT>::GetFlatMesh(void) const
{
    const DelaunayLocalBuilder<PointT> *builder = dynamic_cast<const DelaunayLocalBuilder<PointT> *>(this->builder_.OneDomainBuilder());
    if(builder == nullptr)
    {
        throw MadVoro::Exception::MadVoroException("Voronoi3D::GetFlatMesh: the last build was not a Delaunay build of one domain");
    }
    return builder->GetFlatMesh();
}

template <typename PointT>
const Geometry::FlatVoronoi &Voronoi3D<PointT>::GetFlatVoronoi(void) const
{
    const Geometry::FlatVoronoi *voronoi = this->builder_.OneDomainVoronoi();
    if(voronoi == nullptr)
    {
        throw MadVoro::Exception::MadVoroException("Voronoi3D::GetFlatVoronoi: the last build was not of one domain");
    }
    return *voronoi;
}

template <typename PointT>
vector<PointT> Voronoi3D<PointT>::CreateBoundaryPoints(const std::vector<PointT> &points, vector<std::pair<std::size_t, std::size_t>> const &to_duplicate, vector<vector<size_t>> &past_duplicate)
{
    size_t Ncheck = to_duplicate.size();
    vector<std::pair<std::size_t, std::size_t>> to_add;
    to_add.reserve(Ncheck);
    const bool need_build = box_faces_.empty();
    const vector<Face3D<PointT>> built_box = need_build ? BuildBox(ll_, ur_) : vector<Face3D<PointT>>();
    const vector<Face3D<PointT>> &faces = need_build ? built_box : box_faces_;
    vector<PointT> res;
    bool first_time = past_duplicate.empty();
    if (first_time)
        past_duplicate.resize(faces.size());
    for (std::size_t i = 0; i < Ncheck; ++i)
    {
        if (first_time || !std::binary_search(past_duplicate[to_duplicate[i].first].begin(),
                                                                                    past_duplicate[to_duplicate[i].first].end(), to_duplicate[i].second))
        {
            res.push_back(MirrorPoint(faces[to_duplicate[i].first], points[to_duplicate[i].second]));
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
template <typename PointT>
    vector<vector<std::size_t>> const &Voronoi3D<PointT>::GetGhostIndeces(void) const
    {
        return Nghost_;
    }
#endif // MADVORO_WITH_MPI

/**
 * gets a point index, and returns the maximal radius of the tetrahedra containing that point.
 * @param index the index of the point (within the points list)
*/
template <typename PointT>
double Voronoi3D<PointT>::GetMaxRadius(const size_t &index) const
{
    double res = 0;
    for(const size_t &face : this->FacesInCell_[index])
    {
        for(const size_t &vertex : this->PointsInFace_[face])
        {
            res = std::max(res, fastabs(this->tetra_centers_[vertex] - this->meshPoints_[index]));
        }
    }
    return res;
}

/**
 * gets a point index, and returns the minimal radius of the tetrahedra containing that point.
 * @param index the index of the point (within the points list)
*/
template <typename PointT>
double Voronoi3D<PointT>::GetMinRadius(const size_t &index) const
{
    double res = std::numeric_limits<double>::max();
    for(const size_t &face : this->FacesInCell_[index])
    {
        for(const size_t &vertex : this->PointsInFace_[face])
        {
            res = std::min(res, fastabs(this->tetra_centers_[vertex] - this->meshPoints_[index]));
        }
    }
    return res;
}

/**
 * if the initial box does not exist, builds its faces according to the leftmost and rightmost points.
 * If it does, does not build the faces again.
 * @return the normals to the faces
*/
template <typename PointT>
void Voronoi3D<PointT>::InitialBoxBuild(std::vector<Face3D<PointT>> &box, std::vector<PointT> &normals)
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

/** \brief Initializes internal data structures, for the voronoi build. */
template <typename PointT>
void Voronoi3D<PointT>::BuildInitialize(size_t num_points)
{
    ++buildGeneration_;
    this->periodicImagePhysical_.clear();
    tetra_centers_.clear();
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
template <typename PointT>
bool Voronoi3D<PointT>::PointInMyDomain(const PointT &point) const
{
    size_t containingCell = this->GetContainingCell(point);
    return this->IsPointInCell(point, containingCell);
}

template <typename PointT>
inline int Voronoi3D<PointT>::GetOwner(const PointT &point) const
{
    PointT q = point;
    this->WrapPeriodicPoint(q);
    return this->pointsManager->getEnvironmentAgent()->getOwner(q);
}

template <typename PointT>
inline void Voronoi3D<PointT>::SetImbalanceTolerance(double tolerance)
{
    this->pointsManager->setImbalanceTolerance(tolerance);
}

template <typename PointT>
void Voronoi3D<PointT>::FilterRealGhostPoints()
{
    this->real_duplicated_proc.clear();
    this->real_duplicated_points.clear();

    // std::vector<bool> isNecessaryRecvPoint(this->meshPoints_.size(), false);
    // for(size_t pointIdx = 0; pointIdx < this->Norg_; pointIdx++)
    // {
    //     for(const size_t &neighborIdx : this->GetNeighbors(pointIdx))
    //     {
    //         isNecessaryRecvPoint[neighborIdx] = true;
    //     }
    // }

    // auto ifRecvCopyLambda = [&isNecessaryRecvPoint](const size_t &ghostPointIdx){return isNecessaryRecvPoint[ghostPointIdx];};

    std::vector<size_t> neighbor_buf;
    std::vector<size_t> newSend;
    for(size_t i = 0; i < this->duplicatedprocs_.size(); i++)
    {
        int _rank = this->duplicatedprocs_[i];
        this->real_duplicated_proc.push_back(_rank);
        newSend.clear();
        // check for any original sent point, if it has neighbors that belong to rank `_rank`. If yes, the point is necessary to be sent
        for(const size_t &pointIdxInBuild : this->duplicated_points_[i])
        {
            if(pointIdxInBuild >= this->Norg_)
            {
                // point was not participating in the last built
                continue;
            }
            bool foundNeighbor = false;
            this->GetNeighbors(pointIdxInBuild, neighbor_buf);
            for(const size_t &neighborIdx : neighbor_buf)
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



// void Voronoi3D<PointT>::InitialExchange(const std::vector<PointT> &points, std::vector<int> &sentProc, std::vector<std::vector<size_t>> &sentPoints, const MPI_Comm &comm)
// {
//     const EnvironmentAgent *envAgent = this->pointsManager->getEnvironmentAgent();
//     bool supportsFurthestClosestRanks;
//     std::function<HilbertCurveEnvironmentAgent::DistancesVector(const PointT&)> getFurthestClosestRanks;

//     // check if has 'smartAgent' (an agent that can caluclate distances of ranks as well)            
//     const DistributedOctEnvironmentAgent *distribuedOctEnvAgent = dynamic_cast<const DistributedOctEnvironmentAgent*>(envAgent);
//     if(distribuedOctEnvAgent != nullptr)
//     {
//         supportsFurthestClosestRanks = true;
//         getFurthestClosestRanks = [distribuedOctEnvAgent](const PointT &point){return distribuedOctEnvAgent->getClosestFurthestPointsByRanks(point);};
//     }
//     const HilbertTreeEnvironmentAgent *hilbertTreeEnvAgent = dynamic_cast<const HilbertTreeEnvironmentAgent*>(envAgent);
//     if(hilbertTreeEnvAgent != nullptr)
//     {
//         supportsFurthestClosestRanks = true;
//         getFurthestClosestRanks = [hilbertTreeEnvAgent](const PointT &point){return hilbertTreeEnvAgent->getClosestFurthestPointsByRanks(point);};
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
//     std::vector<std::vector<PointT>> toSend;
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
//                 // its index there will be this->meshPoints_.size() + batchInfo.pointsFromRanks[_rank][i]
//                 this->Nghost_[dupRankIdx].push_back(this->meshPoints_.size() + insertedSoFar + i);
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

//     std::vector<PointT> extraPoints;
//     for(const _3DPoint &_point : almostExtraPoints)
//     {
//         extraPoints.emplace_back(PointT(_point.x, _point.y, _point.z));
//     }

//     this->del_.BuildExtra(extraPoints);

//     this->R_.resize(this->del_.tetras_.size());
//     std::fill(this->R_.begin(), this->R_.end(), RADIUS_UNINITIALIZED);
//     this->tetra_centers_.resize(this->R_.size());
//     this->bigtet_ = SetPointTetras();
// }



/** \brief Rebalances if needed and exchanges the points between the ranks. */
template <typename PointT>
std::vector<PointT> Voronoi3D<PointT>::PrepareToBuildParallel(const std::vector<PointT> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing, bool suppressExchange)
{
    if(this->radiuses.size() < allPoints.size())
    {
        this->radiuses.resize(allPoints.size(), RADIUS_UNINITIALIZED);
    }
    if(this->all_CM.size() < allPoints.size())
    {
        this->all_CM.resize(allPoints.size());
    }
        
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
    {
        this->allPointsWeights.resize(allPoints.size(), 1.0);
    }

    std::vector<VoronoiPayload<PointT>> payloads(allPoints.size());
    for(size_t i = 0; i < allPoints.size(); i++)
    {
        payloads[i] = VoronoiPayload<PointT>(this->radiuses[i], this->all_CM[i]);
    }

    // non-participating points travel too (they are neighbours of participating ones)
    PointsExchangeResult<PointT, VoronoiPayload<PointT>> exchangeResult = this->pointsManager->update(allPoints, allowExchange? allWeights : this->allPointsWeights, payloads, indicesToBuild, allowRebalance, allowExchange, ExchangeScope::AllPointsFlagIndices);

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

    std::vector<PointT> new_points;
    this->indicesInAllMyPoints = AllPointsMap();
    size_t numOfSelfPoints = this->self_index_.size();

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
template <typename PointT>
void Voronoi3D<PointT>::PreparePoints(const std::vector<PointT> &points, const std::vector<size_t> &mask)
{
    if(points.size() != mask.size())
    {
        MadVoro::Exception::MadVoroException eo("In Voronoi3D<PointT>::PreparePoints, mask size is not equal to the points size");
        eo.addEntry("Mask size", mask.size());
        eo.addEntry("Points size", points.size());
        throw eo;
    }
    size_t originalPointsNum = this->allMyPoints.size();
    size_t newPointsNum = points.size();

    std::vector<IVec> oldPoints;
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
        OctTree<IVec> oldPointsTree(this->ll_, this->ur_, oldPoints);
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

template <typename PointT>
bool Voronoi3D<PointT>::DidRebalance(void) const
{
    return this->pointsManager->didRebalance();
}

template <typename PointT>
std::vector<PointT> Voronoi3D<PointT>::BuildPartiallyParallel(const std::vector<PointT> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing, bool suppressExchange)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(const size_t &idx : indicesToBuild)
    {
        if(idx >= allPoints.size())
        {
            MadVoro::Exception::MadVoroException eo("BuildParatiallyParallel: Illegal point was given");
            eo.addEntry("Index", idx);
            eo.addEntry("Size", allPoints.size());
            throw eo;
        }
    }

    std::chrono::high_resolution_clock::time_point start, end;
    START_TIMER_PREEMPTIVE("Prepare to build");

    start = std::chrono::high_resolution_clock::now();

    std::vector<PointT> activePoints = this->PrepareToBuildParallel(allPoints, allWeights, indicesToBuild, suppressRebalancing, suppressExchange);
    end = std::chrono::high_resolution_clock::now();
    if(rank == 0)
    {
        std::cout << "Time for preparing (hilbert tree + load balancing): " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    {
        int myPts = static_cast<int>(activePoints.size());
        int maxPts = myPts, minPts = myPts;
        MPI_Allreduce(MPI_IN_PLACE, &maxPts, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &minPts, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        if (rank == 0)
        {
            std::cout << "Post-LB point distribution: max=" << maxPts << ", min=" << minPts << std::endl;
        }
    }

    this->BuildCells(size > 1);
    return activePoints;
}


inline boost::container::flat_map<size_t, std::pair<rank_t, size_t>> GetGhostInfo(const std::vector<rank_t> &duplicatedProcs, const std::vector<std::vector<size_t>> &duplicatedPoints, const std::vector<std::vector<size_t>> &Nghost, const boost::container::flat_map<size_t, std::pair<rank_t, size_t>> &whereNow)
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
            ghostsInfo.emplace(oldIndex, receivedFromRank[j]);
        }
    }
    return ghostsInfo;
}

inline boost::container::flat_map<size_t, std::pair<rank_t, size_t>> GetRemoteIndices(const std::vector<size_t> &selfIndex, const std::vector<rank_t> &sentProcs, const std::vector<std::vector<size_t>> &sentPoints)
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

struct MockMeshGhostAsk
{
    rank_t receiver = 0;
    size_t pointIndex = 0;
    double imageTx = 0.0;
    double imageTy = 0.0;
    double imageTz = 0.0;

    bool operator==(MockMeshGhostAsk const &other) const
    {
        return receiver == other.receiver and pointIndex == other.pointIndex and imageTx == other.imageTx and imageTy == other.imageTy and imageTz == other.imageTz;
    }

    bool operator<(MockMeshGhostAsk const &other) const
    {
        if(receiver != other.receiver)
        {
            return receiver < other.receiver;
        }
        if(pointIndex != other.pointIndex)
        {
            return pointIndex < other.pointIndex;
        }
        if(imageTx != other.imageTx)
        {
            return imageTx < other.imageTx;
        }
        if(imageTy != other.imageTy)
        {
            return imageTy < other.imageTy;
        }
        return imageTz < other.imageTz;
    }
};

template <typename PointT>
void Voronoi3D<PointT>::MockMesh(void)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const bool periodic = this->HasPeriodicBoundaries();
    // all ranks must agree (the exchanges are collective)
    int useSubdomains = (this->threadSubdomainCount_ > 0 and not periodic) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &useSubdomains, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    if(periodic)
    {
        this->ReportUnusedThreadSubdomains();
    }
    // MADVORO_LB_TIMING=1 prints the phase timings (collective)
    static const bool lbTiming = std::getenv("MADVORO_LB_TIMING") != nullptr;
    auto _prof = [&](const char *label, double localTime)
    {
        if(not lbTiming)
        {
            return;
        }
        struct { double val; int rank; } localIn{localTime, rank}, globalMax;
        MPI_Allreduce(&localIn, &globalMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        double globalSum = 0;
        MPI_Allreduce(&localTime, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if(rank == 0)
        {
            std::cout << label << ": max=" << globalMax.val << "s (rank " << globalMax.rank << "), avg=" << globalSum / size << "s" << std::endl;
        }
    };
    auto _t0 = std::chrono::high_resolution_clock::now();
    auto _t1 = _t0;
    auto _mark = _t0;
    auto _phase = [&](const char *label)
    {
        _prof(label, std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - _mark).count());
        _mark = std::chrono::high_resolution_clock::now();
    };
    
    boost::container::flat_map<size_t, std::pair<rank_t, size_t>> whereNow = GetRemoteIndices(this->self_index_, this->sentprocs_, this->sentpoints_); // previous point -> holder (rank + index)
    assert(whereNow.size() == this->Norg_);
    boost::container::flat_map<size_t, std::pair<rank_t, size_t>> ghostsInfo = GetGhostInfo(this->duplicatedprocs_, this->duplicated_points_, this->Nghost_, whereNow); // previous ghost -> holder (rank + index)
    _phase("MockMesh: RemoteIndices+GhostInfo");

    std::vector<PointT> delaunayPoints = this->meshPoints_;

    std::vector<PointT> new_points;
    this->indicesInAllMyPoints = AllPointsMap();

    size_t allPointsSize = this->allMyPoints.size();
    for(size_t pointIdx = 0; pointIdx < allPointsSize; pointIdx++)
    {
        this->indicesInAllMyPoints[new_points.size()] = pointIdx;
        new_points.push_back(this->allMyPoints[pointIdx]);
    }
    
    // now find what ghosts should be sent
    std::vector<std::vector<MockMeshGhostAsk>> askToSend(size); // [sender][(receiver, point, imageTranslation)...]
    std::vector<std::vector<PointT>> mirrorsToSend(size);
    std::vector<std::vector<PointT>> periodicExtrasToBuild(size);

    auto periodicTranslationAt = [&](size_t delaunayIndex) -> PointT
    {
        if(not periodic or delaunayIndex >= this->periodicImageTranslation_.size())
        {
            return PointT();
        }
        return this->periodicImageTranslation_[delaunayIndex];
    };

    auto addAskToSend = [&](rank_t sender, rank_t receiver, size_t pointIndex, PointT const &imageTranslation)
    {
        if(sender == receiver)
        {
            return;
        }
        askToSend[sender].push_back(MockMeshGhostAsk{receiver, pointIndex, imageTranslation.x, imageTranslation.y, imageTranslation.z});
    };

    auto addPeriodicGhostOrAsk = [&](rank_t sender, rank_t receiver, size_t pointIndex, size_t oldNeighborIndex)
    {
        PointT imageTranslation = periodicTranslationAt(oldNeighborIndex);
        if(periodic)
        {
            if(imageTranslation.x != typename PointT::coord_type(0) or imageTranslation.y != typename PointT::coord_type(0) or imageTranslation.z != typename PointT::coord_type(0))
            {
                if(sender == receiver)
                {
                    periodicExtrasToBuild[receiver].push_back(delaunayPoints[oldNeighborIndex]);
                    return;
                }
            }
            else if(sender == receiver)
            {
                return;
            }
        }
        addAskToSend(sender, receiver, pointIndex, imageTranslation);
    };

    size_t previousN = this->Norg_;

    // pre-image in the previous mesh, verified (the nearest point is arbitrary when the pre-image is not local)
    auto resolveOldPreImage = [&](size_t oldImageIndex) -> size_t
    {
        const size_t preImage = this->ResolvePeriodicImageIndex(oldImageIndex);
        if(preImage >= previousN)
        {
            return previousN;
        }
        PointT expected = delaunayPoints[oldImageIndex] - periodicTranslationAt(oldImageIndex);
        this->WrapPeriodicPoint(expected);
        if(fastabs(delaunayPoints[preImage] - expected) > 1e-8 * fastabs(this->ur_ - this->ll_))
        {
            return previousN;
        }
        return preImage;
    };

    std::vector<size_t> neighbor_buf;
    for(size_t i = 0; i < previousN; i++)
    {
        // looking for previous point `i`
        const auto &[newOwner, newIndex] = whereNow.at(i);
        this->GetNeighbors(i, neighbor_buf);
        for(size_t neighbor : neighbor_buf)
        {
            // looking for previous neighbor `neighbor`
            if(neighbor < previousN)
            {
                // original point and neighbor are both local
                const auto &[neighborNewOwner, neighborNewIndex] = whereNow.at(neighbor);
                addAskToSend(neighborNewOwner, newOwner, neighborNewIndex, periodicTranslationAt(neighbor));
            }
            else
            {
                assert(neighbor >= previousN + 4); // 4 for the big tetra
                // neighbor is either a ghost, a periodic image, or a mirror
                auto it = ghostsInfo.find(neighbor);
                if(it != ghostsInfo.end())
                {
                    // neighbor is a ghost
                    const auto &[neighborNewOwner, neighborNewIndex] = it->second;
                    addPeriodicGhostOrAsk(neighborNewOwner, newOwner, neighborNewIndex, neighbor);
                }
                else if(this->IsPeriodicImage(neighbor))
                {
                    // a periodic image follows its pre-image, not the cell it borders
                    const size_t preImage = resolveOldPreImage(neighbor);
                    if(preImage < previousN)
                    {
                        const auto &[preImageNewOwner, preImageNewIndex] = whereNow.at(preImage);
                        addPeriodicGhostOrAsk(preImageNewOwner, newOwner, preImageNewIndex, neighbor);
                    }
                    else
                    {
                        periodicExtrasToBuild[newOwner].push_back(delaunayPoints[neighbor]);
                    }
                }
                else
                {
                    // Mirror across a non-periodic box face.
                    mirrorsToSend[newOwner].push_back(delaunayPoints[neighbor]);
                }
            }
        }
    }

    // mirrorsToSend = std::vector<std::vector<PointT>>(size, allMirrors); // todo: remove
    
    _phase("MockMesh: CollectGhostAsks");
    _t0 = std::chrono::high_resolution_clock::now();
    for(std::vector<MockMeshGhostAsk> &asksForSender : askToSend)
    {
        if(asksForSender.size() > 1)
        {
            std::sort(asksForSender.begin(), asksForSender.end());
            asksForSender.erase(std::unique(asksForSender.begin(), asksForSender.end()), asksForSender.end());
        }
    }
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("BuildAskToSend", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<rank_t, std::vector<MockMeshGhostAsk>>> whatIshouldSend = MPI_Exchange_sparse_by_rank(askToSend, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 2);
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("AskToSendExchange", std::chrono::duration<double>(_t1 - _t0).count());
    _t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<PointT>> whatIShouldSendToRanks(size);

    std::vector<rank_t> newDuplicatedProcs;
    std::vector<std::vector<size_t>> newDuplicatedPoints;
    boost::container::flat_map<rank_t, boost::container::flat_set<MockMeshGhostAsk>> sentToProcessors;
    boost::container::flat_map<rank_t, size_t> procsToIndices;

    for(const auto &[_requestingRank, sendInfo] : whatIshouldSend)
    {
        for(const MockMeshGhostAsk &ask : sendInfo)
        {
            if(ask.receiver == rank)
            {
                continue;
            }
            auto it = procsToIndices.find(ask.receiver);
            size_t idx = newDuplicatedProcs.size();
            if(it == procsToIndices.end())
            {
                newDuplicatedProcs.push_back(ask.receiver);
                newDuplicatedPoints.emplace_back();
                procsToIndices[ask.receiver] = idx;
            }
            else
            {
                idx = it->second;
            }
            std::vector<size_t> &newDuplicatedPointsOfRank = newDuplicatedPoints[idx];

            boost::container::flat_set<MockMeshGhostAsk> &sentToRank = sentToProcessors[ask.receiver];
            if(sentToRank.insert(ask).second)
            {
                whatIShouldSendToRanks[ask.receiver].push_back(new_points[ask.pointIndex] + PointT(ask.imageTx, ask.imageTy, ask.imageTz));
                newDuplicatedPointsOfRank.push_back(ask.pointIndex);
            }
        }
    }
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("BuildSendToRanks", std::chrono::duration<double>(_t1 - _t0).count());

    _t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<rank_t, std::vector<PointT>>> receivedByRanks = MPI_Exchange_sparse_by_rank(whatIShouldSendToRanks, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 3);
    _t1 = std::chrono::high_resolution_clock::now();
    _prof("FinalPointsExchange", std::chrono::duration<double>(_t1 - _t0).count());

    // the received ghosts, and their mesh indices (from `firstIndex` on) in Nghost
    const auto receivedGhosts = [&](size_t firstIndex, std::vector<PointT> &ghosts, std::vector<std::vector<size_t>> &Nghost)
    {
        size_t currentIndex = firstIndex;
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
                MadVoro::Exception::MadVoroException eo("Voronoi3D<PointT>::MockMesh: received from an unexpected rank");
                eo.addEntry("My rank", rank);
                eo.addEntry("Received from rank", _rank);
                throw eo;
            }
            assert(idx != newDuplicatedProcs.size());
            std::vector<size_t> &NghostOfRank = Nghost[idx];

            for(const PointT &point : receivedFromRank)
            {
                ghosts.push_back(point);
                NghostOfRank.push_back(currentIndex);
                currentIndex++;
            }
        }
    };

    if(useSubdomains == 1)
    {
        // the rank already has every ghost it needs
        this->BuildInitialize(new_points.size());
        std::vector<PointT> ghosts;
        std::vector<std::vector<size_t>> newNghost(newDuplicatedProcs.size());
        receivedGhosts(this->Norg_ + 4, ghosts, newNghost);
        this->duplicatedprocs_ = std::move(newDuplicatedProcs);
        this->duplicated_points_ = std::move(newDuplicatedPoints);
        this->Nghost_ = std::move(newNghost);
        // the old mesh's mirrors (one subdomain needs no search with them)
        std::vector<PointT> mirrors;
        std::set<std::array<typename PointT::coord_type, 3>> addedMirrors;
        for(const auto &[_senderRank, incomingMirrors] : MPI_Exchange_sparse_by_rank(mirrorsToSend, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 1))
        {
            for(const PointT &mirror : incomingMirrors)
            {
                if(addedMirrors.insert({mirror.x, mirror.y, mirror.z}).second)
                {
                    mirrors.push_back(mirror);
                }
            }
        }
        this->BuildCells(false, &ghosts, &mirrors);
        _phase("MockMesh: thread subdomains");
        return;
    }

    _mark = std::chrono::high_resolution_clock::now();
    this->BuildInitialize(new_points.size());
    // the old mesh's periodic images, then its mirrors, from each rank (each once: a point can border several cells), then the ghosts other ranks send
    std::vector<std::vector<PointT>> extras;
    if(periodic)
    {
        std::set<std::array<typename PointT::coord_type, 3>> addedImages;
        for(const auto &[_senderRank, incomingExtras] : MPI_Exchange_sparse_by_rank(periodicExtrasToBuild, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 1))
        {
            std::vector<PointT> newImages;
            for(const PointT &image : incomingExtras)
            {
                if(addedImages.insert({image.x, image.y, image.z}).second)
                {
                    newImages.push_back(image);
                }
            }
            if(this->Norg_ > 0 and not newImages.empty())
            {
                extras.push_back(std::move(newImages));
            }
        }
    }
    std::set<std::array<typename PointT::coord_type, 3>> addedMirrors;
    for(const auto &[_senderRank, incomingMirrors] : MPI_Exchange_sparse_by_rank(mirrorsToSend, MPI_COMM_WORLD, MPI_EXCHANGE_SPARSE_TAG + 1))
    {
        std::vector<PointT> newMirrors;
        for(const PointT &mirror : incomingMirrors)
        {
            if(addedMirrors.insert({mirror.x, mirror.y, mirror.z}).second)
            {
                newMirrors.push_back(mirror);
            }
        }
        if(this->Norg_ > 0 and not newMirrors.empty())
        {
            extras.push_back(std::move(newMirrors));
        }
    }
    size_t ghostsIndex = this->Norg_ + 4;
    for(const std::vector<PointT> &extra : extras)
    {
        ghostsIndex += extra.size();
    }
    std::vector<PointT> buildExtra;
    std::vector<std::vector<size_t>> newNghost(newDuplicatedProcs.size());
    receivedGhosts(ghostsIndex, buildExtra, newNghost);
    extras.push_back(std::move(buildExtra));
    _phase("MockMesh: extras and ghosts");
    std::vector<const std::vector<PointT> *> ghosts;
    for(const std::vector<PointT> &extra : extras)
    {
        ghosts.push_back(&extra);
    }
    std::vector<size_t> cells(new_points.size());
    std::iota(cells.begin(), cells.end(), 0);
    DistributedBuild<PointT> built = this->Output();
    this->builder_.BuildKnown(this->Domain(), this->allMyPoints, cells, this->radiuses, new_points, ghosts, built);
    this->AfterBuild(built);
    this->duplicatedprocs_ = std::move(newDuplicatedProcs);
    this->duplicated_points_ = std::move(newDuplicatedPoints);
    this->Nghost_ = std::move(newNghost);
    this->SyncGhostValues();
    _phase("MockMesh: build");
}

template <typename PointT>
void Voronoi3D<PointT>::SetLoadBalancer(std::shared_ptr<LoadBalancer<PointT>> loadBalancer)
{
    assert(this->pointsManager != nullptr);
    this->pointsManager->setLoadBalancer(loadBalancer);

    std::vector<PointT> previousPoints = this->allMyPoints;
    std::vector<size_t> allIndices(previousPoints.size());
    std::iota(allIndices.begin(), allIndices.end(), 0);

    std::vector<VoronoiPayload<PointT>> payloads(previousPoints.size());
    for(size_t i = 0; i < previousPoints.size(); i++)
    {
        payloads[i] = VoronoiPayload<PointT>(this->radiuses[i], this->all_CM[i]);
    }

    const auto _lbStart = std::chrono::high_resolution_clock::now();
    PointsExchangeResult<PointT, VoronoiPayload<PointT>> exchangeResult = this->pointsManager->update(previousPoints, this->allPointsWeights, payloads, allIndices, false);
    if(std::getenv("MADVORO_LB_TIMING") != nullptr)
    {
        int _rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        double _local = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - _lbStart).count(), _max = 0;
        MPI_Allreduce(&_local, &_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if(_rank == 0)
        {
            std::cout << "SetLoadBalancer: points exchange max=" << _max << "s (" << previousPoints.size() << " points on rank 0 before)" << std::endl;
        }
    }

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

template <typename PointT>
void Voronoi3D<PointT>::PresetLoadBalancer(std::shared_ptr<LoadBalancer<PointT>> loadBalancer)
{
    assert(this->pointsManager != nullptr);
    this->pointsManager->setLoadBalancer(loadBalancer);
}

template <typename PointT>
void Voronoi3D<PointT>::Rebalance(const std::vector<double> &weights)
{
    if(this->pointsManager == nullptr)
    {
        MadVoro::Exception::MadVoroException eo("Voronoi3D<PointT>::Rebalance: pointsManager is nullptr");
        throw eo;
    }

    if(weights.size() != this->allMyPoints.size())
    {
        MadVoro::Exception::MadVoroException eo("Voronoi3D<PointT>::Rebalance: weights.size() != allMyPoints.size()");
        eo.addEntry("weights.size()", weights.size());
        eo.addEntry("allMyPoints.size()", this->allMyPoints.size());
        throw eo;
    }
    
    std::vector<PointT> previousPoints = this->allMyPoints;
    std::vector<size_t> allIndices(previousPoints.size());
    std::iota(allIndices.begin(), allIndices.end(), 0);
    std::vector<VoronoiPayload<PointT>> payloads(previousPoints.size());
    for(size_t i = 0; i < previousPoints.size(); i++)
        payloads[i] = VoronoiPayload<PointT>(this->radiuses[i], this->all_CM[i]);

    PointsExchangeResult<PointT, VoronoiPayload<PointT>> exchangeResult = this->pointsManager->update(previousPoints, weights, payloads, allIndices);

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

#endif // MADVORO_WITH_MPI

template <typename PointT>
void Voronoi3D<PointT>::SetPeriodicBoundaries(bool x, bool y, bool z)
{
    PointT L = this->ur_ - this->ll_;
    if(x and not (L.x > typename PointT::coord_type(0)))
    {
        throw MadVoro::Exception::InvalidArgumentException("SetPeriodicBoundaries: invalid x periodic axis length");
    }
    if(y and not (L.y > typename PointT::coord_type(0)))
    {
        throw MadVoro::Exception::InvalidArgumentException("SetPeriodicBoundaries: invalid y periodic axis length");
    }
    if(z and not (L.z > typename PointT::coord_type(0)))
    {
        throw MadVoro::Exception::InvalidArgumentException("SetPeriodicBoundaries: invalid z periodic axis length");
    }
    this->periodic_[0] = x;
    this->periodic_[1] = y;
    this->periodic_[2] = z;
}

template <typename PointT>
std::array<bool, 3> Voronoi3D<PointT>::GetPeriodicBoundaries() const
{
    return this->periodic_;
}

template <typename PointT>
bool Voronoi3D<PointT>::HasPeriodicBoundaries() const
{
    return this->periodic_[0] or this->periodic_[1] or this->periodic_[2];
}

template <typename PointT>
void Voronoi3D<PointT>::WrapPeriodicPoint(PointT &p) const
{
    if(this->periodic_[0])
    {
        p.x = WrapPeriodicScalar<PointT>(p.x, this->ll_.x, this->ur_.x);
    }
    if(this->periodic_[1])
    {
        p.y = WrapPeriodicScalar<PointT>(p.y, this->ll_.y, this->ur_.y);
    }
    if(this->periodic_[2])
    {
        p.z = WrapPeriodicScalar<PointT>(p.z, this->ll_.z, this->ur_.z);
    }
}

template <typename PointT>
PointT Voronoi3D<PointT>::WrapPeriodicPoint(const PointT &p) const
{
    PointT q = p;
    WrapPeriodicPoint(q);
    return q;
}

template <typename PointT>
bool Voronoi3D<PointT>::IsPeriodicImage(size_t delaunayIndex) const
{
    if(delaunayIndex >= this->periodicImageTranslation_.size())
    {
        return false;
    }
    const PointT &translation = this->periodicImageTranslation_[delaunayIndex];
    return translation.x != typename PointT::coord_type(0) or translation.y != typename PointT::coord_type(0) or translation.z != typename PointT::coord_type(0);
}

template <typename PointT>
const PointT &Voronoi3D<PointT>::GetPeriodicImageTranslation(size_t delaunayIndex) const
{
    static const PointT zero;
    if(delaunayIndex >= this->periodicImageTranslation_.size())
    {
        return zero;
    }
    return this->periodicImageTranslation_[delaunayIndex];
}

template <typename PointT>
size_t Voronoi3D<PointT>::ResolvePeriodicImageIndex(size_t meshIndex) const
{
    if(meshIndex < this->Norg_)
    {
        return meshIndex;
    }
    if(meshIndex >= this->meshPoints_.size() or not this->IsPeriodicImage(meshIndex) or this->myPointsTree == nullptr)
    {
        return meshIndex;
    }
    if(this->periodicImagePhysical_.size() != this->meshPoints_.size())
    {
        this->periodicImagePhysical_.assign(this->meshPoints_.size(), std::numeric_limits<std::size_t>::max());
    }
    if(meshIndex < this->periodicImagePhysical_.size())
    {
        std::size_t &cached = this->periodicImagePhysical_[meshIndex];
        if(cached == std::numeric_limits<std::size_t>::max())
        {
            PointT unwrapped = this->meshPoints_[meshIndex];
            const PointT &translation = this->periodicImageTranslation_[meshIndex];
            unwrapped.x -= translation.x;
            unwrapped.y -= translation.y;
            unwrapped.z -= translation.z;
            PointT wrapped = unwrapped;
            this->WrapPeriodicPoint(wrapped);
            cached = this->myPointsTree->closestPoint(wrapped).getIndex();
        }
        if(cached < this->Norg_)
        {
            return cached;
        }
    }
    return meshIndex;
}

template <typename PointT>
PointT Voronoi3D<PointT>::InferPeriodicTranslation(const PointT &p) const
{
    return PeriodicTranslationOf(p, this->ll_, this->ur_, this->periodic_);
}

template <typename PointT>
BuildDomain<PointT> Voronoi3D<PointT>::Domain(void)
{
    BuildDomain<PointT> domain;
    domain.ll = this->ll_;
    domain.ur = this->ur_;
    domain.boxFaces = this->box_faces_;
    this->InitialBoxBuild(domain.box, domain.normals);
    domain.periodic = this->periodic_;
    return domain;
}

template <typename PointT>
DistributedBuild<PointT> Voronoi3D<PointT>::Output(void)
{
    #ifdef MADVORO_WITH_MPI
        return DistributedBuild<PointT>{this->meshPoints_, this->FacesInCell_, this->FaceNeighbors_, this->PointsInFace_, this->area_, this->Face_CM_, this->tetra_centers_, this->volume_, this->CM_, this->duplicatedprocs_, this->duplicated_points_, this->Nghost_};
    #else // MADVORO_WITH_MPI
        return DistributedBuild<PointT>{this->meshPoints_, this->FacesInCell_, this->FaceNeighbors_, this->PointsInFace_, this->area_, this->Face_CM_, this->tetra_centers_, this->volume_, this->CM_};
    #endif // MADVORO_WITH_MPI
}

template <typename PointT>
void Voronoi3D<PointT>::AfterBuild(const DistributedBuild<PointT> &built)
{
    for(const auto &[meshIndex, pointIdx] : built.ghostRankIndex)
    {
        this->indicesInAllMyPoints[meshIndex] = pointIdx;
    }
    // a periodic image lies outside the box, by the translation from its point
    this->periodicImageTranslation_.assign(this->meshPoints_.size(), PointT());
    if(this->HasPeriodicBoundaries())
    {
        for(size_t i = this->Norg_ + 4; i < this->meshPoints_.size(); ++i)
        {
            this->periodicImageTranslation_[i] = this->InferPeriodicTranslation(this->meshPoints_[i]);
        }
    }
    this->allMyPointsTree = this->builder_.RankTree();
    this->myPointsTree = this->builder_.BuiltTree();
}

template <typename PointT>
void Voronoi3D<PointT>::SyncGhostValues(void)
{
    this->SyncPartialBuildData(this->CM_, this->all_CM);
    this->ShiftPeriodicImageCMs();
    #ifdef MADVORO_WITH_MPI
        // without a points manager (e.g. a thread subdomain) the build talked to no other rank
        if(this->pointsManager == nullptr)
        {
            return;
        }
        // the ghosts' volumes, those of their cells in the ranks that own them
        const std::vector<rank_t> &correspondents = this->GetDuplicatedProcs();
        const std::vector<std::vector<size_t>> &indices = this->GetDuplicatedPoints();
        // indexed by all the rank's points, not by cells (differs in a partial build)
        std::vector<double> allVolumes(this->GetAllPointsNo(), 0.0);
        for(size_t k = 0; k < this->Norg_; ++k)
        {
            allVolumes[this->indicesInAllMyPoints.at(k)] = this->volume_[k];
        }
        std::vector<std::vector<double>> exchange = MPI_exchange_data_indexed(correspondents, allVolumes, indices);
        const std::vector<std::vector<size_t>> &ghost_indices = this->GetGhostIndeces();
        this->volume_.resize(this->GetTotalPointNumber(), 0.0);
        for(size_t i = 0; i < correspondents.size(); ++i)
        {
            const std::vector<double> &data = exchange[i];
            for(size_t j = 0; j < data.size(); ++j)
            {
                this->volume_.at(ghost_indices.at(i).at(j)) = data[j];
            }
        }
    #endif // MADVORO_WITH_MPI
}

/** \brief Shifts the centroids of the periodic images by their translations. */
template <typename PointT>
void Voronoi3D<PointT>::ShiftPeriodicImageCMs(void)
{
    if(this->HasPeriodicBoundaries())
    {
        // images that SyncPartialBuildData did not fill take the centroid of their pre-image
        std::vector<bool> filledBySync(this->CM_.size(), false);
        const AllPointsMap &indicesInAllPoints = this->GetIndicesInAllPoints();
        for(size_t i = this->Norg_; i < this->CM_.size(); ++i)
        {
            filledBySync[i] = (indicesInAllPoints.find(i) != indicesInAllPoints.cend());
        }
#ifdef MADVORO_WITH_MPI
        for(const std::vector<size_t> &ghostIndices : this->Nghost_)
        {
            for(size_t ghostIndex : ghostIndices)
            {
                if(ghostIndex < filledBySync.size())
                {
                    filledBySync[ghostIndex] = true;
                }
            }
        }
#endif // MADVORO_WITH_MPI
        for(size_t i = this->Norg_; i < this->CM_.size(); ++i)
        {
            if(not this->IsPeriodicImage(i))
            {
                continue;
            }
            if(not filledBySync[i])
            {
                const size_t physical = this->ResolvePeriodicImageIndex(i);
                if(physical < this->Norg_)
                {
                    this->CM_[i] = this->CM_[physical];
                }
            }
            this->CM_[i] = this->CM_[i] + this->periodicImageTranslation_[i];
        }
    }
}

#ifdef MADVORO_WITH_MPI
template <typename PointT>
    vector<vector<std::size_t>> &Voronoi3D<PointT>::GetGhostIndeces(void)
    {
        return Nghost_;
    }
#endif // MADVORO_WITH_MPI

template <typename PointT>
std::pair<PointT, PointT> Voronoi3D<PointT>::GetBoxCoordinates(void) const
{
    return std::pair<PointT, PointT>(ll_, ur_);
}

template <typename PointT>
void Voronoi3D<PointT>::ShrinkToFit(void)
{
    meshPoints_.shrink_to_fit();
    tetra_centers_.shrink_to_fit();
    FacesInCell_.shrink_to_fit();
    PointsInFace_.shrink_to_fit();
    FaceNeighbors_.shrink_to_fit();
    all_CM.shrink_to_fit();
    CM_.shrink_to_fit();
    Face_CM_.shrink_to_fit();
    volume_.shrink_to_fit();
    area_.shrink_to_fit();
    allMyPoints.shrink_to_fit();
    allPointsWeights.shrink_to_fit();
    radiuses.shrink_to_fit();
    periodicImageTranslation_.shrink_to_fit();
    periodicImagePhysical_.shrink_to_fit();
    builder_.ShrinkToFit();
#ifdef MADVORO_WITH_MPI
    sentprocs_.shrink_to_fit();
    for(auto &v : sentpoints_)
    {
        v.shrink_to_fit();
    }
    sentpoints_.shrink_to_fit();
    duplicatedprocs_.shrink_to_fit();
    for(auto &v : duplicated_points_)
    {
        v.shrink_to_fit();
    }
    duplicated_points_.shrink_to_fit();
    real_duplicated_proc.shrink_to_fit();
    for(auto &v : real_duplicated_points)
    {
        v.shrink_to_fit();
    }
    real_duplicated_points.shrink_to_fit();
    for(auto &v : Nghost_)
    {
        v.shrink_to_fit();
    }
    Nghost_.shrink_to_fit();
    self_index_.shrink_to_fit();
#endif
}

template <typename PointT>
void Voronoi3D<PointT>::ReleaseMemory(void)
{
    ContainerOps::release_container_memory(meshPoints_);
    ContainerOps::release_container_memory(tetra_centers_);
    ContainerOps::release_container_memory(FacesInCell_);
    ContainerOps::release_container_memory(PointsInFace_);
    ContainerOps::release_container_memory(FaceNeighbors_);
    ContainerOps::release_container_memory(all_CM);
    ContainerOps::release_container_memory(CM_);
    ContainerOps::release_container_memory(Face_CM_);
    ContainerOps::release_container_memory(volume_);
    ContainerOps::release_container_memory(area_);
    ContainerOps::release_container_memory(box_faces_);
    ContainerOps::release_container_memory(allMyPoints);
    ContainerOps::release_container_memory(allPointsWeights);
    ContainerOps::release_container_memory(radiuses);
    ContainerOps::release_container_memory(indicesInAllMyPoints);
    ContainerOps::release_container_memory(periodicImageTranslation_);
    builder_.ReleaseMemory();
#ifdef MADVORO_WITH_MPI
    ContainerOps::release_container_memory(sentprocs_);
    ContainerOps::release_container_memory(sentpoints_);
    ContainerOps::release_container_memory(duplicatedprocs_);
    ContainerOps::release_container_memory(duplicated_points_);
    ContainerOps::release_container_memory(real_duplicated_proc);
    ContainerOps::release_container_memory(real_duplicated_points);
    ContainerOps::release_container_memory(Nghost_);
    ContainerOps::release_container_memory(self_index_);
#endif
}

template <typename PointT>
void Voronoi3D<PointT>::BuildNoBox(vector<PointT> const &points, vector<vector<PointT>> const &ghosts, vector<size_t> toduplicate)
{
    assert(points.size() > 0);
    this->BuildInitialize(points.size());
    // the ghosts, then the points of `toduplicate` mirrored across the box faces
    std::vector<std::vector<PointT>> mirrors;
    vector<vector<size_t>> past_duplicates;
    vector<std::pair<size_t, size_t>> duplicate(6);
    for(const size_t &point : toduplicate)
    {
        for(size_t i = 0; i < 6; ++i)
        {
            duplicate[i] = std::pair<size_t, size_t>(i, point);
        }
        mirrors.push_back(this->CreateBoundaryPoints(points, duplicate, past_duplicates));
    }
    std::vector<const std::vector<PointT> *> extras;
    for(const std::vector<PointT> &group : ghosts)
    {
        extras.push_back(&group);
    }
    for(const std::vector<PointT> &group : mirrors)
    {
        extras.push_back(&group);
    }
    std::vector<double> noRadii;
    DistributedBuild<PointT> built = this->Output();
    this->builder_.BuildKnown(this->Domain(), std::vector<PointT>(), std::vector<size_t>(), noRadii, points, extras, built);
    this->AfterBuild(built);
}

template <typename PointT>
void Voronoi3D<PointT>::BuildPartially(const std::vector<PointT> &allPoints, const std::vector<size_t> &indicesToBuild)
{
    this->PreparePartialBuild(allPoints, indicesToBuild);
    this->BuildCells(false);
}

/** \brief Returns the points of a partial build and clears the tessellation. */
template <typename PointT>
std::vector<PointT> Voronoi3D<PointT>::PreparePartialBuild(const std::vector<PointT> &allPoints, const std::vector<size_t> &indicesToBuild)
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
    std::vector<PointT> activePoints = MadVoro::ContainerOps::VectorValues(allPoints, indicesToBuild);
    
    size_t pointsCounter = 0;
    this->indicesInAllMyPoints = AllPointsMap();
    for(const size_t &pointIdx : indicesToBuild)
    {
        this->indicesInAllMyPoints[pointsCounter] = pointIdx;
        pointsCounter++;
    }

    this->BuildInitialize(activePoints.size());
    return activePoints;
}

/** \brief Says once that a periodic load balance switch rebuilds without thread subdomains. */
template <typename PointT>
void Voronoi3D<PointT>::ReportUnusedThreadSubdomains(void)
{
    if(this->threadSubdomainCount_ <= 1 or this->reportedUnusedSubdomains_)
    {
        return;
    }
    this->reportedUnusedSubdomains_ = true;
    int rank = 0;
    #ifdef MADVORO_WITH_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif // MADVORO_WITH_MPI
    if(rank == 0)
    {
        std::cout << "MadVoro: " << this->threadSubdomainCount_ << " thread subdomains were asked for, but a load balance " << "switch of a periodic box rebuilds each rank as one domain" << std::endl;
    }
}

/** \brief Builds the prepared cells (expects the build prepared up to BuildInitialize). */
template <typename PointT>
void Voronoi3D<PointT>::BuildCells(bool otherRanks, const std::vector<PointT> *knownGhosts, const std::vector<PointT> *knownMirrors)
{
    const size_t pointsNum = this->Norg_;
    std::vector<size_t> cells(pointsNum);
    for(size_t k = 0; k < pointsNum; ++k)
    {
        cells[k] = this->indicesInAllMyPoints.at(k);
    }
    const BuildDomain<PointT> domain = this->Domain();
    if(knownGhosts != nullptr and knownMirrors != nullptr and DistributedBuilder<PointT>::SubdomainsFor(this->threadSubdomainCount_, pointsNum) == 1)
    {
        // nothing to search: the ghosts (at their recorded indices) before the mirrors
        const std::vector<PointT> points = MadVoro::ContainerOps::VectorValues(this->allMyPoints, cells);
        DistributedBuild<PointT> built = this->Output();
        this->builder_.BuildKnown(domain, this->allMyPoints, cells, this->radiuses, points, {knownGhosts, knownMirrors}, built);
        this->AfterBuild(built);
        this->SyncGhostValues();
        return;
    }
    #ifdef MADVORO_WITH_MPI
        std::optional<typename DistributedBuilder<PointT>::ExchangeRequest> request;
        if(otherRanks)
        {
            const std::shared_ptr<EnvironmentAgent<PointT>> envAgent = this->pointsManager->getEnvironmentAgent();
            if(not envAgent)
            {
                throw MadVoro::Exception::MadVoroException("BuildCells: environment agent is null");
            }
            // the last build's real ghosts go first
            request.emplace(typename DistributedBuilder<PointT>::ExchangeRequest{envAgent, this->real_duplicated_proc, this->real_duplicated_points, MPI_COMM_WORLD});
        }
        DistributedBuild<PointT> built = this->Output();
        this->builder_.Build(domain, this->threadSubdomainCount_, this->allMyPoints, cells, this->radiuses, knownGhosts, request ? &*request : nullptr, built);
        this->AfterBuild(built);
    #else // MADVORO_WITH_MPI
        DistributedBuild<PointT> built = this->Output();
        this->builder_.Build(domain, this->threadSubdomainCount_, this->allMyPoints, cells, this->radiuses, knownGhosts, built);
        this->AfterBuild(built);
    #endif // MADVORO_WITH_MPI
    this->SyncGhostValues();
}

/*
template <typename PointT>
void Voronoi3D<PointT>::CalcCellCMVolume(std::size_t index)
{
    volume_[index] = 0;
    CM_[index] = PointT();
    std::size_t Nfaces = FacesInCell_[index].size();
    std::array<PointT, 4> tetra;
    tetra[3] = meshPoints_[index];
    PointT vtemp;
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


template <typename PointT>
size_t Voronoi3D<PointT>::GetContainingCell(const PointT &point) const
{
    if(not this->HasPeriodicBoundaries())
    {
        return this->myPointsTree->closestPoint(point).getIndex();
    }
    PointT q0 = point;
    this->WrapPeriodicPoint(q0);
    PointT L = this->ur_ - this->ll_;
    size_t bestIdx = this->myPointsTree->closestPoint(q0).getIndex();
    typename PointT::coord_type bestDist2 = abs(this->GetMeshPoint(bestIdx) - q0);
    bestDist2 = bestDist2 * bestDist2;
    std::array<int, 3> shiftRanges[3] = {{0}, {0}, {0}};
    if(this->periodic_[0])
    {
        shiftRanges[0] = {-1, 0, 1};
    }
    if(this->periodic_[1])
    {
        shiftRanges[1] = {-1, 0, 1};
    }
    if(this->periodic_[2])
    {
        shiftRanges[2] = {-1, 0, 1};
    }
    for(int sx : shiftRanges[0])
    {
        for(int sy : shiftRanges[1])
        {
            for(int sz : shiftRanges[2])
            {
                if(sx == 0 and sy == 0 and sz == 0)
                {
                    continue;
                }
                PointT qi(q0.x + L.x * sx, q0.y + L.y * sy, q0.z + L.z * sz);
                size_t idx = this->myPointsTree->closestPoint(qi).getIndex();
                typename PointT::coord_type dist2 = abs(this->GetMeshPoint(idx) - qi);
                dist2 = dist2 * dist2;
                if(dist2 < bestDist2)
                {
                    bestDist2 = dist2;
                    bestIdx = idx;
                }
            }
        }
    }
    return bestIdx;
}

template <typename PointT>
std::size_t Voronoi3D<PointT>::GetPointNo(void) const
{
    return Norg_;
}

template <typename PointT>
const PointT &Voronoi3D<PointT>::GetMeshPoint(std::size_t index) const
{
    return meshPoints_[index];
}

template <typename PointT>
double Voronoi3D<PointT>::GetArea(std::size_t index) const
{
    return area_[index];
}

template <typename PointT>
PointT const &Voronoi3D<PointT>::GetCellCM(std::size_t index) const
{
    return this->CM_[index];
}

template <typename PointT>
std::size_t Voronoi3D<PointT>::GetTotalFacesNumber(void) const
{
    return FaceNeighbors_.size();
}

template <typename PointT>
double Voronoi3D<PointT>::GetWidth(std::size_t index) const
{
    return std::pow(3 * volume_[index] * 0.25 / M_PI, 0.3333333333);
}

template <typename PointT>
double Voronoi3D<PointT>::GetVolume(std::size_t index) const
{
    return volume_[index];
}

template <typename PointT>
face_vec const &Voronoi3D<PointT>::GetCellFaces(std::size_t index) const
{
    return FacesInCell_[index];
}

template <typename PointT>
vector<PointT> &Voronoi3D<PointT>::accessMeshPoints(void)
{
    return meshPoints_;
}

template <typename PointT>
const vector<PointT> &Voronoi3D<PointT>::getMeshPoints(void) const
{
    return meshPoints_;
}

template <typename PointT>
const AllPointsMap &Voronoi3D<PointT>::GetIndicesInAllPoints(void) const
{
    return this->indicesInAllMyPoints;
}

template <typename PointT>
const std::vector<PointT> &Voronoi3D<PointT>::getAllPoints(void) const
{
    return this->allMyPoints;
}

template <typename PointT>
std::vector<PointT> &Voronoi3D<PointT>::getAllPoints(void)
{
    return this->allMyPoints;
}

template <typename PointT>
vector<std::size_t> Voronoi3D<PointT>::GetNeighbors(std::size_t index) const
{
    const size_t N = FacesInCell_[index].size();
    vector<size_t> res(N);
    for (size_t i = 0; i < N; ++i)
    {
        size_t face = FacesInCell_[index][i];
        res[i] = FaceNeighbors_[face].first == index ? FaceNeighbors_[face].second : FaceNeighbors_[face].first;
    }
    return res;
}

template <typename PointT>
void Voronoi3D<PointT>::GetNeighbors(size_t index, vector<size_t> &res) const
{
    std::size_t N = FacesInCell_[index].size();
    res.resize(N);
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#pragma ivdep
#endif
    for (std::size_t i = 0; i < N; ++i)
    {
        std::size_t face = FacesInCell_[index][i];
        res[i] = FaceNeighbors_[face].first == index ? FaceNeighbors_[face].second : FaceNeighbors_[face].first;
    }
}


template <typename PointT>
Voronoi3D<PointT>::Voronoi3D(Voronoi3D<PointT> const &other) : ll_(other.ll_), ur_(other.ur_), Norg_(other.Norg_),
                                                set_temp_(other.set_temp_), stack_temp_(other.stack_temp_), meshPoints_(other.meshPoints_),
                                                tetra_centers_(other.tetra_centers_), FacesInCell_(other.FacesInCell_), PointsInFace_(other.PointsInFace_),
                                                FaceNeighbors_(other.FaceNeighbors_), CM_(other.CM_), Face_CM_(other.Face_CM_), volume_(other.volume_), area_(other.area_), builder_(other.builder_), cellBuilder_(other.cellBuilder_),
                                                #ifdef MADVORO_WITH_MPI
                                                    sentprocs_(other.sentprocs_), sentpoints_(other.sentpoints_), duplicatedprocs_(other.duplicatedprocs_), duplicated_points_(other.duplicated_points_),
                                                    Nghost_(other.Nghost_), self_index_(other.self_index_),
                                                #endif // MADVORO_WITH_MPI
                                                temp_points_(std::array<PointT, 4>()), temp_points2_(std::array<PointT, 5>()), box_faces_(other.box_faces_),
                                                periodic_(other.periodic_), periodicImageTranslation_(other.periodicImageTranslation_),
                                                threadSubdomainCount_(other.threadSubdomainCount_),
                                                #ifdef MADVORO_WITH_MPI
                                                    pointsManager(other.pointsManager->clone()), indexingToSave(other.indexingToSave),
                                                    radiuses(other.radiuses), allMyPoints(other.allMyPoints), allPointsWeights(other.allPointsWeights),
                                                #endif // MADVORO_WITH_MPI
                                                indicesInAllMyPoints(other.indicesInAllMyPoints)
                                                {}

template <typename PointT>
bool Voronoi3D<PointT>::NearBoundary(std::size_t index) const
{
    std::size_t N = FacesInCell_[index].size();
    for (std::size_t i = 0; i < N; ++i)
    {
        if (BoundaryFace(FacesInCell_[index][i]))
            return true;
    }
    return false;
}

template <typename PointT>
bool Voronoi3D<PointT>::IsPointInCell(const PointT &point, size_t cellIndex, bool verbose) const
{
    if(cellIndex >= this->Norg_)
    {
        MadVoro::Exception::MadVoroException eo("Voronoi3D<PointT>::IsPointInCell: cell index out of range");
        eo.addEntry("point", point);
        eo.addEntry("cellIndex", cellIndex);
        eo.addEntry("Norg", this->Norg_);
        throw eo;
    }
    std::optional<MadVoro::Exception::MadVoroException> verboseInfo;
    if(verbose)
    {
        verboseInfo.emplace("Voronoi3D<PointT>::IsPointInCell verbose");
        verboseInfo->addEntry("point", point);
        verboseInfo->addEntry("cellIndex", cellIndex);
        verboseInfo->addEntry("cellCenter", this->meshPoints_[cellIndex]);

        PointT cellLL(std::numeric_limits<double>::max());
        PointT cellUR(std::numeric_limits<double>::lowest());
        const std::vector<PointT> &vertices = this->GetFacePoints();
        for(size_t faceIdx : this->FacesInCell_[cellIndex])
        {
            for(size_t vertIdx : this->GetAllPointsInFace()[faceIdx])
            {
                const PointT &v = vertices[vertIdx];
                cellLL.x = std::min(cellLL.x, v.x);
                cellLL.y = std::min(cellLL.y, v.y);
                cellLL.z = std::min(cellLL.z, v.z);
                cellUR.x = std::max(cellUR.x, v.x);
                cellUR.y = std::max(cellUR.y, v.y);
                cellUR.z = std::max(cellUR.z, v.z);
            }
        }
        BoundingBox<PointT> cellBox(cellLL, cellUR);  
        verboseInfo->addEntry("cellBBox", cellBox);
        verboseInfo->addEntry("is point in bounding box of cell?", cellBox.contains(point)? "yes": "no");
    }

    double width = this->GetWidth(cellIndex);
    PointT q = point;
    this->WrapPeriodicPoint(q);
    if(this->HasPeriodicBoundaries())
    {
        PointT L = this->ur_ - this->ll_;
        typename PointT::coord_type deltas[3] = {q.x - this->meshPoints_[cellIndex].x, q.y - this->meshPoints_[cellIndex].y, q.z - this->meshPoints_[cellIndex].z};
        typename PointT::coord_type lengths[3] = {L.x, L.y, L.z};
        typename PointT::coord_type *qCoords[3] = {&q.x, &q.y, &q.z};
        for(size_t d = 0; d < 3; ++d)
        {
            if(this->periodic_[d])
            {
                *qCoords[d] -= lengths[d] * std::nearbyint(deltas[d] / lengths[d]);
            }
        }
    }
    for(size_t faceIdx : this->FacesInCell_[cellIndex])
    {
        const PointT &p1 = this->meshPoints_[this->GetFaceNeighbors(faceIdx).first];
        const PointT &p2 = this->meshPoints_[this->GetFaceNeighbors(faceIdx).second];
        PointT normal = (this->GetFaceNeighbors(faceIdx).second == cellIndex)? p2 - p1 : p1 - p2;
        PointT p = (p1 + p2) * 0.5;
        double dot = ScalarProd(normal, q - p);
        if(verboseInfo)
        {
            size_t neighbor = (this->GetFaceNeighbors(faceIdx).first == cellIndex)? this->GetFaceNeighbors(faceIdx).second : this->GetFaceNeighbors(faceIdx).first;
            verboseInfo->addEntry("face " + std::to_string(faceIdx) + " neighbor", neighbor);
            std::string neighborOwner = "???";
            if(neighbor < this->Norg_)
            {
                neighborOwner = "me";
            }
#ifdef MADVORO_WITH_MPI
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
#endif // MADVORO_WITH_MPI
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

template <typename PointT>
std::pair<bool, size_t> Voronoi3D<PointT>::FindViolatedFaceNeighbor(const PointT &point, size_t cellIndex) const
{
    double width = this->GetWidth(cellIndex);
    PointT q = point;
    this->WrapPeriodicPoint(q);
    if(this->HasPeriodicBoundaries())
    {
        PointT L = this->ur_ - this->ll_;
        typename PointT::coord_type deltas[3] = {q.x - this->meshPoints_[cellIndex].x, q.y - this->meshPoints_[cellIndex].y, q.z - this->meshPoints_[cellIndex].z};
        typename PointT::coord_type lengths[3] = {L.x, L.y, L.z};
        typename PointT::coord_type *qCoords[3] = {&q.x, &q.y, &q.z};
        for(size_t d = 0; d < 3; ++d)
        {
            if(this->periodic_[d])
            {
                *qCoords[d] -= lengths[d] * std::nearbyint(deltas[d] / lengths[d]);
            }
        }
    }
    for(size_t faceIdx : this->FacesInCell_[cellIndex])
    {
        const auto &[n1, n2] = this->GetFaceNeighbors(faceIdx);
        PointT normal = (n2 == cellIndex) ? meshPoints_[n2] - meshPoints_[n1] : meshPoints_[n1] - meshPoints_[n2];
        PointT mid = (meshPoints_[n1] + meshPoints_[n2]) * 0.5;
        double dot = ScalarProd(normal, q - mid);
        if(dot < -1e-12 * width)
        {
            size_t neighbor = (n1 == cellIndex) ? n2 : n1;
            return {false, neighbor};
        }
    }
    return {true, cellIndex};
}

template <typename PointT>
bool Voronoi3D<PointT>::IsPointOutsideBox(size_t index) const
{
    return this->IsPointOutsideBox(this->meshPoints_[index]);
}

template <typename PointT>
bool Voronoi3D<PointT>::IsPointOutsideBox(const PointT &point) const
{
    return this->box_faces_.empty() ? not PointInDomain(this->ll_, this->ur_, point) : not PointInPoly(this->box_faces_, point);
}

template <typename PointT>
bool Voronoi3D<PointT>::BoundaryFace(std::size_t index) const
{
    size_t other = std::max(FaceNeighbors_[index].first, FaceNeighbors_[index].second);
    if(other >= Norg_ and IsPeriodicImage(other))
    {
        return false;
    }
    if (FaceNeighbors_[index].first >= Norg_ or FaceNeighbors_[index].second >= Norg_)
    {
#ifdef MADVORO_WITH_MPI
        if(box_faces_.empty())
        {
            if (PointInDomain(ll_, ur_, meshPoints_[other]))
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        else if(PointInPoly(box_faces_, meshPoints_[other]))
        {
            return false;
        }
        else
#endif
        {
            return true;
        }
    }
    else
    {
        return false;
    }
}

#ifdef MADVORO_WITH_MPI
template <typename PointT>
bool Voronoi3D<PointT>::CheckContinuityOfZone(void) const
{
    std::vector<bool> reached(this->Norg_, false);
    reached[0] = true;
    std::stack<std::size_t> stack;
    stack.push(0);
    std::vector<std::size_t> neighbor_buf;
    while(not stack.empty())
    {
        size_t cell = stack.top();
        stack.pop();
        this->GetNeighbors(cell, neighbor_buf);
        for(const auto &neighbor : neighbor_buf)
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
#endif // MADVORO_WITH_MPI

#ifdef MADVORO_WITH_MPI
template <typename PointT>
    vector<vector<std::size_t>> &Voronoi3D<PointT>::GetDuplicatedPoints(void)
    {
        return duplicated_points_;
    }

template <typename PointT>
    vector<vector<std::size_t>> const &Voronoi3D<PointT>::GetDuplicatedPoints(void) const
    {
        return duplicated_points_;
    }
#endif // MADVORO_WITH_MPI

template <typename PointT>
std::size_t Voronoi3D<PointT>::GetTotalPointNumber(void) const
{
    return meshPoints_.size();
}

template <typename PointT>
vector<PointT> &Voronoi3D<PointT>::GetAllCM(void)
{
    return CM_;
}

template <typename PointT>
vector<PointT> Voronoi3D<PointT>::GetAllCM(void) const
{
    return CM_;
}

template <typename PointT>
void Voronoi3D<PointT>::GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point) const
{
    result.clear();
    result.reserve(70);
    vector<std::size_t> neigh;
    GetNeighbors(point, neigh);
    result = neigh;
    std::size_t N = neigh.size();
    std::sort(neigh.begin(), neigh.end());
    vector<std::size_t> temp;
    for (std::size_t i = 0; i < N; ++i)
    {
        if (neigh[i] < Norg_)
        {
            GetNeighbors(neigh[i], temp);
            result.insert(result.end(), temp.begin(), temp.end());
        }
    }
    std::sort(result.begin(), result.end());
    MadVoro::ContainerOps::unique_inplace(result);
    MadVoro::ContainerOps::RemoveList_inplace(result, neigh);
    MadVoro::ContainerOps::RemoveVal(result, point);
}

template <typename PointT>
vector<boost::container::small_vector<size_t, 8>> &Voronoi3D<PointT>::GetAllPointsInFace(void)
{
    return PointsInFace_;
}

template <typename PointT>
vector<boost::container::small_vector<size_t, 8>> const& Voronoi3D<PointT>::GetAllPointsInFace(void)const
{
    return PointsInFace_;
}

template <typename PointT>
size_t &Voronoi3D<PointT>::GetPointNo(void)
{
    return Norg_;
}

template <typename PointT>
size_t Voronoi3D<PointT>::GetAllPointsNo(void) const
{
    return this->allMyPoints.size();
}

template <typename PointT>
std::vector<std::pair<size_t, size_t>> &Voronoi3D<PointT>::GetAllFaceNeighbors(void)
{
    return FaceNeighbors_;
}

template <typename PointT>
const std::vector<std::pair<size_t, size_t>> &Voronoi3D<PointT>::GetAllFaceNeighbors(void) const
{
    return FaceNeighbors_;
}

template <typename PointT>
vector<double> &Voronoi3D<PointT>::GetAllVolumes(void)
{
    return volume_;
}

template <typename PointT>
vector<double> Voronoi3D<PointT>::GetAllVolumes(void) const
{
    return volume_;
}

template <typename PointT>
PointT Voronoi3D<PointT>::Normal(std::size_t faceindex) const
{
    return meshPoints_[FaceNeighbors_[faceindex].second] - meshPoints_[FaceNeighbors_[faceindex].first];
}

template <typename PointT>
bool Voronoi3D<PointT>::IsGhostPoint(std::size_t index) const
{
    return index >= Norg_;
}

template <typename PointT>
const PointT &Voronoi3D<PointT>::FaceCM(std::size_t index) const
{
    return Face_CM_[index];
}

template <typename PointT>
PointT Voronoi3D<PointT>::CalcFaceVelocity(std::size_t index, PointT const &v0, PointT const &v1) const
{
    std::size_t p0 = FaceNeighbors_[index].first;
    std::size_t p1 = FaceNeighbors_[index].second;
    PointT r0 = GetMeshPoint(p0);
    PointT r1 = GetMeshPoint(p1);
    PointT r_diff = r1 - r0;
    double abs_r_diff = ScalarProd(r_diff, r_diff);

    PointT f = FaceCM(index);
    r1 += r0;
    r1 *= 0.5;
    f -= r1;
    PointT delta_w = ScalarProd((v0 - v1), f) * r_diff / abs_r_diff;
#ifdef RICH_DEBUG
    double dw_abs = fastabs(delta_w);
#endif // RICH_DEBUG
    PointT w = (v0 + v1) * 0.5;
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
        MadVoro::Exception::MadVoroException eo("Bad Face velocity");
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

template <typename PointT>
vector<double> &Voronoi3D<PointT>::GetAllArea(void)
{
    return area_;
}

template <typename PointT>
vector<PointT> &Voronoi3D<PointT>::GetAllFaceCM(void)
{
    return Face_CM_;
}

template <typename PointT>
vector<face_vec> &Voronoi3D<PointT>::GetAllCellFaces(void)
{
    return FacesInCell_;
}

template <typename PointT>
vector<face_vec> const& Voronoi3D<PointT>::GetAllCellFaces(void)const
{
    return FacesInCell_;
}

template <typename PointT>
vector<PointT> &Voronoi3D<PointT>::GetFacePoints(void)
{
    return tetra_centers_;
}

template <typename PointT>
vector<PointT> const &Voronoi3D<PointT>::GetFacePoints(void) const
{
    return tetra_centers_;
}

template <typename PointT>
point_vec const &Voronoi3D<PointT>::GetPointsInFace(std::size_t index) const
{
    return PointsInFace_[index];
}

template <typename PointT>
const std::pair<std::size_t, std::size_t> &Voronoi3D<PointT>::GetFaceNeighbors(std::size_t face_index) const
{
    return FaceNeighbors_[face_index];
}

#ifdef MADVORO_WITH_MPI
template <typename PointT>
    vector<int> Voronoi3D<PointT>::GetDuplicatedProcs(void) const
    {
        return duplicatedprocs_;
    }

template <typename PointT>
    vector<int> Voronoi3D<PointT>::GetSentProcs(void) const
    {
        return sentprocs_;
    }

template <typename PointT>
    vector<vector<std::size_t>> const &Voronoi3D<PointT>::GetSentPoints(void) const
    {
        return sentpoints_;
    }

template <typename PointT>
    vector<std::size_t> const &Voronoi3D<PointT>::GetSelfIndex(void) const
    {
        return self_index_;
    }

template <typename PointT>
    vector<int> &Voronoi3D<PointT>::GetSentProcs(void)
    {
        return sentprocs_;
    }

template <typename PointT>
    vector<vector<std::size_t>> &Voronoi3D<PointT>::GetSentPoints(void)
    {
        return sentpoints_;
    }

template <typename PointT>
    vector<std::size_t> &Voronoi3D<PointT>::GetSelfIndex(void)
    {
        return self_index_;
    }
#endif // MADVORO_WITH_MPI

#ifdef MADVORO_WITH_MPI
template <typename PointT>
void Voronoi3D<PointT>::SetKernel(const std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>> &indexing)
{
    MPI_Barrier(MPI_COMM_WORLD); // everyone should set the kernel
    this->indexingToSave = indexing;
    HilbertPointsManager<PointT, VoronoiPayload<PointT>> *hilbertPointsManager = dynamic_cast<HilbertPointsManager<PointT, VoronoiPayload<PointT>> *>(this->pointsManager.get());
    if(hilbertPointsManager == nullptr)
    {
        // points manager is not a 'HilbertPointsManager', or was not initialized yet
        return;
    }
    hilbertPointsManager->setIndexing(indexing);
}

template <typename PointT>
std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>> Voronoi3D<PointT>::GetKernel() const
{
    return this->indexingToSave;
}
#endif // MADVORO_WITH_MPI

template <typename PointT>
void Voronoi3D<PointT>::SetBox(const PointT &ll, const PointT &ur)
{
    this->ll_ = ll;
    this->ur_ = ur;
    this->box_faces_ = BuildBox(this->ll_, this->ur_);
    #ifdef MADVORO_WITH_MPI
        this->pointsManager = std::make_shared<HilbertPointsManager<PointT, VoronoiPayload<PointT>>>(this->ll_, this->ur_, MPI_COMM_WORLD);
        // this->radiuses.clear();
    #endif // MADVORO_WITH_MPI
}

#ifdef MADVORO_WITH_MPI
template <typename PointT>
void Voronoi3D<PointT>::SetBox(PointT const &ll, PointT const &ur, const std::shared_ptr<const Kernelization3D::IndexingKernel3D<PointT>> &newIndexing)
{
    this->SetBox(ll, ur);
    this->SetKernel(newIndexing);
}

template <typename PointT>
const std::vector<double> &Voronoi3D<PointT>::GetPointsBuildWeights() const
{
    return this->allPointsWeights;
}

template <typename PointT>
const std::shared_ptr<EnvironmentAgent<PointT>> Voronoi3D<PointT>::GetEnvironmentAgent() const
{
    return this->pointsManager->getEnvironmentAgent();
}

#endif // MADVORO_WITH_MPI


template <typename PointT>
template<typename T>
void Voronoi3D<PointT>::SyncPartialBuildData(std::vector<T> &partialBuildData, std::vector<T> &allBuildData) const
{
  size_t Norg = this->GetPointNo();
  if(partialBuildData.size() < Norg)
  {
    MadVoro::Exception::MadVoroException eo("Voronoi3D<PointT>::SyncPartialBuildData: Partial build data has lower size than the number of points");
    eo.addEntry("Partial build data size", partialBuildData.size());
    eo.addEntry("Number of points", Norg);
    throw eo;
  }
  const AllPointsMap &indicesInAllMyPoints = this->GetIndicesInAllPoints();

  allBuildData.resize(this->GetAllPointsNo());

  for(size_t i = 0; i < Norg; i++)
  {
      size_t pointIdx = indicesInAllMyPoints.at(i);
      allBuildData[pointIdx] = partialBuildData[i];
  }

  size_t sizeOfMeshPoints = this->getMeshPoints().size();
  partialBuildData.resize(sizeOfMeshPoints);
  for(size_t i = Norg; i < sizeOfMeshPoints; i++)
  {
      bool pointIsMine = (indicesInAllMyPoints.find(i) != indicesInAllMyPoints.cend());
      if(pointIsMine)
      {
          size_t pointIdx = indicesInAllMyPoints.at(i);
          partialBuildData[i] = allBuildData[pointIdx];
      }
  }

  #ifdef MADVORO_WITH_MPI
      // without a points manager (e.g. a thread subdomain) the build talked to no other rank
      if(this->pointsManager == nullptr)
      {
          return;
      }
      std::vector<std::vector<T>> incoming = MPI_exchange_data_indexed(this->GetDuplicatedProcs(), allBuildData, this->GetDuplicatedPoints());
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

} // namespace MadVoro

#endif // VORONOI3D_HPP
