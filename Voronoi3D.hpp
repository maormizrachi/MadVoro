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
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/container/small_vector.hpp>
#include "utils/debug/SmartTimer.hpp"
#include <omp.h>

#ifdef USE_VCL_VECTORIZATION
  #include <vectorclass.h>
#endif // USE_VCL_VECTORIZATION

#ifdef RICH_MPI
  #include <mpi.h>
#endif // RICH_MPI

#include "VoronoiPayload.hpp"
#include "3D/tessellation/Tessellation3D.hpp"
#include "delaunay/Delaunay3D.hpp"
#include "geometry/Intersections.hpp"
#include "3D/elementary/Mat33.hpp"
#include "utils/Predicates3D.hpp"
#include <MeshDecomposer3D/hilbert/HilbertOrder3D.hpp>
#include "3D/range/SmallRangeAgent.hpp"
#include "3D/range/BigRangeAgent.hpp"
#include "misc/utils.hpp"
#include "misc/io3D.hpp"

#ifdef RICH_MPI
#include <mpi_utils/mpi_commands.hpp>
#endif

// finders
#include "3D/range/finders/BruteForce.hpp"
#include "3D/range/finders/RangeTree.hpp"
#include "3D/range/finders/OctTree.hpp"
#include "3D/range/finders/KDTree.hpp"
#include "3D/range/finders/GroupRangeTree.hpp"

#ifdef RICH_MPI
  // env agents
  #include <MeshDecomposer3D/environment/hilbert/DistributedOctEnvAgent.hpp>
  #include <MeshDecomposer3D/environment/hilbert/HilbertTreeEnvAgent.hpp>
#endif // RICH_MPI

#ifdef RICH_MPI
  #include <MeshDecomposer3D/points_manager/HilbertPointsManager.hpp>
  #define INITIAL_SENDRECV_TAG 1105
#endif 

#define LARGE_POINTS_SHRINK_RADIUS_RATIO 0.95
#define RANGE_MAX_POINTS_TO_GET 15 // 15
#define RADIUSES_GROWING_FACTOR 1.1
#define RADIUS_UNINITIALIZED -1

typedef std::array<std::size_t, 4> b_array_4;
typedef std::array<std::size_t, 3> b_array_3;

//! \brief A three dimensional voronoi tessellation
class Voronoi3D : public Tessellation3D
{
public:
  using Face_T = Face;
  using Point_T = Vector3D;

private:
  Vector3D ll_, ur_;
  std::size_t Norg_, bigtet_;

  std::set<int> set_temp_;
  std::stack<int> stack_temp_;

  /*! \brief Finds intersections for a single cell
    \param box Bounding faces
    \param point Point index
    \param sphere Bounding sphere
    \param intersecting_faces Result
    \param Rtemp TBA
    \param vtemp TBA
   */
  void FindIntersectionsSingle(vector<Face> const& box, std::size_t point, Sphere<Vector3D> &sphere,
			       vector<size_t> &intersecting_faces, std::vector<double> &Rtemp, std::vector<Vector3D> &vtemp);

  std::size_t GetFirstPointToCheck(void)const;

  /*! \brief Get point to check
    \param point Point index
    \param checked Mast of checked points
    \param res Result
   */
  void GetPointToCheck(std::size_t point, vector<unsigned char> const& checked, vector<std::size_t> &res);
  
  void CalcRigidCM(std::size_t face_index);

  void GetTetraCM(std::array<Vector3D, 4> const& points, Vector3D &CM)const;

  double GetTetraVolume(std::array<Vector3D, 4> const& points)const;

  //  void CalcCellCMVolume(std::size_t index);

  double GetRadius(const size_t &index) const;

  void CalcAllCM(void);

  vector<std::pair<std::size_t, std::size_t> > SerialFindIntersections(bool first_run);

  vector<std::pair<std::size_t, std::size_t> > SerialFirstIntersections(void);

  double CalcTetraRadiusCenterHiPrecision(const size_t &index) const;

  double CalcTetraRadiusCenter(const size_t &index) const;

  vector<Vector3D> CreateBoundaryPoints(vector<std::pair<std::size_t, std::size_t> > const& to_duplicate,
					vector<vector<size_t> > &past_duplicate);
  void BuildVoronoi(std::vector<size_t> const& order);

  size_t SetPointTetras(void);

  void InitialBoxBuild(std::vector<Face> &box, std::vector<Vector3D> &normals);

  void BringSelfGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                              BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                              boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                              boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints,
                              RangeFinder::_set<size_t> &selfIgnorePoints);
    
  #ifdef RICH_MPI
    void BringGhostPointsToBuild(const MPI_Comm &comm);
  #else
    void BringGhostPointsToBuild();
  #endif // RICH_MPI

  std::pair<std::vector<SmallRangeQueryData>, std::vector<BigRangeQueryData>> CreateBatches(boost::container::flat_set<size_t> &smallPoints, boost::container::flat_set<size_t> &largePoints, const boost::container::flat_map<size_t, size_t> &firstLargeIteration, std::vector<double> &currentRadiuses, size_t iterations);

  std::pair<boost::container::flat_set<size_t>, boost::container::flat_set<size_t>>
    DetermineNextIterationPoints(size_t iterations,
                                  boost::container::flat_map<size_t, size_t> &firstLargeIteration,
                                  std::vector<double> &currentRadiuses,
                                  const boost::container::flat_map<size_t, size_t> &resultOfSmallPoints,
                                  const boost::container::flat_map<size_t, size_t> &resultOfBigPoints
                                );

  void UpdateRadiuses(const std::vector<Vector3D> &points);

  void UpdateCMs(void);
  
  void UpdateRangeFinder(void);

  void UpdatePointsTree(const std::vector<Vector3D> &activePoints);
  
  #ifdef RICH_MPI
    std::vector<Vector3D> PrepareToBuildParallel(const std::vector<Vector3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing, bool suppressExchange);
    void FilterRealGhostPoints();
    void UpdateDuplicatedPoints(const std::vector<int> &sentProc, const std::vector<std::vector<size_t>> &sentPoints);
    void EnsureSymmetry(const std::vector<int> &sentProc, const std::vector<std::vector<int>> &recvProcLists);
    std::tuple<std::vector<Vector3D>, std::vector<std::vector<size_t>>, std::vector<std::vector<size_t>>> InitialGhostPointsExchange(const MPI_Comm &comm = MPI_COMM_WORLD) const;
    void InitialExchange(const std::vector<Vector3D> &points, std::vector<int> &sentProc, std::vector<std::vector<size_t>> &sentPoints, const MPI_Comm &comm = MPI_COMM_WORLD);
    void SetGhostArray(const std::vector<int> &recvProc, const std::vector<std::vector<size_t>> &recvPoints);  
    void BringRemoteGhostPoints(const std::vector<BigRangeQueryData> &bigQueries, const std::vector<SmallRangeQueryData> &smallQueries,
                                      BigRangeAgent &bigRangeAgent, SmallRangeAgent &smallRangeAgent,
                                      boost::container::flat_map<size_t, size_t> &numOfResultsForBigPoints,
                                      boost::container::flat_map<size_t, size_t> &numOfResultsForSmallPoints);
  #endif // RICH_MPI

  Delaunay3D del_;
  //vector<vector<std::size_t> > PointTetras_; // The tetras containing each point
  vector<tetra_vec > PointTetras_; // The tetras containing each point
  mutable vector<double> R_; // The radius of the sphere of each tetra
  mutable vector<Vector3D> tetra_centers_;
  // Voronoi Data
  //vector<vector<std::size_t> > FacesInCell_;
  vector<face_vec > FacesInCell_;
  std::vector<point_vec > PointsInFace_; // Right hand with regard to first neighbor
  //vector<vector<std::size_t> > PointsInFace_; // Right hand with regard to first neighbor
  vector<std::pair<std::size_t, std::size_t> > FaceNeighbors_;
  vector<Vector3D> all_CM;
  vector<Vector3D> CM_, Face_CM_; // center of masses
  vector<double> volume_; // volumes of each one of the tetrahedra
  vector<double> area_; // surface area of each one of the tetrahedra
  
  #ifdef RICH_MPI
  vector<int> sentprocs_;
  vector<vector<std::size_t>> sentpoints_; // if rank `i` is inside index `j` in `sentprocs_`, then the points in sentpoints_[j] are the points I sent to rank `i` in the initial points exchange in build
  vector<int> duplicatedprocs_; 
  vector<vector<std::size_t>> duplicated_points_;  // if rank `i` is inside index `j` in `duplicatedprocs_`, then Nghost_[j] includes all the points in `i`'s delaunay, which are actually mine
  vector<int> real_duplicated_proc;
  vector<vector<std::size_t>> real_duplicated_points; // indices of points which are a real ghost points
  vector<vector<std::size_t>> Nghost_; // if rank `i` is inside index `j` in `duplicatedprocs_`, then Nghost_[j] includes all the points in my delaunay, which are belongs, originally, to i
  vector<std::size_t> self_index_; // indexes of the points which are truely mine (inside the points list)
  #endif // RICH_MPI

  Voronoi3D();
  Voronoi3D(Voronoi3D const &other);
  std::array<Vector3D, 4> temp_points_;
  std::array<Vector3D, 5> temp_points2_;
  std::vector<Face> box_faces_;

  std::shared_ptr<OctTree<IndexedVector3D>> myPointsTree;
  std::shared_ptr<OctTree<IndexedVector3D>> allMyPointsTree;
  #ifdef RICH_MPI
    std::shared_ptr<PointsManager<Vector3D, VoronoiPayload>> pointsManager;
    std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>> indexingToSave = std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>>();
  #endif // RICH_MPI

  std::shared_ptr<RangeFinder> rangeFinder;
  std::vector<Vector3D> allMyPoints;
  std::vector<double> allPointsWeights;
  std::vector<double> radiuses;

  Tessellation3D::AllPointsMap indicesInAllMyPoints; // the indices of the points in `del_.points_`, in the list of all points

  size_t buildGeneration_ = 0;

public:

  #ifdef RICH_MPI
    const std::vector<double> &GetPointsBuildWeights() const override;
    
    const std::shared_ptr<EnvironmentAgent<Vector3D>> GetEnvironmentAgent() const override;
    
    void SetKernel(const std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>> &indexing = std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>>());
    

    std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>> GetKernel() const;
    
    void SetBox(Vector3D const &ll, Vector3D const &ur, const std::shared_ptr<const Kernelization3D::IndexingKernel3D<Vector3D>> &newIndexing);
  #endif // RICH_MPI

  #ifdef RICH_MPI
    vector<int> &GetSentProcs(void) override;

    vector<vector<size_t>> &GetSentPoints(void) override;

    vector<size_t> &GetSelfIndex(void) override;
  #endif // RICH_MPI

  vector<Vector3D> &GetAllFaceCM(void) override;

  /*! \brief Calculates the centre of mass of a face
    \param index Face index
    \return Face centre of mass
   */
  const Vector3D &FaceCM(std::size_t index) const override;

  /*! \brief class constructor
    \param ll Lower left
    \param ur Upper right
   */
  Voronoi3D(const Vector3D &ll, const Vector3D &ur);

  Voronoi3D(const std::vector<Face> &box_faces);

  void output(std::string const& filename)const override;

  void BuildInitialize(size_t num_points);

  size_t GetBuildGeneration(void) const override { return buildGeneration_; }

  void BuildPartially(const std::vector<Vector3D> &allPoints, const std::vector<size_t> &indicesToBuild) override;

  bool IsPointInCell(const Vector3D &point, size_t cellIndex, bool verbose = false) const override;

#ifdef RICH_MPI

  /*! \brief Output extra build
    \param filename Output file name
   */
  void output_buildextra(std::string const& filename) const;

  void PreparePoints(const std::vector<Vector3D> &points, const std::vector<size_t> &mask) override;

  inline void SetPointsManager(std::shared_ptr<PointsManager<Vector3D, VoronoiPayload>> pointsManager){this->pointsManager = pointsManager;};

  inline std::shared_ptr<PointsManager<Vector3D, VoronoiPayload>> GetPointsManager() const{return this->pointsManager;};

  void MockMesh(void);
  
  void SetLoadBalancer(std::shared_ptr<LoadBalancer<Vector3D>> loadBalancer) override;

  void PresetLoadBalancer(std::shared_ptr<LoadBalancer<Vector3D>> loadBalancer) override;
  
  void Rebalance(const std::vector<double> &weights) override;

  inline bool ShouldRebalance(const std::vector<double> &weights) const override {return this->pointsManager->shouldRebalance(weights);}

  inline bool ShouldRebalance(void) const override {return this->pointsManager->shouldRebalance();}
    
  inline std::shared_ptr<LoadBalancer<Vector3D>> GetLoadBalancer(void) override {return this->pointsManager->getLoadBalancer();};

  inline const std::shared_ptr<LoadBalancer<Vector3D>> GetLoadBalancer(void) const override {return this->pointsManager->getLoadBalancer();};

  void BuildPartiallyParallel(const std::vector<Vector3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToBuild, bool suppressRebalancing = false, bool suppressExchange = false) override;
  
  bool DidRebalance(void) const override;

  bool CheckContinuityOfZone(void) const;
  
  bool PointInMyDomain(const Vector3D &point) const override;

  int GetOwner(const Vector3D &point) const override;

  void SetImbalanceTolerance(double tolerance) override;

#endif // RICH_MPI

  /*! \brief Dump debug information
    \param rank Rank of parallel process
   */
  void BuildDebug(int rank);

  double GetMaxRadius(const size_t &index) const;

  double GetMinRadius(const size_t &index) const;

  size_t GetContainingCell(const Vector3D &point) const override;

  std::size_t GetPointNo(void) const override;

  /*! \brief Get mehs point position
    \param index Index
    \return Position of point
   */
  const Vector3D &GetMeshPoint(std::size_t index) const override;

  /*! \brief Calculate face area
    \param index Face index
    \return Face area
   */
  double GetArea(std::size_t index) const override;

  /*! \brief Get cell centre of mass
    \param index Point index
    \return Centre of mass
   */
  Vector3D const& GetCellCM(std::size_t index) const override;

  std::size_t GetTotalFacesNumber(void) const override;

  /*! \brief Get cell size
    \param index Point index
    \return Cell width
   */
  double GetWidth(std::size_t index) const override;

  /*! \brief Get cell volume
    \param index Point index
    \return Cell volume
   */
  double GetVolume(std::size_t index) const override;

  /*! \brief Get cell faces
    \param index Point index
    \return List of bounding faces
   */ 
  face_vec const& GetCellFaces(std::size_t index) const override;
  
  vector<Vector3D>& accessMeshPoints(void) override;

  const vector<Vector3D>& getMeshPoints(void) const override;

  const Tessellation3D::AllPointsMap &GetIndicesInAllPoints(void) const override;

  const std::vector<Vector3D> &getAllPoints(void) const override;

  std::vector<Vector3D> &getAllPoints(void) override;

  size_t GetAllPointsNo(void) const override;

  /*! \brief Get neighbours
    \param index Point index
    \return List of indices of neighbouring points
   */
  vector<std::size_t> GetNeighbors(std::size_t index)const override;

  Tessellation3D* clone(void) const override;

  /*! \brief Checs if a point is near a boundary
    \param index Point index
    \return True if point is near a boundary
   */
  bool NearBoundary(std::size_t index) const override;

  /*! \brief Checks if a face is on the boundary
    \param index Face index
    \return True if the face is on the boundary
   */
  bool BoundaryFace(std::size_t index) const override;

  #ifdef RICH_MPI
  vector<vector<std::size_t> >& GetDuplicatedPoints(void) override;

  vector<vector<std::size_t> >const& GetDuplicatedPoints(void)const override;

    /*! \brief Get Duplicated processe
      \return List of duplicated points
    */
    vector<int> GetDuplicatedProcs(void)const override;

  /*! \brief Get a list of parallel processes to which points have been sent
    \return List of process numbers
   */
  vector<int> GetSentProcs(void)const override;

  /*! \brief List of point sent to parallel processes, partitioned by processor number
    \return List of list of indices
   */
  vector<vector<std::size_t> > const& GetSentPoints(void)const override;

  /*! \brief Get indices of all real cells
    \return List of indices of all real cells
   */
  vector<std::size_t> const& GetSelfIndex(void) const override;

  #endif // RICH_MPI
  std::size_t GetTotalPointNumber(void)const override;

  vector<Vector3D> & GetAllCM(void) override;

  vector<Vector3D > GetAllCM(void)const override;

  /*! \brief Get neighbours of neighbours
    \param result Result
    \param point Point index
   */
  void GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point)const override;

  /*! \brief Calculate normal vector to face
    \param faceindex Index of face
    \return Vector normal to face
   */
  Vector3D Normal(std::size_t faceindex)const override;

  /*! \brief Checks if a point is a ghost
    \param index Index
    \return bool if the point is a ghost
   */
  bool IsGhostPoint(std::size_t index)const override;

  /*! \brief Calculate the face velocity
    \param index Face index
    \param v0 Velocity of one neighbour
    \param v1 Velocity of second neighbour
    \return Face velocity
   */
  Vector3D CalcFaceVelocity(std::size_t index, Vector3D const& v0, Vector3D const& v1)const override;

  vector<Vector3D>& GetFacePoints(void) override;

  vector<double>& GetAllArea(void) override;

  vector<Vector3D>const& GetFacePoints(void) const override;

  vector<face_vec >& GetAllCellFaces(void) override;

  vector<face_vec >const& GetAllCellFaces(void) const override;

  /*! \brief Get the points in face
    \param index Face index
    \return Indices of points in face
   */
  point_vec const& GetPointsInFace(std::size_t index) const override;

  /*! \brief Get the neighbours across a face
    \param face_index Index of face
    \return Indices of neighbour across face
   */
  const std::pair<std::size_t, std::size_t> &GetFaceNeighbors(std::size_t face_index) const override;

  #ifdef RICH_MPI
    /*! \brief Get the indices of ghost points
      \return List of list of ghost points
    */
    vector<vector<std::size_t> > const& GetGhostIndeces(void) const override;

    /*! \brief Get the indices of ghost points
      \return List of list of ghost points
    */
    vector<vector<std::size_t> >& GetGhostIndeces(void) override;
  #endif // RICH_MPI

  void GetNeighbors(size_t index, vector<size_t> &res)const override;

  /*! \brief Get the positions of opposite corners on the bounding box
    \return Pair of points on opposite corners
   */
  std::pair<Vector3D, Vector3D> GetBoxCoordinates(void)const override;

  /*! \brief Build tessellation without a box
    \param points Mesh generating points
    \param ghosts Ghost points
    \param toduplicate Indices of points to be duplicated
   */
  void BuildNoBox(vector<Vector3D> const& points, vector<vector<Vector3D> > const& ghosts,vector<size_t> toduplicate) override;

  vector<double>& GetAllVolumes(void) override;

  vector<double> GetAllVolumes(void)const override;

  /*! \brief Get all face neighbours
    \return List of pairs of indices to neighbours
   */
  std::vector<std::pair<size_t, size_t>> &GetAllFaceNeighbors(void) override;

  /*! \brief Get all face neighbours
    \return List of pairs of indices to neighbours
   */
  const std::vector<std::pair<size_t, size_t>> &GetAllFaceNeighbors(void) const override;

  /*! \brief List all points in face
    \return List of all points in face
   */
  vector<point_vec > & GetAllPointsInFace(void) override;

  vector<point_vec > const& GetAllPointsInFace(void) const override;

  /*! \brief Get the number of points
    \return The number of points
   */
  size_t& GetPointNo(void) override;

  /*! \brief Check whether a point is inside the computational domain
    \param index Index of point
    \return bool True if point is inside
   */
  bool IsPointOutsideBox(size_t index)const override;

  bool IsPointOutsideBox(const Vector3D &point) const override;

  /*! \brief Adjust position of the boundary
    \param ll Lower left corner
    \param ur Upper right corner
  */
  void SetBox(Vector3D const& ll, Vector3D const& ur) override;

  const std::vector<Face> &GetBoxFaces(void) const override {return this->box_faces_;}

  std::vector<Face> GetBoxFaces(void) override {return this->box_faces_;}

  std::vector<Face>& ModifyBoxFaces(void) override {return this->box_faces_;}
};

/**
 * TODO: following functions should be here?
*/
bool PointInPoly(Tessellation3D const& tess, Vector3D const& point, std::size_t index);

bool PointInPoly(std::vector<Face> const& faces, Vector3D const &point);

#endif // VORONOI3D_HPP
