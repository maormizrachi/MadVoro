/** \file Voronoi3D.hpp
   \brief A 3D Voronoi tessellation, templated on a user-provided vector type.
   \author Elad Steinberg, Maor Mizrachi
*/
#ifndef MADVORO_VORONOI3D_HPP
#define MADVORO_VORONOI3D_HPP

#include <vector>
#include <memory>
#include <utility>
#include "madvoro/types.hpp"
#include "madvoro/Face.hpp"
#include "voronoi/Voronoi3DFull.hpp"

namespace MadVoro
{
  class LoadBalancer;

  /**
   * @brief 3D Voronoi tessellation, templated on a user-provided 3D vector type.
   *
   * @tparam Vec3 A type with public `double x, y, z` members and that is
   *              constructible from three doubles via `Vec3{x, y, z}`.
   *
   * Internally, the tessellation uses its own point type (Point3D).
   * Conversions between Vec3 and Point3D happen automatically at the API boundary.
   */
  template<typename Vec3>
  class Voronoi3D
  {
  public:
    Voronoi3D(const Vec3 &ll, const Vec3 &ur)
      : pImpl(new Voronoi3DFull(Point3D(ll.x, ll.y, ll.z), Point3D(ur.x, ur.y, ur.z))) {}

    Voronoi3D(const std::vector<Face<Vec3>> &box_faces)
      : pImpl(new Voronoi3DFull(toFaces3D(box_faces))) {}

    ~Voronoi3D() { delete pImpl; }

    Voronoi3D(const Voronoi3D&) = delete;
    Voronoi3D& operator=(const Voronoi3D&) = delete;

    // ---------------------- Build ----------------------

    void Build(const std::vector<Vec3> &points)
    {
      pImpl->Build(toPoints(points));
    }

    void BuildPartially(const std::vector<Vec3> &allPoints, const std::vector<std::size_t> &indicesToBuild)
    {
      pImpl->BuildPartially(toPoints(allPoints), indicesToBuild);
    }

    #ifdef MADVORO_WITH_MPI
    const std::vector<double> &GetPointsBuildWeights() const { return pImpl->GetPointsBuildWeights(); }

    std::vector<Vec3> BuildParallel(const std::vector<Vec3> &points, const std::vector<double> &weights, bool suppressRebalancing = false)
    {
      return fromPoints(pImpl->BuildParallel(toPoints(points), weights, suppressRebalancing));
    }

    std::vector<Vec3> BuildParallel(const std::vector<Vec3> &points, bool suppressRebalancing = false)
    {
      return fromPoints(pImpl->BuildParallel(toPoints(points), suppressRebalancing));
    }

    std::vector<Vec3> BuildPartiallyParallel(const std::vector<Vec3> &allPoints, const std::vector<double> &allWeights, const std::vector<std::size_t> &indicesToBuild, bool suppressRebalancing = false)
    {
      return fromPoints(pImpl->BuildPartiallyParallel(toPoints(allPoints), allWeights, indicesToBuild, suppressRebalancing));
    }

    bool PointInMyDomain(const Vec3 &point) const { return pImpl->PointInMyDomain(toPoint(point)); }
    int GetOwner(const Vec3 &point) const { return pImpl->GetOwner(toPoint(point)); }
    void MockMesh(void) { pImpl->MockMesh(); }
    void SetLoadBalancer(std::shared_ptr<LoadBalancer> loadBalancer) { pImpl->SetLoadBalancer(loadBalancer); }
    void Rebalance(const std::vector<double> &weights) { pImpl->Rebalance(weights); }
    void SetImbalanceTolerance(double tolerance) { pImpl->SetImbalanceTolerance(tolerance); }
    #endif // MADVORO_WITH_MPI

    // ---------------------- Query ----------------------

    std::size_t GetContainingCell(const Vec3 &point) const { return pImpl->GetContainingCell(toPoint(point)); }
    Vec3 FaceCM(std::size_t index) const { return fromPoint(pImpl->FaceCM(index)); }
    std::size_t GetPointNo(void) const { return pImpl->GetPointNo(); }
    Vec3 GetMeshPoint(std::size_t index) const { return fromPoint(pImpl->GetMeshPoint(index)); }
    double GetArea(std::size_t faceIndex) const { return pImpl->GetArea(faceIndex); }
    Vec3 GetCellCM(std::size_t index) const { return fromPoint(pImpl->GetCellCM(index)); }
    std::size_t GetTotalFacesNumber(void) const { return pImpl->GetTotalFacesNumber(); }
    double GetWidth(std::size_t index) const { return pImpl->GetWidth(index); }
    double GetVolume(std::size_t index) const { return pImpl->GetVolume(index); }
    const face_vec &GetCellFaces(std::size_t index) const { return pImpl->GetCellFaces(index); }

    std::vector<Vec3> getMeshPoints(void) const { return fromPoints(pImpl->getMeshPoints()); }
    const AllPointsMap &GetIndicesInAllPoints(void) const { return pImpl->GetIndicesInAllPoints(); }
    std::vector<Vec3> getAllPoints(void) const { return fromPoints(pImpl->getAllPoints()); }
    std::size_t GetAllPointsNo(void) const { return pImpl->GetAllPointsNo(); }
    std::vector<std::size_t> GetNeighbors(std::size_t index) const { return pImpl->GetNeighbors(index); }
    bool NearBoundary(std::size_t index) const { return pImpl->NearBoundary(index); }
    bool BoundaryFace(std::size_t index) const { return pImpl->BoundaryFace3D(index); }

    #ifdef MADVORO_WITH_MPI
    const std::vector<std::vector<std::size_t>> &GetDuplicatedPoints(void) const { return pImpl->GetDuplicatedPoints(); }
    std::vector<int> GetDuplicatedProcs(void) const { return pImpl->GetDuplicatedProcs(); }
    std::vector<int> GetSentProcs(void) const { return pImpl->GetSentProcs(); }
    const std::vector<std::vector<std::size_t>> &GetSentPoints(void) const { return pImpl->GetSentPoints(); }
    const std::vector<std::size_t> &GetSelfIndex(void) const { return pImpl->GetSelfIndex(); }
    const std::vector<std::vector<std::size_t>> &GetGhostIndeces(void) const { return pImpl->GetGhostIndeces(); }
    std::vector<std::vector<std::size_t>> &GetGhostIndeces(void) { return pImpl->GetGhostIndeces(); }
    #endif // MADVORO_WITH_MPI

    std::size_t GetTotalPointNumber(void) const { return pImpl->GetTotalPointNumber(); }
    std::vector<Vec3> GetAllCM(void) const { return fromPoints(pImpl->GetAllCM()); }
    std::vector<Vec3> GetAllFaceCM(void) const { return fromPoints(pImpl->GetAllFaceCM()); }
    Vec3 Normal(std::size_t faceindex) const { return fromPoint(pImpl->Normal(faceindex)); }
    bool IsGhostPoint(std::size_t index) const { return pImpl->IsGhostPoint(index); }
    std::vector<Vec3> GetFacePoints(void) const { return fromPoints(pImpl->GetFacePoints()); }
    const std::vector<face_vec> &GetAllCellFaces(void) const { return pImpl->GetAllCellFaces(); }
    const point_vec &GetPointsInFace(std::size_t index) const { return pImpl->GetPointsInFace(index); }
    const std::pair<std::size_t, std::size_t> &GetFaceNeighbors(std::size_t face_index) const { return pImpl->GetFaceNeighbors(face_index); }
    void GetNeighbors(std::size_t index, std::vector<std::size_t> &res) const { pImpl->GetNeighbors(index, res); }

    std::pair<Vec3, Vec3> GetBoxCoordinates(void) const
    {
      auto res = pImpl->GetBoxCoordinates();
      return std::make_pair(fromPoint(res.first), fromPoint(res.second));
    }

    std::vector<double> GetAllVolumes(void) const { return pImpl->GetAllVolumes(); }
    const std::vector<std::pair<std::size_t, std::size_t>> &GetAllFaceNeighbors(void) const { return pImpl->GetAllFaceNeighbors(); }
    const std::vector<point_vec> &GetAllPointsInFace(void) const { return pImpl->GetAllPointsInFace(); }
    bool IsPointOutsideBox(std::size_t index) const { return pImpl->IsPointOutsideBox(index); }

    void SetBox(const Vec3 &ll, const Vec3 &ur)
    {
      pImpl->SetBox(toPoint(ll), toPoint(ur));
    }

    std::vector<Face<Vec3>> GetBoxFaces(void) const
    {
      return fromFaces3D(pImpl->GetBoxFaces());
    }

    void SetVerbosity(bool value) { pImpl->SetVerbosity(value); }

    #ifdef MADVORO_WITH_HDF5
    void ToHDF5(const std::string &fileName, const std::vector<std::string> &fieldNames = {}, const std::vector<std::vector<double>> &fieldValues = {})
    {
      pImpl->ToHDF5(fileName, fieldNames, fieldValues);
    }
    #endif

    #ifdef MADVORO_WITH_VTK
    void ToVTK(const std::string &fileName, const std::vector<std::string> &fieldNames = {}, const std::vector<std::vector<double>> &fieldValues = {})
    {
      pImpl->ToVTK(fileName, fieldNames, fieldValues);
    }
    #endif

  private:
    Voronoi3DFull *pImpl = nullptr;

    // ---------------------- Conversion helpers ----------------------

    static Point3D toPoint(const Vec3 &v) { return Point3D(v.x, v.y, v.z); }

    static Vec3 fromPoint(const Point3D &p) { return Vec3{p.x, p.y, p.z}; }

    static std::vector<Point3D> toPoints(const std::vector<Vec3> &vecs)
    {
      std::vector<Point3D> result;
      result.reserve(vecs.size());
      for (const auto &v : vecs)
        result.emplace_back(v.x, v.y, v.z);
      return result;
    }

    static std::vector<Vec3> fromPoints(const std::vector<Point3D> &pts)
    {
      std::vector<Vec3> result;
      result.reserve(pts.size());
      for (const auto &p : pts)
        result.push_back(Vec3{p.x, p.y, p.z});
      return result;
    }

    static std::vector<Face3D> toFaces3D(const std::vector<Face<Vec3>> &faces)
    {
      std::vector<Face3D> result;
      result.reserve(faces.size());
      for (const auto &f : faces)
      {
        Face3D f3;
        f3.neighbors = f.neighbors;
        for (const auto &v : f.vertices)
          f3.vertices.emplace_back(v.x, v.y, v.z);
        result.push_back(std::move(f3));
      }
      return result;
    }

    static std::vector<Face<Vec3>> fromFaces3D(const std::vector<Face3D> &faces)
    {
      std::vector<Face<Vec3>> result;
      result.reserve(faces.size());
      for (const auto &f : faces)
      {
        Face<Vec3> face;
        face.neighbors = f.neighbors;
        for (const auto &p : f.vertices)
          face.vertices.push_back(Vec3{p.x, p.y, p.z});
        result.push_back(std::move(face));
      }
      return result;
    }
  };
}

#endif // MADVORO_VORONOI3D_HPP
