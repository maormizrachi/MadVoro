/** \file Voronoi3D.hpp
   \brief A 3D Voronoi
   \author Elad Steinberg, Maor Mizrachi
*/
#ifndef VORONOI3D_HPP
#define VORONOI3D_HPP 1

#include <vector>
#include <boost/container/small_vector.hpp>
#include "elementary/Face.hpp"
#include "elementary/Vector3D.hpp"
// #include "exception/InvalidArgumentException.hpp"

namespace MadVoro
{
  //! \brief Container for points defining a face
  typedef boost::container::small_vector<size_t, 24> face_vec;
  //! \brief Container for neighbouring points
  typedef boost::container::small_vector<size_t, 8> point_vec;

  class Voronoi3D
  {
    using AllPointsMap = boost::container::flat_map<std::size_t, std::size_t>;

  public:    
    /**
     * @brief Initializes a 3D Voronoi tessellation with the given computational domain.
     * @param ll The lower left corner of the computational domain.
     * @param ur The upper right corner of the computational domain.
     * @note This constructor initializes the Voronoi tessellation with the given computational domain.
     *       It does not perform any tessellation calculations. Use the `Build` methods for that.
     */
    Voronoi3D(const Vector3D &ll, const Vector3D &ur);

    /**
     * @brief Initializes a 3D Voronoi tessellation with the given computational domain defined by box faces.
     * @param box_faces A vector of Face objects representing the faces of the computational domain.
     * @note This constructor initializes the Voronoi tessellation with the given computational domain.
     *       It does not perform any tessellation calculations. Use the `Build` methods for that.
     * @note The computational domain is defined by the provided box_faces. The faces should be ordered
     *       such that // TODO how
     */
    Voronoi3D(const std::vector<Face> &box_faces);

    ~Voronoi3D();
    
    std::vector<Vector3D> GetAllFaceCM(void) const;

    void BuildPartially(const std::vector<Vector3D> &allPoints, const std::vector<std::size_t> &indicesToBuild);

    void Build(const std::vector<Vector3D> &points);
    #ifdef MADVORO_WITH_MPI
    const std::vector<double> &GetPointsBuildWeights() const;
        
    std::vector<Vector3D> BuildParallel(const std::vector<Vector3D> &points, const std::vector<double> &weights, bool suppressRebalancing = false);

    std::vector<Vector3D> BuildParallel(const std::vector<Vector3D> &points, bool suppressRebalancing = false);

    std::vector<Vector3D> BuildPartiallyParallel(const std::vector<Vector3D> &allPoints, const std::vector<double> &allWeights, const std::vector<std::size_t> &indicesToBuild, bool suppressRebalancing = false);

    /**
     * @brief Gets a point, and return whether this point in my domain or not.
     * @param point a point (three dimensional)
     * @return true if the point is in my domain
     * @return false if the point is not in my domain
     */
    bool PointInMyDomain(const Vector3D &point) const;

    /**
     * @brief Returns the owner rank (MPI rank) of the given point.
     * @param point a point to find its affiliation.
     * @return int the rank number.
     */
    int GetOwner(const Vector3D &point) const;
    #endif // MADVORO_WITH_MPI
    /**
     * @brief Finds the local containing cell for a given point in the 3D Voronoi tessellation.
     * @param point The 3D point for which the containing cell needs to be found.
     * @return The index of the LOCAL cell that contains the given point.
     * @note This function assumes that the Voronoi tessellation has already been 
     * constructed and the mesh points have been set.
     */
    std::size_t GetContainingCell(const Vector3D &point) const;


    /*! \brief Calculates the centre of mass of a face
        \param index Face index
        \return Face centre of mass
    */
    Vector3D FaceCM(std::size_t index) const;

    /**
     * @brief Retrieves the number of points assigned to the current processor in the Voronoi tessellation.
     * @return The total number of points in the Voronoi tessellation.
     */
    std::size_t GetPointNo(void) const;

    /*! \brief Get mehs point position
        \param index Index
        \return Position of point
    */
    Vector3D GetMeshPoint(std::size_t index) const;

    /**
     * @brief Retrieves the area of a specific face in the 3D Voronoi tessellation.
     * @param index The index of the face for which the area needs to be calculated.
     * @return The area of the face with the given index.
     * @note This function assumes that the Voronoi tessellation has already been constructed.
     */
    double GetArea(std::size_t faceIndex) const;

    /*! \brief Get cell centre of mass
        \param index Point index
        \return Centre of mass
    */
    Vector3D GetCellCM(std::size_t index) const;

    /**
     * @brief Retrieves the total number of faces in the 3D Voronoi tessellation.
     * @return The total number of faces in the Voronoi tessellation.
     * @note This function assumes that the Voronoi tessellation has already been constructed.
     */
    std::size_t GetTotalFacesNumber(void) const;

    /**
     * @brief Retrieves the width of a specific cell in the 3D Voronoi tessellation.
     * @param index The index of the cell for which the width needs to be calculated.
     * @return The width of the cell with the given index.
     * @note This function assumes that the Voronoi tessellation has already been constructed.
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
    const face_vec &GetCellFaces(std::size_t index) const;
    
    std::vector<Vector3D> getMeshPoints(void) const;

    const AllPointsMap &GetIndicesInAllPoints(void) const;

    std::vector<Vector3D> getAllPoints(void) const;

    std::size_t GetAllPointsNo(void) const;

    /*! \brief Get neighbours
        \param index Point index
        \return List of indices of neighbouring points
    */
    std::vector<std::size_t> GetNeighbors(std::size_t index) const;

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
    // Communication methods
    const std::vector<std::vector<std::size_t>> &GetDuplicatedPoints(void) const;

    /*! \brief Get Duplicated processe
        \return List of duplicated points
    */
    std::vector<int> GetDuplicatedProcs(void) const;

    /*! \brief Get a list of parallel processes to which points have been sent
        \return List of process numbers
    */
    std::vector<int> GetSentProcs(void) const;

    /*! \brief List of point sent to parallel processes, partitioned by processor number
        \return List of list of indices
    */
    const std::vector<std::vector<std::size_t>> &GetSentPoints(void) const;

    /*! \brief Get indices of all real cells
        \return List of indices of all real cells
    */
    const std::vector<std::size_t> &GetSelfIndex(void) const;

    /*! \brief Get the indices of ghost points
        \return List of list of ghost points
    */
    const std::vector<std::vector<std::size_t>> &GetGhostIndeces(void) const;

    /*! \brief Get the indices of ghost points
        \return List of list of ghost points
    */
    std::vector<std::vector<std::size_t>> &GetGhostIndeces(void);
    #endif // MADVORO_WITH_MPI

    std::size_t GetTotalPointNumber(void) const;

    std::vector<Vector3D> GetAllCM(void) const;

    /*! \brief Calculate normal vector to face
        \param faceindex Index of face
        \return Vector normal to face
    */
    Vector3D Normal(std::size_t faceindex) const;

    /*! \brief Checks if a point is a ghost
        \param index Index
        \return bool if the point is a ghost
    */
    bool IsGhostPoint(std::size_t index) const;

    std::vector<Vector3D> GetFacePoints(void) const;

    const std::vector<face_vec> &GetAllCellFaces(void) const;

    /*! \brief Get the points in face
        \param index Face index
        \return Indices of points in face
    */
    const point_vec &GetPointsInFace(std::size_t index) const;

    /*! \brief Get the neighbours across a face
        \param face_index Index of face
        \return Indices of neighbour across face
    */
    std::pair<std::size_t, std::size_t> GetFaceNeighbors(std::size_t face_index) const;

    /**
     * @brief Retrieves the indices of neighboring cells for a given cell in the 3D Voronoi tessellation.
     * @param index The index of the cell for which the neighbors need to be retrieved.
     * @param res A reference to a vector that will be populated with the indices of the neighboring cells.
     * @return void
     * @note This function assumes that the Voronoi tessellation has already been constructed and the mesh points have been set.
     *       The indices in the 'res' vector correspond to the indices of the cells in the Voronoi tessellation.
     */
    void GetNeighbors(size_t index, std::vector<std::size_t> &res) const;

    /*! \brief Get the positions of opposite corners on the bounding box
        \return Pair of points on opposite corners
    */
    std::pair<Vector3D, Vector3D> GetBoxCoordinates(void) const;

    /**
     * @brief Retrieves the volumes of all cells in the 3D Voronoi tessellation.
     * @return A vector of double values representing the volumes of all cells.
     *         The index of each volume corresponds to the index of the cell in the Voronoi tessellation.
     * @note This function assumes that the Voronoi tessellation has already been constructed.
     *       The returned vector will have the same size as the total number of cells in the tessellation.
     */
    std::vector<double> GetAllVolumes(void) const;

    /*! \brief Get all face neighbours
        \return List of pairs of indices to neighbours
    */
    const std::vector<std::pair<std::size_t, std::size_t>> &GetAllFaceNeighbors(void) const;

    /**
     * @brief Retrieves the indices of all points that form each face in the 3D Voronoi tessellation.
     * @return A const reference to a vector of point_vec objects.
     *         Each point_vec object contains the indices of points that form a specific face.
     *         The index of each point_vec corresponds to the index of the face in the Voronoi tessellation.
     * @note This function assumes that the Voronoi tessellation has already been constructed.
     *       The returned vector will have the same size as the total number of faces in the tessellation.
     */
    const std::vector<point_vec> &GetAllPointsInFace(void) const;

    /**
     * @brief Checks whether a point is inside the computational domain.
     * @param index The index of the point to be checked.
     * @return bool 
     *  - True if the point is outside the computational domain.
     *  - False if the point is inside the computational domain.
     */
    bool IsPointOutsideBox(size_t index) const;

    /**
     * @brief Sets the computational domain box for the Voronoi tessellation.
     * @param ll The lower left corner of the computational domain.
     * @param ur The upper right corner of the computational domain.
     * @return void
     */
    void SetBox(const Vector3D &ll, const Vector3D &ur);

    /**
     * @brief Retrieves the faces of the computational domain.
     * @return A vector of Face objects representing the faces of the computational domain.
     */
    std::vector<Face> GetBoxFaces(void) const;

    void SetVerbosity(bool value);
    
    #ifdef MADVORO_WITH_HDF5
    void ToHDF5(const std::string &fileName, const std::vector<std::string> &fieldNames = std::vector<std::string>(), const std::vector<std::vector<double>> &fieldValues = std::vector<std::vector<double>>());
    #endif // MADVORO_WITH_HDF5
    #ifdef MADVORO_WITH_VTK
    void ToVTK(const std::string &fileName, const std::vector<std::string> &fieldNames = std::vector<std::string>(), const std::vector<std::vector<double>> &fieldValues = std::vector<std::vector<double>>());
    #endif // MADVORO_WITH_VTK

  private:
    class Voronoi3DImpl;
    using impl = Voronoi3DImpl;

    Voronoi3DImpl *pImpl = nullptr;
  };
}

#endif // VORONOI3D_HPP
