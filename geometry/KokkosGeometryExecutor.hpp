#ifndef MADVORO_GEOMETRY_KOKKOS_GEOMETRY_EXECUTOR_HPP
#define MADVORO_GEOMETRY_KOKKOS_GEOMETRY_EXECUTOR_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <Kokkos_Core.hpp>

#include "FlatMesh.hpp"
#include "FlatVoronoi.hpp"
#include "HostGeometryExecutor.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace device_detail
{
    /*! \brief The memory a kernel uses for one host array. Where the execution space
        can access host memory (e.g. the OpenMP backend) that is the host array itself;
        otherwise it is a device copy, kept between builds and grown when needed. */
    template<typename T>
    class DeviceArray
    {
    public:
        static constexpr bool hostAccessible =
            Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace, Kokkos::HostSpace>::accessible;

        //! \brief Kernels read `host`'s current contents.
        void Upload(const std::vector<T> &host)
        {
            this->size_ = host.size();
            if constexpr(hostAccessible)
            {
                this->data_ = const_cast<T *>(host.data());
            }
            else
            {
                this->Reserve(this->size_);
                if(this->size_ > 0)
                    Kokkos::deep_copy(this->DeviceRange(), ConstHostView(host.data(), this->size_));
                this->data_ = this->device_.data();
            }
        }

        //! \brief Kernels write `host`'s new contents (bring them back with Download).
        void Output(std::vector<T> &host)
        {
            this->size_ = host.size();
            if constexpr(hostAccessible)
            {
                this->data_ = host.data();
            }
            else
            {
                this->Reserve(this->size_);
                this->data_ = this->device_.data();
            }
        }

        //! \brief Copies what kernels wrote into `host` (nothing to do when they wrote `host` itself).
        void Download(std::vector<T> &host) const
        {
            if constexpr(!hostAccessible)
            {
                if(this->size_ > 0)
                    Kokkos::deep_copy(HostView(host.data(), this->size_), this->DeviceRange());
            }
        }

        T *Data() const
        {
            return this->data_;
        }

    private:
        using HostView = Kokkos::View<T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        using ConstHostView = Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

        void Reserve(std::size_t size)
        {
            if(this->device_.extent(0) < size)
                this->device_ = Kokkos::View<T *>(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string("MadVoro::DeviceArray")), size);
        }

        auto DeviceRange() const
        {
            return Kokkos::subview(this->device_, std::make_pair(std::size_t(0), this->size_));
        }

        Kokkos::View<T *> device_;
        T *data_ = nullptr;
        std::size_t size_ = 0;
    };
} // namespace device_detail

/*! \brief Runs the geometry kernels with Kokkos (on a device, or on host threads).
    Same calls and results as HostGeometryExecutor: SetMesh, Circumcenters,
    SetVertices, BuildFaces, CellMoments, with inputs and results in the host
    FlatMesh and FlatVoronoi. Work the kernels cannot finish (extended-precision
    circumcenters, exact volumes, faces that do not fit the fixed scratch) is done
    on the host.
    The kernel launches are public member functions because CUDA's extended
    lambdas cannot be defined in private ones. */
class KokkosGeometryExecutor
{
public:
    //! \brief Makes `mesh` the input of the following calls; it must outlive them.
    void SetMesh(const FlatMesh &mesh);

    //! \brief See HostGeometryExecutor::Circumcenters.
    void Circumcenters(const std::vector<std::uint8_t> &compute, FlatVoronoi &out);

    //! \brief Makes out.vertices the face vertices for BuildFaces.
    void SetVertices(const FlatVoronoi &out);

    //! \brief See BuildFacesHost.
    FaceBuildReport BuildFaces(const std::vector<flat_index_t> &order, const std::vector<std::uint8_t> &looseAreaCut,
                               FlatVoronoi &out);

    //! \brief Volumes and centroids of the cells of the last BuildFaces, into `out`.
    void CellMoments(FlatVoronoi &out);

    // Kernel launches (public for CUDA extended lambdas)
    void LaunchFaceCount(std::size_t pointCount);
    void LaunchFaceFill(std::size_t pointCount, const FlatVoronoiView &view);
    void LaunchCellFaceLists(std::size_t pointCount, std::size_t faceCount, FlatVoronoi &hostOut);

private:
    FlatVoronoiView DeviceView(const FlatVoronoi &out) const;

    const FlatMesh *hostMesh_ = nullptr;
    FlatMeshView mesh_;
    device_detail::DeviceArray<FlatPoint> points_;
    device_detail::DeviceArray<FlatTetra> tetras_;
    device_detail::DeviceArray<flat_index_t> pointTetraOffsets_, pointTetras_;

    device_detail::DeviceArray<std::uint8_t> compute_, looseAreaCut_;
    device_detail::DeviceArray<flat_index_t> order_;
    device_detail::DeviceArray<FlatPoint> vertices_;
    device_detail::DeviceArray<double> tetraRadius_;
    device_detail::DeviceArray<GeometryStatus> tetraStatus_;

    device_detail::DeviceArray<flat_index_t> cellFaceOffsets_, cellFaces_, faceVertexOffsets_, faceVertices_;
    device_detail::DeviceArray<flat_index_t> faceNeighbor0_, faceNeighbor1_;
    device_detail::DeviceArray<double> faceArea_;
    device_detail::DeviceArray<FlatPoint> faceCM_;
    device_detail::DeviceArray<double> cellVolume_;
    device_detail::DeviceArray<FlatPoint> cellCM_;
    device_detail::DeviceArray<GeometryStatus> cellStatus_;
    device_detail::DeviceArray<flat_index_t> cellOrder_;

    // Per point counts and offsets of the two-pass face build
    Kokkos::View<FaceBuildResult *> results_;
    Kokkos::View<flat_index_t *> faceOffsets_, vertexOffsets_;
};

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KOKKOS_GEOMETRY_EXECUTOR_HPP
