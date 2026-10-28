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
#include "RawFaces.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace device_detail
{
    //! \brief The host array itself where the execution space can access it, otherwise a reused device copy.
    template<typename T>
    class DeviceArray
    {
    public:
        static constexpr bool hostAccessible = Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace, Kokkos::HostSpace>::accessible;

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
                {
                    Kokkos::deep_copy(this->DeviceRange(), ConstHostView(host.data(), this->size_));
                }
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
            if constexpr(not hostAccessible)
            {
                if(this->size_ > 0)
                {
                    Kokkos::deep_copy(HostView(host.data(), this->size_), this->DeviceRange());
                }
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
            {
                this->device_ = Kokkos::View<T *>(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string("MadVoro::DeviceArray")), size);
            }
        }

        auto DeviceRange() const
        {
            return Kokkos::subview(this->device_, std::make_pair(std::size_t(0), this->size_));
        }

        Kokkos::View<T *> device_;
        T *data_ = nullptr;
        std::size_t size_ = 0;
    };

    //! \brief Each cell's faces (of the `faces` on the device, between cells [0, cells) and others), in increasing face order, into `hostOut`
    void BuildCellFaceLists(flat_index_t cells, flat_index_t faces, const flat_index_t *neighbor0, const flat_index_t *neighbor1, DeviceArray<flat_index_t> &cellFaceOffsets, DeviceArray<flat_index_t> &cellFaces, FlatVoronoi &hostOut);
} // namespace device_detail

//! \brief Same calls and results as HostGeometryExecutor; work the kernels cannot finish is done on the host.
class KokkosGeometryExecutor
{
public:
    //! \brief `mesh` must outlive the following calls.
    void SetMesh(const FlatMesh &mesh);

    void Circumcenters(const std::vector<std::uint8_t> &compute, FlatVoronoi &out);

    //! \brief `out` (its vertices) must outlive the raw faces.
    void SetVertices(const FlatVoronoi &out);

    //! \brief The raw faces of the mesh's points in `order` (which must outlive them), in the default memory space (on the host if a point's ring does not fit the kernels' scratch)
    FaceBuildReport BuildRawFaces(const std::vector<flat_index_t> &order, RawFacesView &raw);

    void Release(void);

    // Kernel launches (public for CUDA extended lambdas)
    void LaunchFaceCount(std::size_t pointCount);
    void LaunchFaceFill(std::size_t pointCount);

private:
    FlatVoronoiView DeviceView(const FlatVoronoi &out) const;

    // Runs the calls made in a SerialGeometryScope
    HostGeometryExecutor host_;

    const FlatMesh *hostMesh_ = nullptr;
    FlatMeshView mesh_;
    device_detail::DeviceArray<FlatPoint> points_;
    device_detail::DeviceArray<FlatTetra> tetras_;
    device_detail::DeviceArray<flat_index_t> pointTetraOffsets_, pointTetras_;

    device_detail::DeviceArray<std::uint8_t> compute_;
    device_detail::DeviceArray<flat_index_t> order_;
    device_detail::DeviceArray<FlatPoint> vertices_;
    device_detail::DeviceArray<double> tetraRadius_;
    device_detail::DeviceArray<GeometryStatus> tetraStatus_;

    const std::vector<FlatPoint> *hostVertices_ = nullptr;

    // Per point counts and offsets of the two-pass raw face build, and the raw faces
    Kokkos::View<FaceBuildResult *> results_;
    Kokkos::View<flat_index_t *> faceOffsets_, vertexOffsets_;
    Kokkos::View<flat_index_t *> rawCell_, rawOther_, rawVertexOffsets_, rawVertices_;
    HostRawFaces hostRaw_;
};

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KOKKOS_GEOMETRY_EXECUTOR_HPP
