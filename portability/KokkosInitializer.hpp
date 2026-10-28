#ifndef MADVORO_KOKKOS_INITIALIZER_HPP
#define MADVORO_KOKKOS_INITIALIZER_HPP

#include <Kokkos_Core.hpp>

namespace MadVoro
{

// Initializes Kokkos on first use and finalizes it at exit, unless another library already did.
class KokkosInitializer
{
public:
    static KokkosInitializer &Instance()
    {
        // Function-local: constructed after Kokkos' static initializers, destroyed before the device runtime's exit handlers
        static KokkosInitializer instance;
        return instance;
    }

    // E.g. before MPI_Finalize when HIP shares resources with an OFI context. Kokkos cannot be initialized again afterwards.
    static void FinalizeIfInitialized()
    {
        if(Created())
        {
            Instance().Finalize();
        }
    }

    KokkosInitializer(const KokkosInitializer &) = delete;
    KokkosInitializer &operator=(const KokkosInitializer &) = delete;

    ~KokkosInitializer()
    {
        this->Finalize();
    }

private:
    KokkosInitializer()
    {
        Created() = true;
        if(not Kokkos::is_initialized() and not Kokkos::is_finalized())
        {
            Kokkos::initialize();
            this->initializedHere_ = true;
        }
    }

    void Finalize()
    {
        if(this->initializedHere_ and Kokkos::is_initialized() and not Kokkos::is_finalized())
        {
            Kokkos::finalize();
        }
        this->initializedHere_ = false;
    }

    static bool &Created()
    {
        static bool created = false;
        return created;
    }

    bool initializedHere_ = false;
};

} // namespace MadVoro

#endif // MADVORO_KOKKOS_INITIALIZER_HPP
