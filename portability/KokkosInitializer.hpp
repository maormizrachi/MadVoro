#ifndef MADVORO_KOKKOS_INITIALIZER_HPP
#define MADVORO_KOKKOS_INITIALIZER_HPP

#include <Kokkos_Core.hpp>

namespace MadVoro
{

// Makes sure Kokkos is initialized when MadVoro first needs it, and finalizes it
// automatically at program exit, so drivers never call Kokkos::initialize/finalize.
//
// Every library using Kokkos (MadVoro, STORM, the host application) may have an
// initializer like this one. Whichever runs first starts Kokkos and owns it;
// the others find it already initialized and never finalize it. Static objects are
// destroyed in reverse order of construction, so the owner, created first, is
// destroyed after every other initializer.
//
// The single instance is a function-local static created on first use:
//  - Kokkos registers its execution spaces from its own static initializers. A
//    namespace-scope global could run Kokkos::initialize before those, since the
//    order of static initialization across translation units is unspecified.
//  - Kokkos::initialize starts the device runtime (CUDA/HIP), which registers its
//    exit handlers while this constructor runs, so this destructor runs before
//    them, while the device runtime is still usable. Objects created in main are
//    destroyed before this, when main returns.
//  - A local static of an inline function is one object for the whole program; a
//    variable in an anonymous namespace would be duplicated in every file that
//    includes this header.
class KokkosInitializer
{
public:
    static KokkosInitializer &Instance()
    {
        static KokkosInitializer instance;
        return instance;
    }

    // Finalize Kokkos before program exit if MadVoro initialized it, e.g. before
    // MPI_Finalize when HIP shares resources with an OFI context. Does nothing when
    // MadVoro never started Kokkos (or another library owns it). Kokkos cannot be
    // initialized again afterwards.
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
