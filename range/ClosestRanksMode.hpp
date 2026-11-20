#ifndef MADVORO_CLOSEST_RANKS_MODE_HPP
#define MADVORO_CLOSEST_RANKS_MODE_HPP

#include <stdexcept>
#include <string>

namespace MadVoro
{
    // How the big-query talk agent evaluates closest/furthest distances to the ranks a query sphere intersects, on a
    // point's first large ghost round (MADVORO_CLOSEST_RANKS_RESTRICTED in RICH):
    //   0  dense: one scan of every routing leaf for all ranks per point (default);
    //   1  restricted: only the leaves of the ranks the query intersects, through a per-owner leaf index;
    //   2  verify: both, selecting with the dense values and counting bit differences; a difference on any rank
    //      makes every rank of the build throw ClosestRanksVerifyError at the end of the ghost loop.
    // Both paths evaluate each leaf with the same non-inlined routine, so modes 0 and 1 give identical talk lists.
    // Process-wide: the caller must set the same value on every rank before the first build that uses it.
    inline int &ClosestRanksRestrictedModeStorage()
    {
        static int Mode = 0;
        return Mode;
    }

    inline void SetClosestRanksRestrictedMode(int Mode)
    {
        if(Mode < 0 || Mode > 2)
        {
            throw std::invalid_argument("MadVoro::SetClosestRanksRestrictedMode: mode must be 0, 1 or 2");
        }
        ClosestRanksRestrictedModeStorage() = Mode;
    }

    inline int GetClosestRanksRestrictedMode()
    {
        return ClosestRanksRestrictedModeStorage();
    }

    namespace Exception
    {
        // Thrown together by every rank of a build when mode 2 found restricted and dense distances that differ. It
        // is not a MadVoroException, so callers that turn build failures into a full-build fallback can let it pass.
        class ClosestRanksVerifyError : public std::runtime_error
        {
        public:
            explicit ClosestRanksVerifyError(const std::string &Message) : std::runtime_error(Message) {}
        };
    }
}

#endif // MADVORO_CLOSEST_RANKS_MODE_HPP
