// Answers a batch of tree queries on host threads, with the answers that one query after the other would give
#ifndef MADVORO_RANGE_PARALLEL_ANSWERS_HPP
#define MADVORO_RANGE_PARALLEL_ANSWERS_HPP

#include <algorithm>
#include <cstddef>
#include <vector>

#include "geometry/HostGeometryExecutor.hpp"

namespace range_detail
{
    // batches smaller than this are answered on one thread
    constexpr std::size_t PARALLEL_BATCH = 64;

    /** \brief The answers to n queries, each excluding the points the queries before it took. find(i, ignore): query i's points outside `ignore` (at most limit(i)); ignoreOf(i): what query i must exclude now; accept(i, points, ignore): records query i's answer (which the queries after it exclude). Every query is first found against what was excluded before the batch, in parallel (the finder only reads); then, in order, the points an earlier query took are dropped, and a query its limit may have cut is found again, so the answers are the serial ones */
    template<typename Find, typename IgnoreOf, typename Limit, typename Accept>
    std::vector<std::vector<std::size_t>> ParallelAnswers(std::size_t n, const Find &find, const IgnoreOf &ignoreOf, const Limit &limit, const Accept &accept)
    {
        std::vector<std::vector<std::size_t>> found(n);
        const auto findOne = [&](std::size_t i)
        {
            found[i] = find(i, ignoreOf(i));
        };
        if(n < PARALLEL_BATCH)
        {
            for(std::size_t i = 0; i < n; ++i)
            {
                findOne(i);
            }
        }
        else
        {
            MadVoro::Geometry::HostParallelFor(n, findOne);
        }
        for(std::size_t i = 0; i < n; ++i)
        {
            auto &ignore = ignoreOf(i);
            const auto taken = [&ignore](std::size_t point)
            {
                return ignore.find(point) != ignore.end();
            };
            if(std::any_of(found[i].begin(), found[i].end(), taken))
            {
                if(found[i].size() >= limit(i))
                {
                    found[i] = find(i, ignore);
                }
                else
                {
                    found[i].erase(std::remove_if(found[i].begin(), found[i].end(), taken), found[i].end());
                }
            }
            accept(i, found[i], ignore);
        }
        return found;
    }
} // namespace range_detail

#endif // MADVORO_RANGE_PARALLEL_ANSWERS_HPP
