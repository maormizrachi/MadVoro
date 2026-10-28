// Brings the ghost points of a rank's meshes (one, or its thread subdomains'), round by round
#ifndef MADVORO_GHOSTS_GHOST_BRINGER_HPP
#define MADVORO_GHOSTS_GHOST_BRINGER_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#include "ghosts/GhostSearch.hpp"
#include "ghosts/GhostExchange.hpp"

namespace MadVoro
{

//! \brief A mesh and its ghost search
template <typename PointT>
struct GhostDomain
{
    GhostSearch<PointT> *search;
    GhostSearchMesh<PointT> *mesh;
};

//! \brief What becomes of the points other ranks send
template <typename PointT>
class GhostReceiver
{
public:
    virtual ~GhostReceiver() = default;

    //! \brief Whether the points are kept for the meshes to take (by index among all the rank's points), rather than added at once
    virtual bool Pooled(void) const = 0;

    //! \brief Points other ranks sent (in one exchange); returns the first one's index among all the rank's points (when pooled)
    virtual std::size_t Receive(const std::vector<PointT> &points) = 0;
};

//! \brief A receiver of functions
template <typename PointT>
class FunctionGhostReceiver : public GhostReceiver<PointT>
{
public:
    using ReceiveFunction = std::function<std::size_t(const std::vector<PointT> &)>;

    FunctionGhostReceiver(bool pooled, ReceiveFunction receive) : pooled(pooled), receive(std::move(receive))
    {
    }

    bool Pooled(void) const override
    {
        return this->pooled;
    }

    std::size_t Receive(const std::vector<PointT> &points) override
    {
        return this->receive(points);
    }

private:
    bool pooled;
    ReceiveFunction receive;
};

namespace ghost_detail
{
    // The end of a domain's round: the points other ranks sent for it, then those it finds among the rank's points, then its mirrors
    template <typename PointT>
    void FinishGhostRound(GhostDomain<PointT> &domain)
    {
        GhostSearch<PointT> &search = *domain.search;
        GhostSearchMesh<PointT> &mesh = *domain.mesh;
        if(search.NothingToSearch())
        {
            return;
        }
        // big searches' points, then small ones' (the order decides the Delaunay of degenerate points)
        for(std::vector<std::size_t> &indices : search.TakeRemoteAnswers())
        {
            std::sort(indices.begin(), indices.end());
            std::vector<std::pair<std::size_t, PointT>> received;
            for(const std::size_t &pointIdxInAll : indices)
            {
                if(search.MarkBrought(pointIdxInAll))
                {
                    received.emplace_back(pointIdxInAll, PointT());
                }
            }
            mesh.AddRankPoints(received);
        }
        if(search.CanFindLocalPoints(mesh.RankPointsNum()))
        {
            mesh.AddRankPoints(search.SearchLocally(LocalSearch::Direct));
        }
        else
        {
            search.CheckNothingLocally(LocalSearch::Direct);
        }
        const std::vector<std::pair<std::size_t, PointT>> images = search.SearchLocally(LocalSearch::Images);
        if(not images.empty())
        {
            mesh.AddRankPoints(images);
        }
        std::vector<std::pair<std::size_t, std::size_t>> mirrorSources;
        const std::vector<PointT> mirrors = search.NewMirrors(mesh, mirrorSources);
        mesh.AddMirrors(mirrors, mirrorSources);
        search.EndRound(mesh);
    }

    #ifdef MADVORO_WITH_MPI
        // The domains' searches that reach other ranks (of one kind), in one exchange
        template <typename PointT, typename Query>
        void ExchangeGhostQueries(std::vector<GhostDomain<PointT>> &domains, GhostExchange<PointT> &exchange, GhostReceiver<PointT> &receiver)
        {
            constexpr bool big = std::is_same_v<Query, BigRangeQueryData<PointT>>;
            const auto queriesOf = [](const GhostSearch<PointT> &search) -> const std::vector<Query> &
            {
                if constexpr(big)
                {
                    return search.BigQueries();
                }
                else
                {
                    return search.SmallQueries();
                }
            };
            std::vector<Query> queries;
            std::vector<std::pair<std::size_t, std::size_t>> origins; // (domain, search)
            for(std::size_t d = 0; d < domains.size(); ++d)
            {
                const GhostSearch<PointT> &search = *domains[d].search;
                for(const std::size_t &j : big ? search.RemoteBig() : search.RemoteSmall())
                {
                    queries.push_back(queriesOf(search)[j]);
                    if(receiver.Pooled())
                    {
                        // a cell's index is unique among the domains as its index among all the rank's points
                        queries.back().pointIdx = domains[d].mesh->AllPointIndex(queries.back().pointIdx);
                    }
                    origins.emplace_back(d, j);
                }
            }
            QueryBatchInfo<Query, PointT> batchInfo = exchange.Exchange(queries);
            const std::size_t first = receiver.Receive(batchInfo.result);
            if(not receiver.Pooled())
            {
                for(const auto &answer : batchInfo.queriesAnswers)
                {
                    const auto &[d, j] = origins[answer.id];
                    domains[d].search->CountFound(big, queriesOf(*domains[d].search)[j].pointIdx, answer.finalResults.size());
                }
                return;
            }
            // each received point is in exactly one answer
            std::map<std::array<double, 3>, std::size_t> received;
            for(std::size_t i = 0; i < batchInfo.result.size(); ++i)
            {
                received.emplace(std::array<double, 3>{batchInfo.result[i].x, batchInfo.result[i].y, batchInfo.result[i].z}, first + i);
            }
            for(const auto &answer : batchInfo.queriesAnswers)
            {
                if(answer.finalResults.empty())
                {
                    continue;
                }
                std::vector<std::size_t> points;
                for(const PointT &point : answer.finalResults)
                {
                    points.push_back(received.at(std::array<double, 3>{point.x, point.y, point.z}));
                }
                const auto &[d, j] = origins[answer.id];
                domains[d].search->AddRemoteAnswer(big, j, std::move(points));
            }
        }
    #endif // MADVORO_WITH_MPI
} // namespace ghost_detail

/** \brief Brings the domains' ghost points until no domain, on any rank, searches; returns the rounds. `forEach(count, task)` runs a task for each domain (e.g. in parallel). `exchange` (or null) searches other ranks; `receiver` takes what they send. */
#ifdef MADVORO_WITH_MPI
    template <typename PointT, typename ForEach>
    std::size_t BringGhosts(std::vector<GhostDomain<PointT>> &domains, GhostExchange<PointT> *exchange, GhostReceiver<PointT> *receiver, const MPI_Comm &comm, const ForEach &forEach)
#else // MADVORO_WITH_MPI
    template <typename PointT, typename ForEach>
    std::size_t BringGhosts(std::vector<GhostDomain<PointT>> &domains, const ForEach &forEach)
#endif // MADVORO_WITH_MPI
{
    int size = 1;
    #ifdef MADVORO_WITH_MPI
        if(exchange != nullptr)
        {
            MPI_Comm_size(comm, &size);
        }
    #endif // MADVORO_WITH_MPI
    std::size_t rounds = 0;
    std::vector<char> searching(domains.size(), 0);
    // MADVORO_GHOST_TIMING=1: each round's phases (the slowest rank's and rank 0's) and searches, printed by rank 0
    static const bool timing = std::getenv("MADVORO_GHOST_TIMING") != nullptr;
    constexpr int PHASES = 5;
    std::vector<double> phaseTimes;
    std::vector<double> searches;
    auto mark = std::chrono::steady_clock::now();
    const auto lap = [&](void)
    {
        const auto now = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(now - mark).count();
        mark = now;
        return seconds;
    };
    const auto report = [&](void)
    {
        if(not timing)
        {
            return;
        }
        std::vector<double> slowest(phaseTimes.size()), total(searches.size());
        int rank = 0;
        #ifdef MADVORO_WITH_MPI
            if(exchange != nullptr)
            {
                MPI_Comm_rank(comm, &rank);
                MPI_Reduce(phaseTimes.data(), slowest.data(), static_cast<int>(phaseTimes.size()), MPI_DOUBLE, MPI_MAX, 0, comm);
                MPI_Reduce(searches.data(), total.data(), static_cast<int>(searches.size()), MPI_DOUBLE, MPI_SUM, 0, comm);
            }
            else
            {
                slowest = phaseTimes;
                total = searches;
            }
        #else // MADVORO_WITH_MPI
            slowest = phaseTimes;
            total = searches;
        #endif // MADVORO_WITH_MPI
        if(rank != 0)
        {
            return;
        }
        double sums[PHASES] = {0, 0, 0, 0, 0};
        for(std::size_t r = 0; r < rounds; ++r)
        {
            const double *t = &slowest[r * PHASES];
            const double *own = &phaseTimes[r * PHASES];
            std::printf("Ghost round %zu: slowest prepare %.4f big %.4f small %.4f finish %.4f wait %.4f | rank 0 %.4f %.4f %.4f %.4f %.4f | searches small %.0f (to other ranks %.0f) big %.0f (to other ranks %.0f)\n", r, t[0], t[1], t[2], t[3], t[4], own[0], own[1], own[2], own[3], own[4], total[r * 4], total[r * 4 + 1], total[r * 4 + 2], total[r * 4 + 3]);
            for(int k = 0; k < PHASES; ++k)
            {
                sums[k] += t[k];
            }
        }
        std::printf("Ghost rounds %zu: slowest prepare %.4f big %.4f small %.4f finish %.4f wait %.4f s\n", rounds, sums[0], sums[1], sums[2], sums[3], sums[4]);
    };
    while(true)
    {
        rounds++;
        lap();
        #ifdef MADVORO_WITH_MPI
            forEach(domains.size(), [&](std::size_t d) { searching[d] = domains[d].search->PrepareRound(*domains[d].mesh, exchange) ? 1 : 0; });
        #else // MADVORO_WITH_MPI
            forEach(domains.size(), [&](std::size_t d) { searching[d] = domains[d].search->PrepareRound(*domains[d].mesh) ? 1 : 0; });
        #endif // MADVORO_WITH_MPI
        double round[PHASES] = {lap(), 0, 0, 0, 0};
        if(timing)
        {
            double small = 0, smallRemote = 0, big = 0, bigRemote = 0;
            for(const GhostDomain<PointT> &domain : domains)
            {
                small += domain.search->SmallQueries().size();
                big += domain.search->BigQueries().size();
                #ifdef MADVORO_WITH_MPI
                    smallRemote += domain.search->RemoteSmall().size();
                    bigRemote += domain.search->RemoteBig().size();
                #endif // MADVORO_WITH_MPI
            }
            searches.insert(searches.end(), {small, smallRemote, big, bigRemote});
        }
        int finished = std::none_of(searching.begin(), searching.end(), [](char c) { return c != 0; }) ? 1 : 0;
        #ifdef MADVORO_WITH_MPI
            MPI_Request finishedReq;
            const int finishedHere = finished;
            if(exchange != nullptr)
            {
                MPI_Iallreduce(&finishedHere, &finished, 1, MPI_INT, MPI_SUM, comm, &finishedReq);
                ghost_detail::ExchangeGhostQueries<PointT, BigRangeQueryData<PointT>>(domains, *exchange, *receiver);
                round[1] = lap();
                ghost_detail::ExchangeGhostQueries<PointT, SmallRangeQueryData<PointT>>(domains, *exchange, *receiver);
                round[2] = lap();
            }
        #endif // MADVORO_WITH_MPI
        forEach(domains.size(), [&](std::size_t d) { ghost_detail::FinishGhostRound(domains[d]); });
        round[3] = lap();
        #ifdef MADVORO_WITH_MPI
            if(exchange != nullptr)
            {
                MPI_Wait(&finishedReq, MPI_STATUS_IGNORE);
            }
        #endif // MADVORO_WITH_MPI
        round[4] = lap();
        if(timing)
        {
            phaseTimes.insert(phaseTimes.end(), round, round + PHASES);
        }
        if(finished == size)
        {
            report();
            return rounds;
        }
    }
}

} // namespace MadVoro

#endif // MADVORO_GHOSTS_GHOST_BRINGER_HPP
