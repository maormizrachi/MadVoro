// A rank's ghost searches of other ranks: what it sent and received, and to and from which ranks
#ifndef MADVORO_GHOSTS_GHOST_EXCHANGE_HPP
#define MADVORO_GHOSTS_GHOST_EXCHANGE_HPP

#ifdef MADVORO_WITH_MPI

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>
#include <mpi.h>

#include "range/finders/OctTreeFinder.hpp"
#include "range/SmallRangeAgent.hpp"
#include "range/BigRangeAgent.hpp"

namespace MadVoro
{

template <typename PointT>
class GhostExchange
{
public:
    // other ranks' searches are answered from the rank's points (in `tree`); the rank sent `alreadySent` and received from `alreadyReceived` before
    GhostExchange(OctTree<IndexedVector<PointT>> *tree, const std::vector<PointT> &points, const std::shared_ptr<EnvironmentAgent<PointT>> &envAgent, const SentPointsContainer &alreadySent, std::vector<int> alreadyReceived, const MPI_Comm &comm) : finder(tree, points), envAgent(envAgent), sent(alreadySent), alreadyReceived(std::move(alreadyReceived)), bigAgent(&this->finder, envAgent, sent, comm), smallAgent(&this->finder, envAgent, sent, comm)
    {
        MPI_Comm_rank(comm, &this->rank);
    }

    GhostExchange(const GhostExchange &) = delete;
    GhostExchange &operator=(const GhostExchange &) = delete;

    //! \brief Whether another rank's domain meets a search's sphere
    template <typename Query>
    bool ReachesOtherRanks(const Query &query) const
    {
        for(int otherRank : this->envAgent->getIntersectingRanks(PointT(query.center.x, query.center.y, query.center.z), query.radius))
        {
            if(otherRank != this->rank)
            {
                return true;
            }
        }
        return false;
    }

    //! \brief One round's searches, answered by the ranks they reach (and other ranks' searches answered)
    QueryBatchInfo<BigRangeQueryData<PointT>, PointT> Exchange(const std::vector<BigRangeQueryData<PointT>> &queries)
    {
        QueryBatchInfo<BigRangeQueryData<PointT>, PointT> batch = this->bigAgent.runBatch(queries);
        this->AddReceived(this->bigAgent.getRecvProc(), this->bigAgent.getRecvPoints(), batch.result.size());
        return batch;
    }

    QueryBatchInfo<SmallRangeQueryData<PointT>, PointT> Exchange(const std::vector<SmallRangeQueryData<PointT>> &queries)
    {
        QueryBatchInfo<SmallRangeQueryData<PointT>, PointT> batch = this->smallAgent.runBatch(queries);
        this->AddReceived(this->smallAgent.getRecvProc(), this->smallAgent.getRecvPoints(), batch.result.size());
        return batch;
    }

    //! \brief The ranks the exchanges received from, and the points of each (positions among all the received points, in the order they came)
    const std::vector<int> &ReceivedFrom(void) const
    {
        return this->receivedFrom;
    }

    const std::vector<std::vector<std::size_t>> &ReceivedPoints(void) const
    {
        return this->receivedPoints;
    }

    //! \brief The ranks the rank sent points to, and the points (indices among the rank's points)
    const std::vector<int> &SentTo(void) const
    {
        return this->sent.getSentProc();
    }

    const std::vector<std::vector<std::size_t>> &SentPoints(void) const
    {
        return this->sent.getSentData();
    }

    //! \brief Every rank the rank received from, before the exchanges and in them
    std::vector<std::vector<int>> AllReceivedFrom(void) const
    {
        return {this->alreadyReceived, this->receivedFrom};
    }

private:
    // a batch's received points (positions in its result), by the ranks that sent them
    void AddReceived(const std::vector<int> &from, const std::vector<std::vector<std::size_t>> &fromPoints, std::size_t batchSize)
    {
        for(std::size_t i = 0; i < from.size(); ++i)
        {
            const std::size_t r = std::find(this->receivedFrom.begin(), this->receivedFrom.end(), from[i]) - this->receivedFrom.begin();
            if(r == this->receivedFrom.size())
            {
                this->receivedFrom.push_back(from[i]);
                this->receivedPoints.emplace_back();
            }
            for(const std::size_t &point : fromPoints[i])
            {
                this->receivedPoints[r].push_back(this->receivedNum + point);
            }
        }
        this->receivedNum += batchSize;
    }

    OctTreeFinder<PointT> finder;
    std::shared_ptr<EnvironmentAgent<PointT>> envAgent;
    int rank = 0;
    SentPointsContainer sent;
    std::vector<int> alreadyReceived;
    std::vector<int> receivedFrom;
    std::vector<std::vector<std::size_t>> receivedPoints;
    std::size_t receivedNum = 0;
    BigRangeAgent<PointT> bigAgent;
    SmallRangeAgent<PointT> smallAgent;
};

} // namespace MadVoro

#endif // MADVORO_WITH_MPI

#endif // MADVORO_GHOSTS_GHOST_EXCHANGE_HPP
