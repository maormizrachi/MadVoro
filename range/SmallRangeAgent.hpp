#ifndef SMALL_RANGE_AGENT_HPP
#define SMALL_RANGE_AGENT_HPP

#include "finders/RangeFinder.hpp"
#include "finders/utils/IndexedVector.hpp"
#include <MeshDecomposer3D/environment/EnvironmentAgent.hpp>
#include <MeshDecomposer3D/environment/hilbert/HilbertTreeEnvAgent.hpp>
#ifdef RICH_MPI
    #include <mpi_utils/queryAgent/BusyWaitQueryAgent.hpp>
    #include <mpi_utils/queryAgent/ThreePhasesQueryAgent.hpp>
    #include <mpi_utils/queryAgent/WaitUntilAnsweredQueryAgent.hpp>
    #include <mpi_utils/queryAgent/BuffersManagerQueryAgent.hpp>
    #include <mpi_utils/queryAgent/thread/ThreadsQueryAgent.hpp>
    #include <MeshDecomposer3D/environment/hilbert/DistributedOctEnvAgent.hpp> 
    #include <mpi_utils/serialize/Serializer.hpp>
    #include "SentPointsContainer.hpp"
#endif // RICH_MPI

#include "RangeQueryData.h"

typedef struct SmallRangeQueryData : public RangeQueryData
{
    size_t maxPointsToGet;

    friend inline std::ostream &operator<<(std::ostream &stream, const SmallRangeQueryData &query)
    {
        return stream << "[SMALL, max points is " << query.maxPointsToGet << ", sphere is (center = " << query.center << ", r = " << query.radius << ")]";
    }

    friend inline std::istream &operator>>(std::istream &stream, SmallRangeQueryData &query)
    {
        return stream >> query.maxPointsToGet >> query.center >> query.radius;
    }

    SmallRangeQueryData(size_t pointIdx, const Vector3D &center, typename Vector3D::coord_type radius, size_t maxPointsToGet):
        RangeQueryData(pointIdx, center, radius), maxPointsToGet(maxPointsToGet)
    {};

    SmallRangeQueryData(): RangeQueryData(), maxPointsToGet(0){};
    
    #ifdef RICH_MPI
        force_inline size_t dump(Serializer *serializer) const override
        {
            size_t bytes = 0;
            bytes += serializer->insert(this->pointIdx);
            bytes += serializer->insert(this->center);
            bytes += serializer->insert(this->radius);
            bytes += serializer->insert(this->maxPointsToGet);
            return bytes;
        }

        force_inline size_t load(const Serializer *serializer, std::size_t byteOffset) override
        {
            size_t bytes = 0;
            bytes += serializer->extract(this->pointIdx, byteOffset);
            bytes += serializer->extract(this->center, byteOffset + bytes);
            bytes += serializer->extract(this->radius, byteOffset + bytes);
            bytes += serializer->extract(this->maxPointsToGet, byteOffset + bytes);
            return bytes;
        }
    #endif // RICH_MPI

} SmallRangeQueryData;

/**
 * The range agent is responsible for running batches of range queries. A batch is a collection of queries, and a range query is an instance of the `SmallRangeQueryData` class, containing a point and a requested radius.
 * The range agent switches between roles - sending queries, receiving answers, and answering for incoming queries. It also supports duplications removal, and returns the results rearranged by processes (what are the points that were received from each one, and what points I sent to each one).
 * In order to answer for incoming requests, a range finder is required. A range finder is an object which holds a list of points, and can answer for range queries.
*/
class SmallRangeAgent
{
private:
    class SmallRangeAnswerAgent
        #ifdef RICH_MPI
            : public AnswerAgent<SmallRangeQueryData, Vector3D>
        #endif // RICH_MPI
    {
        friend class RangeAgent;

    public:
        #ifdef RICH_MPI
            SmallRangeAnswerAgent(const RangeFinder<Vector3D> *rangeFinder, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): rangeFinder(rangeFinder), pointsContainer(pointsContainer)
        #else // RICH_MPI
            SmallRangeAnswerAgent(const RangeFinder<Vector3D> *rangeFinder): rangeFinder(rangeFinder)
        #endif // RICH_MPI
        {}

        std::vector<size_t> selfAnswer(const SmallRangeQueryData &query, std::unordered_set<size_t> &ignore)
        {
            // a small query, bring the requested number of points
            std::vector<size_t> indicesResult = this->rangeFinder->range(Vector3D(query.center.x, query.center.y, query.center.z), query.radius, query.maxPointsToGet, ignore);
            ignore.insert(indicesResult.begin(), indicesResult.end());
            return indicesResult;
        }

        #ifdef RICH_MPI
            std::vector<Vector3D> answer(const SmallRangeQueryData &query, int _rank) override
            {
                std::vector<Vector3D> result;
                std::vector<size_t> indicesResult;

                const SentPointsContainer::PointsSet &ignore = this->pointsContainer.getSentDataSetRank(_rank);

                // a small query, bring the requested number of points
                indicesResult = this->rangeFinder->range(Vector3D(query.center.x, query.center.y, query.center.z), query.radius, query.maxPointsToGet, ignore);
                indicesResult = this->pointsContainer.addPointsAsSent(_rank, indicesResult);

                result.reserve(indicesResult.size());
                for(const size_t &pointIdx : indicesResult)
                {
                    result.push_back(this->rangeFinder->getPoint(pointIdx));
                }
                return result;
            }
        #endif // RICH_MPI

    private:
        const RangeFinder<Vector3D> *rangeFinder;
        #ifdef RICH_MPI
            SentPointsContainer &pointsContainer;
        #endif // RICH_MPI
    };
    
    #ifdef RICH_MPI
        class SmallRangeTalkAgent : public TalkAgent<SmallRangeQueryData>
        {
        public:
            template<typename K, typename V>
            using _map = boost::container::flat_map<K, V>;

            SmallRangeTalkAgent(const std::shared_ptr<EnvironmentAgent<Vector3D>> envAgent,         
                            #ifdef RICH_MPI
                                const MPI_Comm &comm = MPI_COMM_WORLD
                            #endif // RICH_MPI
                            ): envAgent(envAgent)
            {
                #ifdef RICH_MPI
                    MPI_Comm_rank(comm, &this->rank);
                    MPI_Comm_size(comm, &this->size);
                #else
                    this->rank = 0;
                    this->size = 1;
                #endif // RICH_MPI
            };

            inline EnvironmentAgent<Vector3D>::RanksSet getTalkList(const SmallRangeQueryData &query) const override
            {
                // check if has 'smartAgent' (an agent that can caluclate distances of ranks as well)
                EnvironmentAgent<Vector3D>::RanksSet intersectingRanks = this->envAgent->getIntersectingRanks(Vector3D(query.center.x, query.center.y, query.center.z), query.radius);
                if(intersectingRanks.empty())
                {
                    throw UniversalError("In range talk agent, should not reach here: the intersecting ranks list should at least contain the rank itself");
                }
                return intersectingRanks;
            }

        private:
            const std::shared_ptr<EnvironmentAgent<Vector3D>> envAgent;
            int rank, size;
        };
    #endif // RICH_MPI

public:
    template<typename T>
    using _set = std::unordered_set<T>;

    #ifdef RICH_MPI
        SmallRangeAgent(const RangeFinder<Vector3D> *rangeFinder, const std::shared_ptr<EnvironmentAgent<Vector3D>> &envAgent, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): pointsContainer(pointsContainer)
    #else // RICH_MPI
        SmallRangeAgent(const RangeFinder<Vector3D> *rangeFinder)
    #endif // RICH_MPI
    {
        #ifdef RICH_MPI
            this->ansAgent = new SmallRangeAnswerAgent(rangeFinder, pointsContainer, comm);
            this->talkAgent = new SmallRangeTalkAgent(envAgent, comm);
            this->queryAgent = new BuffersManagerQueryAgent<SmallRangeQueryData, Vector3D>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
            // this->queryAgent = new BusyWaitQueryAgent<SmallRangeQueryData, Vector3D>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
        #else // RICH_MPI
            this->ansAgent = new SmallRangeAnswerAgent(rangeFinder);
        #endif // RICH_MPI
    }

    ~SmallRangeAgent()
    {
        #ifdef RICH_MPI
            delete this->queryAgent;
            delete this->talkAgent;
        #endif // RICH_MPI
        delete this->ansAgent;
    }

    std::vector<std::vector<size_t>> selfBatchAnswer(const std::vector<SmallRangeQueryData> &smallQueriesBatch, _set<size_t> &ignore)
    {
        std::vector<std::vector<size_t>> result;
        for(const SmallRangeQueryData &query : smallQueriesBatch)
        {
            result.emplace_back(this->ansAgent->selfAnswer(query, ignore));
        }
        return result;
    }

    #ifdef RICH_MPI
        inline QueryBatchInfo<SmallRangeQueryData, Vector3D> runBatch(const std::vector<SmallRangeQueryData> &queries)
        {
            return this->queryAgent->runBatch(queries);
        };
    #endif // RICH_MPI

    #ifdef RICH_MPI
        inline std::vector<std::vector<std::size_t>> &getSentPoints(){return this->pointsContainer.getSentData();};
        inline std::vector<std::vector<std::size_t>> &getRecvPoints(){return this->queryAgent->getRecvData();};
        inline std::vector<int> &getSentProc(){return this->pointsContainer.getSentProc();};
        inline std::vector<int> &getRecvProc(){return this->queryAgent->getRecvProc();};
    #endif // RICH_MPI

private:
    SmallRangeAnswerAgent *ansAgent;
    #ifdef RICH_MPI
        QueryAgent<SmallRangeQueryData, Vector3D> *queryAgent;
        SmallRangeTalkAgent *talkAgent;
        SentPointsContainer &pointsContainer;
    #endif // RICH_MPI
};

#endif // SMALL_RANGE_AGENT_HPP