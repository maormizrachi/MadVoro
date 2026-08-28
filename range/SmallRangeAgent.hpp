#ifndef SMALL_RANGE_AGENT_HPP
#define SMALL_RANGE_AGENT_HPP

#include "finders/RangeFinder.hpp"
#include "finders/utils/IndexedVector.hpp"
#include "exception/MadVoroException.hpp"
#ifdef MADVORO_WITH_MPI
    #include <MeshDecomposer3D/environment/EnvironmentAgent.hpp>
    #include <MeshDecomposer3D/environment/hilbert/HilbertTreeEnvAgent.hpp>
    #include <mpi_utils/queryAgent/BusyWaitQueryAgent.hpp>
    #include <mpi_utils/queryAgent/ThreePhasesQueryAgent.hpp>
    #include <mpi_utils/queryAgent/WaitUntilAnsweredQueryAgent.hpp>
    #include <mpi_utils/queryAgent/BuffersManagerQueryAgent.hpp>
    #include <mpi_utils/queryAgent/thread/ThreadsQueryAgent.hpp>
    #include <MeshDecomposer3D/environment/hilbert/DistributedOctEnvAgent.hpp> 
    #include <mpi_utils/serialize/Serializer.hpp>
    #include "SentPointsContainer.hpp"
#endif // MADVORO_WITH_MPI

#include "RangeQueryData.h"

template <typename PointT>
struct SmallRangeQueryData : public RangeQueryData<PointT>
{
    using coord_type = typename PointT::coord_type;

    size_t maxPointsToGet;

    friend inline std::ostream &operator<<(std::ostream &stream, const SmallRangeQueryData &query)
    {
        return stream << "[SMALL, max points is " << query.maxPointsToGet << ", sphere is (center = " << query.center << ", r = " << query.radius << ")]";
    }

    friend inline std::istream &operator>>(std::istream &stream, SmallRangeQueryData &query)
    {
        return stream >> query.maxPointsToGet >> query.center >> query.radius;
    }

    SmallRangeQueryData(size_t pointIdx, const PointT &center, coord_type radius, size_t maxPointsToGet):
        RangeQueryData<PointT>(pointIdx, center, radius), maxPointsToGet(maxPointsToGet)
    {};

    SmallRangeQueryData(): RangeQueryData<PointT>(), maxPointsToGet(0){};
    
    #ifdef MADVORO_WITH_MPI
        force_inline size_t dump(Serializer *serializer) const override
        {
            size_t bytes = 0;
            bytes += serializer->insert(this->pointIdx);
            bytes += serializer->insert(this->center);
            bytes += serializer->insert(this->radius);
            bytes += serializer->insert(this->imageTranslation);
            bytes += serializer->insert(this->maxPointsToGet);
            return bytes;
        }

        force_inline size_t load(const Serializer *serializer, std::size_t byteOffset) override
        {
            size_t bytes = 0;
            bytes += serializer->extract(this->pointIdx, byteOffset);
            bytes += serializer->extract(this->center, byteOffset + bytes);
            bytes += serializer->extract(this->radius, byteOffset + bytes);
            bytes += serializer->extract(this->imageTranslation, byteOffset + bytes);
            bytes += serializer->extract(this->maxPointsToGet, byteOffset + bytes);
            return bytes;
        }
    #endif // MADVORO_WITH_MPI
};

/**
 * The range agent is responsible for running batches of range queries. A batch is a collection of queries, and a range query is an instance of the `SmallRangeQueryData` class, containing a point and a requested radius.
 * The range agent switches between roles - sending queries, receiving answers, and answering for incoming queries. It also supports duplications removal, and returns the results rearranged by processes (what are the points that were received from each one, and what points I sent to each one).
 * In order to answer for incoming requests, a range finder is required. A range finder is an object which holds a list of points, and can answer for range queries.
*/
template <typename PointT>
class SmallRangeAgent
{
    using coord_type = typename PointT::coord_type;

private:
    class SmallRangeAnswerAgent
        #ifdef MADVORO_WITH_MPI
            : public AnswerAgent<SmallRangeQueryData<PointT>, PointT>
        #endif // MADVORO_WITH_MPI
    {
        friend class RangeAgent;

    public:
        #ifdef MADVORO_WITH_MPI
            SmallRangeAnswerAgent(const RangeFinder<PointT> *rangeFinder, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): rangeFinder(rangeFinder), pointsContainer(pointsContainer)
        #else // MADVORO_WITH_MPI
            SmallRangeAnswerAgent(const RangeFinder<PointT> *rangeFinder): rangeFinder(rangeFinder)
        #endif // MADVORO_WITH_MPI
        {}

        std::vector<size_t> selfAnswer(const SmallRangeQueryData<PointT> &query, std::array<std::unordered_set<size_t>, NUM_IMAGE_CODES> &ignoreSets)
        {
            int code = ImageCode(query.imageTranslation);
            std::unordered_set<size_t> &ignore = ignoreSets[code];
            std::vector<size_t> indicesResult = this->rangeFinder->range(PointT(query.center.x, query.center.y, query.center.z), query.radius, query.maxPointsToGet, ignore);
            ignore.insert(indicesResult.begin(), indicesResult.end());
            return indicesResult;
        }

        #ifdef MADVORO_WITH_MPI
            std::vector<PointT> answer(const SmallRangeQueryData<PointT> &query, int _rank) override
            {
                std::vector<PointT> result;
                std::vector<size_t> indicesResult;

                int code = ImageCode(query.imageTranslation);
                const SentPointsContainer::PointsSet &ignore = this->pointsContainer.getSentDataSetRank(_rank, code);

                indicesResult = this->rangeFinder->range(PointT(query.center.x, query.center.y, query.center.z), query.radius, query.maxPointsToGet, ignore);
                indicesResult = this->pointsContainer.addPointsAsSent(_rank, indicesResult, code);

                result.reserve(indicesResult.size());
                for(const size_t &pointIdx : indicesResult)
                {
                    result.push_back(this->rangeFinder->getPoint(pointIdx) + query.imageTranslation);
                }
                return result;
            }
        #endif // MADVORO_WITH_MPI

    private:
        const RangeFinder<PointT> *rangeFinder;
        #ifdef MADVORO_WITH_MPI
            SentPointsContainer &pointsContainer;
        #endif // MADVORO_WITH_MPI
    };
    
    #ifdef MADVORO_WITH_MPI
        class SmallRangeTalkAgent : public TalkAgent<SmallRangeQueryData<PointT>>
        {
        public:
            template<typename K, typename V>
            using _map = boost::container::flat_map<K, V>;

            SmallRangeTalkAgent(const std::shared_ptr<EnvironmentAgent<PointT>> envAgent,         
                            #ifdef MADVORO_WITH_MPI
                                const MPI_Comm &comm = MPI_COMM_WORLD
                            #endif // MADVORO_WITH_MPI
                            ): envAgent(envAgent)
            {
                #ifdef MADVORO_WITH_MPI
                    MPI_Comm_rank(comm, &this->rank);
                    MPI_Comm_size(comm, &this->size);
                #else
                    this->rank = 0;
                    this->size = 1;
                #endif // MADVORO_WITH_MPI
            };

            inline typename EnvironmentAgent<PointT>::RanksSet getTalkList(const SmallRangeQueryData<PointT> &query) const override
            {
                typename EnvironmentAgent<PointT>::RanksSet intersectingRanks = this->envAgent->getIntersectingRanks(PointT(query.center.x, query.center.y, query.center.z), query.radius);
                if(intersectingRanks.empty())
                {
                    throw MadVoro::Exception::MadVoroException("In range talk agent, should not reach here: the intersecting ranks list should at least contain the rank itself");
                }
                return intersectingRanks;
            }

        private:
            const std::shared_ptr<EnvironmentAgent<PointT>> envAgent;
            int rank, size;
        };
    #endif // MADVORO_WITH_MPI

public:
    template<typename T>
    using _set = std::unordered_set<T>;

    #ifdef MADVORO_WITH_MPI
        SmallRangeAgent(const RangeFinder<PointT> *rangeFinder, const std::shared_ptr<EnvironmentAgent<PointT>> &envAgent, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): pointsContainer(pointsContainer)
    #else // MADVORO_WITH_MPI
        SmallRangeAgent(const RangeFinder<PointT> *rangeFinder)
    #endif // MADVORO_WITH_MPI
    {
        #ifdef MADVORO_WITH_MPI
            this->ansAgent = new SmallRangeAnswerAgent(rangeFinder, pointsContainer, comm);
            this->talkAgent = new SmallRangeTalkAgent(envAgent, comm);
            this->queryAgent = new BuffersManagerQueryAgent<SmallRangeQueryData<PointT>, PointT>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
            // this->queryAgent = new BusyWaitQueryAgent<SmallRangeQueryData<PointT>, PointT>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
        #else // MADVORO_WITH_MPI
            this->ansAgent = new SmallRangeAnswerAgent(rangeFinder);
        #endif // MADVORO_WITH_MPI
    }

    ~SmallRangeAgent()
    {
        #ifdef MADVORO_WITH_MPI
            delete this->queryAgent;
            delete this->talkAgent;
        #endif // MADVORO_WITH_MPI
        delete this->ansAgent;
    }

    std::vector<std::vector<size_t>> selfBatchAnswer(const std::vector<SmallRangeQueryData<PointT>> &smallQueriesBatch, std::array<std::unordered_set<size_t>, NUM_IMAGE_CODES> &ignoreSets)
    {
        std::vector<std::vector<size_t>> result;
        for(const SmallRangeQueryData<PointT> &query : smallQueriesBatch)
        {
            result.emplace_back(this->ansAgent->selfAnswer(query, ignoreSets));
        }
        return result;
    }

    #ifdef MADVORO_WITH_MPI
        inline QueryBatchInfo<SmallRangeQueryData<PointT>, PointT> runBatch(const std::vector<SmallRangeQueryData<PointT>> &queries)
        {
            return this->queryAgent->runBatch(queries);
        };
    #endif // MADVORO_WITH_MPI

    #ifdef MADVORO_WITH_MPI
        inline std::vector<std::vector<std::size_t>> &getSentPoints(){return this->pointsContainer.getSentData();};
        inline std::vector<std::vector<std::size_t>> &getRecvPoints(){return this->queryAgent->getRecvData();};
        inline std::vector<int> &getSentProc(){return this->pointsContainer.getSentProc();};
        inline std::vector<int> &getRecvProc(){return this->queryAgent->getRecvProc();};
    #endif // MADVORO_WITH_MPI

private:
    SmallRangeAnswerAgent *ansAgent;
    #ifdef MADVORO_WITH_MPI
        QueryAgent<SmallRangeQueryData<PointT>, PointT> *queryAgent;
        SmallRangeTalkAgent *talkAgent;
        SentPointsContainer &pointsContainer;
    #endif // MADVORO_WITH_MPI
};

#endif // SMALL_RANGE_AGENT_HPP
