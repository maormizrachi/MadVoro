#ifndef SMALL_RANGE_AGENT_HPP
#define SMALL_RANGE_AGENT_HPP

#include "range/finders/RangeFinder.hpp"
#include "range/finders/utils/IndexedVector.hpp"
#include "environment/EnvironmentAgent.h"
#include "environment/hilbert/HilbertTreeEnvAgent.hpp"
#ifdef MADVORO_WITH_MPI
    #include "queryAgent/BusyWaitQueryAgent.hpp"
    #include "queryAgent/WaitUntilAnsweredQueryAgent.hpp"
    #include "environment/hilbert/DistributedOctEnvAgent.hpp" 
    #include "SentPointsContainer.hpp"
    #include "mpi/serialize/Serializer.hpp"
#endif // MADVORO_WITH_MPI

#include "RangeQueryData.h"

namespace MadVoro
{
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

        SmallRangeQueryData(size_t pointIdx, const _3DPoint &center, typename _3DPoint::coord_type radius, size_t maxPointsToGet):
            RangeQueryData(pointIdx, center, radius), maxPointsToGet(maxPointsToGet)
        {};

        SmallRangeQueryData(): RangeQueryData(), maxPointsToGet(0){};
        
        #ifdef MADVORO_WITH_MPI
            force_inline size_t dump(MPI::Serializer *serializer) const override
            {
                size_t bytes = 0;
                bytes += serializer->insert(this->pointIdx);
                bytes += serializer->insert(this->center);
                bytes += serializer->insert(this->radius);
                bytes += serializer->insert(this->maxPointsToGet);
                return bytes;
            }

            force_inline size_t load(const MPI::Serializer *serializer, std::size_t byteOffset) override
            {
                size_t bytes = 0;
                bytes += serializer->extract(this->pointIdx, byteOffset);
                bytes += serializer->extract(this->center, byteOffset + bytes);
                bytes += serializer->extract(this->radius, byteOffset + bytes);
                bytes += serializer->extract(this->maxPointsToGet, byteOffset + bytes);
                return bytes;
            }
        #endif // MADVORO_WITH_MPI

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
            #ifdef MADVORO_WITH_MPI
                : public QueryAgent::AnswerAgent<SmallRangeQueryData, _3DPoint>
            #endif // MADVORO_WITH_MPI
        {
        public:
            #ifdef MADVORO_WITH_MPI
                SmallRangeAnswerAgent(const Range::RangeFinder *rangeFinder, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): rangeFinder(rangeFinder), pointsContainer(pointsContainer)
            #else // MADVORO_WITH_MPI
                SmallRangeAnswerAgent(const Range::RangeFinder *rangeFinder): rangeFinder(rangeFinder)
            #endif // MADVORO_WITH_MPI
            {}

            std::vector<size_t> selfAnswer(const SmallRangeQueryData &query, boost::container::flat_set<size_t> &ignore)
            {
                // a small query, bring the requested number of points
                std::vector<size_t> indicesResult = this->rangeFinder->range(Point3D(query.center.x, query.center.y, query.center.z), query.radius, query.maxPointsToGet, ignore);
                ignore.insert(indicesResult.begin(), indicesResult.end());
                return indicesResult;
            }

            #ifdef MADVORO_WITH_MPI
                std::vector<_3DPoint> answer(const SmallRangeQueryData &query, int _rank) override
                {
                    std::vector<_3DPoint> result;
                    std::vector<size_t> indicesResult;

                    const SentPointsContainer::PointsSet &ignore = this->pointsContainer.getSentDataSetRank(_rank);

                    // a small query, bring the requested number of points
                    indicesResult = this->rangeFinder->range(Point3D(query.center.x, query.center.y, query.center.z), query.radius, query.maxPointsToGet, ignore);
                    indicesResult = this->pointsContainer.addPointsAsSent(_rank, indicesResult);

                    result.reserve(indicesResult.size());
                    for(const size_t &pointIdx : indicesResult)
                    {
                        result.push_back(_3DPoint(this->rangeFinder->getPoint(pointIdx)));
                    }
                    return result;
                }
            #endif // MADVORO_WITH_MPI

        private:
            const Range::RangeFinder *rangeFinder;
            #ifdef MADVORO_WITH_MPI
                SentPointsContainer &pointsContainer;
            #endif // MADVORO_WITH_MPI
        };
        
        #ifdef MADVORO_WITH_MPI
            class SmallRangeTalkAgent : public QueryAgent::TalkAgent<SmallRangeQueryData>
            {
            public:
                template<typename K, typename V>
                using _map = boost::container::flat_map<K, V>;

                SmallRangeTalkAgent(const EnvironmentAgent *envAgent,         
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

                inline EnvironmentAgent::RanksSet getTalkList(const SmallRangeQueryData &query) const override
                {
                    // check if has 'smartAgent' (an agent that can caluclate distances of ranks as well)
                    EnvironmentAgent::RanksSet intersectingRanks = this->envAgent->getIntersectingRanks(Point3D(query.center.x, query.center.y, query.center.z), query.radius);
                    if(intersectingRanks.empty())
                    {
                        throw MadVoro::Exception::MadVoroException("In range talk agent, should not reach here: the intersecting ranks list should at least contain the rank itself");
                    }
                    return intersectingRanks;
                }

            private:
                const EnvironmentAgent *envAgent;
                int rank, size;
            };
        #endif // MADVORO_WITH_MPI

    public:
        template<typename T>
        using _set = boost::container::flat_set<T>;

        #ifdef MADVORO_WITH_MPI
            SmallRangeAgent(const Range::RangeFinder *rangeFinder, const EnvironmentAgent *envAgent, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): pointsContainer(pointsContainer)
        #else // MADVORO_WITH_MPI
            SmallRangeAgent(const Range::RangeFinder *rangeFinder)
        #endif // MADVORO_WITH_MPI
        {
            #ifdef MADVORO_WITH_MPI
                this->ansAgent = new SmallRangeAnswerAgent(rangeFinder, pointsContainer, comm);
                this->talkAgent = new SmallRangeTalkAgent(envAgent, comm);
                this->queryAgent = new QueryAgent::BusyWaitQueryAgent<SmallRangeQueryData, _3DPoint>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
                //this->queryAgent = new WaitUntilAnsweredQueryAgent<SmallRangeQueryData, _3DPoint>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
                // this->queryAgent = new ThreadsQueryAgent<SmallRangeQueryData, _3DPoint>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
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

        std::vector<std::vector<size_t>> selfBatchAnswer(const std::vector<SmallRangeQueryData> &smallQueriesBatch, boost::container::flat_set<size_t> &ignore)
        {
            std::vector<std::vector<size_t>> result;
            for(const SmallRangeQueryData &query : smallQueriesBatch)
            {
                result.emplace_back(this->ansAgent->selfAnswer(query, ignore));
            }
            return result;
        }

        #ifdef MADVORO_WITH_MPI
            inline QueryAgent::QueryBatchInfo<SmallRangeQueryData, _3DPoint> runBatch(const std::vector<SmallRangeQueryData> &queries)
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
            QueryAgent::QueryAgent<SmallRangeQueryData, _3DPoint> *queryAgent;
            SmallRangeTalkAgent *talkAgent;
            SentPointsContainer &pointsContainer;
        #endif // MADVORO_WITH_MPI
    };
}

#endif // SMALL_RANGE_AGENT_HPP