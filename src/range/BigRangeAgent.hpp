#ifndef BIG_RANGE_AGENT_HPP
#define BIG_RANGE_AGENT_HPP

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
#include "utils/print/all.hpp"

#include "RangeQueryData.h"

namespace MadVoro
{
    struct BigRangeQueryData : public RangeQueryData
    {
        _3DPoint originalPoint;
        bool askOnlyClose; // in a case of a big query, we can ask all the ranks, or only the close ranks 

        friend inline std::ostream &operator<<(std::ostream &stream, const BigRangeQueryData &query)
        {
            return stream << "[BIG, point is " << query.originalPoint << ", sphere is (center = " << query.center << ", r = " << query.radius << ")]";
        }

        BigRangeQueryData(size_t pointIdx, const _3DPoint &center, typename _3DPoint::coord_type radius, const _3DPoint &originalPoint, bool askOnlyClose): RangeQueryData(pointIdx, center, radius), originalPoint(originalPoint), askOnlyClose(askOnlyClose)
        {}
        
        BigRangeQueryData(): RangeQueryData(), originalPoint(_3DPoint()), askOnlyClose(false)
        {}

        #ifdef MADVORO_WITH_MPI
            force_inline size_t dump(MPI::Serializer *serializer) const override
            {
                size_t bytes = 0;
                bytes += serializer->insert(this->pointIdx);
                bytes += serializer->insert(this->originalPoint);
                bytes += serializer->insert(this->center);
                bytes += serializer->insert(this->radius);
                bytes += serializer->insert(this->askOnlyClose);
                return bytes;
            }

            force_inline size_t load(const MPI::Serializer *serializer, std::size_t byteOffset) override
            {
                size_t bytes = 0;
                bytes += serializer->extract(this->pointIdx, byteOffset);
                bytes += serializer->extract(this->originalPoint, byteOffset + bytes);
                bytes += serializer->extract(this->center, byteOffset + bytes);
                bytes += serializer->extract(this->radius, byteOffset + bytes);
                bytes += serializer->extract(this->askOnlyClose, byteOffset + bytes);
                return bytes;
            }
        #endif // MADVORO_WITH_MPI
    };

    /**
     * The range agent is responsible for running batches of range queries. A batch is a collection of queries, and a range query is an instance of the `RangeQueryData` class, containing a point and a requested radius.
     * The range agent switches between roles - sending queries, receiving answers, and answering for incoming queries. It also supports duplications removal, and returns the results rearranged by processes (what are the points that were received from each one, and what points I sent to each one).
     * In order to answer for incoming requests, a range finder is required. A range finder is an object which holds a list of points, and can answer for range queries.
    */
    class BigRangeAgent
    {
    private:
        class BigRangeAnswerAgent
            #ifdef MADVORO_WITH_MPI
                : public QueryAgent::AnswerAgent<BigRangeQueryData, _3DPoint>
            #endif // MADVORO_WITH_MPI
        {
        public:
            #ifdef MADVORO_WITH_MPI
                BigRangeAnswerAgent(const Range::RangeFinder *rangeFinder, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): rangeFinder(rangeFinder), pointsContainer(pointsContainer)
            #else // MADVORO_WITH_MPI
                BigRangeAnswerAgent(const Range::RangeFinder *rangeFinder): rangeFinder(rangeFinder)
            #endif // MADVORO_WITH_MPI
            {}

            std::vector<size_t> selfAnswer(const BigRangeQueryData &query, boost::container::flat_set<size_t> &ignore)
            {
                // a big query, bring only the closest point
                std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(Point3D(query.center.x, query.center.y, query.center.z), query.radius, Point3D(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
                ignore.insert(indicesResult.begin(), indicesResult.end());
                return indicesResult;
            }

            #ifdef MADVORO_WITH_MPI
                std::vector<_3DPoint> answer(const BigRangeQueryData &query, int _rank) override
                {
                    const SentPointsContainer::PointsSet &ignore = this->pointsContainer.getSentDataSetRank(_rank);

                    // a big query, bring only the closest point
                    std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(Point3D(query.center.x, query.center.y, query.center.z), query.radius, Point3D(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
                    indicesResult = this->pointsContainer.addPointsAsSent(_rank, indicesResult);

                    std::vector<_3DPoint> result;
                    result.reserve(indicesResult.size());
                    for(const size_t &pointIdx : indicesResult)
                    {
                        result.push_back(_3DPoint(this->rangeFinder->getPoint(pointIdx)));
                    }
                    // std::cout << "answering to rank " << _rank << " " << result.size() << " points " << std::endl;            
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
            class BigRangeTalkAgent : public QueryAgent::TalkAgent<BigRangeQueryData>
            {
            public:
                template<typename K, typename V>
                using _map = boost::container::flat_map<K, V>;

                BigRangeTalkAgent(const EnvironmentAgent *envAgent,         
                                #ifdef MADVORO_WITH_MPI
                                    const MPI_Comm &comm = MPI_COMM_WORLD
                                #endif // MADVORO_WITH_MPI
                                ): envAgent(envAgent), supportsFurthestClosestRanks(false)
                {
                    #ifdef MADVORO_WITH_MPI
                        MPI_Initialized(&this->mpiInitialized);

                        if(this->mpiInitialized)
                        {
                            MPI_Comm_rank(comm, &this->rank);
                            MPI_Comm_size(comm, &this->size);
                        }
                        else
                        {
                            this->rank = 0;
                            this->size = 1;
                        }
                    #else
                        this->rank = 0;
                        this->size = 1;
                    #endif // MADVORO_WITH_MPI

                    const DistributedOctEnvironmentAgent *distribuedOctEnvAgent = dynamic_cast<const DistributedOctEnvironmentAgent*>(this->envAgent);
                    if(distribuedOctEnvAgent != nullptr)
                    {
                        this->supportsFurthestClosestRanks = true;
                        this->getFurthestClosestRanks = [distribuedOctEnvAgent](const _3DPoint &point){return distribuedOctEnvAgent->getClosestFurthestPointsByRanks(point);};
                    }
                    const HilbertTreeEnvironmentAgent *hilbertTreeEnvAgent = dynamic_cast<const HilbertTreeEnvironmentAgent*>(this->envAgent);
                    if(hilbertTreeEnvAgent != nullptr)
                    {
                        this->supportsFurthestClosestRanks = true;
                        this->getFurthestClosestRanks = [hilbertTreeEnvAgent](const _3DPoint &point){return hilbertTreeEnvAgent->getClosestFurthestPointsByRanks(point);};
                    }
                };

                inline EnvironmentAgent::RanksSet getTalkList(const BigRangeQueryData &query) const override
                {
                    if(not this->mpiInitialized)
                    {
                        return {this->rank};
                    }
                    if(std::isnan(query.center.x) or std::isnan(query.center.y) or std::isnan(query.center.z))
                    {
                        MadVoro::Exception::MadVoroException eo("In BigRangeTalkAgent, should not reach here, since the query center is NaN");
                        eo.addEntry("Query", query);
                        throw eo;
                    }

                    // std::cout << "rank " << this->rank << " calculates the talk list of query " << query << std::endl;
                    EnvironmentAgent::RanksSet intersectingRanks = this->envAgent->getIntersectingRanks(Point3D(query.center.x, query.center.y, query.center.z), query.radius);
                    if(intersectingRanks.empty())
                    {
                        throw MadVoro::Exception::MadVoroException("In range talk agent, should not reach here: the intersecting ranks list should at least contain the rank itself");
                    }

                    if(intersectingRanks.size() == 1)
                    {
                        return intersectingRanks;
                    }
                    
                    // check if has 'smartAgent' (an agent that can caluclate distances of ranks as well)
                    if(not this->supportsFurthestClosestRanks)
                    {
                        return intersectingRanks;
                    }
                    
                    // if the query requests to ask all the intersecting ranks, return all the intersecting ranks
                    if(not query.askOnlyClose)
                    {
                        return intersectingRanks; // ask all
                    }

                    // otherwise, the queries requests to ask only the close ranks
                    // we calculate the closest distances from the point, to all the other ranks.
                    // maybe the distances were already computed (check in a cache)
                    auto it = this->resultCache.find(query.pointIdx);
                    if(it == this->resultCache.end())
                    {
                        // not in cache, calculate it and insert to the cache
                        bool inserted;
                        std::tie(it, inserted) = this->resultCache.insert(std::make_pair(query.pointIdx, this->getFurthestClosestRanks(query.originalPoint)));
                        assert(inserted);
                        // it = this->resultCache.find(query.pointIdx); // todo: can use previous line
                    }
                    HilbertCurveEnvironmentAgent::DistancesVector &distances = (*it).second;
                    
                    // get the closest rank
                    double minDist = std::numeric_limits<double>::max();
                    int minDistRank = std::numeric_limits<int>::max();
                    for(const int &_rank : intersectingRanks)
                    {
                        if(_rank == this->rank)
                        {
                            continue; // don't count myself
                        }
                        if(distances[_rank].first < minDist)
                        {
                            minDist = distances[_rank].first;
                            minDistRank = _rank;
                        }
                    }
                    if(minDistRank >= this->size)
                    {
                        MadVoro::Exception::MadVoroException eo("In BigRangeTalkAgent, should not reach here, since size of intersectingRanks is > 1");
                        eo.addEntry("Query", query);
                        eo.addEntry("minDistRank", minDistRank);
                        eo.addEntry("Size of intersectingRanks", intersectingRanks.size());
                        eo.addEntry("Distances", distances);
                        throw eo;
                    }
                    // consider the closest rank, and its furthest distance from the point, denoted as `closestDistThreshold`
                    double closestDistThreshold = distances[minDistRank].second;

                    // return all the ranks which their closest point to us is in distance of at most `closestDistThreshold`
                    EnvironmentAgent::RanksSet result;
                    for(const int &_rank : intersectingRanks)
                    {
                        if(distances[_rank].first <= (closestDistThreshold * (1 + EPSILON)))
                        {
                            result.insert(_rank);
                        }
                    }

                    if(result.size() <= 1)
                    {
                        MadVoro::Exception::MadVoroException eo("In BigRangeTalkAgent, should not reach here, since `result` must contain at least one additional rank");
                        eo.addEntry("Query", query);
                        eo.addEntry("Distances", distances);
                        eo.addEntry("closestDistThreshold", closestDistThreshold);
                        eo.addEntry("Size of result", result.size());
                        eo.addEntry("Size of intersectingRanks", intersectingRanks.size());
                        throw eo;
                    }
                    return result;
                }

            private:
                const EnvironmentAgent *envAgent;
                mutable _map<size_t, std::vector<std::pair<double, double>>> resultCache;
                int mpiInitialized;
                int rank, size;
                bool supportsFurthestClosestRanks;
                #ifdef MADVORO_WITH_MPI
                    std::function<HilbertCurveEnvironmentAgent::DistancesVector(const _3DPoint&)> getFurthestClosestRanks;
                #endif // MADVORO_WITH_MPI
            };
        #endif // MADVORO_WITH_MPI

    public:
        template<typename T>
        using _set = boost::container::flat_set<T>;

        #ifdef MADVORO_WITH_MPI
            BigRangeAgent(const Range::RangeFinder *rangeFinder, const EnvironmentAgent *envAgent, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): pointsContainer(pointsContainer)
        #else // MADVORO_WITH_MPI
            BigRangeAgent(const Range::RangeFinder *rangeFinder)
        #endif // MADVORO_WITH_MPI
        {
            #ifdef MADVORO_WITH_MPI
                this->ansAgent = new BigRangeAnswerAgent(rangeFinder, pointsContainer, comm);
                this->talkAgent = new BigRangeTalkAgent(envAgent, comm);
                this->queryAgent = new QueryAgent::BusyWaitQueryAgent<BigRangeQueryData, _3DPoint>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
                //this->queryAgent = new WaitUntilAnsweredQueryAgent<BigRangeQueryData, _3DPoint>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
            #else // MADVORO_WITH_MPI
                this->ansAgent = new BigRangeAnswerAgent(rangeFinder);
            #endif // MADVORO_WITH_MPI
        }

        ~BigRangeAgent()
        {
            #ifdef MADVORO_WITH_MPI
                delete this->queryAgent;
                delete this->talkAgent;
            #endif // MADVORO_WITH_MPI
            delete this->ansAgent;
        }

        #ifdef MADVORO_WITH_MPI
            inline QueryAgent::QueryBatchInfo<BigRangeQueryData, _3DPoint> runBatch(const std::vector<BigRangeQueryData> &queries)
            {
                return this->queryAgent->runBatch(queries);
            };
        #endif // MADVORO_WITH_MPI

        std::vector<std::vector<size_t>> selfBatchAnswer(const std::vector<BigRangeQueryData> &bigQueriesBatch, boost::container::flat_set<size_t> &ignore)
        {
            std::vector<std::vector<size_t>> result;
            for(const BigRangeQueryData &query : bigQueriesBatch)
            {
                result.emplace_back(this->ansAgent->selfAnswer(query, ignore));
            }
            return result;
        }

        #ifdef MADVORO_WITH_MPI
            inline std::vector<std::vector<std::size_t>> &getSentPoints(){return this->pointsContainer.getSentData();};
            inline std::vector<std::vector<std::size_t>> &getRecvPoints(){return this->queryAgent->getRecvData();};
            inline std::vector<int> &getSentProc(){return this->pointsContainer.getSentProc();};
            inline std::vector<int> &getRecvProc(){return this->queryAgent->getRecvProc();};
        #endif // MADVORO_WITH_MPI

    private:
        BigRangeAnswerAgent *ansAgent;
        #ifdef MADVORO_WITH_MPI
            BigRangeTalkAgent *talkAgent;
            QueryAgent::QueryAgent<BigRangeQueryData, _3DPoint> *queryAgent;
            SentPointsContainer &pointsContainer;
        #endif // MADVORO_WITH_MPI
    };
}

#endif // BIG_RANGE_AGENT_HPP