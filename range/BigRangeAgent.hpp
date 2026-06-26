#ifndef BIG_RANGE_AGENT_HPP
#define BIG_RANGE_AGENT_HPP

#include "finders/RangeFinder.hpp"
#include "finders/utils/IndexedVector.hpp"
#include <MeshDecomposer3D/environment/EnvironmentAgent.hpp>
#include <MeshDecomposer3D/environment/hilbert/HilbertTreeEnvAgent.hpp>
#ifdef RICH_MPI
    #include <mpi_utils/queryAgent/BusyWaitQueryAgent.hpp>
    #include <mpi_utils/queryAgent/WaitUntilAnsweredQueryAgent.hpp>
    #include <mpi_utils/queryAgent/BuffersManagerQueryAgent.hpp>
    #include <MeshDecomposer3D/environment/hilbert/DistributedOctEnvAgent.hpp> 
    #include "SentPointsContainer.hpp"
    #include <mpi_utils/serialize/Serializer.hpp>
#endif // RICH_MPI

#include "RangeQueryData.h"

struct BigRangeQueryData : public RangeQueryData
{
    Vector3D originalPoint;
    bool askOnlyClose; // in a case of a big query, we can ask all the ranks, or only the close ranks 

    friend inline std::ostream &operator<<(std::ostream &stream, const BigRangeQueryData &query)
    {
        return stream << "[BIG, point is " << query.originalPoint << ", sphere is (center = " << query.center << ", r = " << query.radius << ")]";
    }

    BigRangeQueryData(size_t pointIdx, const Vector3D &center, typename Vector3D::coord_type radius, const Vector3D &originalPoint, bool askOnlyClose): RangeQueryData(pointIdx, center, radius), originalPoint(originalPoint), askOnlyClose(askOnlyClose)
    {}
    
    BigRangeQueryData(): RangeQueryData(), originalPoint(Vector3D()), askOnlyClose(false)
    {}

    #ifdef RICH_MPI
        force_inline size_t dump(Serializer *serializer) const override
        {
            size_t bytes = 0;
            bytes += serializer->insert(this->pointIdx);
            bytes += serializer->insert(this->originalPoint);
            bytes += serializer->insert(this->center);
            bytes += serializer->insert(this->radius);
            bytes += serializer->insert(this->askOnlyClose);
            return bytes;
        }

        force_inline size_t load(const Serializer *serializer, std::size_t byteOffset) override
        {
            size_t bytes = 0;
            bytes += serializer->extract(this->pointIdx, byteOffset);
            bytes += serializer->extract(this->originalPoint, byteOffset + bytes);
            bytes += serializer->extract(this->center, byteOffset + bytes);
            bytes += serializer->extract(this->radius, byteOffset + bytes);
            bytes += serializer->extract(this->askOnlyClose, byteOffset + bytes);
            return bytes;
        }
    #endif // RICH_MPI
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
        #ifdef RICH_MPI
            : public AnswerAgent<BigRangeQueryData, Vector3D>
        #endif // RICH_MPI
    {
        friend class RangeAgent;

    public:
        #ifdef RICH_MPI
            BigRangeAnswerAgent(const RangeFinder<Vector3D> *rangeFinder, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): rangeFinder(rangeFinder), pointsContainer(pointsContainer)
        #else // RICH_MPI
            BigRangeAnswerAgent(const RangeFinder<Vector3D> *rangeFinder): rangeFinder(rangeFinder)
        #endif // RICH_MPI
        {}

        std::vector<size_t> selfAnswer(const BigRangeQueryData &query, std::unordered_set<size_t> &ignore)
        {
            // a big query, bring only the closest point
            std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(Vector3D(query.center.x, query.center.y, query.center.z), query.radius, Vector3D(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
            ignore.insert(indicesResult.begin(), indicesResult.end());
            return indicesResult;
        }

        #ifdef RICH_MPI
            std::vector<Vector3D> answer(const BigRangeQueryData &query, int _rank) override
            {
                const SentPointsContainer::PointsSet &ignore = this->pointsContainer.getSentDataSetRank(_rank);

                // a big query, bring only the closest point
                std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(Vector3D(query.center.x, query.center.y, query.center.z), query.radius, Vector3D(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
                indicesResult = this->pointsContainer.addPointsAsSent(_rank, indicesResult);

                std::vector<Vector3D> result;
                result.reserve(indicesResult.size());
                for(const size_t &pointIdx : indicesResult)
                {
                    result.push_back(this->rangeFinder->getPoint(pointIdx));
                }
                // std::cout << "answering to rank " << _rank << " " << result.size() << " points " << std::endl;            
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
        class BigRangeTalkAgent : public TalkAgent<BigRangeQueryData>
        {
        public:
            template<typename K, typename V>
            using _map = boost::container::flat_map<K, V>;

            BigRangeTalkAgent(const std::shared_ptr<EnvironmentAgent<Vector3D>> envAgent,         
                            #ifdef RICH_MPI
                                const MPI_Comm &comm = MPI_COMM_WORLD
                            #endif // RICH_MPI
                            ): envAgent(envAgent), supportsFurthestClosestRanks(false)
            {
                #ifdef RICH_MPI
                    MPI_Comm_rank(comm, &this->rank);
                    MPI_Comm_size(comm, &this->size);
                #else
                    this->rank = 0;
                    this->size = 1;
                #endif // RICH_MPI

                const DistributedOctEnvironmentAgent<Vector3D> *distribuedOctEnvAgent = dynamic_cast<const DistributedOctEnvironmentAgent<Vector3D>*>(this->envAgent.get());
                if(distribuedOctEnvAgent != nullptr)
                {
                    this->supportsFurthestClosestRanks = true;
                    this->getFurthestClosestRanks = [distribuedOctEnvAgent](const Vector3D &point){return distribuedOctEnvAgent->getClosestFurthestPointsByRanks(point);};
                }
                const HilbertTreeEnvironmentAgent<Vector3D> *hilbertTreeEnvAgent = dynamic_cast<const HilbertTreeEnvironmentAgent<Vector3D>*>(this->envAgent.get());
                if(hilbertTreeEnvAgent != nullptr)
                {
                    this->supportsFurthestClosestRanks = true;
                    this->getFurthestClosestRanks = [hilbertTreeEnvAgent](const Vector3D &point){return hilbertTreeEnvAgent->getClosestFurthestPointsByRanks(point);};
                }
            };

            inline EnvironmentAgent<Vector3D>::RanksSet getTalkList(const BigRangeQueryData &query) const override
            {
                if(std::isnan(query.center.x) or std::isnan(query.center.y) or std::isnan(query.center.z))
                {
                    UniversalError eo("In BigRangeTalkAgent, should not reach here, since the query center is NaN");
                    eo.addEntry("Query", query);
                    throw eo;
                }
                
                // std::cout << "rank " << this->rank << " calculates the talk list of query " << query << std::endl;
                EnvironmentAgent<Vector3D>::RanksSet intersectingRanks = this->envAgent->getIntersectingRanks(Vector3D(query.center.x, query.center.y, query.center.z), query.radius);
                if(intersectingRanks.empty())
                {
                    throw UniversalError("In range talk agent, should not reach here: the intersecting ranks list should at least contain the rank itself");
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
                    this->resultCache.insert({query.pointIdx, this->getFurthestClosestRanks(query.originalPoint)});
                    it = this->resultCache.find(query.pointIdx); // todo: can use previous line
                }
                HilbertCurveEnvironmentAgent<Vector3D>::DistancesVector &distances = (*it).second;
                
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
                    UniversalError eo("In BigRangeTalkAgent, should not reach here, since size of intersectingRanks is > 1");
                    eo.addEntry("Query", query);
                    eo.addEntry("minDistRank", minDistRank);
                    eo.addEntry("Size of intersectingRanks", intersectingRanks.size());
                    eo.addEntry("Distances", distances);
                    throw eo;
                }
                // consider the closest rank, and its furthest distance from the point, denoted as `closestDistThreshold`
                double closestDistThreshold = distances[minDistRank].second;

                // return all the ranks which their closest point to us is in distance of at most `closestDistThreshold`
                EnvironmentAgent<Vector3D>::RanksSet result;
                for(const int &_rank : intersectingRanks)
                {
                    if(distances[_rank].first <= (closestDistThreshold * (1 + EPSILON)))
                    {
                        result.insert(_rank);
                    }
                }

                if(result.size() <= 1)
                {
                    UniversalError eo("In BigRangeTalkAgent, should not reach here, since `result` must contain at least one additional rank");
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
            const std::shared_ptr<EnvironmentAgent<Vector3D>> envAgent;
            mutable _map<size_t, std::vector<std::pair<double, double>>> resultCache;
            int rank, size;
            bool supportsFurthestClosestRanks;
            #ifdef RICH_MPI
                std::function<HilbertCurveEnvironmentAgent<Vector3D>::DistancesVector(const Vector3D&)> getFurthestClosestRanks;
            #endif // RICH_MPI
        };
    #endif // RICH_MPI

public:
    template<typename T>
    using _set = std::unordered_set<T>;

    #ifdef RICH_MPI
        BigRangeAgent(const RangeFinder<Vector3D> *rangeFinder, const std::shared_ptr<EnvironmentAgent<Vector3D>> &envAgent, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): pointsContainer(pointsContainer)
    #else // RICH_MPI
        BigRangeAgent(const RangeFinder<Vector3D> *rangeFinder)
    #endif // RICH_MPI
    {
        #ifdef RICH_MPI
            this->ansAgent = new BigRangeAnswerAgent(rangeFinder, pointsContainer, comm);
            this->talkAgent = new BigRangeTalkAgent(envAgent, comm);
            this->queryAgent = new BuffersManagerQueryAgent<BigRangeQueryData, Vector3D>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
            // this->queryAgent = new BusyWaitQueryAgent<BigRangeQueryData, Vector3D>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
            //this->queryAgent = new WaitUntilAnsweredQueryAgent<BigRangeQueryData, Vector3D>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
        #else // RICH_MPI
            this->ansAgent = new BigRangeAnswerAgent(rangeFinder);
        #endif // RICH_MPI
    }

    ~BigRangeAgent()
    {
        #ifdef RICH_MPI
            delete this->queryAgent;
            delete this->talkAgent;
        #endif // RICH_MPI
        delete this->ansAgent;
    }

    #ifdef RICH_MPI
        inline QueryBatchInfo<BigRangeQueryData, Vector3D> runBatch(const std::vector<BigRangeQueryData> &queries)
        {
            return this->queryAgent->runBatch(queries);
        };
    #endif // RICH_MPI

    std::vector<std::vector<size_t>> selfBatchAnswer(const std::vector<BigRangeQueryData> &bigQueriesBatch, _set<size_t> &ignore)
    {
        std::vector<std::vector<size_t>> result;
        for(const BigRangeQueryData &query : bigQueriesBatch)
        {
            result.emplace_back(this->ansAgent->selfAnswer(query, ignore));
        }
        return result;
    }

    #ifdef RICH_MPI
        inline std::vector<std::vector<std::size_t>> &getSentPoints(){return this->pointsContainer.getSentData();};
        inline std::vector<std::vector<std::size_t>> &getRecvPoints(){return this->queryAgent->getRecvData();};
        inline std::vector<int> &getSentProc(){return this->pointsContainer.getSentProc();};
        inline std::vector<int> &getRecvProc(){return this->queryAgent->getRecvProc();};
    #endif // RICH_MPI

private:
    BigRangeAnswerAgent *ansAgent;
    #ifdef RICH_MPI
        BigRangeTalkAgent *talkAgent;
        QueryAgent<BigRangeQueryData, Vector3D> *queryAgent;
        SentPointsContainer &pointsContainer;
    #endif // RICH_MPI
};

#endif // BIG_RANGE_AGENT_HPP