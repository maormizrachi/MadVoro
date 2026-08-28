#ifndef BIG_RANGE_AGENT_HPP
#define BIG_RANGE_AGENT_HPP

#include "finders/RangeFinder.hpp"
#include "finders/utils/IndexedVector.hpp"
#include "exception/MadVoroException.hpp"
#ifdef MADVORO_WITH_MPI
    #include <MeshDecomposer3D/environment/EnvironmentAgent.hpp>
    #include <MeshDecomposer3D/environment/hilbert/HilbertTreeEnvAgent.hpp>
    #include <mpi_utils/queryAgent/BusyWaitQueryAgent.hpp>
    #include <mpi_utils/queryAgent/WaitUntilAnsweredQueryAgent.hpp>
    #include <mpi_utils/queryAgent/BuffersManagerQueryAgent.hpp>
    #include <MeshDecomposer3D/environment/hilbert/DistributedOctEnvAgent.hpp> 
    #include "SentPointsContainer.hpp"
    #include <mpi_utils/serialize/Serializer.hpp>
#endif // MADVORO_WITH_MPI

#include "RangeQueryData.h"

template <typename PointT>
struct BigRangeQueryData : public RangeQueryData<PointT>
{
    using coord_type = typename PointT::coord_type;

    PointT originalPoint;
    bool askOnlyClose;

    friend inline std::ostream &operator<<(std::ostream &stream, const BigRangeQueryData &query)
    {
        return stream << "[BIG, point is " << query.originalPoint << ", sphere is (center = " << query.center << ", r = " << query.radius << ")]";
    }

    BigRangeQueryData(size_t pointIdx, const PointT &center, coord_type radius, const PointT &originalPoint, bool askOnlyClose): RangeQueryData<PointT>(pointIdx, center, radius), originalPoint(originalPoint), askOnlyClose(askOnlyClose)
    {}
    
    BigRangeQueryData(): RangeQueryData<PointT>(), originalPoint(PointT()), askOnlyClose(false)
    {}

    #ifdef MADVORO_WITH_MPI
        force_inline size_t dump(Serializer *serializer) const override
        {
            size_t bytes = 0;
            bytes += serializer->insert(this->pointIdx);
            bytes += serializer->insert(this->originalPoint);
            bytes += serializer->insert(this->center);
            bytes += serializer->insert(this->radius);
            bytes += serializer->insert(this->askOnlyClose);
            bytes += serializer->insert(this->imageTranslation);
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
            bytes += serializer->extract(this->imageTranslation, byteOffset + bytes);
            return bytes;
        }
    #endif // MADVORO_WITH_MPI
};

/**
 * The range agent is responsible for running batches of range queries. A batch is a collection of queries, and a range query is an instance of the `RangeQueryData` class, containing a point and a requested radius.
 * The range agent switches between roles - sending queries, receiving answers, and answering for incoming queries. It also supports duplications removal, and returns the results rearranged by processes (what are the points that were received from each one, and what points I sent to each one).
 * In order to answer for incoming requests, a range finder is required. A range finder is an object which holds a list of points, and can answer for range queries.
*/
template <typename PointT>
class BigRangeAgent
{
    using coord_type = typename PointT::coord_type;

private:
    class BigRangeAnswerAgent
        #ifdef MADVORO_WITH_MPI
            : public AnswerAgent<BigRangeQueryData<PointT>, PointT>
        #endif // MADVORO_WITH_MPI
    {
        friend class RangeAgent;

    public:
        #ifdef MADVORO_WITH_MPI
            BigRangeAnswerAgent(const RangeFinder<PointT> *rangeFinder, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): rangeFinder(rangeFinder), pointsContainer(pointsContainer)
        #else // MADVORO_WITH_MPI
            BigRangeAnswerAgent(const RangeFinder<PointT> *rangeFinder): rangeFinder(rangeFinder)
        #endif // MADVORO_WITH_MPI
        {}

        std::vector<size_t> selfAnswer(const BigRangeQueryData<PointT> &query, std::array<std::unordered_set<size_t>, NUM_IMAGE_CODES> &ignoreSets)
        {
            int code = ImageCode(query.imageTranslation);
            std::unordered_set<size_t> &ignore = ignoreSets[code];
            std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(PointT(query.center.x, query.center.y, query.center.z), query.radius, PointT(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
            ignore.insert(indicesResult.begin(), indicesResult.end());
            return indicesResult;
        }

        #ifdef MADVORO_WITH_MPI
            std::vector<PointT> answer(const BigRangeQueryData<PointT> &query, int _rank) override
            {
                int code = ImageCode(query.imageTranslation);
                const SentPointsContainer::PointsSet &ignore = this->pointsContainer.getSentDataSetRank(_rank, code);

                std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(PointT(query.center.x, query.center.y, query.center.z), query.radius, PointT(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
                indicesResult = this->pointsContainer.addPointsAsSent(_rank, indicesResult, code);

                std::vector<PointT> result;
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
        class BigRangeTalkAgent : public TalkAgent<BigRangeQueryData<PointT>>
        {
        public:
            template<typename K, typename V>
            using _map = boost::container::flat_map<K, V>;

            BigRangeTalkAgent(const std::shared_ptr<EnvironmentAgent<PointT>> envAgent,         
                            #ifdef MADVORO_WITH_MPI
                                const MPI_Comm &comm = MPI_COMM_WORLD
                            #endif // MADVORO_WITH_MPI
                            ): envAgent(envAgent), supportsFurthestClosestRanks(false)
            {
                #ifdef MADVORO_WITH_MPI
                    MPI_Comm_rank(comm, &this->rank);
                    MPI_Comm_size(comm, &this->size);
                #else
                    this->rank = 0;
                    this->size = 1;
                #endif // MADVORO_WITH_MPI

                const DistributedOctEnvironmentAgent<PointT> *distribuedOctEnvAgent = dynamic_cast<const DistributedOctEnvironmentAgent<PointT>*>(this->envAgent.get());
                if(distribuedOctEnvAgent != nullptr)
                {
                    this->supportsFurthestClosestRanks = true;
                    this->getFurthestClosestRanks = [distribuedOctEnvAgent](const PointT &point){return distribuedOctEnvAgent->getClosestFurthestPointsByRanks(point);};
                }
                const HilbertTreeEnvironmentAgent<PointT> *hilbertTreeEnvAgent = dynamic_cast<const HilbertTreeEnvironmentAgent<PointT>*>(this->envAgent.get());
                if(hilbertTreeEnvAgent != nullptr)
                {
                    this->supportsFurthestClosestRanks = true;
                    this->getFurthestClosestRanks = [hilbertTreeEnvAgent](const PointT &point){return hilbertTreeEnvAgent->getClosestFurthestPointsByRanks(point);};
                }
            };

            inline typename EnvironmentAgent<PointT>::RanksSet getTalkList(const BigRangeQueryData<PointT> &query) const override
            {
                if(std::isnan(query.center.x) or std::isnan(query.center.y) or std::isnan(query.center.z))
                {
                    MadVoro::Exception::MadVoroException eo("In BigRangeTalkAgent, should not reach here, since the query center is NaN");
                    eo.addEntry("Query", query);
                    throw eo;
                }
                
                typename EnvironmentAgent<PointT>::RanksSet intersectingRanks = this->envAgent->getIntersectingRanks(PointT(query.center.x, query.center.y, query.center.z), query.radius);
                if(intersectingRanks.empty())
                {
                    throw MadVoro::Exception::MadVoroException("In range talk agent, should not reach here: the intersecting ranks list should at least contain the rank itself");
                }

                if(intersectingRanks.size() == 1)
                {
                    return intersectingRanks;
                }
                
                if(not this->supportsFurthestClosestRanks)
                {
                    return intersectingRanks;
                }
                
                if(not query.askOnlyClose)
                {
                    return intersectingRanks;
                }

                auto it = this->resultCache.find(query.pointIdx);
                if(it == this->resultCache.end())
                {
                    this->resultCache.emplace(query.pointIdx, this->getFurthestClosestRanks(query.originalPoint));
                    it = this->resultCache.find(query.pointIdx);
                }
                typename HilbertCurveEnvironmentAgent<PointT>::DistancesVector &distances = (*it).second;
                
                coord_type minDist = std::numeric_limits<coord_type>::max();
                int minDistRank = std::numeric_limits<int>::max();
                for(const int &_rank : intersectingRanks)
                {
                    if(_rank == this->rank)
                    {
                        continue;
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
                coord_type closestDistThreshold = distances[minDistRank].second;

                constexpr coord_type eps = static_cast<coord_type>(1e-12);
                typename EnvironmentAgent<PointT>::RanksSet result;
                for(const int &_rank : intersectingRanks)
                {
                    if(distances[_rank].first <= (closestDistThreshold * (1 + eps)))
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
            const std::shared_ptr<EnvironmentAgent<PointT>> envAgent;
            mutable _map<size_t, std::vector<std::pair<coord_type, coord_type>>> resultCache;
            int rank, size;
            bool supportsFurthestClosestRanks;
            #ifdef MADVORO_WITH_MPI
                std::function<typename HilbertCurveEnvironmentAgent<PointT>::DistancesVector(const PointT&)> getFurthestClosestRanks;
            #endif // MADVORO_WITH_MPI
        };
    #endif // MADVORO_WITH_MPI

public:
    template<typename T>
    using _set = std::unordered_set<T>;

    #ifdef MADVORO_WITH_MPI
        BigRangeAgent(const RangeFinder<PointT> *rangeFinder, const std::shared_ptr<EnvironmentAgent<PointT>> &envAgent, SentPointsContainer &pointsContainer, const MPI_Comm &comm = MPI_COMM_WORLD): pointsContainer(pointsContainer)
    #else // MADVORO_WITH_MPI
        BigRangeAgent(const RangeFinder<PointT> *rangeFinder)
    #endif // MADVORO_WITH_MPI
    {
        #ifdef MADVORO_WITH_MPI
            this->ansAgent = new BigRangeAnswerAgent(rangeFinder, pointsContainer, comm);
            this->talkAgent = new BigRangeTalkAgent(envAgent, comm);
            this->queryAgent = new BuffersManagerQueryAgent<BigRangeQueryData<PointT>, PointT>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
            // this->queryAgent = new BusyWaitQueryAgent<BigRangeQueryData<PointT>, PointT>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
            //this->queryAgent = new WaitUntilAnsweredQueryAgent<BigRangeQueryData<PointT>, PointT>(this->talkAgent, this->ansAgent, false /* dont send messages to self */, comm);
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
        inline QueryBatchInfo<BigRangeQueryData<PointT>, PointT> runBatch(const std::vector<BigRangeQueryData<PointT>> &queries)
        {
            return this->queryAgent->runBatch(queries);
        };
    #endif // MADVORO_WITH_MPI

    std::vector<std::vector<size_t>> selfBatchAnswer(const std::vector<BigRangeQueryData<PointT>> &bigQueriesBatch, std::array<std::unordered_set<size_t>, NUM_IMAGE_CODES> &ignoreSets)
    {
        std::vector<std::vector<size_t>> result;
        for(const BigRangeQueryData<PointT> &query : bigQueriesBatch)
        {
            result.emplace_back(this->ansAgent->selfAnswer(query, ignoreSets));
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
        QueryAgent<BigRangeQueryData<PointT>, PointT> *queryAgent;
        SentPointsContainer &pointsContainer;
    #endif // MADVORO_WITH_MPI
};

#endif // BIG_RANGE_AGENT_HPP
