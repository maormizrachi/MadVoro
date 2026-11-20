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
#include "ClosestRanksMode.hpp"
#include <chrono>
#include <cstring>
#include <sstream>
#include <string>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

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

        std::vector<size_t> selfAnswer(const BigRangeQueryData<PointT> &query, std::unordered_set<size_t> &ignore)
        {
            std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(PointT(query.center.x, query.center.y, query.center.z), query.radius, PointT(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
            ignore.insert(indicesResult.begin(), indicesResult.end());
            return indicesResult;
        }

        #ifdef MADVORO_WITH_MPI
            std::vector<PointT> answer(const BigRangeQueryData<PointT> &query, int _rank) override
            {
                const SentPointsContainer::PointsSet &ignore = this->pointsContainer.getSentDataSetRank(_rank);

                std::vector<size_t> indicesResult = this->rangeFinder->closestPointInSphere(PointT(query.center.x, query.center.y, query.center.z), query.radius, PointT(query.originalPoint.x, query.originalPoint.y, query.originalPoint.z), ignore);
                indicesResult = this->pointsContainer.addPointsAsSent(_rank, indicesResult);

                std::vector<PointT> result;
                result.reserve(indicesResult.size());
                for(const size_t &pointIdx : indicesResult)
                {
                    result.push_back(this->rangeFinder->getPoint(pointIdx));
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
                            ): envAgent(envAgent), supportsFurthestClosestRanks(false), restrictedMode(MadVoro::GetClosestRanksRestrictedMode())
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
                    // The restricted lookup exists for the distributed oct tree only; the Hilbert tree agent keeps the
                    // dense scan in every mode.
                    this->octEnvAgent = distribuedOctEnvAgent;
                }
                const HilbertTreeEnvironmentAgent<PointT> *hilbertTreeEnvAgent = dynamic_cast<const HilbertTreeEnvironmentAgent<PointT>*>(this->envAgent.get());
                if(hilbertTreeEnvAgent != nullptr)
                {
                    this->supportsFurthestClosestRanks = true;
                    this->getFurthestClosestRanks = [hilbertTreeEnvAgent](const PointT &point){return hilbertTreeEnvAgent->getClosestFurthestPointsByRanks(point);};
                    this->octEnvAgent = nullptr;
                }
            };

            // Diagnostics: distinct points that entered the closest-rank lookup in this build (one per point, in every
            // mode; the count only grows).
            inline size_t GetClosestRanksCacheSize() const{return (this->usesRestricted()) ? this->sparseCache.size() : this->resultCache.size();}
            // Diagnostics of the restricted lookup (modes 1 and 2), cumulative over the agent's life: owner pairs
            // evaluated, leaves visited, seconds in the lookup excluding the index build (only while timing), seconds
            // building the per-owner leaf index, and the routing tree's leaf count.
            inline size_t GetOwnerEvaluations() const{return this->ownerEvaluations;}
            inline size_t GetLeafEvaluations() const{return this->leafEvaluations;}
            inline double GetDistanceSeconds() const{return this->distanceSeconds;}
            inline double GetIndexBuildSeconds() const{return this->indexBuildSeconds;}
            inline size_t GetIndexedLeafCount() const{return this->indexedLeafCount;}
            // Mode 2: owner pairs whose restricted and dense values differ in any bit, and the first such sample.
            inline size_t GetVerifyMismatches() const{return this->verifyMismatches;}
            inline const std::string &GetVerifyFirstMismatch() const{return this->verifyFirstMismatch;}
            inline int GetRestrictedMode() const{return this->restrictedMode;}
            inline void SetTiming(bool Timing){this->timing = Timing;}

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

                if(not this->usesRestricted())
                {
                    auto it = this->resultCache.find(query.pointIdx);
                    if(it == this->resultCache.end())
                    {
                        this->resultCache.emplace(query.pointIdx, this->getFurthestClosestRanks(query.originalPoint));
                        it = this->resultCache.find(query.pointIdx);
                    }
                    const typename HilbertCurveEnvironmentAgent<PointT>::DistancesVector &distances = (*it).second;
                    return this->selectCloseRanks(query, intersectingRanks, [&distances](int _rank) -> const std::pair<coord_type, coord_type>& {return distances[_rank];});
                }

                // Restricted (mode 1) or verify (mode 2). The entry of a point keeps the point of its first query, as
                // the dense cache does, and gains the owners of later queries of the same point (other tetrahedra).
                auto sit = this->sparseCache.find(query.pointIdx);
                if(sit == this->sparseCache.end())
                {
                    SparseEntry newEntry;
                    newEntry.point = query.originalPoint;
                    sit = this->sparseCache.emplace(query.pointIdx, std::move(newEntry)).first;
                }
                SparseEntry &entry = (*sit).second;
                this->fillOwners(entry, intersectingRanks);

                if(this->restrictedMode == 2)
                {
                    auto it = this->resultCache.find(query.pointIdx);
                    if(it == this->resultCache.end())
                    {
                        this->resultCache.emplace(query.pointIdx, this->getFurthestClosestRanks(entry.point));
                        it = this->resultCache.find(query.pointIdx);
                    }
                    const typename HilbertCurveEnvironmentAgent<PointT>::DistancesVector &distances = (*it).second;
                    for(const int &_rank : intersectingRanks)
                    {
                        auto checked = entry.verified.find(_rank);
                        if(checked != entry.verified.end())
                        {
                            continue;
                        }
                        entry.verified.insert(_rank);
                        const std::pair<coord_type, coord_type> &restricted = (*entry.owners.find(_rank)).second;
                        const std::pair<coord_type, coord_type> &dense = distances[_rank];
                        if(BigRangeTalkAgent::sameBits(restricted.first, dense.first) and BigRangeTalkAgent::sameBits(restricted.second, dense.second))
                        {
                            continue;
                        }
                        if(this->verifyMismatches == 0)
                        {
                            std::ostringstream sample;
                            sample.precision(17);
                            sample << "point " << query.pointIdx << " " << entry.point << " owner " << _rank << " restricted (" << restricted.first << ", " << restricted.second << ") dense (" << dense.first << ", " << dense.second << ")";
                            this->verifyFirstMismatch = sample.str();
                        }
                        ++this->verifyMismatches;
                    }
                    // Mode 2 selects with the dense values, so its ghost exchange is that of mode 0.
                    return this->selectCloseRanks(query, intersectingRanks, [&distances](int _rank) -> const std::pair<coord_type, coord_type>& {return distances[_rank];});
                }
                return this->selectCloseRanks(query, intersectingRanks, [&entry](int _rank) -> const std::pair<coord_type, coord_type>& {return (*entry.owners.find(_rank)).second;});
            }

        private:
            struct SparseEntry
            {
                PointT point;
                boost::container::flat_map<int, std::pair<coord_type, coord_type>> owners;
                boost::container::flat_set<int> verified;
            };

            inline bool usesRestricted() const{return this->restrictedMode != 0 and this->octEnvAgent != nullptr;}

            static inline bool sameBits(coord_type a, coord_type b)
            {
                return std::memcmp(&a, &b, sizeof(coord_type)) == 0;
            }

            // Evaluates the owners of `ranks` that the entry does not hold yet.
            inline void fillOwners(SparseEntry &entry, const typename EnvironmentAgent<PointT>::RanksSet &ranks) const
            {
                bool missing = false;
                for(const int &_rank : ranks)
                {
                    if(entry.owners.find(_rank) == entry.owners.end())
                    {
                        missing = true;
                        break;
                    }
                }
                if(not missing)
                {
                    return;
                }
                const auto &tree = this->octEnvAgent->getOctTree();
                if(not tree->hasOwnerLeafIndex())
                {
                    const auto indexStart = std::chrono::steady_clock::now();
                    tree->buildOwnerLeafIndex();
                    this->indexBuildSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - indexStart).count();
                }
                this->indexedLeafCount = tree->getIndexedLeafCount();
                std::chrono::steady_clock::time_point lookupStart;
                if(this->timing)
                {
                    lookupStart = std::chrono::steady_clock::now();
                }
                for(const int &_rank : ranks)
                {
                    if(entry.owners.find(_rank) != entry.owners.end())
                    {
                        continue;
                    }
                    entry.owners.emplace(_rank, tree->getClosestFurthestPointOfRank(entry.point, _rank, this->leafEvaluations));
                    ++this->ownerEvaluations;
                }
                if(this->timing)
                {
                    this->distanceSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - lookupStart).count();
                }
            }

            // The selection of the closest ranks: the nearest non-self rank sets the threshold (its furthest distance),
            // and every intersecting rank whose closest distance is within it is asked. `pairOf(rank)` returns the
            // (closest, furthest) squared distances of a rank in `intersectingRanks`.
            template<typename PairOf>
            inline typename EnvironmentAgent<PointT>::RanksSet selectCloseRanks(const BigRangeQueryData<PointT> &query, const typename EnvironmentAgent<PointT>::RanksSet &intersectingRanks, const PairOf &pairOf) const
            {
                coord_type minDist = std::numeric_limits<coord_type>::max();
                int minDistRank = std::numeric_limits<int>::max();
                for(const int &_rank : intersectingRanks)
                {
                    if(_rank == this->rank)
                    {
                        continue;
                    }
                    if(pairOf(_rank).first < minDist)
                    {
                        minDist = pairOf(_rank).first;
                        minDistRank = _rank;
                    }
                }
                if(minDistRank >= this->size)
                {
                    MadVoro::Exception::MadVoroException eo("In BigRangeTalkAgent, should not reach here, since size of intersectingRanks is > 1");
                    eo.addEntry("Query", query);
                    eo.addEntry("minDistRank", minDistRank);
                    eo.addEntry("Size of intersectingRanks", intersectingRanks.size());
                    eo.addEntry("Distances", this->describePairs(intersectingRanks, pairOf));
                    throw eo;
                }
                coord_type closestDistThreshold = pairOf(minDistRank).second;

                constexpr coord_type eps = static_cast<coord_type>(1e-12);
                typename EnvironmentAgent<PointT>::RanksSet result;
                for(const int &_rank : intersectingRanks)
                {
                    if(pairOf(_rank).first <= (closestDistThreshold * (1 + eps)))
                    {
                        result.insert(_rank);
                    }
                }

                if(result.size() <= 1)
                {
                    MadVoro::Exception::MadVoroException eo("In BigRangeTalkAgent, should not reach here, since `result` must contain at least one additional rank");
                    eo.addEntry("Query", query);
                    eo.addEntry("Distances", this->describePairs(intersectingRanks, pairOf));
                    eo.addEntry("closestDistThreshold", closestDistThreshold);
                    eo.addEntry("Size of result", result.size());
                    eo.addEntry("Size of intersectingRanks", intersectingRanks.size());
                    throw eo;
                }
                return result;
            }

            template<typename PairOf>
            inline std::string describePairs(const typename EnvironmentAgent<PointT>::RanksSet &ranks, const PairOf &pairOf) const
            {
                std::ostringstream text;
                text.precision(17);
                for(const int &_rank : ranks)
                {
                    text << "[" << _rank << ": " << pairOf(_rank).first << ", " << pairOf(_rank).second << "] ";
                }
                return text.str();
            }

            const std::shared_ptr<EnvironmentAgent<PointT>> envAgent;
            mutable _map<size_t, std::vector<std::pair<coord_type, coord_type>>> resultCache;
            mutable _map<size_t, SparseEntry> sparseCache;
            int rank, size;
            bool supportsFurthestClosestRanks;
            int restrictedMode;
            const DistributedOctEnvironmentAgent<PointT> *octEnvAgent = nullptr;
            bool timing = false;
            mutable size_t ownerEvaluations = 0;
            mutable size_t leafEvaluations = 0;
            mutable double distanceSeconds = 0;
            mutable double indexBuildSeconds = 0;
            mutable size_t indexedLeafCount = 0;
            mutable size_t verifyMismatches = 0;
            mutable std::string verifyFirstMismatch;
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

    std::vector<std::vector<size_t>> selfBatchAnswer(const std::vector<BigRangeQueryData<PointT>> &bigQueriesBatch, _set<size_t> &ignore)
    {
        std::vector<std::vector<size_t>> result;
        for(const BigRangeQueryData<PointT> &query : bigQueriesBatch)
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

        // Diagnostics only: per-batch record written by runBatch while non-null (see mpi_utils/BatchTrace.hpp).
        inline void SetBatchTrace(QueryBatchTrace *Trace)
        {
            this->BatchTrace = Trace;
            this->queryAgent->SetBatchTrace(Trace);
            this->talkAgent->SetTiming(Trace != nullptr);
        };
        inline const QueryBatchTrace *GetBatchTrace() const{return this->BatchTrace;};
        // Closest-rank computations done so far by the talk agent (its cache only grows): misses of the cache.
        inline size_t GetClosestRanksCacheSize() const{return this->talkAgent->GetClosestRanksCacheSize();};
        // Restricted closest-rank lookup diagnostics and the mode-2 verdict (see BigRangeTalkAgent).
        inline const BigRangeTalkAgent &GetTalkAgent() const{return *this->talkAgent;};
    #endif // MADVORO_WITH_MPI

private:
    BigRangeAnswerAgent *ansAgent;
    #ifdef MADVORO_WITH_MPI
        BigRangeTalkAgent *talkAgent;
        QueryAgent<BigRangeQueryData<PointT>, PointT> *queryAgent;
        SentPointsContainer &pointsContainer;
        QueryBatchTrace *BatchTrace = nullptr;
    #endif // MADVORO_WITH_MPI
};

#endif // BIG_RANGE_AGENT_HPP
