#ifndef SENT_CONTAINER_HPP
#define SENT_CONTAINER_HPP

#ifdef MADVORO_WITH_MPI

#include <vector>
#include <boost/container/flat_set.hpp>
#include <algorithm>
#include "exception/MadVoroException.hpp"

namespace MadVoro
{
    class SentPointsContainer
    {
    public:
        using PointsSet = boost::container::flat_set<size_t>;

        inline SentPointsContainer(const std::vector<int> &sentProc = std::vector<int>(), const std::vector<std::vector<size_t>> &sentData = std::vector<std::vector<size_t>>())
        {
            this->sentProc = sentProc;
            this->sentData = sentData;
            for(size_t i = 0; i < sentProc.size(); i++)
            {
                this->sentDataSet.emplace_back(PointsSet(sentData[i].begin(), sentData[i].end()));
            }
        }

        inline const std::vector<int> &getSentProc() const{return this->sentProc;};

        inline std::vector<int> &getSentProc(){return this->sentProc;};

        inline const std::vector<std::vector<size_t>> &getSentData() const{return this->sentData;};

        inline std::vector<std::vector<size_t>> &getSentData(){return this->sentData;};

        inline const std::vector<size_t> &getSentData(size_t index) const{return this->sentData[index];};

        inline std::vector<size_t> &getSentData(size_t index){return this->sentData[index];};

        inline const std::vector<PointsSet> &getSentDataSet() const{return this->sentDataSet;};

        inline std::vector<PointsSet> &getSentDataSet(){return this->sentDataSet;};

        inline const std::vector<size_t> &getSentDataRank(int rank) const
        {
            size_t index = this->findRankIndex(rank);
            if(index == this->sentProc.size())
            {
                return this->emptyVector;
            }
            return this->sentData[index];
        };

        inline const PointsSet &getSentDataSet(size_t index) const{return this->sentDataSet[index];};

        inline PointsSet &getSentDataSet(size_t index){return this->sentDataSet[index];};

        inline const PointsSet &getSentDataSetRank(int rank) const
        {
            size_t index = this->findRankIndex(rank);
            if(index == this->sentProc.size())
            {
                return this->emptySet;
            }
            return this->sentDataSet[index];
        };

        template<template<typename...> class Container, typename... Ts>
        inline Container<size_t> addPointsAsSent(int rank, const Container<size_t, Ts...> &points)
        {
            Container<size_t> result;
            if(points.empty())
            {
                return result; // `result` is empty
            }
            
            size_t rankIdx = this->findRankIndex(rank);
            if(rankIdx == this->sentProc.size())
            {
                // `_rank` is new
                this->initializeNewRank(rank);
            }

            for(const size_t &dataIdx : points)
            {
                if(this->sentDataSet[rankIdx].find(dataIdx) == this->sentDataSet[rankIdx].end())
                {
                    // `_data` was not sent before
                    result.push_back(dataIdx);
                    this->sentDataSet[rankIdx].insert(dataIdx);
                    this->sentData[rankIdx].push_back(dataIdx);
                }
            }
            return result;
        }

        inline std::vector<size_t> addPointAsSent(int rank, const size_t &point)
        {
            return this->addPointsAsSent(rank, std::vector<size_t>({point}));
        }

        std::vector<int> sentProc;
        std::vector<std::vector<size_t>> sentData;
        std::vector<PointsSet> sentDataSet;
        const std::vector<size_t> emptyVector = std::vector<size_t>();
        const PointsSet emptySet = PointsSet();

    private:

        inline size_t findRankIndex(int rank) const
        {
            return std::distance(this->sentProc.cbegin(), std::find(this->sentProc.cbegin(), this->sentProc.cend(), rank));
        }

        inline void initializeNewRank(int rank)
        {
            if(std::find(this->sentProc.begin(), this->sentProc.end(), rank) != this->sentProc.end())
            {
                MadVoro::Exception::MadVoroException eo("Rank is already in the SentPointsContainer");
                eo.addEntry("Rank", rank);
                throw eo;
            }
            this->sentProc.push_back(rank);
            this->sentData.emplace_back(std::vector<size_t>());
            this->sentDataSet.emplace_back(PointsSet());
        }
    };
}

#endif // MADVORO_WITH_MPI

#endif // SENT_CONTAINER_HPP