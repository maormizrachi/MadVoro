#ifndef SENT_CONTAINER_HPP
#define SENT_CONTAINER_HPP

#ifdef MADVORO_WITH_MPI

#include <vector>
#include <unordered_set>
#include <algorithm>
#include <array>
#include "../exception/MadVoroException.hpp"
#include "RangeQueryData.h"

class SentPointsContainer
{
public:
    using PointsSet = std::unordered_set<size_t>;
    using ImageSets = std::array<PointsSet, NUM_IMAGE_CODES>;

    inline SentPointsContainer(const std::vector<int> &sentProc = std::vector<int>(), const std::vector<std::vector<size_t>> &sentData = std::vector<std::vector<size_t>>())
    {
        this->sentProc = sentProc;
        this->sentData = sentData;
        for(size_t i = 0; i < sentProc.size(); i++)
        {
            this->sentDataSetByImage.emplace_back();
            this->sentDataSetByImage.back()[ZERO_IMAGE_CODE] = PointsSet(sentData[i].begin(), sentData[i].end());
        }
    }

    inline const std::vector<int> &getSentProc() const{return this->sentProc;};

    inline std::vector<int> &getSentProc(){return this->sentProc;};

    inline const std::vector<std::vector<size_t>> &getSentData() const{return this->sentData;};

    inline std::vector<std::vector<size_t>> &getSentData(){return this->sentData;};

    inline const std::vector<size_t> &getSentData(size_t index) const{return this->sentData[index];};

    inline std::vector<size_t> &getSentData(size_t index){return this->sentData[index];};

    inline const std::vector<ImageSets> &getSentDataSetByImage() const{return this->sentDataSetByImage;};

    inline std::vector<ImageSets> &getSentDataSetByImage(){return this->sentDataSetByImage;};

    inline const std::vector<size_t> &getSentDataRank(int rank) const
    {
        size_t index = this->findRankIndex(rank);
        if(index == this->sentProc.size())
        {
            return this->emptyVector;
        }
        return this->sentData[index];
    };

    inline const PointsSet &getSentDataSet(size_t index) const{return this->sentDataSetByImage[index][ZERO_IMAGE_CODE];};

    inline PointsSet &getSentDataSet(size_t index){return this->sentDataSetByImage[index][ZERO_IMAGE_CODE];};

    inline const PointsSet &getSentDataSetRank(int rank) const
    {
        return this->getSentDataSetRank(rank, ZERO_IMAGE_CODE);
    };

    inline const PointsSet &getSentDataSetRank(int rank, int imageCode) const
    {
        size_t index = this->findRankIndex(rank);
        if(index == this->sentProc.size())
        {
            return this->emptySet;
        }
        return this->sentDataSetByImage[index][imageCode];
    };

    template<template<typename...> class Container, typename... Ts>
    inline Container<size_t> addPointsAsSent(int rank, const Container<size_t, Ts...> &points)
    {
        return this->addPointsAsSent(rank, points, ZERO_IMAGE_CODE);
    }

    template<template<typename...> class Container, typename... Ts>
    inline Container<size_t> addPointsAsSent(int rank, const Container<size_t, Ts...> &points, int imageCode)
    {
        Container<size_t> result;
        if(points.empty())
        {
            return result;
        }
        
        size_t rankIdx = this->findRankIndex(rank);
        if(rankIdx == this->sentProc.size())
        {
            this->initializeNewRank(rank);
        }

        PointsSet &ignoreSet = this->sentDataSetByImage[rankIdx][imageCode];
        for(const size_t &dataIdx : points)
        {
            if(ignoreSet.find(dataIdx) == ignoreSet.end())
            {
                result.push_back(dataIdx);
                ignoreSet.insert(dataIdx);
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
    std::vector<ImageSets> sentDataSetByImage;
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
            throw MadVoro::Exception::MadVoroException("Rank is already in the SentPointsContainer");
        }
        this->sentProc.push_back(rank);
        this->sentData.emplace_back(std::vector<size_t>());
        this->sentDataSetByImage.emplace_back();
    }
};

#endif // MADVORO_WITH_MPI

#endif // SENT_CONTAINER_HPP
