#ifndef WEIGHTED_BALANCE_HPP
#define WEIGHTED_BALANCE_HPP

#ifdef MADVORO_WITH_MPI

#include <vector>
#include <algorithm>
#include <optional>
#include <mpi.h>
#include "mpi/serialize/mpi_commands.hpp"
// #include "mpi/mpi_commands.hpp"

#ifndef EPSILON
    #define EPSILON 1e-12
#endif // EPSILON

namespace MadVoro
{
    namespace MPI
    {
        namespace Balance
        {
            template<typename T>
            struct LocationSpecifier : public MadVoro::MPI::Serializable
            {
                T value;
                double weightBefore;
                double weightEqual;
                double weightAfter;
                size_t elementsBefore;
                size_t elementsEqual;
                size_t elementsAfter;

                LocationSpecifier() = default;

                force_inline size_t dump(Serializer *serializer) const
                {
                    size_t bytes = 0;
                    bytes += serializer->insert(this->value);
                    bytes += serializer->insert(this->weightBefore);
                    bytes += serializer->insert(this->weightEqual);
                    bytes += serializer->insert(this->weightAfter);
                    bytes += serializer->insert(this->elementsBefore);
                    bytes += serializer->insert(this->elementsEqual);
                    bytes += serializer->insert(this->elementsAfter);
                    return bytes;
                }

                force_inline size_t load(const Serializer *serializer, std::size_t byteOffset)
                {
                    size_t bytes = 0;
                    bytes += serializer->extract(this->value, byteOffset + bytes);
                    bytes += serializer->extract(this->weightBefore, byteOffset + bytes);
                    bytes += serializer->extract(this->weightEqual, byteOffset + bytes);
                    bytes += serializer->extract(this->weightAfter, byteOffset + bytes);
                    bytes += serializer->extract(this->elementsBefore, byteOffset + bytes);
                    bytes += serializer->extract(this->elementsEqual, byteOffset + bytes);
                    bytes += serializer->extract(this->elementsAfter, byteOffset + bytes);
                    return bytes;
                }

                void syncronizeAll(const MPI_Comm &comm)
                {
                    MPI_Allreduce(MPI_IN_PLACE, &this->weightBefore, 1, MPI_DOUBLE, MPI_SUM, comm);
                    MPI_Allreduce(MPI_IN_PLACE, &this->weightEqual, 1, MPI_DOUBLE, MPI_SUM, comm);
                    MPI_Allreduce(MPI_IN_PLACE, &this->weightAfter, 1, MPI_DOUBLE, MPI_SUM, comm);
                    MPI_Allreduce(MPI_IN_PLACE, &this->elementsBefore, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
                    MPI_Allreduce(MPI_IN_PLACE, &this->elementsEqual, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
                    MPI_Allreduce(MPI_IN_PLACE, &this->elementsAfter, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
                }
            };


            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            struct OrderStatisticsJob
            {
                using ValuesIterator = typename std::vector<T>::const_iterator;
                using WeightsIterator = typename std::vector<double>::const_iterator;
                using OrderStatisticsIterator = typename std::vector<size_t>::const_iterator;

                struct SubJob
                {
                    ValuesIterator valuesVecBegin;
                    ValuesIterator valuesVecEnd;
                    WeightsIterator weightsVecBegin;
                    WeightsIterator weightsVecEnd;
                    OrderStatisticsIterator orderStatisticsBegin;
                    OrderStatisticsIterator orderStatisticsEnd;
                    size_t elementsBefore;
                    size_t elementsAfter;
                    double weightBefore; // weight before the scope of `weightsVecBegin`
                    double weightAfter; // weight after the scope of `weightsVecEnd`
                };

                OrderStatisticsJob(const std::vector<T> &values, const std::vector<double> &weights, const Comparator &comp = std::less<T>{}, const MPI_Comm &comm = MPI_COMM_WORLD):
                    values(values), weights(weights), comp(comp), comm(comm)
                {
                    MPI_Comm_rank(this->comm, &this->rank);
                    MPI_Comm_size(this->comm, &this->size);
                    this->accumulatingWeights.resize(this->values.size());
                    double currentWeight = 0;
                    for(size_t i = 0; i < values.size(); i++)
                    {
                        currentWeight += this->weights[i];
                        this->accumulatingWeights[i] = currentWeight;
                    }
                    this->allLocalWeight = currentWeight;
                    MPI_Allreduce(&this->allLocalWeight, &this->allWeight, 1, MPI_DOUBLE, MPI_SUM, this->comm);
                    this->allElementsNum = this->values.size();
                    MPI_Allreduce(MPI_IN_PLACE, &this->allElementsNum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, this->comm);
                }

                /**
                 * Returns the median weight element of the current watched vector, the number of values smaller or equal to it, and bigger than it
                 */
                std::optional<LocationSpecifier<T>> getCurrentMedian() const
                {
                    // TODO: make more efficient!
                    size_t beginIndex = std::distance(this->values.begin(), this->subjob.valuesVecBegin);
                    size_t endIndex = std::distance(this->values.begin(), this->subjob.valuesVecEnd);
                    if(beginIndex == endIndex)
                    {
                        return std::nullopt;
                    }
                    LocationSpecifier<T> location;
                    location.value = this->values[(beginIndex + endIndex) / 2];

                    location.weightBefore = this->subjob.weightBefore;
                    location.weightAfter = this->subjob.weightAfter;
                    location.weightEqual = 0;
                    location.elementsBefore = this->subjob.elementsBefore;
                    location.elementsAfter = this->subjob.elementsAfter;
                    location.elementsEqual = 0;
                    for(size_t i = beginIndex; i < endIndex; i++)
                    {
                        if(this->values[i] == location.value)
                        {
                            location.weightEqual += this->weights[i];
                            location.elementsEqual++;
                        }
                        else if(this->values[i] < location.value)
                        {
                            location.weightBefore += this->weights[i];
                            location.elementsBefore++;
                        }
                        else
                        {
                            location.weightAfter += this->weights[i];
                            location.elementsAfter++;
                        }
                    }
                    return location;
                }

                /**
                 * Finds the first stat order (between `this->orderStatisticsBegin` and `this->orderStatisticsEnd`) that is larger from the stat order
                 */
                inline OrderStatisticsIterator findFirstLargerStatOrder(size_t statOrder) const noexcept
                {
                    return std::lower_bound(this->subjob.orderStatisticsBegin, this->subjob.orderStatisticsEnd, statOrder);
                }

                inline double getWeightOnBoundaries(size_t index1, size_t index2) const noexcept
                {
                    assert(index1 <= this->values.size());
                    assert(index2 <= this->values.size());
                    assert(index1 <= index2);
                    if(index1 == index2)
                    {
                        return 0;
                    }
                    double w1 = (index1 == 0)? 0 : this->accumulatingWeights[index1 - 1];
                    double w2 = (index2 == 0)? 0 : this->accumulatingWeights[index2 - 1];
                    return w2 - w1;
                }

                LocationSpecifier<T> checkLocalSmallerEqualBiggerElements(const T &value) const;

                const std::vector<T> &values;
                const std::vector<double> &weights;
                const Comparator &comp;
                std::vector<double> accumulatingWeights;
                double allLocalWeight;
                double allWeight;
                size_t allElementsNum;
                rank_t rank, size;
                MPI_Comm comm;
                std::vector<T> result;
                SubJob subjob; // changes recursively
            };

            template<typename T, typename Comparator>
            LocationSpecifier<T> OrderStatisticsJob<T, Comparator>::checkLocalSmallerEqualBiggerElements(const T &value) const
            {
                LocationSpecifier<T> result;
                result.value = value;
                result.elementsBefore = this->subjob.elementsBefore;
                result.elementsAfter = this->subjob.elementsAfter;
                result.weightBefore = this->subjob.weightBefore;
                result.weightAfter = this->subjob.weightAfter;
                result.elementsEqual = 0;
                auto firstEqual = std::lower_bound(this->subjob.valuesVecBegin, this->subjob.valuesVecEnd, value);
                auto firstAboveEqual = firstEqual;
                while((firstAboveEqual < this->subjob.valuesVecEnd) && (*firstAboveEqual) == value)
                {
                    firstAboveEqual++;
                }
                size_t below = std::distance(this->subjob.valuesVecBegin, firstEqual);
                result.elementsBefore += below;
                size_t equal = std::distance(firstEqual, firstAboveEqual);
                result.elementsEqual = equal;
                size_t above = std::distance(firstAboveEqual, this->subjob.valuesVecEnd);
                result.elementsAfter += above;

                assert(result.elementsBefore <= this->values.size());
                assert(result.elementsEqual <= this->values.size());
                assert(result.elementsAfter <= this->values.size());
                assert((result.elementsBefore + result.elementsEqual + result.elementsAfter) <= this->values.size());

                size_t belowIndex = std::distance(this->values.cbegin(), this->subjob.valuesVecBegin);
                size_t equalIndex1 = std::distance(this->values.cbegin(), firstEqual);
                size_t equalIndex2 = std::distance(this->values.cbegin(), firstAboveEqual);
                size_t aboveIndex = std::distance(this->values.cbegin(), this->subjob.valuesVecEnd);
                result.weightBefore += this->getWeightOnBoundaries(belowIndex, equalIndex1);
                result.weightEqual = this->getWeightOnBoundaries(equalIndex1, equalIndex2);
                result.weightAfter += this->getWeightOnBoundaries(equalIndex2, aboveIndex);
                // assert((result.weightBefore + result.weightEqual + result.weightAfter) <= (this->allLocalWeight * (1 + EPSILON)));

                return result;
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            LocationSpecifier<T> getWeightedMedianOfMedians(const OrderStatisticsJob<T> &job, const std::optional<LocationSpecifier<T>> &myMedian)
            {
                LocationSpecifier<T> medianOfMedians;
                std::vector<LocationSpecifier<T>> myMedianToBcast;
                if(myMedian.has_value())
                {
                    myMedianToBcast.push_back(myMedian.value());
                }
                std::vector<LocationSpecifier<T>> medianOfMediansToBcast = MPI_All_cast(myMedianToBcast, job.comm);
                if(medianOfMediansToBcast.empty())
                {
                    MadVoro::Exception::MadVoroException eo("getWeightedMedianOfMedians: median of medians is empty");
                    eo.addEntry("Statistics to find: ", std::vector<size_t>(job.subjob.orderStatisticsBegin, job.subjob.orderStatisticsEnd));
                    eo.addEntry("Number of all elements", job.allElementsNum);
                    throw eo;
                }
                std::sort(medianOfMediansToBcast.begin(), medianOfMediansToBcast.end(), 
                        [&job](const LocationSpecifier<T> &a, const LocationSpecifier<T> &b) -> bool
                        {
                            return job.comp(a.value, b.value);
                        });
                T medianOfMediansValue = medianOfMediansToBcast[medianOfMediansToBcast.size() / 2].value;
                medianOfMedians = job.checkLocalSmallerEqualBiggerElements(medianOfMediansValue);
                medianOfMedians.syncronizeAll(job.comm);
                return medianOfMedians;
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            void recursivelyGetWeightedStatOrder(OrderStatisticsJob<T, Comparator> &job)
            {
                assert(std::distance(job.subjob.valuesVecBegin, job.subjob.valuesVecEnd) == std::distance(job.subjob.weightsVecBegin, job.subjob.weightsVecEnd));
                assert(job.subjob.elementsBefore <= job.values.size());
                assert(job.subjob.elementsAfter <= job.values.size());
                assert(job.subjob.weightBefore <= job.allLocalWeight * (1 + EPSILON));
                assert(job.subjob.weightAfter <= job.allLocalWeight * (1 + EPSILON));
                
                // step 1 - if the list of statistics is empty, return
                if(job.subjob.orderStatisticsBegin == job.subjob.orderStatisticsEnd)
                {
                    return; // no stats to find
                }
                
                // step 2 - get local weighted median of currnet list. If not exist, return
                std::optional<LocationSpecifier<T>> localWeightedMedian = job.getCurrentMedian();
                
                // step 3 (collective) - get the local weighted median of medians
                LocationSpecifier<T> medianOfMedians = getWeightedMedianOfMedians(job, localWeightedMedian);
                
                // step 4 - get the stat order of the element
                // the element might occupity multiple statisticle orders, if there are any duplications
                // so, we calculate the minimal stat order it occupity, and the maximal one
                size_t medianOfMediansMinStatOrder = static_cast<size_t>(std::floor(job.allElementsNum * medianOfMedians.weightBefore / job.allWeight));
                size_t medianOfMediansMaxStatOrder = static_cast<size_t>(std::ceil(job.allElementsNum * (medianOfMedians.weightBefore + medianOfMedians.weightEqual) / job.allWeight));

                // find how many times the median is in the list of the statistics. Cut the list of statistics accordingly
                auto statOrderLastBelow = std::lower_bound(job.subjob.orderStatisticsBegin, job.subjob.orderStatisticsEnd, medianOfMediansMinStatOrder);
                auto statOrderFirstAbove = std::upper_bound(job.subjob.orderStatisticsBegin, job.subjob.orderStatisticsEnd, medianOfMediansMaxStatOrder);
                size_t timesAppearing = std::distance(statOrderLastBelow, statOrderFirstAbove);
                for(size_t times = 0; times < timesAppearing; times++)
                {
                    job.result.push_back(medianOfMedians.value); // found!
                }

                // step 5 - find the elements smaller than the median of medians
                LocationSpecifier<T> medianLocation = job.checkLocalSmallerEqualBiggerElements(medianOfMedians.value);
                
                // steps 6, 7 - generate two sub-jobs, left and right
                typename OrderStatisticsJob<T, Comparator>::SubJob left;
                left.valuesVecBegin = job.subjob.valuesVecBegin;
                left.valuesVecEnd = std::next(job.values.cbegin(), medianLocation.elementsBefore  /* + medianLocation.elementsEqual */);
                left.weightsVecBegin = job.subjob.weightsVecBegin;
                left.weightsVecEnd = std::next(job.weights.cbegin(), medianLocation.elementsBefore  /* + medianLocation.elementsEqual */);
                left.orderStatisticsBegin = job.subjob.orderStatisticsBegin;
                left.orderStatisticsEnd = statOrderLastBelow;
                left.elementsBefore = job.subjob.elementsBefore;
                assert(left.elementsBefore <= job.values.size());
                left.elementsAfter = medianLocation.elementsEqual + medianLocation.elementsAfter; // medianLocation.elementsEqual + medianLocation.elementsAfter
                assert(left.elementsAfter <= job.values.size());
                left.weightBefore = job.subjob.weightBefore;
                left.weightAfter = medianLocation.weightEqual + medianLocation.weightAfter; // medianLocation.weightEqual + medianLocation.weightAfter;

                typename OrderStatisticsJob<T, Comparator>::SubJob right;
                right.valuesVecBegin = std::next(left.valuesVecEnd, medianLocation.elementsEqual);
                right.valuesVecEnd = job.subjob.valuesVecEnd;
                right.weightsVecBegin = std::next(left.weightsVecEnd, medianLocation.elementsEqual);
                right.weightsVecEnd = job.subjob.weightsVecEnd;
                right.orderStatisticsBegin = statOrderFirstAbove;
                right.orderStatisticsEnd = job.subjob.orderStatisticsEnd;
                right.elementsBefore = medianLocation.elementsBefore + medianLocation.elementsEqual; // medianLocation.elementsBefore
                assert(right.elementsBefore <= job.values.size());
                right.elementsAfter = job.subjob.elementsAfter;
                assert(right.elementsAfter <= job.values.size());
                right.weightBefore = medianLocation.weightBefore + medianLocation.weightEqual; // medianLocation.weightBefore
                right.weightAfter = job.subjob.weightAfter;

                // call recursively
                job.subjob = left;
                recursivelyGetWeightedStatOrder(job);
                job.subjob = right;
                recursivelyGetWeightedStatOrder(job);
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            std::vector<T> getWeightedOrderStatistics(const std::vector<T> &values, const std::vector<double> &weights, const std::vector<size_t> &orderStatistics, const Comparator &comp = std::less<T>{}, const MPI_Comm &comm = MPI_COMM_WORLD)
            {
                if(values.size() != weights.size())
                {
                    MadVoro::Exception::MadVoroException eo("Values and weights must have the same size");
                    eo.addEntry("Values size", values.size());
                    eo.addEntry("Weights size", weights.size());
                    throw eo;
                }
                
                size_t N = values.size();
                std::vector<std::pair<T, double>> valuesWeightsSorted;
                for(size_t i = 0; i < N; i++)
                {
                    valuesWeightsSorted.emplace_back(values[i], weights[i]);
                }
                std::sort(valuesWeightsSorted.begin(), valuesWeightsSorted.end(), [&comp](const std::pair<T, double> &lhs, const std::pair<T, double> &rhs){return comp(lhs.first, rhs.first);});
                std::vector<T> valuesSorted;
                std::for_each(valuesWeightsSorted.cbegin(), valuesWeightsSorted.cend(), [&valuesSorted](const std::pair<T, double> &value){valuesSorted.emplace_back(value.first);});
                std::vector<double> weightsSorted;
                std::for_each(valuesWeightsSorted.cbegin(), valuesWeightsSorted.cend(), [&weightsSorted](const std::pair<T, double> &value){weightsSorted.emplace_back(value.second);});
                std::vector<size_t> orderStatisticsSorted = orderStatistics;
                std::sort(orderStatisticsSorted.begin(), orderStatisticsSorted.end());
                OrderStatisticsJob<T, Comparator> job(valuesSorted, weightsSorted, comp, comm);
                job.subjob.valuesVecBegin = valuesSorted.cbegin();
                job.subjob.valuesVecEnd = valuesSorted.cend();
                job.subjob.weightsVecBegin = weightsSorted.cbegin();
                job.subjob.weightsVecEnd = weightsSorted.cend();
                job.subjob.orderStatisticsBegin = orderStatisticsSorted.cbegin();
                job.subjob.orderStatisticsEnd = orderStatisticsSorted.cend();
                job.subjob.elementsBefore = 0;
                job.subjob.elementsAfter = 0;
                job.subjob.weightBefore = 0;
                job.subjob.weightAfter = 0;
                recursivelyGetWeightedStatOrder(job);
                std::sort(job.result.begin(), job.result.end());
                MPI_Barrier(job.comm);
                return job.result;
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            std::vector<T> getWeightedBorders(const std::vector<T> &values, const std::vector<double> &weights, const Comparator &comp = std::less<T>{}, const MPI_Comm &comm = MPI_COMM_WORLD)
            {
                rank_t size;
                MPI_Comm_size(comm, &size);

                size_t totalSize = values.size();
                MPI_Allreduce(MPI_IN_PLACE, &totalSize, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

                size_t factor = totalSize / size; 
                std::vector<size_t> orderStatistics;
                for(int _rank = 0; _rank < size - 1; _rank++)
                {
                    orderStatistics.push_back(factor * (_rank + 1));
                }
                orderStatistics.push_back(totalSize);
                return getWeightedOrderStatistics(values, weights, orderStatistics, comp, comm);
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            std::vector<T> getOrderStatistics(const std::vector<T> &values, const std::vector<size_t> &orderStatistics, const Comparator &comp = std::less<T>{}, const MPI_Comm &comm = MPI_COMM_WORLD)
            {
                return getOrderStatistics(values, std::vector<double>(values.size(), 1.0), orderStatistics, comp, comm);
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            std::vector<T> getBorders(const std::vector<T> &values, const Comparator &comp = std::less<T>{}, const MPI_Comm &comm = MPI_COMM_WORLD)
            {
                return getWeightedBorders(values, std::vector<double>(values.size(), 1.0), comp, comm);
            }
        }
    }
}

#endif // MADVORO_WITH_MPI

#endif // WEIGHTED_BALANCE_HPP