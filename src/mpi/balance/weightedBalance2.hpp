#ifndef WEIGHTED_BALANCE2_HPP
#define WEIGHTED_BALANCE2_HPP

#ifdef MADVORO_WITH_MPI

#include <vector>
#include <functional>
#include <numeric>
#include <iostream>
#include <mpi.h>
#include "mpi/serialize/mpi_commands.hpp"

namespace MadVoro
{
    namespace MPI
    {
        namespace Balance
        {
            template<typename T>
            struct WeightedPair : public MadVoro::MPI::Serializable
            {
                T value;
                double weight;

                inline WeightedPair(const T &value_ = T(), const double &weight_ = 0.0): value(value_), weight(weight_){};

                force_inline size_t dump(Serializer *serializer) const override
                {
                    size_t bytes = 0;
                    bytes += serializer->insert(this->value);
                    bytes += serializer->insert(this->weight);
                    return bytes;
                }

                force_inline size_t load(const Serializer *serializer, std::size_t byteOffset) override
                {
                    size_t bytes = 0;
                    bytes += serializer->extract(this->value, byteOffset + bytes);
                    bytes += serializer->extract(this->weight, byteOffset + bytes);
                    return bytes;
                }
            };

            inline bool isBalanced(double myWeight, double tolerance, const MPI_Comm &comm)
            {
                rank_t rank, size;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &size);

                double totalWeight;
                MPI_Allreduce(&myWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, comm);
                double idealWeight = totalWeight / static_cast<double>(size);

                struct
                {
                    double weight;
                    int rank;
                } maxWeight, myWeight2 = {myWeight, rank};
                MPI_Allreduce(&myWeight2, &maxWeight, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
                if(maxWeight.weight <= (1 + tolerance) * idealWeight)
                {
                    return true;
                }
                return false;
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            std::vector<T> bordersIteration(std::vector<WeightedPair<T>> &valuesAndWeights, const Comparator &comp, const MPI_Comm &comm)
            {
                rank_t rank, size;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &size);

                std::sort(valuesAndWeights.begin(), valuesAndWeights.end(), [&comp](const WeightedPair<T> &lhs, const WeightedPair<T> &rhs){return comp(lhs.value, rhs.value);});
                
                std::vector<double> prefixWeights(valuesAndWeights.size());
                double currentWeight = 0;
                for(size_t i = 0; i < valuesAndWeights.size(); i++)
                {
                    currentWeight += valuesAndWeights[i].weight;
                    prefixWeights[i] = currentWeight;
                }
                double totalWeight = currentWeight;

                std::vector<WeightedPair<T>> localBorders;
                if(valuesAndWeights.empty())
                {
                    localBorders.resize(size - 1, WeightedPair<T>(std::numeric_limits<T>::max(), 0.0));
                }
                else
                {
                    for(rank_t _rank = 1; _rank < size; _rank++)
                    {
                        double target = _rank * totalWeight / static_cast<double>(size);
                        auto it = std::lower_bound(prefixWeights.cbegin(), prefixWeights.cend(), target, comp);
                        size_t index = std::min<size_t>(std::distance(prefixWeights.cbegin(), it), valuesAndWeights.size() - 1);
                        localBorders.emplace_back(valuesAndWeights[index].value, _rank * totalWeight / static_cast<double>(size));
                    }
                }

                std::vector<WeightedPair<T>> allBorders(size * localBorders.size());
                MPI_Allgather(localBorders.data(), sizeof(WeightedPair<T>) * localBorders.size(), MPI_BYTE, allBorders.data(), sizeof(WeightedPair<T>) * localBorders.size(), MPI_BYTE, comm);

                std::sort(allBorders.begin(), allBorders.end(), [&comp](const WeightedPair<T> &lhs, const WeightedPair<T> &rhs){return comp(lhs.value, rhs.value);});
                std::vector<double> bordersSum;
                double totalBordersSum = 0.0;
                for(const WeightedPair<T> &pair : allBorders)
                {
                    totalBordersSum += pair.weight;
                    bordersSum.push_back(totalBordersSum);
                }

                std::vector<T> borders;
                for(rank_t _rank = 1; _rank < size; _rank++)
                {
                    double target = _rank * totalBordersSum / static_cast<double>(size);
                    auto it = std::lower_bound(bordersSum.cbegin(), bordersSum.cend(), target, comp);
                    size_t index = std::min<size_t>(std::distance(bordersSum.cbegin(), it), allBorders.size() - 1);
                    borders.push_back(allBorders[index].value);
                }

                std::function<rank_t(const WeightedPair<T>&)> ownership = [&borders, &comp](const WeightedPair<T> &pair) -> rank_t
                {
                    size_t idx = std::distance(borders.cbegin(), std::upper_bound(borders.cbegin(), borders.cend(), pair.value, comp));
                    return static_cast<rank_t>(idx);
                };

                valuesAndWeights = std::move(MadVoro::MPI::MPI_Exchange_by_ownership(valuesAndWeights, ownership, comm));
                return borders;
            }

            template<typename T, typename Comparator = std::function<bool(const T&, const T&)>>
            std::vector<T> getWeightedBorders2(const std::vector<T> &values, const std::vector<double> &weights, const Comparator &comp = std::less<T>{}, const MPI_Comm &comm = MPI_COMM_WORLD)
            {
                std::vector<WeightedPair<T>> valueWeightPairs;
                valueWeightPairs.reserve(values.size());
                for(size_t i = 0; i < values.size(); i++)
                {
                    valueWeightPairs.emplace_back(values[i], weights[i]);
                }
                std::sort(valueWeightPairs.begin(), valueWeightPairs.end(), [&comp](const WeightedPair<T> &lhs, const WeightedPair<T> &rhs){return comp(lhs.value, rhs.value);});

                size_t iterations = 0;
                std::vector<T> borders;
                double myWeight;
                do
                {
                    borders = bordersIteration(valueWeightPairs, comp, comm);
                    myWeight = std::accumulate(valueWeightPairs.cbegin(), valueWeightPairs.cend(), 0.0, [](double acc, const WeightedPair<T> &pair){return acc + pair.weight;});
                    iterations++;
                } while(not isBalanced(myWeight, 0.05, comm) and iterations < 20);

                rank_t rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                if(rank == 0)
                {
                    std::cout << "Weighted borders determined in " << iterations << " iterations." << std::endl;
                }
                if(!borders.empty())
                {
                    borders.push_back(borders.back() + 1);
                }
                return borders;
            }
        }
    }
}

#endif // MADVORO_WITH_MPI

#endif // WEIGHTED_BALANCE2_HPP
