#ifndef _MADVORO_EXCHANGE_HPP
#define _MADVORO_EXCHANGE_HPP
#ifdef MADVORO_WITH_MPI

#include <mpi.h>
#include <functional>

template<typename T>
struct ExchangeAnswer
{
    std::vector<T> output;
    std::vector<size_t> indicesToMe;
    std::vector<int> processesSend;
    std::vector<std::vector<size_t>> indicesToProcesses;
    std::vector<int> processesRecv;
    std::vector<std::vector<T>> answerByProcesses;
};

namespace MadVoro
{
    namespace MPI
    {
        namespace Exchange
        {
            #define EXCHANGE_DATA_SEND_TAG 2605

            template<typename T>
            void initializeReceive(ExchangeAnswer<T> &answer, const std::vector<size_t> &sizes, std::vector<MPI_Request> &requests, const MPI_Comm &comm)
            {
                int rank, size;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &size);

                for(int _rank = 0; _rank < size; _rank++)
                {
                    if(_rank == rank)
                    {
                        continue;
                    }
                    if(sizes[_rank] == 0)
                    {
                        continue;
                    }

                    // ensure symmetry:
                    if(std::find(answer.processesSend.begin(), answer.processesSend.end(), _rank) == answer.processesSend.end())
                    {
                        answer.processesSend.push_back(_rank);
                        answer.indicesToProcesses.emplace_back(std::vector<size_t>());
                    }

                    answer.processesRecv.push_back(_rank);
                    answer.answerByProcesses.emplace_back(std::vector<T>());
                    std::vector<T> &answerVec = answer.answerByProcesses.back();
                    answerVec.resize(sizes[_rank]);
                    requests.push_back(MPI_REQUEST_NULL);
                    MPI_Irecv(&answerVec[0], sizeof(T) * sizes[_rank], MPI_BYTE, _rank, EXCHANGE_DATA_SEND_TAG, comm, &requests[requests.size() - 1]);
                }
            }

            template<typename T>
            void makeSend(const std::vector<T> &data, ExchangeAnswer<T> &answer, std::vector<MPI_Request> &requests, const MPI_Comm &comm, std::vector<std::vector<T>> &sendVectors)
            {
                int rank, size;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &size);

                sendVectors.reserve(answer.processesSend.size());

                // for each processor to send, send the data by non-
                for(size_t i = 0; i < answer.processesSend.size(); i++)
                {
                    int _rank = answer.processesSend[i];
                    if(_rank == rank)
                    {
                        // irrelevant, and shouldn't happen at all
                        continue;
                    }
                    if(answer.indicesToProcesses[i].size() == 0)
                    {
                        // nothing to send, just move on
                        continue;
                    }
                    sendVectors.push_back(std::vector<T>());
                    std::vector<T> &dataToSend = sendVectors[sendVectors.size() - 1];
                    dataToSend.reserve(answer.indicesToProcesses[i].size());
                    for(const size_t &dataIdx : answer.indicesToProcesses[i])
                    {
                        dataToSend.push_back(data[dataIdx]);
                    }
                    requests.push_back(MPI_REQUEST_NULL);
                    MPI_Isend(&dataToSend[0], sizeof(T) * dataToSend.size(), MPI_BYTE, _rank, EXCHANGE_DATA_SEND_TAG, comm, &requests[requests.size() - 1]);
                }
            }

            template<typename T, typename ExchangeDetermineFunc = std::function<int(const T&)>>
            ExchangeAnswer<T> dataExchange(const std::vector<T> &data, const ExchangeDetermineFunc &getOwner, const MPI_Comm &comm = MPI_COMM_WORLD)
            {
                int rank, size;
                MPI_Comm_rank(comm, &rank);
                MPI_Comm_size(comm, &size);

                std::vector<MPI_Request> requests;
                ExchangeAnswer<T> answer{};

                std::vector<size_t> sizesToSend(size, 0), sizesToMe(size, 0);
                size_t N = data.size();
                for(size_t i = 0; i < N; i++)
                {
                    int _rank = getOwner(data[i]);
                    if(_rank != rank)
                    {
                        if(sizesToSend[_rank] == 0)
                        {
                            // first time
                            answer.processesSend.push_back(_rank);
                            answer.indicesToProcesses.emplace_back(std::vector<size_t>());
                        }
                        // data is not mine
                        size_t index = std::distance(answer.processesSend.begin(), std::find(answer.processesSend.begin(), answer.processesSend.end(), _rank));
                        answer.indicesToProcesses[index].push_back(i);
                        sizesToSend[_rank]++;
                    }
                    else
                    {
                        answer.indicesToMe.push_back(i);
                        answer.output.push_back(data[i]);
                    }
                }
                MPI_Alltoall(&sizesToSend[0], sizeof(size_t), MPI_BYTE, &sizesToMe[0], sizeof(size_t), MPI_BYTE, comm);
                
                requests.reserve(2 * size);
                std::vector<std::vector<T>> sendVectors;

                initializeReceive(answer, sizesToMe, requests, comm);
                makeSend(data, answer, requests, comm, sendVectors);

                if(!requests.empty())
                {
                    MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
                }
                
                // there's a requirement to have all the answers ordered by the order of the sentproc vector
                size_t answerSize = answer.processesSend.size();
                for(size_t i = 0; i < answerSize; i++)
                {
                    int _rank = answer.processesSend[i];
                    if(_rank == rank)
                    {
                        continue; // should not happen, but just to be sure
                    }
                    size_t index = std::distance(answer.processesRecv.begin(), std::find(answer.processesRecv.begin(), answer.processesRecv.end(), _rank));
                    if(index == answer.processesRecv.size())
                    {
                        continue; // no answer
                    }
                    const std::vector<T> &_anserVec = answer.answerByProcesses[index];
                    answer.output.insert(answer.output.end(), _anserVec.cbegin(), _anserVec.cend());
                }
                return answer;
            }
        }
    }
}

#endif // MADVORO_WITH_MPI
#endif // _MADVORO_EXCHANGE_HPP