#ifndef WAIT_UNTIL_ANSWERED_QUERY_AGENT_HPP
#define WAIT_UNTIL_ANSWERED_QUERY_AGENT_HPP

#ifdef MADVORO_WITH_MPI

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <mpi.h>

// set data structure:
#include <boost/container/flat_set.hpp>
#include <unordered_set>

#ifdef TIMING
    #include <chrono>
#endif // TIMING

#include "QueryAgent.hpp"

#define TAG_REQUEST 200
#define TAG_RESPONSE 201
#define TAG_FINISHED 202

#define UNDEFINED_BUFFER_IDX -1
#define FLUSH_QUERIES_NUM 50

namespace MadVoro
{
    namespace QueryAgent
    {

        template<typename QueryData, typename AnswerType>
        class WaitUntilAnsweredQueryAgent : public QueryAgent<QueryData, AnswerType>
        {
        public:
            template<typename T>
            using _set = boost::container::flat_set<T>;

            using _subQueryData = SubQueryData<QueryData>;
            using _queryBatchInfo = QueryBatchInfo<QueryData, AnswerType>;
            using _queryInfo = QueryInfo<QueryData, AnswerType>;

            WaitUntilAnsweredQueryAgent(const TalkAgent<QueryData> *talkAgent, AnswerAgent<QueryData, AnswerType> *answerAgent, bool sendToSelf = false, const MPI_Comm &comm = MPI_COMM_WORLD);

            virtual ~WaitUntilAnsweredQueryAgent() = default;
            
            QueryBatchInfo<QueryData, AnswerType> runBatch(const std::vector<QueryData> &queries) override;
            inline std::vector<std::vector<size_t>> &getRecvData() override{return this->recvData;};
            inline std::vector<int> &getRecvProc() override{return this->recvProcessorsRanks;};

        private:
            std::vector<MPI_Request> requests;
            std::vector<std::vector<char>> buffers; // send buffers, so that they will not be allocated on the stack
            size_t receivedUntilNow; // number of answers I received until now
            size_t shouldReceiveInTotal; // number of answers I have to receive (to know when to finish)
            bool finishedMyQueries; // if finished to answer my queries
            bool answeredAll;
            size_t currentlyAnswered;
            size_t hasToAnswerAtAll;
            std::vector<size_t> numOfSentQueries;
            MPI_Request numOfSentRequest;
            std::vector<int> ones;
            std::vector<char> arriveBuffer;

            std::vector<int> recvProcessorsRanks;  // a vector of ranks, that we received data from
            std::vector<std::vector<size_t>> recvData; // a vector of vectors. The vector in index `i` contains the data indices (relatively to my final answer of the batch) that sent to `this->sentProcessorsRanks[i]`.
            std::vector<size_t> ranksBufferIdx; // maps each rank to an index `i`, where this->buffers[i] contains the prepared buffer for sending to the rank (requests are sent in chunks).

            void receiveQueries(QueryBatchInfo<QueryData, AnswerType> &batch);
            void answerQueries();
            void sendQuery(const QueryInfo<QueryData, AnswerType> &query);
            void rearrangeResult(_queryBatchInfo &queriesBatch);

            void syncAmounts();
            int checkForFinishMessages() const;
            void flushBuffer(int _rank);
            inline void flushAll(){for(int i = 0; i < this->size; i++) this->flushBuffer(i);};
        };

        template<typename QueryData, typename AnswerType>
        WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::WaitUntilAnsweredQueryAgent(const TalkAgent<QueryData> *talkAgent, AnswerAgent<QueryData, AnswerType> *answerAgent, bool sendToSelf, const MPI_Comm &comm):
                QueryAgent<QueryData, AnswerType>(talkAgent, answerAgent, sendToSelf, comm)
        {
            this->ranksBufferIdx.resize(this->size, UNDEFINED_BUFFER_IDX);
            this->ones.resize(this->size, 1);
        }

        template<typename QueryData, typename AnswerType>
        void WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::receiveQueries(_queryBatchInfo &batch)
        {
            if(this->receivedUntilNow >= this->shouldReceiveInTotal)
            {
                return;
            }
            MPI_Status status;
            int receivedAnswer = 0;

            std::vector<_queryInfo> &queries = batch.queriesAnswers;
            MPI_Iprobe(MPI_ANY_SOURCE, TAG_RESPONSE, this->comm, &receivedAnswer, &status);

            std::vector<char> buffer;

            while(receivedAnswer)
            {
                // received a message
                ++this->receivedUntilNow;

                // prepare the reading buffer for receiving
                int count;
                MPI_Get_count(&status, MPI_BYTE, &count);
                if(buffer.size() < static_cast<size_t>(count))
                {
                    buffer.resize(count);
                }

                // receive
                MPI_Recv(&buffer[0], count, MPI_BYTE, status.MPI_SOURCE, TAG_RESPONSE, this->comm, MPI_STATUS_IGNORE);

                // decode the message - id first, then length, then the data itself
                long int id;
                long int length;

                id = *reinterpret_cast<long int*>(buffer.data()); // decode id
                length = *reinterpret_cast<long int*>(buffer.data() + sizeof(long int)); // decode length

                if(id < 0 or static_cast<size_t>(id) >= queries.size())
                {
                    MadVoro::Exception::MadVoroException eo("WaitUntilAnsweredQueryAgent::receiveQueries, id of answered query is illegal");
                    eo.addEntry("Id", id);
                    eo.addEntry("Expected range", "0-" + std::to_string(queries.size()));
                    eo.addEntry("Received Query", *reinterpret_cast<const QueryData*>(buffer.data()));
                    throw eo;
                }
                if(length > 0)
                {
                    // insert the results to the data received by rank `status.MPI_SOURCE` and to the queries result
                    queries[id].finalResults.resize(queries[id].finalResults.size() + length);

                    // base pointers to data
                    AnswerType* base = reinterpret_cast<AnswerType*>(buffer.data() + 2 * sizeof(long int));
                    for(size_t i = 0; i < static_cast<size_t>(length); i++)
                    {
                        queries[id].finalResults[i] = base[i];
                        batch.dataByRanks[status.MPI_SOURCE].emplace_back(base[i]);
                    }
                }
                else
                {
                    if(length < 0)
                    {
                        MadVoro::Exception::MadVoroException eo("In WaitUntilAnsweredQueryAgent::receiveQueries, length of query is negative");
                        eo.addEntry("Length", length);
                        throw eo;
                    }
                }
                MPI_Iprobe(MPI_ANY_SOURCE, TAG_RESPONSE, this->comm, &receivedAnswer, &status);
            }

            if(this->finishedMyQueries and this->shouldReceiveInTotal == this->receivedUntilNow)
            {
                this->syncAmounts();
            }
        }

        template<typename QueryData, typename AnswerType>
        void WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::answerQueries()
        {
            static int totalArrived = 0;

            MPI_Status status;
            int arrivedNew = 0;

            MPI_Iprobe(MPI_ANY_SOURCE, TAG_REQUEST, this->comm, &arrivedNew, &status);
            
            // while arrived new messages, and we should answer until the end, or answer until a bound we haven't reached to
            while(arrivedNew != 0)
            {
                int count;
                MPI_Get_count(&status, MPI_BYTE, &count);
                if(this->arriveBuffer.size() < static_cast<size_t>(count))
                {
                    this->arriveBuffer.resize(count);
                }
                MPI_Recv(this->arriveBuffer.data(), count, MPI_BYTE, status.MPI_SOURCE, TAG_REQUEST, this->comm, MPI_STATUS_IGNORE);
                int subQueries = count / sizeof(_subQueryData);
                totalArrived += subQueries;
                for(int i = 0; i < subQueries; i++)
                {
                    const _subQueryData &query = *reinterpret_cast<_subQueryData*>(&this->arriveBuffer[i * sizeof(_subQueryData)]);
                    // calculate the result
                    std::vector<AnswerType> result = this->answerAgent->answer(query.data, status.MPI_SOURCE);
                    long int resultSize = static_cast<long int>(result.size());

                    this->buffers.push_back(std::vector<char>());
                    std::vector<char> &to_send = this->buffers.back();
                    size_t msg_size = 2 * sizeof(long int) + resultSize * sizeof(AnswerType);
                    to_send.resize(msg_size);

                    long int id = query.parent_id;

                *reinterpret_cast<long int*>(to_send.data()) = id;
                *reinterpret_cast<long int*>(to_send.data() + sizeof(long int)) = resultSize;

                    if(resultSize > 0)
                    {
                        AnswerType *toSendData = reinterpret_cast<AnswerType*>(to_send.data() + sizeof(id) + sizeof(resultSize));
                        std::memcpy(toSendData, result.data(), resultSize * sizeof(AnswerType));
                    }
                    this->requests.push_back(MPI_REQUEST_NULL);
                    MPI_Isend(&to_send[0], msg_size, MPI_BYTE, status.MPI_SOURCE, TAG_RESPONSE, this->comm, &this->requests.back());
                    this->currentlyAnswered++;
                }
                MPI_Iprobe(MPI_ANY_SOURCE, TAG_REQUEST, this->comm, &arrivedNew, &status);
            }
        }

        template<typename QueryData, typename AnswerType>
        void WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::sendQuery(const _queryInfo &query)
        {
            typename TalkAgent<QueryData>::RanksSet talkingRanks = this->talkAgent->getTalkList(query.data);
            for(const int &_rank : talkingRanks)
            {
                if((_rank == this->rank) and (!this->sendToSelf))
                {
                    continue; // unnecessary to send
                }
                int bufferIdx = this->ranksBufferIdx[_rank];
                // check if a flush is needed

                if((bufferIdx != UNDEFINED_BUFFER_IDX) and ((this->buffers[bufferIdx].size() / sizeof(_subQueryData)) >= FLUSH_QUERIES_NUM))
                {
                    // send buffer
                    this->flushBuffer(_rank);
                }
                bufferIdx = this->ranksBufferIdx[_rank];
                if(bufferIdx == UNDEFINED_BUFFER_IDX)
                {
                    this->ranksBufferIdx[_rank] = this->buffers.size();
                    this->buffers.emplace_back(std::vector<char>());
                    this->buffers.back().reserve(sizeof(_subQueryData) * FLUSH_QUERIES_NUM);
                }
                bufferIdx = this->ranksBufferIdx[_rank];
                this->buffers[bufferIdx].resize(this->buffers[bufferIdx].size() + sizeof(_subQueryData));
                _subQueryData &subQuery = *reinterpret_cast<_subQueryData*>(&(*(this->buffers[bufferIdx].end() - sizeof(_subQueryData))));
                subQuery.data = query.data;
                subQuery.parent_id = query.id;
                this->numOfSentQueries[_rank]++;
                ++this->shouldReceiveInTotal;
            }
        }

        template<typename QueryData, typename AnswerType>
        void WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::syncAmounts()
        {
            MPI_Ireduce_scatter(this->numOfSentQueries.data(), &this->hasToAnswerAtAll, this->ones.data(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, this->comm, &this->numOfSentRequest);
        }

        template<typename QueryData, typename AnswerType>
        int WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::checkForFinishMessages() const
        {
            int arrived = 0;
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, TAG_FINISHED, this->comm, &arrived, &status);
            if(arrived)
            {
                int dummy = 0;
                MPI_Recv(&dummy, 1, MPI_BYTE, MPI_ANY_SOURCE, TAG_FINISHED, this->comm, MPI_STATUS_IGNORE);
                return 1;
            }
            return 0;
        }

        template<typename QueryData, typename AnswerType>
        void WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::flushBuffer(int _rank)
        {
            if((_rank < 0) or (_rank >= this->size))
            {
                MadVoro::Exception::MadVoroException eo("WaitUntilAnsweredQueryAgent::flushBuffer: invalid rank buffer flush");
                eo.addEntry("Rank", _rank);
                throw eo;
            }
            int bufferIdx = this->ranksBufferIdx[_rank];
            if(bufferIdx == UNDEFINED_BUFFER_IDX)
            {
                return;
            }

            std::vector<char> &buffer = this->buffers[bufferIdx];
            if(buffer.size() > 0)
            {
                this->requests.push_back(MPI_REQUEST_NULL);
                MPI_Isend(buffer.data(), buffer.size(), MPI_BYTE, _rank, TAG_REQUEST, this->comm, &this->requests.back());
            }
            this->ranksBufferIdx[_rank] = UNDEFINED_BUFFER_IDX;
        }

        template<typename QueryData, typename AnswerType>
        void WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::rearrangeResult(_queryBatchInfo &queriesBatch)
        {
            // make receive, by the order of the recv array
            for(size_t i = 0; i < this->recvProcessorsRanks.size(); i++)
            {
                int _rank = this->recvProcessorsRanks[i];
                std::vector<size_t> &rankRecvData = this->recvData[i]; // indices vector
                const std::vector<AnswerType> &newDataFromRank = queriesBatch.dataByRanks[_rank]; // the data itself
                rankRecvData.reserve(newDataFromRank.size());

                queriesBatch.result.reserve(queriesBatch.result.size() + rankRecvData.size());
                for(const AnswerType &_data : newDataFromRank)
                {
                    rankRecvData.push_back(queriesBatch.result.size()); // the index of the data received by this rank
                    queriesBatch.result.emplace_back(_data);
                }
            }
        }

        template<typename QueryData, typename AnswerType>
        QueryBatchInfo<QueryData, AnswerType> WaitUntilAnsweredQueryAgent<QueryData, AnswerType>::runBatch(const std::vector<QueryData> &queries)
        {
            this->receivedUntilNow = 0; // reset the receive counter
            this->shouldReceiveInTotal = 0; // reset the should-be-received counter
            this->answeredAll = false;
            this->numOfSentQueries = std::vector<size_t>(this->size, 0);
            this->numOfSentRequest = MPI_REQUEST_NULL;

            for(std::vector<size_t> &_receivedDataFromRank : this->recvData)
            {
                _receivedDataFromRank.clear();
            }

            this->buffers.clear();
            size_t originalQueriesNum = queries.size();
            this->buffers.reserve(10 * originalQueriesNum); // heuristic
            this->requests.reserve(10 * originalQueriesNum); // heuristic
            this->requests.clear();
            _queryBatchInfo queriesBatch;
            #ifdef TIMING
                queriesBatch.beginClockTime std::chrono::system_clock::now();
            #endif // TIMING

            queriesBatch.queriesAnswers.reserve(originalQueriesNum);
            std::vector<_queryInfo> &queriesInfo = queriesBatch.queriesAnswers;
            queriesBatch.dataByRanks.resize(this->size);
            this->finishedMyQueries = queries.empty();
            this->currentlyAnswered = 0;
            this->hasToAnswerAtAll = std::numeric_limits<size_t>::max();
            bool updatedHasToAnswerAtAll = false;
            size_t i = 0;

            // if doesn't have any queries, send a finish message
            if(this->finishedMyQueries)
            {
                this->syncAmounts();
                #ifdef TIMING
                    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
                    queriesBatch.finishSubmittingTime = std::chrono::duration_cast<std::chrono::duration<double>>(now - queriesBatch.beginClockTime).count();
                #endif // TIMING
            }
            while((!this->finishedMyQueries) or (this->currentlyAnswered < this->hasToAnswerAtAll) or (this->receivedUntilNow < this->shouldReceiveInTotal))
            {
                if(!this->finishedMyQueries)
                {
                    const QueryData &queryData = queries[i];
                    queriesInfo.push_back({queryData, i, std::vector<AnswerType>()});
                    _queryInfo &query = queriesInfo.back();

                    this->sendQuery(query);
                    if(i == (originalQueriesNum - 1))
                    {
                        // send the rest of the waiting (buffered) requests
                        this->flushAll();
                        #ifdef TIMING
                            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
                            queriesBatch.finishSubmittingTime = std::chrono::duration_cast<std::chrono::duration<double>>(now - queriesBatch.beginClockTime).count();
                        #endif // TIMING

                        this->finishedMyQueries = true;
                        // if had several queries, but no communication was needed, send a finish message
                        if(this->shouldReceiveInTotal == 0)
                        {
                            this->syncAmounts();
                        }
                    }
                }
                
                if(this->receivedUntilNow < this->shouldReceiveInTotal)
                {
                    this->receiveQueries(queriesBatch);
                }
                if(not updatedHasToAnswerAtAll)
                {
                    int arrived;
                    MPI_Test(&this->numOfSentRequest, &arrived, MPI_STATUS_IGNORE);
                    if(arrived)
                    {
                        updatedHasToAnswerAtAll = true;
                    }
                }
                if(not this->answeredAll)
                {
                    this->answerQueries();
                    this->answeredAll = (this->currentlyAnswered >= this->hasToAnswerAtAll);
                }
                ++i;
            }

            if(not this->requests.empty())
            {
                MPI_Waitall(this->requests.size(), this->requests.data(), MPI_STATUSES_IGNORE); // make sure any query was indeed received
            }

            #ifdef TIMING
                std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
                queriesBatch.finishTime = std::chrono::duration_cast<std::chrono::duration<double>>(now - queriesBatch.beginClockTime).count();
            #endif // TIMING

            // add to the list the processors that sent us a message for the first time
            for(int _rank = 0; _rank < this->size; _rank++)
            {
                if(_rank == this->rank)
                {
                    continue; // shouldn't be relevant
                }
                if(queriesBatch.dataByRanks[_rank].empty())
                {
                    continue; // the rank `_rank` did not send us any message
                }
                size_t rankIndex = std::find(this->recvProcessorsRanks.begin(), this->recvProcessorsRanks.end(), _rank) - this->recvProcessorsRanks.begin();
                if(rankIndex == this->recvProcessorsRanks.size())
                {
                    // rank is not inside the recvProcessors rank, add it
                    this->recvProcessorsRanks.push_back(_rank);
                    this->recvData.emplace_back(std::vector<size_t>());
                }
            }

            this->rearrangeResult(queriesBatch);

            MPI_Barrier(this->comm);

            return queriesBatch;
        }
    }
}

#endif // MADVORO_WITH_MPI

#endif // WAIT_UNTIL_ANSWERED_QUERY_AGENT_HPP
