#ifndef QUERY_AGENT_HPP
#define QUERY_AGENT_HPP

#include "AnswerAgent.hpp"
#include "TalkAgent.hpp"
#include "QueryData.hpp"

namespace MadVoro
{
    namespace QueryAgent
    {

        template<typename QueryData, typename AnswerType>
        class QueryAgent
        {
        public:
            QueryAgent(const TalkAgent<QueryData> *talkAgent, AnswerAgent<QueryData, AnswerType> *answerAgent, bool sendToSelf = false, const MPI_Comm &comm = MPI_COMM_WORLD);
            
            virtual ~QueryAgent() = default;

            virtual QueryBatchInfo<QueryData, AnswerType> runBatch(const std::vector<QueryData> &queries) = 0;
            
            virtual std::vector<std::vector<size_t>> &getRecvData() = 0;
            
            virtual std::vector<int> &getRecvProc() = 0;

        protected:
            MPI_Comm comm;
            int rank, size;
            int mpiInitialized;
            const TalkAgent<QueryData> *talkAgent; // answers to whom we should talk to, given a query
            AnswerAgent<QueryData, AnswerType> *answerAgent; // an answer agent
            bool sendToSelf; // should send queries to self
        };

        template<typename QueryData, typename AnswerType>
        QueryAgent<QueryData, AnswerType>::QueryAgent(const TalkAgent<QueryData> *talkAgent, AnswerAgent<QueryData, AnswerType> *answerAgent, bool sendToSelf, const MPI_Comm &comm)
        {
            MPI_Initialized(&this->mpiInitialized);
            if(this->mpiInitialized)
            {
                MPI_Comm_size(comm, &this->size);
                MPI_Comm_rank(comm, &this->rank);
            }
            this->talkAgent = talkAgent;
            this->answerAgent = answerAgent;
            this->sendToSelf = sendToSelf;
            this->comm = comm;
        }
    }
}

#endif // QUERY_AGENT_HPP