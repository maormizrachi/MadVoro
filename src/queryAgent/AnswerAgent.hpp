#ifndef _ANSWER_AGENT_HPP
#define _ANSWER_AGENT_HPP

#include <vector>

namespace MadVoro
{
    namespace QueryAgent
    {
        template<typename QueryData, typename AnswerType>
        class AnswerAgent
        {
        public:
            virtual ~AnswerAgent() = default;
            
            virtual std::vector<AnswerType> answer(const QueryData &query, int _rank) = 0;
        };
    }
}

#endif // _ANSWER_AGENT_HPP