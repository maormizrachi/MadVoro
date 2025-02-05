#ifndef _TALK_AGENT_HPP
#define _TALK_AGENT_HPP

#ifdef MADVORO_WITH_MPI

namespace MadVoro
{
    namespace QueryAgent
    {
        template<typename QueryData>
        class TalkAgent
        {
        public:
            using RanksSet = boost::container::flat_set<int>;

            virtual ~TalkAgent() = default;

            virtual RanksSet getTalkList(const QueryData &query) const = 0;
        };
    }
}

#endif // MADVORO_WITH_MPI

#endif // _TALK_AGENT_HPP