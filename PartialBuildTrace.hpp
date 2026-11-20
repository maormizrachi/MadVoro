#ifndef MADVORO_PARTIAL_BUILD_TRACE_HPP
#define MADVORO_PARTIAL_BUILD_TRACE_HPP

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#ifdef MADVORO_WITH_MPI
  #include <mpi_utils/BatchTrace.hpp>
#endif // MADVORO_WITH_MPI

namespace MadVoro {

// Per-rank diagnostics of Voronoi3D::BuildPartiallyParallel, recorded while the level set by
// Voronoi3D::SetPartialBuildTraceLevel is above 0. Every value is local to this rank: MadVoro neither reduces nor
// prints them. Names and their order are fixed, the same on every rank, attempt and round, so a caller can reduce the
// values by position. Names ending in "_s" are wall seconds; the others are counts, flags or digests, all integers
// stored exactly in a double (below 2^53).
struct PartialBuildTraceRound {
    std::vector<std::pair<std::string, double>> Values;
};

struct PartialBuildTraceAttempt {
    std::vector<std::pair<std::string, double>> Prepare;  // per-call values: preparation, trees, ghost setup, Voronoi
    std::vector<PartialBuildTraceRound> Rounds;           // ghost rounds, at most PartialBuildTraceMaxRounds
};

// Reset at the start of every BuildPartiallyParallel call; holds one attempt (that call) when tracing, none otherwise.
struct PartialBuildTrace {
    std::vector<PartialBuildTraceAttempt> Attempts;
};

// Rounds beyond this count are added into the last entry, whose rounds_aggregated then exceeds 1.
constexpr std::size_t PartialBuildTraceMaxRounds = 32;

namespace PartialBuildTraceDetail {

// X(Identifier, "name"). Every field is additive, so aggregated rounds are plain sums.
#define MADVORO_PBT_PREPARE_FIELDS(X) \
    X(PrepReduceWaitS, "prep_reduce_wait_s") \
    X(PrepPayloadPackS, "prep_payload_pack_s") \
    X(PrepPointsUpdateS, "prep_points_update_s") \
    X(PrepRoutingOctAgent, "prep_routing_oct_agent") \
    X(PrepRoutingRebuilt, "prep_routing_rebuilt") \
    X(PrepRoutingReused, "prep_routing_reused") \
    X(PrepResultUnpackS, "prep_result_unpack_s") \
    X(PrepBuildInitializeS, "prep_build_initialize_s") \
    X(PrepInputPoints, "prep_input_points") \
    X(PrepOwnedPoints, "prep_owned_points") \
    X(PrepTargetPoints, "prep_target_points") \
    X(PrepAllowExchange, "prep_allow_exchange") \
    X(PrepAllowRebalance, "prep_allow_rebalance") \
    X(TreeAllOwnedS, "tree_all_owned_s") \
    X(TreeTargetS, "tree_target_s") \
    X(SeededRadii, "seeded_radii") \
    X(GhostSetupS, "ghost_setup_s") \
    X(GhostEntryVoteWaitS, "ghost_entry_vote_wait_s") \
    X(GhostInitialExchangeS, "ghost_initial_exchange_s") \
    X(GhostInitialPoints, "ghost_initial_points") \
    X(GhostInitialInsertS, "ghost_initial_insert_s") \
    X(GhostLoopS, "ghost_loop_s") \
    X(GhostRounds, "ghost_rounds") \
    X(GhostSmallBatchesSkipped, "ghost_small_batches_skipped") \
    X(GhostBigBatchesSkipped, "ghost_big_batches_skipped") \
    X(GhostFinalPeerMapS, "ghost_final_peer_map_s") \
    X(GhostPeers, "ghost_peers") \
    X(GhostSendOnlyPeers, "ghost_send_only_peers") \
    X(GhostRecvOnlyPeers, "ghost_recv_only_peers") \
    X(GhostPointsReceived, "ghost_points_received") \
    X(GhostPointsSent, "ghost_points_sent") \
    X(PeerOrderDigest, "peer_order_digest") \
    X(GhostIdDigest, "ghost_id_digest") \
    X(PeerListDigest, "peer_list_digest") \
    X(GhostSendIdDigest, "ghost_send_id_digest") \
    X(GhostRecvIdDigest, "ghost_recv_id_digest") \
    X(VorBuildS, "vor_build_s") \
    X(VorCmS, "vor_cm_s") \
    X(VorCentroidSyncS, "vor_centroid_sync_s") \
    X(VorCentroidSyncBarrierS, "vor_centroid_sync_barrier_s") \
    X(VorVolumeFillS, "vor_volume_fill_s") \
    X(VorVolumeExchangeS, "vor_volume_exchange_s") \
    X(VorVolumeExchangeBarrierS, "vor_volume_exchange_barrier_s") \
    X(VorVolumeScatterS, "vor_volume_scatter_s") \
    X(StagePrepareS, "stage_prepare_s") \
    X(StageDelaunayS, "stage_delaunay_s") \
    X(StageTreeS, "stage_tree_s") \
    X(StageRadiiS, "stage_radii_s") \
    X(StageRangeFinderS, "stage_range_finder_s") \
    X(StageGhostsS, "stage_ghosts_s") \
    X(StageVoronoiS, "stage_voronoi_s") \
    X(TraceOverheadS, "trace_overhead_s") \
    X(ClosestRanksMode, "closest_ranks_mode") \
    X(RoutingLeaves, "routing_leaves") \
    X(GhostSetDigest, "ghost_set_digest")

#define MADVORO_PBT_ROUND_FIELDS(X) \
    X(RoundsAggregated, "rounds_aggregated") \
    X(SmallPoints, "small_points") \
    X(LargePoints, "large_points") \
    X(SmallQueries, "small_queries") \
    X(BigQueries, "big_queries") \
    X(SmallBatchRun, "small_batch_run") \
    X(BigBatchRun, "big_batch_run") \
    X(SelfPointsInserted, "self_points_inserted") \
    X(RemoteSmallPointsInserted, "remote_small_points_inserted") \
    X(RemoteBigPointsInserted, "remote_big_points_inserted") \
    X(MirrorPointsInserted, "mirror_points_inserted") \
    X(NetNewTetrahedra, "net_new_tetrahedra") \
    X(Overflows, "overflows") \
    X(OverflowsSelfOnly, "overflows_self_only") \
    X(GrowthSteps, "growth_steps") \
    X(FinishedPoints, "finished_points") \
    X(BigCacheMisses, "big_cache_misses") \
    X(RoundS, "round_s") \
    X(CreateBatchesS, "create_batches_s") \
    X(VoteOverlapS, "vote_overlap_s") \
    X(SelfSmallSearchS, "self_small_search_s") \
    X(SelfBigSearchS, "self_big_search_s") \
    X(SelfInsertS, "self_insert_s") \
    X(VoteWaitS, "vote_wait_s") \
    X(MirrorsS, "mirrors_s") \
    X(PointTetrasS, "point_tetras_s") \
    X(DetermineNextS, "determine_next_s") \
    X(BigOwnerEvals, "big_owner_evals") \
    X(BigLeafEvals, "big_leaf_evals") \
    X(BigDistanceS, "big_distance_s") \
    X(BigIndexBuildS, "big_index_build_s")

// One remote batch; in a round record these appear twice, prefixed "small_" and then "big_".
#define MADVORO_PBT_BATCH_FIELDS(X) \
    X(BatchS, "batch_s") \
    X(PeerMapS, "peer_map_s") \
    X(InsertS, "insert_s") \
    X(ConstructS, "construct_s") \
    X(CtorBarrierS, "ctor_barrier_s") \
    X(SubmitS, "submit_s") \
    X(DrainS, "drain_s") \
    X(TeardownS, "teardown_s") \
    X(FinalBarrierS, "final_barrier_s") \
    X(RearrangeS, "rearrange_s") \
    X(ServeS, "serve_s") \
    X(Fanout, "fanout") \
    X(QueriesAnswered, "queries_answered") \
    X(AnswersReceived, "answers_received") \
    X(SubmitPolls, "submit_polls") \
    X(SubmitPollsEmpty, "submit_polls_empty") \
    X(DrainPolls, "drain_polls") \
    X(DrainPollsEmpty, "drain_polls_empty") \
    X(SubmitSampledEmpty, "submit_sampled_empty") \
    X(SubmitSampledEmptyS, "submit_sampled_empty_s") \
    X(SubmitSampledProductive, "submit_sampled_productive") \
    X(SubmitSampledProductiveS, "submit_sampled_productive_s") \
    X(SubmitSampledQueryS, "submit_sampled_query_s") \
    X(DrainSampledEmpty, "drain_sampled_empty") \
    X(DrainSampledEmptyS, "drain_sampled_empty_s") \
    X(DrainSampledProductive, "drain_sampled_productive") \
    X(DrainSampledProductiveS, "drain_sampled_productive_s") \
    X(QueryMsgsSent, "query_msgs_sent") \
    X(QueryBytesSent, "query_bytes_sent") \
    X(QueryMsgsRecv, "query_msgs_recv") \
    X(QueryBytesRecv, "query_bytes_recv") \
    X(QueryDispatchSize, "query_dispatch_size") \
    X(QueryDispatchTimer, "query_dispatch_timer") \
    X(QueryDispatchFlush, "query_dispatch_flush") \
    X(AnswerMsgsSent, "answer_msgs_sent") \
    X(AnswerBytesSent, "answer_bytes_sent") \
    X(AnswerMsgsRecv, "answer_msgs_recv") \
    X(AnswerBytesRecv, "answer_bytes_recv") \
    X(AnswerDispatchSize, "answer_dispatch_size") \
    X(AnswerDispatchTimer, "answer_dispatch_timer") \
    X(AnswerDispatchFlush, "answer_dispatch_flush")

#define MADVORO_PBT_ENUM_ENTRY(Id, Name) Id,
#define MADVORO_PBT_NAME_ENTRY(Id, Name) Name,

namespace PrepareField {
enum Id : std::size_t { MADVORO_PBT_PREPARE_FIELDS(MADVORO_PBT_ENUM_ENTRY) Count };
}

namespace BatchField {
enum Id : std::size_t { MADVORO_PBT_BATCH_FIELDS(MADVORO_PBT_ENUM_ENTRY) Count };
}

namespace RoundField {
enum Id : std::size_t { MADVORO_PBT_ROUND_FIELDS(MADVORO_PBT_ENUM_ENTRY) GeneralCount };
constexpr std::size_t SmallBatchBase = GeneralCount;
constexpr std::size_t BigBatchBase = GeneralCount + BatchField::Count;
constexpr std::size_t Count = GeneralCount + 2 * BatchField::Count;
}

inline const char *const PrepareFieldNames[PrepareField::Count] = {MADVORO_PBT_PREPARE_FIELDS(MADVORO_PBT_NAME_ENTRY)};
inline const char *const RoundFieldNames[RoundField::GeneralCount] = {MADVORO_PBT_ROUND_FIELDS(MADVORO_PBT_NAME_ENTRY)};
inline const char *const BatchFieldNames[BatchField::Count] = {MADVORO_PBT_BATCH_FIELDS(MADVORO_PBT_NAME_ENTRY)};

#undef MADVORO_PBT_ENUM_ENTRY
#undef MADVORO_PBT_NAME_ENTRY
#undef MADVORO_PBT_PREPARE_FIELDS
#undef MADVORO_PBT_ROUND_FIELDS
#undef MADVORO_PBT_BATCH_FIELDS

using PrepareRecord = std::array<double, PrepareField::Count>;
using RoundRecord = std::array<double, RoundField::Count>;

// The clock of Voronoi3D's existing stage timers, so their time points can be reused instead of read again.
using Clock = std::chrono::high_resolution_clock;

inline double Seconds(Clock::duration Elapsed) {
    return std::chrono::duration<double>(Elapsed).count();
}

// Seconds since Mark; Mark moves to now.
inline double Lap(Clock::time_point &Mark) {
    const Clock::time_point Now = Clock::now();
    const double Elapsed = Seconds(Now - Mark);
    Mark = Now;
    return Elapsed;
}

inline std::string RoundFieldName(std::size_t Index) {
    if(Index < RoundField::GeneralCount) {
        return RoundFieldNames[Index];
    }
    if(Index < RoundField::BigBatchBase) {
        return std::string("small_") + BatchFieldNames[Index - RoundField::SmallBatchBase];
    }
    return std::string("big_") + BatchFieldNames[Index - RoundField::BigBatchBase];
}

inline std::vector<std::pair<std::string, double>> NamedPrepare(const PrepareRecord &Record) {
    std::vector<std::pair<std::string, double>> Values;
    Values.reserve(Record.size());
    for(std::size_t Index = 0; Index < Record.size(); ++Index) {
        Values.emplace_back(PrepareFieldNames[Index], Record[Index]);
    }
    return Values;
}

inline std::vector<std::pair<std::string, double>> NamedRound(const RoundRecord &Record) {
    std::vector<std::pair<std::string, double>> Values;
    Values.reserve(Record.size());
    for(std::size_t Index = 0; Index < Record.size(); ++Index) {
        Values.emplace_back(RoundFieldName(Index), Record[Index]);
    }
    return Values;
}

// Adds a round to the record list, folding rounds beyond PartialBuildTraceMaxRounds into the last entry.
inline void AppendRound(std::vector<RoundRecord> &Rounds, const RoundRecord &Round) {
    if(Rounds.size() < PartialBuildTraceMaxRounds) {
        Rounds.push_back(Round);
        return;
    }
    RoundRecord &Last = Rounds.back();
    for(std::size_t Index = 0; Index < Last.size(); ++Index) {
        Last[Index] += Round[Index];
    }
}

// Order-sensitive 64-bit digest (splitmix64 finalizer per value), reported in its low 53 bits so a double holds it
// exactly. The same input sequence gives the same digest on every rank and run.
inline std::uint64_t DigestMix(std::uint64_t Value) {
    Value += 0x9E3779B97F4A7C15ULL;
    Value = (Value ^ (Value >> 30)) * 0xBF58476D1CE4E5B9ULL;
    Value = (Value ^ (Value >> 27)) * 0x94D049BB133111EBULL;
    return Value ^ (Value >> 31);
}

inline void DigestAdd(std::uint64_t &Digest, std::uint64_t Value) {
    Digest = DigestMix(Digest ^ DigestMix(Value));
}

inline double DigestValue(std::uint64_t Digest) {
    return static_cast<double>(Digest & ((std::uint64_t(1) << 53) - 1));
}

// Order-independent digest of a set of (peer, value) pairs: the sum of a mix of each pair, so the same set gives the
// same digest whatever order its pairs arrive in.
inline void SetDigestAdd(std::uint64_t &Digest, std::uint64_t Peer, std::uint64_t Value) {
    Digest += DigestMix(DigestMix(Peer) ^ Value);
}

#ifdef MADVORO_WITH_MPI
// Adds one batch record of the query agent into a round record at Base (RoundField::SmallBatchBase or BigBatchBase).
inline void AddBatchTrace(RoundRecord &Round, std::size_t Base, const QueryBatchTrace &Batch) {
    using namespace BatchField;
    Round[Base + ConstructS] += Batch.ConstructSeconds;
    Round[Base + CtorBarrierS] += Batch.ConstructorBarrierSeconds;
    Round[Base + SubmitS] += Batch.SubmitSeconds;
    Round[Base + DrainS] += Batch.DrainSeconds;
    Round[Base + TeardownS] += Batch.TeardownSeconds;
    Round[Base + FinalBarrierS] += Batch.FinalBarrierSeconds;
    Round[Base + RearrangeS] += Batch.RearrangeSeconds;
    Round[Base + ServeS] += Batch.ServeSeconds;
    Round[Base + Fanout] += static_cast<double>(Batch.SubQueriesSent);
    Round[Base + QueriesAnswered] += static_cast<double>(Batch.QueriesAnswered);
    Round[Base + AnswersReceived] += static_cast<double>(Batch.AnswersReceived);
    Round[Base + SubmitPolls] += static_cast<double>(Batch.SubmitPolls);
    Round[Base + SubmitPollsEmpty] += static_cast<double>(Batch.SubmitPollsEmpty);
    Round[Base + DrainPolls] += static_cast<double>(Batch.DrainPolls);
    Round[Base + DrainPollsEmpty] += static_cast<double>(Batch.DrainPollsEmpty);
    Round[Base + SubmitSampledEmpty] += static_cast<double>(Batch.SubmitSampled[0]);
    Round[Base + SubmitSampledEmptyS] += Batch.SubmitSampledSeconds[0];
    Round[Base + SubmitSampledProductive] += static_cast<double>(Batch.SubmitSampled[1]);
    Round[Base + SubmitSampledProductiveS] += Batch.SubmitSampledSeconds[1];
    Round[Base + SubmitSampledQueryS] += Batch.SubmitSampledQuerySeconds;
    Round[Base + DrainSampledEmpty] += static_cast<double>(Batch.DrainSampled[0]);
    Round[Base + DrainSampledEmptyS] += Batch.DrainSampledSeconds[0];
    Round[Base + DrainSampledProductive] += static_cast<double>(Batch.DrainSampled[1]);
    Round[Base + DrainSampledProductiveS] += Batch.DrainSampledSeconds[1];
    const BuffersManagerTrace &Queries = Batch.QueryManager;
    Round[Base + QueryMsgsSent] += static_cast<double>(Queries.MessagesSent);
    Round[Base + QueryBytesSent] += static_cast<double>(Queries.BytesSent);
    Round[Base + QueryMsgsRecv] += static_cast<double>(Queries.MessagesReceived);
    Round[Base + QueryBytesRecv] += static_cast<double>(Queries.BytesReceived);
    Round[Base + QueryDispatchSize] += static_cast<double>(Queries.DispatchBySize);
    Round[Base + QueryDispatchTimer] += static_cast<double>(Queries.DispatchByTimer);
    Round[Base + QueryDispatchFlush] += static_cast<double>(Queries.DispatchByFlush);
    const BuffersManagerTrace &Answers = Batch.AnswerManager;
    Round[Base + AnswerMsgsSent] += static_cast<double>(Answers.MessagesSent);
    Round[Base + AnswerBytesSent] += static_cast<double>(Answers.BytesSent);
    Round[Base + AnswerMsgsRecv] += static_cast<double>(Answers.MessagesReceived);
    Round[Base + AnswerBytesRecv] += static_cast<double>(Answers.BytesReceived);
    Round[Base + AnswerDispatchSize] += static_cast<double>(Answers.DispatchBySize);
    Round[Base + AnswerDispatchTimer] += static_cast<double>(Answers.DispatchByTimer);
    Round[Base + AnswerDispatchFlush] += static_cast<double>(Answers.DispatchByFlush);
}
#endif // MADVORO_WITH_MPI

} // namespace PartialBuildTraceDetail

} // namespace MadVoro

#endif // MADVORO_PARTIAL_BUILD_TRACE_HPP
