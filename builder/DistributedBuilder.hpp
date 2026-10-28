// Builds a rank's cells with local builders, bringing the ghosts they need round by round (from the rank's points, its mirrors and periodic images, and other ranks)
#ifndef MADVORO_BUILDER_DISTRIBUTED_BUILDER_HPP
#define MADVORO_BUILDER_DISTRIBUTED_BUILDER_HPP

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>
#ifdef MADVORO_WITH_MPI
    #include <mpi.h>
    #include <mpi_utils/mpi_alltoall.hpp>
#endif // MADVORO_WITH_MPI
#include <MeshDecomposer3D/hilbert/HilbertOrder3D.hpp>

#include "types.hpp"
#include "builder/BuildDomain.hpp"
#include "builder/DistributedBuild.hpp"
#include "builder/LocalBuilder.hpp"
#include "geometry/FaceStage.hpp"
#include "geometry/FlatVoronoi.hpp"
#include "geometry/HostGeometryExecutor.hpp"
#include "geometry/kernels/MirrorCM.hpp"
#include "ghosts/GhostBringer.hpp"
#include "range/finders/OctTreeFinder.hpp"
#include "subdomains/MergeSubdomains.hpp"
#include "subdomains/SubdomainSearch.hpp"
#include "utils/container_utils.hpp"
#include "exception/MadVoroException.hpp"

namespace MadVoro
{

template <typename PointT>
class DistributedBuilder
{
public:
    using IVec = IndexedVector<PointT>;
    // `onHost`: a builder of one of several thread subdomains (geometry on the host, in its thread)
    using LocalBuilderFactory = std::function<std::unique_ptr<LocalBuilder<PointT>>(bool onHost)>;

    #ifdef MADVORO_WITH_MPI
        // A build that searches other ranks: who they are, and what to send them first (the rank's points that were their ghosts in the last build, indices among the rank's points)
        struct ExchangeRequest
        {
            std::shared_ptr<EnvironmentAgent<PointT>> envAgent;
            std::vector<int> carryProcs;
            std::vector<std::vector<std::size_t>> carryPoints;
            MPI_Comm comm;
        };
    #endif // MADVORO_WITH_MPI

    explicit DistributedBuilder(LocalBuilderFactory makeLocal) : makeLocal(std::move(makeLocal))
    {
    }

    // a copy keeps how to build, not the last build
    DistributedBuilder(const DistributedBuilder &other) : makeLocal(other.makeLocal)
    {
    }

    DistributedBuilder &operator=(const DistributedBuilder &other)
    {
        if(this != &other)
        {
            this->makeLocal = other.makeLocal;
            this->single = Part();
            this->parts.clear();
            this->rankTree.reset();
            this->builtTree.reset();
            #ifdef MADVORO_WITH_MPI
                this->exchange.reset();
            #endif // MADVORO_WITH_MPI
        }
        return *this;
    }

    //! \brief How many thread subdomains a build of `cellsNum` cells uses (1: one domain)
    static std::size_t SubdomainsFor(std::size_t requested, std::size_t cellsNum)
    {
        return std::max<std::size_t>(1, std::min({requested, cellsNum, SubdomainGrid<PointT>::MAX_SUBDOMAINS}));
    }

    /** \brief Builds the cells of rankPoints[cells[k]] (cell k), bringing their ghosts. radii: per rank point, in: the search radius to start from (if positive), out: the one to carry. knownGhosts (several subdomains only): other ranks' points, each given at once to the subdomain of its nearest point. `request` (or null): the search of other ranks, which also gives the build's ghost tables. */
    #ifdef MADVORO_WITH_MPI
        void Build(const BuildDomain<PointT> &domain, std::size_t subdomains, const std::vector<PointT> &rankPoints, const std::vector<std::size_t> &cells, std::vector<double> &radii, const std::vector<PointT> *knownGhosts, const ExchangeRequest *request, DistributedBuild<PointT> &out)
    #else // MADVORO_WITH_MPI
        void Build(const BuildDomain<PointT> &domain, std::size_t subdomains, const std::vector<PointT> &rankPoints, const std::vector<std::size_t> &cells, std::vector<double> &radii, const std::vector<PointT> *knownGhosts, DistributedBuild<PointT> &out)
    #endif // MADVORO_WITH_MPI
    {
        const auto start = std::chrono::high_resolution_clock::now();
        const std::size_t count = SubdomainsFor(subdomains, cells.size());
        this->domainLL = domain.ll;
        this->domainUR = domain.ur;
        out.receivedIndex.clear();
        out.rounds = 0;
        #ifdef MADVORO_WITH_MPI
            this->exchange.reset();
            // the last build's ghosts first: what the other ranks had of this rank's points
            std::vector<PointT> carried;
            std::vector<int> carriedSentProcs, carriedReceivedProcs;
            std::vector<std::vector<std::size_t>> carriedSentPoints, carriedReceivedPoints;
            if(request != nullptr)
            {
                // the ghost tables are this search's (a build that searches no other rank keeps the mesh's)
                out.duplicatedProcs.clear();
                out.duplicatedPoints.clear();
                out.ghostIndices.clear();
                this->CarryGhosts(*request, rankPoints, cells.size(), carried, carriedSentProcs, carriedSentPoints, carriedReceivedProcs, carriedReceivedPoints);
            }
            const std::vector<PointT> *carriedGhosts = (request != nullptr) ? &carried : nullptr;
        #else // MADVORO_WITH_MPI
            const std::vector<PointT> *carriedGhosts = nullptr;
        #endif // MADVORO_WITH_MPI
        if(count == 1 and knownGhosts != nullptr)
        {
            throw MadVoro::Exception::MadVoroException("DistributedBuilder::Build: known ghosts need several thread subdomains (BuildKnown builds one domain)");
        }
        // one domain inserts the rank's points in their order, thread subdomains along the Hilbert curve
        const std::vector<std::size_t> builtOrder = this->BuildTrees(domain, rankPoints, cells, count > 1);
        #ifdef MADVORO_WITH_MPI
            if(request != nullptr)
            {
                this->exchange.emplace(this->rankTree.get(), rankPoints, request->envAgent, SentPointsContainer(carriedSentProcs, carriedSentPoints), carriedReceivedProcs, request->comm);
                this->exchangeComm = request->comm;
            }
        #endif // MADVORO_WITH_MPI
        this->treesTime = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
        if(count == 1)
        {
            this->BuildOneDomain(domain, rankPoints, cells, radii, carriedGhosts, out);
        }
        else
        {
            this->BuildSubdomains(domain, count, rankPoints, cells, builtOrder, radii, knownGhosts, carriedGhosts, out);
        }
        #ifdef MADVORO_WITH_MPI
            if(request != nullptr)
            {
                this->GhostTables(carriedReceivedProcs, carriedReceivedPoints, out);
            }
        #endif // MADVORO_WITH_MPI
        int rank = 0;
        #ifdef MADVORO_WITH_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        #endif // MADVORO_WITH_MPI
        if(rank == 0)
        {
            std::cout << "Time for building (" << count << " thread subdomains): " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() << " seconds, " << out.rounds << " rounds" << std::endl;
        }
    }

    /** \brief Builds the cells of `points` as one domain, with every ghost known: the groups of `ghosts`, added in order (after 4 unused points). Builds the rank's trees and the start radii of rankPoints[cells[k]] = points[k]. */
    void BuildKnown(const BuildDomain<PointT> &domain, const std::vector<PointT> &rankPoints, const std::vector<std::size_t> &cells, std::vector<double> &radii, const std::vector<PointT> &points, const std::vector<const std::vector<PointT> *> &ghosts, DistributedBuild<PointT> &out)
    {
        #ifdef MADVORO_WITH_MPI
            this->exchange.reset();
        #endif // MADVORO_WITH_MPI
        this->domainLL = domain.ll;
        this->domainUR = domain.ur;
        out.receivedIndex.clear();
        out.rounds = 0;
        out.carriedGhostsIndex = points.size() + BOUNDING_POINTS;
        if(rankPoints.empty())
        {
            this->rankTree.reset();
            this->builtTree.reset();
        }
        else
        {
            this->BuildTrees(domain, rankPoints, cells, false);
        }
        Part &part = this->SinglePart();
        part.points = points;
        part.cellsNum = points.size();
        part.rankIndex.clear();
        for(std::size_t k = 0; k < cells.size(); ++k)
        {
            part.rankIndex.emplace_hint(part.rankIndex.end(), k, cells[k]);
        }
        if(not cells.empty())
        {
            this->InitialRadii(part, rankPoints, radii);
        }
        part.builder->Build(part.points, domain, HilbertOrder3D(part.points));
        std::unique_ptr<OctTree<IVec>> cellsTree;
        for(const std::vector<PointT> *group : ghosts)
        {
            this->AddKnownGhosts(part, *group, domain, cellsTree);
        }
        this->FinishPart(part, domain);
        this->OneDomainOutput(part, out);
    }

    //! \brief The tree of the rank's points (indices among them), and of the built points (cell indices) of the last build
    const std::shared_ptr<OctTree<IVec>> &RankTree(void) const
    {
        return this->rankTree;
    }

    const std::shared_ptr<OctTree<IVec>> &BuiltTree(void) const
    {
        return this->builtTree;
    }

    #ifdef MADVORO_WITH_MPI
        //! \brief The last build's search of other ranks (null without one)
        GhostExchange<PointT> *Exchange(void)
        {
            return this->exchange ? &*this->exchange : nullptr;
        }
    #endif // MADVORO_WITH_MPI

    //! \brief The builder of the last build as one domain, and its cells (null if it had thread subdomains)
    const LocalBuilder<PointT> *OneDomainBuilder(void) const
    {
        return this->parts.empty() ? this->single.builder.get() : nullptr;
    }

    const Geometry::FlatVoronoi *OneDomainVoronoi(void) const
    {
        return (this->parts.empty() and this->single.builder) ? &this->single.voronoi : nullptr;
    }

    void ShrinkToFit(void)
    {
        if(this->single.builder)
        {
            this->single.builder->ShrinkToFit();
        }
        for(Part &part : this->parts)
        {
            part.builder->ShrinkToFit();
        }
    }

    void ReleaseMemory(void)
    {
        this->single = Part();
        this->parts.clear();
    }

private:
    static constexpr std::size_t BOUNDING_POINTS = 4;
    static constexpr std::size_t NO_INDEX = std::numeric_limits<std::size_t>::max();
    using subdomain_t = typename SubdomainGrid<PointT>::subdomain_t;

    // A local builder, its points and its ghost search
    struct Part
    {
        std::unique_ptr<LocalBuilder<PointT>> builder;
        std::shared_ptr<RangeFinder<PointT>> finder;
        std::unique_ptr<GhostSearch<PointT>> search;
        std::vector<PointT> points; // the cells' points, then the ghosts
        std::size_t cellsNum = 0;
        AllPointsMap rankIndex; // a point -> its index among the rank's points (and the points other ranks sent)
        Geometry::FlatVoronoi voronoi;
        Geometry::FaceStage stage; // finishes the builder's raw faces into `voronoi`
        std::vector<PointT> CM; // the cells', and of a mirror ghost the mirror of its cell's
    };

    // A part as its ghost search sees it
    class PartMesh : public GhostSearchMesh<PointT>
    {
    public:
        // received (or null): the points other ranks sent, after the rank's points
        PartMesh(Part &part, const std::vector<PointT> &rankPoints, const RankGhostPoints<PointT> *received, std::vector<double> &radii) : part(part), rankPoints(rankPoints), received(received), radii(radii)
        {
        }

        const PointT &Point(std::size_t index) const override
        {
            return this->part.points[index];
        }

        std::size_t AllPointIndex(std::size_t cell) const override
        {
            return this->part.rankIndex.at(cell);
        }

        std::size_t RankPointsNum(void) const override
        {
            return this->rankPoints.size() + ((this->received == nullptr) ? 0 : this->received->Points().size());
        }

        void CellVertices(std::size_t cell, std::vector<PointT> &vertices) override
        {
            this->part.builder->CellVertices(cell, vertices);
        }

        double CellRadius(std::size_t cell) override
        {
            return this->part.builder->CellRadius(cell, this->part.points[cell]);
        }

        void FetchVertices(const std::vector<std::size_t> &cells) override
        {
            this->part.builder->FetchVertices(cells);
        }

        void SetCarriedRadius(std::size_t cell, double radius) override
        {
            this->radii[this->part.rankIndex.at(cell)] = radius;
        }

        void AddRankPoints(const std::vector<std::pair<std::size_t, PointT>> &points) override
        {
            std::vector<PointT> newPoints;
            newPoints.reserve(points.size());
            for(const auto &[pointIdxInAll, translation] : points)
            {
                this->part.rankIndex.emplace_hint(this->part.rankIndex.end(), this->part.points.size() + newPoints.size(), pointIdxInAll);
                const PointT &point = (pointIdxInAll < this->rankPoints.size()) ? this->rankPoints[pointIdxInAll] : this->received->Points()[pointIdxInAll - this->rankPoints.size()];
                newPoints.push_back(point + translation);
            }
            this->part.points.insert(this->part.points.end(), newPoints.begin(), newPoints.end());
            this->part.builder->Add(newPoints);
        }

        void AddMirrors(const std::vector<PointT> &points, const std::vector<std::pair<std::size_t, std::size_t>> &sources) override
        {
            std::vector<MirrorSource> mirrorSources(sources.size());
            for(std::size_t i = 0; i < sources.size(); ++i)
            {
                mirrorSources[i].cell = sources[i].second;
                mirrorSources[i].face = sources[i].first;
            }
            this->part.points.insert(this->part.points.end(), points.begin(), points.end());
            this->part.builder->AddMirrors(points, mirrorSources);
        }

    private:
        Part &part;
        const std::vector<PointT> &rankPoints;
        const RankGhostPoints<PointT> *received;
        std::vector<double> &radii;
    };

    Part &SinglePart(void)
    {
        this->parts.clear();
        if(not this->single.builder)
        {
            this->single.builder = this->makeLocal(false);
        }
        return this->single;
    }

    // the tree of the rank's points, and of the built points (a different one in a partial build); `hilbert`: inserted along the Hilbert curve, whose order of the built points it returns
    std::vector<std::size_t> BuildTrees(const BuildDomain<PointT> &domain, const std::vector<PointT> &rankPoints, const std::vector<std::size_t> &cells, bool hilbert)
    {
        const PointT width = domain.ur - domain.ll;
        const IVec lower(domain.ll - width * 0.00001, NO_INDEX), upper(domain.ur + width * 0.00001, NO_INDEX);
        bool builtInOrder = (cells.size() == rankPoints.size());
        for(std::size_t k = 0; k < cells.size() and builtInOrder; ++k)
        {
            builtInOrder = (cells[k] == k);
        }
        this->rankTree = std::make_shared<OctTree<IVec>>(lower, upper);
        const std::vector<std::size_t> order = hilbert ? HilbertOrder3D(rankPoints) : std::vector<std::size_t>();
        for(std::size_t i = 0; i < rankPoints.size(); ++i)
        {
            const std::size_t pointIdx = hilbert ? order[i] : i;
            const PointT &point = rankPoints[pointIdx];
            this->rankTree->insert(IVec(point.x, point.y, point.z, pointIdx));
        }
        if(builtInOrder)
        {
            this->builtTree = this->rankTree;
            return order;
        }
        this->builtTree = std::make_shared<OctTree<IVec>>(lower, upper);
        std::vector<PointT> builtPoints;
        builtPoints.reserve(cells.size());
        for(const std::size_t &pointIdx : cells)
        {
            builtPoints.push_back(rankPoints[pointIdx]);
        }
        const std::vector<std::size_t> builtOrder = hilbert ? HilbertOrder3D(builtPoints) : std::vector<std::size_t>();
        for(std::size_t i = 0; i < cells.size(); ++i)
        {
            const std::size_t k = hilbert ? builtOrder[i] : i;
            this->builtTree->insert(IVec(builtPoints[k].x, builtPoints[k].y, builtPoints[k].z, k));
        }
        return builtOrder;
    }

    #ifdef MADVORO_WITH_MPI
        // sends each rank of `request` its carried points (of the cells), and gives the points received: the ranks and indices sent and received (each received point's index among them)
        void CarryGhosts(const ExchangeRequest &request, const std::vector<PointT> &rankPoints, std::size_t cellsNum, std::vector<PointT> &carried, std::vector<int> &sentProcs, std::vector<std::vector<std::size_t>> &sentPoints, std::vector<int> &receivedProcs, std::vector<std::vector<std::size_t>> &receivedPoints) const
        {
            int size = 1;
            MPI_Comm_size(request.comm, &size);
            std::vector<std::vector<PointT>> toSend(size);
            std::vector<std::vector<std::size_t>> sentIndices(size);
            for(std::size_t k = 0; k < request.carryProcs.size(); ++k)
            {
                const int rank = request.carryProcs[k];
                for(const std::size_t &pointIdx : request.carryPoints[k])
                {
                    if(pointIdx < cellsNum)
                    {
                        sentIndices[rank].push_back(pointIdx);
                        toSend[rank].push_back(rankPoints[pointIdx]);
                    }
                }
            }
            int hasData = std::any_of(toSend.begin(), toSend.end(), [](const std::vector<PointT> &points) { return not points.empty(); }) ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &hasData, 1, MPI_INT, MPI_MAX, request.comm);
            std::vector<std::vector<std::size_t>> receivedIndices(size);
            carried.clear();
            if(hasData != 0)
            {
                const std::vector<std::vector<PointT>> received = MPI_Exchange_all_to_all_sparse(toSend, request.comm);
                for(int rank = 0; rank < size; ++rank)
                {
                    for(std::size_t i = 0; i < received[rank].size(); ++i)
                    {
                        receivedIndices[rank].push_back(carried.size() + i);
                    }
                    carried.insert(carried.end(), received[rank].begin(), received[rank].end());
                }
            }
            for(int rank = 0; rank < size; ++rank)
            {
                if(not sentIndices[rank].empty())
                {
                    sentProcs.push_back(rank);
                    sentPoints.emplace_back(std::move(sentIndices[rank]));
                }
                if(not receivedIndices[rank].empty())
                {
                    receivedProcs.push_back(rank);
                    receivedPoints.emplace_back(std::move(receivedIndices[rank]));
                }
            }
        }

        // the build's ghost tables: the ranks it shares points with (sent iff received), the rank's points each has, and where each one's points are among the ghosts
        void GhostTables(const std::vector<int> &carriedProcs, const std::vector<std::vector<std::size_t>> &carriedPoints, DistributedBuild<PointT> &out) const
        {
            const GhostExchange<PointT> &exchange = *this->exchange;
            const auto rankIndex = [&out](int rank)
            {
                const std::size_t index = std::find(out.duplicatedProcs.begin(), out.duplicatedProcs.end(), rank) - out.duplicatedProcs.begin();
                if(index == out.duplicatedProcs.size())
                {
                    out.duplicatedProcs.push_back(rank);
                    out.duplicatedPoints.emplace_back();
                    out.ghostIndices.emplace_back();
                }
                return index;
            };
            // the ghosts each rank sent: the carried ones, then those the search received
            for(std::size_t k = 0; k < carriedProcs.size(); ++k)
            {
                std::vector<std::size_t> &ghosts = out.ghostIndices[rankIndex(carriedProcs[k])];
                for(const std::size_t &i : carriedPoints[k])
                {
                    ghosts.push_back(out.carriedGhostsIndex + i);
                }
            }
            const std::vector<int> &receivedFrom = exchange.ReceivedFrom();
            const std::vector<std::vector<std::size_t>> &receivedPoints = exchange.ReceivedPoints();
            for(std::size_t k = 0; k < receivedFrom.size(); ++k)
            {
                std::vector<std::size_t> &ghosts = out.ghostIndices[rankIndex(receivedFrom[k])];
                for(const std::size_t &i : receivedPoints[k])
                {
                    ghosts.push_back(out.receivedIndex[i]);
                }
            }
            // the rank's points each rank has
            const std::vector<int> &sentTo = exchange.SentTo();
            const std::vector<std::vector<std::size_t>> &sentPoints = exchange.SentPoints();
            for(std::size_t k = 0; k < sentTo.size(); ++k)
            {
                std::vector<std::size_t> &points = out.duplicatedPoints[rankIndex(sentTo[k])];
                points.insert(points.end(), sentPoints[k].begin(), sentPoints[k].end());
            }
            // only the ranks both sent to and received from
            const std::vector<std::vector<int>> receivedLists = exchange.AllReceivedFrom();
            for(std::size_t i = 0; i < out.duplicatedProcs.size(); ++i)
            {
                const int rank = out.duplicatedProcs[i];
                const bool sent = std::find(sentTo.begin(), sentTo.end(), rank) != sentTo.end();
                const bool received = std::any_of(receivedLists.begin(), receivedLists.end(), [rank](const std::vector<int> &list) { return std::find(list.begin(), list.end(), rank) != list.end(); });
                if(not sent or not received)
                {
                    out.duplicatedProcs.erase(out.duplicatedProcs.begin() + i);
                    out.duplicatedPoints.erase(out.duplicatedPoints.begin() + i);
                    out.ghostIndices.erase(out.ghostIndices.begin() + i);
                    --i;
                }
            }
        }
    #endif // MADVORO_WITH_MPI

    // adds ghosts known before the build, telling the builder which of them mirror a cell (they lie beyond a face that is not periodic, and reflect onto a cell)
    void AddKnownGhosts(Part &part, const std::vector<PointT> &ghosts, const BuildDomain<PointT> &domain, std::unique_ptr<OctTree<IVec>> &cellsTree) const
    {
        std::vector<MirrorSource> sources(ghosts.size());
        bool anyMirror = false;
        const double width = abs(domain.ur - domain.ll);
        for(std::size_t i = 0; i < ghosts.size(); ++i)
        {
            const PointT &ghost = ghosts[i];
            if(part.cellsNum == 0 or not domain.IsOutside(ghost))
            {
                continue;
            }
            if(not cellsTree)
            {
                const PointT margin = (domain.ur - domain.ll) * 0.00001;
                cellsTree = std::make_unique<OctTree<IVec>>(IVec(domain.ll - margin, NO_INDEX), IVec(domain.ur + margin, NO_INDEX));
                for(std::size_t c = 0; c < part.cellsNum; ++c)
                {
                    cellsTree->insert(IVec(part.points[c].x, part.points[c].y, part.points[c].z, c));
                }
            }
            for(std::size_t face = 0; face < domain.box.size() and sources[i].cell == MirrorSource::NOT_A_MIRROR; ++face)
            {
                if(domain.HasPeriodicBoundaries() and domain.periodic[GetBoxFaceAxis(face)])
                {
                    continue;
                }
                const PointT reflected = MirrorPoint(domain.box[face], ghost);
                if(domain.IsOutside(reflected))
                {
                    continue;
                }
                const std::size_t cell = cellsTree->closestPoint(reflected).getIndex();
                // the mirror was made from the cell's point across the face, to round-off
                if(abs(reflected - part.points[cell]) <= 1e-6 * abs(ghost - part.points[cell]) + 1e-13 * width)
                {
                    sources[i].cell = cell;
                    sources[i].face = face;
                    anyMirror = true;
                }
            }
        }
        part.points.insert(part.points.end(), ghosts.begin(), ghosts.end());
        if(anyMirror)
        {
            part.builder->AddMirrors(ghosts, sources);
        }
        else
        {
            part.builder->Add(ghosts);
        }
    }

    // a cell with no radius from the last build starts from the distance to its nearest point
    void InitialRadii(const Part &part, const std::vector<PointT> &rankPoints, std::vector<double> &radii) const
    {
        // three times the nearest neighbor's distance brings a layer of ghosts in the first round (no cell is open towards other ranks when it is first built)
        constexpr double factor = 3;
        // the tree's queries only read it (each thread has its search stack)
        Geometry::HostParallelFor(part.cellsNum, [&](std::size_t j)
        {
            const std::size_t pointIdx = part.rankIndex.find(j)->second;
            if(radii[pointIdx] <= 0)
            {
                radii[pointIdx] = (part.cellsNum > 1) ? factor * this->rankTree->closestPointDistance(rankPoints[pointIdx], false) : abs(this->domainLL - this->domainUR);
            }
        });
    }

    // the part's cells and their moments, and the mirrored centroids of its mirror ghosts
    void FinishPart(Part &part, const BuildDomain<PointT> &domain) const
    {
        const auto isImage = [&](std::size_t p)
        {
            return p >= part.cellsNum and ImageCode(domain.PeriodicTranslation(part.points[p])) != ZERO_IMAGE_CODE;
        };
        std::vector<std::uint8_t> looseAreaCut(part.points.size());
        for(std::size_t p = 0; p < part.points.size(); ++p)
        {
            // a looser small-area cut against points outside the domain (not periodic images)
            looseAreaCut[p] = (domain.IsOutside(part.points[p]) and not isImage(p)) ? 1 : 0;
        }
        // the local builder's raw faces, finished by the shared stage
        Geometry::RawFacesView raw;
        part.builder->Finish(raw);
        const Geometry::FaceStageReport report = part.stage.Finish(raw, looseAreaCut, part.voronoi);
        if(report.Failed())
        {
            MadVoro::Exception::MadVoroException eo((report.result == Geometry::RawFaceResult::NoSecondEdge) ? "DistributedBuilder: a face has no two non degenerate edges" : "DistributedBuilder: a face separates a cell from a point that is not one (the local builder's bounding points)");
            eo.addEntry("cell", report.cell);
            eo.addEntry("point", part.points[report.cell]);
            eo.addEntry("other point index", report.other);
            throw eo;
        }
        const Geometry::FlatVoronoi &voronoi = part.voronoi;
        if(voronoi.CellCount() != part.cellsNum)
        {
            MadVoro::Exception::MadVoroException eo("DistributedBuilder: the local builder built a different number of cells");
            eo.addEntry("cells built", voronoi.CellCount());
            eo.addEntry("cells", part.cellsNum);
            throw eo;
        }
        part.CM.assign(part.points.size(), PointT());
        for(std::size_t c = 0; c < part.cellsNum; ++c)
        {
            part.CM[c] = PointT(voronoi.cellCM[c].x, voronoi.cellCM[c].y, voronoi.cellCM[c].z);
        }
        // the cells behind rigid boundary faces mirror their real neighbor
        for(std::size_t f = 0; f < voronoi.FaceCount(); ++f)
        {
            const std::size_t first = voronoi.faceNeighbor0[f], second = voronoi.faceNeighbor1[f];
            if(first < part.cellsNum and second < part.cellsNum)
            {
                continue;
            }
            const std::size_t real = (first >= part.cellsNum) ? second : first;
            const std::size_t other = (first >= part.cellsNum) ? first : second;
            if(isImage(other))
            {
                continue;
            }
            #ifdef MADVORO_WITH_MPI
                if(not domain.IsOutside(part.points[other]))
                {
                    continue;
                }
            #endif // MADVORO_WITH_MPI
            const Geometry::FlatPoint &vertex = voronoi.vertices[voronoi.faceVertices[voronoi.faceVertexOffsets[f]]];
            part.CM[other] = Geometry::MirroredCM(part.points[first], part.points[second], part.CM[real], PointT(vertex.x, vertex.y, vertex.z));
        }
    }

    // the build of one domain: its points (taken from the part), 4 unused points, then its ghosts
    void OneDomainOutput(Part &part, DistributedBuild<PointT> &out) const
    {
        const Geometry::FlatVoronoi &voronoi = part.voronoi;
        const std::size_t cellsNum = part.cellsNum;
        const std::size_t pointsNum = part.points.size();
        const auto meshIndex = [cellsNum](std::size_t p)
        {
            return (p < cellsNum) ? p : p + BOUNDING_POINTS;
        };
        out.points = std::move(part.points);
        out.points.insert(out.points.begin() + cellsNum, BOUNDING_POINTS, PointT());
        this->UnusedPoints(out.points.begin() + cellsNum);
        out.ghostRankIndex.clear();
        for(const auto &[p, pointIdx] : part.rankIndex)
        {
            if(p >= cellsNum)
            {
                out.ghostRankIndex.emplace_hint(out.ghostRankIndex.end(), meshIndex(p), pointIdx);
            }
        }
        const std::size_t facesNum = voronoi.FaceCount();
        out.cellFaces.resize(cellsNum);
        out.faceNeighbors.resize(facesNum);
        out.faceVertices.resize(facesNum);
        out.faceArea.resize(facesNum);
        out.faceCM.resize(facesNum);
        out.vertices.resize(voronoi.vertices.size());
        out.volume.resize(cellsNum);
        out.CM.assign(out.points.size(), PointT());
        Geometry::HostParallelFor(cellsNum, [&](std::size_t c)
        {
            out.cellFaces[c].assign(voronoi.cellFaces.begin() + voronoi.cellFaceOffsets[c], voronoi.cellFaces.begin() + voronoi.cellFaceOffsets[c + 1]);
            out.volume[c] = voronoi.cellVolume[c];
        });
        Geometry::HostParallelFor(facesNum, [&](std::size_t f)
        {
            out.faceVertices[f].assign(voronoi.faceVertices.begin() + voronoi.faceVertexOffsets[f], voronoi.faceVertices.begin() + voronoi.faceVertexOffsets[f + 1]);
            out.faceNeighbors[f] = std::make_pair(meshIndex(voronoi.faceNeighbor0[f]), meshIndex(voronoi.faceNeighbor1[f]));
            out.faceArea[f] = voronoi.faceArea[f];
            out.faceCM[f] = PointT(voronoi.faceCM[f].x, voronoi.faceCM[f].y, voronoi.faceCM[f].z);
        });
        Geometry::HostParallelFor(voronoi.vertices.size(), [&](std::size_t v)
        {
            out.vertices[v] = PointT(voronoi.vertices[v].x, voronoi.vertices[v].y, voronoi.vertices[v].z);
        });
        for(std::size_t p = 0; p < pointsNum; ++p)
        {
            out.CM[meshIndex(p)] = part.CM[p];
        }
    }

    // the 4 unused points after the cells: far outside the domain, apart from each other
    template <typename Iterator>
    void UnusedPoints(Iterator first) const
    {
        const PointT width = this->domainUR - this->domainLL;
        for(std::size_t k = 0; k < BOUNDING_POINTS; ++k)
        {
            *(first + k) = this->domainLL - width * (1000.0 * static_cast<double>(k + 1));
        }
    }

    void BuildOneDomain(const BuildDomain<PointT> &domain, const std::vector<PointT> &rankPoints, const std::vector<std::size_t> &cells, std::vector<double> &radii, const std::vector<PointT> *carriedGhosts, DistributedBuild<PointT> &out)
    {
        // MADVORO_GHOST_TIMING=1: the build's phases (rank 0)
        static const bool timing = std::getenv("MADVORO_GHOST_TIMING") != nullptr;
        auto mark = std::chrono::high_resolution_clock::now();
        std::vector<double> phases;
        const auto lap = [&](void)
        {
            const auto now = std::chrono::high_resolution_clock::now();
            phases.push_back(std::chrono::duration<double>(now - mark).count());
            mark = now;
        };
        const std::size_t cellsNum = cells.size();
        Part &part = this->SinglePart();
        part.points.resize(cellsNum);
        part.cellsNum = cellsNum;
        part.rankIndex.clear();
        for(std::size_t k = 0; k < cellsNum; ++k)
        {
            part.points[k] = rankPoints[cells[k]];
            part.rankIndex.emplace_hint(part.rankIndex.end(), k, cells[k]);
        }
        this->InitialRadii(part, rankPoints, radii);
        lap();
        part.builder->Build(part.points, domain, HilbertOrder3D(part.points));
        lap();
        out.carriedGhostsIndex = cellsNum + BOUNDING_POINTS;
        if(carriedGhosts != nullptr)
        {
            part.points.insert(part.points.end(), carriedGhosts->begin(), carriedGhosts->end());
            part.builder->Add(*carriedGhosts);
        }

        OctTreeFinder<PointT> finder(this->rankTree.get(), rankPoints);
        GhostSearch<PointT> search(&finder, &finder);
        search.Begin(domain.box, domain.normals, domain.periodic, domain.ll, domain.ur, part.rankIndex, cellsNum, radii);
        PartMesh mesh(part, rankPoints, nullptr, radii);
        std::vector<GhostDomain<PointT>> domains = {{&search, &mesh}};
        const auto oneByOne = [](std::size_t domainsNum, const auto &task)
        {
            for(std::size_t d = 0; d < domainsNum; ++d)
            {
                task(d);
            }
        };
        #ifdef MADVORO_WITH_MPI
            // other ranks' points go into the mesh as they arrive
            FunctionGhostReceiver<PointT> receiver(false, [&](const std::vector<PointT> &points)
            {
                for(std::size_t i = 0; i < points.size(); ++i)
                {
                    out.receivedIndex.push_back(part.points.size() + BOUNDING_POINTS + i);
                }
                part.points.insert(part.points.end(), points.begin(), points.end());
                part.builder->Add(points);
                return std::size_t(0);
            });
            out.rounds = BringGhosts(domains, this->Exchange(), &receiver, this->exchange ? this->exchangeComm : MPI_COMM_SELF, oneByOne);
        #else // MADVORO_WITH_MPI
            out.rounds = BringGhosts(domains, oneByOne);
        #endif // MADVORO_WITH_MPI
        lap();
        this->FinishPart(part, domain);
        lap();
        this->OneDomainOutput(part, out);
        lap();
        int rank = 0;
        #ifdef MADVORO_WITH_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        #endif // MADVORO_WITH_MPI
        if(timing and rank == 0)
        {
            std::printf("Build phases: trees %.4f radii %.4f start %.4f rounds %.4f finish %.4f output %.4f s\n", this->treesTime, phases[0], phases[1], phases[2], phases[3], phases[4]);
        }
    }

    // builtOrder: the built points' Hilbert order
    void BuildSubdomains(const BuildDomain<PointT> &domain, std::size_t count, const std::vector<PointT> &rankPoints, const std::vector<std::size_t> &cells, const std::vector<std::size_t> &builtOrder, std::vector<double> &radii, const std::vector<PointT> *knownGhosts, const std::vector<PointT> *carriedGhosts, DistributedBuild<PointT> &out)
    {
        const std::size_t pointsNum = cells.size(); // the points built
        const std::size_t rankPointsNum = rankPoints.size(); // the rank's points (a partial build does not build all)
        bool builtInOrder = (pointsNum == rankPointsNum);
        for(std::size_t k = 0; k < pointsNum and builtInOrder; ++k)
        {
            builtInOrder = (cells[k] == k);
        }
        std::vector<PointT> builtPoints;
        if(not builtInOrder)
        {
            builtPoints.reserve(pointsNum);
            for(const std::size_t &pointIdx : cells)
            {
                builtPoints.push_back(rankPoints[pointIdx]);
            }
        }
        // subdomain s is the range [s * pointsNum / count, (s + 1) * pointsNum / count) of the built points' Hilbert order
        const auto subdomainBegin = [pointsNum, count](std::size_t s) { return s * pointsNum / count; };
        // the subdomain of each of the rank's points (`count` for points not built)
        std::vector<subdomain_t> subdomainOf(rankPointsNum, static_cast<subdomain_t>(count));
        for(std::size_t s = 0; s < count; ++s)
        {
            for(std::size_t i = subdomainBegin(s); i < subdomainBegin(s + 1); ++i)
            {
                subdomainOf[cells[builtOrder[i]]] = static_cast<subdomain_t>(s);
            }
        }
        const std::shared_ptr<SubdomainGrid<PointT>> grid = std::make_shared<SubdomainGrid<PointT>>(rankPoints, subdomainOf);
        // other ranks also send periodic images, up to a period outside the box
        const PointT period = domain.HasPeriodicBoundaries() ? (domain.ur - domain.ll) : PointT();
        const std::shared_ptr<RankGhostPoints<PointT>> otherRankPoints = std::make_shared<RankGhostPoints<PointT>>(domain.ll - period, domain.ur + period, rankPointsNum);
        const std::shared_ptr<OctTreeFinder<PointT>> pointsFinder = std::make_shared<OctTreeFinder<PointT>>(this->rankTree.get(), rankPoints);
        const auto addOtherRankPoints = [&](const std::vector<PointT> &received)
        {
            for(const PointT &point : received)
            {
                otherRankPoints->Add(point);
            }
            grid->AddOtherRankPoints(received.begin(), received.end());
        };

        this->single = Part();
        this->parts.resize(count);
        std::vector<std::exception_ptr> errors(count);
        const auto forEachSubdomain = [&](const auto &task)
        {
            const auto run = [&](std::size_t s)
            {
                try
                {
                    task(s);
                }
                catch(...)
                {
                    errors[s] = std::current_exception();
                }
            };
            // host threads also in a device build (Kokkos cannot nest parallel loops)
            Geometry::HostParallelFor(count, [&run](std::size_t s)
            {
                const Geometry::SerialGeometryScope serialGeometry;
                run(s);
            });
            for(const std::exception_ptr &error : errors)
            {
                if(error)
                {
                    std::rethrow_exception(error);
                }
            }
        };

        std::vector<std::vector<std::size_t>> subdomainCells(count);
        for(std::size_t s = 0; s < count; ++s)
        {
            // in the order of the build, so one subdomain is the same build as one domain
            std::vector<std::size_t> built(builtOrder.begin() + subdomainBegin(s), builtOrder.begin() + subdomainBegin(s + 1));
            std::sort(built.begin(), built.end());
            for(const std::size_t &k : built)
            {
                subdomainCells[s].push_back(cells[k]);
            }
        }
        forEachSubdomain([&](std::size_t s)
        {
            Part &part = this->parts[s];
            if(not part.builder)
            {
                part.builder = this->makeLocal(true);
            }
            part.finder = std::make_shared<SubdomainRangeFinder<PointT>>(pointsFinder, otherRankPoints, grid, static_cast<subdomain_t>(s));
            part.search = std::make_unique<GhostSearch<PointT>>(part.finder.get(), pointsFinder.get());
            part.cellsNum = subdomainCells[s].size();
            part.points.resize(part.cellsNum);
            part.rankIndex.clear();
            for(std::size_t j = 0; j < part.cellsNum; ++j)
            {
                part.points[j] = rankPoints[subdomainCells[s][j]];
                part.rankIndex.emplace_hint(part.rankIndex.end(), j, subdomainCells[s][j]);
            }
            this->InitialRadii(part, rankPoints, radii);
            part.builder->Build(part.points, domain, HilbertOrder3D(part.points));
            part.search->Begin(domain.box, domain.normals, domain.periodic, domain.ll, domain.ur, part.rankIndex, part.cellsNum, radii);
        });

        std::size_t knownNum = 0;
        if(knownGhosts != nullptr)
        {
            knownNum = knownGhosts->size();
            addOtherRankPoints(*knownGhosts);
            // each known ghost goes at once to the subdomain of its nearest point, so the rounds find little
            std::vector<std::vector<std::size_t>> ghostsOf(count);
            for(std::size_t i = 0; i < knownGhosts->size(); ++i)
            {
                const std::size_t nearest = this->rankTree->closestPoint((*knownGhosts)[i]).getIndex();
                if(subdomainOf[nearest] < count)
                {
                    ghostsOf[subdomainOf[nearest]].push_back(rankPointsNum + i);
                }
            }
            forEachSubdomain([&](std::size_t s)
            {
                Part &part = this->parts[s];
                std::vector<PointT> points;
                for(const std::size_t &pointIdxInAll : ghostsOf[s])
                {
                    if(part.search->MarkBrought(pointIdxInAll))
                    {
                        part.rankIndex.emplace_hint(part.rankIndex.end(), part.points.size() + points.size(), pointIdxInAll);
                        points.push_back(otherRankPoints->Points()[pointIdxInAll - rankPointsNum]);
                    }
                }
                part.points.insert(part.points.end(), points.begin(), points.end());
                part.builder->Add(points);
            });
        }
        out.carriedGhostsIndex = pointsNum + BOUNDING_POINTS + knownNum;
        if(carriedGhosts != nullptr)
        {
            addOtherRankPoints(*carriedGhosts);
        }

        std::vector<std::unique_ptr<PartMesh>> meshes;
        std::vector<GhostDomain<PointT>> domains;
        for(std::size_t s = 0; s < count; ++s)
        {
            meshes.push_back(std::make_unique<PartMesh>(this->parts[s], rankPoints, otherRankPoints.get(), radii));
            domains.push_back({this->parts[s].search.get(), meshes.back().get()});
        }
        const auto inParallel = [&](std::size_t, const auto &task) { forEachSubdomain(task); };
        #ifdef MADVORO_WITH_MPI
            // other ranks' points are kept for all the subdomains
            FunctionGhostReceiver<PointT> receiver(true, [&](const std::vector<PointT> &points)
            {
                const std::size_t first = rankPointsNum + otherRankPoints->Points().size();
                for(std::size_t i = 0; i < points.size(); ++i)
                {
                    out.receivedIndex.push_back(pointsNum + BOUNDING_POINTS + otherRankPoints->Points().size() + i);
                }
                addOtherRankPoints(points);
                return first;
            });
            out.rounds = BringGhosts(domains, this->Exchange(), &receiver, this->exchange ? this->exchangeComm : MPI_COMM_SELF, inParallel);
        #else // MADVORO_WITH_MPI
            out.rounds = BringGhosts(domains, inParallel);
        #endif // MADVORO_WITH_MPI

        forEachSubdomain([&](std::size_t s) { this->FinishPart(this->parts[s], domain); });

        std::vector<SubdomainTessellation<PointT>> tessellations;
        for(const Part &part : this->parts)
        {
            tessellations.push_back({&part.points, part.cellsNum, &part.rankIndex, &part.voronoi, &part.CM});
        }
        // the built points in build order, then the other ranks' points after 4 unused points
        const std::size_t otherRankPointsNum = otherRankPoints->Points().size();
        std::vector<PointT> points = builtInOrder ? rankPoints : builtPoints;
        points.resize(pointsNum + BOUNDING_POINTS);
        this->UnusedPoints(points.begin() + pointsNum);
        points.insert(points.end(), otherRankPoints->Points().begin(), otherRankPoints->Points().end());
        std::vector<std::size_t> mergedIndexOf(rankPointsNum + otherRankPointsNum, NO_INDEX);
        for(std::size_t k = 0; k < pointsNum; ++k)
        {
            mergedIndexOf[cells[k]] = k;
        }
        for(std::size_t i = 0; i < otherRankPointsNum; ++i)
        {
            mergedIndexOf[rankPointsNum + i] = pointsNum + BOUNDING_POINTS + i;
        }
        MergeSubdomains(tessellations, domain, points, pointsNum, mergedIndexOf, rankPointsNum, out);
        // the searches hold the received points and the grid
        for(Part &part : this->parts)
        {
            part.search.reset();
            part.finder.reset();
        }
    }

    LocalBuilderFactory makeLocal;
    double treesTime = 0;
    Part single; // the builder of one domain
    std::vector<Part> parts; // one per thread subdomain
    std::shared_ptr<OctTree<IVec>> rankTree, builtTree;
    PointT domainLL, domainUR;
    #ifdef MADVORO_WITH_MPI
        std::optional<GhostExchange<PointT>> exchange;
        MPI_Comm exchangeComm = MPI_COMM_WORLD;
    #endif // MADVORO_WITH_MPI
};

} // namespace MadVoro

#endif // MADVORO_BUILDER_DISTRIBUTED_BUILDER_HPP
