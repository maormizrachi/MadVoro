// The search for a mesh's ghost points, round by round (arepo paper, section 2.4)
#ifndef MADVORO_GHOSTS_GHOST_SEARCH_HPP
#define MADVORO_GHOSTS_GHOST_SEARCH_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <set>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container_hash/hash.hpp>
#include <boost/unordered/unordered_flat_set.hpp>

#include "types.hpp"
#include "elementary/Face3D.hpp"
#include "elementary/PointOps.hpp"
#include "geometry/Intersections.hpp"
#include "range/SmallRangeAgent.hpp"
#include "range/BigRangeAgent.hpp"
#include "exception/MadVoroException.hpp"

namespace MadVoro
{

template <typename PointT>
using PeriodicIgnoreSets = std::array<std::unordered_set<size_t>, NUM_IMAGE_CODES>;

template <typename PointT>
inline bool SphereIntersectsAABB(const PointT &center, typename PointT::coord_type radius, const PointT &tileLL, const PointT &tileUR)
{
    using coord_type = typename PointT::coord_type;
    coord_type d2 = coord_type(0);
    if(center.x < tileLL.x)
    {
        d2 += (tileLL.x - center.x) * (tileLL.x - center.x);
    }
    else if(center.x > tileUR.x)
    {
        d2 += (center.x - tileUR.x) * (center.x - tileUR.x);
    }
    if(center.y < tileLL.y)
    {
        d2 += (tileLL.y - center.y) * (tileLL.y - center.y);
    }
    else if(center.y > tileUR.y)
    {
        d2 += (center.y - tileUR.y) * (center.y - tileUR.y);
    }
    if(center.z < tileLL.z)
    {
        d2 += (tileLL.z - center.z) * (tileLL.z - center.z);
    }
    else if(center.z > tileUR.z)
    {
        d2 += (center.z - tileUR.z) * (center.z - tileUR.z);
    }
    constexpr coord_type eps = std::numeric_limits<coord_type>::epsilon();
    return d2 <= radius * radius * (coord_type(1) + coord_type(64) * eps);
}

inline int GetBoxFaceAxis(size_t faceIdx)
{
    static const int axisMap[] = {2, 1, 0, 1, 0, 2};
    return axisMap[faceIdx % 6];
}

template <typename QueryDataType, typename PointT>
inline void ExpandPeriodicQueries(const std::vector<QueryDataType> &baseQueries, const std::array<bool, 3> &periodic, const PointT &ll, const PointT &ur, std::vector<QueryDataType> &expanded)
{
    PointT L = ur - ll;
    expanded.clear();
    expanded.reserve(baseQueries.size() * 8);
    for(const QueryDataType &baseQuery : baseQueries)
    {
        expanded.push_back(baseQuery);
        std::array<int, 3> shiftRanges[3] = {{0}, {0}, {0}};
        if(periodic[0])
        {
            shiftRanges[0] = {-1, 0, 1};
        }
        if(periodic[1])
        {
            shiftRanges[1] = {-1, 0, 1};
        }
        if(periodic[2])
        {
            shiftRanges[2] = {-1, 0, 1};
        }
        for(int sx : shiftRanges[0])
        {
            for(int sy : shiftRanges[1])
            {
                for(int sz : shiftRanges[2])
                {
                    if(sx == 0 and sy == 0 and sz == 0)
                    {
                        continue;
                    }
                    PointT translation(L.x * sx, L.y * sy, L.z * sz);
                    PointT tileLL = ll + translation;
                    PointT tileUR = ur + translation;
                    if(not SphereIntersectsAABB(baseQuery.center, baseQuery.radius, tileLL, tileUR))
                    {
                        continue;
                    }
                    QueryDataType periodicQuery = baseQuery;
                    periodicQuery.imageTranslation = translation;
                    periodicQuery.center = baseQuery.center - translation;
                    if constexpr(std::is_same_v<QueryDataType, BigRangeQueryData<PointT>>)
                    {
                        periodicQuery.originalPoint = baseQuery.originalPoint - translation;
                        periodicQuery.askOnlyClose = false;
                    }
                    expanded.push_back(periodicQuery);
                }
            }
        }
    }
}

template <typename PointT>
PointT MirrorPoint(Face3D<PointT> const &face, PointT const &point)
{
    using namespace MadVoro::fallback;
    PointT normal = CrossProduct(face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]);
    normal = normal / abs(normal);
    return point - (2 * ScalarProd(point - face.vertices[0], normal)) * normal;
}

//! \brief The box faces a sphere intersects
template <typename PointT>
void CheckToMirror(const Sphere<PointT> &sphere, const std::vector<Face3D<PointT>> &box, const std::vector<PointT> &normals, std::vector<size_t> &facesItCuts)
{
    facesItCuts.clear();
    for(size_t i = 0; i < box.size(); i++)
    {
        if(FaceSphereIntersections(box[i], sphere, normals[i]))
        {
            facesItCuts.push_back(i);
        }
    }
}

//! \brief The (face, point) pairs to mirror: each query's point across each box face its sphere cuts
template<typename QueryDataType, typename PointT>
std::vector<std::pair<size_t, size_t>> MirrorPoints(const std::vector<QueryDataType> &queries, const std::vector<Face3D<PointT>> &box, const std::vector<PointT> &normals)
{
    static_assert(std::is_convertible<QueryDataType*, RangeQueryData<PointT>*>::value, "MirrorPoints: QueryDataType must inherit 'RangeQueryData'");

    std::vector<std::pair<size_t, size_t>> mirroredPoints;
    std::vector<size_t> facesItCuts;
    for(const QueryDataType &query : queries)
    {
        Sphere<PointT> sphere(PointT(query.center), query.radius);
        CheckToMirror(sphere, box, normals, facesItCuts);
        for(const size_t &faceIdx : facesItCuts)
        {
            mirroredPoints.push_back(std::make_pair(faceIdx, query.pointIdx));
        }
    }
    return mirroredPoints;
}

//! \brief The mesh a ghost search builds: what the search asks of it, and the points it adds to it
template <typename PointT>
class GhostSearchMesh
{
public:
    virtual ~GhostSearchMesh() = default;

    //! \brief A point of the mesh (its cells' points first)
    virtual const PointT &Point(std::size_t index) const = 0;

    //! \brief A cell's point's index among all the rank's points
    virtual std::size_t AllPointIndex(std::size_t cell) const = 0;

    //! \brief How many of the rank's points this mesh can add (by index among all the rank's points)
    virtual std::size_t RankPointsNum(void) const = 0;

    //! \brief The current vertices of a cell
    virtual void CellVertices(std::size_t cell, std::vector<PointT> &vertices) = 0;

    //! \brief The distance from a cell's point to its farthest current vertex
    virtual double CellRadius(std::size_t cell) = 0;

    //! \brief Brings the current vertices of these cells at once, for the CellVertices calls that follow
    virtual void FetchVertices(const std::vector<std::size_t> &)
    {
    }

    //! \brief The search radius a cell starts its next build with
    virtual void SetCarriedRadius(std::size_t cell, double radius) = 0;

    //! \brief Adds points of the rank (index among all the rank's points, and the translation of a periodic image)
    virtual void AddRankPoints(const std::vector<std::pair<std::size_t, PointT>> &points) = 0;

    //! \brief Adds mirrors: points[i] mirrors cell sources[i].second across box face sources[i].first
    virtual void AddMirrors(const std::vector<PointT> &points, const std::vector<std::pair<std::size_t, std::size_t>> &sources) = 0;
};

//! \brief Which of a round's searches to answer from the rank's points
enum class LocalSearch
{
    All,
    Direct, // not periodic images
    Images
};

//! \brief A domain's ghost search: the radii it searches, what it brought and mirrored, and the current round
template <typename PointT>
class GhostSearch
{
public:
    static constexpr std::size_t MAX_POINTS_TO_GET = 15;
    static constexpr double RADIUS_GROWTH = 1.1;
    static constexpr double LARGE_RADIUS_SHRINK = 0.95;

    // `finder` searches the rank's points; `imageFinder`, for the box's periodic images
    GhostSearch(const RangeFinder<PointT> *finder, const RangeFinder<PointT> *imageFinder)
        #ifdef MADVORO_WITH_MPI
            : bigAgent(finder, nullptr, unusedSent, MPI_COMM_SELF), smallAgent(finder, nullptr, unusedSent, MPI_COMM_SELF), imageBigAgent(imageFinder, nullptr, unusedSent, MPI_COMM_SELF), imageSmallAgent(imageFinder, nullptr, unusedSent, MPI_COMM_SELF)
        #else // MADVORO_WITH_MPI
            : bigAgent(finder), smallAgent(finder), imageBigAgent(imageFinder), imageSmallAgent(imageFinder)
        #endif // MADVORO_WITH_MPI
    {
    }

    GhostSearch(const GhostSearch &) = delete;
    GhostSearch &operator=(const GhostSearch &) = delete;

    //! \brief Starts a search for the cells `cells` lists (cell to index among all the rank's points), from their radii
    void Begin(const std::vector<Face3D<PointT>> &boxFaces, const std::vector<PointT> &boxNormals, const std::array<bool, 3> &periodicAxes, const PointT &boxLL, const PointT &boxUR, const AllPointsMap &cells, std::size_t cellsNum, const std::vector<double> &radii)
    {
        this->box = boxFaces;
        this->normals = boxNormals;
        this->periodic = periodicAxes;
        this->ll = boxLL;
        this->ur = boxUR;
        this->smallPoints.clear();
        this->largePoints.clear();
        this->firstLargeIteration.clear();
        this->currentRadiuses.assign(cellsNum, RADIUS_UNINITIALIZED);
        this->ignorePoints = PeriodicIgnoreSets<PointT>();
        for(const std::pair<size_t, size_t> &indices : cells)
        {
            if(indices.first >= cellsNum)
            {
                std::cout << "Error: given point index " << indices.first << ", while built only with " << cellsNum << " points" << std::endl;
            }
            this->smallPoints.insert(indices.first);
            this->currentRadiuses[indices.first] = radii[indices.second];
            this->ignorePoints[ZERO_IMAGE_CODE].insert(indices.second);
        }
        this->mirrored.clear();
        this->checkedVertices.clear();
        this->iterations = 0;
    }

    //! \brief A round's searches and mirrors; returns whether it searches. `reachesOtherRanks` (or null) marks the searches other ranks answer
    template <typename Reach = std::nullptr_t>
    bool PrepareRound(GhostSearchMesh<PointT> &mesh, const Reach *reachesOtherRanks = nullptr)
    {
        this->iterations++;
        this->CreateQueries(mesh);
        this->mirroredPoints = MirrorPoints(this->smallQueries, this->box, this->normals);
        const std::vector<std::pair<size_t, size_t>> moreMirroredPoints = MirrorPoints(this->bigQueries, this->box, this->normals);
        this->mirroredPoints.insert(this->mirroredPoints.end(), moreMirroredPoints.begin(), moreMirroredPoints.end());
        if(this->HasPeriodicBoundaries())
        {
            // no mirrors across periodic faces; the searches also search the periodic images
            this->mirroredPoints.erase(std::remove_if(this->mirroredPoints.begin(), this->mirroredPoints.end(), [this](const std::pair<size_t, size_t> &mirror) { return this->IsPeriodicFace(mirror.first); }), this->mirroredPoints.end());
            std::vector<SmallRangeQueryData<PointT>> expandedSmallQueries;
            std::vector<BigRangeQueryData<PointT>> expandedBigQueries;
            ExpandPeriodicQueries(this->smallQueries, this->periodic, this->ll, this->ur, expandedSmallQueries);
            ExpandPeriodicQueries(this->bigQueries, this->periodic, this->ll, this->ur, expandedBigQueries);
            this->smallQueries = std::move(expandedSmallQueries);
            this->bigQueries = std::move(expandedBigQueries);
        }
        this->smallResults.clear();
        this->bigResults.clear();
        for(const SmallRangeQueryData<PointT> &query : this->smallQueries)
        {
            this->smallResults[query.pointIdx] = 0;
        }
        for(const BigRangeQueryData<PointT> &query : this->bigQueries)
        {
            this->bigResults[query.pointIdx] = 0;
        }
        this->remoteSmall.clear();
        this->remoteBig.clear();
        if constexpr(not std::is_same_v<Reach, std::nullptr_t>)
        {
            if(reachesOtherRanks != nullptr)
            {
                for(size_t j = 0; j < this->smallQueries.size(); ++j)
                {
                    if(reachesOtherRanks->ReachesOtherRanks(this->smallQueries[j]))
                    {
                        this->remoteSmall.push_back(j);
                    }
                }
                for(size_t j = 0; j < this->bigQueries.size(); ++j)
                {
                    if(reachesOtherRanks->ReachesOtherRanks(this->bigQueries[j]))
                    {
                        this->remoteBig.push_back(j);
                    }
                }
            }
        }
        return not this->NothingToSearch();
    }

    bool NothingToSearch(void) const
    {
        return this->smallQueries.empty() and this->bigQueries.empty();
    }

    const std::vector<SmallRangeQueryData<PointT>> &SmallQueries(void) const
    {
        return this->smallQueries;
    }

    const std::vector<BigRangeQueryData<PointT>> &BigQueries(void) const
    {
        return this->bigQueries;
    }

    //! \brief The searches other ranks answer (positions among the round's searches)
    const std::vector<size_t> &RemoteSmall(void) const
    {
        return this->remoteSmall;
    }

    const std::vector<size_t> &RemoteBig(void) const
    {
        return this->remoteBig;
    }

    //! \brief Whether a search in the rank's points (of `allPointsNum`) can still find a point
    bool CanFindLocalPoints(std::size_t allPointsNum) const
    {
        return this->ignorePoints[ZERO_IMAGE_CODE].size() < allPointsNum;
    }

    //! \brief Answers searches from the rank's points: what they found (index among all the rank's points, translation), small searches first
    std::vector<std::pair<std::size_t, PointT>> SearchLocally(LocalSearch which)
    {
        std::vector<std::pair<std::size_t, PointT>> found;
        if(which == LocalSearch::All)
        {
            this->Answer(this->smallQueries, this->smallAgent, this->smallResults, found);
            this->Answer(this->bigQueries, this->bigAgent, this->bigResults, found);
            return found;
        }
        const bool images = (which == LocalSearch::Images);
        std::vector<SmallRangeQueryData<PointT>> smallSelected;
        std::vector<BigRangeQueryData<PointT>> bigSelected;
        this->Select(images, smallSelected, bigSelected);
        if(images and smallSelected.empty() and bigSelected.empty())
        {
            return found;
        }
        this->Answer(smallSelected, images ? this->imageSmallAgent : this->smallAgent, this->smallResults, found);
        this->Answer(bigSelected, images ? this->imageBigAgent : this->bigAgent, this->bigResults, found);
        return found;
    }

    //! \brief Checks (in debug builds) that searches in the rank's points find nothing (`All` or `Direct` ones)
    void CheckNothingLocally(LocalSearch which)
    {
        #ifndef NDEBUG
            std::vector<SmallRangeQueryData<PointT>> smallSelected = this->smallQueries;
            std::vector<BigRangeQueryData<PointT>> bigSelected = this->bigQueries;
            if(which != LocalSearch::All)
            {
                this->Select(false, smallSelected, bigSelected);
            }
            PeriodicIgnoreSets<PointT> ignoreCopy = this->ignorePoints;
            for(const std::vector<size_t> &answer : this->smallAgent.selfBatchAnswer(smallSelected, ignoreCopy))
            {
                assert(answer.empty());
            }
            for(const std::vector<size_t> &answer : this->bigAgent.selfBatchAnswer(bigSelected, ignoreCopy))
            {
                assert(answer.empty());
            }
        #else // NDEBUG
            (void)which;
        #endif // NDEBUG
    }

    //! \brief Counts points other ranks sent for a cell's searches
    void CountFound(bool big, std::size_t cell, std::size_t count)
    {
        (big ? this->bigResults : this->smallResults)[cell] += count;
    }

    //! \brief Points other ranks sent for one of this round's searches (indices among all the rank's points), kept for the round's end
    void AddRemoteAnswer(bool big, std::size_t query, std::vector<std::size_t> points)
    {
        this->remoteAnswers.push_back({big, query, std::move(points)});
    }

    //! \brief The points other ranks sent this round (big, then small searches'), counted
    std::array<std::vector<std::size_t>, 2> TakeRemoteAnswers(void)
    {
        std::array<std::vector<std::size_t>, 2> received;
        for(const RemoteAnswer &answer : this->remoteAnswers)
        {
            std::vector<std::size_t> &into = received[answer.big ? 0 : 1];
            into.insert(into.end(), answer.points.begin(), answer.points.end());
            const std::size_t cell = answer.big ? this->bigQueries[answer.query].pointIdx : this->smallQueries[answer.query].pointIdx;
            this->CountFound(answer.big, cell, answer.points.size());
        }
        this->remoteAnswers.clear();
        return received;
    }

    //! \brief Records a point of the rank as brought; returns whether it was new
    bool MarkBrought(std::size_t pointIdxInAll)
    {
        return this->ignorePoints[ZERO_IMAGE_CODE].insert(pointIdxInAll).second;
    }

    //! \brief The mirrors this round adds (each point once across each face), and the (face, point) each mirrors
    std::vector<PointT> NewMirrors(const GhostSearchMesh<PointT> &mesh, std::vector<std::pair<size_t, size_t>> &sources)
    {
        std::vector<PointT> newPoints;
        newPoints.reserve(this->mirroredPoints.size());
        sources.clear();
        for(const std::pair<size_t, size_t> &pairFacePoint : this->mirroredPoints)
        {
            if(this->mirrored.insert(pairFacePoint).second)
            {
                newPoints.push_back(MirrorPoint(this->box[pairFacePoint.first], mesh.Point(pairFacePoint.second)));
                sources.push_back(pairFacePoint);
            }
        }
        return newPoints;
    }

    //! \brief Ends a round: the cells the next round searches for, and the radii carried to the next build
    void EndRound(GhostSearchMesh<PointT> &mesh)
    {
        boost::container::flat_set<size_t> newSmallPoints, newLargePoints;
        for(const std::pair<size_t, size_t> &cellResult : this->smallResults)
        {
            const size_t &cell = cellResult.first;
            if(cellResult.second > MAX_POINTS_TO_GET)
            {
                // too many points: search for this cell as large from the next round
                newLargePoints.insert(cell);
                this->firstLargeIteration[cell] = this->iterations + 1;
                mesh.SetCarriedRadius(cell, LARGE_RADIUS_SHRINK * this->currentRadiuses[cell]);
            }
            else
            {
                const double maxRadius = mesh.CellRadius(cell);
                if(this->currentRadiuses[cell] < 2 * maxRadius)
                {
                    this->currentRadiuses[cell] *= RADIUS_GROWTH;
                    newSmallPoints.insert(cell);
                }
                else
                {
                    mesh.SetCarriedRadius(cell, RADIUS_GROWTH * (2 * maxRadius));
                }
            }
        }
        for(const std::pair<size_t, size_t> &cellResult : this->bigResults)
        {
            // a large cell is done once its searches find nothing
            if((this->iterations == this->firstLargeIteration.at(cellResult.first)) or (cellResult.second > 0))
            {
                newLargePoints.insert(cellResult.first);
            }
        }
        this->smallPoints = std::move(newSmallPoints);
        this->largePoints = std::move(newLargePoints);
    }

    std::size_t Iterations(void) const
    {
        return this->iterations;
    }

    std::size_t SmallCount(void) const
    {
        return this->smallPoints.size();
    }

    std::size_t LargeCount(void) const
    {
        return this->largePoints.size();
    }

private:
    struct RemoteAnswer
    {
        bool big;
        std::size_t query;
        std::vector<std::size_t> points;
    };

    bool HasPeriodicBoundaries(void) const
    {
        return this->periodic[0] or this->periodic[1] or this->periodic[2];
    }

    bool IsPeriodicFace(size_t faceIdx) const
    {
        if(not this->HasPeriodicBoundaries() or faceIdx >= 6)
        {
            return false;
        }
        return this->periodic[GetBoxFaceAxis(faceIdx)];
    }

    void CreateQueries(GhostSearchMesh<PointT> &mesh)
    {
        this->smallQueries.clear();
        this->bigQueries.clear();
        if(this->iterations > 1 and not this->largePoints.empty())
        {
            // on its first round as large a cell asks only the close ranks
            std::vector<std::pair<std::size_t, bool>> largeCells;
            largeCells.reserve(this->largePoints.size());
            for(const size_t &cell : this->largePoints)
            {
                largeCells.emplace_back(cell, this->iterations == this->firstLargeIteration.at(cell));
            }
            mesh.FetchVertices(std::vector<std::size_t>(this->largePoints.begin(), this->largePoints.end()));
            // a large cell is final once the ball around each of its vertices through its point holds no other point
            std::vector<std::array<double, 3>> searched;
            for(const auto &[cell, askOnlyClose] : largeCells)
            {
                const PointT &point = mesh.Point(cell);
                mesh.CellVertices(cell, this->vertices);
                for(const PointT &vertex : this->vertices)
                {
                    const std::array<double, 3> key = {vertex.x, vertex.y, vertex.z};
                    if(this->checkedVertices.contains(key))
                    {
                        continue;
                    }
                    if(not askOnlyClose)
                    {
                        searched.push_back(key);
                    }
                    this->bigQueries.emplace_back();
                    BigRangeQueryData<PointT> &query = this->bigQueries.back();
                    query.pointIdx = cell;
                    query.center = {vertex.x, vertex.y, vertex.z};
                    query.radius = fastabs(vertex - point);
                    query.originalPoint = {point.x, point.y, point.z};
                    query.askOnlyClose = askOnlyClose;
                }
            }
            // (marked after all the cells: a vertex of several cells is asked about for each)
            this->checkedVertices.insert(searched.begin(), searched.end());
        }
        for(const size_t &cell : this->smallPoints)
        {
            if(this->currentRadiuses[cell] <= 0)
            {
                MadVoro::Exception::MadVoroException eo("Radius for a certain point is <= 0 (in 'GhostSearch::PrepareRound')");
                eo.addEntry("Point Index", cell);
                eo.addEntry("Its radius", this->currentRadiuses[cell]);
                throw eo;
            }
            const PointT &point = mesh.Point(cell);
            this->smallQueries.emplace_back();
            SmallRangeQueryData<PointT> &query = this->smallQueries.back();
            query.pointIdx = cell;
            query.center = {point.x, point.y, point.z};
            query.radius = this->currentRadiuses[cell];
            query.maxPointsToGet = MAX_POINTS_TO_GET + 1;
        }
    }

    // the round's searches of periodic images (or not)
    void Select(bool images, std::vector<SmallRangeQueryData<PointT>> &smallSelected, std::vector<BigRangeQueryData<PointT>> &bigSelected) const
    {
        smallSelected.clear();
        bigSelected.clear();
        for(const SmallRangeQueryData<PointT> &query : this->smallQueries)
        {
            if((ImageCode(query.imageTranslation) != ZERO_IMAGE_CODE) == images)
            {
                smallSelected.push_back(query);
            }
        }
        for(const BigRangeQueryData<PointT> &query : this->bigQueries)
        {
            if((ImageCode(query.imageTranslation) != ZERO_IMAGE_CODE) == images)
            {
                bigSelected.push_back(query);
            }
        }
    }

    template <typename Query, typename Agent>
    void Answer(const std::vector<Query> &queries, Agent &agent, boost::container::flat_map<size_t, size_t> &results, std::vector<std::pair<std::size_t, PointT>> &found)
    {
        const std::vector<std::vector<size_t>> answers = agent.selfBatchAnswer(queries, this->ignorePoints);
        for(size_t i = 0; i < answers.size(); ++i)
        {
            for(const size_t &pointIdxInAll : answers[i])
            {
                found.emplace_back(pointIdxInAll, queries[i].imageTranslation);
            }
            results[queries[i].pointIdx] += answers[i].size();
        }
    }

    std::vector<Face3D<PointT>> box;
    std::vector<PointT> normals;
    std::array<bool, 3> periodic = {false, false, false};
    PointT ll, ur;
    boost::container::flat_set<size_t> smallPoints, largePoints;
    boost::container::flat_map<size_t, size_t> firstLargeIteration;
    std::vector<double> currentRadiuses;
    PeriodicIgnoreSets<PointT> ignorePoints;
    std::set<std::pair<size_t, size_t>> mirrored;
    boost::unordered_flat_set<std::array<double, 3>, boost::hash<std::array<double, 3>>> checkedVertices; // the vertices whose balls were searched on every rank
    std::vector<PointT> vertices;
    size_t iterations = 0;
    // the current round
    std::vector<SmallRangeQueryData<PointT>> smallQueries;
    std::vector<BigRangeQueryData<PointT>> bigQueries;
    std::vector<std::pair<size_t, size_t>> mirroredPoints;
    boost::container::flat_map<size_t, size_t> smallResults, bigResults;
    std::vector<size_t> remoteSmall, remoteBig;
    std::vector<RemoteAnswer> remoteAnswers;
    #ifdef MADVORO_WITH_MPI
        SentPointsContainer unusedSent; // (local searches send nothing)
    #endif // MADVORO_WITH_MPI
    BigRangeAgent<PointT> bigAgent;
    SmallRangeAgent<PointT> smallAgent;
    BigRangeAgent<PointT> imageBigAgent;
    SmallRangeAgent<PointT> imageSmallAgent;
};

} // namespace MadVoro

#endif // MADVORO_GHOSTS_GHOST_SEARCH_HPP
