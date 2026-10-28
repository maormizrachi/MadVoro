// Searches of a thread subdomain for the points of the other subdomains
#ifndef MADVORO_SUBDOMAINS_SUBDOMAIN_SEARCH_HPP
#define MADVORO_SUBDOMAINS_SUBDOMAIN_SEARCH_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include <spatial_ds/OctTree/OctTree.hpp>
#include "../range/finders/RangeFinder.hpp"
#include "../range/finders/utils/IndexedVector.hpp"

namespace MadVoro
{

// Coarse grid of the subdomain range (smallest, largest) and of other ranks' points in each cell
template <typename PointT>
class SubdomainGrid
{
public:
    using subdomain_t = std::uint16_t;
    static constexpr std::size_t MAX_SUBDOMAINS = 0xFFFE;

    struct Reach
    {
        bool otherSubdomains = false;
        bool otherRanks = false;
    };

    SubdomainGrid(const std::vector<PointT> &points, const std::vector<subdomain_t> &subdomainOf)
    {
        // about two points per cell
        const double cells = std::max(1.0, static_cast<double>(points.size()) / 2.0);
        this->side_ = std::clamp(static_cast<int>(std::lround(std::cbrt(cells))), 1, 256);

        // the points' bounding box, with a margin for the points other ranks send
        std::array<double, 3> upper{};
        if(not points.empty())
        {
            this->lower_ = {points[0].x, points[0].y, points[0].z};
            upper = this->lower_;
        }
        for(const PointT &point : points)
        {
            const std::array<double, 3> p = {point.x, point.y, point.z};
            for(int k = 0; k < 3; ++k)
            {
                this->lower_[k] = std::min(this->lower_[k], p[k]);
                upper[k] = std::max(upper[k], p[k]);
            }
        }
        for(int k = 0; k < 3; ++k)
        {
            const double margin = 0.1 * (upper[k] - this->lower_[k]);
            this->lower_[k] -= margin;
            upper[k] += margin;
            const double width = upper[k] - this->lower_[k];
            this->inverseWidth_[k] = (width > 0) ? this->side_ / width : 0.0;
        }

        const std::size_t cellsNum = static_cast<std::size_t>(this->side_) * this->side_ * this->side_;
        this->smallest_.assign(cellsNum, EMPTY_SMALLEST);
        this->largest_.assign(cellsNum, 0);
        this->otherRanks_.assign(cellsNum, 0);
        for(std::size_t i = 0; i < points.size(); ++i)
        {
            const std::size_t cell = this->CellOf(points[i]);
            this->smallest_[cell] = std::min(this->smallest_[cell], subdomainOf[i]);
            this->largest_[cell] = std::max(this->largest_[cell], subdomainOf[i]);
        }
    }

    template<typename Iterator>
    void AddOtherRankPoints(Iterator first, Iterator last)
    {
        for(Iterator it = first; it != last; ++it)
        {
            const std::array<double, 3> p = {it->x, it->y, it->z};
            bool inside = true;
            for(int k = 0; k < 3; ++k)
            {
                const double cell = std::floor((p[k] - this->lower_[k]) * this->inverseWidth_[k]);
                inside = inside and cell >= 0 and cell <= this->side_ - 1;
            }
            if(inside)
            {
                this->otherRanks_[this->CellOf(*it)] = 1;
            }
            else
            {
                this->othersOutside_ = true;
            }
        }
    }

    Reach MayHold(const PointT &center, double radius, subdomain_t subdomain) const
    {
        Reach reach;
        // the searches take points up to slightly beyond the radius
        const double extent = radius * (1 + 1e-9) + 1e-11;
        const std::array<double, 3> c = {center.x, center.y, center.z};
        std::array<int, 3> from{}, to{};
        bool leavesGrid = false;
        for(int k = 0; k < 3; ++k)
        {
            const double low = std::floor((c[k] - extent - this->lower_[k]) * this->inverseWidth_[k]);
            const double high = std::floor((c[k] + extent - this->lower_[k]) * this->inverseWidth_[k]);
            leavesGrid = leavesGrid or low < 0 or high > this->side_ - 1;
            if(high < 0 or low > this->side_ - 1)
            {
                // the box is outside the grid
                reach.otherRanks = this->othersOutside_;
                return reach;
            }
            from[k] = static_cast<int>(std::max(low, 0.0));
            to[k] = static_cast<int>(std::min(high, static_cast<double>(this->side_ - 1)));
        }
        reach.otherRanks = leavesGrid and this->othersOutside_;
        for(int i = from[0]; i <= to[0]; ++i)
        {
            for(int j = from[1]; j <= to[1]; ++j)
            {
                for(int k = from[2]; k <= to[2]; ++k)
                {
                    const std::size_t cell = this->CellIndex(i, j, k);
                    reach.otherSubdomains = reach.otherSubdomains or (this->smallest_[cell] <= this->largest_[cell] and (this->smallest_[cell] != subdomain or this->largest_[cell] != subdomain));
                    reach.otherRanks = reach.otherRanks or (this->otherRanks_[cell] != 0);
                    if(reach.otherSubdomains and reach.otherRanks)
                    {
                        return reach;
                    }
                }
            }
        }
        return reach;
    }

private:
    static constexpr subdomain_t EMPTY_SMALLEST = 0xFFFF;

    template<typename P>
    std::size_t CellOf(const P &point) const
    {
        return this->CellIndex(this->Cell(point.x, 0), this->Cell(point.y, 1), this->Cell(point.z, 2));
    }

    int Cell(double coordinate, int axis) const
    {
        const int cell = static_cast<int>(std::floor((coordinate - this->lower_[axis]) * this->inverseWidth_[axis]));
        return std::clamp(cell, 0, this->side_ - 1);
    }

    std::size_t CellIndex(int i, int j, int k) const
    {
        return (static_cast<std::size_t>(i) * this->side_ + j) * this->side_ + k;
    }

    int side_ = 1;
    std::array<double, 3> lower_{};
    std::array<double, 3> inverseWidth_{};
    std::vector<subdomain_t> smallest_, largest_;
    std::vector<std::uint8_t> otherRanks_;
    bool othersOutside_ = false;
};

// Points other ranks sent, in one tree; the i-th has index FirstIndex() + i. Added only between searches.
template <typename PointT>
class RankGhostPoints
{
public:
    using coord_type = typename PointT::coord_type;
    static constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();

    RankGhostPoints(const PointT &ll, const PointT &ur, std::size_t firstIndex) : tree_(IVec(ll - (ur - ll) * 0.00001, NONE), IVec(ur + (ur - ll) * 0.00001, NONE)), firstIndex_(firstIndex)
    {
    }

    void Add(const PointT &point)
    {
        this->tree_.insert(IVec(point.x, point.y, point.z, this->firstIndex_ + this->points_.size()));
        this->points_.push_back(point);
    }

    const std::vector<PointT> &Points() const
    {
        return this->points_;
    }

    std::size_t FirstIndex() const
    {
        return this->firstIndex_;
    }

    void Range(const PointT &center, coord_type radius, std::size_t N, const std::unordered_set<std::size_t> &ignore, std::vector<std::size_t> &result) const
    {
        if(this->points_.empty() or N == 0)
        {
            return;
        }
        constexpr coord_type eps = static_cast<coord_type>(1e-12);
        for(const IVec &vec : this->tree_.range(Sphere<IVec>(IVec(center.x, center.y, center.z, NONE), radius + eps), N, [&ignore](const IVec &v) { return ignore.find(v.getIndex()) == ignore.cend(); }))
        {
            result.push_back(vec.index);
        }
    }

    // The closest non-ignored point in the sphere and its squared distance (NONE if there is none)
    std::pair<std::size_t, coord_type> ClosestInSphere(const PointT &center, coord_type radius, const PointT &point, const std::unordered_set<std::size_t> &ignore) const
    {
        if(this->points_.empty())
        {
            return {NONE, std::numeric_limits<coord_type>::max()};
        }
        const std::pair<IVec, coord_type> closest = this->tree_.getClosestPointInSphere(Sphere<PointT>(center, radius), point, [&ignore](const IVec &v) { return ignore.find(v.getIndex()) == ignore.cend(); });
        if(closest.second == std::numeric_limits<coord_type>::max())
        {
            return {NONE, closest.second};
        }
        return {closest.first.index, closest.second};
    }

private:
    using IVec = IndexedVector<PointT>;

    OctTree<IVec> tree_;
    std::vector<PointT> points_;
    std::size_t firstIndex_;
};

// Searches other subdomains' points, then other ranks' points, skipping what the grid rules out. Only valid when all of the subdomain's own points are ignored.
template <typename PointT>
class SubdomainRangeFinder : public RangeFinder<PointT>
{
public:
    using coord_type = typename RangeFinder<PointT>::coord_type;
    using subdomain_t = typename SubdomainGrid<PointT>::subdomain_t;

    SubdomainRangeFinder(std::shared_ptr<const RangeFinder<PointT>> points, std::shared_ptr<const RankGhostPoints<PointT>> ghosts, std::shared_ptr<const SubdomainGrid<PointT>> grid, subdomain_t subdomain) : points_(std::move(points)), ghosts_(std::move(ghosts)), grid_(std::move(grid)), subdomain_(subdomain)
    {
    }

    std::vector<size_t> range(const PointT &center, coord_type radius, size_t N, const typename RangeFinder<PointT>::template _set<size_t> &ignore) const override
    {
        const typename SubdomainGrid<PointT>::Reach reach = this->grid_->MayHold(center, radius, this->subdomain_);
        std::vector<size_t> result;
        if(reach.otherSubdomains)
        {
            result = this->points_->range(center, radius, N, ignore);
        }
        if(reach.otherRanks)
        {
            this->ghosts_->Range(center, radius, N - std::min(N, result.size()), ignore, result);
        }
        return result;
    }

    std::vector<size_t> closestPointInSphere(const PointT &center, coord_type radius, const PointT &point, const typename RangeFinder<PointT>::template _set<size_t> &ignore) const override
    {
        const typename SubdomainGrid<PointT>::Reach reach = this->grid_->MayHold(center, radius, this->subdomain_);
        std::vector<size_t> result;
        if(reach.otherSubdomains)
        {
            result = this->points_->closestPointInSphere(center, radius, point, ignore);
        }
        if(not reach.otherRanks)
        {
            return result;
        }
        const std::pair<std::size_t, coord_type> ghost = this->ghosts_->ClosestInSphere(center, radius, point, ignore);
        if(ghost.first != RankGhostPoints<PointT>::NONE)
        {
            const bool closer = result.empty() or ghost.second < SquaredDistance(this->points_->getPoint(result.front()), point);
            if(closer)
            {
                result.assign(1, ghost.first);
            }
        }
        return result;
    }

    const PointT &getPoint(size_t index) const override
    {
        const std::size_t first = this->ghosts_->FirstIndex();
        return (index < first) ? this->points_->getPoint(index) : this->ghosts_->Points()[index - first];
    }

    size_t size() const override
    {
        return this->points_->size() + this->ghosts_->Points().size();
    }

private:
    static coord_type SquaredDistance(const PointT &a, const PointT &b)
    {
        const coord_type dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
        return dx * dx + dy * dy + dz * dz;
    }

    std::shared_ptr<const RangeFinder<PointT>> points_;
    std::shared_ptr<const RankGhostPoints<PointT>> ghosts_;
    std::shared_ptr<const SubdomainGrid<PointT>> grid_;
    subdomain_t subdomain_;
};

} // namespace MadVoro

#endif // MADVORO_SUBDOMAINS_SUBDOMAIN_SEARCH_HPP
