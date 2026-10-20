/** \file SubdomainSearch.hpp
   \brief Searches of a thread subdomain for the points of the other subdomains
   \author Maor Mizrachi
*/
#ifndef MADVORO_SUBDOMAINS_SUBDOMAIN_SEARCH_HPP
#define MADVORO_SUBDOMAINS_SUBDOMAIN_SEARCH_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "../range/finders/RangeFinder.hpp"

namespace MadVoro
{

/*! \brief Which thread subdomains have points in each cell of a coarse grid over a
    rank's points: the smallest and the largest subdomain of the cell's points. */
template <typename PointT>
class SubdomainGrid
{
public:
    using subdomain_t = std::uint16_t;

    /*! \param points The rank's points
        \param subdomainOf The subdomain of each point */
    SubdomainGrid(const std::vector<PointT> &points, const std::vector<subdomain_t> &subdomainOf)
    {
        // about two points per cell
        const double cells = std::max(1.0, static_cast<double>(points.size()) / 2.0);
        this->side_ = std::clamp(static_cast<int>(std::lround(std::cbrt(cells))), 1, 256);

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
            const double width = upper[k] - this->lower_[k];
            this->inverseWidth_[k] = (width > 0) ? this->side_ / width : 0.0;
        }

        const std::size_t cellsNum = static_cast<std::size_t>(this->side_) * this->side_ * this->side_;
        this->smallest_.assign(cellsNum, EMPTY_SMALLEST);
        this->largest_.assign(cellsNum, 0);
        for(std::size_t i = 0; i < points.size(); ++i)
        {
            const std::array<double, 3> p = {points[i].x, points[i].y, points[i].z};
            const std::size_t cell = this->CellIndex(this->Cell(p[0], 0), this->Cell(p[1], 1), this->Cell(p[2], 2));
            this->smallest_[cell] = std::min(this->smallest_[cell], subdomainOf[i]);
            this->largest_[cell] = std::max(this->largest_[cell], subdomainOf[i]);
        }
    }

    /*! \brief Whether points of a subdomain other than `subdomain` may be in the sphere:
        false only when every cell the sphere's bounding box meets has no such point. */
    bool MayHaveOthers(const PointT &center, double radius, subdomain_t subdomain) const
    {
        // the searches take points up to slightly beyond the radius
        const double reach = radius * (1 + 1e-9) + 1e-11;
        const std::array<double, 3> c = {center.x, center.y, center.z};
        std::array<int, 3> from{}, to{};
        for(int k = 0; k < 3; ++k)
        {
            const double low = std::floor((c[k] - reach - this->lower_[k]) * this->inverseWidth_[k]);
            const double high = std::floor((c[k] + reach - this->lower_[k]) * this->inverseWidth_[k]);
            if(high < 0 or low > this->side_ - 1)
            {
                // the box is outside the points' bounding box
                return false;
            }
            from[k] = static_cast<int>(std::max(low, 0.0));
            to[k] = static_cast<int>(std::min(high, static_cast<double>(this->side_ - 1)));
        }
        for(int i = from[0]; i <= to[0]; ++i)
        {
            for(int j = from[1]; j <= to[1]; ++j)
            {
                for(int k = from[2]; k <= to[2]; ++k)
                {
                    const std::size_t cell = this->CellIndex(i, j, k);
                    if(this->smallest_[cell] <= this->largest_[cell] and
                       (this->smallest_[cell] != subdomain or this->largest_[cell] != subdomain))
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

private:
    static constexpr subdomain_t EMPTY_SMALLEST = 0xFFFF;

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
};

/*! \brief A thread subdomain's searches for the points of the other subdomains, in the
    rank's points. A search whose sphere can hold no point of another subdomain (see
    SubdomainGrid) returns nothing without searching; any other search is the rank's.
    Only for searches that ignore all of the subdomain's own points, as the subdomain's
    searches for ghosts do, so that what they find is all of other subdomains. */
template <typename PointT>
class SubdomainRangeFinder : public RangeFinder<PointT>
{
public:
    using coord_type = typename RangeFinder<PointT>::coord_type;
    using subdomain_t = typename SubdomainGrid<PointT>::subdomain_t;

    SubdomainRangeFinder(std::shared_ptr<const RangeFinder<PointT>> points, std::shared_ptr<const SubdomainGrid<PointT>> grid,
                         subdomain_t subdomain)
        : points_(std::move(points)), grid_(std::move(grid)), subdomain_(subdomain)
    {}

    std::vector<size_t> range(const PointT &center, coord_type radius, size_t N,
                              const typename RangeFinder<PointT>::template _set<size_t> &ignore) const override
    {
        if(not this->grid_->MayHaveOthers(center, radius, this->subdomain_))
            return std::vector<size_t>();
        return this->points_->range(center, radius, N, ignore);
    }

    std::vector<size_t> closestPointInSphere(const PointT &center, coord_type radius, const PointT &point,
                                             const typename RangeFinder<PointT>::template _set<size_t> &ignore) const override
    {
        if(not this->grid_->MayHaveOthers(center, radius, this->subdomain_))
            return std::vector<size_t>();
        return this->points_->closestPointInSphere(center, radius, point, ignore);
    }

    const PointT &getPoint(size_t index) const override
    {
        return this->points_->getPoint(index);
    }

    size_t size() const override
    {
        return this->points_->size();
    }

private:
    std::shared_ptr<const RangeFinder<PointT>> points_;
    std::shared_ptr<const SubdomainGrid<PointT>> grid_;
    subdomain_t subdomain_;
};

} // namespace MadVoro

#endif // MADVORO_SUBDOMAINS_SUBDOMAIN_SEARCH_HPP
