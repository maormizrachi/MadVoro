#include "hilbert/rectangular/HilbertRectangularConvertor3D.hpp"

using namespace MadVoro;
using namespace MadVoro::Geometry;

MadVoro::HilbertRectangularConvertor3D::HilbertRectangularConvertor3D(const Point3D &ll, const Point3D &ur, size_t order):
    HilbertConvertor3D(ll, ur, order)
{
    this->spaceBoundingBox = BoundingBox<Point3D>(ll, ur);
    this->changeOrder(order);
}

void MadVoro::HilbertRectangularConvertor3D::changeOrder(size_t order)
{
    this->order = order = std::min<size_t>(MAX_HILBERT_ORDER, order);
    coord_t realWidth = this->ur.x - this->ll.x;
    coord_t realHeight = this->ur.y - this->ll.y;
    coord_t realDepth = this->ur.z - this->ll.z;

    // calculate divisions number in x, y, z axises
    this->div.x = std::ceil(std::pow((realWidth * realWidth) / (realHeight * realDepth), 0.333333333) * std::pow(2, order));
    this->div.y = std::ceil(std::pow((realHeight * realHeight) / (realWidth * realDepth), 0.333333333) * std::pow(2, order));
    this->div.z = std::ceil(std::pow(8, order) / (this->div.x * this->div.y));
    
    this->total_points_num = this->div.x * this->div.y * this->div.z;
    this->step = Point3D(realWidth / div.x, realHeight / div.y, realDepth / div.z);
}

std::vector<MadVoro::HilbertRectangularConvertor3D::RecursionArguments> MadVoro::HilbertRectangularConvertor3D::getRecursionArguments(const RecursionArguments &args) const
{
    const DirectionPoint3D &startPoint = args.startPoint;
    const DirectionPoint3D &a = args.a;
    const DirectionPoint3D &b = args.b;
    const DirectionPoint3D &c = args.c;

    direction_t width = std::abs(a.x + a.y + a.z);
    direction_t height = std::abs(b.x + b.y + b.z);
    direction_t depth = std::abs(c.x + c.y + c.z);

    direction_t dax = SIGN(a.x), day = SIGN(a.y), daz = SIGN(a.z);
    direction_t dbx = SIGN(b.x), dby = SIGN(b.y), dbz = SIGN(b.z);
    direction_t dcx = SIGN(c.x), dcy = SIGN(c.y), dcz = SIGN(c.z);

    DirectionPoint3D a2 = {a.x >> 1, a.y >> 1, a.z >> 1}; /* (a.x//2, a.y//2, a.z//2) */
    DirectionPoint3D b2 = {b.x >> 1, b.y >> 1, b.z >> 1}; /* (b.x//2, b.y//2, b.z//2) */
    DirectionPoint3D c2 = {c.x >> 1, c.y >> 1, c.z >> 1}; /* (c.x//2, c.y//2, c.z//2) */

    direction_t width2 = std::abs(a2.x + a2.y + a2.z);
    direction_t height2 = std::abs(b2.x + b2.y + b2.z);
    direction_t depth2 = std::abs(c2.x + c2.y + c2.z);

    // prefer even steps
    if((width2 % 2) and (width > 2))
    {
        a2.x = a2.x + dax;
        a2.y = a2.y + day;
        a2.z = a2.z + daz;
    }

    if((height2 % 2) and (height > 2))
    {
        b2.x = b2.x + dbx;
        b2.y = b2.y + dby;
        b2.z = b2.z + dbz;
    }

    if((depth2 % 2) and (depth > 2))
    {
        c2.x = c2.x + dcx;
        c2.y = c2.y + dcy;
        c2.z = c2.z + dcz;
    }

    const direction_t &x = startPoint.x;
    const direction_t &y = startPoint.y;
    const direction_t &z = startPoint.z;

    std::vector<RecursionArguments> toReturn;

    if((2 * width > 3 * height) and (2 * width > 3 * depth))
    {
        toReturn.push_back({startPoint, a2, b, c});
        toReturn.push_back({{x + a2.x, y + a2.y, z + a2.z}, {a.x - a2.x, a.y - a2.y, a.z - a2.z}, b, c});
    }
    else if(3 * height > 4 * depth)
    {
        toReturn.push_back({startPoint, b2, c, a2});
        toReturn.push_back({{x + b2.x, y + b2.y, z + b2.z}, a, {b.x - b2.x, b.y - b2.y, b.z - b2.z}, c});
        toReturn.push_back({{x + (a.x - dax) + (b2.x - dbx), y + (a.y - day) + (b2.y - dby), z + (a.z - daz) + (b2.z - dbz)}, {-b2.x, -b2.y, -b2.z}, c, {-(a.x - a2.x), -(a.y - a2.y), -(a.z - a2.z)}});
    }
    else if(3 * depth > 4 * height)
    {
        toReturn.push_back({startPoint, c2, a2, b});
        toReturn.push_back({{x + c2.x, y + c2.y, z + c2.z}, a, b, {c.x - c2.x, c.y - c2.y, c.z - c2.z}});
        toReturn.push_back({{x + (a.x - dax) + (c2.x - dcx), y + (a.y - day) + (c2.y - dcy), z + (a.z - daz) + (c2.z - dcz)}, {-c2.x, -c2.y, -c2.z}, {-(a.x - a2.x), -(a.y - a2.y), -(a.z - a2.z)}, b});
    }
    else
    {
        toReturn.push_back({startPoint, b2, c2, a2});
        toReturn.push_back({{x + b2.x, y + b2.y, z + b2.z}, c, a2, {b.x - b2.x, b.y - b2.y, b.z - b2.z}});
        toReturn.push_back({{x + (b2.x - dbx) + (c.x - dcx), y + (b2.y - dby) + (c.y - dcy), z + (b2.z - dbz) + (c.z - dcz)}, a, {-b2.x, -b2.y, -b2.z}, {-(c.x - c2.x), -(c.y - c2.y), -(c.z - c2.z)}});
        toReturn.push_back({{x + (a.x - dax) + b2.x + (c.x - dcx), y + (a.y - day) + b2.y + (c.y - dcy), z + (a.z - daz) + b2.z + (c.z - dcz)}, {-c.x, -c.y, -c.z}, {-(a.x - a2.x), -(a.y - a2.y), -(a.z - a2.z)}, {b.x - b2.x, b.y - b2.y, b.z - b2.z}});
        toReturn.push_back({{x + (a.x - dax) + (b2.x - dbx), y + (a.y - day) + (b2.y - dby), z + (a.z - daz) + (b2.z - dbz)}, {-b2.x, -b2.y, -b2.z}, c2, {-(a.x - a2.x), -(a.y - a2.y), -(a.z - a2.z)}});
    }
    return toReturn;
}

Point3D MadVoro::HilbertRectangularConvertor3D::WidthHeightDepthToXYZ(direction_t width, direction_t height, direction_t depth) const
{
    coord_t x, y, z;
    x = this->ll[0] + width * this->step[0];
    y = this->ll[1] + height * this->step[1];
    z = this->ll[2] + depth * this->step[2];
    return Point3D(x, y, z);
}

bool MadVoro::HilbertRectangularConvertor3D::d2xyz_helper(const RecursionArguments &args, hilbert_index_t requested_d, hilbert_index_t &current_d, Point3D &result) const
{
    const DirectionPoint3D &startPoint = args.startPoint;
    const DirectionPoint3D &a = args.a;
    const DirectionPoint3D &b = args.b;
    const DirectionPoint3D &c = args.c;
    
    direction_t width = std::abs(a.x + a.y + a.z);
    direction_t height = std::abs(b.x + b.y + b.z);
    direction_t depth = std::abs(c.x + c.y + c.z);

    size_t num_points = width * height * depth;

    direction_t dax = SIGN(a.x), day = SIGN(a.y), daz = SIGN(a.z);
    direction_t dbx = SIGN(b.x), dby = SIGN(b.y), dbz = SIGN(b.z);
    direction_t dcx = SIGN(c.x), dcy = SIGN(c.y), dcz = SIGN(c.z);

    if(requested_d >= current_d + num_points)
    {
        // the rectangle we are iterating over currently is irrelevent
        current_d += num_points;
        return false;
    }

    if(requested_d < current_d)
    {
        throw MadVoro::Exception::MadVoroException("in MadVoro::HilbertRectangularConvertor3D::d2xyz_helper, should not reach here (algorithm failed)");
    }
    hilbert_index_t diff = requested_d - current_d;

    // base cases
    if(height == 1 and depth == 1)
    {
        result = this->WidthHeightDepthToXYZ(startPoint.x + diff * dax, startPoint.y + diff * day, startPoint.z + diff * daz);
        return true;
    }

    if(width == 1 and depth == 1)
    {
        result = this->WidthHeightDepthToXYZ(startPoint.x + diff * dbx, startPoint.y + diff * dby, startPoint.z + diff * dbz);
        return true;
    }

    if(width == 1 and height == 1)
    {
        result = this->WidthHeightDepthToXYZ(startPoint.x + diff * dcx, startPoint.y + diff * dcy, startPoint.z + diff * dcz);
        return true;
    }

    for(const RecursionArguments &nextArgs : this->getRecursionArguments(args))
    {
        if(this->d2xyz_helper(nextArgs, requested_d, current_d, result))
        {
            return true;
        }
    }

    return false;
}

bool MadVoro::HilbertRectangularConvertor3D::xyz2d_helper_base(const DirectionPoint3D &startPoint, size_t steps, const DirectionPoint3D &direction, const DirectionPoint3D &requested_point, hilbert_index_t &current_d) const
{
    direction_t x = startPoint.x, y = startPoint.y, z = startPoint.z;
    for(size_t i = 0; i < steps; i++)
    {
        if((requested_point.x == x) and (requested_point.y == y) and (requested_point.z == z))
        {
            return true;
        }
        x += direction.x;
        y += direction.y;
        z += direction.z;
        current_d++;
    }
    return false;
}

std::pair<typename MadVoro::HilbertRectangularConvertor3D::DirectionPoint3D, typename MadVoro::HilbertRectangularConvertor3D::DirectionPoint3D> MadVoro::HilbertRectangularConvertor3D::getBoundingBox(const RecursionArguments &args) const
{
    const DirectionPoint3D &startPoint = args.startPoint;
    const DirectionPoint3D &a = args.a;
    const DirectionPoint3D &b = args.b;
    const DirectionPoint3D &c = args.c;

    direction_t x_advancing = a.x + b.x + c.x;
    direction_t y_advancing = a.y + b.y + c.y;
    direction_t z_advancing = a.z + b.z + c.z;

    DirectionPoint3D boundary = {startPoint.x + x_advancing + ((x_advancing >= 0)? 1 : 0), startPoint.y + y_advancing + ((y_advancing >= 0)? 1 : 0), startPoint.z + z_advancing + ((z_advancing >= 0)? 1 : 0)};
    return {{std::min(startPoint.x, boundary.x) - 1, std::min(startPoint.y, boundary.y) - 1, std::min(startPoint.z, boundary.z) - 1}, 
            {std::max(startPoint.x, boundary.x) + 1, std::max(startPoint.y, boundary.y) + 1, std::max(startPoint.z, boundary.z) + 1}};    
}

bool MadVoro::HilbertRectangularConvertor3D::xyz2d_helper(const RecursionArguments &args, const DirectionPoint3D &requested_point, hilbert_index_t &current_d) const
{
    const DirectionPoint3D &startPoint = args.startPoint;
    const DirectionPoint3D &a = args.a;
    const DirectionPoint3D &b = args.b;
    const DirectionPoint3D &c = args.c;

    direction_t width = std::abs(a.x + a.y + a.z);
    direction_t height = std::abs(b.x + b.y + b.z);
    direction_t depth = std::abs(c.x + c.y + c.z);

    size_t num_points = width * height * depth;

    std::pair<DirectionPoint3D, DirectionPoint3D> bounding_box = this->getBoundingBox(args);
    if((requested_point.x < bounding_box.first.x) or (requested_point.x > bounding_box.second.x) or
        (requested_point.y < bounding_box.first.y) or (requested_point.y > bounding_box.second.y) or
        (requested_point.z < bounding_box.first.z) or (requested_point.z > bounding_box.second.z))
    {
        // doesn't have a chance to be here
        current_d += num_points;
        return false;
    }    

    direction_t dax = SIGN(a.x), day = SIGN(a.y), daz = SIGN(a.z);
    direction_t dbx = SIGN(b.x), dby = SIGN(b.y), dbz = SIGN(b.z);
    direction_t dcx = SIGN(c.x), dcy = SIGN(c.y), dcz = SIGN(c.z);

    // base cases
    if(height == 1 and depth == 1)
    {
        return this->xyz2d_helper_base(startPoint, width, {dax, day, daz}, requested_point, current_d);
    }

    if(width == 1 and depth == 1)
    {
        return this->xyz2d_helper_base(startPoint, height, {dbx, dby, dbz}, requested_point, current_d);
    }

    if(width == 1 and height == 1)
    {
        return this->xyz2d_helper_base(startPoint, depth, {dcx, dcy, dcz}, requested_point, current_d);
    }

    for(const RecursionArguments &nextArgs : this->getRecursionArguments(args))
    {
        if(this->xyz2d_helper(nextArgs, requested_point, current_d))
        {
            return true;
        }
    }
    return false;
}


Point3D MadVoro::HilbertRectangularConvertor3D::d2xyz(hilbert_index_t d) const
{
    Point3D result;
    hilbert_index_t current_d = 0;
    this->d2xyz_helper({{0, 0, 0}, {this->div.x, 0, 0}, {0, this->div.y, 0}, {0, 0, this->div.z}}, d, current_d, result);
    return result;
}

hilbert_index_t MadVoro::HilbertRectangularConvertor3D::xyz2d(coord_t x, coord_t y, coord_t z) const
{
    // convert (x,y, z) to the integer triple (width, height, width)
    direction_t width = std::floor((x - this->ll.x) / this->step.x);
    direction_t height = std::floor((y - this->ll.y) / this->step.y);
    direction_t depth = std::floor((z - this->ll.z) / this->step.z);

    // if(not this->spaceBoundingBox.contains(Point3D(x, y, z)))
    // {
    //     MadVoro::Exception::MadVoroException eo("MadVoro::HilbertRectangularConvertor3D::xyz2d: Given point is out of the bounding box of the space");
    //     eo.addEntry("Space bounding box", this->spaceBoundingBox);
    //     eo.addEntry("Point", Point3D(x, y, z));
    //     throw eo;
    // }

    if((width < 0) or (height < 0) or (depth < 0) or (width > this->div.x) or (height > this->div.y) or (depth > this->div.z))
    {
        MadVoro::Exception::MadVoroException eo("Should not reach here, overflow (in 3D xyz->d)");
        eo.addEntry("Width", width);
        eo.addEntry("this->div.x", this->div.x);
        eo.addEntry("Height", height);
        eo.addEntry("this->div.y", this->div.y);
        eo.addEntry("Depth", depth);
        eo.addEntry("this->div.z", this->div.z);
        eo.addEntry("x", x);
        eo.addEntry("y", y);
        eo.addEntry("z", z);
        eo.addEntry("step", this->step);
        eo.addEntry("ll", this->ll);
        throw eo;
    }

    hilbert_index_t result = 0;
    if(not this->xyz2d_helper({{0, 0, 0}, {this->div.x, 0, 0}, {0, this->div.y, 0}, {0, 0, this->div.z}}, {width, height, depth}, result))
    {
        MadVoro::Exception::MadVoroException eo("Should not reach here (in 3D xyz->d) (maybe the point is outside the box?)");
        eo.addEntry("x", x);
        eo.addEntry("y", y);
        eo.addEntry("z", z);
        throw eo;
    }
    return result;
}