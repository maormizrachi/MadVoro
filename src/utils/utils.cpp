#include "utils.hpp"
#if defined(_MSC_VER)
/* Microsoft C/C++-compatible compiler */
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include <cfloat>
#include <cmath>

bool MadVoro::Utils::close2zero(const double x)
{
  constexpr double lowest_double = std::numeric_limits<double>::min();
  return std::abs(x)<lowest_double;
}

double MadVoro::Utils::min(vector<double> const& v)
{
  double res = v[0];
  for(size_t i=1;i<v.size();++i)
    res = std::min(res,v[size_t(i)]);
  return res;
}

double MadVoro::Utils::max(vector<double> const& v)
{
  double res = v[0];
  for(size_t i=1;i<v.size();++i)
    res = std::max(res,v[size_t(i)]);
  return res;
}

double MadVoro::Utils::fastsqrt(double x)
{
  if (x<static_cast<double>(FLT_MIN) || x > static_cast<double>(FLT_MAX))
		return std::sqrt(x);
	double res = static_cast<double>(_mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(static_cast<float>(x)))));
	return x*res*(1.5 - 0.5*res*res*x);
}

double MadVoro::Utils::Interpolate2DTable(double const x, double const y, std::vector<double> const& x_vec, std::vector<double> const& y_vec, std::vector<std::vector<double>> const& data,
                                        double const x_vec_high_slope, size_t const slope_length)
{
    // Check if the size of x_vec or y_vec is less than or equal to slope_length
    if(x_vec.size() <= slope_length || y_vec.size() <= slope_length)
        throw MadVoro::Exception::MadVoroException("Interpolate2DTable: x_vec.size() <= slope_length || y__vec.size() <= slope_length");

    // Handle cases where x or y is less than the minimum value in x_vec or y_vec
    if(x < x_vec[0])
    {
        if(y < y_vec[0])
        {
            // Calculate the slope for y and x
            double const y_slope = (data[0][slope_length - 1] -data[0][0]) / (y_vec[slope_length - 1] - y_vec[0]);
            double const x_vec_slope = (data[slope_length - 1][0] -data[0][0]) / (x_vec[slope_length - 1] - x_vec[0]);

            // Perform interpolation using the calculated slopes
            return std::exp(data[0][0] + y_slope * (y - y_vec[0]) + x_vec_slope * (x - x_vec[0]));
        }
        else
        {
            // Calculate the interpolated value at x_vec[0] * 1.00001 and y
            double const data_x0 = BiLinearInterpolation(x_vec, y_vec, data, x_vec[0] * 1.00001, y);

            // Calculate the slope for x
            double const x_vec_slope = (BiLinearInterpolation(x_vec, y_vec, data, x_vec[slope_length - 1], y) - data_x0) / (x_vec[slope_length - 1] - x_vec[0]);

            // Perform interpolation using the calculated slope
            return std::exp(BiLinearInterpolation(x_vec, y_vec, data, x_vec[0] * 1.00001, y) + x_vec_slope * (x - x_vec[0]));
        }
    }

    // Handle cases where x is greater than the maximum value in x_vec
    if(x > x_vec.back())
    {
        if(y < y_vec[0])
        {
            // Calculate the slope for y
            double const y_slope = (data[x_vec.size() - 1][slope_length - 1] -data[x_vec.size() - 1][0]) / (y_vec[slope_length - 1] - y_vec[0]);

            // Perform interpolation using the calculated slope and x_vec_high_slope
            return std::exp(data[x_vec.size() - 1][0] + y_slope * (y - y_vec[0]) + x_vec_high_slope * (x - x_vec.back()));
        }
        else
            // Perform interpolation using the calculated slope and x_vec_high_slope
            return std::exp(BiLinearInterpolation(x_vec, y_vec, data, x_vec.back() * 0.99999, y) + x_vec_high_slope * (x - x_vec.back()));
    }

    // Handle cases where y is less than the minimum value in y_vec
    if(y < y_vec[0])
    {
        // Calculate the interpolated value at x and y_vec[0] * 0.9999
        double const data_y0 = BiLinearInterpolation(x_vec, y_vec, data, x, y_vec[0] * 0.9999);

        // Calculate the slope for y
        double const y_slope =(BiLinearInterpolation(x_vec, y_vec, data, x, y_vec[slope_length - 1]) - data_y0) / (y_vec[slope_length - 1] - y_vec[0]);

        // Perform interpolation using the calculated slope
        return std::exp(BiLinearInterpolation(x_vec, y_vec, data, x, y_vec[0] * 0.9999) + y_slope * (y - y_vec[0]));
    }

    // Perform interpolation using BiLinearInterpolation
    return std::exp(BiLinearInterpolation(x_vec, y_vec, data, x, y));
}