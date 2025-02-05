#ifndef READ_POINTS_HPP
#define READ_POINTS_HPP

#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <madvoro/Vector3D.hpp>

inline std::vector<MadVoro::Vector3D> readPoints(const std::string &filename)
{
    std::vector<MadVoro::Vector3D> points;
    std::ifstream file(filename, std::ios::in);

    if(!file.good())
    {
        std::cerr << "An error in the given file ('" << filename << "')." << std::endl;
        exit(EXIT_FAILURE);
    }
    int lineCounter = 1;
    std::string line;

    while(std::getline(file, line))
    {
        if(*(line.begin()) != '(' or *(line.end() - 1) != ')')
        {
            std::cerr << "Error in line " << lineCounter << ": " << line << std::endl;
            file.close();
            exit(EXIT_FAILURE);
        }
        auto it = std::find(line.begin(), line.end(), ',');
        if(it == line.end())
        {
            std::cerr << "No comma separation was given between x and y in line " << lineCounter << ": " << line << std::endl;
            file.close();
            exit(EXIT_FAILURE);
        }
        double _x = std::stod(std::string(line.begin() + 1, it));
        auto it2 = std::find(it + 1, line.end(), ',');
        if(it2 == line.end())
        {
            std::cerr << "No comma separation was given between y and z in line " << lineCounter << ": " << line << std::endl;
            file.close();
            exit(EXIT_FAILURE);
        }
        double _y = std::stod(std::string(it + 1, it2));
        double _z = std::stod(std::string(it2 + 1, line.end() - 1));

        points.emplace_back(_x, _y, _z);
    }

    file.close();
    return points;
}

#endif // READ_POINTS_HPP