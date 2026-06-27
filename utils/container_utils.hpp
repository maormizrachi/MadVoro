#ifndef MADVORO_CONTAINER_UTILS_HPP
#define MADVORO_CONTAINER_UTILS_HPP

#include <vector>
#include <algorithm>
#include <cstddef>

namespace MadVoro
{
namespace ContainerOps
{

template<typename Container>
inline void conditional_shrink(Container &v)
{
#ifdef RICH_AGGRESSIVE_SHRINK
    if (v.capacity() > v.size() * 1.25)
        v.shrink_to_fit();
#else
    (void)v;
#endif
}

template<typename Container>
inline void release_container_memory(Container &v)
{
    Container empty;
    v.swap(empty);
}

template <class T>
inline void unique_inplace(std::vector<T> &v)
{
    if (v.empty()) return;
    auto last = std::unique(v.begin(), v.end());
    v.erase(last, v.end());
}

template <class T>
inline void RemoveList_inplace(std::vector<T> &v, const std::vector<T> &list)
{
    auto new_end = std::remove_if(v.begin(), v.end(),
        [&list](const T &val){ return std::binary_search(list.begin(), list.end(), val); });
    v.erase(new_end, v.end());
}

template <class T>
inline void RemoveVal(std::vector<T> &vec, T val)
{
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (vec[i] == val)
        {
            vec.erase(vec.begin() + static_cast<long>(i));
            return;
        }
    }
}

template <class T>
inline std::vector<T> VectorValues(const std::vector<T> &v, const std::vector<std::size_t> &index)
{
    if (index.empty() || v.empty())
        return std::vector<T>();
    std::vector<T> result(index.size());
    for (std::size_t i = 0; i < index.size(); ++i)
        result[i] = v[index[i]];
    return result;
}

inline bool close2zero(double x)
{
    return (x < 1e-9 && x > -1e-9);
}

} // namespace ContainerOps
} // namespace MadVoro

#endif // MADVORO_CONTAINER_UTILS_HPP
