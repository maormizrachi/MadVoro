#ifndef SERIALIZABLE_HPP
#define SERIALIZABLE_HPP

namespace MadVoro
{
    namespace MPI
    {
        class Serializer;
        
        class Serializable
        {
        public:
            virtual ~Serializable(void) = default;
        
            virtual size_t dump(Serializer *serializer) const = 0;
        
            virtual size_t load(const Serializer *serializer, std::size_t byteOffset) = 0;
        };
        
        template<typename T>
        using is_serializable = std::is_convertible<T*, Serializable*>;
    }
}

#endif // SERIALIZABLE_HPP