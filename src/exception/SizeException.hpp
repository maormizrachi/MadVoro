#ifndef SIZE_EXCEPTION_HPP
#define SIZE_EXCEPTION_HPP

#include "MadVoroException.hpp"

namespace MadVoro
{
  namespace Exception
  {
    class SizeException : public MadVoroException
    {
    public:
      /*! \brief Class constructor
        \param err_msg Error message
      */
      explicit SizeException(void);

      explicit SizeException(const std::string &err_msg): MadVoroException(err_msg) {};
    };    
  }
}

#endif // SIZE_EXCEPTION_HPP