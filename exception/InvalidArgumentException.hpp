#ifndef ILLEGAL_ARGUMENT_EXCEPTION_HPP
#define ILLEGAL_ARGUMENT_EXCEPTION_HPP

#include "MadVoroException.hpp"

namespace MadVoro
{
  namespace Exception
  {
    class InvalidArgumentException : public MadVoroException
    {
    public:
      /*! \brief Class constructor
        \param err_msg Error message
      */
      explicit InvalidArgumentException(void);

      explicit InvalidArgumentException(const std::string &err_msg): MadVoroException(err_msg) {};
    };
  }
}

#endif // SIZE_EXCEPTION_HPP