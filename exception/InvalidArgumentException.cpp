#include "InvalidArgumentException.hpp"

MadVoro::Exception::InvalidArgumentException::InvalidArgumentException() : MadVoro::Exception::MadVoroException(std::string("Illegal argument was given")) {}