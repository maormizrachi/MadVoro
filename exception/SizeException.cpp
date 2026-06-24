#include "SizeException.hpp"

MadVoro::Exception::SizeException::SizeException() : MadVoro::Exception::MadVoroException(std::string("Wrong/Illegal size")) {}