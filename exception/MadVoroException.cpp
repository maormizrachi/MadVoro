#include <iostream>
#include <algorithm>
#include "MadVoroException.hpp"

using namespace std;

MadVoro::Exception::MadVoroException::MadVoroException(const string& err_msg):
  err_msg_(err_msg),
  fields_() {}

void MadVoro::Exception::MadVoroException::Append2ErrorMessage(string const& msg)
{
  err_msg_ += msg;
}

const string& MadVoro::Exception::MadVoroException::getErrorMessage(void) const
{
  return err_msg_;
}

MadVoro::Exception::MadVoroException::~MadVoroException(void) {}

MadVoro::Exception::MadVoroException::MadVoroException(const MadVoroException& eo):
  err_msg_(eo.getErrorMessage()),
  fields_(eo.fields_) {}


void MadVoro::Exception::reportError(MadVoroException const& eo, std::ostream& os)
{
  os.precision(14);
  os << eo.getErrorMessage() << std::endl;
  for_each(eo.fields_.begin(), eo.fields_.end(),
          [&os](const pair<string, MadVoroException::PrintableAny>& f) {os << f.first << " " << f.second << endl;});
}