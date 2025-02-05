/*! \file universal_error.hpp
  \brief A class for storing error and debug information
  \author Almog Yalinewich, Maor Mizrachi
 */
#ifndef UNIVERSAL_ERROR_HPP
#define UNIVERSAL_ERROR_HPP 1

#include <iostream>
#include <string>
#include <vector>
#include <any>
#include "utils/print/all.hpp"

using std::string;
using std::vector;
using std::pair;

namespace MadVoro
{
  namespace Exception
  {
    /*! \brief Container for error reports
    */
    class MadVoroException
    {
    private:
      struct PrintableAny
      {
        using PrintFunction = void(*)(std::ostream&, const std::any&);

        template<typename T>
        inline PrintableAny(const T &value)
        {
          this->value_ = value;
          this->printFunction_ = [](std::ostream &os, const std::any &value)
                                  {
                                    os << std::any_cast<T>(value);
                                  };
        }

        inline friend std::ostream &operator<<(std::ostream &os, const PrintableAny &p)
        {
          p.printFunction_(os, p.value_);
          return os;
        }

        std::any value_;
        PrintFunction printFunction_;
      };

    public:
      /*! \brief Class constructor
        \param err_msg Error message
      */
      explicit MadVoroException(const string &err_msg);

      /*! \brief Appends std::string to the error message
        \param msg Message to append
      */
      void Append2ErrorMessage(const std::string &msg);

      /*! \brief Returns the error message
        \return Error message
      */
      std::string const& getErrorMessage(void) const;

      ~MadVoroException(void);

      /*! \brief Copy constructor
        \param eo Source
      */
      MadVoroException(const MadVoroException& eo);

      template<typename T>
      inline void addEntry(const std::string &name, const T &value)
      {
        this->fields_.emplace_back(name, value);
      }

      /*! \brief Prints the contents of the error
      \param eo The error object
      */
      friend void reportError(MadVoroException const& eo, std::ostream& os);

    private:

      string err_msg_;

      std::vector<std::pair<std::string, PrintableAny>> fields_;
    };

    void reportError(MadVoroException const& eo, std::ostream& os = std::cout);
  }
}
#endif // UNIVERSAL_ERROR_HPP
