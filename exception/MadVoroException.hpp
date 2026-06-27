/*! \file MadVoroException.hpp
  \brief A class for storing error and debug information
  \author Almog Yalinewich, Maor Mizrachi
 */
#ifndef MADVORO_EXCEPTION_HPP
#define MADVORO_EXCEPTION_HPP 1

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <any>

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
        std::string str_;

        template<typename T>
        static void toStream(std::ostream &os, const T &v) { os << v; }

        template<typename T>
        static void toStream(std::ostream &os, const std::vector<T> &v)
        {
          os << "[";
          for (size_t i = 0; i < v.size(); ++i) { if (i) os << ", "; toStream(os, v[i]); }
          os << "]";
        }

        template<typename A, typename B>
        static void toStream(std::ostream &os, const std::pair<A,B> &p)
        {
          os << "("; toStream(os, p.first); os << ", "; toStream(os, p.second); os << ")";
        }

        template<typename T>
        inline PrintableAny(const T &value)
        {
          std::ostringstream oss;
          toStream(oss, value);
          str_ = oss.str();
        }

        inline friend std::ostream &operator<<(std::ostream &os, const PrintableAny &p)
        {
          return os << p.str_;
        }
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
#endif // MADVORO_EXCEPTION_HPP
