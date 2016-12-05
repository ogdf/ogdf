
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_STRINGIZE_H
#define IGLOO_STRINGIZE_H

#include <cstddef>

namespace snowhouse {
  namespace detail {

    // This type soaks up any implicit conversions and makes the following operator<<
    // less preferred than any other such operator found via ADL.
    struct any
    {
      // Conversion constructor for any type.
      template <class T>
      any(T const&);
    };

    // A tag type returned by operator<< for the any struct in this namespace
    // when T does not support <<.
    struct tag {};

    // Fallback operator<< for types T that don't support <<.
    tag operator<<(std::ostream&, any const&);

    // Two overloads to distinguish whether T supports a certain operator expression.
    // The first overload returns a reference to a two-element character array and is chosen if
    // T does not support the expression, such as <<, whereas the second overload returns a char
    // directly and is chosen if T supports the expression. So using sizeof(check(<expression>))
    // returns 2 for the first overload and 1 for the second overload.
    typedef char yes;
    typedef char (&no)[2];

    no check(tag);

    template <class T>
    yes check(T const&);

    template <class T>
    struct is_output_streamable
    {
      static const T& x;
      static const bool value = sizeof(check(std::cout << x)) == sizeof(yes);
    };

    template<typename T, bool type_is_streamable>
    struct DefaultStringizer
    {
      static std::string ToString(const T& value)
      {
        std::ostringstream buf;
        buf << value;
        return buf.str();
      }
    };

    template<typename T>
    struct DefaultStringizer<T, false>
    {
      static std::string ToString(const T&)
      {
        return "[unsupported type]";
      }
    };
  }

  template<typename T>
  struct Stringizer;

  template<typename T>
  std::string Stringize(const T& value)
  {
    return Stringizer<T>::ToString(value);
  }

  // NOTE: Specialize snowhouse::Stringizer to customize assertion messages
  template<typename T>
  struct Stringizer
  {
    static std::string ToString(const T& value)
    {
      return detail::DefaultStringizer< T, detail::is_output_streamable<T>::value >::ToString(value);
    }
  };

  // We need this because nullptr_t has ambiguous overloads of operator<< in the standard library.
  template<>
  struct Stringizer<std::nullptr_t>
  {
    static std::string ToString(std::nullptr_t)
    {
      return "nullptr";
    }
  };
}

#endif
