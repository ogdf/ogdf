
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_STRINGIZERS_H
#define IGLOO_STRINGIZERS_H

namespace snowhouse
{

  namespace detail
  {

    template<typename Container>
      struct SequentialContainerStringizer
      {
        static std::string
        ToString(const Container& cont)
        {
          std::ostringstream stm;
          typedef typename Container::const_iterator Iterator;

          stm << "[ ";
          for (Iterator it = cont.begin(); it != cont.end();)
          {
            stm << snowhouse::Stringize(*it);

            if (++it != cont.end())
            {
              stm << ", ";
            }
          }
          stm << " ]";
          return stm.str();
        }
      };
  }

  template<typename T>
    struct Stringizer<std::vector<T> > : detail::SequentialContainerStringizer<
        std::vector<T> >
    {
    };

  template<typename T>
    struct Stringizer<std::deque<T> > : detail::SequentialContainerStringizer<
        std::deque<T> >
    {
    };

  template<typename T>
    struct Stringizer<std::list<T> > : detail::SequentialContainerStringizer<
        std::list<T> >
    {
    };
}

#endif
