
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_ASSERTMACRO_H
#define IGLOO_ASSERTMACRO_H

#include "assert.h"

#define SNOWHOUSE_ASSERT_THAT(p1,p2,FAILURE_HANDLER)\
  ::snowhouse::ConfigurableAssert<FAILURE_HANDLER>::That((p1), (p2), __FILE__, __LINE__);\

#ifndef SNOWHOUSE_NO_MACROS

#define AssertThat(p1,p2)\
  SNOWHOUSE_ASSERT_THAT((p1), (p2), ::snowhouse::DefaultFailureHandler);\

#endif // SNOWHOUSE_NO_MACROS

#endif	// IGLOO_ASSERTMACRO_H
