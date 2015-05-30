
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_ASSERTMACRO_H
#define IGLOO_ASSERTMACRO_H

#include "assert.h"

#define AssertThat(p1,p2)\
  Assert::That((p1), (p2), __FILE__, __LINE__);\

#endif	// IGLOO_ASSERTMACRO_H
