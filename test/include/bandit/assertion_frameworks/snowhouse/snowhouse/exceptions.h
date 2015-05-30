
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_EXCEPTIONS_H
#define IGLOO_EXCEPTIONS_H

#include "assert.h"

namespace snowhouse {

  template <typename ExceptionType>
  class ExceptionStorage
  {
  public:
    static std::unique_ptr<ExceptionType>& last_exception()
    {
      static std::unique_ptr<ExceptionType> last;
      return last;
    }

    void compiler_thinks_i_am_unused() {}

    ~ExceptionStorage()
    {
      last_exception().reset(NULL);
    }
  };

  template <typename ExceptionType>
  inline ExceptionType& LastException()
  {
    if(ExceptionStorage<ExceptionType>::last_exception().get() == NULL)
    {
      Assert::Failure("No exception was stored");
    }

    return *(ExceptionStorage<ExceptionType>::last_exception().get());
  }
}

#define IGLOO_CONCAT2(a, b) a##b
#define IGLOO_CONCAT(a, b) IGLOO_CONCAT2(a, b)

#define AssertThrows(EXCEPTION_TYPE, METHOD) \
ExceptionStorage<EXCEPTION_TYPE> IGLOO_CONCAT(IGLOO_storage_, __LINE__); IGLOO_CONCAT(IGLOO_storage_, __LINE__).compiler_thinks_i_am_unused(); \
{ \
  bool wrong_exception = false; \
  bool no_exception = false; \
  try \
  { \
    METHOD; \
    no_exception = true; \
  } \
  catch (const EXCEPTION_TYPE& e) \
  { \
    ExceptionStorage<EXCEPTION_TYPE>::last_exception() = std::unique_ptr<EXCEPTION_TYPE>(new EXCEPTION_TYPE(e)); \
  } \
  catch(...) \
  { \
    wrong_exception = true; \
  } \
  if(no_exception) \
  { \
    std::ostringstream stm; \
    stm << "Expected " << #EXCEPTION_TYPE << ". No exception was thrown."; \
    Assert::Failure(stm.str()); \
  } \
  if(wrong_exception) \
  { \
    std::ostringstream stm; \
    stm << "Expected " << #EXCEPTION_TYPE << ". Wrong exception was thrown."; \
    Assert::Failure(stm.str()); \
  } \
}

#endif


