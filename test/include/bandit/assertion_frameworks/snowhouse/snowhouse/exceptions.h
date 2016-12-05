
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
    static void last_exception(ExceptionType*** e, bool clear=false)
    {
      static ExceptionType* last = NULL;
      if(clear && last)
      {
        delete last;
        return;
      }

      *e = &last;
      silly_warning_about_unused_arg(e);
    }

    static ExceptionType*** silly_warning_about_unused_arg(ExceptionType*** e)
    {
      return e;
    }

    static void store(const ExceptionType& e)
    {
      ExceptionType** last = NULL;
      last_exception(&last);
      if(*last)
      {
        delete *last;
        *last = NULL;
      }

      *last = new ExceptionType(e);
    }

    void compiler_thinks_i_am_unused() {}

    ~ExceptionStorage()
    {
      ExceptionType** e = NULL;
      last_exception(&e);
      if(*e)
      {
        delete *e;
        *e = NULL;
      }
    }
  };

  template <typename ExceptionType>
  inline ExceptionType& LastException()
  {
    ExceptionType** e = NULL;
    ExceptionStorage<ExceptionType>::last_exception(&e);
    if(*e == NULL)
    {
      Assert::Failure("No exception was stored");
    }

    return **e;
  }
}

#define IGLOO_CONCAT2(a, b) a##b
#define IGLOO_CONCAT(a, b) IGLOO_CONCAT2(a, b)

#define SNOWHOUSE_ASSERT_THROWS(EXCEPTION_TYPE, METHOD, FAILURE_HANDLER_TYPE) \
::snowhouse::ExceptionStorage<EXCEPTION_TYPE> IGLOO_CONCAT(IGLOO_storage_, __LINE__); IGLOO_CONCAT(IGLOO_storage_, __LINE__).compiler_thinks_i_am_unused(); \
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
    ::snowhouse::ExceptionStorage<EXCEPTION_TYPE>::store(e); \
  } \
  catch(...) \
  { \
    wrong_exception = true; \
  } \
  if(no_exception) \
  { \
    std::ostringstream stm; \
    stm << "Expected " << #EXCEPTION_TYPE << ". No exception was thrown."; \
    ::snowhouse::ConfigurableAssert<FAILURE_HANDLER_TYPE>::Failure(stm.str()); \
  } \
  if(wrong_exception) \
  { \
    std::ostringstream stm; \
    stm << "Expected " << #EXCEPTION_TYPE << ". Wrong exception was thrown."; \
    ::snowhouse::ConfigurableAssert<FAILURE_HANDLER_TYPE>::Failure(stm.str()); \
  } \
}

#ifndef SNOWHOUSE_NO_MACROS

#define AssertThrows(EXCEPTION_TYPE, METHOD) SNOWHOUSE_ASSERT_THROWS(EXCEPTION_TYPE, (METHOD), ::snowhouse::DefaultFailureHandler)

#endif // SNOWHOUSE_NO_MACROS

#endif
