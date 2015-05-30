#ifndef BANDIT_TEST_RUN_ERROR
#define BANDIT_TEST_RUN_ERROR

namespace bandit { namespace detail {

  struct test_run_error : public std::runtime_error
  {
    test_run_error(const char* message) : std::runtime_error(message) {}
  };
}}

#endif
