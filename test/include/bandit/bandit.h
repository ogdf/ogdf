#ifndef BANDIT_BANDIT_H
#define BANDIT_BANDIT_H

#include <cassert>
#include <functional>
#include <iostream>
#include <list>
#include <deque>
#include <stdexcept>

#define BANDIT_VERSION "1.1.4"

namespace bandit { namespace detail {
  typedef std::function<void ()> voidfunc_t;
}}

#include <bandit/assertion_frameworks/snowhouse/snowhouse/snowhouse.h>
using namespace snowhouse;

#include <bandit/external/optionparser.h>
#include <bandit/options.h>
#include <bandit/test_run_error.h>
#include <bandit/registration/registration.h>
#include <bandit/assertion_exception.h>
#include <bandit/failure_formatters/failure_formatters.h>
#include <bandit/adapters/adapters.h>
#include <bandit/listener.h>
#include <bandit/reporters/reporters.h>
#include <bandit/context.h>
#include <bandit/run_policies/run_policies.h>
#include <bandit/grammar.h>
#include <bandit/runner.h>

#endif
