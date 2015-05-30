#ifndef BANDIT_NEVER_RUN_POLICY_H
#define BANDIT_NEVER_RUN_POLICY_H

namespace bandit { namespace detail {

  struct never_run_policy : public run_policy
  {
    bool should_run(const char* /* it_name */, const contextstack_t& /* contexts */) const
    {
      return false;
    }
  };
}}
#endif
