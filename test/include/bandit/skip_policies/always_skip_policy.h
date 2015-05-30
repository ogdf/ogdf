#ifndef BANDIT_ALWAYS_SKIP_POLICY_H
#define BANDIT_ALWAYS_SKIP_POLICY_H

namespace bandit { namespace detail {

  struct always_skip_policy : public skip_policy
  {
    bool should_skip(const char*) const
    {
      return true;
    }
  };
}}

#endif
