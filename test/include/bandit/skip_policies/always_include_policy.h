#ifndef BANDIT_ALWAYS_INCLUDE_POLICY_H
#define BANDIT_ALWAYS_INCLUDE_POLICY_H

namespace bandit { namespace detail {

  struct always_include_policy : public skip_policy
  {
    bool should_skip(const char*) const
    {
      return false;
    }
  };

}}

#endif
