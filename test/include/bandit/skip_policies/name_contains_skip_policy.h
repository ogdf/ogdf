#ifndef BANDIT_NAME_CONTAINS_SKIP_POLICY_H
#define BANDIT_NAME_CONTAINS_SKIP_POLICY_H

namespace bandit { namespace detail {
  struct name_contains_skip_policy : public skip_policy
  {
    name_contains_skip_policy(const char* pattern)
      : pattern_(pattern)
    {}

    bool should_skip(const char* name) const
    {
      if(pattern_.size() == 0)
      {
        return false;
      }

      std::string n(name);
      bool skip = n.find(pattern_) != std::string::npos;
      return skip;
    }

    private:
    const std::string pattern_;
  };
}}

#endif
