#ifndef BANDIT_SKIP_POLICY_H
#define BANDIT_SKIP_POLICY_H

namespace bandit {

  struct skip_policy
  {
    virtual bool should_skip(const char* name) const = 0;
  };
  typedef std::unique_ptr<skip_policy> skip_policy_ptr;

  namespace detail {

    inline skip_policy& registered_skip_policy(skip_policy* policy = NULL)
    {
      static struct skip_policy* policy_;

      if(policy)
      {
        policy_ = policy;
      }

      return *policy_;
    }
  }

}

#endif
