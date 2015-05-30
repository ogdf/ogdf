#ifndef BANDIT_RUN_POLICY_H
#define BANDIT_RUN_POLICY_H

namespace bandit { namespace detail {

  struct run_policy
  {
    virtual ~run_policy() {}
    virtual bool should_run(const char* it_name, const contextstack_t& contexts) const = 0;
  };
  typedef std::unique_ptr<run_policy> run_policy_ptr;

  inline run_policy& registered_run_policy(run_policy* policy = NULL)
  {
    static struct run_policy* policy_;

    if(policy)
    {
      policy_ = policy;
    }

    return *policy_;
  }

}}

#endif
