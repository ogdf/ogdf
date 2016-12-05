#ifndef BANDIT_RUN_POLICY_H
#define BANDIT_RUN_POLICY_H

namespace bandit { namespace detail {

  struct run_policy
  {
    run_policy() : encountered_failure_(false) {}
    run_policy(const run_policy& other) = default;

#ifndef _MSC_VER
    run_policy(run_policy&&) = default;
#else
	run_policy(run_policy&& other) : encountered_failure_(other.encountered_failure_) {}
#endif

    virtual ~run_policy() {}

    virtual bool should_run(const char* it_name, const contextstack_t& contexts) const = 0;

    virtual void encountered_failure()
    {
        encountered_failure_ = true;
    }

    virtual bool has_encountered_failure() const
    {
        return encountered_failure_;
    }

  private:
    bool encountered_failure_;
  };
  typedef std::unique_ptr<run_policy> run_policy_ptr;

  inline run_policy& registered_run_policy(run_policy* policy = nullptr)
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
