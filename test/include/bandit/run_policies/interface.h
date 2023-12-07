#ifndef BANDIT_RUN_POLICIES_INTERFACE_H
#define BANDIT_RUN_POLICIES_INTERFACE_H

#include <stdexcept>
#include <bandit/context.h>

namespace bandit {
  namespace run_policy {
    struct interface {
      interface() : encountered_failure_(false) {}
      interface(const interface& other) = default;

#ifndef _MSC_VER
      interface(interface&&) = default;
#else
      interface(interface&& other) : encountered_failure_(other.encountered_failure_) {}
#endif

      virtual ~interface() = default;

      virtual bool should_run(const std::string& it_name, const context::stack_t& contexts) const = 0;

      virtual void encountered_failure() {
        encountered_failure_ = true;
      }

      virtual bool has_encountered_failure() const {
        return encountered_failure_;
      }

    private:
      bool encountered_failure_;
    };
  }

  namespace detail {
    using run_policy_t = run_policy::interface;
  }
}
#endif
