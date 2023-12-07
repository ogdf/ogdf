#ifndef BANDIT_RUN_POLICIES_NEVER_H
#define BANDIT_RUN_POLICIES_NEVER_H

#include <bandit/run_policies/interface.h>

namespace bandit {
  namespace run_policy {
    struct never : public interface {
      bool should_run(const std::string&, const context::stack_t&) const override {
        return false;
      }
    };
  }
}
#endif
