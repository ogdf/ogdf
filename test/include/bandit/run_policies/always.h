#ifndef BANDIT_RUN_POLICIES_ALWAYS_H
#define BANDIT_RUN_POLICIES_ALWAYS_H

#include <bandit/run_policies/interface.h>

namespace bandit {
  namespace run_policy {
    struct always : public interface {
      bool should_run(const std::string&, const context::stack_t&) const override {
        return true;
      }
    };
  }
}
#endif
