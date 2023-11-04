#ifndef BANDIT_REGISTRATION_SPEC_REGISTRY_H
#define BANDIT_REGISTRATION_SPEC_REGISTRY_H

#include <functional>
#include <vector>

namespace bandit {
  namespace detail {
    using spec_registry = std::vector<std::function<void()>>;

    inline detail::spec_registry& specs() {
      static detail::spec_registry registry;
      return registry;
    }
  }
}
#endif
