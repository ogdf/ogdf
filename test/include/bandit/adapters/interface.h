#ifndef BANDIT_ADAPTERS_INTERFACE_H
#define BANDIT_ADAPTERS_INTERFACE_H

#include <functional>

namespace bandit {
  namespace adapter {
    struct interface {
      virtual ~interface() = default;

      virtual void adapt_exceptions(std::function<void()>) = 0;
    };
  }

  namespace detail {
    using assertion_adapter_t = adapter::interface;
  }
}
#endif
