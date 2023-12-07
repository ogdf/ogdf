#ifndef BANDIT_FAILURE_FORMATTERS_INTERFACE_H
#define BANDIT_FAILURE_FORMATTERS_INTERFACE_H

#include <bandit/assertion_exception.h>

namespace bandit {
  namespace failure_formatter {
    struct interface {
      virtual ~interface() = default;

      virtual std::string format(const detail::assertion_exception&) const = 0;
    };
  }

  namespace detail {
    using failure_formatter_t = ::bandit::failure_formatter::interface;
  }
}
#endif
