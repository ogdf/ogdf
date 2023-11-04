#ifndef BANDIT_COLORIZERS_INTERFACE_H
#define BANDIT_COLORIZERS_INTERFACE_H

#include <string>

namespace bandit {
  namespace colorizer {
    struct interface {
      virtual ~interface() = default;

      virtual const std::string good() const = 0;
      virtual const std::string neutral() const = 0;
      virtual const std::string info() const = 0;
      virtual const std::string bad() const = 0;
      virtual const std::string emphasize() const = 0;
      virtual const std::string reset() const = 0;
    };
  }

  namespace detail {
    using colorizer_t = ::bandit::colorizer::interface;
  }
}
#endif
