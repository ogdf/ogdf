#ifndef BANDIT_COLORIZERS_DARK_H
#define BANDIT_COLORIZERS_DARK_H

#include <bandit/colorizers/backend.h>
#include <bandit/colorizers/interface.h>

namespace bandit {
  namespace colorizer {
    struct dark : interface, private backend {
      const std::string good() const override {
        return dark_green();
      }

      const std::string neutral() const override {
        return brown();
      }

      const std::string info() const override {
        return dark_blue();
      }

      const std::string bad() const override {
        return dark_red();
      }

      const std::string emphasize() const override {
        return dark_gray();
      }

      const std::string reset() const override {
        return backend::reset();
      }
    };
  }
}
#endif
