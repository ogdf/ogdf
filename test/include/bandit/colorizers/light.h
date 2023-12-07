#ifndef BANDIT_COLORIZERS_LIGHT_H
#define BANDIT_COLORIZERS_LIGHT_H

#include <bandit/colorizers/backend.h>
#include <bandit/colorizers/interface.h>

namespace bandit {
  namespace colorizer {
    struct light : interface, private backend {
      const std::string good() const override {
        return green();
      }

      const std::string neutral() const override {
        return yellow();
      }

      const std::string info() const override {
        return blue();
      }

      const std::string bad() const override {
        return red();
      }

      const std::string emphasize() const override {
        return white();
      }

      const std::string reset() const override {
        return backend::reset();
      }
    };
  }
}
#endif
