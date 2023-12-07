#ifndef BANDIT_COLORIZERS_OFF_H
#define BANDIT_COLORIZERS_OFF_H

#include <bandit/colorizers/interface.h>

namespace bandit {
  namespace colorizer {
    struct off : interface {
      const std::string good() const override {
        return "";
      }

      const std::string neutral() const override {
        return "";
      }

      const std::string info() const override {
        return "";
      }

      const std::string bad() const override {
        return "";
      }

      const std::string emphasize() const override {
        return "";
      }

      const std::string reset() const override {
        return "";
      }
    };
  }
}
#endif
