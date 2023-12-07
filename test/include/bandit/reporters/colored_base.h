#ifndef BANDIT_REPORTERS_COLORED_BASE_H
#define BANDIT_REPORTERS_COLORED_BASE_H

#include <ostream>
#include <bandit/colorizers/interface.h>
#include <bandit/reporters/progress_base.h>

namespace bandit {
  namespace reporter {
    struct colored_base : public progress_base {
      colored_base(std::ostream& stm,
          const detail::failure_formatter_t& formatter,
          const detail::colorizer_t& colorizer)
          : progress_base(formatter), stm_(stm), colorizer_(colorizer) {}

      ~colored_base() override {
        stm_ << colorizer_.reset();
      }

    protected:
      std::ostream& stm_;
      const detail::colorizer_t& colorizer_;
    };
  }
}
#endif
