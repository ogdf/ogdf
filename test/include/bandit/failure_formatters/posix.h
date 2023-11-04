#ifndef BANDIT_FAILURE_FORMATTERS_POSIX_H
#define BANDIT_FAILURE_FORMATTERS_POSIX_H

#include <bandit/failure_formatters/generic.h>

namespace bandit {
  namespace failure_formatter {
    struct posix : public generic {
      std::string format(const detail::assertion_exception& err) const override {
        return generic_format(err, "", ":", "", ": ", "");
      }
    };
  }
}
#endif
