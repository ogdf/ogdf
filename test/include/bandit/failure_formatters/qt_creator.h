#ifndef BANDIT_FAILURE_FORMATTERS_QT_CREATOR_H
#define BANDIT_FAILURE_FORMATTERS_QT_CREATOR_H

#include <bandit/failure_formatters/generic.h>

namespace bandit {
  namespace failure_formatter {
    struct qt_creator : public generic {
      std::string format(const detail::assertion_exception& err) const override {
        return generic_format(err, "file://", ":", "", ": ", "");
      }
    };
  }
}
#endif
