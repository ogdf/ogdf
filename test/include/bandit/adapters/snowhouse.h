#ifndef BANDIT_ADAPTERS_SNOWHOUSE_H
#define BANDIT_ADAPTERS_SNOWHOUSE_H

#include <bandit/adapters/interface.h>
#include <bandit/assertion_exception.h>
#include <bandit/assertion_frameworks/snowhouse/snowhouse.h>

namespace bandit {
  namespace adapter {
    struct snowhouse : public interface {
      void adapt_exceptions(std::function<void()> func) override {
        try {
          func();
        } catch (const ::snowhouse::AssertionException& ex) {
          throw detail::assertion_exception(ex.what(), ex.file(), ex.line());
        }
      }
    };
  }
}
#endif
