#ifndef BANDIT_REPORTERS_INTERFACE_H
#define BANDIT_REPORTERS_INTERFACE_H

#include <bandit/assertion_exception.h>
#include <bandit/test_run_error.h>

namespace bandit {
  namespace reporter {
    struct interface {
      virtual ~interface() = default;

      virtual void test_run_starting() = 0;
      virtual void test_run_complete() = 0;

      virtual void context_starting(const std::string& desc) = 0;
      virtual void context_ended(const std::string& desc) = 0;
      virtual void test_run_error(const std::string& desc, const detail::test_run_error& error) = 0;

      virtual void it_starting(const std::string& desc) = 0;
      virtual void it_succeeded(const std::string& desc) = 0;
      virtual void it_failed(const std::string& desc, const detail::assertion_exception& ex) = 0;
      virtual void it_unknown_error(const std::string& desc) = 0;
      virtual void it_skip(const std::string& desc) = 0;

      virtual bool did_we_pass() const = 0;
    };
  }

  namespace detail {
    using reporter_t = reporter::interface;
  }
}
#endif
