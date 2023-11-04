#ifndef BANDIT_GRAMMAR_H
#define BANDIT_GRAMMAR_H

#include <bandit/controller.h>

namespace bandit {
  inline void describe(const std::string& desc, const std::function<void()>& func,
      bool hard_skip = false,
      detail::controller_t& controller = detail::registered_controller()) {
    controller.describe(desc, func, hard_skip);
  }

  inline void describe_skip(const std::string& desc, const std::function<void()>& func,
      detail::controller_t& controller = detail::registered_controller()) {
    describe(desc, func, true, controller);
  }

  inline void xdescribe(const std::string& desc, const std::function<void()>& func,
      detail::controller_t& controller = detail::registered_controller()) {
    describe(desc, func, true, controller);
  }

  inline void before_each(const std::function<void()>& func,
      detail::controller_t& controller = detail::registered_controller()) {
    controller.before_each(func);
  }

  inline void after_each(const std::function<void()>& func,
      detail::controller_t& controller = detail::registered_controller()) {
    controller.after_each(func);
  }

  inline void it(const std::string& desc, const std::function<void()>& func,
      bool hard_skip = false, detail::controller_t& controller = detail::registered_controller()) {
    controller.it(desc, func, hard_skip);
  }

  inline void it_skip(const std::string& desc, const std::function<void()>& func,
      detail::controller_t& controller = detail::registered_controller()) {
    it(desc, func, true, controller);
  }

  inline void xit(const std::string& desc, const std::function<void()>& func,
      detail::controller_t& controller = detail::registered_controller()) {
    it(desc, func, true, controller);
  }
}
#endif
