#ifndef BANDIT_REGISTRATION_REGISTRAR_H
#define BANDIT_REGISTRATION_REGISTRAR_H

#include <bandit/registration/spec_registry.h>

namespace bandit {
  namespace detail {
    struct spec_registrar {
      spec_registrar(std::function<void()> func) {
        bandit::detail::specs().push_back(func);
      }
    };
  }
}

#define BANDIT_CONCAT2(a, b) a##b
#define BANDIT_CONCAT(a, b) BANDIT_CONCAT2(a, b)
#define BANDIT_ADD_COUNTER(a) BANDIT_CONCAT(a, __COUNTER__)

#define go_bandit \
  static bandit::detail::spec_registrar BANDIT_ADD_COUNTER(bandit_registrar_)

#define SPEC_BEGIN(name) \
  go_bandit([]{
#define SPEC_END \
  });
#endif
