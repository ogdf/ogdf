#ifndef BANDIT_REGISTRAR_H
#define BANDIT_REGISTRAR_H

namespace bandit { namespace detail {

    struct spec_registrar
    {
      spec_registrar( bandit::detail::voidfunc_t func)
      {
        bandit::detail::specs().push_back(func);
      }
    };

}}

#define go_bandit \
  static bandit::detail::spec_registrar bandit_registrar

#endif
