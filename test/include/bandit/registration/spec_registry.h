#ifndef BANDIT_SPEC_REGISTRY_H
#define BANDIT_SPEC_REGISTRY_H

namespace bandit {
  namespace detail {
    typedef std::list<voidfunc_t> spec_registry;

    inline detail::spec_registry& specs()
    {
      static detail::spec_registry registry;
      return registry;
    }
  }

}

#endif
