#ifndef BANDIT_RUN_POLICIES_FILTER_CHAIN_H
#define BANDIT_RUN_POLICIES_FILTER_CHAIN_H

#include <string>
#include <vector>

namespace bandit {
  namespace run_policy {
    struct filter_chain_element {
      std::string pattern;
      bool skip;
    };

    using filter_chain_t = std::vector<filter_chain_element>;
  }
}
#endif
