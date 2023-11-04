#ifndef BANDIT_REPORTERS_SUMMARY_H
#define BANDIT_REPORTERS_SUMMARY_H

#include <algorithm>
#include <vector>
#include <iostream>
#include <bandit/colorizers/interface.h>

namespace bandit {
  namespace reporter {
    struct summary {
      summary() = delete;

      static void write(std::ostream& stm,
          int specs_run, int specs_failed, int specs_succeeded, int specs_skipped,
          const std::vector<std::string>& failures, const std::vector<std::string>& errors,
          const detail::colorizer_t& colorizer) {
        if (specs_run == 0 && errors.size() == 0) {
          stm << colorizer.bad()
              << "Could not find any tests."
              << colorizer.reset()
              << std::endl;
          return;
        }

        if (specs_failed == 0 && errors.size() == 0) {
          stm << colorizer.good()
              << "Success!"
              << colorizer.reset()
              << std::endl;
        }

        if (errors.size() > 0) {
          for (const auto& error : errors) {
            stm << error << std::endl;
          }
        }

        if (specs_failed > 0) {
          stm << colorizer.bad()
              << "There were failures!"
              << colorizer.reset() << std::endl;
          for (const auto& failure : failures) {
            stm << failure << std::endl;
          }
        }

        stm << "Test run complete. "
            << specs_run << " tests run. "
            << specs_succeeded << " succeeded.";

        if (specs_skipped > 0) {
          stm << " " << specs_skipped << " skipped.";
        }

        if (specs_failed > 0) {
          stm << " " << specs_failed << " failed.";
        }

        if (errors.size() > 0) {
          stm << " " << errors.size() << " test run errors.";
        }

        stm << std::endl;
      }
    };
  }
}
#endif
