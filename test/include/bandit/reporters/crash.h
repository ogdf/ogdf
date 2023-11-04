#ifndef BANDIT_REPORTERS_CRASH_H
#define BANDIT_REPORTERS_CRASH_H

#include <iostream>
#include <vector>
#include <bandit/reporters/progress_base.h>

namespace bandit {
  namespace reporter {
    struct crash : public progress_base {
      crash(std::ostream& stm, const detail::failure_formatter_t& formatter)
          : progress_base(formatter), stm_(stm) {}

      crash(const detail::failure_formatter_t& formatter)
          : crash(std::cout, formatter) {}

      crash& operator=(const crash&) {
        return *this;
      }

      void test_run_complete() override {
        progress_base::test_run_complete();
        for (auto failure : failures_) {
          stm_ << std::endl
               << "# FAILED " << failure;
        }
        for (auto error : test_run_errors_) {
          stm_ << std::endl
               << "# ERROR " << error;
        }
        stm_.flush();
      }

      void test_run_error(const std::string& desc, const detail::test_run_error& err) override {
        progress_base::test_run_error(desc, err);
        std::stringstream ss;
        ss << current_context_name() << ": " << desc << ": " << err.what() << std::endl;
        test_run_errors_.push_back(ss.str());
      }

      void it_skip(const std::string& desc) override {
        progress_base::it_skip(desc);
      }

      void it_starting(const std::string& desc) override {
        progress_base::it_starting(desc);
        for (auto context : contexts_) {
          stm_ << context << " | ";
        }
        stm_ << desc << std::endl;
      }

      void it_succeeded(const std::string& desc) override {
        progress_base::it_succeeded(desc);
      }

      void it_failed(const std::string& desc, const detail::assertion_exception& ex) override {
        ++specs_failed_;

        std::stringstream ss;
        ss << current_context_name() << " " << desc << ":" << std::endl
           << failure_formatter_.format(ex);
        failures_.push_back(ss.str());

        stm_ << "FAILED" << std::endl;
      }

      void it_unknown_error(const std::string& desc) override {
        ++specs_failed_;

        std::stringstream ss;
        ss << current_context_name() << " " << desc << ": Unknown exception" << std::endl;
        failures_.push_back(ss.str());

        stm_ << "UNKNOWN EXCEPTION" << std::endl;
      }

    private:
      std::ostream& stm_;
    };
  }
}
#endif
