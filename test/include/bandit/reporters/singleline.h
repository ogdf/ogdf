#ifndef BANDIT_REPORTERS_SINGLELINE_H
#define BANDIT_REPORTERS_SINGLELINE_H

#include <bandit/reporters/colored_base.h>
#include <bandit/reporters/summary.h>

namespace bandit {
  namespace reporter {
    struct singleline : public colored_base {
      singleline(std::ostream& stm, const detail::failure_formatter_t& formatter,
          const detail::colorizer_t& colorizer)
          : colored_base(stm, formatter, colorizer) {}

      singleline(const detail::failure_formatter_t& formatter,
          const detail::colorizer_t& colorizer)
          : singleline(std::cout, formatter, colorizer) {}

      singleline& operator=(const singleline&) {
        return *this;
      }

      void test_run_complete() override {
        progress_base::test_run_complete();

        stm_ << std::endl;

        summary::write(stm_, specs_run_, specs_failed_, specs_succeeded_, specs_skipped_, failures_,
            test_run_errors_, colorizer_);
      }

      void test_run_error(const std::string& desc, const detail::test_run_error& err) override {
        progress_base::test_run_error(desc, err);

        std::stringstream ss;
        ss << std::endl;
        ss << "Failed to run \"" << current_context_name() << "\": error \"" << err.what() << "\"" << std::endl;

        test_run_errors_.push_back(ss.str());
      }

      void it_starting(const std::string& desc) override {
        print_status_line();
        progress_base::it_starting(desc);
      }

      void it_succeeded(const std::string& desc) override {
        progress_base::it_succeeded(desc);
        print_status_line();
      }

      void it_failed(const std::string& desc, const detail::assertion_exception& ex) override {
        progress_base::it_failed(desc, ex);
        print_status_line();
      }

      void it_unknown_error(const std::string& desc) override {
        progress_base::it_unknown_error(desc);
        print_status_line();
      }

    private:
      void print_status_line() {
        stm_ << '\r';
        stm_ << "Executed " << specs_run_ << " tests.";

        if (specs_failed_) {
          stm_
              << " " << specs_succeeded_ << " succeeded. "
              << colorizer_.bad() << specs_failed_ << " failed."
              << colorizer_.reset();
        }
        stm_.flush();
      }
    };
  }
}
#endif
