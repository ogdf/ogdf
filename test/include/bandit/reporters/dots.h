#ifndef BANDIT_REPORTERS_DOTS_H
#define BANDIT_REPORTERS_DOTS_H

#include <bandit/reporters/colored_base.h>
#include <bandit/reporters/summary.h>

namespace bandit {
  namespace reporter {
    struct dots : public colored_base {
      dots(std::ostream& stm, const detail::failure_formatter_t& formatter,
          const detail::colorizer_t& colorizer)
          : colored_base(stm, formatter, colorizer) {}

      dots(const detail::failure_formatter_t& formatter, const detail::colorizer_t& colorizer)
          : dots(std::cout, formatter, colorizer) {}

      dots& operator=(const dots&) {
        return *this;
      }

      void test_run_complete() override {
        progress_base::test_run_complete();

        stm_ << std::endl;

        summary::write(stm_, specs_run_, specs_failed_, specs_succeeded_, specs_skipped_, failures_,
            test_run_errors_, colorizer_);
        stm_.flush();
      }

      void test_run_error(const std::string& desc, const detail::test_run_error& err) override {
        progress_base::test_run_error(desc, err);

        std::stringstream ss;
        ss << std::endl;
        ss << "Failed to run \"" << current_context_name() << "\": error \"" << err.what() << "\"" << std::endl;

        test_run_errors_.push_back(ss.str());
      }

      void it_succeeded(const std::string& desc) override {
        progress_base::it_succeeded(desc);
        stm_ << colorizer_.good() << "." << colorizer_.reset();
        stm_.flush();
      }

      void it_failed(const std::string& desc, const detail::assertion_exception& ex) override {
        progress_base::it_failed(desc, ex);
        stm_ << colorizer_.bad() << "F" << colorizer_.reset();
        stm_.flush();
      }

      void it_skip(const std::string& desc) override {
        progress_base::it_skip(desc);
        stm_ << colorizer_.neutral() << 'S' << colorizer_.reset();
        stm_.flush();
      }

      void it_unknown_error(const std::string& desc) override {
        progress_base::it_unknown_error(desc);
        stm_ << colorizer_.bad() << "E" << colorizer_.reset();
        stm_.flush();
      }
    };
  }
}
#endif
