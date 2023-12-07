#ifndef BANDIT_REPORTERS_SPEC_H
#define BANDIT_REPORTERS_SPEC_H

#include <bandit/reporters/colored_base.h>
#include <bandit/reporters/summary.h>

namespace bandit {
  namespace reporter {
    struct spec : public colored_base {
      spec(std::ostream& stm, const detail::failure_formatter_t& formatter,
          const detail::colorizer_t& colorizer)
          : colored_base(stm, formatter, colorizer), indentation_(0) {}

      spec(const detail::failure_formatter_t& formatter, const detail::colorizer_t& colorizer)
          : spec(std::cout, formatter, colorizer) {}

      spec& operator=(const spec&) {
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

      void context_starting(const std::string& desc) override {
        progress_base::context_starting(desc);

        stm_ << indent();
        stm_ << "describe " << desc << std::endl;
        increase_indent();
        stm_.flush();
      }

      void context_ended(const std::string& desc) override {
        progress_base::context_ended(desc);
        decrease_indent();
      }

      void it_starting(const std::string& desc) override {
        progress_base::it_starting(desc);
        stm_ << indent() << "- it " << desc << " ... ";
        stm_.flush();
      }

      void it_succeeded(const std::string& desc) override {
        progress_base::it_succeeded(desc);
        stm_ << colorizer_.good();
        stm_ << "OK";
        stm_ << colorizer_.reset();
        stm_ << std::endl;
        stm_.flush();
      }

      void it_failed(const std::string& desc, const detail::assertion_exception& ex) override {
        progress_base::it_failed(desc, ex);
        stm_ << colorizer_.bad();
        stm_ << "FAILED";
        stm_ << colorizer_.reset();
        stm_ << std::endl;
        stm_.flush();
      }

      void it_unknown_error(const std::string& desc) override {
        progress_base::it_unknown_error(desc);
        stm_ << colorizer_.bad();
        stm_ << "ERROR";
        stm_ << colorizer_.reset();
        stm_ << std::endl;
        stm_.flush();
      }

      void it_skip(const std::string& desc) override {
        progress_base::it_skip(desc);
        stm_ << indent() << "- it " << desc << " ... SKIPPED" << std::endl;
        stm_.flush();
      }

    private:
      void increase_indent() {
        indentation_++;
      }

      void decrease_indent() {
        indentation_--;
      }

      std::string indent() {
        return std::string(indentation_, '\t');
      }

    private:
      int indentation_;
    };
  }
}
#endif
