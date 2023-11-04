#ifndef BANDIT_REPORTERS_INFO_H
#define BANDIT_REPORTERS_INFO_H

#include <iostream>
#include <bandit/reporters/colored_base.h>

namespace bandit {
  namespace reporter {
    struct info : public colored_base {
      info(std::ostream& stm, const detail::failure_formatter_t& formatter,
          const detail::colorizer_t& colorizer,
          bool report_timing)
          : colored_base(stm, formatter, colorizer),
            active_context_index_(0), report_timing_(report_timing) {}

      info(const detail::failure_formatter_t& formatter,
          const detail::colorizer_t& colorizer,
          bool report_timing)
          : info(std::cout, formatter, colorizer, report_timing) {}

      info& operator=(const info&) {
        return *this;
      }

      void test_run_complete() override {
        colored_base::test_run_complete();
        stm_ << std::endl;
        list_failures_and_errors();
        summary();
        stm_.flush();
      }

      void test_run_error(const std::string& desc, const detail::test_run_error& err) override {
        colored_base::test_run_error(desc, err);

        std::stringstream ss;
        ss << "Failed to run \"" << current_context_name() << "\": error \"" << err.what() << "\"" << std::endl;
        test_run_errors_.push_back(ss.str());
      }

      void context_starting(const std::string& desc) override {
        colored_base::context_starting(desc);
        context_stack_.emplace_back();
        if (context_stack_.size() == 1) {
          output_context_start_message();
        }
      }

      void context_ended(const std::string& desc) override {
        if (context_stack_.size() == 1 || !context_stack_.back().skipped_all()) {
          output_context_end_message();
        }

        colored_base::context_ended(desc);
        const auto context = context_stack_.back();
        context_stack_.pop_back();
        if (!context_stack_.empty()) {
          context_stack_.back().merge(context);
        }
      }

      void it_skip(const std::string& desc) override {
        colored_base::it_skip(desc);
        context_stack_.back().count_skip();
      }

      void it_starting(const std::string& desc) override {
        if (context_stack_.size() > 1 && context_stack_.back().skipped_all()) {
          output_not_yet_shown_context_start_messages();
        }

        colored_base::it_starting(desc);
        stm_
            << indent()
            << colorizer_.neutral()
            << "[ TEST ]"
            << colorizer_.reset()
            << " it " << desc;
        ++active_context_index_;
        stm_.flush();
      }

      void it_succeeded(const std::string& desc) override {
        colored_base::it_succeeded(desc);
        context_stack_.back().count_success();
        --active_context_index_;
        stm_
            << "\r" << indent()
            << colorizer_.good()
            << "[ PASS ]"
            << colorizer_.reset()
            << " it " << desc;
        print_timing();
        stm_ << std::endl;
        stm_.flush();
      }

      void it_failed(const std::string& desc, const detail::assertion_exception& ex) override {
        colored_base::it_failed(desc, ex);
        context_stack_.back().count_failure();
        --active_context_index_;
        stm_
            << "\r" << indent()
            << colorizer_.bad()
            << "[ FAIL ]"
            << colorizer_.reset()
            << " it " << desc;
        print_timing();
        stm_ << std::endl;
        stm_.flush();
      }

      void it_unknown_error(const std::string& desc) override {
        colored_base::it_unknown_error(desc);
        context_stack_.back().count_failure();
        --active_context_index_;
        stm_
            << "\r" << indent()
            << colorizer_.bad()
            << "-ERROR->"
            << colorizer_.reset()
            << " it " << desc;
        print_timing();
        stm_ << std::endl;
        stm_.flush();
      }

    protected:
      struct context_info {
        context_info() : total_(0), skipped_(0), failed_(0) {}

        void merge(const context_info& ci) {
          total_ += ci.total();
          skipped_ += ci.skipped();
          failed_ += ci.failed();
        }

        void count_success() {
          ++total_;
        }

        void count_skip() {
          ++total_;
          ++skipped_;
        }

        void count_failure() {
          ++total_;
          ++failed_;
        }

        bool skipped_all() const {
          return total_ == skipped_;
        }

        int total() const { return total_; }
        int skipped() const { return skipped_; }
        int failed() const { return failed_; }

      private:
        int total_;
        int skipped_;
        int failed_;
      };

      std::string indent() {
        return std::string(2 * active_context_index_, ' ');
      }

      void print_timing() {
          if (report_timing_) {
              std::chrono::duration<double> dur_in_sec(testcase_duration_);
                  stm_
                  << colorizer_.neutral()
                  << " ("
                  << time_to_string(dur_in_sec.count())
                  << " seconds)";
          }
      }

      void list_failures_and_errors() {
        if (specs_failed_ > 0) {
          stm_
              << colorizer_.bad()
              << "List of failures:"
              << std::endl;
          for (const auto& failure : failures_) {
            stm_
                << colorizer_.emphasize()
                << " (*) "
                << colorizer_.bad()
                << failure << std::endl;
          }
        }
        if (test_run_errors_.size() > 0) {
          stm_
              << colorizer_.bad()
              << "List of run errors:"
              << std::endl;
          for (const auto& error : test_run_errors_) {
            stm_
                << colorizer_.emphasize()
                << " (*) "
                << colorizer_.bad()
                << error << std::endl;
          }
        }
      }

      void summary() {
        stm_
            << colorizer_.emphasize()
            << "Tests run: " << specs_run_
            << std::endl;
        if (specs_skipped_ > 0) {
          stm_
              << colorizer_.neutral()
              << "Skipped: " << specs_skipped_
              << std::endl;
        }
        if (specs_succeeded_ > 0) {
          stm_
              << colorizer_.good()
              << "Passed: " << specs_succeeded_
              << std::endl;
        }
        if (specs_failed_ > 0) {
          stm_
              << colorizer_.bad()
              << "Failed: " << specs_failed_
              << std::endl;
        }
        if (test_run_errors_.size() > 0) {
          stm_
              << colorizer_.bad()
              << "Errors: " << test_run_errors_.size()
              << std::endl;
        }
        stm_
            << colorizer_.reset()
            << std::endl;
      }

      void output_context_start_message() {
        stm_
            << indent()
            << colorizer_.info()
            << "begin "
            << colorizer_.emphasize()
            << contexts_[active_context_index_]
            << colorizer_.reset()
            << std::endl;
        ++active_context_index_;
        stm_.flush();
      }

      void output_not_yet_shown_context_start_messages() {
        for (auto i = active_context_index_; i < context_stack_.size(); ++i) {
          output_context_start_message();
        }
      }

      void output_context_end_message() {
        --active_context_index_;
        stm_
            << indent()
            << colorizer_.info()
            << "end "
            << colorizer_.reset()
            << contexts_.back();

        const auto& context = context_stack_.back();
        if (context.total() > 0) {
          stm_
              << colorizer_.emphasize()
              << " " << context.total() << " total";
        }
        if (context.skipped() > 0) {
          stm_
              << colorizer_.neutral()
              << " " << context.skipped() << " skipped";
        }
        if (context.failed() > 0) {
          stm_
              << colorizer_.bad()
              << " " << context.failed() << " failed";
        }
        stm_ << colorizer_.reset() << std::endl;
      }

      std::vector<context_info> context_stack_;
      decltype(context_stack_)::size_type active_context_index_;
      bool report_timing_;
    };
  }
}
#endif
