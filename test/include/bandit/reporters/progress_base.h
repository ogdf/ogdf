#ifndef BANDIT_REPORTERS_PROGRESS_BASE_H
#define BANDIT_REPORTERS_PROGRESS_BASE_H

#include <algorithm>
#include <sstream>
#include <vector>
#include <chrono>
#include <bandit/failure_formatters/interface.h>
#include <bandit/reporters/interface.h>

namespace bandit {
  namespace reporter {
    struct progress_base : public interface {
      progress_base(const detail::failure_formatter_t& formatter)
          : specs_run_(0), specs_succeeded_(0), specs_failed_(0), specs_skipped_(0),
            failure_formatter_(formatter),
            testcase_start_time_point_(), testcase_duration_(), testsuite_runtime_() {}

      progress_base& operator=(const progress_base&) {
        return *this;
      }

      void test_run_starting() override {
        specs_run_ = 0;
        specs_succeeded_ = 0;
        specs_failed_ = 0;
        specs_skipped_ = 0;
        failures_.clear();
        contexts_.clear();
        testsuite_runtime_ = std::chrono::nanoseconds(0);
      }

      void test_run_complete() override {}

      void context_starting(const std::string& desc) override {
        contexts_.push_back(std::string(desc));
      }

      void context_ended(const std::string&) override {
        contexts_.pop_back();
      }

      void test_run_error(const std::string&, const detail::test_run_error&) override {}

      void it_starting(const std::string&) override {
        specs_run_++;
        testcase_start_time_point_ = std::chrono::high_resolution_clock::now();
      }

      void it_succeeded(const std::string&) override {
        specs_succeeded_++;
        update_test_duration();
      }

      void it_failed(const std::string& desc, const detail::assertion_exception& ex) override {
        specs_failed_++;
        update_test_duration();

        std::stringstream ss;
        ss << current_context_name() << " " << desc << ":" << std::endl;
        ss << failure_formatter_.format(ex);

        failures_.push_back(ss.str());
      }

      void it_unknown_error(const std::string& desc) override {
        specs_failed_++;
        update_test_duration();

        std::stringstream ss;
        ss << current_context_name() << " " << desc << ":" << std::endl;
        ss << "Unknown exception";
        ss << std::endl;

        failures_.push_back(ss.str());
      }

      void it_skip(const std::string&) override {
        specs_skipped_++;
      }

      bool did_we_pass() const override final {
        return specs_failed_ == 0 && test_run_errors_.size() == 0;
      }

    protected:
      std::string current_context_name() {
        std::string name;

        for (const auto& context : contexts_) {
          if (name.size() > 0) {
            name += " ";
          }

          name += context;
        }

        return name;
      }

    protected:
      void update_test_duration() {
        testcase_duration_ = std::chrono::high_resolution_clock::now() - testcase_start_time_point_;
        testsuite_runtime_ += std::chrono::duration_cast<std::chrono::nanoseconds>(testcase_duration_);
      }

      std::string time_to_string(double t) {
        std::stringstream ss;
        ss.precision(3);
        ss << std::fixed << t;
        return ss.str();
      }

      int specs_run_;
      int specs_succeeded_;
      int specs_failed_;
      int specs_skipped_;
      const detail::failure_formatter_t& failure_formatter_;
      std::vector<std::string> contexts_;
      std::vector<std::string> failures_;
      std::vector<std::string> test_run_errors_;
      std::chrono::high_resolution_clock::time_point testcase_start_time_point_;
      std::chrono::duration<double> testcase_duration_;
      std::chrono::nanoseconds testsuite_runtime_;
    };
  }
}
#endif
