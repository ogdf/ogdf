#ifndef BANDIT_REPORTERS_XUNIT_H
#define BANDIT_REPORTERS_XUNIT_H

#include <iostream>
#include <bandit/reporters/progress_base.h>

namespace bandit {
  namespace reporter {
    struct xunit : public progress_base {
      xunit(std::ostream& stm, const detail::failure_formatter_t& formatter, bool report_timing)
          : progress_base(formatter), stm_(stm), report_timing_(report_timing) {
      }

      xunit(const detail::failure_formatter_t& formatter, bool report_timing)
          : xunit(std::cout, formatter, report_timing) {}

      void it_starting(const std::string& desc) override {
        progress_base::it_starting(desc);
        work_stm_ << "\t<testcase classname=\"" << escape(current_context_name()) << "\" ";
        work_stm_ << "name=\"" << escape(desc) << "\"";
      }

      void it_succeeded(const std::string& desc) override {
        progress_base::it_succeeded(desc);
        print_remaining_header_with_time();
        work_stm_ << "\t</testcase>\n";
      }

      void it_failed(const std::string& desc, const detail::assertion_exception& ex) override {
        progress_base::it_failed(desc, ex);
        print_remaining_header_with_time();
        work_stm_ << "\t\t<failure message=\"" << escape(failure_formatter_.format(ex)) << "\" />\n";
        work_stm_ << "\t</testcase>\n";
      }

      void it_unknown_error(const std::string& desc) override {
        progress_base::it_unknown_error(desc);
        print_remaining_header_with_time();
        work_stm_ << "\t\t<failure message=\"Unknown exception\" />\n";
        work_stm_ << "\t</testcase>\n";
      }

      void it_skip(const std::string& desc) override {
        progress_base::it_skip(desc);
        work_stm_ << "\t<testcase classname=\"" << escape(current_context_name()) << "\" ";
        work_stm_ << "name=\"" << escape(desc) << "\" time=\"0.0\">\n";
        work_stm_ << "\t\t<skipped />\n";
        work_stm_ << "\t</testcase>\n";
      }

      void test_run_complete() override {
        stm_ << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
        stm_ << "<testsuite name=\"bandit\" tests=\"" << specs_run_ << "\" errors=\"0\" failures=\""
             << specs_failed_ << "\"";

        if (specs_skipped_ > 0) {
          stm_ << " skipped=\"" << specs_skipped_ << "\"";
        }

        std::chrono::duration<double> dur_in_sec(testsuite_runtime_);
        stm_ << " time=\""
             << time_to_string(report_timing_ ? dur_in_sec.count() : 0.000)
             << "\">\n";

        stm_ << work_stm_.str();

        stm_ << "</testsuite>\n";
      }

    private:
      std::string escape(const std::string& str) {
        std::stringstream stm;

        for (char c : str) {
          switch (c) {
          case '&':
            stm << "&amp;";
            break;
          case '<':
            stm << "&lt;";
            break;
          case '>':
            stm << "&gt;";
            break;
          case '\\':
            stm << "&apos;";
            break;
          case '\"':
            stm << "&quot;";
            break;
          default:
            stm << c;
          }
        }

        return stm.str();
      }

      void print_remaining_header_with_time() {
        work_stm_ << " time=\""
                  << time_to_string(report_timing_ ? testcase_duration_.count() : 0.000)
                  << "\">\n";
      }

      std::ostream& stm_;
      std::stringstream work_stm_;
      bool report_timing_;
    };
  }
}
#endif
