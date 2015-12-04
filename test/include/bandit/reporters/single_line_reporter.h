#ifndef BANDIT_REPORTERS_SINGLE_LINE_REPORTER_H
#define BANDIT_REPORTERS_SINGLE_LINE_REPORTER_H

namespace bandit { namespace detail {

  struct single_line_reporter : public progress_reporter
  {
    single_line_reporter(std::ostream& stm, const failure_formatter& failure_formatter,
        const detail::colorizer& colorizer)
      : progress_reporter(failure_formatter), stm_(stm), colorizer_(colorizer)
    {}

    single_line_reporter(const failure_formatter& failure_formatter,
        const detail::colorizer& colorizer)
      : progress_reporter(failure_formatter), stm_(std::cout), colorizer_(colorizer)
    {}

	single_line_reporter& operator=(const single_line_reporter&) { return *this; }

    void test_run_complete()
    {
      progress_reporter::test_run_complete();

      stm_ << std::endl;

      test_run_summary summary(specs_run_, specs_failed_, specs_succeeded_, specs_skipped_, failures_,
          test_run_errors_, colorizer_);
      summary.write(stm_);
    }

    void test_run_error(const char* desc, const struct test_run_error& err)
    {
      progress_reporter::test_run_error(desc, err);

      std::stringstream ss;
      ss << std::endl;
      ss << "Failed to run \"" << current_context_name() << "\": error \"" << err.what() << "\"" << std::endl;

      test_run_errors_.push_back(ss.str());
    }

    void it_starting(const char* desc)
    {
      print_status_line();
      progress_reporter::it_starting(desc);
    }

    void it_succeeded(const char* desc)
    {
      progress_reporter::it_succeeded(desc);
      print_status_line();
    }

    void it_failed(const char* desc, const assertion_exception& ex)
    {
      progress_reporter::it_failed(desc, ex);
      print_status_line();
    }

    void it_unknown_error(const char* desc)
    {
      progress_reporter::it_unknown_error(desc);
      print_status_line();
    }

    void it_list(const char* desc)
    {
      progress_reporter::it_list(desc);
      stm_ << current_context_name() << " " << desc << std::endl;
    }

    private:
    void print_status_line()
    {
      stm_ << '\r';
      stm_ << "Executed " << specs_run_ << " tests.";

      if(specs_failed_)
      {
        stm_ << " " << specs_succeeded_ << " succeeded. " << colorizer_.red() << specs_failed_ <<
          " failed." << colorizer_.reset();
      }
      stm_.flush();
    }

    private:
    std::ostream& stm_;
    const detail::colorizer& colorizer_;
  };
}}

#endif
