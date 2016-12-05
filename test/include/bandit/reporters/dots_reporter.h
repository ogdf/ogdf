#ifndef BANDIT_DOTS_REPORTER_H
#define BANDIT_DOTS_REPORTER_H

namespace bandit { namespace detail {

  struct dots_reporter : public progress_reporter
  {
    dots_reporter(std::ostream& stm, const failure_formatter& failure_formatter,
        const detail::colorizer& colorizer)
      : progress_reporter(failure_formatter), stm_(stm), colorizer_(colorizer)
    {}

    dots_reporter(const failure_formatter& failure_formatter, const detail::colorizer& colorizer)
      : progress_reporter(failure_formatter), stm_(std::cout), colorizer_(colorizer)
    {}

	dots_reporter& operator=(const dots_reporter&) { return *this; }

    void test_run_complete()
    {
      progress_reporter::test_run_complete();

      stm_ << std::endl;

      test_run_summary summary(specs_run_, specs_failed_, specs_succeeded_, specs_skipped_, failures_,
          test_run_errors_, colorizer_);
      summary.write(stm_);
      stm_.flush();
    }

    void test_run_error(const char* desc, const struct test_run_error& err)
    {
      progress_reporter::test_run_error(desc, err);

      std::stringstream ss;
      ss << std::endl;
      ss << "Failed to run \"" << current_context_name() << "\": error \"" << err.what() << "\"" << std::endl;

      test_run_errors_.push_back(ss.str());
    }

    void it_succeeded(const char* desc)
    {
      progress_reporter::it_succeeded(desc);
      stm_ << colorizer_.green() << "." << colorizer_.reset();
      stm_.flush();
    }

    void it_failed(const char* desc, const assertion_exception& ex)
    {
      progress_reporter::it_failed(desc, ex);
      stm_ << colorizer_.red() << "F" << colorizer_.reset();
      stm_.flush();
    }

    void it_skip(const char* desc)
    {
        progress_reporter::it_skip(desc);
        stm_ << colorizer_.yellow() << 'S' << colorizer_.reset();
        stm_.flush();
    }

    void it_unknown_error(const char* desc)
    {
      progress_reporter::it_unknown_error(desc);
      stm_ << colorizer_.red() << "E" << colorizer_.reset();
      stm_.flush();
    }

    void it_list(const char* desc)
    {
      progress_reporter::it_list(desc);
      stm_ << current_context_name() << " " << desc << std::endl;
    }

    private:
    std::ostream& stm_;
    const detail::colorizer& colorizer_;
  };
}}

#endif
