#ifndef BANDIT_SPEC_REPORTER_H
#define BANDIT_SPEC_REPORTER_H

namespace bandit { namespace detail {

  struct spec_reporter : public progress_reporter
  {
    spec_reporter(std::ostream& stm, const failure_formatter& failure_formatter,
        const detail::colorizer& colorizer)
      : progress_reporter(failure_formatter),  stm_(stm), colorizer_(colorizer), indentation_(0)
    {}

    spec_reporter(const failure_formatter& failure_formatter, const detail::colorizer& colorizer)
      : progress_reporter(failure_formatter), stm_(std::cout), colorizer_(colorizer), indentation_(0)
    {}

	spec_reporter& operator=(const spec_reporter&) { return *this; }

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

    virtual void context_starting(const char* desc)
    {
      progress_reporter::context_starting(desc);

      stm_ << indent();
      stm_ << "describe " << desc << std::endl;
      increase_indent();
      stm_.flush();

    }

    virtual void context_ended(const char* desc)
    {
      progress_reporter::context_ended(desc);
      decrease_indent();
    }

    virtual void it_starting(const char* desc)
    {
      progress_reporter::it_starting(desc);
      stm_ << indent() << "- it " << desc << " ... ";
      stm_.flush();
    }

    virtual void it_succeeded(const char* desc)
    {
      progress_reporter::it_succeeded(desc);
      stm_ << colorizer_.green();
	  stm_ << "OK";
	  stm_ << colorizer_.reset();
	  stm_ << std::endl;
      stm_.flush();
    }

    virtual void it_failed(const char* desc, const assertion_exception& ex)
    {
      progress_reporter::it_failed(desc, ex);
      stm_ << colorizer_.red();
	  stm_ << "FAILED";
	  stm_ << colorizer_.reset();
	  stm_ << std::endl;
      stm_.flush();
    }

    virtual void it_unknown_error(const char* desc)
    {
      progress_reporter::it_unknown_error(desc);
      stm_ << colorizer_.red();
	  stm_ << "ERROR";
	  stm_ << colorizer_.reset();
	  stm_ << std::endl;
      stm_.flush();
    }

    virtual void it_skip(const char* desc)
    {
      progress_reporter::it_skip(desc);
      stm_ << indent() << "- it " << desc << " ... SKIPPED" << std::endl;
      stm_.flush();
    }

    virtual void it_list(const char* desc)
    {
      progress_reporter::it_skip(desc);
      stm_ << indent() << "- it " << desc << std::endl;
      stm_.flush();
    }

    private:
    void increase_indent()
    {
      indentation_++;
    }

    void decrease_indent()
    {
      indentation_--;
    }

    std::string indent()
    {
      return std::string(indentation_, '\t');
    }

    private:
    std::ostream& stm_;
    const detail::colorizer& colorizer_;
    int indentation_;
  };
}}

#endif
