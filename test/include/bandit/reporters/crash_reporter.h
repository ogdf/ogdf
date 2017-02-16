#ifndef BANDIT_CRASH_REPORTER_H
#define BANDIT_CRASH_REPORTER_H

namespace bandit {
namespace detail {

struct crash_reporter : public progress_reporter
{
	crash_reporter(std::ostream &stm, const failure_formatter &failure_formatter, const detail::colorizer &)
	  : progress_reporter(failure_formatter)
	  , stm_(stm)
	  , contexts_()
	{
	}

	crash_reporter(const failure_formatter &failure_formatter, const detail::colorizer &colorizer)
	  : crash_reporter(std::cout, failure_formatter, colorizer)
	{
	}

	crash_reporter &operator=(const crash_reporter &)
	{
		return *this;
	}

	void test_run_complete() override
	{
		progress_reporter::test_run_complete();
		for (auto failure : failures_) {
			stm_ << std::endl << "# FAILED " << failure;
		}
		for (auto error : test_run_errors_) {
			stm_ << std::endl << "# ERROR " << error;
		}
		stm_.flush();
	}

	void test_run_error(const char *desc, const struct test_run_error &err) override
	{
		progress_reporter::test_run_error(desc, err);
		std::stringstream ss;
		ss << current_context_name() << " " << desc << ":" << err.what() << std::endl;
		test_run_errors_.push_back(ss.str());
	}

	void context_starting(const char *desc) override
	{
		progress_reporter::context_starting(desc);
		contexts_.emplace_back(desc);
	}

	void context_ended(const char *desc) override
	{
		progress_reporter::context_ended(desc);
		contexts_.pop_back();
	}

	void it_skip(const char *desc) override
	{
		progress_reporter::it_skip(desc);
	}

	void it_starting(const char *desc) override
	{
		progress_reporter::it_starting(desc);
		for (auto context : contexts_) {
			stm_ << context << " | ";
		}
		stm_ << desc << std::endl;
	}

	void it_succeeded(const char *desc) override
	{
		progress_reporter::it_succeeded(desc);
	}

	void it_failed(const char *desc, const assertion_exception &ex) override
	{
		++specs_failed_;

		std::stringstream ss;
		ss << current_context_name() << " " << desc << ":" << std::endl << failure_formatter_.format(ex);
		failures_.push_back(ss.str());

		stm_ << "FAILED" << std::endl;
	}

	void it_unknown_error(const char *desc) override
	{
		++specs_failed_;

		std::stringstream ss;
		ss << current_context_name() << " " << desc << ": Unknown exception" << std::endl;
		failures_.push_back(ss.str());

		stm_ << "UNKNOWN EXCEPTION" << std::endl;
	}

	bool did_we_pass() const override
	{
		return specs_failed_ == 0 && test_run_errors_.size() == 0;
	}

private:
	std::ostream &stm_;
	std::vector<const char *> contexts_;
};
}
}

#endif
