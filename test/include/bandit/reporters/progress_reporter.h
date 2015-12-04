#ifndef BANDIT_PROGRESS_REPORTER_H
#define BANDIT_PROGRESS_REPORTER_H

namespace bandit { namespace detail {

  struct progress_reporter : public listener
  {
    progress_reporter(const detail::failure_formatter& failure_formatter)
      : specs_run_(0), specs_succeeded_(0), specs_failed_(0), specs_skipped_(0),
      failure_formatter_(failure_formatter)
    {}

	progress_reporter& operator=(const progress_reporter&) { return *this; }

    virtual void test_run_starting()
    {
      specs_run_ = 0;
      specs_succeeded_ = 0;
      specs_failed_ = 0;
      specs_skipped_ = 0;
      failures_.clear();
      contexts_.clear();
    }

    virtual void test_run_complete()
    {
    }

    virtual void context_starting(const char* desc)
    {
      contexts_.push_back(std::string(desc));
    }

    virtual void context_ended(const char*)
    {
      contexts_.pop_back();
    }

    virtual void test_run_error(const char*, const struct test_run_error&)
    {}

    void it_starting(const char*)
    {
      specs_run_++;
    }

    void it_succeeded(const char*)
    {
      specs_succeeded_++;
    }

    void it_failed(const char* desc, const assertion_exception& ex)
    {
      specs_failed_++;

      std::stringstream ss;
      ss << std::endl;
      ss << current_context_name() << " " << desc << ":" << std::endl;
      ss << failure_formatter_.format(ex);

      failures_.push_back(ss.str());
    }

    void it_unknown_error(const char* desc)
    {
      specs_failed_++;

      std::stringstream ss;
      ss << std::endl;
      ss << current_context_name() << " " << desc << ":" << std::endl;
      ss << "Unknown exception";
      ss << std::endl;

      failures_.push_back(ss.str());
    }

    void it_skip(const char* /* desc */)
    {
      specs_skipped_++;
    }

    void it_list(const char* /* desc */)
    {
    }

    bool did_we_pass() const
    {
      return specs_run_ > 0 && specs_failed_ == 0 && test_run_errors_.size() == 0;
    }

    protected:
    std::string current_context_name()
    {
      std::string name;

      std::for_each(contexts_.begin(), contexts_.end(), [&](const std::string context){
          if(name.size() > 0)
          {
          name += " ";
          }

          name += context;
          });

      return name;
    }

    protected:
    int specs_run_;
    int specs_succeeded_;
    int specs_failed_;
    int specs_skipped_;
    const detail::failure_formatter& failure_formatter_;
    std::list<std::string> contexts_;
    std::list<std::string> failures_;
    std::list<std::string> test_run_errors_;
  };
}}

#endif
