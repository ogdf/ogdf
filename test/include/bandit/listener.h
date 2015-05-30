#ifndef BANDIT_LISTENER_H
#define BANDIT_LISTENER_H

namespace bandit { namespace detail {
  struct listener
  {
    virtual ~listener() {}

    virtual void test_run_starting() = 0;
    virtual void test_run_complete() = 0;

    virtual void context_starting(const char* desc) = 0;
    virtual void context_ended(const char* desc) = 0;
    virtual void test_run_error(const char* desc, const test_run_error& error) = 0;

    virtual void it_starting(const char* desc) = 0;
    virtual void it_succeeded(const char* desc) = 0;
    virtual void it_failed(const char* desc, const detail::assertion_exception& ex) = 0;
    virtual void it_unknown_error(const char* desc) = 0;
    virtual void it_skip(const char* desc) = 0;

    virtual bool did_we_pass() const = 0;
  };
  typedef std::unique_ptr<listener> listener_ptr;
}}

#endif
