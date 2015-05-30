#ifndef BANDIT_ASSERTION_EXCEPTION_H
#define BANDIT_ASSERTION_EXCEPTION_H

namespace bandit { namespace detail {

  struct assertion_exception : public std::runtime_error
  {
    assertion_exception(const std::string& message,
        const std::string& file_name, const unsigned int line_number)
      : std::runtime_error(message), file_name_(file_name), line_number_(line_number)
    {}

    assertion_exception(const std::string& message)
      : std::runtime_error(message), line_number_(0)
    {}

    //
    // To make gcc < 4.7 happy.
    //
    virtual ~assertion_exception() throw()
    {}

    const std::string& file_name() const
    {
      return file_name_;
    }

    unsigned int line_number() const
    {
      return line_number_;
    }

    private:
    std::string file_name_;
    unsigned int line_number_;
  };
}}

#endif
