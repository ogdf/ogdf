#ifndef BANDIT_ASSERTION_EXCEPTION_H
#define BANDIT_ASSERTION_EXCEPTION_H

namespace bandit { namespace detail {

  struct assertion_exception : public std::runtime_error
  {
    assertion_exception(const std::string& message,
        const std::string& filename, const unsigned int linenumber)
      : std::runtime_error(message), file_name_(filename), line_number_(linenumber)
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
