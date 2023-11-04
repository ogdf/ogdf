#ifndef BANDIT_FAILURE_FORMATTERS_GENERIC_H
#define BANDIT_FAILURE_FORMATTERS_GENERIC_H

#include <sstream>
#include <bandit/failure_formatters/interface.h>

namespace bandit {
  namespace failure_formatter {
    struct generic : public interface {
    protected:
      std::string generic_format(const detail::assertion_exception& err,
          std::string&& before_file,
          std::string&& before_line,
          std::string&& after_line,
          std::string&& after_file,
          std::string&& file_replacement) const {
        std::stringstream ss;
        if (err.file_name().size() > 0) {
          ss << before_file << err.file_name();

          if (err.line_number() > 0) {
            ss << before_line << err.line_number() << after_line;
          }

          ss << after_file;
        } else {
          ss << file_replacement;
        }

        ss << err.what();

        return ss.str();
      }
    };
  }
}
#endif
