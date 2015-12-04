#ifndef BANDIT_OPTIONS_H
#define BANDIT_OPTIONS_H

namespace bandit { namespace detail {

    // TODO: print any unknown options
    // TODO: check for parser errors
    struct options
    {

      options(int argc, char* argv[])
      {
        argc -= (argc>0); argv += (argc>0); // Skip program name (argv[0]) if present
        option::Stats stats(usage(), argc, argv);
        options_.resize(stats.options_max);
        std::vector<option::Option> buffer(stats.buffer_max);
        option::Parser parse(usage(), argc, argv, options_.data(), buffer.data());
        parsed_ok_ = !parse.error();
      }

      bool help() const
      {
        return options_[HELP] != nullptr;
      }

      void print_usage() const
      {
        option::printUsage(std::cout, usage());
      }

      bool version() const
      {
        return options_[VERSION] != nullptr;
      }

      bool list_tests() const
      {
        return options_[LIST_TESTS] != nullptr;
      }

      const char* reporter() const
      {
        return options_[REPORTER].arg;
      }

      bool no_color() const
      {
        return options_[NO_COLOR] != nullptr;
      }

      typedef enum
      {
        FORMATTER_DEFAULT,
        FORMATTER_VS,
        FORMATTER_UNKNOWN
      } formatters;

      formatters formatter() const
      {
        std::string arg = options_[FORMATTER].arg ? options_[FORMATTER].arg : "";
        if(arg == "vs")
        {
          return formatters::FORMATTER_VS;
        }

        return formatters::FORMATTER_DEFAULT;
      }

      const char* skip() const
      {
        return options_[SKIP].arg ? options_[SKIP].arg : "";
      }

      const char* only() const
      {
        return options_[ONLY].arg ? options_[ONLY].arg : "";
      }

      bool break_on_failure() const
      {
          return options_[BREAK_ON_FAILURE] != nullptr;
      }

      private:
        enum option_index { UNKNOWN, VERSION, HELP, REPORTER, NO_COLOR,
          FORMATTER, SKIP, ONLY, BREAK_ON_FAILURE, LIST_TESTS };

        static const option::Descriptor* usage()
        {
          static const option::Descriptor usage[] =
          {
            {UNKNOWN, 0, "", "", option::Arg::None,         "USAGE: <executable> [options]\n\n"
                                                            "Options:" },
            {VERSION, 0, "", "version", option::Arg::None, "  --version, \tPrint version of bandit"},
            {HELP, 0, "", "help", option::Arg::None,        "  --help, \tPrint usage and exit."},
            {REPORTER, 0, "", "reporter", option::Arg::Optional, "  --reporter=<reporter>, \tSelect reporter (dots, singleline, xunit, info, spec)"},
            {NO_COLOR, 0, "", "no-color", option::Arg::None,     "  --no-color, \tSuppress colors in output"},
            {FORMATTER, 0, "", "formatter", option::Arg::Optional, "  --formatter=<formatter>, \tSelect formatting of errors (default, vs)"},
            {SKIP, 0, "", "skip", option::Arg::Optional, "  --skip=<substring>, \tskip all 'describe' and 'it' containing substring"},
            {ONLY, 0, "", "only", option::Arg::Optional, "  --only=<substring>, \tonly run 'describe' and 'it' containing substring"},
            {BREAK_ON_FAILURE, 0, "", "break-on-failure", option::Arg::Optional, "  --break-on-failure, \tstop test run on first failing test"},
            {LIST_TESTS, 0, "", "list-tests", option::Arg::None, "  --list-tests, \tLists all tests."},
            {0, 0, 0, 0, 0, 0}
          };

          return usage;
        }

      private:
        std::vector<option::Option> options_;
        bool parsed_ok_;

    };

}}

#endif
