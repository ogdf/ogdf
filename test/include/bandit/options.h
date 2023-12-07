#ifndef BANDIT_OPTIONS_H
#define BANDIT_OPTIONS_H

#include <algorithm>
#include <functional>
#include <map>
#include <vector>
#include <iostream>

#include <bandit/external/optionparser.h>
#include <bandit/run_policies/filter_chain.h>
#include <bandit/controller.h>

namespace bandit {
  namespace detail {
    using controller_func_t = std::function<void(detail::controller_t&)>;

    struct option_map {
      void add(const std::string& name, controller_func_t func, bool is_default = false) {
        auto it = map.emplace(std::make_pair(name, func));
        if (is_default || map.size() == 1) {
          default_it = it.first;
        }
      }

      std::string comma_separated_list() const {
        std::string csl;
        auto it = map.begin();
        if (it != map.end()) {
          csl += it->first;
          std::for_each(++it, map.end(), [&](const std::pair<std::string, controller_func_t>& val) {
            csl += ", " + val.first;
          });
        } else {
          csl = "[no registered entities]";
        }
        return csl;
      }

      using map_type = std::map<std::string, controller_func_t>;
      map_type map;
      map_type::iterator default_it;
    };

    struct choice_options {
      option_map colorizers;
      option_map formatters;
      option_map reporters;
    };

    struct options {
      struct argument : public option::Arg {
        static std::string name(const option::Option& option) {
          std::string copy(option.name);
          return copy.substr(0, option.namelen);
        }

        static option::ArgStatus Required(const option::Option& option, bool msg) {
          if (option.arg != nullptr) {
            return option::ARG_OK;
          }
          if (msg) {
            std::cerr << "Option '" << name(option) << "' requires an argument\n";
          }
          return option::ARG_ILLEGAL;
        }
      };

      static const choice_options& empty_choice_options() {
        static choice_options choices;
        return choices;
      }

      options(int argc, char* argv[], const choice_options& choices = empty_choice_options())
          : choices_(choices),
            usage_{
                {UNKNOWN, 0, "", "", argument::None,
                    "USAGE: <executable> [options]\n\n"
                    "Options:"},
                {VERSION, 0, "", "version", argument::None,
                    "  --version, \tPrint version of bandit"},
                {HELP, 0, "", "help", argument::None,
                    "  --help, \tPrint usage and exit."},
                {SKIP, 0, "", "skip", argument::Required,
                    "  --skip=<substring>, "
                    "\tSkip all 'describe' and 'it' containing substring"},
                {ONLY, 0, "", "only", argument::Required,
                    "  --only=<substring>, "
                    "\tRun only 'describe' and 'it' containing substring"},
                {BREAK_ON_FAILURE, 0, "", "break-on-failure", argument::None,
                    "  --break-on-failure, "
                    "\tStop test run on first failing test"},
                {DRY_RUN, 0, "", "dry-run", argument::None,
                    "  --dry-run, "
                    "\tSkip all tests. Use to list available tests"},
                {REPORT_TIMING, 0, "", "report-timing", argument::None,
                    "  --report-timing, "
                    "\tInstruct reporter to report timing information"},
            },
            reporter_help_("  --reporter=<reporter>, "
                "\tSelect reporter: " + choices_.reporters.comma_separated_list()),
            colorizer_help_("  --colorizer=<colorizer>, "
                "\tSelect color theme: " + choices_.colorizers.comma_separated_list()),
            formatter_help_("  --formatter=<formatter>, "
                "\tSelect error formatter: " + choices_.formatters.comma_separated_list()) {
        usage_.push_back(option::Descriptor{REPORTER, 0, "", "reporter", argument::Required, reporter_help_.c_str()});
        usage_.push_back(option::Descriptor{COLORIZER, 0, "", "colorizer", argument::Required, colorizer_help_.c_str()});
        usage_.push_back(option::Descriptor{FORMATTER, 0, "", "formatter", argument::Required, formatter_help_.c_str()});
        usage_.push_back(option::Descriptor{0, 0, nullptr, nullptr, nullptr, nullptr});

        argc -= (argc > 0);
        argv += (argc > 0); // Skip program name (argv[0]) if present
        option::Stats stats(usage(), argc, argv);
        options_.resize(stats.options_max);
        std::vector<option::Option> buffer(stats.buffer_max);
        option::Parser parse(usage(), argc, argv, options_.data(), buffer.data());
        parsed_ok_ = !parse.error();
        has_further_arguments_ = (parse.nonOptionsCount() != 0);
        has_unknown_options_ = (options_[UNKNOWN] != nullptr);

        for (int i = 0; i < parse.optionsCount(); ++i) {
          option::Option& opt = buffer[i];
          switch (opt.index()) {
          case SKIP:
            filter_chain_.push_back({opt.arg, true});
            break;
          case ONLY:
            filter_chain_.push_back({opt.arg, false});
            break;
          }
        }
      }

      bool help() const {
        return options_[HELP] != nullptr;
      }

      bool parsed_ok() const {
        return parsed_ok_;
      }

      bool has_further_arguments() const {
        return has_further_arguments_;
      }

      bool has_unknown_options() const {
        return has_unknown_options_;
      }

      void print_usage() const {
        option::printUsage(std::cout, usage());
      }

      bool version() const {
        return options_[VERSION] != nullptr;
      }

      bool update_controller_settings(controller_t& controller) {
        struct chooser_t {
          const char* name;
          const option_map& choice;
          option_index index;
        };
        for (auto chooser : {
                 chooser_t{"colorizer", choices_.colorizers, COLORIZER},
                 chooser_t{"formatter", choices_.formatters, FORMATTER},
                 chooser_t{"reporter", choices_.reporters, REPORTER},
             }) {
          if (chooser.choice.map.empty()) {
            throw std::runtime_error(std::string("No ") + chooser.name + " set.");
          }
          if (!apply(controller, chooser.choice, chooser.index)) {
            print_usage(chooser.name, chooser.choice, chooser.index);
            parsed_ok_ = false;
            return false;
          }
        }
        return true;
      }

      const run_policy::filter_chain_t& filter_chain() const {
        return filter_chain_;
      }

      bool break_on_failure() const {
        return options_[BREAK_ON_FAILURE] != nullptr;
      }

      bool dry_run() const {
        return options_[DRY_RUN] != nullptr;
      }

      bool report_timing() const {
        return options_[REPORT_TIMING] != nullptr;
      }

    private:
      enum option_index {
        UNKNOWN,
        VERSION,
        HELP,
        REPORTER,
        COLORIZER,
        FORMATTER,
        SKIP,
        ONLY,
        BREAK_ON_FAILURE,
        DRY_RUN,
        REPORT_TIMING,
      };

      bool apply(controller_t& controller, const option_map& choice, enum option_index index) const {
        option_map::map_type::const_iterator it{choice.default_it};
        if (options_[index] != nullptr) {
          it = choice.map.find(std::string(options_[index].arg));
          if (it == choice.map.end()) {
            return false;
          }
        }
        it->second(controller);
        return true;
      }

      void print_usage(std::string&& what, const option_map& choice, option_index index) const {
        std::cerr << "Unknown " << what << " '" << options_[index].arg << "'." << std::endl;
        std::cerr << "Option argument of '--" << what << "' must be one of: " << choice.comma_separated_list() << std::endl;
      }

      const option::Descriptor* usage() const {
        return &usage_[0];
      }

      const choice_options& choices_;

      std::vector<option::Descriptor> usage_;

      const std::string reporter_help_;
      const std::string colorizer_help_;
      const std::string formatter_help_;

      std::vector<option::Option> options_;
      run_policy::filter_chain_t filter_chain_;
      bool parsed_ok_;
      bool has_further_arguments_;
      bool has_unknown_options_;
    };
  }
}
#endif
