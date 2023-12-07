#ifndef BANDIT_RUNNER_H
#define BANDIT_RUNNER_H

#include <bandit/registration/registrar.h>
#include <bandit/options.h>
#include <bandit/controller.h>
#include <bandit/version.h>

#include <bandit/colorizers.h>
#include <bandit/failure_formatters.h>
#include <bandit/reporters.h>
#include <bandit/run_policies.h>

namespace bandit {
  inline int run(const detail::options& opt, const detail::spec_registry& specs,
      detail::controller_t& controller = detail::registered_controller()) {
    if (opt.help() || !opt.parsed_ok()) {
      opt.print_usage();
      return !opt.parsed_ok();
    }

    if (opt.version()) {
      std::cout << "bandit version " << BANDIT_VERSION << std::endl;
      return 0;
    }

    controller.get_reporter().test_run_starting();

    bool hard_skip = false;
    context::bandit global_context("", hard_skip);
    controller.get_contexts().push_back(&global_context);

    for (auto func : specs) {
      func();
    };

    controller.get_reporter().test_run_complete();

    return controller.get_reporter().did_we_pass() ? 0 : 1;
  }

  inline void use_default_colorizers(detail::choice_options& choices) {
    choices.colorizers.add("off", [&](detail::controller_t& controller) {
      controller.set_colorizer(new colorizer::off());
    });
    choices.colorizers.add("dark", [&](detail::controller_t& controller) {
      controller.set_colorizer(new colorizer::dark());
    });
    choices.colorizers.add("light", [&](detail::controller_t& controller) {
      controller.set_colorizer(new colorizer::light());
    }, true);
  }

  inline void use_default_formatters(detail::choice_options& choices) {
    choices.formatters.add("qt", [&](detail::controller_t& controller) {
      controller.set_formatter(new failure_formatter::qt_creator());
    });
    choices.formatters.add("vs", [&](detail::controller_t& controller) {
      controller.set_formatter(new failure_formatter::visual_studio());
    });
    choices.formatters.add("posix", [&](detail::controller_t& controller) {
      controller.set_formatter(new failure_formatter::posix());
    }, true);
  }

  inline void use_default_reporters(detail::choice_options& choices) {
    choices.reporters.add("singleline", [&](detail::controller_t& controller) {
      controller.set_reporter(new bandit::reporter::singleline(controller.get_formatter(), controller.get_colorizer()));
    });
    choices.reporters.add("xunit", [&](detail::controller_t& controller) {
      controller.set_reporter(new bandit::reporter::xunit(controller.get_formatter(), controller.get_report_timing()));
    });
    choices.reporters.add("info", [&](detail::controller_t& controller) {
      controller.set_reporter(new bandit::reporter::info(controller.get_formatter(), controller.get_colorizer(), controller.get_report_timing()));
    }, true);
    choices.reporters.add("spec", [&](detail::controller_t& controller) {
      controller.set_reporter(new bandit::reporter::spec(controller.get_formatter(), controller.get_colorizer()));
    });
    choices.reporters.add("crash", [&](detail::controller_t& controller) {
      controller.set_reporter(new bandit::reporter::crash(controller.get_formatter()));
    });
    choices.reporters.add("dots", [&](detail::controller_t& controller) {
      controller.set_reporter(new bandit::reporter::dots(controller.get_formatter(), controller.get_colorizer()));
    });
  }

  inline void use_defaults(detail::choice_options& choices) {
    use_default_colorizers(choices);
    use_default_formatters(choices);
    use_default_reporters(choices);
  }

  inline int run(int argc, char* argv[], const detail::choice_options& choices, bool allow_further = true) {
    detail::options opt(argc, argv, choices);

    if (!allow_further &&
        (opt.has_further_arguments() || opt.has_unknown_options())) {
      opt.print_usage();
      return 1;
    }

    detail::controller_t controller;

    controller.set_report_timing(opt.report_timing());

    if (!opt.update_controller_settings(controller)) {
      return 1;
    }

    controller.set_policy(new run_policy::bandit(opt.filter_chain(), opt.break_on_failure(), opt.dry_run()));

    detail::register_controller(&controller);
    return run(opt, detail::specs());
  }

  inline int run(int argc, char* argv[], bool allow_further = true) {
    detail::choice_options choices;

    use_defaults(choices);

    return run(argc, argv, choices, allow_further);
  }
}
#endif
