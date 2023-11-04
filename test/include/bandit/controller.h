#ifndef BANDIT_CONTROLLER_H
#define BANDIT_CONTROLLER_H

#include <memory>
#include <bandit/adapters.h>
#include <bandit/colorizers/interface.h>
#include <bandit/failure_formatters/interface.h>
#include <bandit/reporters/interface.h>
#include <bandit/run_policies/interface.h>

namespace bandit {
  namespace detail {
    struct controller_t {
      controller_t() : adapter(new adapter::snowhouse), report_timing(false) {}

      context::stack_t& get_contexts() {
        return context_stack;
      }

      assertion_adapter_t& get_adapter() {
        throw_if_nullptr(adapter.get(), "assertion adapter", "bandit::detail::controller_t::set_adapter()");
        return *adapter;
      }

      void set_adapter(assertion_adapter_t* adapter_) {
        adapter.reset(adapter_);
      }

      colorizer_t& get_colorizer() {
        throw_if_nullptr(colorizer.get(), "colorizer", "bandit::detail::controller_t::set_colorizer()");
        return *colorizer;
      }

      void set_colorizer(colorizer_t* colorizer_) {
        colorizer.reset(colorizer_);
      }

      failure_formatter_t& get_formatter() {
        throw_if_nullptr(formatter.get(), "formatter", "bandit::detail::controller_t::set_formatter()");
        return *formatter;
      }

      void set_formatter(failure_formatter_t* formatter_) {
        formatter.reset(formatter_);
      }

      reporter_t& get_reporter() {
        throw_if_nullptr(reporter.get(), "reporter", "bandit::detail::controller_t::set_reporter()");
        return *reporter;
      }

      void set_reporter(reporter_t* reporter_) {
        reporter.reset(reporter_);
      }

      run_policy_t& get_policy() {
        throw_if_nullptr(run_policy.get(), "run policy", "bandit::detail::controller_t::set_policy()");
        return *run_policy;
      }

      void set_policy(run_policy_t* run_policy_) {
        run_policy.reset(run_policy_);
      }

      bool get_report_timing() {
        return report_timing;
      }

      void set_report_timing(bool report_timing_) {
        report_timing = report_timing_;
      }

      void describe(const std::string& desc, const std::function<void()>& func, bool hard_skip) {
        reporter->context_starting(desc);

        context_stack.back()->execution_is_starting();

        context::bandit ctxt(desc, hard_skip);
        context_stack.push_back(&ctxt);

        try {
          func();
        } catch (const bandit::detail::test_run_error& error) {
          reporter->test_run_error(desc, error);
        }

        context_stack.pop_back();

        reporter->context_ended(desc);
      }

      void before_each(const std::function<void()>& func) {
        context_stack.throw_if_empty("before_each");
        context_stack.back()->register_before_each(func);
      }

      void after_each(const std::function<void()>& func) {
        context_stack.throw_if_empty("after_each");
        context_stack.back()->register_after_each(func);
      }

      void it(const std::string& desc, const std::function<void()>& func, bool hard_skip) {
        context_stack.throw_if_empty("it");
        if (hard_skip || !run_policy->should_run(desc, context_stack)) {
          reporter->it_skip(desc);
        } else {
          reporter->it_starting(desc);
          context_stack.back()->execution_is_starting();

          bool success = false;
          context::interface* last_successful_before_each_context = nullptr;
          try_with_adapter(desc, true, [&] {
            for (auto context : context_stack) {
              context->run_before_eaches();
              last_successful_before_each_context = context;
            }

            func();
            success = true;
          });

          try_with_adapter(desc, success, [&] {
            bool do_run_after_each = false;
            std::for_each(context_stack.rbegin(), context_stack.rend(), [&](context::interface* context) {
              if (context == last_successful_before_each_context) {
                do_run_after_each = true;
              }
              if (do_run_after_each) {
                context->run_after_eaches();
              }
            });

            if (success) {
              reporter->it_succeeded(desc);
            }
          });
        }
      }

      // A function is required to initialize a static controller variable in a header file
      // and this struct aims at encapsulating this function
      static void register_controller(controller_t* controller) {
        if (controller == nullptr) {
          throw std::runtime_error("Invalid null controller passed to bandit::detail::register_controller()");
        }
        get_controller_address() = controller;
      }

      static controller_t& registered_controller() {
        auto controller = get_controller_address();
        throw_if_nullptr(controller, "controller", "bandit::detail::register_controller()");
        return *controller;
      }

    private:
      static controller_t*& get_controller_address() {
        static controller_t* controller_ = nullptr;
        return controller_;
      }

      static void throw_if_nullptr(const void* ptr, const std::string& name, const std::string& setter) {
        if (ptr == nullptr) {
          throw std::runtime_error("No " + name + " set. Please set it using " + setter);
        }
      }

      void try_with_adapter(const std::string& desc, bool allow_fail, const std::function<void()>& do_it) {
        if (allow_fail) {
          try {
            adapter->adapt_exceptions([&] { do_it(); });
          } catch (const bandit::detail::assertion_exception& ex) {
            reporter->it_failed(desc, ex);
            run_policy->encountered_failure();
          } catch (const std::exception& ex) {
            std::string err = std::string("exception: ") + ex.what();
            reporter->it_failed(desc, bandit::detail::assertion_exception(err));
            run_policy->encountered_failure();
          } catch (...) {
            reporter->it_unknown_error(desc);
            run_policy->encountered_failure();
          }
        } else {
          try {
            do_it();
          } catch (...) {
            /* ignore */
          }
        }
      }

      context::stack_t context_stack;
      std::unique_ptr<assertion_adapter_t> adapter;
      std::unique_ptr<colorizer_t> colorizer;
      std::unique_ptr<failure_formatter_t> formatter;
      std::unique_ptr<reporter_t> reporter;
      std::unique_ptr<run_policy_t> run_policy;
      bool report_timing;
    };

    inline void register_controller(controller_t* controller) {
      controller_t::register_controller(controller);
    }

    inline controller_t& registered_controller() {
      return controller_t::registered_controller();
    }
  }
}
#endif
