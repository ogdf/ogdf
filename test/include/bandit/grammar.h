#ifndef BANDIT_GRAMMAR_H
#define BANDIT_GRAMMAR_H

namespace bandit {

  inline void describe(const char* desc, detail::voidfunc_t func,
      detail::listener& listener, detail::contextstack_t& context_stack,
      bool hard_skip = false)
  {
    listener.context_starting(desc);

    context_stack.back()->execution_is_starting();

    detail::bandit_context ctxt(desc, hard_skip);

    context_stack.push_back(&ctxt);
    try
    {
      func();
    }
    catch(const bandit::detail::test_run_error& error)
    {
      listener.test_run_error(desc, error);
    }

    context_stack.pop_back();

    listener.context_ended(desc);
  }

  inline void describe(const char* desc, detail::voidfunc_t func)
  {
    describe(desc, func, detail::registered_listener(), detail::context_stack());
  }

  inline void describe_skip(const char* desc, detail::voidfunc_t func,
      detail::listener& listener, detail::contextstack_t& context_stack)
  {
    bool skip = true;
    describe(desc, func, listener, context_stack, skip);
  }

  inline void describe_skip(const char* desc, detail::voidfunc_t func)
  {
    describe_skip(desc, func, detail::registered_listener(),
        detail::context_stack());
  }

  inline void before_each(detail::voidfunc_t func,
      detail::contextstack_t& context_stack)
  {
    context_stack.back()->register_before_each(func);
  }

  inline void before_each(detail::voidfunc_t func)
  {
    before_each(func, detail::context_stack());
  }

  inline void after_each(detail::voidfunc_t func,
      detail::contextstack_t& context_stack)
  {
    context_stack.back()->register_after_each(func);
  }

  inline void after_each(detail::voidfunc_t func)
  {
    after_each(func, detail::context_stack());
  }

  inline void it_skip(const char* desc, detail::voidfunc_t, detail::listener& listener)
  {
    listener.it_skip(desc);
  }

  inline void it_skip(const char* desc, detail::voidfunc_t func)
  {
    it_skip(desc, func, detail::registered_listener());
  }

  inline void it(const char* desc, detail::voidfunc_t func, detail::listener& listener,
      detail::contextstack_t& context_stack,
      bandit::adapters::assertion_adapter& assertion_adapter,
      const detail::run_policy& run_policy)
  {
    if(!run_policy.should_run(desc, context_stack))
    {
      it_skip(desc, func, listener);
      return;
    }

    listener.it_starting(desc);

    context_stack.back()->execution_is_starting();

    auto run_before_eaches = [&](){
      for_each(context_stack.begin(), context_stack.end(), [](detail::context* ctxt){
          ctxt->run_before_eaches();
      });
    };

    auto run_after_eaches = [&](){
      for_each(context_stack.begin(), context_stack.end(), [](detail::context* ctxt){
          ctxt->run_after_eaches();
      });
    };

    try
    {
      assertion_adapter.adapt_exceptions([&](){
          run_before_eaches();

          func();
          listener.it_succeeded(desc);
      });
    }
    catch(const bandit::detail::assertion_exception& ex)
    {
      listener.it_failed(desc, ex);
    }
    catch(...)
    {
      listener.it_unknown_error(desc);
    }

    try
    {
      run_after_eaches();
    }
    catch(...)
    {
      listener.it_unknown_error(desc);
    }
  }

  inline void it(const char* desc, detail::voidfunc_t func)
  {
    it(desc, func, detail::registered_listener(), detail::context_stack(),
        detail::registered_adapter(), detail::registered_run_policy());
  }


}

#endif
