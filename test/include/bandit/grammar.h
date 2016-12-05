#ifndef BANDIT_GRAMMAR_H
#define BANDIT_GRAMMAR_H

namespace bandit {

  inline bool has_context_with_list_tests(const detail::contextstack_t& contexts)
  {
    detail::contextstack_t::const_iterator it;
    for(it = contexts.begin(); it != contexts.end(); it++)
    {
      if((*it)->list_tests())
      {
        return true;
      }
    }

    return false;
  }

  inline void describe(const std::string &desc, detail::voidfunc_t func,
      detail::listener& listener, detail::contextstack_t& context_stack,
      bool hard_skip = false, bool list_tests = false)
  {
    listener.context_starting(desc.c_str());

    context_stack.back()->execution_is_starting();

    detail::bandit_context ctxt(desc.c_str(), hard_skip, list_tests);

    context_stack.push_back(&ctxt);
    try
    {
      func();
    }
    catch(const bandit::detail::test_run_error& error)
    {
      listener.test_run_error(desc.c_str(), error);
    }

    context_stack.pop_back();

    listener.context_ended(desc.c_str());
  }

  inline void describe(const std::string &desc, detail::voidfunc_t func, bool skip = false)
  {
    describe(desc, func, detail::registered_listener(), detail::context_stack(), skip);
  }

  inline void describe_skip(const std::string &desc, detail::voidfunc_t func,
      detail::listener& listener, detail::contextstack_t& context_stack)
  {
    describe(desc, func, listener, context_stack, true);
  }

  inline void describe_skip(const std::string &desc, detail::voidfunc_t func)
  {
    describe_skip(desc, func, detail::registered_listener(), detail::context_stack());
  }

  inline void before_each(detail::voidfunc_t func, detail::contextstack_t& context_stack)
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

  inline void it_skip(const std::string &desc, detail::voidfunc_t, detail::listener& listener)
  {
    listener.it_skip(desc.c_str());
  }

  inline void it_skip(const std::string &desc, detail::voidfunc_t func)
  {
    it_skip(desc, func, detail::registered_listener());
  }

  inline void it_list(const std::string &desc, detail::voidfunc_t, detail::listener& listener)
  {
    listener.it_list(desc.c_str());
  }

  inline void it_list(const std::string &desc, detail::voidfunc_t func)
  {
    it_list(desc, func, detail::registered_listener());
  }

  inline void it(const std::string &desc, detail::voidfunc_t func, detail::listener& listener,
      detail::contextstack_t& context_stack,
      bandit::adapters::assertion_adapter& assertion_adapter,
      detail::run_policy& run_policy)
  {
    if(has_context_with_list_tests(context_stack))
    {
      it_list(desc, func, listener);
      return;
    }

    if(!run_policy.should_run(desc.c_str(), context_stack))
    {
      it_skip(desc, func, listener);
      return;
    }

    listener.it_starting(desc.c_str());

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

    bool we_have_been_successful_so_far = false;
    try
    {
      assertion_adapter.adapt_exceptions([&](){
          run_before_eaches();

          func();
          we_have_been_successful_so_far = true;
      });
    }
    catch(const bandit::detail::assertion_exception& ex)
    {
      listener.it_failed(desc.c_str(), ex);
      run_policy.encountered_failure();
    }
    catch(const std::exception& ex)
    {
      std::string err = std::string("exception: ") + ex.what();
      listener.it_failed(desc.c_str(), bandit::detail::assertion_exception(err));
      run_policy.encountered_failure();
    }
    catch(...)
    {
      listener.it_unknown_error(desc.c_str());
      run_policy.encountered_failure();
    }

    try
    {
      assertion_adapter.adapt_exceptions([&](){
          run_after_eaches();

          if(we_have_been_successful_so_far)
          {
            listener.it_succeeded(desc.c_str());
          }
      });
    }
    catch(const bandit::detail::assertion_exception& ex)
    {
      listener.it_failed(desc.c_str(), ex);
      run_policy.encountered_failure();
    }
    catch(const std::exception& ex)
    {
      std::string err = std::string("exception: ") + ex.what();
      listener.it_failed(desc.c_str(), bandit::detail::assertion_exception(err));
      run_policy.encountered_failure();
    }
    catch(...)
    {
      listener.it_unknown_error(desc.c_str());
      run_policy.encountered_failure();
    }
  }

  inline void it(const std::string &desc, detail::voidfunc_t func, bool skip = false)
  {
    if (skip) {
      it_skip(desc, func);
    } else {
      it(desc, func, detail::registered_listener(), detail::context_stack(),
         detail::registered_adapter(), detail::registered_run_policy());
    }
  }
}

#endif
