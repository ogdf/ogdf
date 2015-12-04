#ifndef BANDIT_CONTEXT_H
#define BANDIT_CONTEXT_H

namespace bandit {
  namespace detail {

    class context
    {
      public:
        virtual ~context() {}
        virtual const std::string& name() = 0;
        virtual void execution_is_starting() = 0;
        virtual void register_before_each(voidfunc_t func) = 0;
        virtual void register_after_each(voidfunc_t func) = 0;
        virtual void run_before_eaches() = 0;
        virtual void run_after_eaches() = 0;
        virtual bool hard_skip() = 0;
        virtual bool list_tests() = 0;
    };

    class bandit_context : public context
    {
      public:
        bandit_context(const char* desc, bool hard_skip, bool list_tests)
          : desc_(desc), hard_skip_(hard_skip), list_tests_(list_tests), is_executing_(false)
        {}

        const std::string& name()
        {
          return desc_;
        }

        void execution_is_starting()
        {
          is_executing_ = true;
        }

        void register_before_each(voidfunc_t func)
        {
          if(is_executing_)
          {
            throw test_run_error("before_each was called after 'describe' or 'it'");
          }

          before_eaches_.push_back(func);
        }

        void register_after_each(voidfunc_t func)
        {
          if(is_executing_)
          {
            throw test_run_error("after_each was called after 'describe' or 'it'");
          }

          after_eaches_.push_back(func);
        }

        void run_before_eaches()
        {
          run_all(before_eaches_);
        }

        void run_after_eaches()
        {
          run_all(after_eaches_);
        }

        bool hard_skip()
        {
          return hard_skip_;
        }

        bool list_tests()
        {
          return list_tests_;
        }

      private:
        void run_all(const std::list<voidfunc_t>& funcs)
        {
          auto call_func = [](voidfunc_t f){ f(); };

          for_each(funcs.begin(), funcs.end(), call_func);
        }

      private:
        std::string desc_;
        bool hard_skip_;
        bool list_tests_;
        bool is_executing_;
        std::list<voidfunc_t> before_eaches_;
        std::list<voidfunc_t> after_eaches_;
    };
    typedef std::deque<context*> contextstack_t;

    inline contextstack_t& context_stack()
    {
      static contextstack_t contexts;
      return contexts;
    }
  }
}

#endif
