#ifndef BANDIT_BANDIT_RUN_POLICY_H
#define BANDIT_BANDIT_RUN_POLICY_H

namespace bandit { namespace detail {

  struct bandit_run_policy : public run_policy
  {
    bandit_run_policy(const char* skip_pattern, const char* only_pattern)
      : skip_pattern_(skip_pattern), only_pattern_(only_pattern)
    {}

    bool should_run(const char* it_name, const contextstack_t& contexts) const
    {
      //
      // Never run if a context has been marked as skip
      // using 'describe_skip'
      //
      if(has_context_with_hard_skip(contexts))
      {
        return false;
      }

      //
      // Always run if no patterns have been specifed
      //
      if(!has_skip_pattern() && !has_only_pattern())
      {
        return true;
      }

      if(has_only_pattern() && !has_skip_pattern())
      {
        return context_matches_only_pattern(contexts)
          || matches_only_pattern(it_name);
      }

      if(has_skip_pattern() && !has_only_pattern())
      {
        bool skip = context_matches_skip_pattern(contexts) ||
          matches_skip_pattern(it_name);
        return !skip;
      }

      //
      // If we've come this far, both 'skip' and 'only'
      // have been specified.
      //
      // If our contexts match 'only' we're still good
      // regardless of whether there's a 'skip' somewhere
      // in the context stack as well.
      if(context_matches_only_pattern(contexts))
      {
        //
        // We can still mark the current 'it' as 'skip'
        // and ignore it. We check that here.
        //
        return !matches_skip_pattern(it_name);
      }

      //
      // If we've gotten this far, the context matches 'skip'
      // We can still run this spec if it is specifically marked
      // as 'only'.
      //
      return matches_only_pattern(it_name);
    }

    private:
    bool has_context_with_hard_skip(const contextstack_t& contexts) const
    {
      contextstack_t::const_iterator it;
      for(it = contexts.begin(); it != contexts.end(); it++)
      {
        if((*it)->hard_skip())
        {
          return true;
        }
      }

      return false;
    }

    bool has_only_pattern() const
    {
      return only_pattern_.size() > 0;
    }

    bool has_skip_pattern() const
    {
      return skip_pattern_.size() > 0;
    }

    bool context_matches_only_pattern(const contextstack_t& contexts) const
    {
      contextstack_t::const_iterator it;
      for(it = contexts.begin(); it != contexts.end(); it++)
      {
        if(matches_only_pattern((*it)->name()))
        {
          return true;
        }
      }

      return false;
    }

    bool context_matches_skip_pattern(const contextstack_t& contexts) const
    {
      contextstack_t::const_iterator it;
      for(it = contexts.begin(); it != contexts.end(); it++)
      {
        if(matches_skip_pattern((*it)->name()))
        {
          return true;
        }
      }

      return false;
    }

    bool matches_only_pattern(const char* name) const
    {
      std::string n(name);
      return matches_only_pattern(n);
    }

    bool matches_only_pattern(const std::string& name) const
    {
      return matches_pattern(name, only_pattern_);
    }

    bool matches_skip_pattern(const char* name) const
    {
      std::string n(name);
      return matches_skip_pattern(n);
    }

    bool matches_skip_pattern(const std::string& name) const
    {
      return matches_pattern(name, skip_pattern_);
    }

    bool matches_pattern(const std::string& name, const std::string& pattern) const
    {
      return name.find(pattern) != std::string::npos;
    }

    private:
    std::string skip_pattern_;
    std::string only_pattern_;
  };

}}

#endif
