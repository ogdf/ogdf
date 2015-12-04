#ifndef BANDIT_REPORTERS_XUNIT_REPORTER_H
#define BANDIT_REPORTERS_XUNIT_REPORTER_H

namespace bandit { namespace detail {

  struct xunit_reporter : public progress_reporter
  {
    xunit_reporter(std::ostream& stm, const failure_formatter& formatter)
      : progress_reporter(formatter), stm_(stm)
    {
    }

    xunit_reporter(const failure_formatter& formatter)
      : progress_reporter(formatter), stm_(std::cout)
    {
    }

    void it_starting(const char* desc)
    {
      progress_reporter::it_starting(desc);
      work_stm_ << "\t<testcase classname=\"" << escape(current_context_name()) << "\" ";
      work_stm_ << "name=\"" << escape(desc) << "\" time=\"0\">\n";
    }

    void it_succeeded(const char* desc)
    {
      progress_reporter::it_succeeded(desc);
      work_stm_ << "\t</testcase>\n";
    }

    void it_failed(const char* desc, const assertion_exception& ex)
    {
      progress_reporter::it_failed(desc, ex);
      work_stm_ << "\t\t<failure message=\"" << escape(failure_formatter_.format(ex)) << "\" />\n";
      work_stm_ << "\t</testcase>\n";
    }

    void it_unknown_error(const char* desc)
    {
      progress_reporter::it_unknown_error(desc);
      work_stm_ << "\t\t<failure message=\"Unknown exception\" />\n";
      work_stm_ << "\t</testcase>\n";
    }

    void it_skip(const char* desc)
    {
      progress_reporter::it_skip(desc);
      work_stm_ << "\t<testcase classname=\"" << escape(current_context_name()) << "\" ";
      work_stm_ << "name=\"" << escape(desc) << "\" time=\"0\">\n";
      work_stm_ << "\t\t<skipped />\n";
      work_stm_ << "\t</testcase>\n";
    }

    void it_list(const char* desc)
    {
      progress_reporter::it_list(desc);
      work_stm_ << "\t<testcase classname=\"" << escape(current_context_name()) << "\" ";
      work_stm_ << "name=\"" << escape(desc) << "\" time=\"0\"/>\n";
    }

    void test_run_complete()
    {
      stm_ << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
      stm_ << "<testsuite name=\"bandit\" tests=\"" << specs_run_ << "\" errors=\"0\" failures=\""
           << specs_failed_ << "\"";

      if(specs_skipped_ > 0)
      {
        stm_ << " skipped=\"" << specs_skipped_ << "\"";
      }

      stm_ << ">\n";

      stm_ << work_stm_.str();

      stm_ << "</testsuite>\n";
    }

    private:
    std::string escape(const std::string& str)
    {
      std::stringstream stm;

      std::for_each(str.begin(), str.end(), [&](char c){
            switch(c)
            {
              case '&':
                stm << "&amp;";
                break;
              case '<':
                stm << "&lt;";
                break;
              case '>':
                stm << "&gt;";
                break;
              case '\\':
                stm << "&apos;";
                break;
              case '\"':
                stm << "&quot;";
                break;
              default:
                stm << c;
            }
          });

      return stm.str();
    }

    private:
    std::ostream& stm_;
    std::stringstream work_stm_;
  };
}}

#endif
