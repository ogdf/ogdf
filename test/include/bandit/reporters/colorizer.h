#ifndef BANDIT_REPORTERS_COLORIZER_H
#define BANDIT_REPORTERS_COLORIZER_H

#ifdef _WIN32
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif

  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#endif

namespace bandit { namespace detail {

#ifdef _WIN32
  struct colorizer
  {
    colorizer(bool colors_enabled = true)
      : colors_enabled_(colors_enabled),
		stdout_handle_(GetStdHandle(STD_OUTPUT_HANDLE))
    {
		original_color_ = get_console_color();
	}

    const char* green() const
    {
      if(colors_enabled_)
	  {
		  set_console_color(FOREGROUND_GREEN);
	  }
	  return "";
    }

    const char* yellow() const
    {
      if(colors_enabled_)
	  {
		  set_console_color(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
	  }
	  return "";
    }

    const char* blue() const
    {
      if(colors_enabled_)
	  {
		  set_console_color(FOREGROUND_BLUE);
	  }
	  return "";
    }

    const char* red() const
    {
      if(colors_enabled_)
	  {
		  set_console_color(FOREGROUND_RED);
	  }
	  return "";
    }

    const char* white() const
    {
      if(colors_enabled_)
	  {
		  set_console_color(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
	  }
	  return "";
    }

    const char* reset() const
    {
      if(colors_enabled_)
	  {
		  set_console_color(original_color_);
	  }
	  return "";
    }

  private:
	  WORD get_console_color() const
	  {
		  CONSOLE_SCREEN_BUFFER_INFO info{};
		  GetConsoleScreenBufferInfo(stdout_handle_, &info);
		  return info.wAttributes;
	  }

	  void set_console_color(WORD color) const
	  {
		  SetConsoleTextAttribute(stdout_handle_, color);
	  }

    private:
    bool colors_enabled_;
	HANDLE stdout_handle_;
	WORD original_color_;
  };

#else
  struct colorizer
  {
    colorizer(bool colors_enabled = true)
      : colors_enabled_(colors_enabled)
    {}

    const char* green() const
    {
      return colors_enabled_ ? "\033[1;32m" : "";
    }

    const char* yellow() const
    {
      return colors_enabled_ ? "\033[1;33m" : "";
    }

    const char* blue() const
    {
      return colors_enabled_ ? "\033[1;34m" : "";
    }

    const char* red() const
    {
      return colors_enabled_ ? "\033[1;31m" : "";
    }

    const char* white() const
    {
      return colors_enabled_ ? "\033[1;37m" : "";
    }

    const char* reset() const
    {
      return colors_enabled_ ? "\033[0m" : "";
    }

    private:
    bool colors_enabled_;
  };
#endif
}}

#endif
