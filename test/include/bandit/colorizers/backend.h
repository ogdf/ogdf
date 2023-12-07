#ifndef BANDIT_COLORIZERS_BACKEND_H
#define BANDIT_COLORIZERS_BACKEND_H

#include <string>

#if defined(_WIN32) && !defined(BANDIT_CONFIG_COLOR_ANSI)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

namespace bandit {
  namespace colorizer {
    struct backend {
#if defined(_WIN32) && !defined(BANDIT_CONFIG_COLOR_ANSI)
      backend() : stdout_handle_(GetStdHandle(STD_OUTPUT_HANDLE)) {
        background_color_ = original_color_ = get_console_color();
        background_color_ &= BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED | BACKGROUND_INTENSITY;
      }

      const std::string dark_green() const {
        set_console_color(FOREGROUND_GREEN | background_color_);
        return "";
      }

      const std::string green() const {
        set_console_color(FOREGROUND_GREEN | FOREGROUND_INTENSITY | background_color_);
        return "";
      }

      const std::string brown() const {
        set_console_color(FOREGROUND_RED | FOREGROUND_GREEN | background_color_);
        return "";
      }

      const std::string yellow() const {
        set_console_color(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY | background_color_);
        return "";
      }

      const std::string dark_blue() const {
        set_console_color(FOREGROUND_BLUE | FOREGROUND_INTENSITY | background_color_);
        return "";
      }

      const std::string blue() const {
        set_console_color(FOREGROUND_BLUE | FOREGROUND_INTENSITY | background_color_);
        return "";
      }

      const std::string dark_red() const {
        set_console_color(FOREGROUND_RED | background_color_);
        return "";
      }

      const std::string red() const {
        set_console_color(FOREGROUND_RED | FOREGROUND_INTENSITY | background_color_);
        return "";
      }

      const std::string dark_gray() const {
        set_console_color(FOREGROUND_INTENSITY | background_color_);
        return "";
      }

      const std::string white() const {
        set_console_color(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE |
            FOREGROUND_INTENSITY | background_color_);
        return "";
      }

      const std::string reset() const {
        set_console_color(original_color_);
        return "";
      }

    private:
      WORD get_console_color() const {
        CONSOLE_SCREEN_BUFFER_INFO info{};
        GetConsoleScreenBufferInfo(stdout_handle_, &info);
        return info.wAttributes;
      }

      void set_console_color(WORD color) const {
        SetConsoleTextAttribute(stdout_handle_, color);
      }

    private:
      HANDLE stdout_handle_;
      WORD original_color_;
      WORD background_color_;
#else
      const std::string dark_green() const {
        return "\033[0;32m";
      }

      const std::string green() const {
        return "\033[1;32m";
      }

      const std::string brown() const {
        return "\033[0;33m";
      }

      const std::string yellow() const {
        return "\033[1;33m";
      }

      const std::string dark_blue() const {
        return "\033[0;34m";
      }

      const std::string blue() const {
        return "\033[1;34m";
      }

      const std::string dark_red() const {
        return "\033[0;31m";
      }

      const std::string red() const {
        return "\033[1;31m";
      }

      const std::string dark_gray() const {
        return "\033[1;30m";
      }

      const std::string white() const {
        return "\033[1;37m";
      }

      const std::string reset() const {
        return "\033[0m";
      }
#endif
    };
  }
}
#endif
