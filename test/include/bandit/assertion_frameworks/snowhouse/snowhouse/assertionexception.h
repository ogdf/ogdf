
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_ASSERTIONEXCEPTION_H
#define IGLOO_ASSERTIONEXCEPTION_H

namespace snowhouse {
  class AssertionException : public std::exception
  {
    public:
      AssertionException(const std::string& message)
        : m_message(message), m_fileName(""), m_line(0)
      {}

      AssertionException(const std::string& message, const std::string& fileName, unsigned int line)
        : m_message(message), m_fileName(fileName), m_line(line)
      {}

      virtual ~AssertionException() throw()
      {
      }

      std::string GetMessage() const
      {
        return m_message;
      }

      std::string GetFilename() const
      {
        return m_fileName;
      }

      unsigned int GetLineNumber() const
      {
        return m_line;
      }

    private:
      std::string m_message;
      std::string m_fileName;
      unsigned int m_line;
  };
}

#endif // IGLOO_ASSERTIONEXCEPTION_H
