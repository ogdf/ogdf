
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_ASSERT_H
#define IGLOO_ASSERT_H

#include "stringize.h"
#include "stringizers.h"

namespace snowhouse {

   struct DefaultFailureHandler
   {
      template <class ExpectedType, class ActualType>
      static void Handle(const ExpectedType& expected, const ActualType& actual, const char* file_name, int line_number)
      {
         std::ostringstream str;

         str << "Expected: " << snowhouse::Stringize(expected) << std::endl;
         str << "Actual: " << snowhouse::Stringize(actual) << std::endl;

         throw AssertionException(str.str(), file_name, line_number);
      }

      static void Handle(const std::string& message)
      {
         throw AssertionException(message);
      }
   };

   template<typename FailureHandler>
   class ConfigurableAssert
   {
   public:

      template <typename ActualType, typename ConstraintListType>
      static void That(const ActualType& actual, ExpressionBuilder<ConstraintListType> expression)
      {
        const char* no_file = "";
        int line_number = 0;

        ConfigurableAssert<FailureHandler>::That(actual, expression, no_file, line_number);
      }

      template <typename ActualType, typename ConstraintListType>
      static void That(const ActualType& actual, ExpressionBuilder<ConstraintListType> expression, const char* file_name, int line_number)
      {
         try
         {
            ResultStack result;
            OperatorStack operators;
            expression.Evaluate(result, operators, actual);

            while (!operators.empty())
            {
               ConstraintOperator* op = operators.top();
               op->PerformOperation(result);
               operators.pop();
            }

            if (result.empty())
            {
               throw InvalidExpressionException("The expression did not yield any result");
            }

            if (!result.top())
            {
               FailureHandler::Handle(expression, actual, file_name, line_number);
            }
         }
         catch (const InvalidExpressionException& e)
         {
            FailureHandler::Handle("Malformed expression: \"" + snowhouse::Stringize(expression) + "\"\n" + e.Message());
         }
      }

      template <typename ConstraintListType>
      static void That(const char* actual, ExpressionBuilder<ConstraintListType> expression)
      {
         return That(std::string(actual), expression);
      }

      template <typename ActualType, typename ExpressionType>
      static void That(const ActualType& actual, const ExpressionType& expression)
      {
        const char* no_file = "";
        int no_line = 0;
        That(actual, expression, no_file, no_line);
      }

      template <typename ActualType, typename ExpressionType>
      static void That(const ActualType& actual, const ExpressionType& expression, const char* file_name, int line_number)
      {
         if (!expression(actual))
         {
            FailureHandler::Handle(expression, actual, file_name, line_number);
         }
      }

      template <typename ExpressionType>
      static void That(const char* actual, const ExpressionType& expression)
      {
         return That(std::string(actual), expression);
      }

      static void That(bool actual)
      {
         if (!actual)
         {
            FailureHandler::Handle("Expected: true\nActual: false");
         }
      }

      static void Failure(const std::string& message)
      {
         FailureHandler::Handle(message);
      }
   };

   typedef ConfigurableAssert<DefaultFailureHandler> Assert;
}

#endif	// IGLOO_ASSERT_H
