
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_CONTRAINTOPERATOR_H
#define IGLOO_CONTRAINTOPERATOR_H

namespace snowhouse {

  struct InvalidExpressionException
  {
    InvalidExpressionException(const std::string& message) : m_message(message)
    {
    }

    const std::string& Message() const
    {
      return m_message;
    }

    std::string m_message;
  };

  struct ConstraintOperator
  {
    virtual ~ConstraintOperator() {}

    virtual void PerformOperation(ResultStack& result) = 0;
    virtual int Precedence() const = 0;

    template <typename ConstraintListType, typename ActualType>
    static bool EvaluateElementAgainstRestOfExpression(ConstraintListType& list, const ActualType& actual)
    {
       ResultStack innerResult;
       OperatorStack innerOperators;

       EvaluateConstraintList(list.m_tail, innerResult, innerOperators, actual);
       EvaluateAllOperatorsOnStack(innerOperators, innerResult);

       if(innerResult.empty())
       {
         throw InvalidExpressionException("The expression after \"" + snowhouse::Stringize(list.m_head) + "\" operator does not yield any result");
       }

       return innerResult.top();
    }

    static void EvaluateOperatorsWithLessOrEqualPrecedence(const ConstraintOperator& op, OperatorStack& operators, ResultStack& result)
    {
      while(!operators.empty())
      {
        ConstraintOperator* op_from_stack = operators.top();

        if(op_from_stack->Precedence() > op.Precedence())
        {
          break;
        }

        op_from_stack->PerformOperation(result);
        operators.pop();
      }
    }

    static void EvaluateAllOperatorsOnStack(OperatorStack& operators, ResultStack& result)
    {
      while(!operators.empty())
      {
        ConstraintOperator* op = operators.top();
        op->PerformOperation(result);
        operators.pop();
      }
    }
  };

}

#endif
