
//          Copyright Joakim Karlsson & Kim Gräsman 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef IGLOO_FLUENT_H
#define IGLOO_FLUENT_H

#include "constraintlist.h"
#include "constraintadapter.h"
#include "operators/constraintoperator.h"
#include "operators/andoperator.h"
#include "operators/oroperator.h"
#include "operators/collections/collectionconstraintevaluator.h"
#include "operators/collections/alloperator.h"
#include "operators/collections/noneoperator.h"
#include "operators/collections/atleastoperator.h"
#include "operators/collections/exactlyoperator.h"
#include "operators/collections/atmostoperator.h"
#include "operators/notoperator.h"
#include "expressionbuilder.h"

namespace snowhouse {

  inline ExpressionBuilder<Nil> Is()
  {
    return ExpressionBuilder<Nil>(Nil());
  }

  inline ExpressionBuilder<Nil> Has()
  {
     return ExpressionBuilder<Nil>(Nil());
  }

}

#endif
