// $Id: CglCutGenerator.cpp 948 2011-01-04 23:28:13Z lou $
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cstdlib>
#include <cassert>
//#include <cfloat>
//#include <iostream>

#include "CglCutGenerator.hpp"
#include "CoinHelperFunctions.hpp"


//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglCutGenerator::CglCutGenerator ()
  : aggressive_(0)
{
  // nothing to do here
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglCutGenerator::CglCutGenerator (
                  const CglCutGenerator & source)
  : aggressive_(source.aggressive_)
{
  // nothing to do here
}


//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglCutGenerator::~CglCutGenerator ()
{
  // nothing to do here
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglCutGenerator &
CglCutGenerator::operator=(
                   const CglCutGenerator& rhs)
{
  if (this != &rhs) {
    aggressive_ = rhs.aggressive_;
  }
  return *this;
}
bool
CglCutGenerator::mayGenerateRowCutsInTree() const
{
  return true;
}
// Return true if needs optimal basis to do cuts
bool
CglCutGenerator::needsOptimalBasis() const
{
  return false;
}

#ifdef NDEBUG
#undef NDEBUG
#endif

#if 0
//--------------------------------------------------------------------------
// test EKKsolution methods.
//--------------------------------------------------------------------------
void
CglCutGenerator::unitTest()
{
}
#endif
