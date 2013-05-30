/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Callable         */
/* Library.                                                                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2004-2006 Ted Ralphs and Lehigh University.                 */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* The authors of this file are Menal Guzelsoy and Ted Ralphs                */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include "SymWarmStart.hpp"

#include <iostream>

#include "symphony.h"

//#############################################################################

SymWarmStart::SymWarmStart(warm_start_desc * ws)
{
   warmStart_ = sym_create_copy_warm_start(ws);
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::SymWarmStart(char * fileName)
{

   warmStart_ = sym_read_warm_start(fileName);
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::SymWarmStart(const SymWarmStart & symWS)
{
   warm_start_desc * wsCopy;
   SymWarmStart * sWS = const_cast<SymWarmStart *>(&symWS);
   wsCopy = const_cast<warm_start_desc *>(sWS->getCopyOfWarmStartDesc());

   warmStart_ = wsCopy;
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::~SymWarmStart()
{
   sym_delete_warm_start(warmStart_);
}

/*===========================================================================*/
/*===========================================================================*/

CoinWarmStart * SymWarmStart::clone () const
{
   return new SymWarmStart(*this);
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc * SymWarmStart::getCopyOfWarmStartDesc()
{

   if(warmStart_){
      return(sym_create_copy_warm_start(warmStart_));
   }
   else{
      std::cout<<"getWarmStart(): No loaded warm start desc. to return!"<<std::endl;
      return 0;
   }
}

/*===========================================================================*/
/*===========================================================================*/

int SymWarmStart::writeToFile(char * fileName)
{
   return(sym_write_warm_start_desc(warmStart_, fileName));
}

/*===========================================================================*/
/*===========================================================================*/
