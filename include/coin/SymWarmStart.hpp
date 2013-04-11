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

#ifndef SymWarmStart_H
#define SymWarmStart_H

#include "CoinWarmStart.hpp"

typedef struct WARM_START_DESC warm_start_desc;

//#############################################################################

class SymWarmStart : public CoinWarmStart 
{

public:

   /* Default constructor. Will do nothing! */
   SymWarmStart(){}
   
   /* Initialize the warmStart_ using the given warm start. If dominate
      WarmStart is set, then, SymWarmStart will take the control of the 
      given description, otherwise, will copy everything.
   */
   SymWarmStart(warm_start_desc * ws);
   
   /*Get the warmStart info from a file*/
   SymWarmStart(char *f);
   
   /* Copy constructor */
   SymWarmStart(const SymWarmStart & symWS);

   /* Destructor */
   virtual ~SymWarmStart();

   /* Clone the warmstart */
   virtual CoinWarmStart * clone() const; 

   /* Get the pointer to the loaded warmStart_ */
   virtual warm_start_desc * getCopyOfWarmStartDesc();

   /* Move the pointer to the rootnode of the warmStart to another
      node which will change the underlying tree 
   */
   // virtual void setRoot(bc_node *root) {} //FIX_ME! Ask Prof. Ralphs.

   /* Write the current warm start info to a file */
   virtual int writeToFile(char * f);

private:

   /* Private warm start desc. to keep everything */
   warm_start_desc *warmStart_;

};

#endif
