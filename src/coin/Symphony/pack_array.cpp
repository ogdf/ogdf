/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2011 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>         /* has malloc() in AIX ... */

#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_pack_array.h"
#include "sym_proccomm.h"

/*===========================================================================*/

void pack_array_desc(array_desc *adesc)
{
   send_char_array((char *)adesc, sizeof(array_desc));
   if (adesc->type != NO_DATA_STORED && adesc->size > 0)
      send_int_array(adesc->list, adesc->size);
}

/*===========================================================================*/

array_desc *unpack_array_desc(array_desc *padesc)
{
   array_desc *adesc =
      padesc ? padesc : (array_desc *) malloc( sizeof(array_desc) );
   receive_char_array((char *)adesc, sizeof(array_desc));
   if (adesc->type != NO_DATA_STORED && adesc->size > 0){
      adesc->list = (int *) malloc(adesc->size * ISIZE);
      receive_int_array(adesc->list, adesc->size);
   }else{
      adesc->list = NULL;
   }
   if (adesc->type == EXPLICIT_LIST)
      adesc->added = adesc->size;
   return(adesc);
}

/*===========================================================================*/

void pack_double_array_desc(double_array_desc *dad, char explicit_packing)
{
   send_char_array(&dad->type, 1);
   send_int_array(&dad->size, 1);
   if (dad->size > 0){
      if (!explicit_packing && dad->type == WRT_PARENT)
	 send_int_array(dad->list, dad->size);
      send_int_array(dad->stat, dad->size);
   }
}

/*===========================================================================*/

void unpack_double_array_desc(double_array_desc *dad, char explicit_packing)
{
   receive_char_array(&dad->type, 1);
   receive_int_array(&dad->size, 1);
   if (dad->size > 0){
      if (!explicit_packing && dad->type == WRT_PARENT){
	 dad->list = (int *) malloc(dad->size * ISIZE);
	 receive_int_array(dad->list, dad->size);
      }else{
	 dad->list = NULL;
      }
      dad->stat = (int *) malloc(dad->size * ISIZE);
      receive_int_array(dad->stat, dad->size);
   }else{
      dad->list = NULL;
      dad->stat = NULL;
   }
}

/*===========================================================================*/

void pack_basis(basis_desc *basis, char explicit_packing)
{
   send_char_array(&basis->basis_exists, 1);
   if (basis->basis_exists){
      pack_double_array_desc(&basis->baserows, explicit_packing);
      pack_double_array_desc(&basis->extrarows, explicit_packing);
      pack_double_array_desc(&basis->basevars, explicit_packing);
      pack_double_array_desc(&basis->extravars, explicit_packing);
   }
}

/*===========================================================================*/

basis_desc *unpack_basis(basis_desc *pbasis, char explicit_packing)
{
   basis_desc *basis =
      pbasis ? pbasis : (basis_desc *) calloc(1, sizeof(basis_desc) );
   receive_char_array(&basis->basis_exists, 1);
   if (basis->basis_exists){
      unpack_double_array_desc(&basis->baserows, explicit_packing);
      unpack_double_array_desc(&basis->extrarows, explicit_packing);
      unpack_double_array_desc(&basis->basevars, explicit_packing);
      unpack_double_array_desc(&basis->extravars, explicit_packing);
   }else{
      basis->baserows.list = basis->baserows.stat = NULL;
      basis->extrarows.list = basis->extrarows.stat = NULL;
      basis->basevars.list = basis->basevars.stat = NULL;
      basis->extravars.list = basis->extravars.stat = NULL;
   }

   return(basis);
}

