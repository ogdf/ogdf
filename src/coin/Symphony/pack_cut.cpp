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

#include <stdlib.h>        /* has malloc() in AIX ... */

#include "sym_constants.h"
#include "sym_pack_cut.h"
#include "sym_proccomm.h"

/*===========================================================================*/

void pack_cut(cut_data *cut)
{
   send_char_array((char *)cut, sizeof(cut_data));
   if (cut->size > 0)
      send_char_array(cut->coef, cut->size);
}

/*===========================================================================*/

cut_data *unpack_cut(cut_data *pcut)
{
   cut_data *cut = pcut ? pcut : (cut_data *) malloc(sizeof(cut_data));
   char *coef = pcut ? pcut->coef : NULL;

   receive_char_array((char *)cut, sizeof(cut_data));
   cut->coef = coef;
   if (cut->size > 0){
      if (!cut->coef)
	 cut->coef = (char *) malloc(cut->size);
      receive_char_array(cut->coef, cut->size);
   }
   return(cut);
}
