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

#include <stdio.h>
#include <stdlib.h>

#include "sym_types.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_proccomm.h"
#include "sym_messages.h"
#include "sym_cp.h"
#include "sym_cp_u.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains CP wrapper functions that interface with the user.
\*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_receive_cp_data that
 * receives the initial data from the Master process. Returns TRUE if
 * succeeded, FALSE otherwise.
\*===========================================================================*/

int receive_cp_data_u(cut_pool *cp)
{
   int r_bufid;

   r_bufid = receive_msg(cp->master, CP_DATA);
   receive_char_array((char *)&cp->par, sizeof(cp_params));
#ifdef USE_SYM_APPLICATION
   switch( user_receive_cp_data(&cp->user) ){
    case USER_SUCCESS:
    case USER_NO_PP:
    case USER_AND_PP:
      /* User function terminated without problems. No post-processing. */
    case USER_DEFAULT:
      freebuf(r_bufid);
      return(TRUE);
    case USER_ERROR:
    default:
      freebuf(r_bufid);
      /* Unexpected return value. Do something!! */
      return(FALSE);
   }
#else
   freebuf(r_bufid);
   return(TRUE);
#endif
}

/*===========================================================================*/

int receive_lp_solution_cp_u(cut_pool *cp)
{
   int termcode = 0;
#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_receive_lp_solution_cp(&cp->user) );
#endif
   return(termcode);
}

/*===========================================================================*/

int check_cuts_u(cut_pool *cp, lp_sol *cur_sol)
{
   int num_cuts = 0, i, violated;
   cp_cut_data **pcp_cut;
   double quality;
   int cuts_to_check = MIN(cp->cut_num, cp->par.cuts_to_check);

#ifdef USE_SYM_APPLICATION
   if (user_prepare_to_check_cuts(cp->user, cur_sol->xlength, cur_sol->xind,
				  cur_sol->xval) == USER_ERROR){
      return(0);
   }
#endif

   switch(cp->par.check_which){ /* decide which cuts to check for violation */

    case CHECK_ALL_CUTS: /* check all cuts in the pool */
      for (i = 0, pcp_cut = cp->cuts; i < cuts_to_check; i++, pcp_cut++){
	 if (check_cut_u(cp, cur_sol, &(*pcp_cut)->cut,
			 &violated, &quality) == USER_ERROR)
	    break;
	 (*pcp_cut)->quality =
	    ((*pcp_cut)->quality*(double)((*pcp_cut)->check_num) + quality)/
	    (double)((*pcp_cut)->check_num+1);
	 (*pcp_cut)->check_num++;
	 if ( violated ){
	    num_cuts++;
	    (*pcp_cut)->touches = 0;
	    cut_pool_send_cut(cp, &(*pcp_cut)->cut, cur_sol->lp);
	 }else{
	    (*pcp_cut)->touches++;
	 }
      }
      break;

    case CHECK_LEVEL: /* only check cuts generated at a level higher
			 than the current level. This prevents checking
			 cuts generated in other parts of the tree which
			 are not as likely to be violated */
      for (i = 0, pcp_cut = cp->cuts; i < cuts_to_check; i++, pcp_cut++){
	 if ((*pcp_cut)->level >= cur_sol->xlevel)
	    continue;
	 if (check_cut_u(cp, cur_sol, &(*pcp_cut)->cut,
			 &violated, &quality) == USER_ERROR)
	    break;
	 (*pcp_cut)->quality =
	    ((*pcp_cut)->quality*(double)((*pcp_cut)->check_num) + quality)/
	    (double)((*pcp_cut)->check_num+1);
	 (*pcp_cut)->check_num++;
	 if ( violated ){
	    num_cuts++;
	    (*pcp_cut)->touches = 0;
	    cut_pool_send_cut(cp, &(*pcp_cut)->cut, cur_sol->lp);
	 }else{
	    (*pcp_cut)->touches++;
	 }
      }
      break;

    case CHECK_TOUCHES: /* only check cuts which have been recently
			   violated */
      for (i = 0, pcp_cut = cp->cuts; i < cuts_to_check; i++, pcp_cut++){
	 if ((*pcp_cut)->touches > cp->par.touches_until_deletion)
	    continue;
	 if (check_cut_u(cp, cur_sol, &(*pcp_cut)->cut,
			 &violated, &quality) == USER_ERROR)
	    break;
	 (*pcp_cut)->quality =
	    ((*pcp_cut)->quality*(double)((*pcp_cut)->check_num) + quality)/
	    (double)((*pcp_cut)->check_num+1);
	 (*pcp_cut)->check_num++;
	 if ( violated ){
	    num_cuts++;
	    (*pcp_cut)->touches = 0;
	    cut_pool_send_cut(cp, &(*pcp_cut)->cut, cur_sol->lp);
	 }else{
	    (*pcp_cut)->touches++;
	 }
      }
      break;

    case CHECK_LEVEL_AND_TOUCHES: /* a combination of the above two
				     options */
      for (i = 0, pcp_cut = cp->cuts; i < cuts_to_check; i++, pcp_cut++){
	 if ((*pcp_cut)->touches > cp->par.touches_until_deletion ||
	     (*pcp_cut)->level > cur_sol->xlevel)
	    continue;
	 if (check_cut_u(cp, cur_sol, &(*pcp_cut)->cut,
			 &violated, &quality) == USER_ERROR)
	    break;
	 (*pcp_cut)->quality =
	    ((*pcp_cut)->quality*(double)((*pcp_cut)->check_num) + quality)/
	    (double)((*pcp_cut)->check_num+1);
	 (*pcp_cut)->check_num++;
	 if ( violated ){
	    num_cuts++;
	    (*pcp_cut)->touches = 0;
	    cut_pool_send_cut(cp, &(*pcp_cut)->cut, cur_sol->lp);
	 }else{
	    (*pcp_cut)->touches++;
	 }
      }
      break;

   default:
      printf("Unknown rule for checking cuts \n\n");
      break;
   }
#ifdef USE_SYM_APPLICATION
   user_finished_checking_cuts(cp->user);
#endif

   return(num_cuts);
}

/*===========================================================================*/

int check_cut_u(cut_pool *cp, lp_sol *cur_sol, cut_data *cut, int *is_violated,
		double *quality)
{
   int varnum = cur_sol->xlength, nzcnt;
   int *indices = cur_sol->xind, *matind;
   double *values = cur_sol->xval, *matval;
   double lhs = 0, etol = cur_sol->lpetol;
   int i, j;

   switch (cut->type){

    case EXPLICIT_ROW:
      nzcnt = ((int *) (cut->coef))[0];
      matval = (double *) (cut->coef + DSIZE);
      matind = (int *) (cut->coef + (nzcnt+ 1) * DSIZE);
      for (i = 0, j = 0; i < nzcnt && j < varnum; ){
	 if (matind[i] == indices[j]){
	    lhs += matval[i++]*values[j++];
	 }else if (matind[i] < indices[j]){
	    i++;
	 }else if (matind[i] > indices[j]){
	    j++;
	 }
      }

      switch (cut->sense){

       case 'G':

	 *is_violated = (lhs < cut->rhs - etol);
	 *quality = cut->rhs - lhs;

       case 'L':

	 *is_violated = (lhs > cut->rhs + etol);
	 *quality = lhs - cut->rhs;

       case 'R':

	 if (cut->range > 0){
	    *is_violated = ((lhs < cut->rhs - etol) ||
			    (lhs > cut->rhs + cut->range + etol));
	    *quality = lhs < cut->rhs - etol ? cut->rhs - lhs :
	       lhs - cut->rhs + cut->range;
	 }else{
	    *is_violated = ((lhs > cut->rhs + etol) ||
			    (lhs < cut->rhs + cut->range - etol));
	    *quality = lhs > cut->rhs + etol ? lhs - cut->rhs :
	       cut->rhs + cut->range - lhs;
	 }
      }
      return(0);

    default:
#ifdef USE_SYM_APPLICATION
      return(user_check_cut(cp->user, cur_sol->lpetol, varnum, indices, values,
			    cut, is_violated, quality));
#else
      return(USER_DEFAULT);
#endif
   }
}


/*===========================================================================*/

void free_cut_pool_u(cut_pool *cp)
{
   int i;
#ifdef USE_SYM_APPLICATION
   user_free_cp(&cp->user);
#endif
   for (i = cp->cut_num - 1; i >= 0; i--){
      FREE(cp->cuts[i]->cut.coef);
      FREE(cp->cuts[i]);
   }
   FREE(cp->cuts);
   FREE(cp->cur_sol.xind);
   FREE(cp->cur_sol.xval);
#ifdef COMPILE_IN_CP
   FREE(cp->cuts_to_add);
#endif
   FREE(cp);
}

/*===========================================================================*/

