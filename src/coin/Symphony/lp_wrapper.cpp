/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2011 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* The OSI interface in this file was written by Menal Guzelsoy.             */
/* The OSL interface was written by Ondrej Medek.                            */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "sym_lp.h"
#include "sym_master.h"
#include "sym_proccomm.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_types.h"
#include "sym_lp_solver.h"
#include "sym_primal_heuristics.h"
#ifdef USE_CGL_CUTS
#include "sym_cg.h"
#endif
#if defined (COMPILE_IN_LP) && defined (COMPILE_IN_TM)
#include "sym_master_u.h"
#endif
#ifdef COMPILE_IN_CP
#include "sym_cp.h"
#endif

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP wrapper functions that interface with the user.
\*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_receive_lp_data that
 * receives the initial data from the Master process. Returns TRUE if
 * succeeded, FALSE otherwise.
\*===========================================================================*/

int receive_lp_data_u(lp_prob *p)
{
   int r_bufid;
   char has_desc;
   char has_colnames;
   int i;

   r_bufid = receive_msg(p->master, LP_DATA);
   receive_char_array((char *)(&p->par), sizeof(lp_params));
   receive_int_array(&p->has_ub, 1);
   if (p->has_ub){
      receive_dbl_array(&p->ub, 1);
   }else{
      p->ub = - (MAXDOUBLE / 2);
   }
   if(p->par.multi_criteria){
      receive_int_array(&p->has_mc_ub, 1);
      if (p->has_mc_ub){
	 receive_dbl_array(&p->mc_ub, 1);
	 receive_dbl_array(p->obj, 2);
      }else{
	 p->mc_ub = - (MAXDOUBLE / 2);
      }
      receive_dbl_array(p->utopia, 2);
   }
   receive_int_array(&p->draw_graph, 1);
   receive_int_array(&p->base.varnum, 1);
   if (p->base.varnum > 0){
      p->base.userind = (int *) malloc(p->base.varnum * ISIZE);
      receive_int_array(p->base.userind, p->base.varnum);
   }
   receive_int_array(&p->base.cutnum, 1);
   MIPdesc *mip = p->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));
   receive_int_array(&(mip->m), 1);
   receive_int_array(&(mip->n), 1);
   receive_int_array(&(mip->nz), 1);
   receive_char_array(&(mip->obj_sense), 1);
   receive_dbl_array(&(mip->obj_offset), 1);
   receive_char_array(&has_desc, 1);

   if (has_desc){
      /* Allocate memory */
      mip->matbeg = (int *) malloc(ISIZE * (mip->n + 1));
      mip->matind = (int *)    malloc(ISIZE * mip->nz);
      mip->matval = (double *) malloc(DSIZE * mip->nz);
      mip->obj    = (double *) malloc(DSIZE * mip->n);
      if (p->par.multi_criteria){
	 mip->obj1    = (double *) malloc(DSIZE * mip->n);
	 mip->obj2    = (double *) malloc(DSIZE * mip->n);
      }
      mip->rhs    = (double *) malloc(DSIZE * mip->m);
      mip->sense  = (char *)   malloc(CSIZE * mip->m);
      mip->rngval = (double *) malloc(DSIZE * mip->m);
      mip->ub     = (double *) malloc(DSIZE * mip->n);
      mip->lb     = (double *) malloc(DSIZE * mip->n);
      mip->is_int = (char *)   calloc(CSIZE, mip->n);

      /* Receive the problem description */
      receive_int_array(mip->matbeg, mip->n+1);
      receive_int_array(mip->matind, mip->nz);
      receive_dbl_array(mip->matval, mip->nz);
      receive_dbl_array(mip->obj, mip->n);
      if (p->par.multi_criteria){
	 receive_dbl_array(mip->obj1, mip->n);
	 receive_dbl_array(mip->obj2, mip->n);
      }
      receive_dbl_array(mip->rhs, mip->m);
      receive_char_array(mip->sense, mip->m);
      receive_dbl_array(mip->rngval, mip->m);
      receive_dbl_array(mip->ub, mip->n);
      receive_dbl_array(mip->lb, mip->n);
      receive_char_array(mip->is_int, mip->n);
      receive_char_array(&has_colnames, 1);
      if (has_colnames){
	 mip->colname = (char **) malloc(sizeof(char *) * mip->n);
	 for (i = 0; i < mip->n; i++){
	    mip->colname[i] = (char *) malloc(CSIZE * 9);
	    receive_char_array(mip->colname[i], 8);
	    mip->colname[i][8] = 0;
	 }
      }
   }

#ifdef USE_SYM_APPLICATION
   switch( user_receive_lp_data(&p->user)){
    case USER_ERROR:
      freebuf(r_bufid);
      return(ERROR__USER);
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
      break;
    default:
      freebuf(r_bufid);
      /* Unexpected return value. Do something!! */
      return(ERROR__USER);
   }
#endif

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int comp_cut_name(const void *c0, const void *c1)
{
   return((*((cut_data **)c0))->name - (*((cut_data **)c1))->name);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_create_subproblem that
 * creates the problem matrix.
\*===========================================================================*/

int create_subproblem_u(lp_prob *p)
{
   node_desc *desc = p->desc;

   LPdata *lp_data = p->lp_data;
   int i, j, k, maxm, maxn, maxnz;
   row_data *row, *rows;

   int bvarnum = p->base.varnum;
   int bcutnum = p->base.cutnum;
   var_desc **vars;

   int *d_uind = NULL, *d_cind = NULL; /* just to keep gcc quiet */

   double *rhs, *rngval, *darray;
   char *sense;
   char *status;
   cut_data *cut;
   branch_desc *bobj;

   int new_row_num;
   waiting_row **new_rows;

   int user_res;
   int *userind;

   double *lb, *ub;
   char *is_int;
   MIPdesc *lp_data_mip, *p_mip, *tmp_mip;
   char *lp_data_mip_is_int, *p_mip_is_int;
   double *lp_data_mip_ub, *p_mip_ub;
   double *lp_data_mip_lb, *p_mip_lb;
   double *lp_data_mip_obj, *p_mip_obj;
   int *lp_data_mip_matbeg, *p_mip_matbeg;
   int *lp_data_mip_matind, *p_mip_matind;
   double *lp_data_mip_matval, *p_mip_matval;
   char *lp_data_mip_sense, *p_mip_sense;
   double *lp_data_mip_rhs, *p_mip_rhs;
   double *lp_data_mip_rngval, *p_mip_rngval;

   p->par.lp_data_mip_is_copied = TRUE;
   p_mip = p->mip;
   tmp_mip = lp_data->mip;

   lp_data->n = bvarnum + desc->uind.size;
   lp_data->m = bcutnum + desc->cutind.size;

   maxm = lp_data->maxm;
   maxn = lp_data->maxn;
   maxnz = lp_data->maxnz;

   lp_data->nf_status = desc->nf_status;
   if (desc->nf_status == NF_CHECK_AFTER_LAST ||
       desc->nf_status == NF_CHECK_UNTIL_LAST){
      lp_data->not_fixed_num = desc->not_fixed.size;
      memcpy(lp_data->not_fixed, desc->not_fixed.list,
	     lp_data->not_fixed_num * ISIZE);
   }

   if (desc->uind.size > 0){ /* fill up the rest of lp_data->vars */
      if (MAX(maxn, bvarnum) < lp_data->n){
	 vars = lp_data->vars = (var_desc **)
	    realloc(lp_data->vars, lp_data->n * sizeof(var_desc *));
	 for (i = lp_data->n - 1; i >= MAX(maxn, bvarnum); i--){
	    vars[i] = (var_desc *) malloc(sizeof(var_desc) );
	 }
      }
      vars = lp_data->vars;
      d_uind = desc->uind.list;
      for (i = desc->uind.size - 1; i >= 0; i--){
	 vars[i + bvarnum]->userind = d_uind[i];
	 vars[i + bvarnum]->colind = bvarnum + i;
      }

      if (p->par.multi_criteria && !p->par.mc_find_supported_solutions){
	 vars[lp_data->n - 1]->userind = p_mip->n;
	 vars[lp_data->n - 1]->colind  = lp_data->n - 1;
      }
   }
   lp_data->ordering = COLIND_AND_USERIND_ORDERED;

   lp_data_mip    = lp_data->mip;
   lp_data_mip->n = lp_data->n;
   lp_data_mip->m = lp_data->m;

   /* Create the list of indices to pass to the user */
   userind = (int *) malloc(lp_data->n * ISIZE);
   vars = lp_data->vars;
   p->par.is_userind_in_order = TRUE;
   for (i = lp_data->n - 1; i >= 0; --i){
      userind[i] = vars[i]->userind;
      if (userind[i]!=i) {
         p->par.is_userind_in_order = FALSE;
      }
   }
   if (lp_data->n != p_mip->n) {
      p->par.is_userind_in_order = FALSE;
   }

#ifdef USE_SYM_APPLICATION
   user_res = user_create_subproblem(p->user,
       /* list of base and extra variables */
       userind,
       /* description of the LP relaxation to be filled out by the user */
       lp_data_mip,
       /* max sizes (estimated by the user) */
       &maxn, &maxm, &maxnz);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){

    case USER_DEFAULT:

      if (!p_mip->matbeg){
	 printf("Illegal return code.\n");
	 printf("Trying to use default user_create_subproblem without");
	 printf("reading an MPS or AMPL file. Exiting...\n\n");
	 return(ERROR__ILLEGAL_RETURN_CODE);
      }

      lp_data_mip->nz = p_mip->nz;
      if (p->par.multi_criteria && !p->par.mc_find_supported_solutions){
	 lp_data_mip->nz += 2 * lp_data->n;
      }

      if (p->par.is_userind_in_order == FALSE || p->bc_index == 0) {
         /* Allocate the arrays.*/
         lp_data_mip->matbeg  = (int *) malloc((lp_data_mip->n + 1) * ISIZE);
         lp_data_mip->matind  = (int *) malloc((lp_data_mip->nz) * ISIZE);
         lp_data_mip->matval  = (double *) malloc((lp_data_mip->nz) * DSIZE);
         lp_data_mip->obj     = (double *) malloc(lp_data_mip->n * DSIZE);
         lp_data_mip->ub      = (double *) malloc(lp_data_mip->n * DSIZE);
         lp_data_mip->lb      = (double *) calloc(lp_data_mip->n, DSIZE);
         lp_data_mip->rhs     = (double *) malloc(lp_data_mip->m * DSIZE);
         lp_data_mip->sense   = (char *)   malloc(lp_data_mip->m * CSIZE);
         lp_data_mip->rngval  = (double *) calloc(lp_data_mip->m, DSIZE);
         lp_data_mip->is_int  = (char *)   calloc(lp_data_mip->n, CSIZE);

         lp_data_mip_is_int        = lp_data_mip->is_int;
         lp_data_mip_ub            = lp_data_mip->ub;
         lp_data_mip_lb            = lp_data_mip->lb;
         lp_data_mip_obj           = lp_data_mip->obj;
         lp_data_mip_matbeg        = lp_data_mip->matbeg;
         lp_data_mip_matind        = lp_data_mip->matind;
         lp_data_mip_matval        = lp_data_mip->matval;
         lp_data_mip_sense         = lp_data_mip->sense;
         lp_data_mip_rhs           = lp_data_mip->rhs;
         lp_data_mip_rngval        = lp_data_mip->rngval;

         p_mip_is_int              = p_mip->is_int;
         p_mip_ub                  = p_mip->ub;
         p_mip_lb                  = p_mip->lb;
         p_mip_obj                 = p_mip->obj;
         p_mip_matbeg              = p_mip->matbeg;
         p_mip_matind              = p_mip->matind;
         p_mip_matval              = p_mip->matval;
         p_mip_sense               = p_mip->sense;
         p_mip_rhs                 = p_mip->rhs;
         p_mip_rngval              = p_mip->rngval;
         /* Fill out the appropriate data structures*/
         lp_data_mip_matbeg[0]    = 0;
         for (j = 0, i = 0; i < lp_data_mip->n; i++){
            if (userind[i] == p_mip->n){
               /* We should only be here with multi-criteria problems. */
               /* This is the artifical variable added for finding nondominated */
               /* solutions. */
               lp_data_mip_is_int[i]    = FALSE;
               lp_data_mip_ub[i]        = MAXINT;
               lp_data_mip_lb[i]        = -MAXINT;
               lp_data_mip_obj[i]       = 1.0;
               lp_data_mip_matval[j]    = -1.0;
               lp_data_mip_matind[j++]  = bcutnum - 2;
               lp_data_mip_matval[j]    = -1.0;
               lp_data_mip_matind[j++]  = bcutnum - 1;
               lp_data_mip_matbeg[i+1]  = j;
               continue;
            }
            lp_data_mip_ub[i] = p_mip_ub[userind[i]];
            lp_data_mip_lb[i] = p_mip_lb[userind[i]];
            lp_data_mip_is_int[i] = p_mip_is_int[userind[i]];
            for (k = p_mip_matbeg[userind[i]]; k < p_mip_matbeg[userind[i]+1];
                  k++){
               lp_data_mip_matind[j]   = p_mip_matind[k];
               lp_data_mip_matval[j++] = p_mip_matval[k];
            }
            if (p->par.multi_criteria && !p->par.mc_find_supported_solutions){
               lp_data_mip_obj[i] = p->par.mc_rho*(p_mip->obj1[userind[i]] +
                     p_mip->obj2[userind[i]]);
               lp_data_mip_matval[j] = p->par.mc_gamma*p_mip->obj1[userind[i]];
               lp_data_mip_matind[j++] = bcutnum - 2;
               lp_data_mip_matval[j] = p->par.mc_tau*p_mip->obj2[userind[i]];
               lp_data_mip_matind[j++] = bcutnum - 1;
            }else{
               lp_data_mip_obj[i] = p_mip_obj[userind[i]];
            }
            lp_data_mip_matbeg[i+1] = j;
         }
         lp_data_mip->nz = j;
         for (i = 0; i < p_mip->m; i++){
            lp_data_mip_rhs[i] = p_mip_rhs[i];
            lp_data_mip_sense[i] = p_mip_sense[i];
            lp_data_mip_rngval[i] = p_mip_rngval[i];
         }
         if (p->par.multi_criteria && !p->par.mc_find_supported_solutions){
            lp_data_mip_rhs[bcutnum - 2] = p->par.mc_gamma * p->utopia[0];
            lp_data_mip_sense[bcutnum - 2] = 'L';
            lp_data_mip_rhs[bcutnum - 1] = p->par.mc_tau * p->utopia[1];
            lp_data_mip_sense[bcutnum - 1] = 'L';
         }
      } else {
         p->par.lp_data_mip_is_copied = FALSE;
         lp_data->mip = p->mip;
         lp_data_mip = p->mip;
      }
      maxm = lp_data->m;
      maxn = lp_data->n;
      maxnz = lp_data->nz;
      lp_data->m = bcutnum;
      lp_data->nz = lp_data_mip->nz;
      break;
    case USER_SUCCESS:
       /* Fall through to next case */
    case USER_AND_PP:
    case USER_NO_PP:

      /* User function terminated without problems. In the post-processing
       * the extra cuts are added. HOWEVER, this might not be done until the
       * problem is loaded into the lp solver (for cplex it is not possible).
       * So for now just reset lp_data->m, do everything to load in the
       * stuff into the lp solver then come back to adding the cuts. */

       /* change obj coeff only if the obj funct. was created through
	  user_create_subproblem() */

      if (p_mip->obj_sense == SYM_MAXIMIZE){
         lp_data_mip_obj = lp_data_mip->obj;
         for (i = 0; i < lp_data_mip->n; i++){
	    lp_data_mip_obj[i] *= -1.0;
	 }
      }

      lp_data->m = bcutnum;
      lp_data->nz = lp_data_mip->nz;
      break;

    case USER_ERROR:

      /* Error. The search tree node will not be processed. */
      FREE(userind);
      return(ERROR__USER);

    default:

      /* Unexpected return value. Do something!! */
      FREE(userind);
      return(ERROR__USER);
   }

   FREE(userind); /* No longer needed */

   /*------------------------------------------------------------------------*\
    * Let's see about reallocing...
   \*----------------------------------------------------------------------- */

   if (maxm  < lp_data->m)  maxm  = lp_data->m;
   if (maxn  < lp_data->n)  maxn  = lp_data->n;
   if (maxnz < lp_data->nz) maxnz = lp_data->nz;

   size_lp_arrays(lp_data, FALSE, TRUE, maxm, maxn, maxnz);

   /* generate the random hash. useful for checking duplicacy of cuts and
    * solutions from feasibility pump
    */

   if (p->par.is_userind_in_order == FALSE || p->bc_index == 0) {
      darray = lp_data->random_hash;
      for (i=0; i<lp_data->n; i++) {
         darray[i] = CoinDrand48();
      }
   }

   if (lp_data->maxn > lp_data->n){
      vars = lp_data->vars = (var_desc **)
	 realloc(lp_data->vars, lp_data->maxn * sizeof(var_desc *));
      for (i = lp_data->n; i < lp_data->maxn; i++){
	 vars[i] = (var_desc *) malloc( sizeof(var_desc) );
      }
   }

   // TODO: fix char vs int
   /* Default status of every variable is NOT_FIXED */
   status = lp_data->status;
   if (bvarnum > 0) {
      //memset(lp_data->status, NOT_FIXED | BASE_VARIABLE, bvarnum);
      for (i=0; i<bvarnum; i++) {
         status[i] = NOT_FIXED | BASE_VARIABLE;
      }
   }
   if (bvarnum < lp_data->n) {
      //memset(lp_data->status + bvarnum, NOT_FIXED, lp_data->n - bvarnum);
      for (i=bvarnum; i<lp_data->n; i++) {
         status[i] = NOT_FIXED;
      }
   }

   /*------------------------------------------------------------------------*\
    * Set the necessary fields in rows
   \*----------------------------------------------------------------------- */

   rows = lp_data->rows;
   rhs = lp_data_mip->rhs;
   rngval = lp_data_mip->rngval;
   sense = lp_data_mip->sense;
   for (i = bcutnum - 1; i >= 0; i--){
      row = rows + i;
      cut = row->cut;
      cut->rhs = rhs[i];
      cut->range = rngval[i];
      cut->branch = (((cut->sense = sense[i]) != 'E') ?
		     ALLOWED_TO_BRANCH_ON : DO_NOT_BRANCH_ON_THIS_ROW);
      cut->size = 0;
      row->eff_cnt = 1;
      row->free = FALSE;
      cut->name = BASE_CONSTRAINT;
      cut->type = ORIGINAL_CONSTRAINT;
   }

   /*------------------------------------------------------------------------*\
    * Set the upper and lower bounds and integer status
   \*----------------------------------------------------------------------- */

   lb = lp_data_mip->lb;
   ub = lp_data_mip->ub;
   is_int = lp_data_mip->is_int;

   vars = lp_data->vars;
   for (i = lp_data->n - 1; i >= 0; i--){
      vars[i]->lb = vars[i]->new_lb = lb[i];
      vars[i]->ub = vars[i]->new_ub = ub[i];
      vars[i]->is_int = is_int[i];
   }

   /*------------------------------------------------------------------------*\
    * Load the lp problem (load_lp is an lp solver dependent routine).
   \*----------------------------------------------------------------------- */

   if (p->bc_index == 0 || p->par.should_reuse_lp == FALSE) {
      load_lp_prob(lp_data, p->par.scaling, p->par.fastmip); //load new
   } else {
      reset_lp_prob(lp_data, p->par.scaling, p->par.fastmip); //use old
   }


   /* Free the user's description */
   if (p->par.lp_data_mip_is_copied == TRUE) {
      free_mip_desc(lp_data_mip);
   }
   lp_data->mip = tmp_mip;

   if (desc->cutind.size > 0){
      unpack_cuts_u(p, CUT_FROM_TM, UNPACK_CUTS_SINGLE,
		    desc->cutind.size, desc->cuts, &new_row_num, &new_rows);
      add_row_set(p, new_rows, new_row_num);
      FREE(new_rows);
   }

   /* We don't need the cuts anymore. Free them. */
   if (desc->cutind.size > 0){
#ifndef COMPILE_IN_LP /*If we are using shared memory, we don't need to free*/
      free_cuts(desc->cuts, desc->cutind.size);
#endif
      FREE(desc->cuts);
   }else{
      desc->cuts = NULL;
   }
   lp_data->cgl = p->par.cgl;

#ifdef COMPILE_IN_LP
   /* reliability branching */
   /* pseudo costs and reliability measures */
   if (p->tm->pcost_down==NULL) {
      p->pcost_down = (double *)calloc(p->mip->n, DSIZE);
      p->pcost_up = (double *)calloc(p->mip->n, DSIZE);
      p->br_rel_down = (int *)calloc(p->mip->n, ISIZE);
      p->br_rel_up = (int *)calloc(p->mip->n, ISIZE);
      p->br_rel_cand_list = (int *)calloc(p->mip->n, ISIZE);
      p->br_rel_down_min_level = (int *)malloc(p->mip->n* ISIZE);
      p->br_rel_up_min_level = (int *)malloc(p->mip->n* ISIZE);
      for(i = 0; i <p->mip->n; i++){
	 p->br_rel_down_min_level[i] =
	    p->br_rel_up_min_level[i] = (int)1e7;
      }
      p->tm->pcost_down = p->pcost_down;
      p->tm->pcost_up = p->pcost_up;
      p->tm->br_rel_down = p->br_rel_down;
      p->tm->br_rel_up = p->br_rel_up;
      p->tm->br_rel_cand_list = p->br_rel_cand_list;
      p->tm->br_rel_down_min_level = p->br_rel_down_min_level;
      p->tm->br_rel_up_min_level = p->br_rel_up_min_level;
   } else {
      p->pcost_down = p->tm->pcost_down;
      p->pcost_up = p->tm->pcost_up;
      p->br_rel_down = p->tm->br_rel_down;
      p->br_rel_up = p->tm->br_rel_up;
      p->br_rel_down_min_level = p->tm->br_rel_down_min_level;
      p->br_rel_up_min_level = p->tm->br_rel_up_min_level;
   }
#endif

   p->str_br_check = TRUE;

   /*------------------------------------------------------------------------*\
    * Now go through the branching stuff
   \*----------------------------------------------------------------------- */
   if (p->par.lp_data_mip_is_copied == FALSE) {
      /* first reset all bounds */
      for (j=0; j<lp_data->n; j++) {
         change_ub(lp_data, j, p->mip->ub[j]);
         change_lb(lp_data, j, p->mip->lb[j]);
      }
   }

   d_cind = desc->cutind.list;
   vars = lp_data->vars;
   rows = lp_data->rows;
   if (p->bc_level){
      status = lp_data->status;
      for (i = 0; i < p->bc_level; i++){
	 bobj = p->bdesc + i;
         //bd_change = p->bnd_change + i;
	 if (bobj->type == BRANCHING_VARIABLE){
	    j = bobj->name < 0 ? /* base variable : extra variable */
	       -bobj->name-1 :
	       bfind(bobj->name, d_uind, desc->uind.size) + bvarnum;
	    switch (bobj->sense){
	     case 'E':
	       change_lbub(lp_data, j, bobj->rhs, bobj->rhs);
	       vars[j]->lb = vars[j]->ub = bobj->rhs;
	       vars[j]->new_lb = vars[j]->new_ub = bobj->rhs;
	       break;
	     case 'L':
	       change_ub(lp_data, j, bobj->rhs);
	       vars[j]->ub = bobj->rhs;
	       vars[j]->new_ub = bobj->rhs;
	       break;
	     case 'G':
	       change_lb(lp_data, j, bobj->rhs);
	       vars[j]->lb = bobj->rhs;
	       vars[j]->new_lb = bobj->rhs;
	       break;
	     case 'R':
	       change_lbub(lp_data, j, bobj->rhs, bobj->rhs + bobj->range);
	       vars[j]->lb = bobj->rhs;
	       vars[j]->new_lb = bobj->rhs;
	       vars[j]->ub = bobj->rhs + bobj->range;
	       vars[j]->new_ub = bobj->rhs + bobj->range;
	       break;
	    }
	    status[j] |= VARIABLE_BRANCHED_ON;
	 }else{ /* BRANCHING_CUT */
	    j = bobj->name < 0 ? /* base constraint : extra constraint */
	       -bobj->name-1 :
	       bfind(bobj->name, d_cind, desc->cutind.size) + bcutnum;
	    change_row(lp_data, j, bobj->sense, bobj->rhs, bobj->range);
#ifdef COMPILE_IN_LP
	    /* Because these cuts are shared with the treemanager, we have to
	       make a copy before changing them if the LP is compiled in */
	    cut = (cut_data *) malloc(sizeof(cut_data));
	    memcpy((char *)cut, (char *)rows[j].cut, sizeof(cut_data));
	    if (cut->size){
	       cut->coef = (char *) malloc(cut->size);
	       memcpy((char *)cut->coef, (char *)rows[j].cut->coef,
		      cut->size);
	    }
	    rows[j].cut = cut;
#else
	    cut = rows[j].cut;
#endif
	    cut->rhs = bobj->rhs;
	    cut->range = bobj->range;
	    cut->sense = bobj->sense;
	    cut->branch |= CUT_BRANCHED_ON;
	 }
      }
   }

   /*------------------------------------------------------------------------*\
    * Change bounds of variables that got changed in previous nodes
   \*----------------------------------------------------------------------- */
   /*
   for (i=0; i<lp_data->n; i++) {
      if (vars[i]->lb != vars[i]->new_lb) {
         printf("new lb of %d = %f, old = %f\n", i, vars[i]->new_lb,
         vars[i]->lb);
      }
      if (vars[i]->ub != vars[i]->new_ub) {
         printf("new ub of %d = %f, old = %f\n", i, vars[i]->new_ub,
         vars[i]->ub);
      }
   }
   */
#ifdef COMPILE_IN_LP
   if (p->desc->bnd_change) {
      bounds_change_desc *bnd_change = p->desc->bnd_change;
      int *index = bnd_change->index;
      char *lbub = bnd_change->lbub;
      double *value = bnd_change->value;
      int tmp_index = -1;
      for (i=0; i<bnd_change->num_changes; i++) {
         tmp_index = -1;
         if (vars[index[i]]->userind == index[i]) {
            tmp_index = index[i];
         } else {
            for (j=0; j<lp_data->n; j++) {
               if (vars[j]->userind==index[i]) {
                  tmp_index = j;
               }
            }
         }
         if (tmp_index<0) {
            /*
             * the variable with userind index[i] does not exist in this
             * formulation
             */
            continue;
         }
         if (lbub[i] == 'L') {
            if (vars[tmp_index]->lb<value[i]) {
               vars[tmp_index]->lb = value[i];
               vars[tmp_index]->new_lb = value[i];
               change_lb(lp_data, tmp_index, value[i]);
            }
         }
         if (lbub[i] == 'U') {
            if (vars[tmp_index]->ub>value[i]) {
               vars[tmp_index]->ub = value[i];
               vars[tmp_index]->new_ub = value[i];
               change_ub(lp_data, tmp_index, value[i]);
            }
         }
      }
      /* p->desc->bnd_change_desc no longer needed. free it */
      FREE(bnd_change->index);
      FREE(bnd_change->lbub);
      FREE(bnd_change->value);
      FREE(p->desc->bnd_change);
   }
#else
   p->desc->bnd_change = NULL;
#endif

   /*------------------------------------------------------------------------*\
    * The final step: load in the basis.
    * This is cplex style. sorry about it... Still, it
    * might be ok if {VAR,SLACK}_{B,LB,UB} are properly defined
   \*----------------------------------------------------------------------- */

   if (p->par.should_warmstart_chain == TRUE &&
         desc->basis.basis_exists == TRUE){
      int *rstat, *cstat;
      if (desc->basis.extravars.size == 0){
	 cstat = desc->basis.basevars.stat;
      }else if (desc->basis.basevars.size == 0){
	 cstat = desc->basis.extravars.stat;
      }else{ /* neither is zero */
	 cstat = lp_data->tmp.i1; /* n */
	 memcpy(cstat,
		desc->basis.basevars.stat, desc->basis.basevars.size *ISIZE);
	 memcpy(cstat + desc->basis.basevars.size,
		desc->basis.extravars.stat, desc->basis.extravars.size *ISIZE);
      }
      if (desc->basis.extrarows.size == 0){
	 rstat = desc->basis.baserows.stat;
      }else if (desc->basis.baserows.size == 0){
	 rstat = desc->basis.extrarows.stat;
      }else{ /* neither is zero */
	 rstat = lp_data->tmp.i2; /* m */
	 memcpy(rstat,
		desc->basis.baserows.stat, desc->basis.baserows.size *ISIZE);
	 memcpy(rstat + desc->basis.baserows.size,
		desc->basis.extrarows.stat, desc->basis.extrarows.size *ISIZE);
      }
      load_basis(lp_data, cstat, rstat);
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/* is_last_iter == TRUE ==> it is the last iteration before branching. However
 * if a solutiion is found by the following heuristics, then this may no
 * longer be a last call. */
int is_feasible_u(lp_prob *p, char branching, char is_last_iter)
{
#ifndef COMPILE_IN_LP
   int s_bufid;
#endif
   int user_res;
   int feasible = IP_INFEASIBLE;
   double true_objval = p->lp_data->objval;
   LPdata *lp_data = p->lp_data;
   double lpetol = lp_data->lpetol;
   double lpetol100 = lpetol*100, lpetol1 = 1 - lpetol100;
   int *indices;
   double *values, valuesi, *heur_solution = NULL, *col_sol = NULL,
          new_obj_val;
   int cnt, i, termcode;
   var_desc **vars = lp_data->vars;
   char found_better_solution;
   int should_call_fp = FALSE;
   double *x;
   int n = lp_data->n;
   double gran_round;

   get_x(lp_data); /* maybe just fractional -- parameter ??? */

   indices = lp_data->tmp.i1; /* n */
   values = lp_data->tmp.d; /* n */

   char do_local_search = TRUE;
   double d_gap;
   heur_solution = p->lp_data->heur_solution;

#ifdef USE_SYM_APPLICATION
   cnt = collect_nonzeros(p, lp_data->x, indices, values);
   user_res = user_is_feasible(p->user, lpetol, cnt, indices, values,
			       &feasible, &true_objval, branching,
			       heur_solution);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){
    case USER_ERROR: /* Error. Consider as feasibility not recognized. */
      return(FALSE);
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    case USER_DEFAULT: /* set the default */
      user_res = TEST_INTEGRALITY;
      if (feasible != IP_INFEASIBLE){
	 printf("Warning: User set feasibility status of solution, but\n");
	 printf("SYMPHONY to check feasibility. Ignoring request.");
	 user_res = USER_SUCCESS;
      }
      break;
    default:
      break;
   }

   switch (user_res){
    case TEST_ZERO_ONE: /* User wants us to test 0/1 -ness. */
      cnt = collect_nonzeros(p, lp_data->x, indices, values);
       for (i=cnt-1; i>=0; i--){
	 if (!vars[indices[i]]->is_int)
	    continue; /* Not an integer variable */
	 if (values[i] < lpetol1) break;
       }
      feasible = i < 0 ? IP_FEASIBLE : IP_INFEASIBLE;
      break;
    case TEST_INTEGRALITY:
      x = lp_data->x;
      for (i = n - 1; i >= 0; i--){
	 if (!vars[i]->is_int)
	    continue; /* Not an integer variable */
	 valuesi = x[i];
	 if (valuesi-floor(valuesi) > lpetol100 &&
	     ceil(valuesi)-valuesi > lpetol100 &&
             valuesi>vars[i]->lb-lpetol && valuesi<vars[i]->ub+lpetol){
	    break;
	 }
      }
      feasible = i < 0 ? IP_FEASIBLE : IP_INFEASIBLE;
      break;
    default:
      break;
   }

#ifdef COMPILE_IN_LP
   /* try rounding first */
   if (user_res == TEST_INTEGRALITY && feasible != IP_FEASIBLE && feasible != IP_HEUR_FEASIBLE &&
       p->par.do_primal_heuristic && !p->par.multi_criteria){

      if (feasible == IP_INFEASIBLE){
	 true_objval = SYM_INFINITY;
      }

      if (p->has_ub){
	 d_gap = (p->ub-p->lp_data->objval)/(fabs(p->ub)+0.0001)*100;
	 if(d_gap > 0.0001){
	    true_objval = p->ub;
	    if (round_solution(p, &true_objval, heur_solution)){
	       feasible = IP_HEUR_FEASIBLE;
	    }
	 }
	 do_local_search = FALSE;
	 if(do_local_search && p->mip->n - p->mip->mip_inf->cont_var_num < 12500){// &&
	    //	    p->bc_level <=10){
	    if(feasible == IP_HEUR_FEASIBLE){
	       if((true_objval - p->lp_data->objval)/
		  (fabs(true_objval)+0.0001)*100 - d_gap < 0.0133){
		  col_sol = (double *)calloc(DSIZE, lp_data->n);
		  memcpy(col_sol, heur_solution, DSIZE*lp_data->n);
		  do_local_search = TRUE;
	       }
	    }else if(d_gap > p->par.fp_min_gap){
	       col_sol = (double *)calloc(DSIZE, lp_data->n);
	       for(i = 0; i< p->best_sol.xlength; i++) {
		  col_sol[p->best_sol.xind[i]] = p->best_sol.xval[i];
	       }
	       do_local_search = TRUE;
	    }
	 }
	 if(do_local_search){
	    printf("callin ls-1\n");
	    if (local_search(p, &true_objval, col_sol, heur_solution)){
	       feasible = IP_HEUR_FEASIBLE;
	       printf("found ls -1\n");
	    }
	 }
      }else{
	 if (round_solution(p, &true_objval, heur_solution)){
	    feasible = IP_HEUR_FEASIBLE;
	 }
      }
   }

   if (user_res == TEST_INTEGRALITY && feasible != IP_FEASIBLE && feasible != IP_HEUR_FEASIBLE) {
      fp_should_call_fp(p,branching,&should_call_fp,is_last_iter);
      if (should_call_fp==TRUE) {
         termcode    = feasibility_pump (p, &found_better_solution,
					 new_obj_val, heur_solution);

	 if (termcode!=FUNCTION_TERMINATED_NORMALLY) {
            PRINT(p->par.verbosity,0,("warning: feasibility pump faced some "
                     "difficulties.\n"));
         } else if (found_better_solution) {
            feasible    = IP_HEUR_FEASIBLE;
            true_objval = new_obj_val;
         }
      }
   }
#endif

   if (feasible == IP_FEASIBLE && p->par.multi_criteria){
      cnt = collect_nonzeros(p, lp_data->x, indices, values);
      if (analyze_multicriteria_solution(p, indices, values, cnt,
					 &true_objval, lpetol, branching) > 0){
	 if(feasible == IP_FEASIBLE){
	    if (p->par.mc_add_optimality_cuts || branching){
	       feasible = IP_FEASIBLE_BUT_CONTINUE;
	    }else{
	       feasible = IP_FEASIBLE;
	    }
	 }else{
	    feasible = IP_FEASIBLE;
	 }
      }
   }

   if (feasible == IP_FEASIBLE || feasible == IP_FEASIBLE_BUT_CONTINUE ||
       feasible == IP_HEUR_FEASIBLE){
      if (feasible == IP_HEUR_FEASIBLE) {
         cnt = collect_nonzeros(p, heur_solution, indices, values);
      } else {
         cnt = collect_nonzeros(p, lp_data->x, indices, values);
      }
      gran_round = p->par.granularity;
      gran_round = floor(gran_round + 0.5);
      if (p->par.granularity > lpetol100 &&
            fabs(gran_round-p->par.granularity) < lpetol100) {
         /* we have granularity. symphony now uses granularity to set ub on
          * lp-solver using granularity. so we round the solution to the
          * nearest integer so that this tighter ub does not cut off other
          * good solutions.
          */
         true_objval = floor(true_objval+0.5);
      }
      /* Send the solution value to the treemanager */
      if (p->has_ub && true_objval >= p->ub - p->par.granularity){
	 //FREE(heur_solution);
	 //FREE(col_sol);
	 if (!p->par.multi_criteria){
	    PRINT(p->par.verbosity, 0,
		  ("\n* Found Another Feasible Solution.\n"));
	    if (p->mip->obj_sense == SYM_MAXIMIZE){
	       PRINT(p->par.verbosity, 0, ("* Cost: %f\n\n", -true_objval
			+ p->mip->obj_offset));
	    }else{
	       PRINT(p->par.verbosity, 0, ("****** Cost: %f\n\n", true_objval
			+ p->mip->obj_offset));
	    }
	 }
	 return(feasible);
      }
      p->has_ub = TRUE;
      p->ub = true_objval;
      if (p->par.set_obj_upper_lim) {
	 set_obj_upper_lim(p->lp_data, p->ub - p->par.granularity + lpetol);
      }
      if (!p->par.multi_criteria){
	 p->best_sol.xlevel = p->bc_level;
	 p->best_sol.xindex = p->bc_index;
	 p->best_sol.xiter_num = p->iter_num;
	 p->best_sol.xlength = cnt;
	 p->best_sol.lpetol = lpetol;
	 p->best_sol.objval = true_objval;
	 FREE(p->best_sol.xind);
	 FREE(p->best_sol.xval);
	 if(cnt){
	    p->best_sol.xind = (int *) malloc(cnt*ISIZE);
	    p->best_sol.xval = (double *) malloc(cnt*DSIZE);
	    memcpy((char *)p->best_sol.xind, (char *)indices, cnt*ISIZE);
	    memcpy((char *)p->best_sol.xval, (char *)values, cnt*DSIZE);
	 }
	 if(!p->best_sol.has_sol)
	    p->best_sol.has_sol = TRUE;
	 PRINT(p->par.verbosity, 0,
	       ("\n****** Found Better Feasible Solution !\n"));
	 if (feasible == IP_HEUR_FEASIBLE){
	    PRINT(p->par.verbosity, 2,
		  ("****** After Calling Heuristics !\n"));
	 }
	 if (p->mip->obj_sense == SYM_MAXIMIZE){
	    PRINT(p->par.verbosity, 1, ("****** Cost: %f\n\n", -true_objval
					+ p->mip->obj_offset));
	 }else{
	    PRINT(p->par.verbosity, 1, ("****** Cost: %f\n\n", true_objval
					+ p->mip->obj_offset));
	 }
      }
#ifdef COMPILE_IN_LP
#pragma omp critical (new_ub)
      {
	 install_new_ub(p->tm, p->ub, p->proc_index, p->bc_index, branching,
			feasible);
	 if (p->bc_index>0) {
	    tighten_root_bounds(p);
	 }
      }
      if (!p->par.multi_criteria){
	 display_lp_solution_u(p, DISP_FEAS_SOLUTION);
      }
#else
      s_bufid = init_send(DataInPlace);
      send_dbl_array(&true_objval, 1);
      send_int_array(&(p->bc_index), 1);
      send_int_array(&feasible, 1);
      send_char_array(&branching, 1);
      send_msg(p->tree_manager, UPPER_BOUND);
      freebuf(s_bufid);
#endif
#if !defined(COMPILE_IN_LP) || !defined(COMPILE_IN_TM)
      send_feasible_solution_u(p, p->bc_level, p->bc_index, p->iter_num,
			       lpetol, true_objval, cnt, indices, values);
#endif
   }

   if (feasible == IP_FEASIBLE){
      lp_data->termcode = LP_OPT_FEASIBLE;
      p->lp_stat.lp_sols++;

#ifdef COMPILE_IN_LP
      sp_add_solution(p,cnt,indices,values,true_objval+p->mip->obj_offset,
            p->bc_index);
#endif
   }

#if 0
   if(is_last_iter){
      for (i=p->lp_data->mip->n-1; i>=0; i--){
	 if (vars[i]->is_int)
	    p->lp_data->si->setInteger(i);
      }
      write_mps(p->lp_data, "test");
   }
#endif
   //printf("feasible: solution = %f\n", lp_data->objval);
   FREE(col_sol);
   return(feasible);
}

/*===========================================================================*/

void send_feasible_solution_u(lp_prob *p, int xlevel, int xindex,
			      int xiter_num, double lpetol, double new_ub,
			      int cnt, int *xind, double *xval)
{
   int s_bufid, msgtag, user_res;

   /* Send to solution to the master */
   s_bufid = init_send(DataInPlace);
   send_int_array(&xlevel, 1);
   send_int_array(&xindex, 1);
   send_int_array(&xiter_num, 1);
   send_dbl_array(&lpetol, 1);
   send_dbl_array(&new_ub, 1);
   send_int_array(&cnt, 1);
   if (cnt > 0){
      send_int_array(xind, cnt);
      send_dbl_array(xval, cnt);
   }
#ifdef USE_SYM_APPLICATION
   user_res = user_send_feasible_solution(p->user, lpetol, cnt, xind, xval);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    case USER_ERROR: /* Error. Do the default */
    case USER_DEFAULT: /* set the default */
	 user_res = p->par.send_feasible_solution_default;
      break;
   }
   switch (user_res){
    case SEND_NONZEROS:
      msgtag = FEASIBLE_SOLUTION_NONZEROS;
      break;
    default: /* Otherwise the user packed it */
      msgtag = FEASIBLE_SOLUTION_USER;
      break;
   }
   send_msg(p->master, msgtag);
   freebuf(s_bufid);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_display_solution
 * that (graphically) displays the current solution.
\*===========================================================================*/

void display_lp_solution_u(lp_prob *p, int which_sol)
{
   int user_res;
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol = lp_data->lpetol;

   int number = 0;
   int i, *xind = lp_data->tmp.i1; /* n */
   double tmpd, *xval = lp_data->tmp.d; /* n */

   if (p->par.verbosity < 0) return;

   number = collect_nonzeros(p, x, xind, xval);

   /* Invoke user written function. */
#ifdef USE_SYM_APPLICATION
   user_res = user_display_lp_solution(p->user, which_sol, number, xind, xval);
#else
   user_res = USER_DEFAULT;
#endif

   switch(user_res){
    case USER_ERROR:
      /* SYMPHONY ignores error message. */
      return;
    case USER_AND_PP:
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
      return;
    case USER_DEFAULT:
      user_res = p->par.display_solution_default;
      break;
    default:
      break;
   }

   switch(user_res){
    case DISP_NOTHING:
      break;
    case DISP_NZ_INT:
      if (p->mip->colname){
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" Column names and values of nonzeros in the solution\n");
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; i++){
	    if (xind[i] == p->mip->n) continue; /* For multi-criteria */
	    printf("%8s %10.7f\n", p->mip->colname[xind[i]], xval[i]);
	 }
	 printf("\n");
      }else{
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" User indices and values of nonzeros in the solution\n");
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; i++){
	    if (xind[i] == p->mip->n) continue; /* For multi-criteria */
	    printf("%7d %10.7f\n", xind[i], xval[i]);
	 }
	 printf("\n");
      }
      break;
    case DISP_NZ_HEXA:
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf(" User indices (hexa) and values of nonzeros in the solution\n");
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      for (i = 0; i < number; i++){
	 if (xind[i] == p->mip->n) continue; /* For multi-criteria */
	 printf("%7x %10.7f ", xind[i], xval[i]);
	 if (!(++i & 3)) printf("\n"); /* new line after every four pair*/
      }
      printf("\n");
      break;
    case DISP_FRAC_INT:
      if (p->mip->colname){
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" Column names and values of fractional vars in solution\n");
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; i++){
	    if (xind[i] == p->mip->n) continue; /* For multi-criteria */
	    tmpd = xval[i];
	    if ((tmpd > floor(tmpd)+lpetol) && (tmpd < ceil(tmpd)-lpetol)){
	       printf("%8s %10.7f\n", p->mip->colname[xind[i]], tmpd);
	    }
	 }
	 printf("\n");
      }else{
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 printf(" User indices and values of fractional vars in solution\n");
	 printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 for (i = 0; i < number; i++){
	    if (xind[i] == p->mip->n) continue; /* For multi-criteria */
	    tmpd = xval[i];
	    if ((tmpd > floor(tmpd)+lpetol) && (tmpd < ceil(tmpd)-lpetol)){
	       printf("%7d %10.7f ", xind[i], tmpd);
	       if (!(++i & 3)) printf("\n"); /* new line after every four*/
	    }
	 }
      }
      printf("\n");
      break;
    case DISP_FRAC_HEXA:
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf(" User indices (hexa) and values of frac vars in the solution\n");
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      for (i = 0; i < number; i++){
	 if (xind[i] == p->mip->n) continue; /* For multi-criteria */
	 tmpd = xval[i];
	 if ((tmpd > floor(tmpd)+lpetol) && (tmpd < ceil(tmpd)-lpetol)){
	    printf("%7x %10.7f ", xind[i], tmpd);
	    if (!(++i & 3)) printf("\n"); /* new line after every four pair*/
	 }
      }
      printf("\n");
      break;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_branch that selects
 * candidates to branch on. It receives a number of arguments:
 *  sim_num : slacks in matrix number (the
\*===========================================================================*/

int select_candidates_u(lp_prob *p, int *cuts, int *new_vars,
			int *cand_num, branch_obj ***candidates)
{
   int user_res, action = USER__BRANCH_IF_MUST;
   LPdata *lp_data = p->lp_data;
   row_data *rows = lp_data->rows;
   int i, j = 0, m = lp_data->m;
   int *candidate_rows;
   branch_obj *can;
   cut_data **slacks_in_matrix = NULL; /* just to keep gcc quiet */

   /* If the user might need to generate rows, we better have the
    * columns COLIND_ORDERED */
   colind_sort_extra(p);

   candidate_rows = lp_data->tmp.i2; /* m */
   if (p->par.branch_on_cuts){
      slacks_in_matrix = (cut_data **)lp_data->tmp.p2; /* m */
      /* get a list of row that are candidates for branching */
      for (i=0; i<m; i++){ /* can't branch on original rows */
	 if ((rows[i].cut->branch & CANDIDATE_FOR_BRANCH)){
	    slacks_in_matrix[j] = rows[i].cut;
	    candidate_rows[j++] = i;
	 }
      }
   }

   /* First decide if we are going to branch or not */
#ifdef USE_SYM_APPLICATION
   user_res = user_shall_we_branch(p->user, lp_data->lpetol, *cuts, j,
				   slacks_in_matrix, p->slack_cut_num,
				   p->slack_cuts, lp_data->n, lp_data->vars,
				   lp_data->x, lp_data->status, cand_num,
				   candidates, &action);
   check_tailoff(p);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    case USER_ERROR:   /* In case of error, default is used. */
    case USER_DEFAULT:
      action = p->par.shall_we_branch_default;
      break;
   }

   if (p->bc_level <= p->par.load_balance_level &&
       p->node_iter_num >= p->par.load_balance_iterations)
      action = USER__DO_BRANCH;

   if ((action == USER__DO_NOT_BRANCH) ||
       (p->bound_changes_in_iter>0) ||
       (action == USER__BRANCH_IF_TAILOFF && *cuts > 0 && p->has_tailoff==FALSE) ||
       (action == USER__BRANCH_IF_MUST && *cuts > 0))
      return(DO_NOT_BRANCH);

   p->comp_times.strong_branching += used_time(&p->tt);
   {
      /* it seems we are going to branch. Before doing that, we should invoke
       * heuristics. */
      double oldobj = (p->has_ub ? p->ub : SYM_INFINITY);
      int feas_status = is_feasible_u(p, FALSE, TRUE);
      p->comp_times.primal_heur += used_time(&p->tt);
      if (feas_status == IP_FEASIBLE || (feas_status==IP_HEUR_FEASIBLE &&
					 p->ub < oldobj - lp_data->lpetol)){// && //){
	  //					 *cuts > 0)){
         return(DO_NOT_BRANCH__FEAS_SOL);
      }
   }


   action = col_gen_before_branch(p, new_vars);
   /* vars might have been added, so tmp arrays might be freed/malloc'd,
      but only those where maxn plays any role in the size. Therefore tmp.i2
      and tmp.p2 do NOT change. Phew... */

   if (action == DO_NOT_BRANCH__FATHOMED)
      return(DO_NOT_BRANCH__FATHOMED);

   /* In the other two cases we may have to re-generate the rows
      corresponding to slacks not in the matrix (whether violated
      slacks or branching candidates), depending on new_vars */
   if (*new_vars > 0 && *cand_num > 0){
      cut_data **regen_cuts = (cut_data **) malloc(*cand_num*sizeof(cut_data));
      for (j = 0, i = 0; i < *cand_num; i++){
	 can = (*candidates)[i];
	 if (can->type == VIOLATED_SLACK ||
	     can->type == CANDIDATE_CUT_NOT_IN_MATRIX){
	    regen_cuts[j++] = can->row->cut;
	 }
      }
      if (j > 0){
	 int new_row_num;
	 waiting_row **new_rows;
	 unpack_cuts_u(p, CUT_FROM_TM, UNPACK_CUTS_SINGLE,
		       j, regen_cuts, &new_row_num, &new_rows);
	 for (j = 0, i = 0; i < *cand_num; i++){
	    can = (*candidates)[i];
	    if (can->type == VIOLATED_SLACK ||
		can->type == CANDIDATE_CUT_NOT_IN_MATRIX){
	       free_waiting_row(&can->row);
	       can->row = new_rows[j++];
	    }
	 }
	 FREE(new_rows);
      }
      FREE(regen_cuts);
   }

   if (action == DO_NOT_BRANCH)
      return(DO_NOT_BRANCH);

   /* So the action from col_gen_before_branch is DO_BRANCH */

   action = USER__DO_BRANCH;

   /* before branching, update the control parameters for cut generation
    * --asm4
    */
   //   if (p->bc_level==0) {
   //update_cut_parameters(p);
   // }

   /* OK, so we got to branch */
#ifdef USE_SYM_APPLICATION
   user_res = user_select_candidates(p->user, lp_data->lpetol, *cuts, j,
				     slacks_in_matrix, p->slack_cut_num,
				     p->slack_cuts, lp_data->n, lp_data->vars,
				     lp_data->x, lp_data->status, cand_num,
				     candidates, &action, p->bc_level);
#else
   user_res = USER_DEFAULT;
#endif

   /* Get rid of any contsraint from slack_cuts which is listed in candidates
    * and rewrite the position of the CANDIDATE_CUT_IN_MATRIX ones */
   if (p->par.branch_on_cuts){
      for (i = 0; i < *cand_num; ){
	 can = (*candidates)[i];
	 switch (can->type){
	  case CANDIDATE_VARIABLE:
	    i++;
	    break;
	  case CANDIDATE_CUT_IN_MATRIX:
	    can->position = candidate_rows[can->position];
	    i++;
	    break;
	  case VIOLATED_SLACK:
	  case CANDIDATE_CUT_NOT_IN_MATRIX:
	    free_cut(p->slack_cuts + can->position);
	    i++;
	    break;
	  case SLACK_TO_BE_DISCARDED:
	    free_cut(p->slack_cuts + can->position);
	    free_candidate(*candidates + i);
	    (*candidates)[i] = (*candidates)[--(*cand_num)];
	    break;
	 }
      }
      compress_slack_cuts(p);
   }

   if (action == USER__DO_NOT_BRANCH)
      return(DO_NOT_BRANCH);

   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      if (! *cand_num){
	 printf("Error! User didn't select branching candidates!\n");
	 return(ERROR__NO_BRANCHING_CANDIDATE);
      }
      return(DO_BRANCH);
    case USER_ERROR:    /* In case of error, default is used. */
    case USER_DEFAULT:
      user_res = p->par.select_candidates_default;
      break;
    default:
      break;
   }

   i = (int) (p->par.strong_branching_cand_num_max -
	      p->par.strong_branching_red_ratio * p->bc_level);
   i = MAX(i, p->par.strong_branching_cand_num_min);

   switch(user_res){
    case USER__CLOSE_TO_HALF:
      branch_close_to_half(p, i, cand_num, candidates);
      break;
    case USER__CLOSE_TO_HALF_AND_EXPENSIVE:
      branch_close_to_half_and_expensive(p, i, cand_num, candidates);
      break;
    case USER__CLOSE_TO_ONE_AND_CHEAP:
      branch_close_to_one_and_cheap(p, i, cand_num, candidates);
      break;

    default:
      /* Unexpected return value. Do something!! */
      break;
   }

   if (! *cand_num){
      PRINT(p->par.verbosity, 2,
	    ("No branching candidates found using default rule...\n"));
      return(DO_NOT_BRANCH);
   }
   return(DO_BRANCH);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_compare_candidates
 * that compares to branching candidates.
\*===========================================================================*/

int compare_candidates_u(lp_prob *p, double oldobjval,
			 branch_obj *best, branch_obj *can)
{
   int user_res;
   int i;
   double low0, low1, high0, high1;
   double lpetol = p->lp_data->lpetol;
   const double ub_minus_gran = p->ub - p->par.granularity;
   const double alpha = p->par.strong_branching_high_low_weight;
   const double infinity = SYM_INFINITY;
#ifdef COMPILE_FRAC_BRANCHING
   int frl0, frl1, frh0, frh1;
#endif
   for (i = can->child_num-1; i >= 0; i--){
      switch (can->termcode[i]){
       case LP_OPTIMAL:
       case LP_OPT_FEASIBLE_BUT_CONTINUE:
#ifdef DO_TESTS
	 if (can->objval[i] < oldobjval - .01){
	    printf("#####Error: Branching candidate has lower objval ");
	    printf("(%.3f) than parent (%.3f)\n", can->objval[i],  oldobjval);
	 }
#endif
	 break;
       case LP_OPT_FEASIBLE:
       case LP_D_UNBOUNDED:
       case LP_D_OBJLIM:
	 can->objval[i] = MAXDOUBLE / 2;
	 break;
       case LP_D_ITLIM:
	 can->objval[i] = MAX(can->objval[i], oldobjval);
	 break;
       case LP_D_INFEASIBLE:
       case LP_ABANDONED:
	 can->objval[i] = oldobjval;
	 break;
      }
   }

   /*------------------------------------------------------------------------*\
    * If ALL descendants in cand terminated with primal infeasibility
    * or high cost, that proves that the current node can be fathomed,
    * so we select cand and force branching on it.
    *
    * MAYBE THIS SHOULD BE LEFT TO THE USER ?????????????????
   \*------------------------------------------------------------------------*/

   for (i = can->child_num-1; i >= 0; i--){
      if (! (can->termcode[i] == LP_D_UNBOUNDED ||
	     can->termcode[i] == LP_D_OBJLIM ||
	     can->termcode[i] == LP_OPT_FEASIBLE ||
	     can->termcode[i] == LP_OPT_FEASIBLE_BUT_CONTINUE ||
	     (can->termcode[i] == LP_OPTIMAL && p->has_ub &&
	      can->objval[i] > ub_minus_gran))){
	 break;
      }
   }

   if (i < 0){
      /* i.e., we did not break, i.e., we'll select this cand */
      return(SECOND_CANDIDATE_BETTER_AND_BRANCH_ON_IT);
   }

   /* If this is the first, keep it */
   if (best == NULL){
      return(SECOND_CANDIDATE_BETTER);
   }

   /* Otherwise, first give the choice to the user */
#ifdef USE_SYM_APPLICATION
   user_res = user_compare_candidates(p->user, best, can, p->ub,
				      p->par.granularity, &i);

#else
   user_res = USER_DEFAULT;
#endif

   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
       /* User function terminated without problems. No post-processing. */
      return(i);
    case USER_ERROR:
      /* In case of error, default is used. */
    case USER_DEFAULT:
      user_res = p->par.compare_candidates_default;
      break;
    default:
      break;
   }

   /* Well, the user let us make the choice.
    *
    * If something had gone wrong with at least one descendant in
    * can, then prefer to choose something else. */
   for (i = can->child_num-1; i >= 0; i--)
      if (can->termcode[i] == LP_ABANDONED)
	 return(FIRST_CANDIDATE_BETTER);

   /* OK, so all descendants in can finished fine. Just do whatever
    * built-in was asked */
#ifdef COMPILE_FRAC_BRANCHING
   for (frl0 = frh0 = best->frac_num[0], i = best->child_num-1; i; i--){
      frl0 = MIN(frl0, best->frac_num[i]);
      frh0 = MAX(frh0, best->frac_num[i]);
   }
   for (frl1 = frh1 = can->frac_num[0], i = can->child_num-1; i; i--){
      frl1 = MIN(frl1, can->frac_num[i]);
      frh1 = MAX(frh1, can->frac_num[i]);
   }
#endif
   for (low0 = high0 = best->objval[0], i = best->child_num-1; i; i--){
      low0 = MIN(low0, best->objval[i]);
      high0 = MAX(high0, best->objval[i]);
   }
   for (low1 = high1 = can->objval[0], i = can->child_num-1; i; i--){
      low1 = MIN(low1, can->objval[i]);
      high1 = MAX(high1, can->objval[i]);
   }

   switch(user_res){
    case HIGH_LOW_COMBINATION:
      if (high0 > ub_minus_gran) {
         high0 = infinity;
      }
      if (low0 > ub_minus_gran) {
         low0 = infinity;
      }
      if (high1 > ub_minus_gran) {
         high1 = infinity;
      }
      if (low1 > ub_minus_gran) {
         low1 = infinity;
      }
      i = (alpha*low0 + (1 - alpha)*high0 > alpha*low1 + (1 - alpha)*high1) ?
         0 : 1;
      break;
    case BIGGEST_DIFFERENCE_OBJ:
      i = (high0 - low0 >= high1 - low1) ? 0 : 1;
      break;
    case LOWEST_LOW_OBJ:
      i = (fabs(low0-low1)<lpetol) ? (high0 <= high1 ? 0 : 1) : (low0 < low1 ? 0 : 1);
      break;
    case HIGHEST_LOW_OBJ:
      i = (fabs(low0-low1)<lpetol) ? (high0 >= high1 ? 0 : 1) : (low0 > low1 ? 0 : 1);
      break;
    case LOWEST_HIGH_OBJ:
      i = (fabs(high0-high1)<lpetol) ? (low0 <= low1 ? 0 : 1) : (high0 < high1 ? 0 : 1);
      break;
    case HIGHEST_HIGH_OBJ:
      i = (fabs(high0-high1)<lpetol) ? (low0 >= low1 ? 0 : 1) : (high0 > high1 ? 0 : 1);
      break;
#ifdef COMPILE_FRAC_BRANCHING
    case HIGHEST_LOW_FRAC:
      i = (frl0 == frl1) ? (frh0 >= frh1 ? 0 : 1) : (frl0 > frl1 ? 0 : 1);
      break;
    case LOWEST_LOW_FRAC:
      i = (frl0 == frl1) ? (frh0 <= frh1 ? 0 : 1) : (frl0 < frl1 ? 0 : 1);
      break;
    case HIGHEST_HIGH_FRAC:
      i = (frh0 == frh1) ? (frl0 >= frl1 ? 0 : 1) : (frh0 > frh1 ? 0 : 1);
      break;
    case LOWEST_HIGH_FRAC:
      i = (frh0 == frh1) ? (frl0 <= frl1 ? 0 : 1) : (frh0 < frh1 ? 0 : 1);
      break;
#endif
    default: /* Unexpected return value. Do something!! */
      break;
   }
   return(i == 0 ? FIRST_CANDIDATE_BETTER : SECOND_CANDIDATE_BETTER);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function invokes the user written function user_select_child that
 * selects one of the candidates after branching for further processing.
\*===========================================================================*/

int select_child_u(lp_prob *p, branch_obj *can, char *action)
{
   int user_res;
   int ind, i;

#ifdef DO_TESTS
   char sense;
   for (i = can->child_num-1; i >= 0; i--){
      sense = can->sense[i];
      if (sense != 'E' && sense != 'L' && sense != 'G' && sense != 'R'){
	 printf("Error! The sense of a child doesn't make sense!");
	 printf("(nonexistent)\n\n");
	 return(ERROR__ILLEGAL_BRANCHING);
      }
   }
#endif

   for (ind = -1, i = 0; i < can->child_num; i++){
      action[i] = RETURN_THIS_CHILD;
      if (p->lp_data->nf_status == NF_CHECK_NOTHING){
	 /*see which one is infeasible!*/
	 if (can->termcode[i] == LP_OPTIMAL ||
	     can->termcode[i] == LP_D_ITLIM){
	    if (p->has_ub &&
		can->objval[i] > p->ub - p->par.granularity){
	       action[i] = PRUNE_THIS_CHILD_FATHOMABLE;
	    }
	 }else if (can->termcode[i] == LP_OPT_FEASIBLE ||
		   can->termcode[i] == LP_OPT_FEASIBLE_BUT_CONTINUE){
	    action[i] = PRUNE_THIS_CHILD_FATHOMABLE;
	 }else{
	    action[i] = PRUNE_THIS_CHILD_INFEASIBLE;
	 }
      }
   }

#ifdef USE_SYM_APPLICATION
   user_res = user_select_child(p->user, p->ub, can, action);
#else
   user_res = USER_DEFAULT;
#endif

   switch(user_res){
    case USER_NO_PP:
    case USER_AND_PP:
      /* User function terminated without problems. Skip post-processing. */
      break;
    case USER_ERROR:
      /* In case of error, default is used. */
    case USER_DEFAULT:
      user_res = p->par.select_child_default;
      break;
    default:
      break;
   }

   switch(user_res){
    case PREFER_LOWER_OBJ_VALUE:
      for (ind = 0, i = can->child_num-1; i; i--){
	 if (can->objval[i] < can->objval[ind])
	    ind = i;
      }
      if (!p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that if the lowest objval child is fathomed then everything is */
      break;

    case PREFER_HIGHER_OBJ_VALUE:
      for (ind = 0, i = can->child_num-1; i; i--){
	 if ((can->objval[i] > can->objval[ind]) &&
	     (! p->has_ub ||
	      (p->has_ub && can->objval[i] < p->ub - p->par.granularity)))
	    ind = i;
      }
      if (! p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that this selects the highest objval child NOT FATHOMED, thus
       * if the highest objval child is fathomed then so is everything */
      break;

#ifdef COMPILE_FRAC_BRANCHING
    case PREFER_MORE_FRACTIONAL:
      for (ind = 0, i = can->child_num-1; i; i--){
	 if ((can->frac_num[i] > can->frac_num[ind]) &&
	     (! p->has_ub ||
	      (p->has_ub && can->objval[i] < p->ub - p->par.granularity)))
	    ind = i;
      }
      if (! p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that this selects the most fractional child NOT FATHOMED, thus
       * if that child is fathomed then so is everything */
      break;

    case PREFER_LESS_FRACTIONAL:
      for (ind = 0, i = can->child_num-1; i; i--){
	 if ((can->frac_num[i] < can->frac_num[ind]) &&
	     (! p->has_ub ||
	      (p->has_ub && can->objval[i] < p->ub - p->par.granularity)))
	    ind = i;
      }
      if (! p->has_ub ||
	  (p->has_ub && can->objval[ind] < p->ub - p->par.granularity))
	 action[ind] = KEEP_THIS_CHILD;
      /* Note that this selects the least fractional child NOT FATHOMED, thus
       * if that child is fathomed then so is everything */
      break;
#endif

    case USER_NO_PP:
    case USER_AND_PP:
      break;

    default:
      /* Unexpected return value. Do something!! */
      break;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function prints whatever statistics we want on branching
\*===========================================================================*/

void print_branch_stat_u(lp_prob *p, branch_obj *can, char *action)
{
   int i;

   if (can->type == CANDIDATE_VARIABLE){
      if (p->mip){
	 if (p->mip->colname){
	    printf("Branching on variable %s \n   children: ",
		   p->mip->colname[p->lp_data->vars[can->position]->userind]);
	 }
      }else{
	 printf("Branching on variable %i ( %i )\n   children: ",
		can->position, p->lp_data->vars[can->position]->userind);
      }
   }else{ /* must be CANDIDATE_CUT_IN_MATRIX */
      printf("Branching on a cut %i\n   children: ",
	     p->lp_data->rows[can->position].cut->name);
   }
   for (i=0; i<can->child_num; i++){
      if (can->objval[i] != MAXDOUBLE / 2){
	 if (p->mip->obj_sense == SYM_MAXIMIZE){
	    printf("[%.3f, %i,%i]  ", -can->objval[i] + p->mip->obj_offset,
		   can->termcode[i], can->iterd[i]);
	 }else{
	    printf("[%.3f, %i,%i]  ", can->objval[i] + p->mip->obj_offset,
		   can->termcode[i], can->iterd[i]);
	 }
      }else{
	 printf("[*, %i,%i]  ", can->termcode[i], can->iterd[i]);
      }
   }
   printf("\n");

#ifdef USE_SYM_APPLICATION
   if (can->type == CANDIDATE_VARIABLE){
      user_print_branch_stat(p->user, can, NULL,
			     p->lp_data->n, p->lp_data->vars, action);
   }else{
      user_print_branch_stat(p->user, can,
			     p->lp_data->rows[can->position].cut,
			     p->lp_data->n, p->lp_data->vars, action);
   }
#endif
}

/*===========================================================================*/

/*===========================================================================*\
 * Append additional information to the description of an active node
 * before it is sent back to the tree manager.
\*===========================================================================*/

void add_to_desc_u(lp_prob *p, node_desc *desc)
{
   desc->desc_size = 0;
   desc->desc = NULL;
#ifdef USE_SYM_APPLICATION
   user_add_to_desc(p->user, &desc->desc_size,
		    &desc->desc);
#endif
}

/*===========================================================================*/

int same_cuts_u(lp_prob *p, waiting_row *wrow1, waiting_row *wrow2)
{
   int user_res;
   int same_cuts = DIFFERENT_CUTS;
   cut_data *rcut1 = NULL, *rcut2 = NULL;

#ifdef USE_SYM_APPLICATION
   user_res = user_same_cuts(p->user, wrow1->cut, wrow2->cut, &same_cuts);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){
    case USER_SUCCESS:
    case USER_NO_PP:
    case USER_AND_PP:
      break;
    case USER_ERROR: /* Error. Use the default */
    case USER_DEFAULT: /* the only default is to compare byte by byte */
      rcut1 = wrow1->cut;
      rcut2 = wrow2->cut;
      if (rcut1->type != rcut2->type || rcut1->sense != rcut2->sense ||
	  rcut1->size != rcut2->size ||
	  memcmp(rcut1->coef, rcut2->coef, rcut1->size))
	 break; /* if LHS is different, then just break out. */

      /* Otherwise the two cuts have the same left hand side. Test which
       * one is stronger */
      /********* something should be done about ranged constraints ***********/
      /* FIXME! */
      if (rcut1->sense == 'L'){
	 same_cuts = rcut1->rhs > rcut2->rhs - p->lp_data->lpetol ?
	    SECOND_CUT_BETTER : FIRST_CUT_BETTER;
	 break;
      }else if (rcut1->sense == 'G'){
	 same_cuts = rcut1->rhs < rcut2->rhs + p->lp_data->lpetol ?
	    SECOND_CUT_BETTER : FIRST_CUT_BETTER;
	 break;
      }
      same_cuts = wrow1->source_pid < wrow2->source_pid ?
	 SECOND_CUT_BETTER : FIRST_CUT_BETTER;
      break;
   }

   switch(same_cuts){
    case SECOND_CUT_BETTER: /* effective replace the old with the new, then..*/
      same_cuts = SAME_CUTS;
      wrow1->violation += fabs(rcut1->rhs - rcut2->rhs);
      rcut1->rhs = rcut2->rhs;
      rcut1->name = rcut2->name;
    case SAME_CUTS:
    case FIRST_CUT_BETTER:  /* delete the new */
      FREE(rcut2->coef);
      break;

    case DIFFERENT_CUTS:
      break;
   }

   return(same_cuts);
}

/*===========================================================================*/

void unpack_cuts_u(lp_prob *p, int from, int type,
		   int cut_num, cut_data **cuts,
		   int *new_row_num, waiting_row ***new_rows)
{
   LPdata       *lp_data = p->lp_data;
   int           user_res;
   int           i, j, k, l = 0, nzcnt, real_nzcnt, explicit_row_num = 0;
   const int     n = lp_data->n;
   int          *matind, *row_matind;
   double       *matval, *row_matval;
   waiting_row **row_list = NULL;
   double       *obj1 = p->mip->obj1;
   double       *obj2 = p->mip->obj2;
   var_desc    **vars = lp_data->vars;
   const int     is_userind_in_order = p->par.is_userind_in_order;

   colind_sort_extra(p);

   if (cut_num > 0)
      row_list = (waiting_row **) calloc (cut_num, sizeof(waiting_row *));

   /* First SYMPHONY looks for cut types that it knows */
   for (i = 0; i < cut_num; i++){

      switch (cuts[i]->type){

      case EXPLICIT_ROW:
	 real_nzcnt = 0;
	 row_list[explicit_row_num] =
	    (waiting_row *) malloc(sizeof(waiting_row));
	 row_list[explicit_row_num]->cut = cuts[i];
	 nzcnt = ((int *) (cuts[i]->coef))[0];
	 matval = (double *) (cuts[i]->coef + DSIZE);
	 matind = (int *) (cuts[i]->coef + (nzcnt + 1)*DSIZE);
	 row_matval = row_list[explicit_row_num]->matval =
            (double *) malloc(nzcnt * DSIZE);
	 row_matind = row_list[explicit_row_num]->matind =
            (int *) malloc(nzcnt * ISIZE);
         if (is_userind_in_order) {
            memcpy(row_matind, matind, nzcnt*ISIZE);
            memcpy(row_matval, matval, nzcnt*DSIZE);
            real_nzcnt = nzcnt;
         } else {
            for (j = 0; j < n; j++){
               for (k = 0; k < nzcnt; k++){
                  if (matind[k] == vars[j]->userind){
                     row_matind[real_nzcnt]   = j;
                     row_matval[real_nzcnt++] = matval[k];
                  }
               }
            }
         }
	 row_list[explicit_row_num++]->nzcnt = real_nzcnt;
	 cuts[i] = NULL;
	 break;

      case OPTIMALITY_CUT_FIRST:
	 row_list[explicit_row_num] =
	    (waiting_row *) malloc(sizeof(waiting_row));
	 row_matind = row_list[explicit_row_num]->matind =
            (int *) malloc (lp_data->n * ISIZE);
	 row_matval = row_list[explicit_row_num]->matval =
	    (double *) malloc (lp_data->n * DSIZE);
	 row_list[explicit_row_num]->cut = cuts[i];
	 for (nzcnt = 0, j = 0; j < lp_data->n; j++){
	    if (vars[j]->userind == p->mip->n)
	       continue;
	    row_matind[nzcnt] = j;
	    row_matval[nzcnt++] = obj1[vars[j]->userind];
	 }
	 cuts[i]->sense = 'L';
	 cuts[i]->deletable = FALSE;
	 cuts[i]->branch = DO_NOT_BRANCH_ON_THIS_ROW;
         row_list[explicit_row_num++]->nzcnt = nzcnt;
	 cuts[i] = NULL;
	 break;

      case OPTIMALITY_CUT_SECOND:
	 row_list[explicit_row_num] =
	    (waiting_row *) malloc(sizeof(waiting_row));
	 row_list[explicit_row_num]->matind =
	    (int *) malloc (lp_data->n * ISIZE);
	 row_list[explicit_row_num]->matval =
	    (double *) malloc (lp_data->n * DSIZE);
	 row_list[explicit_row_num]->cut = cuts[i];
	 for (nzcnt = 0, j = 0; j < lp_data->n; j++){
	    if (vars[j]->userind == p->mip->n)
	       continue;
	    row_list[explicit_row_num]->matind[nzcnt] = j;
	    row_list[explicit_row_num]->matval[nzcnt++] =
	       obj2[vars[j]->userind];
	 }
	 cuts[i]->sense = 'L';
	 cuts[i]->deletable = FALSE;
	 cuts[i]->branch = DO_NOT_BRANCH_ON_THIS_ROW;
         row_list[explicit_row_num++]->nzcnt = nzcnt;
	 cuts[i] = NULL;
	 break;

       default: /* A user cut type */
#if defined(COMPILE_IN_CG) && defined(CHECK_CUT_VALIDITY)
	  check_validity_of_cut_u(p->cgp, cuts[i]);
#endif
	  if (l != i){
	     cuts[l++] = cuts[i];
	     cuts[i] = NULL;
	  }else{
	     l++;
	  }
	 break;
      }
   }

   *new_row_num = 0;

#ifdef USE_SYM_APPLICATION
   user_res = user_unpack_cuts(p->user, from, type,
			       lp_data->n, lp_data->vars,
			       l, cuts, new_row_num, new_rows);
#else
   user_res = USER_DEFAULT;
#endif

   for (i = 0; i < l; i++){
      if (cuts[i]){
	 (*new_rows)[i]->cut = cuts[i];
	 cuts[i] = NULL;
      }
   }

   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
    case USER_DEFAULT:

      /* Combine the user's rows with SYMPHONY's rows */
      if (*new_row_num == 0 && explicit_row_num == 0){
	 FREE(row_list);
	 *new_row_num = 0;
	 *new_rows = NULL;
      }else if (*new_row_num == 0 && explicit_row_num > 0){
	 *new_row_num = explicit_row_num;
	 *new_rows = row_list;
      }else if (*new_row_num > 0 && explicit_row_num > 0){
	 if (*new_row_num + explicit_row_num > cut_num){
	    row_list = (waiting_row **) realloc(row_list, *new_row_num +
						explicit_row_num);
	 }
	 for (i = explicit_row_num; i < *new_row_num + explicit_row_num; i++){
	    row_list[i] = (*new_rows)[i - explicit_row_num];
	 }

	 FREE(*new_rows);
	 *new_row_num += explicit_row_num;
	 *new_rows = row_list;
      }else{
	 FREE(row_list);
      }

      break;

    case USER_ERROR: /* Error. ??? what will happen ??? */
      *new_row_num = 0;
      FREE(*new_rows);

      break;

    default: /* No builtin possibility. Counts as ERROR. */
      break;
   }

   free_cuts(cuts, cut_num);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user packs together and sends a message to the cut generator or
 * cut pool process to obtain violated cuts.
 * Default options: SEND_NONZEROS, SEND_FRACTIONS.
 * The function return 1 or 0, depending on whether the sending of the
 * lp solution was successful or not.
\*===========================================================================*/

int send_lp_solution_u(lp_prob *p, int tid)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   int user_res, nzcnt, s_bufid, msgtag = ANYTHING;
   int *xind = lp_data->tmp.i1; /* n */
   double *xval = lp_data->tmp.d; /* n */

   s_bufid = init_send(DataInPlace);
   send_int_array(&p->bc_level, 1);
   send_int_array(&p->bc_index, 1);
   send_int_array(&p->iter_num, 1);
   send_dbl_array(&lp_data->lpetol, 1);
   if (tid == p->cut_gen){
      send_dbl_array(&lp_data->objval, 1);
      send_int_array(&p->has_ub, 1);
      if (p->has_ub)
	 send_dbl_array(&p->ub, 1);
   }
   colind_sort_extra(p);
#ifdef USE_SYM_APPLICATION
   user_res = user_send_lp_solution(p->user, lp_data->n, lp_data->vars, x,
				    tid == p->cut_gen ?
				    LP_SOL_TO_CG : LP_SOL_TO_CP);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){
    case USER_ERROR: /* Error. Consider as couldn't send to cut_gen, i.e.,
		   equivalent to NO_MORE_CUTS_FOUND */
      freebuf(s_bufid);
      return(0);
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      msgtag = LP_SOLUTION_USER;
      break;
    case USER_DEFAULT: /* set the default */
      user_res = p->par.pack_lp_solution_default; /* SEND_NONZEROS */
      break;
   }

   if (msgtag == LP_SOLUTION_USER){
      send_msg(tid, LP_SOLUTION_USER);
      freebuf(s_bufid);
      return(1);
   }

   switch(user_res){
    case SEND_NONZEROS:
      nzcnt = collect_nonzeros(p, x, xind, xval);
      msgtag = LP_SOLUTION_NONZEROS;
      break;
    case SEND_FRACTIONS:
      nzcnt = collect_fractions(p, x, xind, xval);
      msgtag = LP_SOLUTION_FRACTIONS;
      break;
   }
   /* send the data */
   send_int_array(&nzcnt, 1);
   send_int_array(xind, nzcnt);
   send_dbl_array(xval, nzcnt);
   send_msg(tid, msgtag);
   freebuf(s_bufid);

   return(1);
}

/*===========================================================================*/

void logical_fixing_u(lp_prob *p)
{
   char *status = p->lp_data->tmp.c; /* n */
   char *lpstatus = p->lp_data->status;
   char *laststat = status + p->lp_data->n;
   int  fixed_num = 0, user_res;

   colind_sort_extra(p);
   //memcpy(status, lpstatus, p->lp_data->n);
   memcpy(status, lpstatus, CSIZE*p->lp_data->n);

#ifdef USE_SYM_APPLICATION
   user_res = user_logical_fixing(p->user, p->lp_data->n, p->lp_data->vars,
				  p->lp_data->x, status, &fixed_num);
#else
   user_res = USER_DEFAULT;
#endif

   switch(user_res){
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      if (fixed_num > 0){
	 while (status != laststat) {
	    *lpstatus &= NOT_REMOVABLE;
	    *lpstatus++ |= (*status++ & (NOT_FIXED |
					 TEMP_FIXED_TO_LB | TEMP_FIXED_TO_UB |
					 PERM_FIXED_TO_LB | PERM_FIXED_TO_UB));
	 }
      }
    case USER_DEFAULT:
      break;
   }
}

/*===========================================================================*/

int generate_column_u(lp_prob *p, int lpcutnum, cut_data **cuts,
		      int prevind, int nextind, int generate_what,
		      double *colval, int *colind, int *collen, double *obj,
		      double *lb, double *ub)
{
   int real_nextind = nextind;
#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_generate_column(p->user, generate_what,
					    p->lp_data->m - p->base.cutnum,
					    cuts, prevind, nextind,
					    &real_nextind, colval, colind,
					    collen, obj, lb, ub) );
#endif
   return(real_nextind);
}

/*===========================================================================*/

int generate_cuts_in_lp_u(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   int user_res, new_row_num = 0;
   waiting_row **new_rows = NULL;
   char deleted_cut;
   cut_data **cuts = NULL;
   int i, j;

   colind_sort_extra(p);

#if defined(COMPILE_IN_CG) || defined(COMPILE_IN_CP)
   {
#ifdef COMPILE_IN_CP
      int cp_new_row_num = 0;
      waiting_row **cp_new_rows = NULL;
#endif
#ifdef COMPILE_IN_CG
      int cg_new_row_num = 0;
      waiting_row **cg_new_rows = NULL;
#endif
      int user_res2, xlength = 0, *xind = NULL;
      lp_sol *cur_sol = &(p->cgp->cur_sol);
      double *xval = NULL, lpetol = 0;

#ifdef USE_SYM_APPLICATION
      user_res2 = user_send_lp_solution(p->user,
					lp_data->n, lp_data->vars, x,
					LP_SOL_WITHIN_LP);
#else
      user_res2 = USER_DEFAULT;
#endif

      if (user_res2 == USER_DEFAULT)
	 user_res2 = p->par.pack_lp_solution_default;

      switch (user_res2){
       case USER_ERROR:
	 return(ERROR__USER);
       case USER_SUCCESS:
       case USER_AND_PP:
       case USER_NO_PP:
	 break;
       case SEND_NONZEROS:
       case SEND_FRACTIONS:
	 cur_sol->xind = xind = lp_data->tmp.i1; /* n */
	 cur_sol->xval = xval = lp_data->tmp.d; /* n */
	 cur_sol->lpetol = lpetol = lp_data->lpetol;
	 cur_sol->xlevel = p->bc_level;
	 cur_sol->xindex = p->bc_index;
	 cur_sol->xiter_num = p->iter_num;
	 cur_sol->objval = lp_data->objval;
	 if (p->has_ub)
	    p->cgp->ub = p->ub;
	 cur_sol->xlength = xlength = user_res2 == SEND_NONZEROS ?
	                             collect_nonzeros(p, x, xind, xval) :
		                     collect_fractions(p, x, xind, xval);
	 break;
      }
#ifdef COMPILE_IN_CG
      if (p->cgp->par.do_findcuts && !new_row_num)
	 find_cuts_u(p->cgp, p->lp_data, &cg_new_row_num);
#endif

      if (p->cgp->cuts_to_add_num){
	 unpack_cuts_u(p, CUT_FROM_CG, UNPACK_CUTS_MULTIPLE,
		       p->cgp->cuts_to_add_num, p->cgp->cuts_to_add,
		       &cg_new_row_num, &cg_new_rows);
	 p->cgp->cuts_to_add_num = 0;
	 if (cg_new_row_num){
	    for (i = 0; i < cg_new_row_num; i++){
	       if (cg_new_rows[i]->cut->name != CUT__SEND_TO_CP)
		  cg_new_rows[i]->cut->name = CUT__DO_NOT_SEND_TO_CP;
	       cg_new_rows[i]->source_pid = INTERNAL_CUT_GEN;
	       for (j = p->waiting_row_num - 1; j >= 0; j--){
		  if (same_cuts_u(p, p->waiting_rows[j],
				  cg_new_rows[i]) !=
		      DIFFERENT_CUTS){
		     free_waiting_row(cg_new_rows+i);
		     break;
		  }
	       }
	       if (j < 0){
		  add_new_rows_to_waiting_rows(p, cg_new_rows+i, 1);
	       }
	    }
	    FREE(cg_new_rows);
	 }
      }
#if defined(COMPILE_IN_CP) && defined(COMPILE_IN_LP)

      if ((p->iter_num == 1 && (p->bc_level > 0 || p->phase==1)) ||
	  (p->iter_num % p->par.cut_pool_check_freq == 0) ||
	  (!cg_new_row_num)){
	 cut_pool *cp = p->tm->cpp[p->cut_pool];
	 p->comp_times.separation += used_time(&p->tt);
	 cur_sol->lp = 0;
#pragma omp critical(cut_pool)
	 if (cp){
	    cp_new_row_num = check_cuts_u(cp, cur_sol);
	    if (++cp->reorder_count % 10 == 0){
	       delete_duplicate_cuts(cp);
	       order_cuts_by_quality(cp);
	       cp->reorder_count = 0;
	    }
	    if (cp_new_row_num){
	       unpack_cuts_u(p, CUT_FROM_CG, UNPACK_CUTS_MULTIPLE,
			     cp->cuts_to_add_num, cp->cuts_to_add,
			     &cp_new_row_num, &cp_new_rows);
	       cp->cuts_to_add_num = 0;
	    }
	 }
	 if (cp_new_row_num){
	    for (i = 0; i < cp_new_row_num; i++){
	       if (cp_new_rows[i]->cut->name != CUT__SEND_TO_CP)
		  cp_new_rows[i]->cut->name = CUT__DO_NOT_SEND_TO_CP;
	       cp_new_rows[i]->source_pid = INTERNAL_CUT_POOL;
	       for (j = p->waiting_row_num - 1; j >= 0; j--){
		  if (same_cuts_u(p, p->waiting_rows[j], cp_new_rows[i]) !=
		      DIFFERENT_CUTS){
		     free_waiting_row(cp_new_rows+i);
		     break;
		  }
	       }
	       if (j < 0){
		  add_new_rows_to_waiting_rows(p, cp_new_rows+i, 1);
	       }
	    }
	    FREE(cp_new_rows);
	 }
	 p->comp_times.cut_pool += used_time(&p->tt);
      }
#endif
   }
#endif

#ifdef USE_SYM_APPLICATION
   user_res = user_generate_cuts_in_lp(p->user, lp_data, lp_data->n,
				       lp_data->vars, x,
				       &new_row_num, &cuts);
#else
   user_res = GENERATE_CGL_CUTS;
#endif

   switch(user_res){
    case USER_ERROR:
      FREE(cuts);
      return(ERROR__USER);
    case GENERATE_CGL_CUTS:
    case USER_DEFAULT:
      /* Add to the user's list of cuts */
#ifdef USE_CGL_CUTS
      if (p->par.cgl.generate_cgl_cuts){
         int bound_changes = 0;
         /*
         double ub = p->has_ub ? p->ub : SYM_INFINITY;
	 generate_cgl_cuts(lp_data, &new_row_num, &cuts, FALSE,
			   p->bc_index, p->bc_level, p->node_iter_num,
                            p->par.max_cut_num_per_iter_root,
                           ub, &bound_changes, &(p->lp_stat), &(p->comp_times),
                           p->par.verbosity);
                           */
	 generate_cgl_cuts_new(p, &new_row_num, &cuts, FALSE,
			       &bound_changes);
	 if (bound_changes>0) {
	    p->bound_changes_in_iter += bound_changes;
	 }
	 // if(p->bc_index < 1 && p->iter_num == 1 ){
	 //  p->par.cgl = lp_data->cgl;
	 // }
      }
#endif
      /* Fall through to next case */

    case DO_NOT_GENERATE_CGL_CUTS:
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      /* Process the generated cuts */
      if (new_row_num){
	 unpack_cuts_u(p, CUT_FROM_CG, UNPACK_CUTS_MULTIPLE,
		       new_row_num, cuts, &new_row_num, &new_rows);
	 for (i = 0; i < new_row_num; i++){
	    if (new_rows[i]->cut->name != CUT__SEND_TO_CP)
	       new_rows[i]->cut->name = CUT__DO_NOT_SEND_TO_CP;
	    new_rows[i]->source_pid = INTERNAL_CUT_GEN;
	 }
      }
      /* Test whether any of the new cuts are identical to any of
         the old ones. */
      if (p->waiting_row_num && new_row_num){
	 for (i = 0, deleted_cut = FALSE; i < new_row_num;
	      deleted_cut = FALSE){
	    for (j = p->waiting_row_num - 1; j >= 0; j--){
	       if (same_cuts_u(p, p->waiting_rows[j], new_rows[i]) !=
		   DIFFERENT_CUTS){
		  free_waiting_row(new_rows+i);
		  new_rows[i] = new_rows[--new_row_num];
		  deleted_cut = TRUE;
		  break;
	       }
	    }
	    if (!deleted_cut) i++;
	 }
      }
      if (new_row_num){
	 add_new_rows_to_waiting_rows(p, new_rows, new_row_num);
	 FREE(new_rows);
      }
      FREE(cuts);
      return(FUNCTION_TERMINATED_NORMALLY);
    default:
      /* Unexpected return value. Do something!! */
      FREE(cuts);
      return(ERROR__USER);
   }
}

/*===========================================================================*/

void print_stat_on_cuts_added_u(lp_prob *p, int added_rows)
{
   int user_res;

#ifdef USE_SYM_APPLICATION
   user_res = user_print_stat_on_cuts_added(p->user, added_rows,
					    p->waiting_rows);
#else
   user_res = USER_DEFAULT;
#endif

   switch(user_res){
    case USER_ERROR:
    case USER_DEFAULT:
      /* print out how many cuts have been added */
      PRINT(p->par.verbosity, 5,
	    ("Number of cuts added to the problem: %i\n", added_rows));
      break;
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      break;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }
}

/*===========================================================================*/

void purge_waiting_rows_u(lp_prob *p)
{
   int user_res, i, j;
   waiting_row **wrows = p->waiting_rows;
   int wrow_num = p->waiting_row_num;
   char *delete_rows;
   int   max_cut_num_per_iter;

   REMALLOC(p->lp_data->tmp.cv, char, p->lp_data->tmp.cv_size, wrow_num,
	    BB_BUNCH);
   delete_rows = p->lp_data->tmp.cv; /* wrow_num */

   memset(delete_rows, 0, wrow_num);

#ifdef USE_SYM_APPLICATION
   user_res = user_purge_waiting_rows(p->user, wrow_num, wrows, delete_rows);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){
    case USER_ERROR: /* purge all */
      free_waiting_rows(wrows, wrow_num);
      p->waiting_row_num = 0;
      break;
    case USER_DEFAULT: /* the only default is to keep enough for one
			  iteration */
      max_cut_num_per_iter = (p->bc_level<1) ? p->par.max_cut_num_per_iter_root
                                             : p->par.max_cut_num_per_iter;
      if (wrow_num - max_cut_num_per_iter > 0){
	 free_waiting_rows(wrows + max_cut_num_per_iter,
			   wrow_num-max_cut_num_per_iter);
	 p->waiting_row_num = max_cut_num_per_iter;
      }
      break;
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      for (i = j = 0; i < wrow_num; i++){
	 if (delete_rows[i]){
	    free_waiting_row(wrows + i);
	 }else{
	    wrows[j++] = wrows[i];
	 }
      }
      p->waiting_row_num = j;
      break;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }
}

/*===========================================================================*/
/*===========================================================================*\
 * This function invokes the user written function user_free_prob_dependent
 * that deallocates the user defined part of the data structure. Returns TRUE
 * if succeeded, FALSE otherwise.
\*===========================================================================*/

void free_prob_dependent_u(lp_prob *p)
{

#ifdef USE_SYM_APPLICATION
   switch (user_free_lp(&p->user)){
    case USER_ERROR:
      /* SYMPHONY ignores error message */
    case USER_SUCCESS:
    case USER_AND_PP:
    case USER_NO_PP:
      /* User function terminated without problems. No post-processing. */
      return;
    default:
      /* Unexpected return value. Do something!! */
      break;
   }
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int analyze_multicriteria_solution(lp_prob *p, int *indices, double *values,
				    int length, double *true_objval,
				    double etol, char branching)
{
  double obj[2] = {0.0, 0.0};
  int i;
  char new_solution = FALSE;
  int continue_with_node = FALSE;
  bool has_artificial = false;

  for (i = 0; i < length; i++){
     if (indices[i] == p->mip->n){
	has_artificial = true;
	continue;
     }
     obj[0] += p->mip->obj1[indices[i]]*values[i];
     obj[1] += p->mip->obj2[indices[i]]*values[i];
  }
  if (has_artificial) length--;

  if (p->has_mc_ub && *true_objval-p->par.mc_rho*(obj[0]+obj[1]) >
      p->mc_ub + etol + MAX(0, MIN(p->par.mc_gamma, p->par.mc_tau))){
     return(FALSE);
  }

  if (p->par.mc_gamma == 1.0){
     if (!p->has_mc_ub || obj[0] < p->obj[0] + etol){
	if (!p->has_mc_ub || (obj[0] < p->obj[0] - etol ||
			      (obj[0] >= p->obj[0] - etol
			       && obj[1] < p->obj[1] - etol))){
	    if (p->par.verbosity >= 1){
	       printf("\nBetter Solution Found:\n");
	       if(p->mip->obj_sense == SYM_MAXIMIZE){
		  printf("First Objective Cost: %.1f\n", -obj[0]);
		  printf("Second Objective Cost: %.1f\n", -obj[1]);
	       }else{
		  printf("First Objective Cost: %.1f\n", obj[0]);
		  printf("Second Objective Cost: %.1f\n", obj[1]);
	       }
	    }
	    p->obj[1] = obj[1];
	    p->obj[0] = obj[0];
	    p->mc_ub = *true_objval-p->par.mc_rho*(obj[0]+obj[1]);
	    p->has_mc_ub = TRUE;
	    new_solution = TRUE;
	 }
	/* Add an optimality cut for the second objective */
	 if (!branching && p->par.mc_add_optimality_cuts){
	    cut_data *new_cut = (cut_data *) calloc(1, sizeof(cut_data));
	    new_cut->coef = NULL;
	    new_cut->rhs = obj[1] - 1 + etol;
	    new_cut->size = 0;
	    new_cut->type = OPTIMALITY_CUT_SECOND;
	    new_cut->name = CUT__DO_NOT_SEND_TO_CP;
	    continue_with_node = cg_add_user_cut(new_cut,
						 &p->cgp->cuts_to_add_num,
						 &p->cgp->cuts_to_add_size,
						 &p->cgp->cuts_to_add);
	    FREE(new_cut);
	 }else{
	    continue_with_node = TRUE;
	 }
     }
  }else if (p->par.mc_tau == 1.0){
     if (!p->has_mc_ub || obj[1] < p->obj[1] + etol){
	if (!p->has_mc_ub || (obj[1] < p->obj[1] - etol ||
			      (obj[1] >= p->obj[1] - etol
			       && obj[0] < p->obj[0] - etol))){
	   if (p->par.verbosity >= 1){
	      printf("\nBetter Solution Found:\n");
	      if(p->mip->obj_sense == SYM_MAXIMIZE){
		 printf("First Objective Cost: %.1f\n", -obj[0]);
		 printf("Second Objective Cost: %.1f\n", -obj[1]);
	      }else{
		 printf("First Objective Cost: %.1f\n", obj[0]);
		 printf("Second Objective Cost: %.1f\n", obj[1]);
	      }
	   }
	   p->obj[1] = obj[1];
	   p->obj[0] = obj[0];
	   p->mc_ub = *true_objval-p->par.mc_rho*(obj[0]+obj[1]);
	   p->has_mc_ub = TRUE;
	   new_solution = TRUE;
	}
	/* Add an optimality cut for the second objective */
	if (!branching && p->par.mc_add_optimality_cuts){
	   cut_data *new_cut = (cut_data *) calloc(1, sizeof(cut_data));
	   new_cut->coef = NULL;
	   new_cut->rhs = obj[0] - 1 + etol;
	   new_cut->size = 0;
	   new_cut->type = OPTIMALITY_CUT_FIRST;
	   new_cut->name = CUT__DO_NOT_SEND_TO_CP;
	   continue_with_node = cg_add_user_cut(new_cut,
						&p->cgp->cuts_to_add_num,
						&p->cgp->cuts_to_add_size,
						&p->cgp->cuts_to_add);
	   FREE(new_cut);
	}else{
	   continue_with_node = TRUE;
	}
     }
  }else{
     if (!p->has_mc_ub ||
	 (p->has_mc_ub && *true_objval-p->par.mc_rho*(obj[0]+obj[1]) <
	  p->mc_ub - MIN(p->par.mc_gamma, p->par.mc_tau) + 100*etol) ||
	 (obj[0] < p->obj[0] - etol &&
	  obj[1] < p->obj[1] + etol + MIN(p->par.mc_gamma, p->par.mc_tau)) ||
	 (obj[1] < p->obj[1] - etol &&
	  obj[0] < p->obj[0] + etol + MIN(p->par.mc_gamma, p->par.mc_tau))){
	if (p->par.verbosity >= 1){
	   printf("\nBetter Solution Found:\n");
	   if(p->mip->obj_sense == SYM_MAXIMIZE){
	      printf("First Objective Cost: %.1f\n", -obj[0]);
	      printf("Second Objective Cost: %.1f\n", -obj[1]);
	   }else{
	      printf("First Objective Cost: %.1f\n", obj[0]);
	      printf("Second Objective Cost: %.1f\n", obj[1]);
	   }
	}
	p->obj[1] = obj[1];
	p->obj[0] = obj[0];
	p->mc_ub = *true_objval-p->par.mc_rho*(obj[0]+obj[1]);
	p->has_mc_ub = TRUE;
	new_solution = TRUE;
     }
     if (!branching && !p->par.mc_find_supported_solutions &&
	 p->par.mc_add_optimality_cuts){
	if (p->par.mc_gamma*(obj[0] - p->utopia[0]) >
	    *true_objval-p->par.mc_rho*(obj[0]+obj[1])-etol){
	   /* Add an optimality cut for the second objective */
	   cut_data *new_cut = (cut_data *) calloc(1, sizeof(cut_data));
	   new_cut->coef = NULL;
	   new_cut->rhs = obj[1] - 1 + etol;
	   new_cut->size = 0;
	   new_cut->type = OPTIMALITY_CUT_SECOND;
	   new_cut->name = CUT__DO_NOT_SEND_TO_CP;
	   continue_with_node = cg_add_user_cut(new_cut,
						&p->cgp->cuts_to_add_num,
						&p->cgp->cuts_to_add_size,
						&p->cgp->cuts_to_add);
	   FREE(new_cut);
	}else if (!p->par.mc_find_supported_solutions &&
		  p->par.mc_add_optimality_cuts){
	   /* Add an optimality cut for the second objective */
	   cut_data *new_cut = (cut_data *) calloc(1, sizeof(cut_data));
	   new_cut->coef = NULL;
	   new_cut->rhs = obj[0] - 1 + etol;
	   new_cut->size = 0;
	   new_cut->type = OPTIMALITY_CUT_FIRST;
	   new_cut->name = CUT__DO_NOT_SEND_TO_CP;
	   continue_with_node = cg_add_user_cut(new_cut,
						&p->cgp->cuts_to_add_num,
						&p->cgp->cuts_to_add_size,
						&p->cgp->cuts_to_add);
	   FREE(new_cut);
	}
     }else if (branching){
	continue_with_node = TRUE;
     }
  }

  if (!new_solution){
     return(continue_with_node);
  }

  p->best_sol.xlevel = p->bc_level;
  p->best_sol.xindex = p->bc_index;
  p->best_sol.xiter_num = p->iter_num;
  p->best_sol.xlength = length;
  p->best_sol.lpetol = etol;
  p->best_sol.objval = *true_objval-p->par.mc_rho*(obj[0]+obj[1]);
  FREE(p->best_sol.xind);
  FREE(p->best_sol.xval);
  p->best_sol.xind = (int *) malloc(length*ISIZE);
  p->best_sol.xval = (double *) malloc(length*DSIZE);
  memcpy((char *)p->best_sol.xind, (char *)indices, length * ISIZE);
  memcpy((char *)p->best_sol.xval, (char *)values, length * DSIZE);

  if(!p->best_sol.has_sol){
     p->best_sol.has_sol = TRUE;
  }

#ifndef COMPILE_IN_TM
  send_feasible_solution_u(p, p->bc_level, p->bc_index, p->iter_num,
			   lpetol, *true_objval-p->par.mc_rho*(obj[0]+obj[1]),
			   length, indices, values);
#endif
  display_lp_solution_u(p, DISP_FEAS_SOLUTION);

  return(continue_with_node);

}
