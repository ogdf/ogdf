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

#include <math.h>
#include <memory.h>
#include <stdlib.h>

#include "sym_lp.h"
#include "sym_qsort.h"
#include "sym_proccomm.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_types.h"
#include "sym_lp_solver.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP functions related to branching.
\*===========================================================================*/

void add_slacks_to_matrix(lp_prob *p, int cand_num, branch_obj **candidates)
{
   LPdata *lp_data = p->lp_data;
   int *index;
   int m = p->lp_data->m;
   int j, k;
   branch_obj *can;
   row_data *newrows;
   waiting_row **wrows;

   for (j=cand_num-1; j >= 0; j--)
      if (candidates[j]->type == CANDIDATE_CUT_NOT_IN_MATRIX)
	 break;

   if (j < 0) /* there is nothing to add */
      return;

   /* We'll create a waiting row for each cut, add them to the matrix,
      then set their status to be free */
   wrows = (waiting_row **) malloc(cand_num * sizeof(waiting_row *));
   /* can't use tmp.p, because that might get resized in add_row_set */
   for (k=0; j >= 0; j--){
      can = candidates[j];
      if (can->type == CANDIDATE_CUT_NOT_IN_MATRIX){
	 wrows[k] = can->row;
	 can->row = NULL;
	 can->position = m + k;
	 can->type = CANDIDATE_CUT_IN_MATRIX;
	 k++;
      }
   }
   add_row_set(p, wrows, k);
   /* To satisfy the size requirements in free_row_set, the following sizes
    * are needed: tmp.c:2*m   tmp.i1:3*m   tmp.d:m */
   FREE(wrows);
   index = lp_data->tmp.i1;
   for (j = 0; j < k; j++)
      index[j] = m + j;
   free_row_set(lp_data, k, index);
   newrows = lp_data->rows + m; /* m is still the old one! */
   for (j = 0; j < k; j++){
      newrows[j].ineff_cnt = (MAXINT) >> 1; /* it is slack... */
      newrows[j].free = TRUE;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Ok. So there were violated cuts, either in waiting_rows or among the
 * slacks (or both). We just have to add those among the slacks to the
 * waiting_rows; add the appropriate number of cuts to the problem and
 * continue the node (this second part is just like the end of receive_cuts().
\*===========================================================================*/

int add_violated_slacks(lp_prob *p, int cand_num, branch_obj **candidates)
{
   LPdata *lp_data = p->lp_data;
   waiting_row **new_rows;
   int i, new_row_num = 0;

   /* If there are any violated (former) slack, unpack them and add them
    * to the set of waiting rows. */
   if (cand_num > 0){
      new_rows = (waiting_row **) lp_data->tmp.p1; /* m (actually, candnum<m */
      for (i=0; i<cand_num; i++){
	 if (candidates[i]->type == VIOLATED_SLACK){
	    new_rows[new_row_num++] = candidates[i]->row;
	    candidates[i]->row = NULL;
	 }
      }
      if (new_row_num > 0)
	 add_new_rows_to_waiting_rows(p, new_rows, new_row_num);
   }

   return( p->waiting_row_num == 0 ? 0 : add_best_waiting_rows(p) );
}

/*===========================================================================*/

int select_branching_object(lp_prob *p, int *cuts, branch_obj **candidate)
{

   LPdata *lp_data = p->lp_data;
   var_desc **vars;
   row_data *rows;
   int m;
#ifndef MAX_CHILDREN_NUM
   int maxnum;
   double *objval, *pobj;
   int *termcode, *pterm, *feasible, *pfeas, *iterd, *piter;
#ifdef COMPILE_FRAC_BRANCHING
   int *frnum, *pfrnum, **frind, **pfrind;
   double **frval, **pfrval;
#endif
#endif
   int i, j, k, l, branch_var, branch_row, min_ind;
   double lb, ub, oldobjval, min_obj;
   cut_data *cut;
   branch_obj *can, *best_can = NULL;
#ifdef COMPILE_FRAC_BRANCHING
   int *xind;
   double *xval;
#endif

   /* These are the return values from select_candidates_u() */
   int cand_num = 0, new_vars = 0, *indices;
   double *values;
   branch_obj **candidates = NULL;
#ifdef STATISTICS
   int itlim = 0, cnum = 0;
#endif
   int should_use_hot_starts = FALSE;
   double total_time = 0;
   double st_time = 0;
   int total_iters, should_continue;//, max_iter_num;
   int should_use_rel_br = p->par.should_use_rel_br;
   double high, low, down_obj, up_obj, best_var_score;
   int *br_rel_down = p->br_rel_down;
   int *br_rel_up = p->br_rel_up;
   double *pcost_down = p->pcost_down;
   double *pcost_up = p->pcost_up;
   int rel_threshold = p->par.rel_br_threshold;
   int best_var;
   // double lp_time, scaled_by;
   int max_presolve_iter;
   const int bc_level = p->bc_level;
   int strong_br_min_level = p->par.strong_br_min_level;

   /*------------------------------------------------------------------------*\
    * First we call select_candidates_u() to select candidates. It can
    * -- return with DO_BRANCH and a bunch of candidates, or
    * -- return with DO_NOT_BRANCH along with a bunch of violated cuts
    *    in the matrix and/or among the slack_cuts, or
    * -- return with DO_NOT_BRANCH__FATHOMED, i.e., the node can be fathomed.
   \*------------------------------------------------------------------------*/

   j = select_candidates_u(p, cuts, &new_vars, &cand_num, &candidates);

   switch (j){
    case DO_NOT_BRANCH__FATHOMED:
      *candidate = NULL;
      return(DO_NOT_BRANCH__FATHOMED);
    case DO_NOT_BRANCH__FEAS_SOL:
      *candidate = NULL;
      return(DO_NOT_BRANCH__FEAS_SOL);
    case DO_NOT_BRANCH:
      if (cand_num)
	 *cuts += add_violated_slacks(p, cand_num, candidates);
#ifdef DO_TESTS
      if (*cuts == 0 && new_vars == 0){
	 printf("Error! Told not to branch, but there are no new cuts or ");
	 printf("variables!\n");
	 *candidate = NULL;
	 return(ERROR__NO_BRANCHING_CANDIDATE);
      }
#endif
      /* Free the candidates */
      if (candidates){
	 for (i=0; i<cand_num; i++){
	    free_candidate(candidates + i);
	 }
	 FREE(candidates);
      }
      *candidate = NULL;
      return(DO_NOT_BRANCH);

    case DO_BRANCH:
      break;

   case ERROR__NO_BRANCHING_CANDIDATE:
      *candidate = NULL;
      return(ERROR__NO_BRANCHING_CANDIDATE);
   }

   /* OK, now we have to branch. */

   /* First of all, send everything to the cutpool that hasn't been sent
      before and send the current node description to the TM. */
   p->comp_times.strong_branching += used_time(&p->tt);
#pragma omp critical(cut_pool)
   send_cuts_to_pool(p, -1);
   send_node_desc(p, NODE_BRANCHED_ON);
   p->comp_times.communication += used_time(&p->tt);

   /* Add all the branching cuts */
   if (p->par.branch_on_cuts)
      add_slacks_to_matrix(p, cand_num, candidates);
   m = lp_data->m;
   rows = lp_data->rows;


#ifndef MAX_CHILDREN_NUM
   /* The part below is not needed when we have MAX_CHILDREN_NUM specified */
   /* Count how many objval/termcode/feasible entry we might need
      and allocate space for it */
   for (maxnum = candidates[0]->child_num, j=0, i=1; i<cand_num; i++){
      if (maxnum < candidates[i]->child_num)
	 maxnum = candidates[i]->child_num;
   }

   objval   = (double *) malloc(maxnum * DSIZE);
   termcode = (int *) malloc(maxnum * ISIZE);
   feasible = (int *) malloc(maxnum * ISIZE);
   iterd    = (int *) malloc(maxnum * ISIZE);
#ifdef COMPILE_FRAC_BRANCHING
   frval = (double **) malloc(maxnum * sizeof(double *));
   pfrval = (double **) malloc(maxnum * sizeof(double *));
   frind = (int **) malloc(maxnum * sizeof(int *));
   pfrind = (int **) malloc(maxnum * sizeof(int *));
   frnum = (int *) malloc(maxnum * ISIZE);
   pfrnum = (int *) malloc(maxnum * ISIZE);
#endif
   pobj  = (double *) malloc(maxnum * DSIZE);
   pterm = (int *) malloc(maxnum * ISIZE);
   pfeas = (int *) malloc(maxnum * ISIZE);
   piter = (int *) malloc(maxnum * ISIZE);
#endif

   /* Look at the candidates one-by-one and presolve them. */
   vars = lp_data->vars;
   oldobjval   = lp_data->objval;
   st_time     = used_time(&total_time);
   total_iters = 0;

   if (should_use_rel_br==TRUE) {

      const double lpetol100 = lp_data->lpetol*100;
      double *x = (double *)malloc (lp_data->n*DSIZE);
      double *bnd_val = (double *)malloc (2*lp_data->n*DSIZE);
      int *bnd_ind = (int *)malloc (2*lp_data->n*ISIZE);
      char *bnd_sense = (char *)malloc (2*lp_data->n*CSIZE);
      int num_bnd_changes = 0;
      double xval, floorx, ceilx, var_score;
      int full_solves = 0, down_is_est, up_is_est,
          max_solves_since_impr = p->par.rel_br_cand_threshold,
          stop_solving = FALSE, both_children_inf = FALSE, rel_up,
          rel_down, solves_since_impr = 0;
      int max_solves = p->par.rel_br_max_solves;
      double alpha = p->par.strong_branching_high_low_weight;
      double one_m_alpha = 1.0 - alpha;
      int check_first = FALSE;
      int check_level = 0;
      int num_up_iters = 0, num_down_iters = 0;
      int up_status = -1, down_status = -1;

      best_var = -1;
      best_var_score = -SYM_INFINITY;
      memcpy(x, lp_data->x, lp_data->n*DSIZE);

#ifdef COMPILE_IN_LP

      int hotstart_reg_cand = 0;
      int hotstart_ch_cand = 0;
      int check_off = TRUE;

      if(p->par.rel_br_override_default && p->mip->mip_inf){


	 int weighted_iter =
	    p->lp_stat.lp_total_iter_num/(p->lp_stat.lp_calls -
					  p->lp_stat.str_br_lp_calls -
					  p->lp_stat.fp_lp_calls + 1);
	 if(p->mip->nz > 5e3){
	    weighted_iter = (int)
	       ((1.0*weighted_iter * p->mip->nz) / 5e3);
	 }

	 if(p->mip->nz > 5e4){
	    rel_threshold = MAX(2, (int)(1.0 * rel_threshold * 5e4/p->mip->nz));
	 }

	 if(p->bc_level < 1){
	    if(p->iter_num > 2 && weighted_iter <= 1000){
	       if(p->mip->mip_inf){
		  if(p->mip->mip_inf->prob_type == BINARY_TYPE){
		     strong_br_min_level =
			MIN(p->par.strong_br_min_level,
			    (int)((p->mip->mip_inf->binary_var_num)/10.0) + 1);
		  }
		  if(p->mip->mip_inf->prob_type == BIN_CONT_TYPE){
		     if(p->mip->mip_inf->bin_var_ratio < 0.1){
			strong_br_min_level =
			   MIN(MAX(p->par.strong_br_min_level,
				   (int)((p->mip->mip_inf->binary_var_num)/10.0)
				   + 1), 10);
			if(p->mip->nz < 5e4){
			   strong_br_min_level = MAX(p->par.strong_br_min_level, 8);
			}
		     }
		  }
	       }
	    }
	 }

	 if(weighted_iter * p->bc_index < 5e7){
	    check_off = FALSE;
	 }


	 if(p->mip->mip_inf && p->mip->mip_inf->bin_cont_row_num > 0 &&
	    (p->mip->mip_inf->bin_cont_row_num >= p->mip->m ||
	     (p->mip->mip_inf->bin_var_ratio < 0.2) ||
	     p->mip->n - p->mip->mip_inf->cont_var_num <= 100 ||
	     (p->mip->mip_inf->sos_bin_row_ratio < 0.05 &&
	      p->mip->mip_inf->bin_var_ratio < 0.6) ||
	     (p->mip->mip_inf->max_row_ratio < 0.01 &&
	      p->mip->mip_inf->prob_type != BIN_CONT_TYPE))){
	    /* -either we have all continuos rows
	       -less number of bin vars
	       -small number of int vars
	       -large bin but less sos rows
	       -small max_row_size - if rest are all binary, skip to latter
	    */

	    if(p->mip->mip_inf->bin_cont_row_num >= p->mip->m){
	       max_solves = MIN(max_solves, 2*cand_num);

	    }else if(p->mip->mip_inf->bin_var_ratio < 0.2){
	       max_solves = MIN(max_solves, 2*cand_num);
	       if(p->mip->mip_inf->bin_var_ratio > 0.05){
		  strong_br_min_level = (int)strong_br_min_level/2;
	       }
	    }else{// if(p->mip->mip_inf->max_row_ratio < 0.01){
	       max_solves = MIN(2*max_solves, 2*cand_num);
	       if(p->mip->mip_inf->sos_bin_row_ratio > 0.05){
		  //  max_solves = MIN(2*max_solves, 2*cand_num);
		  strong_br_min_level = (int)(2.0*strong_br_min_level);
		  rel_threshold = 2*rel_threshold;
	       }//else if(p->mip->mip_inf->bin_var_ratio > 0.2){
		//  max_solves /= 2;
	       // }
	    }///else{
	     //  max_solves = MIN(2*max_solves, 2*cand_num);
	    //if(p->mip->n - p->mip->mip_inf->cont_var_num <= 100 &&
	    //  p->mip->nz < 1e4){
	    //printf("here");
	    //max_solves = 2*cand_num;
	    //rel_threshold = MAX(2, rel_threshold/2);
	    //printf("rel_the %i\n",rel_threshold);
	    // }
	    // }
	 }else{
	    double imp_avg = 0.0;
	    int backtrack = 0;

	    bc_node *node = p->tm->active_nodes[p->proc_index];
	    node = p->tm->active_nodes[p->proc_index];
	    /*
	    if(p->bc_level > 0){

		  //if(node->start_objval > node->parent->end_objval + lp_data->lpetol){

	       if(weighted_iter * p->bc_index < 5e7){
		  check_off = FALSE;
	       }

	       if(p->mip->mip_inf && (p->mip->mip_inf->prob_type == BIN_CONT_TYPE ||
				      p->mip->mip_inf->prob_type == BINARY_TYPE)){
		  if(!check_off) check_first = TRUE;
	       }
	    }
	    */
	    if(p->bc_level >= 1){
	       while(node->parent){
		  if(node->start_objval > node->parent->end_objval){
		     imp_avg +=
			fabs(node->start_objval/(node->parent->end_objval +0.0001) - 1.0);
		  }
		  node = node->parent;
		  if(backtrack++ > p->par.rel_br_chain_backtrack) break;
	       }
	    }

	    if(backtrack > 0){
	       imp_avg /= backtrack;
	    }

	    if(imp_avg > p->par.rel_br_min_imp &&
	       imp_avg < p->par.rel_br_max_imp){
	       if(bc_level <= strong_br_min_level ){
		  max_solves = MIN(3*max_solves, 2*cand_num);
	       }else{
		  max_solves = MIN(2*max_solves, 2*cand_num);
	       }
	    }else{

	       int c_cnt = 0;
	       double d_avg = 0.0;

	       for (i=0; i<cand_num; i++) {
		  branch_var = p->br_rel_cand_list[i];
		  xval = x[branch_var];
		  floorx = floor(xval);
		  ceilx = ceil(xval);
		  rel_down = br_rel_down[branch_var];
		  rel_up = br_rel_up[branch_var];
		  if(xval - floorx > 0.5){
		     d_avg += ceilx - xval;
		  }else{
		     d_avg += xval - floorx;
		  }
	       }

	       d_avg /= cand_num;

	       for (i=0; i<cand_num; i++) {
		  branch_var = p->br_rel_cand_list[i];
		  xval = x[branch_var];
		  if(xval - floor(xval) > d_avg && ceil(xval) - xval > d_avg) c_cnt++;
		  else break;
	       }

	       if(bc_level < 1){
		  max_solves = cand_num < p->par.rel_br_override_max_solves ?
		     MIN(p->par.rel_br_override_max_solves/2, cand_num) :
		     p->par.rel_br_override_max_solves;
	       }else if(bc_level < 4){
		  max_solves = MIN((int)(0.75*c_cnt), (int)(0.3 * cand_num) + 1);
	       }else if(bc_level < 20){
		  max_solves = MIN(c_cnt/2, (int)(0.25 * cand_num) + 1);
	       }else if(bc_level < 40){
		  max_solves = MIN(c_cnt/3, (int)(0.20 * cand_num) + 1);
	       }else{
		  max_solves = MAX(c_cnt/4, (int)(0.15 * cand_num) + 1);
	       }
	    }

	    max_solves_since_impr  = 5;

	    //printf("level - set to : %i %i\n", p->bc_level, max_solves);
	    //printf("c_cnt - cand num - max_solves : %i %i %i\n\n",
	    //c_cnt,cand_num, max_solves);
	    //if(check_off){
	    int int_num = p->mip->n - p->mip->mip_inf->cont_var_num;
	    int max_level = (p->mip->mip_inf == 0) ? 500 :
	       (int_num)/2;
	    max_level = MIN(500, MAX(100, max_level));
	    // if(cand_num > 0.25*(p->mip->n -
	    //		   p->mip->mip_inf->cont_var_num))
	    //  max_level = MIN(50, max_level);
	    //if(cand_num > 500) max_level = MIN(100, max_level);
	    //max_level = 50;
	    if(cand_num > 100 && int_num > 500){
	       max_level = MIN(100, max_level);
	       if((p->mip->mip_inf->prob_type == BINARY_TYPE ||
		   p->mip->mip_inf->prob_type == BIN_CONT_TYPE) &&
		  cand_num > 0.05*int_num){
		  max_level /= 2;
	       }
	    }

	    if(bc_level > max_level){
	       rel_threshold = max_solves = 0;
	       strong_br_min_level = 1;
	       //cand_num = 1;
	    }
	 }
	 //}
	 /*
	 for (i=0; i<cand_num; i++) {
	    branch_var = p->br_rel_cand_list[i];
	    printf("x - xval %s %f\n", p->mip->colname[branch_var], lp_data->x[branch_var]);
	 }
	 */

	 max_solves = MIN(p->par.rel_br_override_max_solves, max_solves);

	 double rel_limit = 0.05;
	 if((p->mip->mip_inf && ((p->mip->mip_inf->mat_density < rel_limit &&
				 p->mip->mip_inf->int_var_ratio > rel_limit &&
				 (p->mip->mip_inf->max_col_ratio > rel_limit ||
				  p->mip->mip_inf->max_row_ratio > rel_limit))||
	     (p->mip->nz > 1e5 && p->mip->mip_inf->mat_density > rel_limit/50) ||
	     (p->mip->mip_inf->max_row_ratio < rel_limit/5 &&
	      p->mip->mip_inf->prob_type != BIN_CONT_TYPE)))){
#ifdef __OSI_CLP__
	    lp_data->si->setupForRepeatedUse(2,0);
#endif
	 }

	 if(p->mip->mip_inf && !check_off &&
	    (p->mip->mip_inf->prob_type == BINARY_TYPE ||
	     p->mip->mip_inf->prob_type == BIN_CONT_TYPE) &&
	    (p->mip->n - p->mip->mip_inf->cont_var_num < 100 ||
	     (p->mip->mip_inf->int_var_ratio > rel_limit &&
	      p->mip->mip_inf->row_density/(p->mip->n + 1) > rel_limit/5))){//
	    //	      && p->mip->mip_inf->int_var_ratio < 1-rel_limit))){
	    check_first = TRUE;
	 }

	 //  printf("c_cnt - cand num - max_solves : %i %i %i\n\n",
	 //c_cnt,cand_num, max_solves);
	 //printf("cand num - max_solves : %i %i \n\n", cand_num, max_solves);
	 //    printf("str_level - thresh: %i %i \n",strong_br_min_level, rel_threshold);

	 if(p->mip->mip_inf->binary_sos_row_num > 0){
	    double bin_den = (1.0*p->mip->mip_inf->binary_sos_row_num)/
	       (p->mip->m + 1);
	    //1.0*p->mip->mip_inf->binary_sos_row_num/(p->mip->m+1));
	    if( bin_den > rel_limit && ((bin_den < 10*rel_limit &&
					 p->mip->mip_inf->prob_type != BINARY_TYPE &&
					 p->mip->mip_inf->bin_var_ratio > 10*rel_limit) ||
					(bin_den < 2*rel_limit &&
					 p->mip->mip_inf->prob_type == BINARY_TYPE))){
	       /* give priority to vars appear in sos rows */
	       int *sos_ind = (int *)(malloc)(ISIZE*cand_num);
	       int *sos_tot_var = (int *)(malloc)(ISIZE*cand_num);
	       int sos_cnt = 0;
	       for (i=0; i<cand_num; i++) {
		  branch_var = p->br_rel_cand_list[i];
		  //printf("%i %i\n", branch_var, p->mip->mip_inf->cols[branch_var].sos_num);
		  // if(p->mip->mip_inf->cols[branch_var].sos_num > 0.1*p->mip->n){
		  if(p->mip->mip_inf->cols[branch_var].sos_num >= (1.0*p->mip->nz)/(p->mip->m + 1)){
		     sos_tot_var[sos_cnt] = p->mip->n -
			p->mip->mip_inf->cols[branch_var].sos_num;
		     sos_ind[sos_cnt] = i;
		     sos_cnt++;
		  }
	       }
	       //printf("sos_cnt %i\n", sos_cnt);
	       if(sos_cnt > 0){
		  qsort_ii(sos_tot_var, sos_ind, sos_cnt);
		  int *sos_chosen = (int *)(calloc)(ISIZE,cand_num);
		  int *new_ord = (int *)(malloc)(ISIZE*cand_num);

		  for (i=0; i<MIN(max_solves/2 + 1, sos_cnt); i++) {
		     new_ord[i] = p->br_rel_cand_list[sos_ind[i]];
		     sos_chosen[sos_ind[i]] = TRUE;
		  }

		  if(i < cand_num){
		     int rest_cnt = 0;
		  for(j = 0; j < cand_num; j++){
		     if(sos_chosen[j]) continue;
		     else new_ord[i+rest_cnt++] = p->br_rel_cand_list[j];
		  }
		  }

		  memcpy(p->br_rel_cand_list, new_ord, ISIZE*cand_num);

		  FREE(sos_chosen);
		  FREE(new_ord);
	       }

	       FREE(sos_ind);
	       FREE(sos_tot_var);
	    }
	 }
      }

      for (i=0; i<cand_num; i++) {
	 branch_var = p->br_rel_cand_list[i];
	 rel_down = br_rel_down[branch_var];
	 rel_up = br_rel_up[branch_var];

	 if (rel_down <= rel_threshold ||
	     bc_level <= strong_br_min_level){
	    hotstart_reg_cand++;
	 }
	 if(rel_up <= rel_threshold ||
	    bc_level <= strong_br_min_level){
	    hotstart_reg_cand++;
	 }

	 if(check_first && i < check_level + 1){
	    hotstart_ch_cand += 2;
	 }
      }

      if(p->par.use_hot_starts && !p->par.branch_on_cuts &&
	 (hotstart_reg_cand > 0 || hotstart_ch_cand > 0)){
	 should_use_hot_starts = TRUE;
      }else{
	 should_use_hot_starts = FALSE;
      }
#else
      if(p->par.use_hot_starts && !p->par.branch_on_cuts){
	 should_use_hot_starts = TRUE;
      }else{
	 should_use_hot_starts = FALSE;
      }
#endif

      if (should_use_hot_starts) {
	 mark_hotstart(lp_data);
      }

      /* Set the iteration limit */

      if (p->par.max_presolve_iter > 0) {
	 max_presolve_iter = p->par.max_presolve_iter - bc_level;

#ifdef COMPILE_IN_LP
	 if(p->mip->nz > 5e4){
	    max_presolve_iter = (int)(1.0 * max_presolve_iter * 5e4/p->mip->nz);
	 }
#endif
	 if (max_presolve_iter < 5) {
	    max_presolve_iter = 5;
	 }

	 if(should_use_hot_starts){
	   set_itlim_hotstart(lp_data, max_presolve_iter);
	 }
	 set_itlim(lp_data, max_presolve_iter);
      }

      for (i=0; i<cand_num; i++) {
         branch_var = p->br_rel_cand_list[i];
         lb = vars[branch_var]->new_lb;
         ub = vars[branch_var]->new_ub;
         xval = x[branch_var];
	 floorx = floor(xval);
         ceilx = ceil(xval);
         rel_down = br_rel_down[branch_var];
         rel_up = br_rel_up[branch_var];

	 if ((rel_down > rel_threshold &&
	      bc_level > strong_br_min_level) &&
	     (i > check_level || (i < check_level + 1 && !check_first ))){
	    down_obj = oldobjval + pcost_down[branch_var] * (xval - floorx);
            down_is_est = TRUE;
	    p->lp_stat.rel_br_pc_down_num++;
         } else {
            if (stop_solving == TRUE){
               continue;
            }
            if (strong_branch(p, branch_var, lb, ub, lb, floorx, &down_obj,
			      should_use_hot_starts, &down_status, &num_down_iters)) {
               // lp was abandoned
               continue;
            }
            down_is_est = FALSE;
	    if(p->bc_level < p->br_rel_down_min_level[branch_var]){
	       p->br_rel_down_min_level[branch_var] =  p->bc_level;
	    }
            if (down_status == LP_D_INFEASIBLE || down_status == LP_D_OBJLIM ||
                down_status == LP_D_UNBOUNDED) {
               // update bounds
               bnd_val[num_bnd_changes] = ceilx;
               bnd_sense[num_bnd_changes] = 'G';
               bnd_ind[num_bnd_changes] = branch_var;
               num_bnd_changes++;
               change_lbub(lp_data, branch_var, ceilx, ub);
               vars[branch_var]->new_lb = ceilx;
               vars[branch_var]->lb = ceilx;
               lb = ceilx;
            } else {
              // update pcost
	       pcost_down[branch_var] = (pcost_down[branch_var]*
		  rel_down + (down_obj - oldobjval)/(xval-floorx))/
		  (rel_down + 1);
               br_rel_down[branch_var]++;
	       p->lp_stat.rel_br_down_update++;
            }
            full_solves++;
            solves_since_impr++;
	    p->lp_stat.rel_br_full_solve_num++;
         }
	 if ((br_rel_up[branch_var] > rel_threshold &&
	      bc_level > strong_br_min_level) &&
	     (i > check_level || (i < check_level + 1 && !check_first ))){
	    up_obj   = oldobjval + pcost_up[branch_var] * (ceilx - xval);
            up_is_est = TRUE;
	    p->lp_stat.rel_br_pc_up_num++;
         } else {
            if (stop_solving == TRUE){
               continue;
            }

            if (strong_branch(p, branch_var, lb, ub, ceilx, ub, &up_obj,
			      should_use_hot_starts, &up_status, &num_up_iters)) {
               // lp was abandoned
               continue;
            }
	    if(p->bc_level < p->br_rel_up_min_level[branch_var]){
	       p->br_rel_up_min_level[branch_var] =  p->bc_level;
	    }
            up_is_est = FALSE;
            if (up_status == LP_D_INFEASIBLE || up_status == LP_D_OBJLIM ||
                up_status == LP_D_UNBOUNDED) {
               // update bounds
               bnd_val[num_bnd_changes] = floorx;
               bnd_sense[num_bnd_changes] = 'L';
               bnd_ind[num_bnd_changes] = branch_var;
               num_bnd_changes++;
               change_lbub(lp_data, branch_var, lb, floorx);
               vars[branch_var]->new_ub = floorx;
               vars[branch_var]->ub = floorx;
               ub = floorx;
            } else {
              // update pcost
               pcost_up[branch_var] = (pcost_up[branch_var]*
	         rel_up + (up_obj - oldobjval)/(ceilx-xval))/
	         (rel_up + 1);
               br_rel_up[branch_var]++;
	       p->lp_stat.rel_br_up_update++;
            }
            full_solves++;
            solves_since_impr++;
	    p->lp_stat.rel_br_full_solve_num++;
         }

         if (down_obj > SYM_INFINITY/10 && up_obj > SYM_INFINITY/10) {
            both_children_inf = TRUE;
            best_can = candidates[0];
            break;
         }

         if (down_obj < up_obj) {
            low = down_obj;
            high = up_obj;
         } else {
            low = up_obj;
            high = down_obj;
         }
	  var_score = alpha * low + one_m_alpha * high;
         if (var_score > best_var_score + lpetol100 || best_can == NULL) {
	    if ( var_score > best_var_score + 0.1 &&(down_is_est != TRUE ||
							   up_is_est != TRUE)) {
               solves_since_impr = 0;
	       if(best_can!= NULL){
		  p->lp_stat.rel_br_impr_num++;
	       }
	    }

            best_var_score = var_score;
            best_var = branch_var;
            best_can = candidates[0];
            best_can->position = branch_var;
            best_can->solutions = NULL;
            best_can->sol_inds = NULL;
            best_can->sol_sizes = NULL;
            best_can->sense[0] = 'L';
            best_can->sense[1] = 'G';
            if (down_is_est==TRUE) {
              best_can->objval[0] = oldobjval;
              best_can->iterd[0] = 0;
              best_can->termcode[0] = LP_D_ITLIM;
            } else {
              best_can->objval[0] = down_obj;
              best_can->iterd[0] = num_down_iters;
              best_can->termcode[0] = down_status;
              // added by asm4 because  hot starts dont generate a reliable
              // bound.
              //if (should_use_hot_starts && down_status==LP_D_ITLIM) {
              //  down_is_est = TRUE;
              //  best_can->objval[0] = oldobjval;
              //}
            }
            if (up_is_est==TRUE) {
              best_can->objval[1] = oldobjval;
              best_can->iterd[1] = 0;
              best_can->termcode[1] = LP_D_ITLIM;
            } else {
              best_can->objval[1] = up_obj;
              best_can->iterd[1] = num_up_iters;
              best_can->termcode[1] = up_status;
              // added by asm4 because  hot starts dont generate a reliable
              // bound.
              //if (should_use_hot_starts && up_status==LP_D_ITLIM) {
              //  up_is_est = TRUE;
              //  best_can->objval[1] = oldobjval;
              //}
            }
            best_can->is_est[0] = down_is_est;
            best_can->is_est[1] = up_is_est;
            best_can->rhs[0] = floorx;
            best_can->rhs[1] = ceilx;
            best_can->value = xval;
         }
         if (solves_since_impr > max_solves_since_impr ||
	     full_solves >= max_solves) {
            //printf("breaking because of no gain at iter %d\n", i);
            stop_solving = TRUE;
         }
      }
      //printf("reliability branching: selected var %d with score %f\n", best_var, best_var_score);
#ifdef COMPILE_IN_LP
      if (num_bnd_changes > 0) {
         str_br_bound_changes(p, num_bnd_changes, bnd_val, bnd_ind, bnd_sense);
      }
#endif

      cand_num = 1;
      FREE(x);
      FREE(bnd_val);
      FREE(bnd_ind);
      FREE(bnd_sense);
      if (both_children_inf == TRUE) {
         FREE(best_can);
         FREE(candidates);
         *candidate = NULL;
         p->lp_stat.str_br_nodes_pruned++;
         if (should_use_hot_starts){
            unmark_hotstart(lp_data);
            set_itlim_hotstart(lp_data, -1);
         }
         set_itlim(lp_data, -1); //both limits should be set for hotstarts
         return (DO_NOT_BRANCH__FATHOMED);
      }
   } else {
      /* do the default symphony branching */

      /*
       * see if hot-starts should be used. in theory if strong branching is used
       * and only variable bounds are changed then, hot-starts should be faster
       */

      if (p->par.use_hot_starts && !p->par.branch_on_cuts) {
	 should_use_hot_starts = TRUE;
      } else {
	 should_use_hot_starts = FALSE;
      }

      /* Set the iteration limit */
      if (should_use_hot_starts) {
	 mark_hotstart(lp_data);
      }

      if (p->par.max_presolve_iter > 0) {
	 max_presolve_iter = p->par.max_presolve_iter - bc_level;

	 if (max_presolve_iter < 5) {
	    max_presolve_iter = 5;
	 }
	 if(should_use_hot_starts){
	   set_itlim_hotstart(lp_data, max_presolve_iter);
	 }
	 set_itlim(lp_data, max_presolve_iter);
      }

      for (i=0; i<cand_num; i++){
         can = candidates[i];

#ifndef MAX_CHILDREN_NUM
         can->objval = pobj;
         can->termcode = pterm;
         can->feasible = pfeas;
         can->iterd = piter;
         if (p->tm->par.keep_description_of_pruned == KEEP_IN_MEMORY){
            can->solutions = (double **) calloc(maxnum, sizeof(double *));
            can->sol_inds = (int **) calloc(maxnum, sizeof(int *));
            can->sol_size = (int *) calloc(maxnum, ISIZE);
         }

#ifdef SENSITIVITY_ANALYSIS
         if (p->tm->par.sensitivity_analysis){
            can->duals = (double **) calloc(maxnum, sizeof(double *));
         }else{
            can->duals = NULL;
         }
#endif
#ifdef COMPILE_FRAC_BRANCHING
         can->frac_num = pfrnum;
         can->frac_ind = pfrind;
         can->frac_val = pfrval;
#endif

#else
         if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
            can->solutions = (double **) calloc (MAX_CHILDREN_NUM,
                  sizeof(double *));
            can->sol_inds = (int **) calloc(MAX_CHILDREN_NUM,
                  sizeof(int *));
            can->sol_sizes = (int *) calloc(MAX_CHILDREN_NUM, ISIZE);
         }
#ifdef SENSITIVITY_ANALYSIS
         if (p->par.sensitivity_analysis){
            can->duals = (double **) calloc (MAX_CHILDREN_NUM, sizeof(double *));
         }else{
            can->duals = NULL;
         }
#endif
#endif

#ifdef STATISTICS
         cnum += can->child_num;
#endif

         /* Now depending on the type, adjust ub/lb or rhs/range/sense */
         switch (can->type){
          case CANDIDATE_VARIABLE:
            branch_var = can->position;
#if 0
            if (lp_data->status[branch_var] & PERM_FIXED_TO_LB ||
                  lp_data->status[branch_var] & PERM_FIXED_TO_UB){
               if (vars[branch_var]->lb == vars[branch_var]->ub){
                  printf("Warning -- branching candidate is already fixed. \n");
                  printf("SYMPHONY has encountered numerical difficulties \n");
                  printf("With the LP solver. Exiting...\n\n");
               }
               /* } to unconfuse vi*/
#endif
            lb = vars[branch_var]->new_lb;
            ub = vars[branch_var]->new_ub;
            for (j = 0; j < can->child_num; j++){
               switch (can->sense[j]){
                case 'E':
                  change_lbub(lp_data, branch_var, can->rhs[j], can->rhs[j]);
                  break;
                case 'R':
                  change_lbub(lp_data, branch_var, can->rhs[j],
                        can->rhs[j] + can->range[j]);
                  break;
                case 'L':
                  change_lbub(lp_data, branch_var, lb, can->rhs[j]);
                  break;
                case 'G':
                  change_lbub(lp_data, branch_var, can->rhs[j], ub);
                  break;
               }
               check_ub(p);
               /* The original basis is in lp_data->lpbas */
               if (should_use_hot_starts) {
                  can->termcode[j] = solve_hotstart(lp_data, can->iterd+j);
                  total_iters+=*(can->iterd+j);
               } else {
                  can->termcode[j] = dual_simplex(lp_data, can->iterd+j);
                  total_iters+=*(can->iterd+j);
               }
               p->lp_stat.lp_calls++;
               p->lp_stat.str_br_lp_calls++;
	       p->lp_stat.str_br_total_iter_num += *(can->iterd+j);
               can->objval[j] = lp_data->objval;
               get_x(lp_data);

#ifdef SENSITIVITY_ANALYSIS
               if (p->par.sensitivity_analysis){
                  get_dj_pi(lp_data);
                  can->duals[j] = (double *) malloc (DSIZE*p->base.cutnum);
                  memcpy(can->duals[j], lp_data->dualsol, DSIZE*p->base.cutnum);
               }
#endif

               if (can->termcode[j] == LP_OPTIMAL){
                  /* is_feasible_u() fills up lp_data->x, too!! */
                  switch (is_feasible_u(p, TRUE, FALSE)){

                     /*NOTE: This is confusing but not all that citical...*/
                     /*The "feasible" field is only filled out for the
                       purposes of display (in vbctool) to keep track of
                       where in the tree the feasible solutions were
                       found. Since this may not be the actual candidate
                       branched on, we need to pass this info on to whatever
                       candidate does get branched on so the that the fact that
                       a feasible solution was found in presolve can be recorded*/

                   case IP_FEASIBLE:
                     can->termcode[j] = LP_OPT_FEASIBLE;
                     can->feasible[j] = TRUE;
                     if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
                        can->solutions[j] = (double *) malloc (DSIZE*lp_data->n);
                        memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
                     }
                     break;

                   case IP_FEASIBLE_BUT_CONTINUE:
                     can->termcode[j] = LP_OPT_FEASIBLE_BUT_CONTINUE;
                     can->feasible[j] = TRUE;
                     if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
                        can->solutions[j] = (double *) malloc (DSIZE*lp_data->n);
                        memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
                     }
                     break;

                   default:
                     break;
                  }
               } else if (can->termcode[j] == LP_D_OBJLIM ||
                     can->termcode[j] == LP_D_UNBOUNDED ||
                     can->termcode[j] == LP_D_INFEASIBLE){
                  //p->bound_changes_in_iter++;
                  switch (can->sense[j]){
                   case 'L':
                     /* decreasing the ub made the problem inf, so change lb */
                     //lb = can->rhs[j] + 1;
                     //vars[can->position]->new_lb = lb;
                     break;
                   case 'G':
                     //ub = can->rhs[j] - 1;
                     //vars[can->position]->new_ub = ub;
                     break;
                   case 'E':
                     /* problem becomes infeasible */
                     /* dont know what to do */
                     break;
                  }
               }
#ifdef COMPILE_FRAC_BRANCHING
               else{
                  if (can->termcode[j] != LP_ABANDONED){
                     get_x(lp_data);
                  }
               }
               if (can->termcode[j] != LP_ABANDONED){
                  xind = lp_data->tmp.i1; /* n */
                  xval = lp_data->tmp.d; /* n */
                  can->frac_num[j] = collect_fractions(p, lp_data->x, xind, xval);
                  if (can->frac_num[j] > 0){
                     can->frac_ind[j] = (int *) malloc(can->frac_num[j] * ISIZE);
                     can->frac_val[j] = (double *) malloc(can->frac_num[j]*DSIZE);
                     memcpy(can->frac_ind[j], xind, can->frac_num[j] * ISIZE);
                     memcpy(can->frac_val[j], xval, can->frac_num[j] * DSIZE);
                  }
               }else{
                  can->frac_num[j] = 0;
               }
#endif
#ifdef STATISTICS
               if (can->termcode[j] == LP_D_ITLIM)
                  itlim++;
#endif
            }
            change_lbub(lp_data, branch_var, lb, ub);
            break;

          case CANDIDATE_CUT_IN_MATRIX:
            branch_row = can->position;
            for (j = 0; j < can->child_num; j++){
               change_row(lp_data, branch_row,
                     can->sense[j], can->rhs[j], can->range[j]);
               check_ub(p);
               /* The original basis is in lp_data->lpbas */
               can->termcode[j] = dual_simplex(lp_data, can->iterd+j);
               p->lp_stat.lp_calls++;
               p->lp_stat.str_br_lp_calls++;
	       p->lp_stat.str_br_total_iter_num += *(can->iterd+j);
               can->objval[j] = lp_data->objval;


               get_x(lp_data);

#ifdef SENSITIVITY_ANALYSIS
               if (p->par.sensitivity_analysis){
                  get_dj_pi(lp_data);
                  can->duals[j] = (double *) malloc (DSIZE*p->base.cutnum);
                  memcpy(can->duals[j], lp_data->dualsol, DSIZE*p->base.cutnum);
               }
#endif

               if (can->termcode[j] == LP_OPTIMAL){
                  /* is_feasible_u() fills up lp_data->x, too!! */
                  switch (is_feasible_u(p, TRUE, FALSE)){

                     /*NOTE: This is confusing but not all that citical...*/
                     /*The "feasible" field is only filled out for the
                       purposes of display (in vbctool) to keep track of
                       where in the tree the feasible solutions were
                       found. Since this may not be the actual candidate
                       branched on, we need to pass this info on to whatever
                       candidate does get branched on so the that the fact that
                       a feasible solution was found in presolve can be recorded*/

                   case IP_FEASIBLE:
                     can->termcode[j] = LP_OPT_FEASIBLE;
                     can->feasible[j] = TRUE;
                     if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
                        can->solutions[j] = (double *) malloc (DSIZE*
                              lp_data->n);
                        memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
                     }
                     break;

                   case IP_FEASIBLE_BUT_CONTINUE:
                     can->termcode[j] = LP_OPT_FEASIBLE_BUT_CONTINUE;
                     can->feasible[j] = TRUE;
                     if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
                        can->solutions[j] = (double *) malloc (DSIZE*
                              lp_data->n);
                        memcpy(can->solutions[j], lp_data->x, DSIZE*lp_data->n);
                     }
                     break;

                   default:
                     break;
                  }
               }
#ifdef COMPILE_FRAC_BRANCHING
               else{
                  if (can->termcode[j] != LP_ABANDONED)
                     get_x(lp_data);
               }
               if (can->termcode[j] != LP_ABANDONED){
                  xind = lp_data->tmp.i1; /* n */
                  xval = lp_data->tmp.d; /* n */
                  can->frac_num[j] = collect_fractions(p, lp_data->x, xind, xval);
                  if (can->frac_num[j] > 0){
                     can->frac_ind[j] = (int *) malloc(can->frac_num[j] * ISIZE);
                     can->frac_val[j] = (double *) malloc(can->frac_num[j]*DSIZE);
                     memcpy(can->frac_ind[j], xind, can->frac_num[j] * ISIZE);
                     memcpy(can->frac_val[j], xval, can->frac_num[j] * DSIZE);
                  }
               }else{
                  can->frac_num[j] = 0;
               }
#endif
#ifdef STATISTICS
               if (can->termcode[j] == LP_D_ITLIM)
                  itlim++;
#endif
            }
            cut = rows[branch_row].cut;
            change_row(lp_data, branch_row, cut->sense, cut->rhs, cut->range);
            free_row_set(lp_data, 1, &branch_row);
            break;
         }

         switch ((j = compare_candidates_u(p, oldobjval, best_can, can))){
          case FIRST_CANDIDATE_BETTER:
          case FIRST_CANDIDATE_BETTER_AND_BRANCH_ON_IT:
            if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
               min_ind = -1;
               for (k = can->child_num - 1; k >= 0; k--){
                  if (can->feasible[k]){
                     if (min_ind < 0){
                        min_obj = SYM_INFINITY;
                        for (l = best_can->child_num - 1; l >= 0; l--){
                           if (best_can->feasible[l] && best_can->objval[k] <
                                 min_obj){
                              min_obj = best_can->objval[l];
                              min_ind = l;
                           }
                        }
                     }
                     if (min_ind > -1){
                        if(can->objval[k] > best_can->objval[min_ind]){
                           best_can->feasible[k] = TRUE;
                           best_can->solutions[k] = can->solutions[k];
                           can->solutions[k] = 0;
                           min_ind = -1;
                        }
                     }
                  }
               }
            } else{
               for (k = best_can->child_num - 1; k >= 0; k--){
                  /* Again, this is only for tracking that there was a feasible
                     solution discovered in presolve for display purposes */
                  if (can->feasible[k]){
                     best_can->feasible[k] = TRUE;
                  }
               }
            }
            free_candidate(candidates + i);
            break;
          case SECOND_CANDIDATE_BETTER:
          case SECOND_CANDIDATE_BETTER_AND_BRANCH_ON_IT:
#ifndef MAX_CHILDREN_NUM
            if (best_can == NULL){
               pobj  = objval;
               pterm = termcode;
               pfeas = feasible;
               piter = iterd;
#ifdef COMPILE_FRAC_BRANCHING
               pfrnum = frnum;
               pfrind = frind;
               pfrval = frval;
#endif
            }else{
               pobj  = best_can->objval;
               pterm = best_can->termcode;
               pfeas = best_can->feasible;
               piter = best_can->iterd;
#ifdef COMPILE_FRAC_BRANCHING
               pfrnum = best_can->frac_num;
               pfrind = best_can->frac_ind;
               pfrval = best_can->frac_val;
#endif
            }
#endif
            if (best_can){
               if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
                  min_ind = -1;
                  for (k = best_can->child_num - 1; k >= 0; k--){
                     if (best_can->feasible[k]){
                        if (min_ind < 0){
                           min_obj = SYM_INFINITY;
                           for (l = can->child_num - 1; l >= 0; l--){
                              if (can->feasible[l] && can->objval[k] <
                                    min_obj){
                                 min_obj = can->objval[l];
                                 min_ind = l;
                              }
                           }
                        }
                        if (min_ind > -1){
                           if(best_can->objval[k] > can->objval[min_ind]){
                              can->feasible[k] = TRUE;
                              can->solutions[k] = best_can->solutions[k];
                              best_can->solutions[k] = 0;
                              min_ind = -1;
                           }
                        }
                     }
                  }
               }
               else{
                  for (k = can->child_num - 1; k >= 0; k--){
                     /* Again, this is only for tracking that there was a feasible
                        solution discovered in presolve for display purposes */
                     if (best_can->feasible[k]){
                        can->feasible[k] = TRUE;
                     }
                  }
               }
               free_candidate(&best_can);
            }
            best_can = can;
            candidates[i] = NULL;
            break;
         }
         if ((j & BRANCH_ON_IT)){
            break;
         }
         st_time += used_time(&total_time);
         should_continue_strong_branching(p,i,cand_num,st_time,total_iters,
               &should_continue);
         if (should_continue==FALSE) {
            PRINT(p->par.verbosity, 0,
                  ("too much time in strong branching, breaking\n"));
            break;
         }
      }
   }
   //printf ("total_iters = %d \n",total_iters);
   //printf ("candidates evaluated = %d \n",i);

   if (should_use_hot_starts) {
      unmark_hotstart(lp_data);
   }

#if 0
   if (best_can->type == CANDIDATE_VARIABLE &&
       vars[best_can->position]->lb == vars[best_can->position]->ub){
      printf("Error -- branching variable is already fixed. \n");
      printf("SYMPHONY has encountered numerical difficulties \n");
      printf("with the LP solver. Exiting...\n\n");
   }

   if (best_can->type == CANDIDATE_VARIABLE &&
       best_can->rhs[0] < vars[best_can->position]->lb ||
       best_can->rhs[1] > vars[best_can->position]->ub){
      printf("Warning -- possible illegal branching. \n");
      printf("SYMPHONY has encountered possible numerical difficulties \n");
      printf("with the LP solver. Exiting...\n\n");
   }
#endif

#ifndef MAX_CHILDREN_NUM
   FREE(pobj); FREE(pterm); FREE(pfeas); FREE(piter);
#  ifdef COMPILE_FRAC_BRANCHING
   FREE(pfrnum); FREE(pfrind); FREE(pfrval);
#  endif
#endif

   if (p->par.max_presolve_iter > 0) {
      if (should_use_hot_starts == TRUE) {
         set_itlim_hotstart(lp_data, -1);
      }
      set_itlim(lp_data, -1); // both limits should be set if using hotstarts
   }

#ifdef STATISTICS
   PRINT(p->par.verbosity, 5,
	 ("Itlim reached %i times out of %i .\n\n", itlim, cnum));
#endif

   if (should_use_rel_br != TRUE) {
      for (i++; i<cand_num; i++){
         /* Free the remaining candidates */
         free_candidate(candidates + i);
      }
   }
   FREE(candidates);

   if (p->par.keep_description_of_pruned == KEEP_IN_MEMORY){
      indices = lp_data->tmp.i1;
      values = lp_data->tmp.d;
      for (k = best_can->child_num - 1; k >= 0; k--){
	 if (best_can->feasible[k]){
	    best_can->sol_sizes[k] = collect_nonzeros(p,
						      best_can->solutions[k],
						     indices, values);
	    FREE(best_can->solutions[k]);
	    best_can->sol_inds[k] = (int *) malloc(best_can->sol_sizes[k]*
						   ISIZE);
	    best_can->solutions[k] = (double *) malloc(best_can->sol_sizes[k]*
						       DSIZE);
	    memcpy(best_can->sol_inds[k], indices, best_can->sol_sizes[k] *
		   ISIZE);
	    memcpy(best_can->solutions[k], values, best_can->sol_sizes[k]*
		   DSIZE);
	    break;
	 }
      }
   }

   *candidate = best_can;

   return(DO_BRANCH);
}

/*===========================================================================*/

int branch(lp_prob *p, int cuts)
{
   LPdata *lp_data = p->lp_data;
   branch_obj *can;
   char *action;
   int branch_var, branch_row, keep;
   var_desc *var;
   cut_data *cut;
   node_desc *desc;
   int termcode;

   termcode = select_branching_object(p, &cuts, &can);

   if (termcode == ERROR__NO_BRANCHING_CANDIDATE){
      return(termcode);
   }

   if (can == NULL){
      if (termcode == DO_NOT_BRANCH__FEAS_SOL) {
         /* a better feasible solution was found. return to do reduced cost
          * fixing etc.
          */
         return(FEAS_SOL_FOUND);
      }
      /* We were either able to fathom the node or found violated cuts
       * In any case, send the qualifying cuts to the cutpool */
      p->comp_times.strong_branching += used_time(&p->tt);
#pragma omp critical(cut_pool)
      send_cuts_to_pool(p, p->par.eff_cnt_before_cutpool);
      p->comp_times.communication += used_time(&p->tt);
      return (termcode == DO_NOT_BRANCH__FATHOMED ? BRANCHING_INF_NODE : cuts);
   }

   /*------------------------------------------------------------------------*\
    * Now we evaluate can, the best of the candidates.
   \*------------------------------------------------------------------------*/
   action = lp_data->tmp.c; /* n (good estimate... can->child_num) */
   if ((termcode = select_child_u(p, can, action)) < 0)
      return(termcode);
   if (p->par.verbosity > 4)
      print_branch_stat_u(p, can, action);

   for (keep = can->child_num-1; keep >= 0; keep--)
      if (action[keep] == KEEP_THIS_CHILD) break;

   /* Send the branching information to the TM and inquire whether we
      should dive */
   p->comp_times.strong_branching += used_time(&p->tt);
   send_branching_info(p, can, action, &keep);
   p->comp_times.communication += used_time(&p->tt);

   /* If we don't dive then return quickly */
   if (keep < 0 || p->dive == DO_NOT_DIVE){
      free_candidate_completely(&can);
      return(FATHOMED_NODE);
   }

   desc = p->desc;
   switch (can->type){
    case CANDIDATE_VARIABLE:
      var = lp_data->vars[branch_var = can->position];
      switch (can->sense[keep]){
       case 'E':
	 var->new_lb = var->new_ub = can->rhs[keep];
	 var->lb = var->ub = can->rhs[keep];                             break;
       case 'R':
	 var->new_lb = can->rhs[keep];
         var->new_ub = var->lb + can->range[keep];
	 var->lb = can->rhs[keep]; var->ub = var->lb + can->range[keep]; break;
       case 'L':
	 var->new_ub = can->rhs[keep];
	 var->ub = can->rhs[keep];                                       break;
       case 'G':
	 var->new_lb = can->rhs[keep];
	 var->lb = can->rhs[keep];                                       break;
      }
      change_col(lp_data, branch_var, can->sense[keep], var->lb, var->ub);
      lp_data->status[branch_var] |= VARIABLE_BRANCHED_ON;
      break;
    case CANDIDATE_CUT_IN_MATRIX:
      branch_row = can->position;
      cut = lp_data->rows[branch_row].cut;
      /* To maintain consistency with TM we have to fix a few more things if
	 we had a non-base, new branching cut */
      if (branch_row >= p->base.cutnum && !(cut->branch & CUT_BRANCHED_ON)){
	 /* insert cut->name into p->desc.cutind.list, and insert SLACK_BASIC
	    to the same position in p->desc.basis.extrarows.stat */
#ifdef DO_TESTS
	 if (desc->cutind.size != desc->basis.extrarows.size){
	    printf("Oops! desc.cutind.size != desc.basis.extrarows.size! \n");
	    exit(-123);
	 }
#endif
#ifdef COMPILE_IN_LP
	 /* Because these cuts are shared with the treemanager, we have to
	    make a copy before changing them if the LP is compiled in */
	 cut = (cut_data *) malloc(sizeof(cut_data));
	 memcpy((char *)cut, (char *)lp_data->rows[branch_row].cut,
		sizeof(cut_data));
	 if (cut->size){
	    cut->coef = (char *) malloc(cut->size);
	    memcpy((char *)cut->coef,
		   (char *)lp_data->rows[branch_row].cut->coef, cut->size);
	 }
	 lp_data->rows[branch_row].cut = cut;
#endif
	 if (desc->cutind.size == 0){
	    desc->cutind.size = 1;
	    desc->cutind.list = (int *) malloc(ISIZE);
	    desc->cutind.list[0] = cut->name;
	    desc->basis.extrarows.size = 1; /* this must have been 0, too */
	    desc->basis.extrarows.stat = (int *) malloc(ISIZE);
	    desc->basis.extrarows.stat[0] = SLACK_BASIC;
	 }else{
	    int i, name = cut->name;
	    int *list;
	    int *stat;
	    /* most of the time the one extra element will fit into the
	       already allocated memory chunk, so it's worth to realloc */
	    desc->cutind.size++;
	    list = desc->cutind.list =
	       (int *) realloc(desc->cutind.list, desc->cutind.size * ISIZE);
	    desc->basis.extrarows.size++;
	    stat = desc->basis.extrarows.stat =
	       (int *) realloc(desc->basis.extrarows.stat,
			       desc->cutind.size * ISIZE);
	    for (i = desc->cutind.size - 1; i > 0; i--){
#ifdef DO_TESTS
	       if (name == list[i-1]){
		  printf("Oops! name == desc.cutind.list[i] !\n");
		  exit(-124);
	       }
#endif
	       if (name < list[i-1]){
		  list[i] = list[i-1];
		  stat[i] = stat[i-1];
	       }else{
		  break;
	       }
	    }
	    list[i] = name;
	    stat[i] = SLACK_BASIC;
	 }
      }
      cut->rhs = can->rhs[keep];
      if ((cut->sense = can->sense[keep]) == 'R')
	 cut->range = can->range[keep];
      cut->branch = CUT_BRANCHED_ON | can->branch[keep];
      constrain_row_set(lp_data, 1, &branch_row);
      lp_data->rows[branch_row].free = FALSE;
      break;
   }

   /* Since this is a child we dived into, we know that TM stores the stati of
      extra vars/rows wrt the parent */
   p->desc->basis.extravars.type = WRT_PARENT;
   p->desc->basis.extrarows.type = WRT_PARENT;

   free_candidate_completely(&can);

   /* the new p->bc_index is received in send_branching_info() */
   p->bc_level++;
   /*p->iter_num = 0;*/

   return(NEW_NODE);
}

/*===========================================================================*/

int col_gen_before_branch(lp_prob *p, int *new_vars)
{
   our_col_set *new_cols;
   int dual_feas;

   check_ub(p);
   if (! p->has_ub ||
       (p->colgen_strategy & BEFORE_BRANCH__DO_NOT_GENERATE_COLS) ||
       (p->lp_data->nf_status & NF_CHECK_NOTHING))
      return(DO_BRANCH);

   PRINT(p->par.verbosity, 2, ("Generating cols before branching.\n"));
   p->comp_times.strong_branching += used_time(&p->tt);
   new_cols = price_all_vars(p);
   p->comp_times.pricing += used_time(&p->tt);
   /*price_all_vars sorts by user_ind. We need things sorted by colind */
   colind_sort_extra(p);
   *new_vars = new_cols->num_vars + new_cols->rel_ub + new_cols->rel_lb;
   dual_feas = new_cols->dual_feas;
   free_col_set(&new_cols);
   check_ub(p);
   if (dual_feas == NOT_TDF){
      return(DO_NOT_BRANCH);
   }else{
      if (p->ub - p->par.granularity < p->lp_data->objval ||
	  p->lp_data->termcode == LP_D_OBJLIM ||
	  p->lp_data->termcode == LP_OPT_FEASIBLE){
	 /* If total dual feas and high cost or feasibility ==> fathomable */
	 PRINT(p->par.verbosity, 1, ("Managed to fathom the node.\n"));
	 send_node_desc(p, p->lp_data->termcode == LP_OPT_FEASIBLE ?
			FEASIBLE_PRUNED : OVER_UB_PRUNED);
	 p->comp_times.communication += used_time(&p->tt);
	 return(DO_NOT_BRANCH__FATHOMED);
      }else{
	 return(DO_BRANCH); /* if we got here, then DO_BRANCH */
      }
   }
   return(DO_BRANCH); /* fake return */
}

/*===========================================================================*/

/*****************************************************************************/
/* This is a generic function                                                */
/*****************************************************************************/

void branch_close_to_half(lp_prob *p, int max_cand_num, int *cand_num,
			  branch_obj ***candidates)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol100 = lp_data->lpetol*100, lpetol1 = 1 - lpetol100;
   int *xind = lp_data->tmp.i1; /* n */
   double fracx, *xval = lp_data->tmp.d; /* n */
   branch_obj *cand;
   int i, j, cnt = 0;
   double lim[7] = {.1, .15, .20, .233333, .266667, .3, 1};
   var_desc **vars = lp_data->vars;
   int should_use_rel_br = p->par.should_use_rel_br;

   /* first get the fractional values */
   if (should_use_rel_br == TRUE) {
      xind = p->br_rel_cand_list;
   }
   for (i = lp_data->n-1; i >= 0; i--){
      if (vars[i]->is_int){
	 if (x[i] > vars[i]->new_lb && x[i] < vars[i]->new_ub){
	    fracx = x[i] - floor(x[i]);
	    if (fracx > lpetol100 && fracx < lpetol1){
	       xind[cnt] = i;
	       xval[cnt++] = fabs(fracx - .5);
	    }
	 }
      }
      *cand_num = cnt;
   }

   if (should_use_rel_br == TRUE) {
      *candidates = (branch_obj **) malloc(1 * sizeof(branch_obj *));
      cand = (*candidates)[0] = (branch_obj *) calloc(1, sizeof(branch_obj) );
      cand->type = CANDIDATE_VARIABLE;
      cand->child_num = 2;
      cand->sense[0] = 'L';
      cand->sense[1] = 'G';
      cand->range[0] = cand->range[1] = 0;
      qsort_di(xval, xind, cnt);
   } else {
      qsort_di(xval, xind, cnt);
      if (p->bc_level>p->par.strong_br_all_candidates_level ||
            p->par.user_set_strong_branching_cand_num) {
         for (j = 0, i = 0; i < cnt;){
            if (xval[i] > lim[j]){
               if (i == 0){
                  j++; continue;
               }else{
                  break;
               }
            }else{
               i++;
            }
         }
         cnt = i;
         *cand_num = MIN(max_cand_num, cnt);
      } else {
         *cand_num = cnt;
      }

      if (!*candidates)
         *candidates = (branch_obj **) malloc(*cand_num * sizeof(branch_obj *));
      for (i=*cand_num-1; i>=0; i--){
         cand = (*candidates)[i] = (branch_obj *) calloc(1, sizeof(branch_obj) );
         cand->type = CANDIDATE_VARIABLE;
         cand->child_num = 2;
         cand->position = xind[i];
         cand->sense[0] = 'L';
         cand->sense[1] = 'G';
         cand->rhs[0] = floor(x[xind[i]]);
         cand->rhs[1] = cand->rhs[0] + 1;
         cand->range[0] = cand->range[1] = 0;
      }
   }
}

/*===========================================================================*/

/*****************************************************************************/
/* This is a generic function                                                */
/*****************************************************************************/

void branch_close_to_half_and_expensive(lp_prob *p, int max_cand_num,
					int *cand_num, branch_obj ***candidates)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol = lp_data->lpetol, lpetol1 = 1 - lpetol;
   int *xind = lp_data->tmp.i1; /* n */
   double fracx, *xval = lp_data->tmp.d; /* n */
   branch_obj *cand;
   int i, j, cnt = 0;
   double lim[7] = {.1, .15, .20, .233333, .266667, .3, 1};

   /* first get the fractional values */
   for (i = lp_data->n-1; i >= 0; i--){
      fracx = x[i] - floor(x[i]);
      if (fracx > lpetol && fracx < lpetol1){
	 xind[cnt] = i;
	 xval[cnt++] = fabs(fracx - .5);
       }
   }

   qsort_di(xval, xind, cnt);

   for (j=0, i=0; i<cnt; i++){
      if (xval[i] > lim[j]){
	 if (i == 0){
	    j++; continue;
	 }else{
	    break;
	 }
      }
   }
   cnt = i;

   if (max_cand_num >= cnt){
      *cand_num = cnt;
   }else{
      for (i=cnt-1; i>=0; i--){
	 get_objcoef(p->lp_data, xind[i], xval+i);
	 xval[i] *= -1;
      }
      qsort_di(xval, xind, cnt);
      *cand_num = max_cand_num;
   }

   if (!*candidates)
      *candidates = (branch_obj **) malloc(*cand_num * sizeof(branch_obj *));
   for (i=*cand_num-1; i>=0; i--){
      cand = (*candidates)[i] = (branch_obj *) calloc(1, sizeof(branch_obj) );
      cand->type = CANDIDATE_VARIABLE;
      cand->child_num = 2;
      cand->position = xind[i];
      cand->sense[0] = 'L';
      cand->sense[1] = 'G';
      cand->rhs[0] = floor(x[xind[i]]);
      cand->rhs[1] = cand->rhs[0] + 1;
      cand->range[0] = cand->range[1] = 0;
   }
}

/*===========================================================================*/

/*****************************************************************************/
/* This works only for 0/1 problems!!!                                       */
/*****************************************************************************/

void branch_close_to_one_and_cheap(lp_prob *p, int max_cand_num, int *cand_num,
				   branch_obj ***candidates)
{
   LPdata *lp_data = p->lp_data;
   double *x = lp_data->x;
   double lpetol = lp_data->lpetol, lpetol1 = 1 - lpetol;
   int *xind = lp_data->tmp.i1; /* n */
   double *xval = lp_data->tmp.d; /* n */
   branch_obj *cand;
   int i, j, cnt = 0;
   double lim[8] = {.1, .2, .25, .3, .333333, .366667, .4, 1};

   /* first get the fractional values */
   for (i = lp_data->n-1; i >= 0; i--){
      if (x[i] > lpetol && x[i] < lpetol1){
	 xind[cnt] = i;
	 xval[cnt++] = 1 - x[i];
      }
   }
   qsort_di(xval, xind, cnt);

   for (j=0, i=0; i<cnt; i++){
      if (xval[i] > lim[j]){
	 if (i == 0){
	    j++; continue;
	 }else{
	    break;
	 }
      }
   }
   cnt = i;

   if (max_cand_num >= cnt){
      *cand_num = cnt;
   }else{
      for (i=cnt-1; i>=0; i--){
	 get_objcoef(p->lp_data, xind[i], xval+i);
      }
      qsort_di(xval, xind, cnt);
      *cand_num = max_cand_num;
   }

   if (!*candidates)
      *candidates = (branch_obj **) malloc(*cand_num * sizeof(branch_obj *));
   for (i=*cand_num-1; i>=0; i--){
      cand = (*candidates)[i] = (branch_obj *) calloc(1, sizeof(branch_obj) );
      cand->type = CANDIDATE_VARIABLE;
      cand->child_num = 2;
      cand->position = xind[i];
      cand->sense[0] = 'L';
      cand->sense[1] = 'G';
      cand->rhs[0] = floor(x[xind[i]]);
      cand->rhs[1] = cand->rhs[0] + 1;
      cand->range[0] = cand->range[1] = 0;
   }
}

/*===========================================================================*/
int should_continue_strong_branching(lp_prob *p, int i, int cand_num,
                                     double st_time, int total_iters,
                                     int *should_continue)
{
   double allowed_time = 0;
   *should_continue = TRUE;
   int min_cands;
   int verbosity = p->par.verbosity;
   //verbosity = 20;
   if (p->bc_level<1) {
      allowed_time = 20*p->comp_times.lp/p->iter_num;
      //allowed_iter = 20*p->lp_stat.lp_total_iter_num/(p->iter_num + 1);
      if (allowed_time < 2) {
	 allowed_time = 2;
      }
      //allowed_iter = MAX(allowed_iter, 1000);
      min_cands = MIN(cand_num,p->par.strong_branching_cand_num_max);
   } else {
      allowed_time = p->comp_times.lp/2 - p->comp_times.strong_branching;
      //allowed_iter = (int)((p->lp_stat.lp_total_iter_num -
      //		    p->lp_stat.str_br_total_iter_num)/2.0);
      min_cands = MIN(cand_num,p->par.strong_branching_cand_num_min);
   }
   PRINT(verbosity,10,("allowed_time = %f\n",allowed_time));
   if (st_time/(i+1)*cand_num < allowed_time) {
      /* all cands can be evaluated in given time */
      *should_continue = TRUE;
   } else if (i >= min_cands-1 && st_time>allowed_time) {
      /* time is up and min required candidates have been evaluated */
      *should_continue = FALSE;
   } else if (p->par.user_set_max_presolve_iter == TRUE) {
      /* user specified a limit and we wont change it */
      *should_continue = TRUE;
   } else {
      /* we will not be able to evaluate all candidates in given time. we
       * reduce the number of iterations */
      double min_iters =
         (allowed_time-st_time)*total_iters/st_time/(cand_num-i+1);
      if (min_iters<10) {
         /*
          * cant evaluate all candidates in given time with just 10 iters
          * as well. we have a choice: increase iters and do min_cands
          * or do ten iters and try to evaluate max possible num of cands.
          * we like the second option more.
          */
         min_iters = 10;
      }
      if (p->par.use_hot_starts && !p->par.branch_on_cuts) {
         set_itlim_hotstart(p->lp_data, (int) min_iters);
         set_itlim(p->lp_data, (int) min_iters);
      } else {
         set_itlim(p->lp_data, (int) min_iters);
      }
      PRINT(verbosity,6, ("iteration limit set to %d\n", (int )min_iters));
      *should_continue = TRUE;
   }
   PRINT(verbosity,29, ("strong branching i = %d\n",i));
   return 0;
}

/*===========================================================================*/
int strong_branch(lp_prob *p, int branch_var, double lb, double ub,
		  double new_lb, double new_ub, double *obj, int should_use_hot_starts,
                  int *termstatus, int *iterd)
{
   int status = 0;
   LPdata *lp_data = p->lp_data;

   // TODO: LP_ABANDONED
   /* change the lb and ub */
   change_lbub(lp_data, branch_var, new_lb, new_ub);


   //   if (p->par.use_hot_starts && !p->par.branch_on_cuts) {
   if (should_use_hot_starts) {
      *termstatus = solve_hotstart(lp_data, iterd);
   } else {
      *termstatus = dual_simplex(lp_data, iterd);
   }

   if (*termstatus == LP_D_INFEASIBLE || *termstatus == LP_D_OBJLIM ||
         *termstatus == LP_D_UNBOUNDED) {
      *obj = SYM_INFINITY;
      p->lp_stat.str_br_bnd_changes++;
   } else {
      *obj = lp_data->objval;
      if (*termstatus == LP_OPTIMAL) {
         if (!p->has_ub || *obj < p->ub - lp_data->lpetol) {
            is_feasible_u(p, TRUE, FALSE);
         } else {
            *obj = SYM_INFINITY;
            p->lp_stat.str_br_bnd_changes++;
         }
      } else if (*termstatus == LP_ABANDONED) {
         status = LP_ABANDONED;
      }
   }
   p->lp_stat.lp_calls++;
   p->lp_stat.str_br_lp_calls++;
   p->lp_stat.str_br_total_iter_num += *iterd;
   p->lp_stat.num_str_br_cands_in_path++;

   change_lbub(lp_data, branch_var, lb, ub);
   return status;
}
/*===========================================================================*/
/*===========================================================================*/
