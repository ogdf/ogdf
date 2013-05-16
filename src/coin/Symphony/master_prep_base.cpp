/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Menal Guzelsoy                                 */
/*                                                                           */
/* (c) Copyright 2006-2011 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include <memory.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sym_types.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_master.h"
#include "sym_prep.h"

/*===========================================================================*/
/*
  This function is the internal master of the preprocessing part.
*/
/*===========================================================================*/
int prep_basic(PREPdesc *P)
{

   /* initialization */
   prep_stats *stats = &(P->stats);
   prep_params params = P->params;

   const int p_level = params.level;
   const int iter_cnt_limit = params.iteration_limit;
   const int verbosity = params.verbosity;
   const int do_sr_rlx = params.do_single_row_rlx;
   const int dive_level = params.dive_level;
   const int impl_dive_level = params.impl_dive_level;
   const double etol = params.etol;
   const double time_limit = params.time_limit;

   int termcode;     /* return status of functions called herein */
   int iter_cnt = 0;
   double a_val, new_bound = 0.0;// min_ub, max_lb;
   int i, j, m, n, nz, *r_matbeg, *r_matind, *matbeg, *matind;
   double *obj, *rhs, *r_matval, *matval, *ub, *lb;
   char *sense;
   int col_ind, row_ind, fix_type;

   char can_impl = FALSE, bin_type = FALSE;
   int old_changes_cnt, changes_diff = 0, new_changes_cnt, init_changes_cnt = 0;
   int old_others_cnt, new_others_cnt, mark_others_cnt = 0;
   double start_impl_time, mark_time, impl_time = 0.0;

   MIPdesc *mip = P->mip;
   MIPinfo *mip_inf = mip->mip_inf;
   COLinfo * cols = mip_inf->cols;
   ROWinfo *rows = mip_inf->rows;

   stats->row_infeas_ind = stats->col_infeas_ind = -1;

   if (mip_inf->prob_type == CONTINUOUS_TYPE){
      /* no need for prep, just quit */
      return PREP_UNMODIFIED;
   }

   /*
    * round up/down the bounds on integer variables.
    * if p_level>2, also round up/down the bounds on vars that are implied to
    * be integers.
    */
   termcode = prep_integerize_bounds(P);
   if (PREP_QUIT(termcode)){
      return termcode;
   }

   m = mip->m;
   n = mip->n;
   nz = mip->nz;

   r_matbeg = mip->row_matbeg;
   r_matind = mip->row_matind;
   r_matval = mip->row_matval;

   matbeg = mip->matbeg;
   matind = mip->matind;
   matval = mip->matval;

   ub = mip->ub;
   lb = mip->lb;

   obj = mip->obj;
   sense = mip->sense;
   rhs = mip->rhs;

   char need_reset;

   /* check if we have binary columns */
   if (mip_inf->prob_type == BINARY_TYPE ||
      mip_inf->prob_type == BIN_CONT_TYPE ||
      mip_inf->prob_type == BIN_INT_TYPE ||
      mip_inf->prob_type == ALL_MIXED_TYPE){

      bin_type = TRUE;
   }

   /* first check duplicate rows, cols */
   termcode = prep_delete_duplicate_rows_cols(P, TRUE, TRUE);
   if (PREP_QUIT(termcode)){
      return termcode;
   }

   init_changes_cnt = stats->vars_fixed + stats->rows_deleted;

   /*initialize implications data structure */
   if (p_level >= 5){ /* disabled now */
      /* probably doesnt worth it for this version if n or nz is too large */
      if (p_level >= 10 || (n < 1e4 && nz < 1e5)){
	 if (bin_type){
	    /* for now, just between binary variables */
	    P->impl_rows = (ROWinfo *)malloc(sizeof(ROWinfo)*m);
	    P->impl_cols = (COLinfo *)malloc(sizeof(COLinfo)*n);
	    P->impl_ub = (double *) malloc(DSIZE*n);
	    P->impl_lb = (double *) malloc(DSIZE*n);

	    P->ulist_checked = (char *)malloc(CSIZE * n);
	    P->llist_checked = (char *)malloc(CSIZE * n);

	    P->impl_limit = params.impl_limit;
	    can_impl = TRUE;

	    /* get the list of columns to apply impl on */
	    /* not effective */
	    P->impl_vars = (char *) calloc(CSIZE,n);

	    for (i = n - 1; i >= 0; i--){
	       P->impl_vars[i] = TRUE;
	    }
	 }
      }
   }

   /* main preprocessing loop */
   old_changes_cnt = new_changes_cnt = 0;
   old_others_cnt = new_others_cnt = 0;
   while(iter_cnt < iter_cnt_limit){

      iter_cnt++;

      PRINT(verbosity, 2, ("Basic iteration number: %d\n", iter_cnt));

      /* check the updated bounds and cols to iterate on*/
      /* for each updated row and column do bla bla...*/
      /* while iterating update the new cols and rows indices*/

      /*=====================================================================*/
      /*=====================================================================*/

      for (col_ind = 0; col_ind < n; col_ind++){
	 if (cols[col_ind].var_type == 'F'){
	    continue;
	 }
	 /* can we fix it? first check implications */
	 /* disabled now */
	 if (can_impl && impl_time < time_limit){
	    if (cols[col_ind].var_type == 'B' &&
	       P->impl_vars[col_ind] && (iter_cnt < 2 ||
					 (iter_cnt > 1 &&
					  changes_diff > 0))){
	       /* fist copy initial info */
	       /* do once for each variable */
	       start_impl_time = wall_clock(NULL);
	       fix_type = FIX_NO_BOUND;

	       need_reset = FALSE;
	       mark_time = wall_clock(NULL);
	       memcpy(P->impl_rows, rows, sizeof(ROWinfo)*m);
	       memcpy(P->impl_cols, cols, sizeof(COLinfo)*n);
	       memcpy(P->impl_ub, ub, DSIZE*n);
	       memcpy(P->impl_lb, lb, DSIZE*n);
	       P->impl_stats = P->stats;
	       P->alloc_time += wall_clock(NULL) - mark_time;

	       if (cols[col_ind].sign_type != ALL_NEG_VEC){
		  need_reset = TRUE;
		  //if (!cols[col_ind].ulist){
		     // P->impl_cols[col_ind].ulist = (cols[col_ind].ulist =
		     //	(IMPlist *)calloc(sizeof(IMPlist),1));
		  //}

		  P->list = cols[col_ind].ulist;
		  /* fix it to 1.0 and see if that causes any infeasibility
		     otherwise get the impllist and continue*/
		  /* get the implication list */

		  /* fix this column, update row bounds of this column
		     check for redundancy, */
		  P->impl_col_ind = col_ind;
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind, -1,
							    impl_dive_level,
							    1.0,
							    FIX_BINARY, TRUE,
							    TRUE);
		  free_imp_list(&(cols[col_ind].ulist));
		  P->impl_cols[col_ind].ulist = 0;
		  P->list = 0;
		  if (termcode == PREP_INFEAS){
		     /*then this column is fixable to its lower bound! */
		     new_bound = 0.0;
		     fix_type = FIX_BINARY;
		  }
	       }

	       if (fix_type != FIX_BINARY &&
		  cols[col_ind].sign_type != ALL_POS_VEC){
		  mark_time = wall_clock(NULL);
		  /* reset what we had */
		  if (need_reset){
		     memcpy(rows, P->impl_rows,sizeof(ROWinfo)*mip->m);
		     memcpy(cols, P->impl_cols, sizeof(COLinfo)*mip->n);

		     memcpy(ub, P->impl_ub, DSIZE*mip->n);
		     memcpy(lb, P->impl_lb, DSIZE*mip->n);
		     P->stats = P->impl_stats;
		     P->alloc_time += wall_clock(NULL) - mark_time;
		  }
		  //if (!cols[col_ind].llist){
		     //P->impl_cols[col_ind].llist = (cols[col_ind].llist =
		     //	(IMPlist *)calloc(sizeof(IMPlist),1));
		  //}

		  P->list = cols[col_ind].llist;
		  P->impl_col_ind = col_ind;

		  termcode = prep_modified_cols_update_info(P, 1, &col_ind, -1,
							    impl_dive_level,
							    0.0, FIX_BINARY,
							    TRUE, TRUE);
		  free_imp_list(&(cols[col_ind].llist));
		  P->impl_cols[col_ind].llist = 0;
		  P->list = 0;
		  if (termcode == PREP_INFEAS){
		     new_bound = 1.0;
		     fix_type = FIX_BINARY;
		  }
	       }

	       /* now get back */
	       mark_time = wall_clock(NULL);
	       memcpy(rows, P->impl_rows,sizeof(ROWinfo)*mip->m);
	       memcpy(cols, P->impl_cols, sizeof(COLinfo)*mip->n);
	       memcpy(ub, P->impl_ub, DSIZE*mip->n);
	       memcpy(lb, P->impl_lb, DSIZE*mip->n);
	       P->stats = P->impl_stats;
	       P->alloc_time += wall_clock(NULL) - mark_time;
	       /* and now check if we can fix anything */

	       impl_time += wall_clock(NULL) - start_impl_time;
	       if (fix_type == FIX_BINARY){
		  termcode = prep_modified_cols_update_info(P, 1,
							    &col_ind, -1,
							    dive_level,
							    new_bound,
							    fix_type, TRUE,
							    FALSE);
		  if (PREP_QUIT(termcode)){
		     return termcode;
		  }
		  continue;
	       }
	    }
	 }

	 /* couldnt fix it, continue */
	 /* for each coefficient, and the corresponding row and column
	    check if we can:
	    -detect row redundancy, or tighten the rhs
	    -tighten the variable bounds
	    -and finally, tighten the coefficient itself
	 */
	 for (j = matbeg[col_ind]; j < matbeg[col_ind+1]; j++){
	    row_ind = matind[j];

	    if (rows[row_ind].is_redundant){
	       continue;
	    } else {
	       if (rows[row_ind].ub >= INF && rows[row_ind].lb <= -INF){
		  continue;
	       }
	       /* here check redundancy */
	       termcode = prep_check_redundancy(P, row_ind, FALSE, 0.0, 0.0,
						FALSE, 0);
	       if (PREP_QUIT(termcode)){
		  return termcode;
	       }
	       if (rows[row_ind].is_redundant){
		  continue;
	       }
	    }

	    /* so the row is not redundant, continue... */
	    a_val = matval[j];

	    /* check whether it is an unbounded column */
	    if (rows[row_ind].ub_inf_var_num <= 1 ||
	       rows[row_ind].lb_inf_var_num <= 1){

	       if (ub[col_ind] >= INF || lb[col_ind] <= -INF){
		  if ((a_val > etol && ub[col_ind] >= INF) ||
		     (a_val < -etol && lb[col_ind] <= -INF) ||
		     sense[row_ind] == 'E'){

		     /* if so, check if we can put a bound */
		     termcode =
			prep_force_row_bounds(P, row_ind, col_ind, j);

		     if (PREP_QUIT(termcode)){
			return termcode;
		     }

		     if (rows[row_ind].is_redundant){
			continue;
		     }
		  }
	       }
	    }

	    if ((a_val > etol || a_val < -etol) &&
	       cols[col_ind].var_type != 'F'){

	       /* then, try to improve the bounds or the coefficient */
	       termcode = prep_improve_variable(P, col_ind,
						row_ind, j,
						dive_level, TRUE,
						FALSE, FALSE, 0.0,
						0.0, MAT_COL_ORDERED);
	       if (PREP_QUIT(termcode)){
		  return termcode;
	       }
	    }
	 }
      }

      /* now check the changes and decide to iterate again */
      new_changes_cnt = stats->rows_deleted +
	 stats->vars_fixed;

      new_others_cnt = stats->coeffs_changed +
	 stats->bounds_tightened;

      changes_diff = new_changes_cnt - old_changes_cnt;

      if (changes_diff > 0){
	 old_changes_cnt = new_changes_cnt;
	 old_others_cnt = new_others_cnt;
      } else {
	 if (new_others_cnt > old_others_cnt){
	    old_others_cnt = new_others_cnt;
	    mark_others_cnt++;
	 } else {
	    break;
	 }
	 if (mark_others_cnt > 3){
	    break;
	 }
      }
   }


   /* disabled now*/
   /* call single row relaxation presolver */
   if (do_sr_rlx){
      for (row_ind = 0; row_ind < m; row_ind++){
	 if (!rows[row_ind].is_redundant){
	    termcode = prep_solve_sr_rlx(P, 1, &row_ind);
	    if (PREP_QUIT(termcode)){
	       return termcode;
	    }
	 }
      }
   }

   /* if we have deleted some rows, the column sizes might have changed,
      so try to improve those variables one more time */
   if (stats->rows_deleted > 0){
      for (col_ind = 0; col_ind < n; col_ind++){
	 if (cols[col_ind].var_type != 'F' && cols[col_ind].col_size == 0){
	    termcode = prep_improve_variable(P, col_ind, -1, 0,
					     dive_level, TRUE, FALSE,
					     FALSE,
					     0.0,0.0, MAT_COL_ORDERED);
	    if (PREP_QUIT(termcode)){
	       return termcode;
	    }
	 }
      }
   }

   /* similary deleting columns might have ended up with 0 sized rows,
      eliminate them */
   if (stats->vars_fixed > 0){
      for (row_ind = 0; row_ind < m; row_ind++){
	 if (!rows[row_ind].is_redundant &&
	    rows[row_ind].fixed_var_num  >= rows[row_ind].size - 1){
	    termcode = prep_check_redundancy(P, row_ind, FALSE, 0.0, 0.0,
					     FALSE, dive_level);
	    if (PREP_QUIT(termcode)){
	       return termcode;
	    }
	 }
      }
   }

   /* finally, if we have modified the matrix, check duplicacy again */

   if (new_changes_cnt > init_changes_cnt){
      termcode = prep_delete_duplicate_rows_cols(P, TRUE, TRUE);
   }

#if 0
   if (verbosity >= 2){
      printf("total alloc time: %f\n", P->alloc_time);
      printf("total alloc time2: %f\n", P->alloc2_time);
      printf("total impl time2: %f\n", impl_time);
      printf("total impl_cols_time: %f\n", P->impl_cols_time);
      printf("total impl_rows_time: %f\n", P->impl_rows_time);
   }
#endif

   if (new_changes_cnt + stats->coeffs_changed + mip_inf->fixed_var_num > 0){
      termcode = prep_cleanup_desc(P);
   }

   if (PREP_QUIT(termcode)){
      return termcode;
   }

   if (P->mip->mip_inf->binary_sos_row_num){
      prep_sos_fill_var_cnt(P);
   }

   if (stats->rows_deleted +
      stats->vars_fixed +
      stats->bounds_integerized +
      stats->coeffs_changed +
      stats->bounds_tightened > 0){
      return PREP_MODIFIED;
   }

   /* exit preprocessor */

   return PREP_UNMODIFIED;
}

/*===========================================================================*/
/* We check duplicacy here. */
/*===========================================================================*/

int prep_delete_duplicate_rows_cols(PREPdesc *P, char check_rows,
				    char check_cols)
{

   /* start - initialization */
   int termcode = PREP_UNMODIFIED;

   if ((!check_cols && !check_rows) || P->mip->n < 1){
     return termcode;
   } else{
     if(P->mip->m < 2) check_rows = FALSE;
     if(P->mip->n < 2) check_cols = FALSE;
   }

   const double etol = P->params.etol;
   const int verbosity = P->params.verbosity;

   int i, j, k, l, delete_ind, l_ind, r_ind, cr_ind, cl_ind;
   int obj_ind, col_ind, row_ind, end, obj_size, row_size, delete_row_ind;
   char can_iterate, in_conflict, delete_row;
   int fix_type, diff_cnt, diff_ind;
   double new_bound, diff_val, diff_obj_val, diff_row_val, rhs_obj, rhs_row;
   double delete_val;

   MIPdesc *mip = P->mip;
   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   int dive_level = 0;
   prep_stats *stats = &(P->stats);

   int m = mip->m;
   int n = mip->n;

   int *matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;

   double *ub = mip->ub, new_row_ub;
   double *lb = mip->lb, new_row_lb;

   char *sense = mip->sense;
   double *rhs = mip->rhs;
   double *obj = mip->obj;

   double *col_sum = NULL, *col_factor = NULL;
   double *row_sum = NULL, *row_factor = NULL;
   int last_lloc, last_rloc, *r_loc = NULL, *c_loc = NULL;

   int * col_del_ind = NULL, col_del_cnt = 0;
   int * col_fix_type = NULL, dup_type;
   double *col_fix_val = NULL;
   char *col_orig_type = NULL, type_l, type_r, bin_type;
   double obj_l, obj_r, obj_diff;

   /* end - initialization */

   /* to detect duplicate columns we use a very simple but efficient
      hash function: we generate random numbers of size #of rows,
      multiply them with the coefficients of the column vector and
      add them up. So that, if the sum of a column is same with
      one other, then there is a very high chance that they have
      the same column vectors.
   */

   /* same for detecting the duplicacy of rows */

   /* generate the hash table */
   if (check_rows){
      col_factor = (double *)malloc(n*DSIZE);
      row_sum = (double *)calloc(m,DSIZE);
      for (i = 0; i < n; i++){
	 col_factor[i] = 1 + CoinDrand48();
	 if (CoinDrand48() < 0.5) col_factor[i] *= -1.0;
      }

      r_loc = (int *)malloc(m*ISIZE);
      memcpy(r_loc, P->user_row_ind, ISIZE*m);
  }

   if (check_cols){
      col_del_ind = (int *)malloc(n*ISIZE);
      col_fix_type = (int *)malloc(n*ISIZE);
      col_fix_val = (double *)malloc(n*DSIZE);
      col_orig_type = (char *)malloc(n*CSIZE);
      row_factor = (double *)malloc(m*DSIZE);
      col_sum = (double *)calloc(n,DSIZE);

      for (i = 0; i < m; i++){
	 row_factor[i] = 1 + CoinDrand48();
	 if (CoinDrand48() < 0.5) row_factor[i] *= -1.0;
      }

      c_loc = (int *)malloc(n*ISIZE);
      memcpy(c_loc, P->user_col_ind, ISIZE*n);
   }

   if (check_rows && check_cols){
      for (col_ind = 0; col_ind < n; col_ind++){
	 end = matbeg[col_ind + 1];
	 for (j = matbeg[col_ind]; j < end; j++){
	    row_ind = matind[j];
	    row_sum[row_ind] += matval[j]*col_factor[col_ind];
	    col_sum[col_ind] += matval[j]*row_factor[row_ind];
	 }
      }
   } else if (check_rows){
      for (col_ind = 0; col_ind < n; col_ind++){
	 end = matbeg[col_ind + 1];
	 for (j = matbeg[col_ind]; j < end; j++){
	    row_ind = matind[j];
	    row_sum[row_ind] += matval[j]*col_factor[col_ind];
	 }
      }
   } else {
      for (col_ind = 0; col_ind < n; col_ind++){
	 end = matbeg[col_ind + 1];
	 for (j = matbeg[col_ind]; j < end; j++){
	    row_ind = matind[j];
	    col_sum[col_ind] += matval[j]*row_factor[row_ind];
	 }
      }
   }

   /* now, we start with checking the columns */
   /* for same cols,
      -if objs are same, then aggregate them (though for binary,
            check if they conflict first)
      -otherwise,
      -for binary cols, check if they conflict
      -for others, we can delete the duplicate col only under some
      specific requirements
    */

   if (check_cols){
      //qsort_di(col_sum, c_loc, n);
      CoinSort_2(col_sum, col_sum+n, c_loc);
      last_lloc = last_rloc = 0;
      while(TRUE){
	 if (last_lloc == n - 1){
	    break;
	 }
	 /* search for same cols */
	 for (i = last_lloc; i < n - 1; i++){
            /* if fixed or fixable to a bound, continue */
	    if (cols[c_loc[i]].var_type == 'F' ||
	       cols[c_loc[i]].var_type == 'U' ||
	       cols[c_loc[i]].var_type == 'L'){
	       continue;
	    }
	    if (prep_is_equal(col_sum[i], col_sum[i+1], etol)){
	       last_rloc = i+1;
	       if ( i < n - 2 ){
		  for (j = i+2; j < n; j++){
		     if (!prep_is_equal(col_sum[i], col_sum[j], etol)){
			last_rloc = j;
			break;
		     }
		  }
	       }
	       break;
	    }
	 }

	 if (i == n - 1){
	    break;
	 }

	 /* now we got 2 candidate cols*/
	 l_ind = i;
	 r_ind = l_ind + 1;
	 last_lloc = last_rloc;

	 //printf("starting while loop - cols section \n");
	 while(l_ind < last_rloc - 1){

	    cl_ind = c_loc[l_ind];
	    cr_ind = c_loc[r_ind];

	    //printf("processing cl_ind, %i cr_ind %i\n", cl_ind, cr_ind);

	    if (r_ind == last_rloc || cols[cl_ind].var_type == 'F' ||
	       cols[cl_ind].var_type == 'U' ||
	       cols[cl_ind].var_type == 'L'){
	       //cols[cl_ind].var_type != 'B')
	       l_ind++;
	       r_ind = l_ind + 1;
	       continue;
	    }

	    //if (cols[cr_ind].var_type != 'B'){ }
	    if (cols[cr_ind].var_type == 'F' ||
	       cols[cl_ind].var_type == 'U' ||
	       cols[cl_ind].var_type == 'L'){
	       r_ind++;
	       continue;
	    }

	    if (cols[cl_ind].col_size == 0){
	       /* fix this here */
	       col_orig_type[col_del_cnt] = cols[cl_ind].var_type;
	       cols[cl_ind].var_type = 'F';
	       col_del_ind[col_del_cnt] = cl_ind;
	       col_fix_type[col_del_cnt] = FIX_OTHER;
	       if (obj[cl_ind] >= 0){
		  col_fix_val[col_del_cnt] = lb[cl_ind];
	       } else {
		  col_fix_val[col_del_cnt] = ub[cl_ind];
	       }
	       col_del_cnt++;
	       l_ind++;
	       r_ind = l_ind + 1;
	       continue;
	    }

	    if (cols[cr_ind].col_size == 0){
	       /* fix this here */
	       col_orig_type[col_del_cnt] = cols[cr_ind].var_type;
	       cols[cr_ind].var_type = 'F';
	       col_del_ind[col_del_cnt] = cr_ind;
	       col_fix_type[col_del_cnt] = FIX_OTHER;
	       if (obj[cr_ind] >= 0){
		  col_fix_val[col_del_cnt] = lb[cr_ind];
	       } else {
		  col_fix_val[col_del_cnt] = ub[cr_ind];
	       }
	       col_del_cnt++;
	       r_ind++;
	       continue;
	    }

	    /* if sizes are different, then they are definitely not
	       same cols
	       -remember, col_size = real_size - redundant_rows
	    */
	    /* also I dont want to mess with diff type of columns now,
	       so skip if one is int and the other is cont column */
	    type_l = cols[cl_ind].var_type;
	    type_r = cols[cr_ind].var_type;

	    if ((cols[cl_ind].col_size != cols[cr_ind].col_size) ||
	       (type_l == 'C' && type_r != 'C') ||
	       (type_l != 'C' && type_r == 'C')){
	       r_ind++;
	       continue;
	    }

	    /* now we have two cols with same size, but this is not enough,
	       we can iterate only if we have two binary cols or one of
	       the following conditions are satisfied */

	    obj_l = obj[cl_ind];
	    obj_r = obj[cr_ind];
	    obj_diff = obj_l - obj_r;

	    if (type_l == 'B' && type_r == 'B'){
	       bin_type = TRUE;
	    } else {
	       bin_type = FALSE;
	    }

	    dup_type = -1;
    /*
	       dup_type: -1 -> others,
	                  0 -> obj_l = obj_r, obj_diff = 0,

	                  1 -> obj_l, obj_r >= 0, obj_diff > 0, l->lb
                               lb[cl_ind] > -INF, ub[cr_ind] = INF

	                  2 -> obj_l, obj_r <= 0, obj_diff > 0, r->lb
                               lb[cl_ind] = -INF, ub[cr_ind] < INF

	                  3 -> obj_1, obj_r >= 0, obj_diff < 0, r->ub
                               ub[cl_ind] = INF, lb[cr_ind] > -INF

	                  4 -> obj_l, obj_r <= 0, obj_diff < 0, l->ub
                               ub[cl_ind] < INF, lb[cr_ind] = -INF

			      --- unbounded cases --

			  5 -> ub[cl_ind] = INF, lb[cr_ind] = -INF
                               obj_l <= 0, obj_r >= 0

		          6 -> ub[cl_ind] = INF, lb[cr_ind] = -INF
			       obj_l >= 0, obj_r >= 0, obj_diff < 0

		          7 -> ub[cl_ind] = INF, lb[cr_ind] = -INF
			       obj_l <= 0, obj_r <= 0, obj_diff < 0


		          8 -> lb[cl_ind] = -INF, ub[cr_ind] = INF
                               obj_l >= 0, obj_r <= 0

		          9 -> lb[cl_ind] = -INF, ub[cr_ind] = INF
                               obj_l >= 0, obj_r >= 0, obj_diff > 0

		         10 -> lb[cl_ind] = -INF, ub[cr_ind] = INF
                               obj_l <= 0, obj_r <= 0, obj_diff > 0

	    */
	    /* get dup type */

	    if (obj_l == obj_r){
	       /* fixme, make this more efficient */
	       //if (prep_is_equal(obj_r, obj_l, etol)){ }
	       dup_type = 0;
	    } else if (!bin_type){
	       /* at least one inf */
	       /* fix this ugly thing here */
	       if (ub[cl_ind] >= INF || lb[cl_ind] <= -INF ||
		  ub[cr_ind] >= INF || lb[cl_ind] <= -INF){
		  if (obj_diff > 0.0){
		     if (lb[cl_ind] > -INF && ub[cr_ind] >= INF){
			if (obj_l >= 0 && obj_r >= 0) dup_type = 1;
		     } else if (lb[cl_ind] <= -INF && ub[cr_ind] < INF){
			if (obj_l <= 0.0 && obj_r <= 0.0) dup_type = 2;
		     } else if (lb[cl_ind] <= -INF && ub[cr_ind] >= INF){
			if (obj_l >= 0.0 && obj_r <= 0.0) dup_type = 8;
			else if (obj_l >= 0.0 && obj_r >= 0.0) dup_type = 9;
			else if (obj_l <= 0.0 && obj_r <= 0.0) dup_type = 10;
		     }
		  } else {
		     if (ub[cl_ind] >= INF && lb[cr_ind] > -INF){
			if (obj_l >= 0 && obj_r >= 0) dup_type = 3;
		     } else if (ub[cl_ind] < INF && lb[cr_ind] <= -INF){
			if (obj_l <= 0.0 && obj_r <= 0.0) dup_type = 4;
		     } else if (ub[cl_ind] >= INF && lb[cr_ind] <= -INF){
			if (obj_l <= 0.0 && obj_r >= 0.0) dup_type = 5;
			else if (obj_l >= 0.0 && obj_r >= 0.0) dup_type = 6;
			else if (obj_l <= 0.0 && obj_r <= 0.0) dup_type = 7;
		     }
		  }
	       }
	    }

	    /* fixme - check other cases!
	       cant iterate in this case for now */
	    if (dup_type < 0 && !bin_type){
	       r_ind++;
	       continue;
	    }

	    /* so now, we have a case*/
	    /* first check if we have same cols here */
	    /* also check conflict if it is binary */

	    in_conflict = FALSE;
	    can_iterate = TRUE;

	    for ( k = matbeg[cl_ind], l = matbeg[cr_ind];;){
	       if (k < matbeg[cl_ind + 1] &&
		  (matind[k] < matind[l] ||
		   l >= matbeg[cr_ind + 1])){
		  if (!rows[matind[k]].is_redundant){
		     can_iterate = FALSE;
		     break;
		  }
		  k++;
	       } else if (l < matbeg[cr_ind + 1] &&
			(matind[k] > matind[l] ||
			 k >= matbeg[cl_ind+1])){
		  if (!rows[matind[l]].is_redundant){
		     can_iterate = FALSE;
		     break;
		  }
		  l++;
	       } else {
		  if (!rows[matind[l]].is_redundant){
		     if (!prep_is_equal(matval[l], matval[k], etol)){
			can_iterate = FALSE;
			break;
		     }
		     if (bin_type){
			if (!in_conflict){
			   new_row_lb = rows[matind[l]].lb;
			   new_row_ub = rows[matind[l]].ub;
			   if (matval[l] > 0.0){
			      new_row_lb += matval[l];
			   } else {
			      new_row_ub += matval[l];
			   }
			   if (matval[k] > 0.0){
			      new_row_lb += matval[k];
			   } else {
			      new_row_ub += matval[k];
			   }

			   switch(sense[matind[l]]){
			    case 'E':
			       if (new_row_lb > rhs[matind[l]] + etol ||
				  new_row_ub < rhs[matind[l]] - etol){
				  in_conflict = TRUE;
			       }
			       break;
			    case 'L':
			       if (new_row_lb > rhs[matind[l]] + etol){
				  in_conflict = TRUE;
			       }
			       break;
			   }
			}
		     }
		  }
		  k++;
		  l++;
	       }
	       if ((k == matbeg[cl_ind + 1] && l == matbeg[cr_ind + 1])){
		  break;
	       }
	    }

	    if (!can_iterate){
	       /* so columns are not same, go back */
	       r_ind++;
	       continue;
	    }

	    /* so we have same columns, check if we can aggregate or
	       delete one of them */

	    /* first check if we have conflict in binary case */
	    delete_ind = -1;
	    if (bin_type && in_conflict){
	       if (obj_diff > 0.0){
		  delete_ind = cl_ind;
		  delete_val = 0.0;
	       } else {
		  delete_ind = cr_ind;
		  delete_val = 0.0;
	       }
	       fix_type = FIX_BINARY;
	    } else {
	       /* so not binary or are not in conflict */
	       /*check if we can aggregate first*/
	       if (dup_type == 0){
		  /* same obj, same col */
		  /* just merge it to the one on the left and
		     make the one on the right invisible
		     row bounds wont change in this case */

		  /* first merge */
		  /* just the bounds */
		  if (lb[cl_ind] > -INF){
		     if (lb[cr_ind] <= -INF){
			lb[cl_ind] = -INF;
		     } else {
			lb[cl_ind] += lb[cr_ind];
		     }
		  }
		  if (ub[cl_ind] < INF){
		     if (ub[cr_ind] >= INF){
			ub[cl_ind] = INF;
		     } else {
			ub[cl_ind] += ub[cr_ind];
		     }
		  }

		  if (type_l == 'B'){
		     cols[cl_ind].var_type = 'I';
		  }

		  stats->vars_aggregated++;
		  /* now cleanup*/
		  /* well just fix it to 0.0 so that we
		     will take care of other attr with this row*/

		  delete_ind = cr_ind;
		  delete_val = 0.0;
		  fix_type = FIX_AGGREGATE;
	       } else {
		  if (bin_type){
		     /* cant do anything with this pair*/
		     r_ind++;
		     continue;
		  }

		  /* now check dup_type and see what we can do */
		  fix_type = FIX_OTHER;
		  switch(dup_type){
		   case 1:
		      delete_ind = cl_ind;
		      delete_val = lb[cl_ind];
		      break;
		   case 2:
		      delete_ind = cr_ind;
		      delete_val = ub[cr_ind];
		      break;
		   case 3:
		      delete_ind = cr_ind;
		      delete_val = lb[cr_ind];
		      break;
		   case 4:
		      delete_ind = cl_ind;
		      delete_val = ub[cl_ind];
		      break;
		   case 5:
		   case 6:
		   case 10:
		      stats->col_unbound_ind = cr_ind;
		      return PREP_UNBOUNDED;
		   case 7:
		   case 8:
		   case 9:
		      stats->col_unbound_ind = cl_ind;
		      return PREP_UNBOUNDED;
		   default:
		      printf("error in prep_delete_duplicate_row_cols()\n");
		      return PREP_OTHER_ERROR;
		  }
	       }
	    }

	    if (delete_ind >= 0){
	       if (delete_ind == cl_ind){
		  l_ind++;
		  r_ind = l_ind + 1;
	       } else {
		  r_ind++;
	       }
	       col_orig_type[col_del_cnt] = cols[delete_ind].var_type;
	       cols[delete_ind].var_type = 'F';
	       col_del_ind[col_del_cnt] = delete_ind;
	       col_fix_val[col_del_cnt] = delete_val;
	       col_fix_type[col_del_cnt++] = fix_type;
	    }
	 }
      }

      /* ok, now fix each of these duplicate cols */
      for (i = 0; i < col_del_cnt; i++){
	 cols[col_del_ind[i]].var_type = col_orig_type[i];
	 termcode = prep_modified_cols_update_info(P, 1, &col_del_ind[i],
						   -1, dive_level,
						   col_fix_val[i],
						   col_fix_type[i],
						   TRUE, FALSE);
	 if (PREP_QUIT(termcode)){
	    return termcode;
	 }
      }
   }

#if 0
   mark_time = wall_clock(NULL);
   printf("Total duplicate cols Time: %f...\n\n",
	  mark_time - start_time);
#endif

   /* now same for rows */

   if (check_rows){
      //qsort_di(row_sum, r_loc, m);
      CoinSort_2(row_sum, row_sum+m, r_loc);
      last_lloc = last_rloc = 0;
      while(TRUE && check_rows){

	 if (last_lloc == m - 1){
	    break;
	 }

	 for (i = last_lloc; i < m - 1; i++){
	    if (prep_is_equal(row_sum[i], row_sum[i+1], etol)){
	       last_rloc = i+1;
	       if ( i < m - 2 ){
		  for (j = i+2; j < m; j++){
		     if (!prep_is_equal(row_sum[i], row_sum[j], etol)){
			last_rloc = j;
			break;
		     }
		  }
	       }
	       break;
	    }
	 }

	 if (i == m - 1){
	    break;
	 }

	 l_ind = i;
	 r_ind = l_ind + 1;
	 last_lloc = last_rloc;
	 while(l_ind < last_rloc - 1){

	    obj_ind = r_loc[l_ind];
	    row_ind = r_loc[r_ind];

	    if (r_ind == last_rloc ||
	       rows[obj_ind].is_redundant){
	       l_ind++;
	       r_ind = l_ind + 1;
	       continue;
	    }

	    if (rows[row_ind].is_redundant){
	       r_ind++;
	       continue;
	    }

	    obj_size = rows[obj_ind].size - rows[obj_ind].fixed_var_num;
	    row_size = rows[row_ind].size - rows[row_ind].fixed_var_num;

	    if (obj_size - row_size > 2 ||
	       obj_size - row_size < -2){
	       r_ind++;
	       continue;
	    }

	    delete_row = FALSE;

	    if (obj_size == 0){
	       delete_row = TRUE;
	       delete_row_ind = obj_ind;
	       l_ind++;
	       r_ind = l_ind + 1;
	    }

	    if (!delete_row){
	       if (row_size == 0){
		  delete_row = TRUE;
		  delete_row_ind = row_ind;
		  r_ind++;
	       }
	    }

	    if (!delete_row){

	       /* now check if rows are same */
	       diff_cnt = 0;
	       diff_ind = 0;
	       diff_obj_val = 0;
	       diff_row_val = 0;

	       for ( k = r_matbeg[obj_ind], l = r_matbeg[row_ind];;){
		  if (k < r_matbeg[obj_ind + 1] &&
		     (r_matind[k] < r_matind[l] ||
		      l >= r_matbeg[row_ind + 1])){

		     if (cols[r_matind[k]].var_type != 'F'){
			diff_ind = r_matind[k];
			diff_obj_val = r_matval[k];
			diff_cnt++;
		     }
		     k++;
		  } else if (l < r_matbeg[row_ind + 1] &&
			   (r_matind[k] > r_matind[l] ||
			    k >= r_matbeg[obj_ind+1])){
		     if (cols[r_matind[l]].var_type != 'F'){
			diff_ind = r_matind[l];
			diff_row_val = r_matval[l];
			diff_cnt++;
		     }
		     l++;
		  } else {
		     if (cols[r_matind[l]].var_type != 'F'){
			if (!prep_is_equal(r_matval[l], r_matval[k], etol)){
			   diff_ind = r_matind[k];
			   diff_obj_val = r_matval[k];
			   diff_row_val = r_matval[l];
			   diff_cnt++;
			}
		     }
		     k++;
		     l++;
		  }
		  if (diff_cnt > 1 ||
		     (k == r_matbeg[obj_ind + 1] &&
		      l == r_matbeg[row_ind + 1])){
		     break;
		  }
	       }

	       if (diff_cnt < 2 &&(diff_cnt == 0 ||
				  prep_is_equal(diff_obj_val -
						diff_row_val,
						0.0,
						etol))){
		  rhs_obj = rhs[obj_ind] - rows[obj_ind].fixed_lhs_offset;
		  rhs_row = rhs[row_ind] + rows[row_ind].fixed_lhs_offset;
		  delete_row = TRUE;
		  if (sense[obj_ind] == 'E'){
		     if (sense[row_ind] == 'E'){
			if (!prep_is_equal(rhs_obj, rhs_row, etol)){
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}
		     } else {
			if (rhs_row < rhs_obj - etol){
			   stats->row_infeas_ind = obj_ind;
			   return PREP_INFEAS;
			}
		     }
		     delete_row_ind = row_ind;
		     r_ind++;
		  } else {
		     if (sense[row_ind] == 'E'){
			if (rhs_row > rhs_obj + etol){
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}
			delete_row_ind = obj_ind;
			l_ind++;
			r_ind = l_ind + 1;
		     } else {
			if (rhs_row < rhs_obj - etol){
			   delete_row_ind = obj_ind;
			   l_ind++;
			   r_ind = l_ind + 1;
			} else {
			   delete_row_ind = row_ind;
			   r_ind++;
			}
		     }
		  }
	       } else if (diff_cnt == 1){
		  rhs_obj = rhs[obj_ind] - rows[obj_ind].fixed_lhs_offset;
		  rhs_row = rhs[row_ind] + rows[row_ind].fixed_lhs_offset;
		  diff_val = diff_obj_val - diff_row_val;
		  new_bound = (rhs_obj - rhs_row)/diff_val;

		  if (sense[obj_ind] == 'E'){
		     if (sense[row_ind] == 'E'){
			if (obj_size > row_size){
			   delete_row_ind = row_ind;
			   r_ind++;
			} else {
			   delete_row_ind = obj_ind;
			   l_ind++;
			   r_ind = l_ind + 1;
			}
			termcode = prep_modified_cols_update_info(P, 1,
								  &diff_ind,
								  -1,
								  dive_level,
								  new_bound,
								  FIX_OTHER,
								  TRUE, FALSE);
			if (PREP_QUIT(termcode)){
			   return termcode;
			}
			delete_row = TRUE;
		     } else {
			fix_type = FIX_NO_BOUND;
			if (obj_size > row_size){
			   if (diff_val > 0.0){
			      if (lb[diff_ind] < new_bound - etol){
				 fix_type = IMPROVE_LB;
			      }
			   } else {
			      if (ub[diff_ind] > new_bound + etol){
				 fix_type = IMPROVE_UB;
			      }
			   }
			} else {
			   if (diff_val > 0.0){
			      fix_type = IMPROVE_UB;
			   } else {
			      fix_type = IMPROVE_LB;
			   }
			}
			if (fix_type == IMPROVE_UB){
			   if (ub[diff_ind] < new_bound - etol){
			      fix_type = FIX_NO_BOUND;
			   }
			} else if (fix_type == IMPROVE_LB){
			   if (lb[diff_ind] > new_bound + etol){
			      fix_type = FIX_NO_BOUND;
			   }
			}
			if (fix_type != FIX_NO_BOUND){
			   termcode =
			      prep_modified_cols_update_info(P, 1,
							     &diff_ind,
							     -1,
							     dive_level,
							     new_bound,
							     fix_type,
							     TRUE, FALSE);
			   if (PREP_QUIT(termcode)){
			      return termcode;
			   }
			}
			r_ind++;
		     }
		  } else {
		     if (sense[obj_ind] != 'E'){
			fix_type = FIX_NO_BOUND;
			if (obj_size > row_size){
			   if (diff_val > 0.0){
			      if (ub[diff_ind] > new_bound + etol){
				 fix_type = IMPROVE_UB;
			      }
			   } else {
			      if (lb[diff_ind] < new_bound - etol){
				 fix_type = IMPROVE_LB;
			      }
			   }
			} else {
			   if (diff_val > 0.0){
			      fix_type = IMPROVE_LB;
			   } else {
			      fix_type = IMPROVE_UB;
			   }
			}

			if (fix_type == IMPROVE_UB){
			   if (ub[diff_ind] < new_bound - etol){
			      fix_type = FIX_NO_BOUND;
			   }
			} else if (fix_type == IMPROVE_LB){
			   if (lb[diff_ind] > new_bound + etol){
			      fix_type = FIX_NO_BOUND;
			   }
			}
			if (fix_type != FIX_NO_BOUND){
			   termcode =
			      prep_modified_cols_update_info(P, 1,
							     &diff_ind,
							     -1,
							     dive_level,
							     new_bound,
							     fix_type,
							     TRUE, FALSE);
			   if (PREP_QUIT(termcode)){
			      return termcode;
			   }
			}
		     }
		     r_ind++;
		  }
	       } else {
		  r_ind++;
	       }

	    }

	    if (delete_row){
	       stats->rows_deleted++;
	       if (verbosity >= 13){
		  prep_declare_redundant_row(rows[delete_row_ind],
					     delete_row_ind,
					     sense[delete_row_ind],
					     rhs[delete_row_ind]);
	       }
	       termcode = prep_deleted_row_update_info(mip, delete_row_ind);
	       if (PREP_QUIT(termcode)){
		  return termcode;
	       }
	    }
	 }
      }
   }

#if 0
   printf("Total duplicate rows Time: %f...\n\n",
	  wall_clock(NULL) - mark_time);
#endif

   /* free memory */

   if (check_cols){
      FREE(col_orig_type);
      FREE(col_del_ind);
      FREE(col_fix_type);
      FREE(col_fix_val);

      FREE(col_sum);
      FREE(row_factor);
      FREE(c_loc);
   }
   if (check_rows){
      FREE(row_sum);
      FREE(col_factor);
      FREE(r_loc);
   }

   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
/* check if we can bound an unbounded variable
   only used when this row has only 1 lb_inf_var_num or 1 ub_inf_var_num */

int prep_force_row_bounds(PREPdesc *P, int row_ind, int col_ind, int a_loc)
{

   int termcode;
   MIPdesc *mip = P->mip;
   ROWinfo *rows = mip->mip_inf->rows;

   double *rhs = mip->rhs;
   double *ub = mip->ub;
   double *lb = mip->lb;
   char sense = mip->sense[row_ind];

   double new_bound;

   char row_lb_improved = FALSE;
   char row_ub_improved = FALSE;
   int fix_type = FIX_NO_BOUND;

   double etol = P->params.etol;
   int dive_level = 0;
   if (rows[row_ind].lb <= -INF && rows[row_ind].ub >= INF){
      return PREP_UNMODIFIED;
   }

   double a_val = mip->matval[a_loc];

   /* debug check */

   if (sense != 'E'){
      if (!((a_val > 0.0 && ub[col_ind] >= INF) ||
	   (a_val < 0.0 && lb[col_ind] <= -INF))){
	 printf("error in prep_force_row_bounds()\n");
	 return PREP_OTHER_ERROR;
      }
   } else {
      if (!((a_val > 0.0 && ub[col_ind] >= INF) ||
	   (a_val < 0.0 && lb[col_ind] <= -INF)||
	   (a_val < 0.0 && ub[col_ind] >= INF) ||
	   (a_val > 0.0 && lb[col_ind] <= -INF))){
	 printf("error -1 in prep_force_row_bounds()\n");
	 return PREP_OTHER_ERROR;
      }
   }

   /* note that, if this row has only one unbounded variable,
      we can try to bound this variable using the rhs of this row
      and the bounds of the other present variables
   */

   if (rows[row_ind].ub_inf_var_num <= 1){
      if (a_val > etol && ub[col_ind] >= INF){
	 if (rows[row_ind].lb > -INF){
	    new_bound = (double)((rhs[row_ind] - rows[row_ind].lb +
				  lb[col_ind]*a_val)/a_val);

	    fix_type = IMPROVE_UB;
	    row_ub_improved = TRUE;
	 }
      } else if (a_val < -etol && lb[col_ind] <= -INF){
	 if (rows[row_ind].lb > -INF){

	    new_bound = (double)((rhs[row_ind] - rows[row_ind].lb +
				  ub[col_ind]*a_val)/a_val);

	    fix_type = IMPROVE_LB;
	    row_ub_improved = TRUE;
	 }
      }
   } else if (sense == 'E'){
      if (a_val > etol && lb[col_ind] <= -INF){
	 if (rows[row_ind].ub < INF){

	    new_bound = (double)((rhs[row_ind] - rows[row_ind].ub +
				  ub[col_ind]*a_val)/a_val);

	    fix_type = IMPROVE_LB;
	    row_lb_improved = TRUE;
	 }
      } else if (a_val < -etol && ub[col_ind] >= INF){
	 if (rows[row_ind].ub < INF){

	    new_bound = (double)((rhs[row_ind] - rows[row_ind].ub +
				  lb[col_ind]*a_val)/a_val);

	    fix_type = IMPROVE_UB;
	    row_lb_improved = TRUE;
	 }
      }
   }

   if (row_lb_improved || row_ub_improved){
      //ub[col_ind] = old_ub;
      //lb[col_ind] = old_lb;

      /* we have already obtained this rows bound, however we need to
	 apply this bound to other rows */

      termcode = prep_modified_cols_update_info(P, 1, &col_ind, row_ind,
					       dive_level,
					       new_bound,
					       fix_type, TRUE, FALSE);

      if (PREP_QUIT(termcode)){
	 return termcode;
      }
      return PREP_MODIFIED;
   }

   return PREP_UNMODIFIED;

}

/*===========================================================================*/
/*===========================================================================*/
/* check if we can improve the variable bounds or the rhs or the coefficient
   itself
*/

int prep_improve_variable(PREPdesc *P, int col_ind, int row_ind, int a_loc,
			  int dive_level, char check_improve, char impl_mode,
			  char use_sr_bounds, double sr_ub, double sr_lb,
			  int use_mip)

{

   int i, fix_type = FIX_NO_BOUND, termcode = PREP_UNMODIFIED;

   MIPdesc *mip = P->mip;

   COLinfo *cols = mip->mip_inf->cols;
   ROWinfo *rows = mip->mip_inf->rows;
   double *maj_matval;

   if (use_mip == MAT_COL_ORDERED){
      maj_matval = mip->matval;
   } else {
      maj_matval = mip->row_matval;
   }

   double *ub = mip->ub;
   double *lb = mip->lb;
   double *obj = mip->obj;

   char is_int = mip->is_int[col_ind];
   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double a_val = maj_matval[a_loc];

   char fix_to_lb, fix_to_ub, improve_coef;
   char col_lb_unbounded, col_ub_unbounded;

   double new_bound, new_ub, new_lb;

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   double coeff_etol = 1e-15;
   prep_stats *stats = &(P->stats);

   /* first we check if we can fix this column under the condition that
      it has only one nonzero element */

   if (cols[col_ind].var_type != 'U' &&
      cols[col_ind].var_type != 'L'){
      if (cols[col_ind].col_size <= 1){
	 if (obj[col_ind] >= 0.0 && ( (sense == 'G' && a_val < -etol) ||
				     (sense == 'L' && a_val > etol) ||
				     (cols[col_ind].col_size == 0) ) ){
	   if (lb[col_ind] <= -INF){
	      if(obj[col_ind] > coeff_etol){
		 stats->col_unbound_ind = col_ind;
	         return PREP_UNBOUNDED;
	      }
	   } else {
	      new_bound = lb[col_ind];
	       if (cols[col_ind].var_type == 'B'){
		  fix_type = FIX_BINARY;
	       } else {
		  fix_type = FIX_OTHER;
	       }
 	   }
	 } else if (obj[col_ind] <= 0.0 && ( (sense == 'G' && a_val > etol) ||
					   (sense == 'L' && a_val < -etol) ||
					   (cols[col_ind].col_size == 0)) ){
	   if (ub[col_ind] >= INF){
	      if(obj[col_ind] < -coeff_etol){
	         stats->col_unbound_ind = col_ind;
	         return PREP_UNBOUNDED;
	      }
	    } else{
	       new_bound = ub[col_ind];
	       if (cols[col_ind].var_type == 'B'){
	  	  fix_type = FIX_BINARY;
	       } else {
		  fix_type = FIX_OTHER;
	       }
	    }
	 }
	 if (fix_type != FIX_NO_BOUND){
 	    termcode = PREP_MODIFIED;
	 }
      }
   }

   if (termcode == PREP_UNMODIFIED){
      /* continue */
      /* firt case: it is a binary variable */
      if (cols[col_ind].var_type == 'B'){

	 /* see if we can fix this variable */
	 /* or improve the coefficients */

	 fix_to_lb = FALSE;
	 fix_to_ub = FALSE;
	 improve_coef = FALSE;

	 /* if the coefficient is positive: */
	 if (a_val > etol){
	    switch(sense){
	     case 'L':
		/* try to fix */
		if (rows[row_ind].lb > -INF){
		   if (use_sr_bounds){
		      new_lb = sr_lb;
		   } else {
		      new_lb = rows[row_ind].lb + a_val;
		   }
		   if (new_lb > rhs + etol)
		      fix_to_lb = TRUE;
		}
		/* try to improve */
		if (!fix_to_lb && check_improve && !impl_mode){
		   if (rows[row_ind].ub < INF){
		      /* this row might have been redundant all the
			 way up here */
		      if (use_sr_bounds){
			 new_ub = sr_ub;
			 if (sr_ub < rhs - etol){
			    new_ub = sr_ub - rhs;
			    maj_matval[a_loc] -= new_ub;
			    mip->rhs[row_ind] -= new_ub;

			    if (prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
			    }
			    rows[row_ind].ub +=
			       (maj_matval[a_loc] - a_val) *
			       ub[col_ind];

			    improve_coef = TRUE;
			 }
		      } else {
			 new_ub = rows[row_ind].ub - a_val;
			 if (new_ub < rhs - etol){

			    /* update coefs */
			    maj_matval[a_loc] = rows[row_ind].ub - rhs;
			    mip->rhs[row_ind] = new_ub;

			    /* debug */
			    if (maj_matval[a_loc] < -etol){
			       printf("error -0 in prep_improve_variable()\n");
			       return PREP_OTHER_ERROR;
			    }

			    if (prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
			    }

			    /* update bounds */
			    rows[row_ind].ub +=
			       (maj_matval[a_loc] - a_val) *
			       ub[col_ind];

			    improve_coef = TRUE;
			 }
		      }
		   }
		}
		break;
	     case 'G':
		/* debug */
		/* cant have 'G' row */
		printf("error -2 in prep_improve_variable()\n");
		return PREP_OTHER_ERROR;
	     case 'E':
		if (rows[row_ind].lb > -INF){
		   if (use_sr_bounds){
		      new_lb = sr_lb;
		   } else {
		      new_lb = rows[row_ind].lb + a_val;
		   }
		   if (new_lb > rhs){
		      fix_to_lb = TRUE;
		   }
		}
		if (rows[row_ind].ub < INF){
		   if (use_sr_bounds){
		      new_ub = sr_ub;
		   } else {
		      new_ub = rows[row_ind].ub - a_val;
		   }
		   if (new_ub < rhs){
		      fix_to_ub = TRUE;
		   }
		}
		if (fix_to_lb && fix_to_ub){
		   stats->col_infeas_ind = col_ind;
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
		break;
	    }
	 } else if (a_val < -etol){
	    /* if the coefficient is negative: */
	    switch(sense){
	     case 'L':
		if (rows[row_ind].lb > -INF){
		   if (use_sr_bounds){
		      new_lb = sr_lb;
		   } else {
		      new_lb = rows[row_ind].lb - a_val;
		   }
		   if (new_lb > rhs)
		      fix_to_ub = TRUE;
		}
		if (!fix_to_ub && check_improve && !impl_mode){
		   if (rows[row_ind].ub < INF){
		      if (use_sr_bounds){
			 new_ub = sr_ub;
			 if (sr_ub < rhs - etol){
			    new_ub = sr_ub - rhs;
			    maj_matval[a_loc] -= new_ub;

			    if (prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
			    }
			    rows[row_ind].lb +=
			       (maj_matval[a_loc] - a_val) *
			       ub[col_ind];

			    improve_coef = TRUE;
			 }
		      } else {
			 new_ub = rows[row_ind].ub + a_val;
			 if (new_ub < rhs - etol){
			    /* update coef*/
			    maj_matval[a_loc] -= new_ub - rhs;

			    /* debug */
			    if (maj_matval[a_loc] > etol){
			       printf("error -3 in prep_improve_variable()\n");
			       return PREP_OTHER_ERROR;
			    }

			    if (prep_is_equal(maj_matval[a_loc], 0.0, etol)){
			       maj_matval[a_loc] = 0.0;
			    }

			    /* update bounds */
			    if (rows[row_ind].lb > -INF){
			       rows[row_ind].lb +=
				  (maj_matval[a_loc] - a_val) *
				  ub[col_ind];
			    }

			    improve_coef = TRUE;
			 }
		      }
		   }
		}
		break;
	     case 'G':
		/* debug */
		/* cant have 'G' row */
		printf("error -5 in prep_improve_variable()\n");
		return PREP_OTHER_ERROR;
		break;
	     case 'E':
		if (rows[row_ind].lb > -INF){
		   if (use_sr_bounds){
		      new_lb = sr_lb;
		   } else {
		      new_lb = rows[row_ind].lb - a_val;
		   }
		   if (new_lb > rhs){
		      fix_to_ub = TRUE;
		   }
		}
		if (rows[row_ind].ub < INF){
		   if (use_sr_bounds){
		      new_ub = sr_ub;
		   } else {
		      new_ub = rows[row_ind].ub + a_val;
		   }
		   if (new_ub < rhs){
		      fix_to_lb = TRUE;
		   }
		}
		if (fix_to_lb && fix_to_ub){
		   stats->col_infeas_ind = col_ind;
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
		break;
	    }
	 }

	 if (fix_to_lb || fix_to_ub){
	    /* we have managed to fix it to either its lower or upper bound */
	    if (fix_to_lb){
	       new_bound = 0.0;
	    } else {
	       new_bound = 1.0;
	    }
	    fix_type = FIX_BINARY;
	    termcode = PREP_MODIFIED;
	 } else if (improve_coef){
	    /* so we can improve a_val and rhs */

	    /* need to update row bounds again here */
	    /* debug -fixme */
	    /* i really dont like this brute forcing here, try to fix it*/

	    if (use_mip == MAT_COL_ORDERED){
	       for (i = mip->row_matbeg[row_ind]; i <
		      mip->row_matbeg[row_ind + 1]; i++){
		  if (mip->row_matind[i] == col_ind){
		     mip->row_matval[i] = maj_matval[a_loc];
		     break;
		  }
	       }
	       /* debug */
	       if (i == mip->row_matbeg[row_ind + 1]){
		  printf("error -1 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;
	       }
	    } else {
	       for (i = mip->matbeg[col_ind]; i <
		      mip->matbeg[col_ind + 1]; i++){
		  if (mip->matind[i] == row_ind){
		     mip->matval[i] = maj_matval[a_loc];
		     break;
		  }
	       }
	       /* debug */
	       if (i == mip->matbeg[col_ind + 1]){
		  printf("error -6 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;
	       }
	    }

	    if (verbosity >=14){
	       if (mip->colname){
		  prep_declare_coef_change(row_ind, col_ind,
					   mip->colname[col_ind],
					   maj_matval[a_loc],
					   mip->rhs[row_ind]);
	       } else {
		  prep_declare_coef_change(row_ind, col_ind,
					   NULL,
					   maj_matval[a_loc],
					   mip->rhs[row_ind]);
	       }
	    }
	    if (!(stats->nz_coeff_changed[a_loc])){
	       stats->nz_coeff_changed[a_loc] = TRUE;
	       stats->coeffs_changed++;
	    }

	    /* since row bound(s) has been updated, check for redundancy */
	    /* usually this shouldnt be here, but we have updated row bounds
	       here */

	    termcode = prep_check_redundancy(P, row_ind, FALSE, 0.0, 0.0,
					     impl_mode, dive_level);

	    if (PREP_QUIT(termcode)){
	       return termcode;
	    } else if (cols[col_ind].var_type == 'F'){
	       return PREP_MODIFIED;
	    }

	    fix_type = IMPROVE_COEF;
	    termcode = PREP_MODIFIED;
	 } else if (FALSE && !impl_mode && //disabled --
		  ((a_val > etol && !P->ulist_checked[col_ind]) ||
		   (a_val < -etol && !P->llist_checked[col_ind]))){

	    /* for now - just among the binary variables*/
	    /* so cant fix it, cant improve it and binary
	       try logical fixing - with dive-level = 1 */

	    /* fist copy initial info
	       fixme! - work on this!
	    */

	    /* do once for each variable */
	    int impl_dive_level = 2;
	    memcpy(P->impl_rows, rows, sizeof(ROWinfo)*mip->m);
	    memcpy(P->impl_cols, cols, sizeof(COLinfo)*mip->n);
	    memcpy(P->impl_ub, ub, DSIZE*mip->n);
	    memcpy(P->impl_lb, lb, DSIZE*mip->n);
	    P->impl_stats = P->stats;
	    if (a_val > etol){

	       //if (!cols[col_ind].ulist){
		  // cols[col_ind].ulist = (IMPlist *)calloc(sizeof(IMPlist),1);
	       //}

	       P->list = cols[col_ind].ulist;
	       P->ulist_checked[col_ind] = TRUE;

	       /* fix it to 1.0 and see if that causes any infeasibility
		  otherwise get the impllist and continue*/
	       /* get the implication list */

	       /* fix this column, update row bounds of this column
		  check for redundancy, */

	       termcode = prep_modified_cols_update_info(P, 1, &col_ind,
							row_ind,
							 impl_dive_level,
							 1.0,
							 IMPROVE_LB,
							 TRUE, TRUE);
	       if (termcode == PREP_INFEAS){
		  printf("infeasibility detected!\n");
		  /*then this column is fixable to its lower bound! */
		  new_bound = 0.0;
		  fix_type = FIX_BINARY;
		  termcode = PREP_MODIFIED;
	       } else {
		  termcode = PREP_UNMODIFIED;
	       }
	    } else if (a_val < etol){

	       //if (!cols[col_ind].llist){
		  //  cols[col_ind].llist = (IMPlist *)calloc(sizeof(IMPlist),1);
	       //}

	       P->list = cols[col_ind].llist;
	       P->llist_checked[col_ind] = TRUE;

	       /* fix it to 1.0 and see if that causes any infeasibility
		  otherwise get the impllist and continue*/

	       termcode = prep_modified_cols_update_info(P, 1, &col_ind,
							 row_ind,
							 impl_dive_level,
							 0.0, IMPROVE_UB,
							 TRUE, TRUE);

	       if (termcode == PREP_INFEAS){
		  printf("infeasibility detected!\n");
		  /*then this column is fixable to its lower bound! */
		  new_bound = 1.0;
		  fix_type = FIX_BINARY;
		  termcode = PREP_MODIFIED;
	       } else {
		  termcode = PREP_UNMODIFIED;
	       }
	    }

	    /* now get back */
	    memcpy(rows, P->impl_rows,sizeof(ROWinfo)*mip->m);
	    memcpy(cols, P->impl_cols, sizeof(COLinfo)*mip->n);
	    memcpy(ub, P->impl_ub, DSIZE*mip->n);
	    memcpy(lb, P->impl_lb, DSIZE*mip->n);
	    P->stats = P->impl_stats;
	 }
      } else if (cols[col_ind].var_type == 'U' ||
	       cols[col_ind].var_type == 'L'){
	 if (cols[col_ind].var_type == 'U'){
	    new_bound = ub[col_ind];
	    if (is_int){
	       new_bound = prep_rnd_integral(new_bound, etol, RND_FLOOR);
	    }
	 } else {
	    new_bound = lb[col_ind];
	    if (is_int){
	       new_bound = prep_rnd_integral(new_bound, etol, RND_CEIL);
	    }
	 }

	 fix_type = FIX_OTHER;
	 termcode = PREP_MODIFIED;

      } else {

	 /* not binary, not fixable etc. */
	 /* now try bounds improvement */

	 col_lb_unbounded = FALSE;
	 col_ub_unbounded = FALSE;

	 if (a_val > etol){
	    if (lb[col_ind] <= -INF){

	       /* debug */
	       if (rows[row_ind].lb > -INF){
		  printf("error -7 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;
	       }

	       col_lb_unbounded = TRUE;

	       /* can we fix it? */
	       /* is fixable if sense = 'E' and ub < INF*/

	       if (sense == 'E' && rows[row_ind].ub < INF){
		  new_bound = (double)((rhs - rows[row_ind].ub +
					a_val*ub[col_ind])/a_val);
		  if (cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_CEIL);
		  }
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind,
							   row_ind, dive_level,
							   new_bound,
							   IMPROVE_LB, TRUE,
							    impl_mode);
		  if (PREP_QUIT(termcode)){
		     return termcode;
		  } else if (rows[row_ind].is_redundant){
		     return PREP_MODIFIED;
		  }
		  termcode = PREP_UNMODIFIED;
		  col_lb_unbounded = FALSE;
	       }
	    }
	    if (!col_lb_unbounded){
	       if (rows[row_ind].lb > -INF){
		  new_bound = (double)((rhs - rows[row_ind].lb +
					a_val*lb[col_ind])/a_val);
		  if (cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_FLOOR);
		  }

		  if (new_bound < ub[col_ind]){
		     termcode = PREP_MODIFIED;
		     fix_type = IMPROVE_UB;
		  }
	       }
	    }
	 } else if (a_val < -etol){
	    if (ub[col_ind] >= INF){

	       /* debug */
	       if (rows[row_ind].lb > -INF){
		  printf("error -2 in prep_improve_variable()\n");
		  return PREP_OTHER_ERROR;
	       }

	       col_ub_unbounded = TRUE;
	       /* can we fix it? */
	       /* is fixable if sense = 'E' and ub < INF*/
	       if (sense == 'E' && rows[row_ind].ub < INF){
		  new_bound = (double)((rhs - rows[row_ind].ub +
					a_val*lb[col_ind])/a_val);
		  if (cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_FLOOR);
		  }
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind,
							   row_ind, dive_level,
							   new_bound,
							   IMPROVE_UB, TRUE,
							    impl_mode);
		  if (PREP_QUIT(termcode)){
		     return termcode;
		  } else {
		     if (rows[row_ind].is_redundant){
			return PREP_MODIFIED;
		     }
		  }
		  termcode = PREP_UNMODIFIED;
		  col_ub_unbounded = FALSE;
	       }
	    }
	    if (!col_ub_unbounded){
	       if (rows[row_ind].lb > -INF){
		  new_bound = (double)((rhs - rows[row_ind].lb +
					a_val*ub[col_ind])/a_val);
		  if (cols[col_ind].var_type != 'C'){
		     new_bound = prep_rnd_integral(new_bound, etol, RND_CEIL);
		  }

		  if (new_bound > lb[col_ind]){
		     termcode = PREP_MODIFIED;
		     fix_type = IMPROVE_LB;
		  }
	       }
	    }
	 }
      }
   }

   /* now check if we need to update row bounds */
   if (termcode == PREP_MODIFIED && fix_type != IMPROVE_COEF){

      /* set col type to 'T', set it to F after you have visited all
	 other rows? */
      /* have col.fix_row_ind to mark on which row you have fixed it? */
      /* isnt worth it, update row bounds here */
      //cols[col_ind].fix_row_ind = row_ind;
      /* this col might have been improved above */
      if (cols[col_ind].var_type != 'F'){

	 termcode = prep_modified_cols_update_info(P, 1, &col_ind, row_ind,
						  dive_level,
						  new_bound,
						  fix_type, TRUE, impl_mode);
	 if (PREP_QUIT(termcode)){
	    return termcode;
	 } else {
	    return PREP_MODIFIED;
	 }
      }
   }

   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
/* Under Development */
int prep_add_to_impl_list(IMPlist *list, int ind, int fix_type,
			  double val){

   if (!list){
      printf("error in prep_add_to_impl_list\n");
      exit(0);
   }

   IMPvar * var = (IMPvar *)calloc(sizeof(IMPvar),1);

   var->ind = ind;
   var->fix_type = fix_type;
   var->val = val;

   if (!list->head){
      list->head = list->tail = var;
   } else {
      list->tail->right = var;
      list->tail = var;
   }

   list->size++;
   return 0;

}
/*===========================================================================*/
/*===========================================================================*/
/* Imply the modified column to the presolve description */

int prep_modified_cols_update_info(PREPdesc *P, int col_cnt, int *col_start,
				   int row_ind, int dive_level,
				   double fixed_bound,  int intl_fix_type,
				  char check_redundancy, char impl_mode)
{
   /* fix_type
      0 FIX_NO_BOUND
      1 FIX_BINARY
      2 FIX_FIXABLE
      3 FIX_OTHER
      4 IMPROVE_UB
      5 IMPROVE_LB
      6 IMPROVE_COEF
      7 FIX_ALL_LB
      8 FIX_ALL_UB
   */

   int termcode = PREP_UNMODIFIED, i, j, k, r_ind, end;
   int col_ind, a_loc_ref = 0;
   MIPdesc *mip = P->mip;
   int *matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;
   char *is_int = mip->is_int;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;
   //double * obj = mip->obj;

   double old_ub;
   double old_lb;

   double *ub = mip->ub;
   double *lb = mip->lb;

   double a_val;
   char get_row_ubounds;
   char get_row_lbounds;

   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   prep_stats *stats = &(P->stats);

   int fix_type = 0, row_cnt = 0;
   char row_updated;
   char *is_row_updated = NULL;
   int *row_updated_ind = NULL;
   IMPlist *imp_list;
   IMPvar *imp_var;

    /*first, if in impl mode, add these variables to current impl_list */
   /* do this above */
   int can_iterate = TRUE;

   double mark_time = wall_clock(NULL);
   //if (intl_fix_type != FIX_AGGREGATE){
   is_row_updated = (char *)calloc(CSIZE,mip->m);
   row_updated_ind = (int *)malloc(ISIZE*mip->m);
   //}

   P->alloc2_time += wall_clock(NULL) - mark_time;

   if (intl_fix_type == FIX_ROW_LB ||
      intl_fix_type == FIX_ROW_UB){
      a_loc_ref = r_matbeg[row_ind];
   }

   mark_time = wall_clock(NULL);

   for (j = 0; j < col_cnt; j++){

      col_ind = col_start[j];

      if (cols[col_ind].var_type == 'F'){
	 if (intl_fix_type != FIX_ROW_LB &&
	    intl_fix_type != FIX_ROW_UB){
	    if (!prep_is_equal(ub[col_ind], fixed_bound, etol)){
	       if (!impl_mode){
		  stats->col_infeas_ind = col_ind;
	       }
	       termcode = PREP_INFEAS;
	       can_iterate = FALSE;
	       break;
	    }
	 }
	 continue;
      }

      fix_type = intl_fix_type;

      old_ub = ub[col_ind];
      old_lb = lb[col_ind];

      if (fix_type == FIX_ROW_LB || fix_type == FIX_ROW_UB){
	 a_val = r_matval[a_loc_ref + j];
	 if (fix_type == FIX_ROW_LB){
	    if (a_val > etol){
	       fixed_bound = ub[col_ind] = lb[col_ind];
	    } else if (a_val < -etol){
	       fixed_bound = lb[col_ind] = ub[col_ind];
	    } else {
	       continue;
	    }
	    if (cols[col_ind].var_type == 'B'){
	       fix_type = FIX_BINARY;
	    } else {
	       fix_type = FIX_OTHER;
	    }
	 } else {
	    if (a_val > etol){
	       fixed_bound = lb[col_ind] = ub[col_ind];
	    } else if (a_val < -etol){
	       fixed_bound = ub[col_ind] = lb[col_ind];
	    } else {
	       continue;
	    }
	    if (cols[col_ind].var_type == 'B'){
	       fix_type = FIX_BINARY;
	    } else {
	       fix_type = FIX_OTHER;
	    }
	 }
      } else {
	 if (fix_type != IMPROVE_LB){
	    if (fix_type != FIX_AGGREGATE){
	      if (ub[col_ind] < fixed_bound - etol){
		if (!impl_mode){
		  stats->col_infeas_ind = col_ind;
		  stats->row_infeas_ind = row_ind;
		}
		termcode = PREP_INFEAS;
		can_iterate = FALSE;
		break;
	      }
	    }
	    ub[col_ind] = fixed_bound;
	    //(stats->bounds_tightened)++;
	 }
	 if (fix_type != IMPROVE_UB){
	    if (fix_type != FIX_AGGREGATE){
	      if (lb[col_ind] > fixed_bound + etol){
		if (!impl_mode){
		  stats->col_infeas_ind = col_ind;
		  stats->row_infeas_ind = row_ind;
		}
		termcode = PREP_INFEAS;
		can_iterate = FALSE;
		break;
	      }
	    }
	    lb[col_ind] = fixed_bound;
	 }
      }

      if (verbosity >= 12){
	 if (fix_type == FIX_AGGREGATE){
	    if (mip->colname){
	       printf("var %s [%i] is aggregated: \n", mip->colname[col_ind],
		      col_ind);
	    } else {
	       printf("var [%i] is aggregated: \n", col_ind);
	    }
	 } else {
	    if (mip->colname){
	       printf("var %s [%i] bounds are improved: ",
		      mip->colname[col_ind], col_ind);
	    } else {
	       printf("var [%i] bounds are improved: ", col_ind);
	    }
	    if (lb[col_ind] > -INF){
	       printf("\t lb:%f", lb[col_ind]);
	    }
	    if (ub[col_ind] < INF){
	       printf("\t ub:%f ", ub[col_ind]);
	    }
	    printf("\n");
	 }
      }

      if (fix_type != FIX_BINARY){
	 if (fix_type != FIX_AGGREGATE){
	    if (prep_is_equal(ub[col_ind], lb[col_ind], etol)){
	       if (cols[col_ind].var_type == 'B'){
		  fix_type = FIX_BINARY;
	       } else {
		  fix_type = FIX_OTHER;
	       }
	       cols[col_ind].var_type = 'F';
	    }
	 } else {
	    cols[col_ind].var_type = 'F';
	 }
      } else {
	 cols[col_ind].var_type = 'F';
      }

      /* now add to impl list if in impl_mode */

      if (fix_type != FIX_AGGREGATE && FALSE){ //disabled now
	 if (impl_mode && P->impl_col_ind != col_ind){
	    if (P->list->size < P->impl_limit){
	       prep_add_to_impl_list(P->list, col_ind, fix_type,
				     fixed_bound);
	    }
	 }
      }

      if (cols[col_ind].var_type == 'F'){
	 if (fix_type != FIX_AGGREGATE){

	 /* first see if you can fix any other variables from the
	    impl list of this variable */
	    if (fix_type == FIX_BINARY && FALSE){ //disabled now
	       if (lb[col_ind] >= 1.0 - etol){
		  imp_list = cols[col_ind].ulist;
	       } else {
		  imp_list = cols[col_ind].llist;
	       }
	       if (imp_list){
		  if (imp_list->size > 0){
		     for (imp_var = imp_list->head; imp_var != 0;
			 imp_var = imp_var->right){
			termcode = prep_modified_cols_update_info(P, 1,
								  &imp_var->ind,
								  -1, 0,
								  imp_var->val,
								  imp_var->fix_type,//FIX_BINARY,
								  FALSE, FALSE);
			if (PREP_QUIT(termcode)){
			   can_iterate = FALSE;
			   break;
			}
		     }

		     if (!can_iterate) break;
		  }
	       }
	    }

	    if (verbosity >= 12){
	       if (mip->colname){
		  prep_declare_fixed_var(col_ind, mip->colname[col_ind],
					 ub[col_ind]);
	       } else {
		  prep_declare_fixed_var(col_ind, NULL,
					 ub[col_ind]);
	       }
	    }

	    if (!impl_mode){
	       mip->mip_inf->sum_obj_offset += mip->obj[col_ind]*ub[col_ind];
	    }
	 }

	 (stats->vars_fixed)++;

      } else {
	 (stats->bounds_tightened)++;
      }

      if (cols[col_ind].col_size == 0 ){
	 continue;
      } else if (cols[col_ind].col_size < 0){
	 printf("error -00 in prep_fixed_col_update_info()\n");
	 termcode = PREP_OTHER_ERROR;
	 can_iterate = FALSE;
	 break;
      }

      end = matbeg[col_ind + 1];

      for (i = matbeg[col_ind]; i < end; i++){

	 if (!(rows[matind[i]].is_redundant)){

	    a_val = matval[i];
	    get_row_ubounds = FALSE;
	    get_row_lbounds = FALSE;
	    r_ind = matind[i];
	    row_updated = FALSE;
	    if (fix_type != IMPROVE_UB && fix_type != IMPROVE_LB){
	       rows[r_ind].fixed_var_num++;
	       if (!is_int[col_ind]){
		  (rows[r_ind].cont_var_num)--;
	       }

	       if (!prep_is_integral(a_val, etol)){
		  (rows[r_ind].frac_coef_num)--;
	       }

	       if (fix_type == FIX_BINARY){
		  rows[r_ind].bin_var_num--;
	       }
	       /* debug */
	       if (rows[r_ind].bin_var_num < 0 ||
		  rows[r_ind].fixable_var_num < 0){
		  printf("error -0 in prep_fixed_col_update_info()\n");
		  termcode = PREP_OTHER_ERROR;
		  break;
	       }
	       if (fix_type != FIX_AGGREGATE){
		  rows[r_ind].fixed_lhs_offset += a_val * fixed_bound;
	       }
	    }

	    if (old_ub >= INF && fix_type != IMPROVE_LB){
	       if (a_val > etol){
		  rows[r_ind].ub_inf_var_num--;
		  if (rows[r_ind].ub_inf_var_num == 0){
		     get_row_ubounds = TRUE;
		  }
	       } else if (a_val < -etol){
		  rows[r_ind].lb_inf_var_num--;
		  if (rows[r_ind].lb_inf_var_num == 0){
		     get_row_ubounds = TRUE;
		  }
	       }
	    }

	    if (old_lb <= -INF && fix_type != IMPROVE_UB){
	       if (a_val > etol){
		  rows[r_ind].lb_inf_var_num--;
		  if (rows[r_ind].lb_inf_var_num == 0){
		     get_row_lbounds = TRUE;
		  }
	       } else if (a_val < -etol){
		  rows[r_ind].ub_inf_var_num--;
		  if (rows[r_ind].lb_inf_var_num == 0){
		     get_row_lbounds = TRUE;
		  }
	       }
	    }

	    if (fix_type != FIX_AGGREGATE){
	       if (a_val > etol){
		  if (fix_type != IMPROVE_LB){
		     if (rows[r_ind].ub < INF){
			/* debug */
			if (old_ub >= INF){
			   printf("error -1 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if (fixed_bound != old_ub){
			   rows[r_ind].ub += a_val*(fixed_bound - old_ub);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }
		  if (fix_type != IMPROVE_UB && fix_type != FIX_AGGREGATE){
		     if (rows[r_ind].lb > -INF){
			/* debug */
			if (old_lb <= -INF){
			   printf("error -2 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if (fixed_bound != old_lb){
			   rows[r_ind].lb += a_val*(fixed_bound - old_lb);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }
	       } else if (a_val < -etol){
		  if (fix_type != IMPROVE_UB){
		     if (rows[r_ind].ub < INF){
			/* debug */
			if (old_lb <= -INF){
			   printf("error -3 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if (fixed_bound != old_lb){
			   rows[r_ind].ub += a_val*(fixed_bound - old_lb);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }

		  if (fix_type != IMPROVE_LB){
		     if (rows[r_ind].lb > -INF){
			/* debug */
			if (old_ub >= INF){
			   printf("error -4 in prep_fixed_col_update_info()\n");
			   termcode = PREP_OTHER_ERROR;
			   break;
			}
			if (fixed_bound != old_ub){
			   rows[r_ind].lb += a_val*(fixed_bound - old_ub);
			   rows[r_ind].is_updated = TRUE;
			   row_updated = TRUE;
			}
		     }
		  }
	       }
	    }

	    /* debug */
	    if (rows[r_ind].fixed_var_num +
	       rows[r_ind].fixable_var_num > rows[r_ind].size){
	       printf("error in fixing vars 2, prep_fix_variable()\n");
	       termcode = PREP_OTHER_ERROR;
	       break;
	    }

	    if (get_row_lbounds || get_row_ubounds){
	       rows[r_ind].is_updated = TRUE;
	       row_updated = TRUE;
	       prep_get_row_bounds(mip, r_ind, etol);
	    }
	    if (row_updated){
	       if (!is_row_updated[r_ind]){
		  is_row_updated[r_ind] = TRUE;
		  row_updated_ind[row_cnt++] = r_ind;
	       }
	    }
	 }
      }
   }
   if (impl_mode){
      P->impl_cols_time += wall_clock(NULL) - mark_time;
   }

   /* if row_updated_cnt > 0 also just rows updated? this is inefficient*/

   mark_time = wall_clock(NULL);

   if (!PREP_QUIT(termcode) && check_redundancy &&
      fix_type != IMPROVE_LB &&
      fix_type != IMPROVE_UB && fix_type){

      for (i = 0; i < row_cnt; i++){
	 r_ind = row_updated_ind[i];
	 if (rows[r_ind].is_updated && !rows[r_ind].is_redundant){
	    //if (impl_mode){
	    //  printf("processing row %i\n", r_ind);
	    //}

	    termcode = prep_check_redundancy(P, r_ind, FALSE,
					     0.0, 0.0, impl_mode,
					     dive_level);
	    if (PREP_QUIT(termcode)){
	       break;
	    }

	    rows[r_ind].is_updated = FALSE;

	    /* now do we want to dive on variables of the rows those share
	       a comman variable
	       with these fixed column(s)? */
	    if (dive_level > 0){
	       for (k = r_matbeg[r_ind]; k < r_matbeg[r_ind + 1]; k++){
		  if (rows[r_ind].is_redundant){
		     break;
		  }
		  col_ind = r_matind[k];
		  //if (rows[r_ind].vars_checked){
		  // break;
		  //}

		  if (cols[col_ind].var_type != 'F'){
		     termcode = prep_improve_variable(P, col_ind,
						      r_ind, k,
						      (dive_level - 1), TRUE,
                                                      impl_mode, FALSE, 0.0,
                                                      0.0, MAT_ROW_ORDERED);
		     if (PREP_QUIT(termcode)){
			break;
		     }
		  }
	       }
	    }
	 }
	 if (PREP_QUIT(termcode)){
	    break;
	 }
      }
   }

   if (impl_mode){
      P->impl_rows_time += wall_clock(NULL) - mark_time;
   }

   FREE(row_updated_ind);
   FREE(is_row_updated);


   if (PREP_QUIT(termcode)){
      return termcode;
   }

   return PREP_MODIFIED;
}


/*===========================================================================*/
/*===========================================================================*/
/* get the upper and lower bound of the corresponding row */

int prep_get_row_bounds(MIPdesc *mip, int r_ind, double etol)
{

   ROWinfo *rows = mip->mip_inf->rows;

   int j, c_ind, *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;

   double *ub = mip->ub;
   double *lb = mip->lb;

   double a_val;

   rows[r_ind].ub = rows[r_ind].lb = 0.0;

   for (j = r_matbeg[r_ind]; j < r_matbeg[r_ind + 1]; j++){
      a_val = r_matval[j];
      c_ind = r_matind[j];
      if (a_val > etol){
	 if (rows[r_ind].ub < INF){
	    if (ub[c_ind] >= INF){
	       rows[r_ind].ub = INF;
	    } else {
	       rows[r_ind].ub += a_val * ub[c_ind];
	    }
	 }
	 if (rows[r_ind].lb > -INF){
	    if (lb[c_ind] <= -INF){
	       rows[r_ind].lb = -INF;
	    } else {
	       rows[r_ind].lb += a_val * lb[c_ind];
	    }
	 }
      } else if (a_val < -etol){
	 if (rows[r_ind].ub < INF){
	    if (lb[c_ind] <= -INF){
	       rows[r_ind].ub = INF;
	    } else {
	       rows[r_ind].ub += a_val * lb[c_ind];
	    }
	 }
	 if (rows[r_ind].lb > -INF){
	    if (ub[c_ind] >= INF){
	       rows[r_ind].lb = -INF;
	    } else {
	       rows[r_ind].lb += a_val * ub[c_ind];
	    }
	 }
      }
   }

   return 0;

}

/*===========================================================================*/
/*===========================================================================*/

double prep_rnd_integral(double val, double etol, char rnd_type)
{

   double new_bound = 0.0;

   if (rnd_type == RND_FLOOR){
      new_bound = ceil(val);
      if (val < new_bound - etol){
	 new_bound = floor(val);
      }
   } else {
      new_bound = floor(val);
      if (val > new_bound + etol){
	 new_bound = ceil(val);
      }
   }

   return new_bound;
}


/*===========================================================================*/
/* check if this row is redundant */

int prep_check_redundancy(PREPdesc *P, int row_ind,
			  char use_sr_bounds,
			  double sr_ub, double sr_lb, char impl_mode,
			  int dive_level)
{
   int i, termcode = PREP_UNMODIFIED, fix_type = FIX_NO_BOUND;
   int fixed_row = FALSE, fix_all_lb = FALSE, fix_all_ub = FALSE, col_ind;
   int debug_cnt = 0;
   double a_val, ub, lb, new_bound = 0.0, rnd_floor, rnd_ceil;

   MIPdesc *mip = P->mip;
   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;

   char sense = mip->sense[row_ind];
   double rhs = mip->rhs[row_ind];
   double *c_ub = mip->ub;
   double *c_lb = mip->lb;
   double *obj = mip->obj;

   int *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   double *r_matval = mip->row_matval;

   int verbosity = P->params.verbosity;
   double etol = P->params.etol;
   prep_stats *stats = &(P->stats);
   //int dive_level = 1;

  /*
      ratio_type =
      0 c_val >0, a_val>0
      1 c_val >= 0, a_val <= 0
      2 c_val <= 0, a_val >= 0
      3 c_val < 0, a_val < 0
   */

   //   for (i = 0; i < row_cnt; i++){
   //      row_ind = row_start[i];

      if (!use_sr_bounds && rows[row_ind].fixed_var_num  >= rows[row_ind].size){
	 if ((sense == 'L' &&
	      rows[row_ind].fixed_lhs_offset > rhs + etol) ||
	     (sense == 'E' &&
	      !prep_is_equal(rows[row_ind].fixed_lhs_offset, rhs, etol)) ||
	     (sense == 'R' &&
	      (rows[row_ind].fixed_lhs_offset > rhs + etol ||
	       rows[row_ind].fixed_lhs_offset <
	       rhs - mip->rngval[row_ind] - etol))){
	   stats->row_infeas_ind = row_ind;
	   return PREP_INFEAS;
	 }
	 termcode = PREP_MODIFIED;
      } else if (sense != 'R' && !use_sr_bounds &&
	       rows[row_ind].fixed_var_num >= rows[row_ind].size - 1){
	 for (i = r_matbeg[row_ind]; i < r_matbeg[row_ind + 1]; i++){
	    col_ind = r_matind[i];
	    a_val = r_matval[i];
	    if (cols[col_ind].var_type != 'F'){
	       debug_cnt++;
	       if (a_val > etol || a_val < -etol){
		  new_bound = (double)
		     ((rhs - rows[row_ind].fixed_lhs_offset)/a_val);

		  if (sense == 'E'){
		     if (new_bound > c_ub[col_ind] + etol ||
			new_bound < c_lb[col_ind] - etol){
			stats->col_infeas_ind = col_ind;
			stats->row_infeas_ind = row_ind;
			return PREP_INFEAS;
		     }
		     if (cols[col_ind].var_type != 'C'){
			rnd_floor = floor(new_bound);
			rnd_ceil = ceil(new_bound);
			if (new_bound >= rnd_floor + etol &&
			   new_bound <= rnd_ceil - etol){
			   stats->col_infeas_ind = col_ind;
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			} else {
			   if (new_bound < rnd_floor + etol){
			      new_bound = rnd_floor;
			   } else {
			      new_bound = rnd_ceil;
			   }
			}
		     }
		     fix_type = FIX_OTHER;
		  } else {
		     if (cols[col_ind].col_size > 1){
			if ((sense == 'G' && a_val > etol) ||
			   (sense == 'L' && a_val < -etol)){
			   if (new_bound > c_ub[col_ind] + etol){
			      stats->col_infeas_ind = col_ind;
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if (new_bound > c_lb[col_ind] + etol ){
			      if (cols[col_ind].var_type != 'C'){
				 new_bound = prep_rnd_integral(new_bound, etol,
							       RND_CEIL);
			      }
			      fix_type = IMPROVE_LB;
			   }
			} else {
			   if (new_bound < c_lb[col_ind] - etol){
			      stats->col_infeas_ind = col_ind;
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if (new_bound < c_ub[col_ind] - etol){
			      if (cols[col_ind].var_type != 'C'){
				 new_bound= prep_rnd_integral(new_bound, etol,
							      RND_FLOOR);
			      }
			      fix_type = IMPROVE_UB;
			   }
			}
		     } else {
			/*so we can fix this column and delete this row*/
			/* and sense == 'L' */
			if (a_val > etol){
			   if (new_bound < c_lb[col_ind] - etol){
			      stats->col_infeas_ind = col_ind;
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if (obj[col_ind] >= 0.0){
			      new_bound = c_lb[col_ind];
			   } else {
			      if (new_bound > c_ub[col_ind] + etol){
				 new_bound = c_ub[col_ind];
			      } else {
				 if (cols[col_ind].var_type != 'C'){
				    new_bound= prep_rnd_integral(new_bound,
								 etol,
								 RND_FLOOR);
				 }
			      }
			   }
			} else if (a_val < -etol){
			   if (new_bound > c_ub[col_ind] + etol){
			      stats->col_infeas_ind = col_ind;
			      stats->row_infeas_ind = row_ind;
			      return PREP_INFEAS;
			   }
			   if (obj[col_ind] <= 0.0){
			      new_bound = c_ub[col_ind];
			   } else {
			      if (new_bound <  c_lb[col_ind] - etol){
				 new_bound = c_lb[col_ind];
			      } else {
				 if (cols[col_ind].var_type != 'C'){
				    new_bound= prep_rnd_integral(new_bound,
								 etol,
								 RND_CEIL);
				 }
			      }
			   }
			}

			fix_type = FIX_OTHER;
		     }
		  }
	       } else {
		  /* so a_val has been fixed to 0.0 */
		  /* this row is redundant */
		  /* check if we can fix this column*/
		  if (cols[col_ind].col_size == 1){
		     if (sense == 'E'){
			if (!prep_is_equal(rows[row_ind].fixed_lhs_offset, rhs,
					  etol)){
			   stats->col_infeas_ind = col_ind;
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}
		     } else { /* sense is 'L' */
			if (rows[row_ind].fixed_lhs_offset > rhs + etol){
			   stats->col_infeas_ind = col_ind;
			   stats->row_infeas_ind = row_ind;
			   return PREP_INFEAS;
			}

			if (obj[col_ind] >= 0.0){
			   new_bound = c_lb[col_ind];
			} else {
			   new_bound = c_ub[col_ind];
			}
		     }
		     fix_type = FIX_OTHER;
		  } else {
		     /* just declare this row to be redundant
			here*/
		     //cols[col_ind].col_size--;
		     termcode = PREP_MODIFIED;
		  }
	       }

	       if (fix_type != FIX_NO_BOUND){
		  if (cols[col_ind].var_type == 'B'){
		     fix_type = FIX_BINARY;
		  }
		  termcode = prep_modified_cols_update_info(P, 1, &col_ind,
							    row_ind,
							    dive_level,
							    new_bound,
							    fix_type, TRUE,
							    impl_mode);
		  if (PREP_QUIT(termcode)){
		     return termcode;
		  } else if (rows[row_ind].is_redundant){
		     return PREP_MODIFIED;
		  }
	       }
	       termcode = PREP_MODIFIED;
	       /* debug - need to break here */
	       break;
	    }
	 }
	 /* debug */
	 if (debug_cnt > 1){
	    printf("error in prep_check_redundancy()\n");
	    return PREP_OTHER_ERROR;
	 } else if (debug_cnt == 0){
	    /* end point of recursive procedure */
	    /*
	      probably many variables have been fixed during recursive
	      procedure, so just declare redundancy here
	      no need to do anything */
	    termcode = PREP_MODIFIED;
	 }
      }

      if (termcode == PREP_UNMODIFIED){
	 if (use_sr_bounds){
	    ub = sr_ub;
	    lb = sr_lb;
	 } else {
	    ub = rows[row_ind].ub;
	    lb = rows[row_ind].lb;
	 }

	 /*check redundancy and infeasibility*/
	 if (lb > ub + etol){
	    stats->row_infeas_ind = row_ind;
	    /* debug, can this happen? */
	    return PREP_INFEAS;
	 } else if (lb > ub - etol){
	    fixed_row = TRUE;
	    fix_all_ub = TRUE;
	    /* check infeasibility here */
	    if (lb > rhs + etol){
	       stats->row_infeas_ind = row_ind;
	       return PREP_INFEAS;
	    }
	    if (sense == 'E'){
	       if (ub < rhs - etol){
		  stats->row_infeas_ind = row_ind;
		  return PREP_INFEAS;
	       }
	    }
	 }

	 if (!fixed_row){
	    switch(sense){
	     case 'L':
		if (lb > rhs + etol){
		   stats->row_infeas_ind = row_ind;
		   /* prob infeasible */
		   return PREP_INFEAS;
		}
		if (ub < rhs - etol){
		   //rows[row_ind].is_redundant = TRUE;
		   termcode = PREP_MODIFIED;
		}
		if (lb > rhs - etol){
		   fix_all_lb = TRUE;
		}
		break;
	     case 'E':
		if (lb > rhs + etol ||
		   ub < rhs - etol){
		   /* prob infeasible */
		   stats->row_infeas_ind = row_ind;
		   return PREP_INFEAS;
		}
		if (prep_is_equal(lb, rhs, etol*1e-5)){
		   fix_all_lb = TRUE;
		}
		if (prep_is_equal(ub, rhs, etol*1e-5)){
 		   fix_all_ub = TRUE;
		   if(fix_all_lb){
		     if(fabs(ub - rhs) < fabs(rhs - lb)){
		       fix_all_lb = FALSE;
		     }else{
		       fix_all_ub = FALSE;
		     }
		   }
		}
		break;
	     default:
		/* debug */
		break;
	    }
	 }
      }

      if (fix_all_lb || fix_all_ub){
	 if (!use_sr_bounds){

	    /* first temporarily fix the columns
	       this is to detect infeasibilities during
	       recursive procedure */

	    /* keep bounds for recursive infeasibility detection */
	    /* not very efficient */
	    rows[row_ind].is_redundant = TRUE;

	    if (fix_all_lb){
	       fix_type = FIX_ROW_LB;
	    } else {
	       fix_type = FIX_ROW_UB;
	    }

	    //fix_type = FIX_OTHER;

	    /* now go ahead and imply them */
	    termcode = prep_modified_cols_update_info(P, r_matbeg[row_ind + 1]
						      - r_matbeg[row_ind],
						      &r_matind[r_matbeg[row_ind]],
						      row_ind,
						      dive_level, 0.0,
						      fix_type, TRUE, impl_mode);
	    if (PREP_QUIT(termcode)){
	       return termcode;
	    }

	    termcode = PREP_MODIFIED;

	 } else {
	    if (fix_all_lb && fix_all_ub){
	       printf("sr bounds are equal to rhs - row redundant!\n");
	       termcode =  PREP_MODIFIED;
	    }
	 }
      }

      if (termcode == PREP_MODIFIED){
	 stats->rows_deleted++;
	 if (verbosity >= 13){
	    prep_declare_redundant_row(rows[row_ind], row_ind, sense,
				       rhs);
	 }
	 termcode = prep_deleted_row_update_info(mip, row_ind);
	 if (PREP_QUIT(termcode)){
	    return termcode;
	 } else {
	    return PREP_MODIFIED;
	 }

      }
      //   }

   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
/* make these more robust */
int prep_declare_redundant_row(ROWinfo row, int row_ind, char sense,
			       double rhs)
{
   printf("row [%i] is redundant: ", row_ind);
   printf("ub: ");
   if (row.ub < INF){
      printf("%f",row.ub);
   } else {
      printf("INF");
   }
   printf("\t lb: ");
   if (row.lb > -INF){
      printf("%f",row.lb);
   } else {
      printf("-INF");
   }
   printf("\t sense: %c \t rhs: %f\n", sense, rhs);

   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_declare_fixed_var(int col_ind, char *name, double fixed_bound){

   if (name){
      printf("var %s [%i] is fixed to %f\n",
	     name, col_ind, fixed_bound);
   } else {
      printf("var [%i] is fixed to %f\n",
	     col_ind, fixed_bound);
   }
   return 0;

}
/*===========================================================================*/
/*===========================================================================*/
int prep_declare_coef_change(int row_ind, int col_ind,
			     char *name, double a_val,
			     double rhs){
   if (name){
      printf("row [%i] with rhs %f: col %s [%i]: coeff improved to %f\n",
	     row_ind, rhs, name, col_ind, a_val);
   } else {
      printf("row [%i] with rhs %f: col [%i]: coeff improved to %f\n",
	     row_ind, rhs, col_ind, a_val);
   }
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
int prep_deleted_row_update_info(MIPdesc *mip, int row_ind)
{
   /* nothing fancy to do here */
   mip->mip_inf->rows[row_ind].is_redundant = TRUE;

   /* and update col sizes */
   COLinfo *cols = mip->mip_inf->cols;

   int i, end, *r_matbeg = mip->row_matbeg;
   int *r_matind = mip->row_matind;
   //   double *r_matval = mip->row_matval;

   end = r_matbeg[row_ind + 1];
   for (i = r_matbeg[row_ind]; i < end; i++){
      if (cols[r_matind[i]].var_type != 'F'){
	 (cols[r_matind[i]].col_size)--;
	 /* debug */
	 if (cols[r_matind[i]].col_size < 0){
	    printf("error in prep_deleted_row_update_info()\n");
	    return PREP_OTHER_ERROR;
	 }
      }
   }
   return 0;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_initialize_mipinfo(PREPdesc *P)

{
   int i, j;
   double coef_val, fixed_obj_offset;
   int row_ind, cont_var_cnt = 0, bin_var_cnt = 0, fixed_var_cnt = 0;
   int  row_unbounded_cnt, max_row_size, max_col_size, bin_var_nz_cnt = 0;
   int * row_coef_bin_cnt = NULL, *row_sign_pos_cnt = NULL, bin_row_cnt;
   char is_binary, is_bounded, unbounded_below, unbounded_above;
   int gen_type; /* 0 fractional, 1 integer, 2 binary */
   int col_size, col_coef_bin_cnt, col_coef_frac_cnt, col_sign_pos_cnt;
   int *rows_integerized_var_ind = NULL, integerizable_var_num;
   char is_opt_val_integral = TRUE;
   int is_col_all_neg; /* if we convert all constraints to 'L' */
   int is_col_all_pos; /* if we convert all constraints to 'L' */
   int obj_size; /* number of nonzeros in objective function */
   int bin_sos_row_cnt = 0, bin_cont_row_cnt, cont_row_cnt = 0;

   MIPdesc *mip = P->mip;
   prep_stats *stats = &(P->stats);
   prep_params params = P->params;

   double etol = params.etol;
   double coeff_etol = 1e-15;
   int verbosity = params.verbosity;
   int p_level = params.level;
   //   int termcode;

   /* fixme! objsense min max issue!!! will always assume that it is a
      min problem here!!!!
   */

   if (!mip){
      if (verbosity >= 1){
	 printf("prep_initialize_mipinfocollect_mipinfo():"
		"Empty mip description...\n");
      }
      return(PREP_OTHER_ERROR);
   } else{
     if (mip->n < 1 && p_level > 2 ){
       if (verbosity >= 1){
	 printf("Empty problem...\n");
       }
       return PREP_SOLVED;
     }
   }

   int n = mip->n;
   int m = mip->m;
   int * matbeg = mip->matbeg;
   int *matind = mip->matind;
   double *matval = mip->matval;
   double * ub = mip->ub;
   double *lb = mip->lb;
   char *is_int = mip->is_int;
   double *obj = mip->obj;
   char *sense = mip->sense;
   double *rhs = mip->rhs;
   //   char obj_sense = env->mip->obj_sense;

   MIPinfo *mip_inf = (MIPinfo *)calloc (1, sizeof(MIPinfo));
   COLinfo *cols = NULL;
   ROWinfo *rows = NULL;

   if (m > 0){
      rows = (ROWinfo *)calloc(m, sizeof(ROWinfo));
      row_coef_bin_cnt = (int *)calloc(ISIZE,m);
      row_sign_pos_cnt = (int *)calloc(ISIZE,m);
      rows_integerized_var_ind = (int *)malloc(ISIZE*m);
   }
   if (n > 0){
      cols = (COLinfo *)calloc(n, sizeof(COLinfo));
   }

   /*
      Number of continuous variables that will always have an integer value in
      a solution. For now, we check rows with only one cont var
      this might be further improved in advanced option of preprocessor
   */
   integerizable_var_num = 0;

   max_col_size = 0;
   fixed_obj_offset = 0;
   obj_size = 0;

   for (i = 0; i < n; i++){
      is_binary = FALSE;
      is_bounded = FALSE;
      unbounded_below = FALSE;
      unbounded_above = FALSE;

      if (fabs(obj[i]) > etol) {
         obj_size++;
      }

      cols[i].var_type = 'I';
      /* check if a variable has bad bounds */
      if (lb[i] >= ub[i] + etol && p_level > 2) {
	 stats->col_infeas_ind = i;
         /* fixme: who will free the above allocs? */
	 return(PREP_INFEAS);
      }

      /* check for unboundedness */
      if ((lb[i] >= INF || ub[i] <= -INF) && p_level > 2) {
        stats->col_numeric_ind = i;
        /* fixme: who will free the above allocs? */
        return(PREP_NUMERIC_ERROR);
      }

      /* check if a variable can be fixed because of its bounds */
      if (lb[i] > (ub[i] - etol)) {
	 cols[i].var_type = 'F';
	 fixed_obj_offset += obj[i]*ub[i];
	 fixed_var_cnt++;
      } else if (is_int[i]) {
	 if (lb[i] > (-1.0 + etol) &&
	    ub[i] < (2.0 - etol)){
	    is_binary = is_bounded = TRUE; //is_lb_zero = TRUE;
	    cols[i].var_type = 'B';
	    bin_var_cnt++;
	    bin_var_nz_cnt += matbeg[i+1] - matbeg[i];
	 } else if (lb[i] > (-2.0 + etol) && ub[i] < (1.0 - etol)) {
           /* fixme: move this to later. there may be other variables as well
            * that can be reduced to binary variables. */
	    is_binary = is_bounded = TRUE; //is_lb_zero = TRUE;
	    cols[i].var_type = 'R';
	    bin_var_cnt++;
	    bin_var_nz_cnt += matbeg[i+1] - matbeg[i];
	 }
      } else {
	 cols[i].var_type = 'C';
	 cont_var_cnt++;
      }

      if (!is_bounded) {
	 if (lb[i] <= -INF) {
	    unbounded_below = TRUE;
	 }
	 if (ub[i] >= INF) {
	    unbounded_above = TRUE;
	 }
	 if (!unbounded_below && !unbounded_above) {
	    is_bounded = TRUE;
	 }
      }

      col_coef_bin_cnt = 0;
      col_coef_frac_cnt = 0;
      col_sign_pos_cnt = 0;
      is_col_all_neg = TRUE;
      is_col_all_pos = TRUE;

      for (j = matbeg[i]; j < matbeg[i+1]; j++) {
	 row_ind = matind[j];
	 coef_val = matval[j];
	 rows[row_ind].size++;

	 /* for types of variables in a row */
	 if (cols[i].var_type == 'F'){
	    rows[row_ind].fixed_var_num++;
         } else if (is_int[i]) {
            if (is_binary){
              rows[row_ind].bin_var_num++;
            }
         } else {
           rows[row_ind].cont_var_num++;
           if (rows[row_ind].cont_var_num < 2){
             rows_integerized_var_ind[row_ind] = i;
           }
         }

	 /* for bounds on the activity of a row */
	 if (!is_bounded) {
	    if (unbounded_above) {
	       if (coef_val > 0.0) {
		  rows[row_ind].ub_inf_var_num++;
	       } else {
		  rows[row_ind].lb_inf_var_num++;
	       }
	    }
	    if (unbounded_below) {
	       if (coef_val > 0.0) {
		  rows[row_ind].lb_inf_var_num++;
	       } else {
		  rows[row_ind].ub_inf_var_num++;
	       }
	    }
	 }

	 /* for coef types */
	 if (fabs(coef_val-floor(coef_val+0.5)) > coeff_etol) {
	    rows[row_ind].frac_coef_num++;
	    col_coef_frac_cnt++;
	 } else if (fabs(coef_val - 1.0) < coeff_etol ||
	            fabs(coef_val + 1.0) < coeff_etol) {
            row_coef_bin_cnt[row_ind]++;
            col_coef_bin_cnt++;
         }

	 /* for sign types and update bounds */
	 if (coef_val > 0.0) {
	    row_sign_pos_cnt[row_ind]++;
	    col_sign_pos_cnt++;
	    if (rows[row_ind].ub < INF) {
	       if (ub[i] >= INF) {
		  rows[row_ind].ub = INF;
	       } else {
		  rows[row_ind].ub += coef_val * ub[i];
	       }
	    }
	    if (rows[row_ind].lb > -INF) {
	       if (lb[i] <= -INF) {
		  rows[row_ind].lb = -INF;
	       } else {
		  rows[row_ind].lb += coef_val * lb[i];
	       }
	    }

	    if (is_col_all_neg) {
               /* fixme: this is never going to be the case, since we
                * eliminated >= constraints? */
	       if (sense[row_ind] != 'G') {
		  is_col_all_neg = FALSE;
	       }
	    }
	    if (is_col_all_pos) {
	       if (sense[row_ind] != 'L') {
		  is_col_all_pos = FALSE;
	       }
	    }

	 } else if (coef_val < 0.0) {
	    if (rows[row_ind].ub < INF) {
	       if (lb[i] <= -INF) {
		  rows[row_ind].ub = INF;
	       } else {
		  rows[row_ind].ub += coef_val * lb[i];
	       }
	    }
	    if (rows[row_ind].lb > -INF){
	       if (ub[i] >= INF){
		  rows[row_ind].lb = -INF;
	       } else {
		  rows[row_ind].lb += coef_val * ub[i];
	       }
	    }

	    if (is_col_all_neg) {
	       if (sense[row_ind] != 'L') {
		  is_col_all_neg = FALSE;
	       } else {
                 /* do nothing */
	       }
	    }
	    if (is_col_all_pos){
              /* fixme: this is never going to be the case, since we
               * eliminated >= constraints? */
	       if (sense[row_ind] != 'G'){
		  is_col_all_pos = FALSE;
	       }
	    }
	 }

	 if (cols[i].var_type == 'F'){
	    rows[row_ind].fixed_obj_offset += obj[i]*ub[i];
	    rows[row_ind].fixed_lhs_offset += coef_val * ub[i];
	 }
      }


      /* check if this column is fixable */
      /* if so set type, update the obj_offset etc */

      col_size = cols[i].col_size = matbeg[i+1] - matbeg[i];

      if (is_col_all_neg || col_size <= 0){
	 //if ((obj[i] > 0.0 && obj_sense == SYM_MAXIMIZE) ||
	 //(obj[i] < 0.0 && obj_sense == SYM_MINIMIZE)){
	 if (obj[i] < 0.0){
	    if (ub[i] >= INF && p_level > 2){
	       stats->col_unbound_ind = i;
	       return(PREP_UNBOUNDED); /* fixme: unbounded return code */
	    } else {
	       /* total obj offset here is for fixed variables
		  later(if prep.is used) will include 'U' and 'L'
		  variables */
	       //total_obj_offset += obj[i] * ub[i];
	       if (verbosity >= 12){
		  if (mip->colname){
		     printf("var %s [%i] is fixable to its upper bound: %f\n",
			    mip->colname[i], i, ub[i]);
		  } else {
		     printf("var [%i] is fixable to its upper bound: %f\n",
			    i, ub[i]);
		  }
		  cols[i].var_type = 'U';
	       }
	    }
	 }
      }
      if (is_col_all_pos || col_size <= 0){
	 if (obj[i] > 0.0){
	    if (lb[i] <= -INF && p_level > 2){
	       stats->col_unbound_ind = i;
	       return(PREP_UNBOUNDED);
	    } else {
	       //total_obj_offset += obj[i] * lb[i];
	       if (verbosity >= 12){
		  if (mip->colname){
		     printf("var %s [%i] is fixable to its lower bound: %f\n",
			    mip->colname[i], i, lb[i]);
		  } else {
		     printf("var [%i] is fixable to its lower bound: %f\n",
			    i, lb[i]);
		  }
		  cols[i].var_type = 'L';
	       }
	    }
	 }
      }

      /* fill in col info - if not a fixed variable */
      if (col_size && cols[i].var_type != 'F'){

	 if (col_size > max_col_size){
	    max_col_size = col_size;
	 }

	 gen_type = 0;
	 if (col_coef_frac_cnt > 0){
	    gen_type = FRACTIONAL_VEC;
	 } else {
	    if (col_coef_bin_cnt < col_size){
	       gen_type = ALL_INTEGER_VEC;
	    } else {
	       gen_type = ALL_BINARY_VEC;
	    }
	 }
	 cols[i].coef_type = gen_type;

	 gen_type = 0;
	 if (col_sign_pos_cnt > 0){
	    if (col_sign_pos_cnt < col_size){
	       gen_type = MIXED_TYPE_VEC;
	    } else {
	       gen_type = ALL_POS_VEC;
	    }
	 } else {
	    gen_type = ALL_NEG_VEC;
	 }
	 cols[i].sign_type = gen_type;
      }
   }

   /* fill in row info */

   max_row_size = 0;
   bin_row_cnt = 0;
   cont_row_cnt = 0;
   bin_cont_row_cnt = 0;

   for (j = 0; j < m; j++){

      if (rows[j].size > max_row_size){
	 max_row_size = rows[j].size;
      }

      gen_type = 0;
      if (rows[j].cont_var_num > 0 ){
	 if (rows[j].bin_var_num > 0){
	    if ( rows[j].cont_var_num + rows[j].bin_var_num +
		rows[j].fixed_var_num
		>= rows[j].size){
	       gen_type = BIN_CONT_TYPE;
	    } else {
	       gen_type = ALL_MIXED_TYPE;
	    }
	    bin_row_cnt++;
	    bin_cont_row_cnt++;
	 } else {
	    if (rows[j].cont_var_num + rows[j].fixed_var_num < rows[j].size){
	       gen_type = INT_CONT_TYPE;
	    } else {
	       gen_type = CONTINUOUS_TYPE;
	    }
	 }
	 cont_row_cnt++;
      } else {
	 if (rows[j].bin_var_num > 0){
	    if (rows[j].bin_var_num + rows[j].fixed_var_num < rows[j].size){
	       gen_type = BIN_INT_TYPE;
	    } else {
	       gen_type = BINARY_TYPE;
	    }
	    bin_row_cnt++;
	 } else {
	    gen_type = INTEGER_TYPE;
	 }
      }

      rows[j].type = gen_type;
      gen_type = 0;

      row_unbounded_cnt = rows[j].lb_inf_var_num +
	 rows[j].ub_inf_var_num;

      if (row_unbounded_cnt == 0){
	 gen_type = ALL_BOUNDED_ROW;
      } else if (row_unbounded_cnt + rows[j].fixed_var_num < rows[j].size){
	 gen_type = MIXED_BOUNDED_ROW;
      } else {
	 gen_type = OPEN_ROW;
      }

      rows[j].bound_type = gen_type;
      gen_type = 0;

      if (rows[j].frac_coef_num > 0){
	 gen_type = FRACTIONAL_VEC;
      } else {
	 if (row_coef_bin_cnt[j] +  rows[j].fixed_var_num < rows[j].size){
	    gen_type = ALL_INTEGER_VEC;
	 } else {
	    gen_type = ALL_BINARY_VEC;
	 }
      }

      rows[j].coef_type = gen_type;
      gen_type = 0;

      if (row_sign_pos_cnt[j] > 0){
	 if (row_sign_pos_cnt[j] + rows[j].fixed_var_num < rows[j].size){
	    gen_type = MIXED_TYPE_VEC;
	 } else {
	    gen_type = ALL_POS_VEC;
	 }
      } else {
	 gen_type = ALL_NEG_VEC;
      }

      rows[j].sign_type = gen_type;


      if (rows[j].coef_type == ALL_BINARY_VEC &&
	 ((rows[j].sign_type == ALL_POS_VEC &&
	   (sense[j] == 'E' || sense[j] == 'L') &&
	   rhs[j] > 0 && rhs[j] < 2.0) ||
	  (rows[j].sign_type == ALL_NEG_VEC &&
	   (sense[j] == 'E' || sense[j] == 'G') &&
	   rhs[j] < 0 && rhs[j] > -2.0))){
	 rows[j].is_sos_row = TRUE;
	 bin_sos_row_cnt++;
      }

      /* work on integerizable vars */

      if (sense[j] == 'E' && rows[j].cont_var_num == 1 &&
	 rows[j].coef_type != FRACTIONAL_VEC &&
	 prep_is_integral(rhs[j], coeff_etol) &&
	 prep_is_integral(rows[j].fixed_lhs_offset, coeff_etol)){
	 if (cols[rows_integerized_var_ind[j]].var_type != 'Z'){
	    cols[rows_integerized_var_ind[j]].var_type = 'Z';
	    integerizable_var_num++;
	 }
      }

      rows[j].sr_ub = rows[j].ub;
      rows[j].sr_lb = rows[j].lb;
   }

  /* work on obj */

   if (!(cont_var_cnt - integerizable_var_num)){
      for (i = 0; i < n; i++){
	 coef_val = obj[i];
	 if (!prep_is_integral(coef_val, coeff_etol)){
	    if (cols[i].var_type == 'F'){
	       if (ub[i] < etol && ub[i] > -etol){
		  continue;
	       }
	    }
	    is_opt_val_integral = FALSE;
	    break;
	 }
      }
   } else {
      is_opt_val_integral = FALSE;
   }

   /* finally prob type */
   gen_type = 0;
   if (cont_var_cnt > 0 ){
      if (bin_var_cnt > 0){
	 if (cont_var_cnt + bin_var_cnt + fixed_var_cnt < n){
	    gen_type = ALL_MIXED_TYPE;
	 } else {
	    gen_type = BIN_CONT_TYPE;
	 }
      } else {
	 if (cont_var_cnt + fixed_var_cnt < n){
	    gen_type = INT_CONT_TYPE;
	 } else {
	    gen_type = CONTINUOUS_TYPE;
	 }
      }
   } else {
      if (bin_var_cnt > 0){
	 if (bin_var_cnt + fixed_var_cnt < n){
	    gen_type  = BIN_INT_TYPE;
	 } else {
	    gen_type = BINARY_TYPE;
	 }
      } else {
	 gen_type = INTEGER_TYPE;
      }
   }

   mip_inf->prob_type = gen_type;
   mip_inf->cont_var_num = cont_var_cnt;
   mip_inf->binary_var_num = bin_var_cnt;
   mip_inf->binary_var_nz = bin_var_nz_cnt;
   mip_inf->int_var_ratio = (1.0*(n - cont_var_cnt))/(n + 1);
   mip_inf->cont_var_ratio = 1.0*cont_var_cnt/(n + 1);
   mip_inf->bin_var_ratio = 1.0*bin_var_cnt/(n +1);
   mip_inf->fixed_var_num = fixed_var_cnt;
   mip_inf->max_row_size = max_row_size;
   mip_inf->max_col_size = max_col_size;
   mip_inf->max_row_ratio = 1.0*max_row_size/(n+1);
   mip_inf->max_col_ratio = 1.0*max_col_size/(m+1);
   mip_inf->obj_size = obj_size;
   mip_inf->mat_density = 1.0*mip->nz/(n*m + 1);
   mip_inf->row_density = 1.0*mip->nz/(m+1);
   mip_inf->col_density = 1.0*mip->nz/(n+1);
   mip_inf->integerizable_var_num = integerizable_var_num;
   mip_inf->is_opt_val_integral = is_opt_val_integral;
   mip_inf->sum_obj_offset = fixed_obj_offset;
   mip_inf->binary_row_num = bin_row_cnt;
   mip_inf->binary_sos_row_num = bin_sos_row_cnt;
   mip_inf->cont_row_num = cont_row_cnt;
   mip_inf->bin_cont_row_num = bin_cont_row_cnt;

   mip_inf->sos_bin_row_ratio = 1.0*bin_sos_row_cnt/(bin_row_cnt +1);
   mip_inf->bin_row_ratio = 1.0*bin_row_cnt/(m +1);

   //   printf("m, cont_row_num: %i\t%i\n", m, cont_row_cnt);
   if (bin_var_cnt){
      mip_inf->row_bin_den = (int)
	 (bin_var_nz_cnt/bin_row_cnt) + 1;

      if (bin_var_cnt < n){
      	 mip_inf->row_bin_den = (int)(mip_inf->row_bin_den*
				      n/bin_var_cnt) + 1;
      }

      mip_inf->col_bin_den = (int)
	 (bin_var_nz_cnt/bin_var_cnt) + 1;

      mip_inf->row_bin_den_mean =
	 (int)2*mip_inf->row_bin_den * max_row_size/
	 (mip_inf->row_bin_den + max_row_size) + 1;

      mip_inf->col_bin_den_mean =
	 (int)2*mip_inf->col_bin_den * max_row_size/
	 (mip_inf->col_bin_den + max_col_size) + 1;


   }

   mip_inf->rows = rows;
   mip_inf->cols = cols;

   if (mip->mip_inf){
      FREE(mip->mip_inf->rows);
      FREE(mip->mip_inf->cols);
      FREE(mip->mip_inf);
   }

   mip->mip_inf = mip_inf;
   /*
   double mat_den = (1.0)* mip->nz/(mip->m * mip->n + 1);
   double int_den = (1.0*(mip->n-mip->mip_inf->cont_var_num))/(mip->n + 1);
   double max_row_col_den = (1.0*mip->mip_inf->max_col_size *
			     mip->mip_inf->max_row_size)/(mip->n * mip->m);
   double max_col_den = 1.0*mip->mip_inf->max_col_size/(mip->m + 1);
   double max_row_den = 1.0*mip->mip_inf->max_row_size/(mip->n + 1);
   printf("mat_den - int_den  - row_col_den: %f %f %f \n",
	  mat_den, int_den, max_row_col_den);
   printf("col_den - row_den %f %f \n", max_col_den, max_row_den);

   if ((mat_den < 0.05 && int_den > 0.05 && (max_col_den > 0.05 ||
					    max_row_den > 0.05))
      || mip->nz > 1e5){
      printf("TRUE");
   }
   */

   FREE(row_coef_bin_cnt);
   FREE(row_sign_pos_cnt);
   FREE(rows_integerized_var_ind);

   return(PREP_MODIFIED);
}

/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/
int prep_integerize_bounds(PREPdesc *P)
{
   /* Change the bounds of integer variables to floor/ceiling appropriately */
   int termcode = 0;
   MIPdesc *mip = P->mip;
   MIPinfo *mip_inf = mip->mip_inf;
   COLinfo *cols = mip_inf->cols;
   int i, b_cnt = 0, n = mip->n;
   double *ub = mip->ub;
   double *lb = mip->lb;
   double temp_fl, temp_cl;
   double diff_ub, diff_lb;
   double etol = P->params.etol;
   int verbosity = P->params.verbosity;

   //   int * bounds_updated = (int *)calloc(ISIZE,n);

   if (P->params.level >= 6 && mip_inf->integerizable_var_num){
      for (i = 0; i < n; i++) {
	 if (cols[i].var_type == 'Z'){
	    termcode = prep_integerize_var(P, i);
	    if (PREP_QUIT(termcode)){
	       return termcode;
	    }
	 }
      }
   }

   for (i = 0; i < n; i++) {
     if (cols[i].var_type != 'F' &&
	 cols[i].var_type != 'C'){
       if(!(mip->is_int[i] || cols[i].var_type == 'Z')){
	   continue;
         }
	 diff_ub = diff_lb = 0.0;
	 if (ub[i] < INF) {
	    temp_fl = floor(ub[i]);
	    temp_cl = ceil(ub[i]);
	    if (temp_cl - ub[i] < etol ){
	       ub[i] = temp_cl;
	    } else {
	       diff_ub = ub[i] - temp_fl;
	       ub[i] = temp_fl;
	    }
	 }
	 if (lb[i] > -INF){
	    temp_fl = floor(lb[i]);
	    temp_cl = ceil(lb[i]);
	    if (lb[i] - temp_fl < etol){
	       lb[i] = temp_fl;
	    } else {
	       diff_lb = temp_cl - lb[i];
	       lb[i] = temp_cl;
	    }
	 }
	 if (diff_ub >= etol || diff_lb >= etol ){
	    if (ub[i] > lb[i] - etol && ub[i] < lb[i] + etol){
	       if (cols[i].var_type =='B'){
		  mip_inf->binary_var_num--;
		  mip_inf->binary_var_nz -= mip->matbeg[i+1] - mip->matbeg[i];
	       }
	       mip_inf->fixed_var_num++;
	       cols[i].var_type = 'F';
	    }
	    b_cnt++;
	    if (verbosity>=11) {
	       if (mip->colname){
		  printf("integerized bounds [lb-ub] of variable %s:"
			 "%f - %f\n",
			 mip->colname[i],lb[i],ub[i]);
	       } else {
		  printf("integerized bounds [lb-ub] of variable: "
			 "%f - %f\n",
			 lb[i],ub[i]);
	       }
	    }
	 }
      }
   }

#if 0
   if (keeptrack){
      P->mip_diff->bounds_integerized_num = b_cnt;
      P->mip_diff->bounds_integerized_ind = bounds_updated;
   } else {
      FREE(bounds_updated);
   }
#endif
   P->stats.bounds_integerized = b_cnt;
   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/

int prep_integerize_var(PREPdesc *P, int col_ind) {

   int j, k, row_ind, termcode = PREP_MODIFIED;
   MIPdesc *mip = P->mip;
   ROWinfo *rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;
   double etol = P->params.etol;
   double coeff_etol = 1e-15;

   if (P->params.verbosity >= 11 ){
      printf("col %i is integerized\n", col_ind);
   }

   (P->stats.vars_integerized)++;
   mip->is_int[col_ind] = TRUE;
   cols[col_ind].var_type = 'I';
   if (mip->lb[col_ind] > (-1.0 + etol) &&
      mip->ub[col_ind] < (2.0 - etol)){
      cols[col_ind].var_type = 'B';
   }
   for (j = mip->matbeg[col_ind];
       j < mip->matbeg[col_ind+1]; j++){
      row_ind = mip->matind[j];
      if (cols[col_ind].var_type == 'B'){
	 rows[row_ind].bin_var_num++;
      }
      rows[row_ind].cont_var_num--;
      if (rows[row_ind].cont_var_num < 0){
	 printf("error: prep_integerize_var()\n");
	 return PREP_OTHER_ERROR;
      } else if (rows[row_ind].cont_var_num < 1){
	 if (rows[row_ind].bin_var_num){
	    if (rows[row_ind].bin_var_num +
	       rows[row_ind].fixed_var_num
	       >= rows[row_ind].size){
	       rows[row_ind].type = BINARY_TYPE;
	    } else {
	       rows[row_ind].type = BIN_INT_TYPE;
	    }
	 } else {
	    rows[row_ind].type = INTEGER_TYPE;
	 }
      } else if (rows[row_ind].cont_var_num == 1){
	 if (mip->sense[row_ind] == 'E' &&
	    rows[row_ind].coef_type != FRACTIONAL_VEC &&
	    prep_is_integral(mip->rhs[row_ind], coeff_etol) &&
	    prep_is_integral(rows[row_ind].fixed_lhs_offset, coeff_etol)){
	    for (k = mip->row_matbeg[row_ind];
		k < mip->row_matbeg[row_ind + 1]; k++){
	       if (cols[mip->row_matind[k]].var_type == 'C'){
		  termcode = prep_integerize_var(P, mip->row_matind[k]);
		  break;
	       }
	    }
	 }
      }
      if (PREP_QUIT(termcode)){
	 break;
      }
   }
   return termcode;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_fill_row_ordered(PREPdesc *P)
{
   /*
     recreates 'A' matrix using three matrices just like the standard
     notation. However, matrices contain row representations rather than
     column representation.

     Also replace any >= constraints by <=.
   */

   int i, j, *o_ind, *c_lengths;
   int row_ind, elem_ind, *matind, *matbeg, *r_matind, *r_matbeg, *r_lengths;
   double * matval, *r_matval, *rhs;
   MIPdesc *mip = P->mip;
   int n = mip->n;
   int m = mip->m;
   int nz = mip->nz;
   char *sense, *o_sense;
   int *u_col_ind, *u_row_ind;

   matval = mip->matval;
   matind = mip->matind;
   matbeg = mip->matbeg;
   rhs = mip->rhs;
   sense = mip->sense;

   /* allocate space for different arrays */

   r_matval = (mip->row_matval = (double *)malloc(nz*DSIZE));
   r_matind = (mip->row_matind = (int *)malloc(nz*ISIZE));
   r_matbeg = (mip->row_matbeg = (int *)malloc((m+1)*ISIZE));
   r_lengths = (mip->row_lengths = (int *)calloc(m,ISIZE));
   o_sense = (mip->orig_sense = (char *)malloc(m *CSIZE));
   o_ind = (mip->orig_ind = (int *)malloc(n*ISIZE));
   u_col_ind = (P->user_col_ind) = (int *)malloc(n*ISIZE);
   u_row_ind = (P->user_row_ind) = (int *)malloc(m*ISIZE);
   c_lengths = (mip->col_lengths = (int *)calloc(n,ISIZE));
   /* these are initialized here, we have to visit this function anyway */

   /* first get row legths */
   for (i = 0; i < n; i++){
      /* set orig indices here */
      o_ind[i] = u_col_ind[i] = i;
      for (j = matbeg[i]; j < matbeg[i+1]; j++){
	 r_lengths[matind[j]]++;
      }
      c_lengths[i] = matbeg[i+1] - matbeg[i];
   }

   r_matbeg[0] = 0;

   /* fill in matbegs */
   for (i = 0; i < m; i++){
      u_row_ind[i] = i;
      r_matbeg[i + 1] = r_matbeg[i] + r_lengths[i];
   }

   /* get matrix, change 'G' rows to 'L'*/
   for (i = 0; i < n; i++){
      for (j = matbeg[i]; j < matbeg[i+1]; j++){
	 row_ind = matind[j];
	 elem_ind = r_matbeg[row_ind];
	 r_matind[elem_ind] = i;
	 if (sense[row_ind] == 'G'){
	    matval[j] = -matval[j];
	 }
	 r_matval[elem_ind] = matval[j];
	 r_matbeg[row_ind] = elem_ind + 1;
      }
   }
   /* and update matbegs, rhs, and rows with 'G' to 'L'*/
   memcpy(o_sense, sense, CSIZE*m);

   for (i = 0; i < m; i++){
      r_matbeg[i] -= r_lengths[i];
      if (sense[i] == 'G'){
	 sense[i] = 'L';
	 rhs[i] = -rhs[i];
      }
   }

   return PREP_UNMODIFIED;
}
/*===========================================================================*/
/*===========================================================================*/
int prep_cleanup_desc(PREPdesc *P)
{

   int i, j, col_nz, col_num, row_num, fixed_nz, *fixed_ind, *o_ind;
   int row_ind, elem_ind, *matind, *matbeg, *r_matind, *r_matbeg, *r_lengths;
   double *ub, *lb, *matval, *r_matval, *obj, *rhs, *rngval, *fixed_val;
   double obj_offset, debug_offset;
   int new_del_cnt, *c_lengths, binary_var_nz, binary_var_num, bin_row_cnt,
      cont_row_cnt, bin_cont_row_cnt;
   int sos_row_cnt;
   int termcode = PREP_UNMODIFIED;

   MIPdesc *mip = P->mip;
   int n = mip->n;
   int r_ind, m = mip->m;
   int nz = mip->nz;
   char *is_int, *sense, *o_sense, **colnames;
   ROWinfo * rows = mip->mip_inf->rows;
   COLinfo *cols = mip->mip_inf->cols;
   //  if (!rows){
      /* debug */
   //  printf("error in prep_fill_row_ordered - 1");
   //  return PREP_OTHER_ERROR;
   // }

   prep_stats * stats = &(P->stats);
   prep_params params = P->params;

   double etol = params.etol;
   double coeff_etol = 1e-15;
   int deleted_row_cnt = stats->rows_deleted;
   //int vars_fixed = stats->vars_fixed;
   //int keep_row_ordered = params.keep_row_ordered;
   //int reduce_mip = params.reduce_mip;

   int max_row_size, max_col_size;
   int old_start, *row_new_inds = NULL;
   int fixed_num = stats->vars_fixed + mip->mip_inf->fixed_var_num;

   obj_offset = 0.0;
   debug_offset = 0.0;
   sense = mip->sense;
   rngval = mip->rngval;
   o_sense = mip->orig_sense;
   rhs = mip->rhs;
   obj = mip->obj;
   ub = mip->ub;
   lb = mip->lb;
   is_int = mip->is_int;

   if (!fixed_num && !stats->rows_deleted){
     return termcode;
   }

   fixed_nz = 0;
   fixed_ind = mip->fixed_ind = (int *)malloc(n*ISIZE);
   fixed_val = mip->fixed_val = (double *)malloc(n*DSIZE);
   o_ind = mip->orig_ind;


   if (!params.reduce_mip || fixed_num == n || stats->rows_deleted == m){
      if (fixed_num == n || stats->rows_deleted == m){
	 /* get fixed nz vals */
	 for (i = 0; i < n; i++){
	    if (cols[i].var_type == 'F'){
	       if (!prep_is_equal(mip->ub[i], 0.0, etol)){
		  fixed_ind[fixed_nz] = i;
		  fixed_val[fixed_nz++] = mip->ub[i];
	       }
	    } else {
	       if (obj[i] > 0.0){
		  if (lb[i] <= -INF){
		    termcode = PREP_UNBOUNDED;
		    break;
		  } else {
		     if (!prep_is_equal(mip->lb[i], 0.0, etol)){
			fixed_ind[fixed_nz] = i;
			fixed_val[fixed_nz++] = mip->lb[i];
			obj_offset += obj[i]*mip->lb[i];
		     }
		  }
	       } else if (obj[i] < 0.0){
		  if (ub[i] >= INF){
 		     termcode = PREP_UNBOUNDED;
		     break;
		  } else {
		     if (!prep_is_equal(mip->ub[i], 0.0, etol)){
			fixed_ind[fixed_nz] = i;
			fixed_val[fixed_nz++] = mip->ub[i];
			obj_offset += obj[i]*mip->ub[i];
		     }
		  }
	       }
	    }
	 }

	 if(termcode == PREP_UNMODIFIED){
	   mip->fixed_n = fixed_nz;
	   mip->obj_offset = mip->mip_inf->sum_obj_offset + obj_offset;
	   return PREP_SOLVED;
	 }else{
	   mip->fixed_n = 0;
	   FREE(mip->fixed_ind);
	   FREE(mip->fixed_val);
	   return termcode;
	 }
      }
      return termcode;
   }

   row_new_inds = (int *)calloc(m,ISIZE);

   mip->alloc_n = n;
   mip->alloc_m = m;
   mip->alloc_nz = nz;

   matval = mip->matval;
   matind = mip->matind;
   matbeg = mip->matbeg;

   colnames = mip->colname;
   binary_var_nz = 0;
   binary_var_num = 0;

   /* first get new row indices */

   row_num = 0;

   for (i = 0; i < m; i++){
      if (!(rows[i].is_redundant)){
	 row_new_inds[i] = row_num;
	 row_num++;
      }
      rows[i].size = 0; /* might have zero matval vals that we want
			   to exclude, so we will re-evaluate */
   }

   /* debug */
   /* ------ */
   if (row_num != m - deleted_row_cnt){
      printf("error: missing rows \n");
      return PREP_OTHER_ERROR;
   }
   /* ------ */

   /* first fill in col ordered */

   col_nz = col_num = 0;
   old_start = 0;
   for (i = 0; i < n; i++){
       if (cols[i].var_type != 'F'){
	 for (j = old_start; j < matbeg[i+1]; j++){
	    r_ind = matind[j];
	    if (!(rows[r_ind].is_redundant)){
	       if (!prep_is_equal(matval[j], 0.0, coeff_etol)){
		  matind[col_nz] = row_new_inds[r_ind];
		  matval[col_nz] = matval[j];
		  rows[r_ind].size++;
		  col_nz++;
	       }
	    }
	 }

	 /* check if this column has any nonzeros, otherwise
	    fix it*/

	 if (col_nz == matbeg[col_num]){
	    cols[i].var_type = 'F';
	    if (obj[i] >= 0.0){
	       obj_offset += obj[i]*lb[i];
	       fixed_val[fixed_nz] = lb[i];
	    } else {
	       obj_offset += obj[i]*ub[i];
	       fixed_val[fixed_nz] = ub[i];
	    }
	    if (!prep_is_equal(fixed_val[fixed_nz], 0.0, etol)){
	       fixed_ind[fixed_nz++] = i;
	    }
	    stats->vars_fixed++;
	 } else {
	    o_ind[col_num] = i;
	    obj[col_num] = obj[i];
	    ub[col_num] = ub[i];
	    lb[col_num] = lb[i];
	    is_int[col_num] = is_int[i];
	    if (col_num != i){
	       cols[col_num] = cols[i];
	       cols[i].ulist = cols[i].llist = 0;
	       if (colnames){
		  strcpy(colnames[col_num], colnames[i]);
	       }
	    }
	    cols[col_num].col_size = col_nz - matbeg[col_num];
	    old_start = matbeg[i+1];
	    matbeg[(++col_num)] = col_nz;
	    /* debug */
	    /* ----------- */
	    if (cols[col_num - 1].col_size <= 0){
	       printf("error: empty size column \n");
	       return PREP_OTHER_ERROR;
	    }
	    /* ---------- */
	 }

	 if (cols[i].var_type == 'B'){
	    binary_var_num++;
	    binary_var_nz += matbeg[col_num] - matbeg[col_num -1];
	 }
      } else {

	 /* debug */
	 /*---- -*/
	 if (stats->vars_aggregated <= 0){
	    if (!prep_is_equal(ub[i], lb[i], etol)){
	       printf("error: not fixed column? \n");
	       return PREP_OTHER_ERROR;
	    }
	 }
	 /* ----- */
	 debug_offset += obj[i]*ub[i];
	 old_start = matbeg[i+1];
	 if (!prep_is_equal(ub[i], 0.0, etol)){
	    fixed_ind[fixed_nz] = i;
	    fixed_val[fixed_nz++] = ub[i];
	 }
      }
   }

   /* debug */
   /* -------------- */
   if (col_num != n - (stats->vars_fixed + mip->mip_inf->fixed_var_num)){
      printf("error: missing cols \n");
      return PREP_OTHER_ERROR;
   }
   /* ---------- */

   /* now update row_info */
   row_num = 0;
   new_del_cnt = 0;
   for (i = 0; i < m; i++){
      if (!(rows[i].is_redundant) && rows[i].size >= 0){
	 row_new_inds[row_num + new_del_cnt] = row_num;
	 if (rows[i].size == 0){
	    rows[i].is_redundant = TRUE;
	    new_del_cnt++;
	 } else {
	    if (i != row_num){
	       rows[row_num] = rows[i];
	       sense[row_num] = sense[i];
	       if (sense[row_num] == 'R'){
		  rngval[row_num] = rngval[i];
	       }
	    }
	    rhs[row_num] = rhs[i] - rows[i].fixed_lhs_offset;
	    row_num++;
	 }
      }
   }

   stats->rows_deleted += new_del_cnt;
   stats->vars_fixed += mip->mip_inf->fixed_var_num;

   /* now convert it row ordered if asked*/
   /* get row_lengths and fill in r_matbeg, sense, rhs etc */

   r_matbeg = mip->row_matbeg;
   r_matind = mip->row_matind;
   r_matval = mip->row_matval;
   r_lengths = mip->row_lengths;
   c_lengths = mip->col_lengths;

   for (i = 0; i < row_num; i++){
      r_lengths[i] = rows[i].size;
      r_matbeg[i+1] = r_matbeg[i] +
	 rows[i].size;
      rows[i].bin_var_num = 0;
   }

   /* debug */
   if (r_matbeg[row_num ] != col_nz){
      printf("error; missing nonzeros\n");
      return PREP_OTHER_ERROR;
   }

   /* */
   max_row_size = max_col_size = 0;

   for (i = 0; i < col_num; i++){
      for (j = matbeg[i]; j < matbeg[i+1]; j++){
	 matind[j] = row_new_inds[matind[j]];
	 row_ind = matind[j];
	 elem_ind = r_matbeg[row_ind];
	 r_matind[elem_ind] = i;
	 r_matval[elem_ind] = matval[j];
	 r_matbeg[row_ind] = elem_ind + 1;

	 if (cols[i].var_type == 'B'){
	    rows[row_ind].bin_var_num++;
	 }
      }
      c_lengths[i] = matbeg[i+1] - matbeg[i];
       if (max_col_size < c_lengths[i]){
	  max_col_size = c_lengths[i];
       }
   }

   bin_row_cnt = 0;
   cont_row_cnt = 0;
   bin_cont_row_cnt = 0;
   sos_row_cnt = 0;
   for (i = 0; i < row_num; i++){
      r_matbeg[i] -= r_lengths[i];
      if (max_row_size < r_lengths[i]){
	 max_row_size = r_lengths[i];
      }
      if (rows[i].bin_var_num){
	 bin_row_cnt++;
	 if (rows[i].cont_var_num){
	    bin_cont_row_cnt++;
	 }
      }

      if (rows[i].cont_var_num){
	 cont_row_cnt++;
      }

      if (rows[i].is_sos_row){
	 sos_row_cnt++;
      }
   }

   MIPinfo * mip_inf = mip->mip_inf;

   mip->n = col_num;
   mip->m = row_num;
   mip->nz = col_nz;
   mip->obj_offset = mip_inf->sum_obj_offset + obj_offset;
   mip->fixed_n = fixed_nz;
   mip_inf->binary_var_num = binary_var_num;
   mip_inf->binary_var_nz = binary_var_nz;
   mip_inf->max_row_size = max_row_size;
   mip_inf->max_col_size = max_col_size;
   mip_inf->binary_row_num = bin_row_cnt;
   mip_inf->cont_row_num = cont_row_cnt;
   mip_inf->bin_cont_row_num = bin_cont_row_cnt;
   mip_inf->binary_sos_row_num = sos_row_cnt;
   mip_inf->int_var_ratio = (1.0*(col_num - mip_inf->cont_var_num))/(col_num + 1);
   mip_inf->cont_var_ratio = 1.0*mip_inf->cont_var_num/(col_num + 1);
   mip_inf->bin_var_ratio = 1.0*binary_var_num/(col_num +1);
   mip_inf->max_row_ratio = 1.0*max_row_size/(col_num+1);
   mip_inf->max_col_ratio = 1.0*max_col_size/(row_num+1);
   mip_inf->mat_density = 1.0*col_nz/(col_num*row_num + 1);
   mip_inf->row_density = 1.0*col_nz/(row_num+1);
   mip_inf->col_density = 1.0*col_nz/(col_num+1);
   mip_inf->sos_bin_row_ratio = 1.0*sos_row_cnt/(bin_row_cnt +1);
   mip_inf->bin_row_ratio = 1.0*bin_row_cnt/(row_num +1);

   if (binary_var_num){
      mip_inf->row_bin_den = (int)
	 (binary_var_nz/bin_row_cnt) + 1;

      mip_inf->col_bin_den = (int)
	 (binary_var_nz/binary_var_num) + 1;

      if (binary_var_num < col_num){
      	 mip_inf->row_bin_den = (int)(mip_inf->row_bin_den*
					   col_num/
					   binary_var_num) + 1;
      }

      mip_inf->row_bin_den_mean =
	 (int)2*mip_inf->row_bin_den * max_row_size/
	 (mip_inf->row_bin_den + max_row_size) + 1;

      mip_inf->col_bin_den_mean =
	 (int)2*mip_inf->col_bin_den * max_row_size/
	 (mip_inf->col_bin_den + max_col_size) + 1;
   }

   FREE(row_new_inds);

   if (mip->n <= 0 || mip->m <= 0){
      return PREP_SOLVED;
   }

   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/
char prep_is_equal(double lval, double rval, double etol)
{

   double diff = lval - rval;

   if (diff < etol &&
      diff > -etol){
      return TRUE;
   } else {
      return FALSE;
   }
}
/*===========================================================================*/
/*===========================================================================*/

char prep_is_integral(double val, double etol)
{

   if (val - floor(val) < etol ||
      ceil(val) - val < etol){
      return TRUE;
   } else {
      return FALSE;
   }

}

/*===========================================================================*/
/*===========================================================================*/

int prep_report(PREPdesc *P, int termcode)
{

   MIPdesc *mip = P->mip;
   int i;
   prep_stats stats = P->stats;
   int p_level = P->params.level;
   char report_input  = FALSE;
   if (p_level > 2){
      switch(termcode){
       case PREP_INFEAS:
	 printf("Preprocessing detected infeasibility...");
	 if (stats.col_infeas_ind >= 0 ||
	    stats.row_infeas_ind >= 0){
	    printf("while improving bounds of \n\t");
	    if (stats.col_infeas_ind >= 0){
	       printf("variable ");
	       if (mip->colname){
		  printf("%s ", mip->colname[stats.col_infeas_ind]);
	       }
	       printf("[%i]", stats.col_infeas_ind);
	       if (stats.row_infeas_ind >= 0){
		  printf(" on the ");
	       }
	    }
	    if (stats.row_infeas_ind >= 0){
	       printf("row [%i]", stats.row_infeas_ind);
	    }
	    printf("\n");
	 }
	 break;
       case PREP_UNBOUNDED:
	 printf("Preprocessing detected unbounded problem...");
	 if (stats.col_unbound_ind >= 0){
	    printf("while improving bounds on \n");
	    if (mip->colname){
	       printf("variable %s [%i]\n",
		      mip->colname[stats.col_unbound_ind],
		      stats.col_unbound_ind);
	    } else {
	       printf("variable [%i]\n",
		      stats.col_unbound_ind);
	    }
	 }
	 break;
       case PREP_NUMERIC_ERROR:
	 printf("Preprocessing detected numerical problems ");
	 if (stats.col_numeric_ind >= 0){
	    printf("while improving bounds on \n");
	    if (mip->colname){
	       printf("variable %s [%i]\n",
		      mip->colname[stats.col_numeric_ind],
		      stats.col_numeric_ind);
	    } else {
	       printf("variable [%i]\n",
		      stats.col_numeric_ind);
	    }
	 }
	 break;
       case PREP_OTHER_ERROR:
	 printf("Preprocessing - unknown error.. ignoring presolve...\n");
	 break;
       case PREP_SOLVED:
	 printf("Preprocessing found the optimum:\n");
	 printf("Solution Cost: %.10f\n:",
		mip->obj_sense == SYM_MAXIMIZE ? -(mip->obj_offset) :
		mip->obj_offset);
	 if (mip->colname){
	    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	    printf("Column names and values of nonzeros in the solution\n");
	    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	    for (i = 0; i < mip->fixed_n; i++){
	       printf("%8s %10.10f\n", P->orig_mip->colname[mip->fixed_ind[i]],
		      mip->fixed_val[i]);
	    }
	    printf("\n");
	 } else {
	    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	    printf("User indices and values of nonzeros in the solution\n");
	    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	    for (i = 0; i < mip->fixed_n; i++){
	       printf("%7d %10.10f\n", mip->fixed_ind[i], mip->fixed_val[i]);
	    }
	    printf("\n");
	 }
	 break;
       default:
	 printf("Preprocessing finished...\n ");
	 report_input = TRUE;
	 if (stats.coeffs_changed +
	    stats.bounds_tightened +
	    stats.rows_deleted +
	    stats.vars_fixed +
	    stats.vars_aggregated +
	    stats.vars_integerized > 0){
	    if (stats.coeffs_changed > 0){
	       printf("\t coefficients modified: %i\n",
		      stats.coeffs_changed);
	    }
	    if (stats.bounds_tightened > 0){
	       printf("\t bounds improved: %i\n",
		      stats.bounds_tightened);
	    }
	    if (stats.rows_deleted +
	       stats.vars_fixed > 0){
	       if (stats.rows_deleted > 0){
		  printf("\t constraints removed: %i\n",
			 stats.rows_deleted);
		  //printf("\t %i remained\n", mip->m);
	       }
	       if (stats.vars_fixed > 0){
		  printf("\t variables fixed: %i\n", stats.vars_fixed);
		  //printf("\t %i remained\n", mip->n);
	       }
	    }
	    if (stats.vars_aggregated > 0){
	       printf("\t variables aggregated: %i\n",
		      stats.vars_aggregated);
	    }
	    if (stats.vars_integerized > 0){
	       printf("\t variables integerized: %i\n",
		      stats.vars_integerized);
	    }

	 } else {
	    printf("\t with no modifications...\n");
	 }
      }
   }else{
     report_input = TRUE;
   }

   if (report_input && P->params.verbosity >= 0){
     printf("Problem has \n"
	    "\t %i constraints \n"
	    "\t %i variables \n"
	    "\t %i nonzero coefficients\n",
	    mip->m, mip->n, mip->nz);
   }
   printf("\n");
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/
void prep_sos_fill_row(ROWinfo *row, int alloc_size, int size,
			   int *ind)
{
   /* check CHAR_BIT? */

   int i, sos_size = (alloc_size >> 3) + 1;
   if (row->sos_rep){
      memset(row->sos_rep, 0, CSIZE*sos_size);
   } else {
      row->sos_rep = (char *)calloc(CSIZE,sos_size);
   }

   for (i = 0; i < size; i++){
      (row->sos_rep[ind[i] >> 3]) |= (1 << (ind[i] & 7));
   }
}
/*===========================================================================*/
/*===========================================================================*/
void prep_sos_fill_var_cnt(PREPdesc *P)
{
   /* for now, just count the obvious variables those to be fixed
      if col_ind is fixed */

   ROWinfo * rows = P->mip->mip_inf->rows;
   COLinfo *cols = P->mip->mip_inf->cols;
   int n = P->mip->n;
   int m = P->mip->m;
   int sos_row_size = (n >> 3) + 1;

   int i, j, k;
   char * sos_final = (char *)malloc(CSIZE*sos_row_size);
   int sos_cnt = 0;

   int *matbeg = P->mip->matbeg;
   int *matind = P->mip->matind;
   int *r_matbeg = P->mip->row_matbeg;
   int *r_matind = P->mip->row_matind;
   int row_ind;

   for (i = 0; i < m; i++){
      if (rows[i].is_sos_row){
	 prep_sos_fill_row(&rows[i], n, r_matbeg[i+1] - r_matbeg[i],
			   &r_matind[i]);
      }
   }

   for (i = 0; i < n; i++){
      memset(sos_final, 0, CSIZE*sos_row_size);
      sos_cnt = 0;
      for (j = matbeg[i]; j < matbeg[i + 1]; j++){
	 row_ind = matind[j];
	 if (rows[row_ind].is_sos_row){
	    for (k = 0; k < sos_row_size; k++){
	       sos_final[k] |= rows[row_ind].sos_rep[k];
	    }
	 }
      }

      for (j = 0; j < sos_row_size; j++){
	 for (k = 7; k >= 0; k--){
	    sos_cnt += (sos_final[j] & (1 << k)) ? 1 : 0;
	 }
      }

      cols[i].sos_num = sos_cnt;
      //printf("col %i - sos_vars: %i\n", i, sos_cnt);
   }

   /* since we wont use these for now, delete sos row representations */
   for (i = 0; i < m; i++){
      if (rows[i].is_sos_row){
	 FREE(rows[i].sos_rep);
	 rows[i].sos_rep = 0;
      }
   }

   FREE(sos_final);

}
/*===========================================================================*/
/*===========================================================================*/
void free_prep_desc(PREPdesc *P)
{
   if (P){
      if (P->sr){
	 free_sr_desc(P->sr);
      }
      if (P->d_sr){
	 free_sr_desc(P->d_sr);
      }

      if (P->mip){
	 free_mip_desc(P->mip);
      }
      /* fixme - add impl stuff here - disabled now*/
      FREE(P->impl_vars);
      FREE(P->impl_ub);
      FREE(P->impl_lb);
      FREE(P->ulist_checked);
      FREE(P->llist_checked);
      FREE(P->rows_checked);

      /* since used to keep only static row and col info,
      */

      FREE(P->impl_cols);
      FREE(P->impl_rows);
      FREE(P->user_col_ind);
      FREE(P->user_row_ind);
      FREE(P->stats.nz_coeff_changed);
      FREE(P);
   }
}
/*===========================================================================*/
/*===========================================================================*/

void free_imp_list(IMPlist **list)
{

   IMPvar *imp_var;
   IMPvar *tmp_var;

   if (*list){
      for (imp_var = (*list)->head; imp_var != 0;){
	 tmp_var = imp_var->right;
	 FREE(imp_var);
	 imp_var = tmp_var;
      }
      FREE(*list);
      *list = 0;
   }
}

/*===========================================================================*/
/*===========================================================================*/
