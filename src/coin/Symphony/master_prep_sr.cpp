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
/* last modified: June 09, menal*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "sym_qsort.h"
#include "sym_macros.h"
#include "sym_constants.h"
#include "sym_prep.h"

/* Under Development */

/*===========================================================================*/
/*===========================================================================*/
void sr_initialize(SRdesc **sr, int n){

   int do_clean = FALSE;

   if(!(*sr)){
      *sr = (SRdesc *)calloc(1, sizeof(SRdesc));
      do_clean = TRUE;
   }

   if(!do_clean){
      (*sr)->prob_type = 0;
      (*sr)->max_n = (*sr)->min_n = 0;
      (*sr)->ub = (*sr)->lb = 0.0;
      (*sr)->ub_offset = (*sr)->lb_offset = 0.0;
      (*sr)->sum_a_max = (*sr)->sum_a_min = 0.0;
      (*sr)->sum_c_max = (*sr)->sum_c_min = 0.0;
      (*sr)->ub_updated = (*sr)->lb_updated = FALSE;
      (*sr)->rhs = (*sr)->rhs_max = (*sr)->rhs_min = 0.0;
      (*sr)->sense = ' ';
      if((*sr)->obj_max){
	 memset((*sr)->reversed_max, FALSE, CSIZE*n);
	 memset((*sr)->reversed_min, FALSE, CSIZE*n);
	 memset((*sr)->var_stat_max, SR_VAR_IN, ISIZE*n);
	 memset((*sr)->var_stat_min, SR_VAR_IN, ISIZE*n);
      }
   }
}



/*===========================================================================*/
/*===========================================================================*/

void sr_allocate(SRdesc **sr, int n){

   int k;
   (*sr)->obj_max = (double *)malloc(DSIZE*n);
   (*sr)->matval_max = (double *)malloc(DSIZE*n);
   (*sr)->matind_max = (int *)malloc(ISIZE*n);
   (*sr)->ratio_max = (double *)malloc(DSIZE*n);
   (*sr)->reversed_max = (char *)malloc(CSIZE*n);

   (*sr)->obj_min = (double *)malloc(DSIZE*n);
   (*sr)->matval_min = (double *)malloc(DSIZE*n);
   (*sr)->matind_min = (int *)malloc(ISIZE*n);
   (*sr)->ratio_min = (double *)malloc(DSIZE*n);
   (*sr)->reversed_min = (char *)malloc(CSIZE*n);

   /* for variable fixing, tightening etc... */

   (*sr)->var_max_opt = (double *)malloc(n* DSIZE);
   (*sr)->var_min_opt = (double *)malloc(n* DSIZE);
   (*sr)->var_stat_max = (int *)malloc(ISIZE*n);
   (*sr)->var_stat_min = (int *)malloc(n* ISIZE);
   (*sr)->var_obj_max = (double *)malloc(n* DSIZE);
   (*sr)->var_obj_min = (double *)malloc(n* DSIZE);
   (*sr)->var_matval_max = (double *)malloc(n* DSIZE);
   (*sr)->var_matval_min = (double *)malloc(n* DSIZE);

   /* debug, get something smart instead of these */
   (*sr)->tmp_ind = (int *)malloc(ISIZE*n);
   (*sr)->fixed_ind = (int *)malloc(ISIZE*n);

   for(k = 0; k < n; k++){
      (*sr)->fixed_ind[k] = k;
   }
}

/*===========================================================================*/
/*===========================================================================*/
int prep_solve_sr_rlx(PREPdesc *P, int row_cnt, int *row_indices)
{

   int i, j, k, l;
   int termcode = SR_NO_UPDATES;

   MIPdesc * mip = P->mip;
   prep_params params = P->params;
   MIPinfo *mip_inf = mip->mip_inf;

   COLinfo *cols = mip_inf->cols;
   ROWinfo *rows = mip_inf->rows;

   int n = mip->n, m = mip->m;
   int *c_matbeg = mip->matbeg;
   int *c_matind = mip->matind;

   int * r_matbeg = mip->row_matbeg;
   int * r_matind = mip->row_matind;
   double * r_matval = mip->row_matval;
   double *rhs = mip->rhs;
   char *sense = mip->sense;

   double *ub = mip->ub;
   double *lb = mip->lb;

   int max_sr_cnt, max_aggr_cnt, verbosity; //max_aggr_row_num, verbosity;
   int p_level, do_sr_rlx, do_aggr_row_rlx;
   int obj_ind, tot_sub_pr;
   char can_iterate = TRUE;
   double etol;

   do_sr_rlx = params.do_single_row_rlx;
   do_aggr_row_rlx = params.do_aggregate_row_rlx;
   p_level = params.level;
   etol = params.etol;
   verbosity = params.verbosity;

   /* get the max iteration limits */
   max_sr_cnt = params.max_sr_cnt;
   max_aggr_cnt = 0;

   /* initialize arrays to be used for each subproblem*/


   SRdesc * sr, *d_sr;

   if(!(P->rows_checked)){
      P->rows_checked = (char *)malloc(m* CSIZE);
   }

   char *rows_checked = P->rows_checked;
   double old_bound;
   //char no_upper, no_lower;
   int row_ind;// const_row_ind;
   //int updated_lb_cnt = 0, updated_ub_cnt = 0;
   int last_col_loc, last_row_loc;

   /*max_sr_cnt should be up to some ratio!!!! may be too many to handle*/
   tot_sub_pr = max_sr_cnt + max_aggr_cnt;

   /* do this only for all unbounded or all bounded rows for now...*/
   /* fixme... extend this if results happen to be good...*/

   for(i = 0; i < row_cnt; i++){

      obj_ind = row_indices[i];

      if(rows[obj_ind].bound_type == MIXED_BOUNDED_ROW ||
	 rows[obj_ind].is_redundant){
	 continue;
      }

      rows[obj_ind].orig_ub = rows[obj_ind].sr_ub = rows[obj_ind].ub;
      rows[obj_ind].orig_lb = rows[obj_ind].sr_lb = rows[obj_ind].lb;

      if(verbosity >=4){
	 printf("init bounds: row: %i", i);
	 printf("\told_lb:");
	 if(rows[obj_ind].sr_lb > -INF){
	 printf("%f", rows[obj_ind].sr_lb);
	 }else{
	    printf("-inf");
	 }
	 printf("\told_ub:");
	 if(rows[obj_ind].sr_ub < INF){
	    printf("%f", rows[obj_ind].sr_ub);
	 }else{
	    printf("inf");
	 }
	 printf("\n");
      }


      //     srows[i] = (SRrlx *)calloc(tot_sub_pr, sizeof(SRrlx));
      memset(rows_checked, FALSE, CSIZE*m);
      last_col_loc = r_matbeg[obj_ind];
      last_row_loc = c_matbeg[r_matind[last_col_loc]];

      for(j = 0; j < tot_sub_pr; j++){
	 /* first do single row stuff if can*/
	 //	 is_open_prob = 1;
	 row_ind = -1;
	 can_iterate = TRUE;
	 if(j < max_sr_cnt){
	    //int max_shared_row_ind = -1;
	    //int max_shared_size = 0;
	    //int row_search_iter_cnt = 0;
	    /* get something smarter here */
	    /*find a row that has the most common shared vars with this
	      one */
	    /*find a row to be used as a constraint*/
	    for(k = last_col_loc; k < r_matbeg[obj_ind+1]; k++){
	       for(l = last_row_loc; l < c_matbeg[r_matind[k]+1];
		   l++){
		  if(!rows[c_matind[l]].is_redundant &&
		     !rows_checked[c_matind[l]]){
		     rows_checked[c_matind[l]] = TRUE;
		     if(rows[obj_ind].bound_type ==
			rows[c_matind[l]].bound_type
			&& c_matind[l] != obj_ind) {
			row_ind = c_matind[l];
			break;
		     }
		  }
	       }

	       if(row_ind >= 0){
		  last_col_loc = k;
		  last_row_loc = l;
		  break;
	       }
	    }

	    if(row_ind >= 0){

	       sr_initialize(&(P->sr), n);
	       sr = P->sr;
	       sr->prob_type = rows[obj_ind].bound_type;
	       sr->rhs = rhs[row_ind];
	       sr->sense = sense[row_ind];

	       /* convert the problem to <= constraint*/
	       /* or solve it if it is unbounded */

	       switch(rows[obj_ind].bound_type){
		case OPEN_ROW: /* easiest case */

		   sr->rhs_max = sr->rhs_min = sr->rhs;

                   sr_solve_open_prob(P, sr, obj_ind, row_ind, r_matbeg,
				      r_matind, r_matval, cols, ub, lb, etol);

		   break;
		case ALL_BOUNDED_ROW:
		   if(rows[obj_ind].ub_inf_var_num +
		      rows[obj_ind].lb_inf_var_num +
		      rows[obj_ind].free_var_num > 0 ||
		      rows[row_ind].ub_inf_var_num +
		      rows[row_ind].lb_inf_var_num +
		      rows[row_ind].free_var_num > 0){
		      /* debug */
		      /* fixme, get rid of this bug*/
		      printf("something is wrong -case all_bounded_row-"
			     "prep_solve_sr_rlx(), exiting...\n");
		      return PREP_OTHER_ERROR;
		   }

		   /* always convert the problem into ax <= b ineq for max
		      ax >= b for min */

		   /* so,

		      _max arrays for solving [max cx st. ax <= b]
		      _min arrays for solving [min cx st. ax >= b]

		      if sense = E;

		      d_max arrays for solving [max cx st. -ax <= -b]
		      d_min arrays for solving [min cx st. -ax >= -b]

		   */

		   if(!sr->obj_max && rows[obj_ind].bound_type != OPEN_ROW){
		      sr_allocate(&sr, n);
		   }
		   switch(sr->sense){
		    case 'G':
		       sr->rhs_max = -sr->rhs;
		       sr->rhs_min = sr->rhs;
		       break;
		    case 'L':
		       sr->rhs_max = sr->rhs;
		       sr->rhs_min = -sr->rhs;
		       break;
		    case 'E':
		       sr->rhs_max = sr->rhs;
		       sr->rhs_min = -sr->rhs;

		       sr_initialize(&(P->d_sr), n);
		       d_sr = P->d_sr;
		       d_sr->prob_type = rows[obj_ind].bound_type;
		       d_sr->rhs = rhs[row_ind];
		       d_sr->sense = sense[row_ind];

		       d_sr->rhs_max = -d_sr->rhs;
		       d_sr->rhs_min = d_sr->rhs;

		       if(!d_sr->obj_max){
			  sr_allocate(&d_sr, n);
		       }

			break;
		   }

	           sr_solve_bounded_prob(P, sr, d_sr, obj_ind, row_ind,
					 r_matbeg, r_matind, r_matval,
					 cols, ub, lb, etol);
		   if(!rows[obj_ind].is_redundant){
		      if(sr->sense == 'E'){
			 if(sr->ub > d_sr->ub){
			    sr->ub = d_sr->ub;
			 }
			 if(sr->lb < d_sr->lb){
			    sr->lb = d_sr->lb;
			 }
		      }

		      sr->lb_updated = sr->ub_updated = TRUE;
		   }else{
		      break;
		   }
	       }

	       /* check for any progress */
	       if(sr->lb_updated){
		  if(rows[obj_ind].sr_lb < sr->lb){
		     old_bound = rows[obj_ind].sr_lb;
		     rows[obj_ind].sr_lb = sr->lb;
		     /* debug */
		     if(termcode != SR_BOUNDS_UPDATED){
			termcode = SR_BOUNDS_UPDATED;
		     }

		     if(verbosity >=5){
			printf("lb improved, "
			       "row: %i \told_lb:%f \tnew_lb:%f\n",
			       obj_ind, old_bound <= -INF ? 1 : old_bound, sr->lb);
		     }
		  }else if (rows[obj_ind].orig_lb > sr->lb + etol){
		     /* debug */
		     printf("error-lb, row: %i \told_lb:%f \tnew_lb:%f\n",
			    obj_ind, rows[obj_ind].orig_lb, sr->lb);
		  }
	       }
	       if(sr->ub_updated){
		  if(rows[obj_ind].sr_ub > sr->ub){
		     old_bound = rows[obj_ind].sr_ub;
		     rows[obj_ind].sr_ub = sr->ub;
		     /* debug */
		     if(termcode != SR_BOUNDS_UPDATED){
			termcode = SR_BOUNDS_UPDATED;
		     }
		     if(verbosity >=5){
			printf("ub improved, "
			       "row: %i \told_ub:%f \tnew_ub:%f\n",
			       obj_ind, old_bound >= INF ? -1 : old_bound, sr->ub);
		     }
		  }else if(rows[obj_ind].orig_ub < sr->ub - etol){
		     /*debug*/
		     //		     if(verbosity >=5){
		     printf("error-ub, row: %i \told_ub:%f \tnew_ub:%f\n",
			    obj_ind, rows[obj_ind].orig_ub, sr->ub);
			//		     }
		  }
		  if(sr->lb_updated){
		     if(sr->ub < sr->lb - etol){
			/* debug */
			printf("bounds err : "
			       "row: %i \tnew_ub:%f \tnew_lb:%f\n",
			       obj_ind, sr->ub, sr->lb);
			termcode = SR_INFEAS;
			break;
		     }
		  }
	       }

	    }
	 }
      }

      /* debug */
      if(termcode == SR_INFEAS){
	 break;
      }

      if(verbosity >=4){
	 printf("finl bounds: row: %i", i);
	 printf("\tnew_lb:");
	 if(rows[obj_ind].sr_lb > -INF){
	    printf("%f", rows[obj_ind].sr_lb);
	 }else{
	    printf("-inf");
      }
	 printf("\tnew_ub:");
	 if(rows[obj_ind].sr_ub < INF){
	 printf("%f", rows[obj_ind].sr_ub);
	 }else{
	    printf("inf");
	 }
	 printf("\n\n");

      }
   }

   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/

int sr_solve_bounded_prob(PREPdesc *P, SRdesc *sr, SRdesc *d_sr,
			  int obj_ind, int row_ind,
			  int *r_matbeg, int *r_matind, double *r_matval,
			  COLinfo *cols, double *ub, double *lb, double etol)
{

   int k, l, col_ind;
   double c_val, a_val;

   for( k = r_matbeg[obj_ind], l = r_matbeg[row_ind];;){
      if(k < r_matbeg[obj_ind + 1] &&
	 (r_matind[k] < r_matind[l] ||
	  l >= r_matbeg[row_ind + 1])){
	 c_val = r_matval[k];
	 col_ind = r_matind[k];
	 sr_add_new_col(sr, d_sr, c_val, 0.0, col_ind,
			cols[col_ind].var_type, ub[col_ind], lb[col_ind],
			sr->sense, 1, 1);
	 k++;
      }else if(l < r_matbeg[row_ind + 1] &&
	       (r_matind[k] > r_matind[l] ||
		k >= r_matbeg[obj_ind+1])){
	 a_val = r_matval[l];
	 col_ind = r_matind[l];

	 sr_add_new_col(sr, d_sr, 0.0, a_val, col_ind,
			cols[col_ind].var_type, ub[col_ind], lb[col_ind],
			sr->sense, 0, 1);
	 l++;
      }else{
	 /* now the indices are equal, fill in the arrays */
	 c_val = r_matval[k];
	 a_val = r_matval[l];
	 col_ind = r_matind[k];

	 if(c_val == 0.0 || a_val == 0.0){
	    printf("not nonzero???"
		   "numerical issues -case bounded row-"
		   "sr_solve_bounded_prob(), exiting...\n");
	    return PREP_OTHER_ERROR;
	 }

	 sr_add_new_col(sr, d_sr, c_val, a_val, col_ind,
			cols[col_ind].var_type, ub[col_ind], lb[col_ind],
			sr->sense, 2, 1);
	 k++;
	 l++;
      }
      if(k == r_matbeg[obj_ind + 1] && l == r_matbeg[row_ind + 1]){
	 break;
      }
   }

   /* now solve the problem */
   if(!P->mip->mip_inf->rows[obj_ind].is_redundant){
      sr_find_opt_bounded(P, sr, obj_ind, ub, lb);
   }

   if(!P->mip->mip_inf->rows[obj_ind].is_redundant){
      if(sr->sense == 'E'){
	 sr_find_opt_bounded(P, d_sr, obj_ind, ub, lb);
      }
   }


   int termcode = 0;
   ROWinfo *rows = P->mip->mip_inf->rows;
   double min_ub = sr->ub;
   double max_lb = sr->lb;

   if(sr->sense == 'E'){
      if(!P->mip->mip_inf->rows[obj_ind].is_redundant){
	 if(min_ub > d_sr->ub){
	    min_ub = d_sr->ub;
	 }

	 if(max_lb < d_sr->lb){
	    max_lb = d_sr->lb;
	 }
      }
   }
   if(rows[obj_ind].ub > min_ub || rows[obj_ind].lb < max_lb){
      termcode = prep_check_redundancy(P, obj_ind, TRUE, min_ub, max_lb,
				       FALSE, 0);
   }

   if(PREP_QUIT(termcode)){
      return termcode;
   }

   return(0);

}

/*===========================================================================*/
/*===========================================================================*/

int sr_find_opt_bounded(PREPdesc *P, SRdesc *sr, int obj_ind,
			 double *ub, double *lb)

{
   int i, last_ind, col_loc, col_ind, *var_stat; //,j, var_ind;
   char max_solved = FALSE, min_solved = FALSE;
   double lhs, ax, var_frac_val;
   /* get opt for each column (col fixed ub in min solved and
      lb in max solved - check also a_vals)*/
   double bound;

   int * tmp_ind = sr->tmp_ind;
   double etol = P->params.etol;

   if(sr->sum_a_max < sr->rhs_max +etol || sr->max_n <= 0){
      sr->ub += sr->sum_c_max + sr->ub_offset;
      max_solved = TRUE;
   }

   if(sr->sum_a_min > sr->rhs_min - etol|| sr->min_n <= 0){
      sr->lb += sr->sum_c_min + sr->lb_offset;
      min_solved = TRUE;
   }

   if(max_solved && min_solved){
      /* check redundancy */
      /* redundant and useless row*/
      return PREP_UNMODIFIED;
   }

   if(!max_solved){ /* otherwise, the row is redundant and useless */

      var_stat = sr->var_stat_max;
      memcpy(tmp_ind, sr->fixed_ind, ISIZE*sr->max_n);
      qsort_di(sr->ratio_max, tmp_ind, sr->max_n);
      //CoinSort_2(sr->ratio_max, sr->ratio_max + sr->max_n, tmp_ind);

      /* now fill in knapsack */
      lhs = 0;
      for(i = sr->max_n - 1; i >=0; i--){
	 col_loc = tmp_ind[i];
	 col_ind = sr->matind_max[col_loc];
	 bound = ub[col_ind] - lb[col_ind];

	 if(lhs > sr->rhs_max - etol){
	    break;
	 }

	 ax = sr->matval_max[col_loc] * bound;

	 if(lhs + ax < sr->rhs_max - etol){
	    sr->ub += bound *
	       sr->obj_max[col_loc];
	    lhs += ax;
	    var_stat[col_ind] = SR_VAR_IN_FIXED_UB;
	 }else{
	    var_frac_val = sr->obj_max[col_loc] *
	       (sr->rhs_max - lhs)/sr->matval_max[col_loc];
	    sr->ub += var_frac_val; //sr->obj_max[col_loc] * var_frac_val;
	    var_stat[col_ind] = SR_VAR_IN_FRAC;
	    last_ind = i;
	    break;
	 }
      }
      sr->ub += sr->ub_offset;
   }

   if(!min_solved){ /* otherwise this row is redundant and useless */
      memcpy(tmp_ind, sr->fixed_ind, ISIZE*sr->min_n);
      qsort_di(sr->ratio_min, tmp_ind, sr->min_n);
      //CoinSort_2(sr->ratio_min, sr->ratio_min + sr->min_n, tmp_ind);
      /* now fill in knapsack */
      lhs = 0;
      var_stat = sr->var_stat_min;
      for(i = 0; i < sr->min_n; i++){
	 col_loc = tmp_ind[i];
	 col_ind = sr->matind_min[col_loc];
	 bound = ub[col_ind] - lb[col_ind];
	 ax = sr->matval_min[col_loc] * bound;

	 if(lhs > sr->rhs_min - etol){
	    break;
	 }

	 if(lhs + ax < sr->rhs_min - etol){
	    sr->lb += bound *
	       sr->obj_min[col_loc];
	    lhs += ax;
	    var_stat[col_ind] = SR_VAR_IN_FIXED_UB;
	 }else{
	    //	    ax = (sr->rhs_max - lhs)/sr->matval_max[col_loc];
	    sr->lb += sr->obj_min[col_loc] *
	       (sr->rhs_min - lhs)/sr->matval_min[col_loc];
	    var_stat[col_ind] = SR_VAR_IN_FIXED_UB;
	    last_ind = i;
	    break;
	 }
      }
      sr->lb += sr->lb_offset;
   }

   return 0;

}

/*===========================================================================*/
/*===========================================================================*/

/* will add the column to problem if necessary */

int sr_add_new_col(SRdesc *sr, SRdesc *d_sr, double c_val, double a_val,
		   int col_ind, char var_type, double col_ub,
		   double col_lb, char sense,
		   int col_type, int col_bound_type)
{
   /* col_type =
      0 => c_val = 0, a_val != 0
      1 => c_val != 0, a_val = 0
      2 => c_val != 0, a_val != 0
   */

   /* col_bound_type =
      0 => open row
      1 => all bounded row
      2 => mixed bounded row
   */

   double rhs_ub_offset = a_val * col_ub;
   double rhs_lb_offset = a_val * col_lb;

   double obj_ub_offset = c_val * col_ub;
   double obj_lb_offset = c_val * col_lb;

   if(col_bound_type == 1){
      if(col_type >= 0){
	 if(var_type != 'F'){
	    switch(sense){
	     case 'L':
		sr_add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MAX, var_type);
		sr_add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MIN, var_type);
		break;
	     case 'G':
		sr_add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MAX, var_type);
		sr_add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MIN, var_type);
		break;
	     case 'E':
		sr_add_new_bounded_col(sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MAX, var_type);
		sr_add_new_bounded_col(sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MIN, var_type);
		sr_add_new_bounded_col(d_sr, c_val, -a_val, col_ind,
				    -rhs_ub_offset, -rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MAX, var_type);
		sr_add_new_bounded_col(d_sr, c_val, a_val, col_ind,
				    rhs_ub_offset, rhs_lb_offset,
				    obj_ub_offset, obj_lb_offset,
				    col_ub, col_lb,
				    SR_MIN, var_type);
		break;
	    }
	 }else{

	    sr->ub_offset += obj_ub_offset;
	    sr->lb_offset += obj_ub_offset;
	    sr->rhs_max -= rhs_ub_offset;
	    sr->rhs_min -= rhs_ub_offset;

	    if(sense == 'E'){
	       d_sr->ub_offset += obj_ub_offset;
	       d_sr->lb_offset += obj_ub_offset;
	       d_sr->rhs_max -= rhs_ub_offset;
	       d_sr->rhs_min -= rhs_ub_offset;
	    }
	 }
      }
   }
   return 0;
}

/*===========================================================================*/
/*===========================================================================*/

/* will add the column to problem if necessary, otherwise will update the
   offset values.
   will assume the sense is L for max and G for min, so
   the a_val has to be sent after being updated.
   For E, this function will be called twice each for max and min*/

int sr_add_new_bounded_col(SRdesc *sr, double c_val, double a_val,
			   int col_ind,
			   double rhs_ub_offset, double rhs_lb_offset,
			   double obj_ub_offset, double obj_lb_offset,
			   double col_ub, double col_lb, int obj_sense,
			   char var_type)
{
   /*
      ratio_type =
      0 c_val >0, a_val>0
      1 c_val >= 0, a_val <= 0
      2 c_val <= 0, a_val >= 0
      3 c_val < 0, a_val < 0
   */

   //  int n;// = sr->max_n, min_n = sr->min_n;

   /* we will convert the vars so that u-l >= x >= 0 */


   int ratio_type = 0;

   if(c_val > 0.0){
      if(a_val <= 0.0){
	 ratio_type = 1;
      }
   }else if(c_val < 0.0){
      if(a_val >= 0.0){
	 ratio_type = 2;
      }else{
	 ratio_type = 3;
      }
   }else{
      if(a_val <= 0.0){
	 ratio_type = 1;
      }else{
	 ratio_type = 2;
      }
   }

   int *n, *matind, *var_stat;
   double *obj, *matval, *rhs, *obj_offset, *sum, *obj_sum, *ratios;
   double *var_matval, *var_obj;
   char *is_reversed;
   if(obj_sense == SR_MAX){
      n = &(sr->max_n);
      obj_offset = &(sr->ub_offset);
      sum = &(sr->sum_a_max);
      obj_sum = &(sr->sum_c_max);
      rhs = &(sr->rhs_max);
      obj = sr->obj_max;
      matind = sr->matind_max;
      matval = sr->matval_max;
      ratios = sr->ratio_max;
      is_reversed = sr->reversed_max;
      var_stat = sr->var_stat_max;
      var_matval = sr->var_matval_max;
      var_obj = sr->var_obj_max;
      //var_lhs_offset = sr->opt_ub_var_offset;
   }else{
      n = &(sr->min_n);
      obj_offset = &(sr->lb_offset);
      sum = &(sr->sum_a_min);
      obj_sum = &(sr->sum_c_min);
      rhs = &(sr->rhs_min);
      obj = sr->obj_min;
      matind = sr->matind_min;
      matval = sr->matval_min;
      ratios = sr->ratio_min;
      is_reversed = sr->reversed_min;
      var_stat = sr->var_stat_min;
      var_matval = sr->var_matval_min;
      var_obj = sr->var_obj_min;
   }

#if 0
   sr->var_ub_lhs_offset[col_ind] = -rhs_ub_offset;
   sr->var_lb_lhs_offset[col_ind] = -rhs_lb_offset;
   sr->var_lb_obj_offset[col_ind] = obj_lb_offset;
   sr->var_ub_obj_offset[col_ind] = obj_ub_offset;
#endif

   if(ratio_type == 0){
      obj[*n] = c_val;
      matval[*n] = a_val;
      matind[*n] = col_ind;
      ratios[*n] = c_val/a_val;
      if(obj_sense == SR_MAX){
	 *sum += (rhs_ub_offset - rhs_lb_offset);
	 *obj_sum += (obj_ub_offset - obj_ub_offset);
      }else{
	 /* since bounds are converted to be
	    u - l >= x >= 0 */
	 *sum += 0.0;//rhs_lb_offset;
	 *obj_sum += 0.0;//obj_lb_offset;
      }
      (*n)++;
      /* to solve bounds problem */
      *rhs += -(rhs_lb_offset); /* conversion by x = y + l */
      *obj_offset += obj_lb_offset;

   }else if((ratio_type == 1 && obj_sense == SR_MAX) ||
	    (ratio_type == 2 && obj_sense == SR_MIN)){
      *rhs += -rhs_ub_offset;
      *obj_offset += obj_ub_offset;
      var_stat[col_ind] = SR_VAR_FIXED_UB;
      var_matval[col_ind] = a_val;
      var_obj[col_ind] = c_val;
   }else if((ratio_type == 1 && obj_sense == SR_MIN) ||
	    (ratio_type == 2 && obj_sense == SR_MAX)){
      *rhs += -rhs_lb_offset;
      *obj_offset += obj_lb_offset;
      var_stat[col_ind] = SR_VAR_FIXED_LB;
      var_matval[col_ind] = a_val;
      var_obj[col_ind] = c_val;
   }else{
      obj[*n] = -c_val;
      matval[*n] = -a_val;
      matind[*n] = col_ind;
      ratios[*n] = c_val/a_val;
      is_reversed[*n] = TRUE;
      if(obj_sense == SR_MAX){
	 *sum += -rhs_ub_offset + +rhs_lb_offset;
	 *obj_sum += -obj_ub_offset + rhs_lb_offset;
      }else{
	 *sum += 0.0; //-rhs_lb_offset;
	 *obj_sum += 0.0; //-obj_lb_offset;
      }
      (*n)++;
      /* to solve bounds problem */
      *rhs += -(rhs_ub_offset); /* conversion by x = -y + u */
      *obj_offset += obj_ub_offset;
   }

   return 0;
}
/*===========================================================================*/
/*===========================================================================*/

/* will modify the constraint to E and solve it */

int sr_solve_open_prob(PREPdesc *P, SRdesc *sr, int obj_ind,
		       int row_ind, int *r_matbeg,
		       int *r_matind, double *r_matval, COLinfo *cols,
		       double *ub, double *lb, double etol)
{

   int l, k, col_ind;

   double max_dual_ub = INF, min_dual_ub = INF;
   double max_dual_lb = -INF, min_dual_lb = -INF;
   double d_ratio, obj_val, a_val;

   char no_upper = FALSE, no_lower = FALSE, is_free_column;
   char is_fixed_column = FALSE;
   char can_iterate = TRUE, prob_infeasible = FALSE, is_null_obj;

   double *ub_offset = &(sr->ub_offset);
   double *lb_offset = &(sr->lb_offset);
   double rhs = sr->rhs;
   char sense = sr->sense;

   double obj_ub_offset;
   double obj_lb_offset;


   //  sr->prob_type = OPEN_PROB;

   for( k = r_matbeg[obj_ind], l = r_matbeg[row_ind];;){
      if(k < r_matbeg[obj_ind + 1] &&
	 (r_matind[k] < r_matind[l] ||
	  l >= r_matbeg[row_ind + 1])){
	 if(r_matval[k] > 0.0){
	    if(!no_upper){
	       if(ub[r_matind[k]] >= INF){
		  no_upper = TRUE;
	       }else{
		  *ub_offset += ub[r_matind[k]] * r_matval[k];
	       }
	    }
	    if(!no_lower){
	       if(lb[r_matind[k]] <= -INF){
		  no_lower = TRUE;
	       }else{
		  *lb_offset += lb[r_matind[k]] * r_matval[k];
	       }
	    }
	 }else if (r_matval[k] < 0.0){
	    if(!no_lower){
	       if(ub[r_matind[k]] >= INF){
		  no_lower = TRUE;
	       }else{
		  *lb_offset += ub[r_matind[k]] * r_matval[k];
	       }
	    }
	    if(!no_upper){
	       if(lb[r_matind[k]] <= -INF){
		  no_upper = TRUE;
	       }else{
		  *ub_offset += lb[r_matind[k]] * r_matval[k];
	       }
	    }
	 }
	 k++;
      }else{
	 if(l < r_matbeg[row_ind + 1] &&
	    (r_matind[k] > r_matind[l] ||
	     k >= r_matbeg[obj_ind+1])) {
	    is_null_obj = TRUE;
	    obj_val = 0.0;
	 }else{
	    is_null_obj = FALSE;
	    obj_val = r_matval[k];
	 }

	 /* first convert the column into 0 <= x <= inf */

	 a_val = r_matval[l];
	 col_ind = r_matind[l];
	 is_free_column = FALSE;
	 is_fixed_column = FALSE;
	 if(ub[col_ind] < INF && lb[col_ind] > -INF){
	    /* debug - get vars.type here */
	    if(ub[col_ind] > lb[col_ind] + etol){
		  /* debug */
	       printf("bounded column -case all open row-"
		      "sr_solve_open_prob(), exiting...\n");
	       return PREP_OTHER_ERROR;
	    }else{
	       /* fix column, take care of it here */
	       if(!is_null_obj){
		  obj_lb_offset = obj_val * lb[col_ind];
		  if(!no_upper){
		     *ub_offset += obj_lb_offset;
		  }
		  if(!no_lower){
		     *lb_offset += obj_lb_offset;
		  }
	       }
	       rhs += -(a_val *lb[col_ind]);
	       is_fixed_column = TRUE;
	    }
	 }else if(ub[col_ind] >= INF){
	    if(lb[col_ind] > -INF){
	       if(!is_null_obj){
		  obj_lb_offset = obj_val * lb[col_ind];
		  if(!no_upper){
		     *ub_offset += obj_lb_offset;
		  }
		  if(!no_lower){
		     *lb_offset += obj_lb_offset;
		  }
	       }
	       rhs += -(a_val * lb[col_ind]);
	    }else{
	       is_free_column = TRUE;
	    }
	 }else{
	    if(!is_null_obj){
	       obj_ub_offset = obj_val * ub[col_ind];
	       if(!no_upper){
		  *ub_offset += obj_ub_offset;
	       }
	       if(!no_lower){
		  *lb_offset += obj_ub_offset;
	       }
	    }
	    rhs += -(a_val * ub[col_ind]);
	    obj_val = -obj_val;
	    a_val = -a_val;
	 }

	 /* now get dual bounds */
	 if(!is_fixed_column){
	    if(a_val > 0.0 || a_val < 0.0){
	       d_ratio = obj_val/a_val;
	       if(a_val > 0.0){
		  if(d_ratio < min_dual_ub){
		     min_dual_ub = d_ratio;
		  }
		  if(-d_ratio < max_dual_ub){
		     max_dual_ub = -d_ratio;
		  }
		  if(is_free_column){
		     if(d_ratio > min_dual_lb){
			min_dual_lb = d_ratio;
		     }
		     if(-d_ratio > max_dual_lb){
			max_dual_lb = -d_ratio;
		     }
		  }
	       }else{
		  if(d_ratio > min_dual_lb){
		     min_dual_lb = d_ratio;
		  }
		  if(-d_ratio > max_dual_lb){
		     max_dual_lb = -d_ratio;
		  }

		  if(is_free_column){
		     if(d_ratio < min_dual_ub){
			min_dual_ub = d_ratio;
		     }
		     if(-d_ratio < max_dual_ub){
			max_dual_ub = -d_ratio;
		     }
		  }
	       }

	       if(min_dual_lb > min_dual_ub){
		  no_lower = TRUE;
	       }
	       if(max_dual_lb > max_dual_ub){
		  no_upper = TRUE;
		  /* debug */
		  //printf("unbounded or infeasible problem?"
		  // "-case all open row-"
		  // "prep_solve_sr_rlx(), exiting...\n");
		  //return PREP_OTHER_ERROR;
	       }
	    }else{
	       /* debug */
	       printf("not nonzero???"
		      "numerical issues -case all open row-"
		      "prep_solve_sr_rlx(), exiting...\n");
	       return PREP_OTHER_ERROR;
	    }
	 }
	 l++;
	 if(!is_null_obj){
	    k++;
	 }
      }

      if((no_upper && no_lower)){
	 can_iterate = FALSE;
	 break;
      }
      if(k == r_matbeg[obj_ind + 1] && l == r_matbeg[row_ind + 1]){
	 break;
      }
   }

   if(can_iterate){
      /* update the bounds for this row */

      switch(sense){
       case 'L':
	  if(max_dual_ub > 0.0){
	     max_dual_ub = 0.0;
	  }
	  if(min_dual_ub > 0.0){
	     min_dual_ub = 0.0;
	  }
	  break;
       case 'G':
	  if(max_dual_lb < 0.0){
	     max_dual_lb = 0.0;
	  }
	  if(min_dual_lb < 0.0){
	     min_dual_lb = 0.0;
	  }
	  break;
      }

      /* check again */
      //   if(min_dual_lb > min_dual_ub ||/
      //	 max_dual_lb > max_dual_ub){/
      //	 printf("unbounded or infeasible problem?"
      //		"-case all open row-"
      //		"prep_solve_sr_rlx(), exiting...\n");
      //	 return PREP_OTHER_ERROR;
      //  }

      if(!no_lower){
	 if(rhs >= 0){
	    if(min_dual_ub < INF){
	       sr->lb = min_dual_ub * rhs;
	    }else{
	       prob_infeasible = TRUE;
	    }
	 }else{
	    if(min_dual_lb > -INF){
	       sr->lb = min_dual_lb *rhs;
	    }else{
	       prob_infeasible = TRUE;
	    }
	 }
	 if(!prob_infeasible){
	    sr->lb += *lb_offset;
	    sr->lb_updated = TRUE;
	 }

      }

      if(!prob_infeasible){
	 if(!no_upper){
	    if(rhs >= 0){
	       if(max_dual_ub < INF){
		  sr->ub = -(max_dual_ub * rhs);
	       }else{
		  prob_infeasible = TRUE;
	       }
	    }else{
	       if(max_dual_lb > -INF){
		  sr->ub = -(max_dual_lb *rhs);
	       }else{
		  prob_infeasible = TRUE;
	       }
	    }
	    if(!prob_infeasible){
	       sr->ub += *ub_offset;
	       sr->ub_updated = TRUE;
	    }

	 }
      }
   }

   //   return(sr->lb_updated || sr->ub_updated);
   return(prob_infeasible);
}

/*===========================================================================*/
/*===========================================================================*/
void free_sr_desc(SRdesc *sr)
{
   if(sr){
      FREE(sr->obj_max);
      FREE(sr->matval_max);
      FREE(sr->matind_max);
      FREE(sr->ratio_max);

      FREE(sr->obj_min);
      FREE(sr->matval_min);
      FREE(sr->matind_min);
      FREE(sr->ratio_min);

      FREE(sr->fixed_ind);
      FREE(sr->tmp_ind);

      FREE(sr);
   }
}
/*===========================================================================*/
/*===========================================================================*/
