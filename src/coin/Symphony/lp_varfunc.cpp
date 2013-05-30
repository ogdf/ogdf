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

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <math.h>

#include "sym_lp.h"
#include "sym_proccomm.h"
#include "sym_types.h"
#include "sym_macros.h"
#include "sym_constants.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP functions related to column operations.
\*===========================================================================*/

void add_col_set(lp_prob *p, our_col_set *new_cols)
{
   LPdata *lp_data = p->lp_data;
   var_desc *evar, **extra, **vars = lp_data->vars;

   char *status = lp_data->status;

   int new_vars = new_cols->num_vars;
   int i, j, oldn;
   char *where_to_move;

   int cnt = 0;
   int *index;
   char *lu;
   double *bd;

   int to_lb_num, *to_lb_ind, to_ub_num, *to_ub_ind;

   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;

   colind_sort_extra(p);

   if (new_cols->dual_feas == NOT_TDF){
      to_lb_num = new_cols->rel_ub;
      to_lb_ind = new_cols->rel_ub_ind;
      to_ub_num = new_cols->rel_lb;
      to_ub_ind = new_cols->rel_lb_ind;
   }else{
      to_ub_num = new_cols->rel_ub;
      to_ub_ind = new_cols->rel_ub_ind;
      to_lb_num = new_cols->rel_lb;
      to_lb_ind = new_cols->rel_lb_ind;
   }

   if (new_vars){
      size_lp_arrays(lp_data, TRUE, FALSE, 0, new_vars, new_cols->nzcnt);
   }

   lu = lp_data->tmp.c; /* n (max(n, new_vars), but already resized) */
   index = lp_data->tmp.i1; /* n */
   bd = lp_data->tmp.d; /* 2*n (MAX(n, 2*new_vars), but already resized) */

   if (to_ub_num > 0){
      memset(lu, 'U', to_ub_num);
      /* Put the branching variables and base variable to the end
       * of the list and the extra variables to the beginning */
      for (i = to_ub_num - 1; i >= 0; i--){
	 j = to_ub_ind[i];
	 release_var(lp_data, j, MOVE_TO_UB); /* empty function for cplex */
	 status[j] = NOT_FIXED | (status[j] & NOT_REMOVABLE);
	 bd[cnt] = vars[j]->ub;
	 index[cnt++] = j;
      }
   }

   if (to_lb_num > 0){
      memset(lu + cnt, 'L', to_lb_num);
      for (i = to_lb_num - 1; i >= 0; i--){
	 j = to_lb_ind[i];
	 release_var(lp_data, j, MOVE_TO_LB); /* empty function for cplex */
	 status[j] = NOT_FIXED | (status[j] & NOT_REMOVABLE);
	 bd[cnt] = vars[j]->lb;
	 index[cnt++] = j;
      }
   }

   if (cnt > 0)
      change_bounds(lp_data, cnt, index, lu, bd);

   if (! new_vars)
      return;

   where_to_move = lp_data->tmp.c; /* new_vars */
   /* In the current implementation everything not in the matrix was at its
    * lower bound (0), therefore to restore dual feasibility they have to be
    * moved to their upper bounds, while when we just add them to the problem
    * they have to be moved to their lower bound */
   memset(where_to_move,
	  new_cols->dual_feas == NOT_TDF ? MOVE_TO_UB : MOVE_TO_LB, new_vars);

   add_cols(lp_data, new_vars, new_cols->nzcnt, new_cols->objx,
	    new_cols->matbeg, new_cols->matind, new_cols->matval,
	    new_cols->lb, new_cols->ub, where_to_move);
   lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
   lp_data->col_set_changed = TRUE;
   p->colset_changed = TRUE;
   lp_data->ordering = COLIND_ORDERED;

   /* update the extra parts of vars */
   oldn = lp_data->n - new_vars;
   extra = lp_data->vars + oldn;
   for (i = new_vars - 1; i >= 0; i--){
      evar = extra[i];
      evar->userind = new_cols->userind[i];
      evar->colind = oldn + i;
      evar->lb = new_cols->lb[i];
      evar->ub = new_cols->ub[i];
   }

   /* zero out x, i.e., set it to the LB */
   memset(lp_data->x + oldn, 0, new_vars * DSIZE);
   /* set status of the new vars to NOT_FIXED */
   //memset(lp_data->status + oldn, NOT_FIXED, new_vars);
   //TODO: char vs int
   for (i = oldn; i<oldn+new_vars; i++) {
      lp_data->status[i] = NOT_FIXED;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Try to tighten bounds based on reduced cost and logical fixing
\*===========================================================================*/

void tighten_bounds(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   double *dj = lp_data->dj;
   //double *x = lp_data->x;
   char *status = lp_data->status;
   var_desc **vars = lp_data->vars;
   int n = lp_data->n;
   double lpetol = lp_data->lpetol;

   double gap = 0.0, max_change;
   int i, vars_recently_fixed_to_ub = 0;
   int did_logical_fixing = FALSE,  did_reduced_cost_fixing = FALSE;
   int lb_vars = 0, perm_lb_vars = 0, ub_vars = 0, perm_ub_vars = 0;
   int del_vars = 0, *delstat = NULL;

   //char not_fixed__lb__switch, not_fixed__ub__switch;
   int *ind = 0;
   char *lu = 0;
   double *bd = 0, *ub = 0, *lb = 0;
   int cnt = 0;

   colind_sort_extra(p);

   check_ub(p);
   if (p->has_ub){
      gap = p->ub - lp_data->objval - p->par.granularity;
   }

   /*========================================================================*\
    *                   Here is the reduced cost fixing.
    *
    * If the gap is negative that means that we are above the limit, so don't
    * do anything.
    * Otherwise we do regular rc tightening if one of the following holds:
    * -- if we have done rc fixing before then the gap must have decreased
    *    significantly
    * -- if we haven't done rc tightening before, then the gap must be relatively
    * small compared to the upper bound
    \*=======================================================================*/
   if (p->par.do_reduced_cost_fixing && p->has_ub && gap > 0){
      if (p->last_gap == 0.0 ?
	  (gap < p->par.gap_as_ub_frac * p->ub) :
	  (gap < p->par.gap_as_last_gap_frac * p->last_gap)){
	 /* Tighten the upper/lower bounds on the variables,
	    prepare to delete them and do some statistics. */
	 delstat = lp_data->tmp.i1;   /* 2*n */
	 ind = lp_data->tmp.i1 + n;
	 lu = lp_data->tmp.c;         /* n */
	 bd = lp_data->tmp.d;         /* n */

	 get_bounds(lp_data);
	 ub = lp_data->ub;
	 lb = lp_data->lb;

	 p->vars_deletable = 0;
	 memset((char *)delstat, 0, n * ISIZE);
	 lb_vars = perm_lb_vars = ub_vars = perm_ub_vars = 0;
	 for (cnt = 0, i = n-1; i >= 0; i--){
	    if (fabs(dj[i]) < lpetol || !vars[i]->is_int){
	       continue;
	    }
	    max_change = gap/dj[i];
	    if (max_change > 0 && max_change < ub[i] - lb[i]){
	       if (lp_data->nf_status & NF_CHECK_NOTHING){
		  status[i] ^= NOT_FIXED__PERM_LB__SWITCH;
		  perm_lb_vars++;
	       }else{
		  status[i] ^= NOT_FIXED__TEMP_LB__SWITCH;
		  lb_vars++;
	       }
	       ind[cnt] = i;
	       lu[cnt] = 'U';
	       bd[cnt++] = vars[i]->is_int ? floor(lb[i] + max_change) :
		  lb[i] + max_change;
               vars[i]->new_ub = bd[cnt-1];
               p->bound_changes_in_iter++;
	       if (! (status[i] & NOT_REMOVABLE) && lb[i] == 0 &&
		   lb[i] == ub[i]){
		  p->vars_deletable++;
		  delstat[i] = 1;
	       }
	    }else if (max_change < 0 && max_change > lb[i] - ub[i]){
	       if (lp_data->nf_status & NF_CHECK_NOTHING){
		  status[i] ^= NOT_FIXED__PERM_UB__SWITCH;
		  perm_ub_vars++;
	       }else{
		  status[i] ^= NOT_FIXED__TEMP_UB__SWITCH;
		  ub_vars++;
	       }
	       ind[cnt] = i;
	       lu[cnt] = 'L';
	       bd[cnt++] = vars[i]->is_int ? ceil(ub[i] + max_change) :
		  ub[i] + max_change;
               vars[i]->new_lb = bd[cnt-1];
               p->bound_changes_in_iter++;
	       if (! (status[i] & NOT_REMOVABLE) && lb[i] == 0 &&
		   lb[i] == ub[i]){
		  p->vars_deletable++;
		  delstat[i] = 1;
	       }
	       vars_recently_fixed_to_ub++;
	    }
	    did_reduced_cost_fixing = TRUE;
	 }
	 p->vars_recently_fixed_to_ub += vars_recently_fixed_to_ub;
      }
   }

#ifdef COMPILE_IN_LP
   if (p->bc_level==0 && p->par.do_reduced_cost_fixing) {
      /* we are root node. we will save the reduced costs after each round of
       * cuts. whenever ub is updated, we can come back and update bounds in
       * the root
       */
      save_root_reduced_costs(p);
   }
#endif

   if (cnt > 0){
      change_bounds(lp_data, cnt, ind, lu, bd);
   }

   /*========================================================================*\
    * Logical fixing is done only if the number of variables recently fixed
    * to upper bound reaches a given constant AND is at least a certain
    * fraction of the total number of variables.
    \*=======================================================================*/

   if ((p->par.do_logical_fixing) &&
       (p->vars_recently_fixed_to_ub >
	p->par.fixed_to_ub_before_logical_fixing) &&
       (p->vars_recently_fixed_to_ub >
	p->par.fixed_to_ub_frac_before_logical_fixing * n)){
      logical_fixing_u(p);
      did_logical_fixing = TRUE;
   }

   if (! did_reduced_cost_fixing && ! did_logical_fixing)
      return;

   if (did_reduced_cost_fixing)
      p->last_gap = gap;
   if (did_logical_fixing)
      p->vars_recently_fixed_to_ub = 0;

   if (p->par.verbosity > 3){
      if (ub_vars)
	 printf("total of %i variables with temp adjusted UB ...\n",ub_vars);
      if (perm_ub_vars)
	 printf("total of %i variables with perm adjusted UB ...\n",perm_ub_vars);
      if (lb_vars)
	 printf("total of %i variables with temp adjusted LB ...\n",lb_vars);
      if (perm_lb_vars)
	 printf("total of %i variables with perm adjusted LB ...\n",perm_lb_vars);
   }
   p->vars_at_lb = lb_vars;
   p->vars_at_ub = ub_vars;

   /* if enough variables have been fixed, then physically compress the matrix,
    * eliminating the columns that are fixed to zero */
   if (p->vars_deletable > p->par.mat_col_compress_num &&
       p->vars_deletable > n * p->par.mat_col_compress_ratio){

      PRINT(p->par.verbosity,3, ("Compressing constraint matrix (col) ...\n"));
      del_vars = delete_cols(lp_data, p->vars_deletable, delstat);
      if (del_vars > 0){
	 lp_data->lp_is_modified = LP_HAS_BEEN_MODIFIED;
	 lp_data->col_set_changed = TRUE;
      }
      if (p->vars_deletable > del_vars){
	 PRINT(p->par.verbosity, 3,
	       ("%i vars were not removed because they were basic ...\n",
		p->vars_deletable - del_vars));
      }
      if (del_vars > 0){
	 p->vars_deletable -= del_vars;
	 PRINT(p->par.verbosity, 3,
	       ("%i vars successfully removed from the problem ...\n",
		del_vars));
	 for (i = p->base.varnum; i < n; i++){
	    if (delstat[i] != -1){
	       *(vars[delstat[i]]) = *(vars[i]);
	       vars[delstat[i]]->colind = delstat[i];
	    }
	 }
      }
   }
}

/*===========================================================================*\
 *===========================================================================*
 *
 * IMPORTANT:
 * No matter whether this routine was called with a primal feasible tableau
 * or not, if everything prices out, that means that if we were to add every
 * variable right now, we would still have a dual feasible tableau, therefore
 * the dual simplex could be continued, with the dual obj value only
 * increasing.
 *
 * NOTE: The tableau we know of is always dual feasible.
 *
 * This routine starts by collecting ALL (known and not known)
 * variables having reduced cost between 0 and 'gap' into
 * 'new_cols'. As soon as a non-dual-feasible variable is encountered,
 * the routine switches to collect the non-dual-feasibles into 'new_cols'.
 *
 * At the end, the result is:
 * -- if dual_feas is TDF_HAS_ALL, then new_cols contains the set of
 *    non-fixable variables not in the matrix.
 * -- if dual_feas is NOT_TDF, then new_cols contains all (well, at most
 *    max_ndf_vars) not dual feasible variables, which could be added.
 * -- if dual_feas is TDF_NOT_ALL then the whole tableau is dual feas, just
 *    we have ran out of space for storing the non-fixable variables
 *
\*===========================================================================*/

/*===========================================================================*\
 * THIS FUNCTION CAN BE CALLED ONLY IF
 *       p->lp_data->nf_status[0] != NF_CHECK_NOTHING !!!
\*===========================================================================*/

our_col_set *price_all_vars(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   double lpetol = lp_data->lpetol;
   int m = lp_data->m, n = lp_data->n;
   char *status = lp_data->status;
   double *dj = lp_data->dj;
   double *dual = lp_data->dualsol;

   int bvarnum = p->base.varnum;
   int extranum = lp_data->n - bvarnum;
   var_desc **vars = lp_data->vars;
   var_desc **extra = vars + bvarnum;

   int next_not_fixed, not_fixed_num = lp_data->not_fixed_num;
   int *not_fixed = lp_data->not_fixed;
   char new_nf_status = NF_CHECK_UNTIL_LAST;
   int  nf_status = lp_data->nf_status;
   int tmp_not_fixed_num = 0, *tmp_not_fixed, *itmp;

   int cutnum;
   cut_data **cuts;
   row_data *rows;

   char basevar = TRUE; /* just to keep gcc quiet */
   int ind, prevind, curind, nextind = -1; /* just to keep gcc quiet */
   double gap, red_cost;

   char must_add;
   int dual_feas;
   int termcode = p->lp_data->termcode;

   our_col_set *new_cols = (our_col_set *) calloc(1, sizeof(our_col_set));
   int max_ndf_vars, max_nfix_vars,  max_new_vars, max_new_nzcnt;
   int new_vars=0, new_nzcnt=0, rel_lb=0, rel_ub=0;

   int collen, *colind;
   double obj, lb, ub, *colval;

   int i, j, k;

#ifdef STATISTICS
   int nfix = 0;
#endif

   /* Compute how many non-dual-feasible cols we are willing to add.
    * Also compute how many non-fixable cols we are willing to add.
    * Then, to start with, set max_new_vars to be the latter. */
   max_ndf_vars = (int) (n * p->par.max_non_dual_feas_to_add_frac);
   if (max_ndf_vars < p->par.max_non_dual_feas_to_add_min)
      max_ndf_vars = p->par.max_non_dual_feas_to_add_min;
   if (max_ndf_vars > p->par.max_non_dual_feas_to_add_max)
      max_ndf_vars = p->par.max_non_dual_feas_to_add_max;

   if (termcode != LP_D_UNBOUNDED){
      max_nfix_vars = 0;
   }else{
      max_nfix_vars = (int) (n * p->par.max_not_fixable_to_add_frac);
      if (max_nfix_vars < p->par.max_not_fixable_to_add_min)
	 max_nfix_vars = p->par.max_not_fixable_to_add_min;
      if (max_nfix_vars > p->par.max_not_fixable_to_add_max)
	 max_nfix_vars = p->par.max_not_fixable_to_add_max;
   }

   tmp_not_fixed = (int *) malloc(p->par.not_fixed_storage_size * ISIZE);

   max_new_vars = MAX(max_nfix_vars, max_ndf_vars);
   max_new_nzcnt = m*max_new_vars;

   new_cols->rel_lb_ind = p->vars_at_lb ?
      (int *) malloc(p->vars_at_lb * ISIZE) : NULL;
   new_cols->rel_ub_ind = p->vars_at_ub ?
      (int *) malloc(p->vars_at_ub * ISIZE) : NULL;
   new_cols->objx = (double *) malloc(max_new_vars * DSIZE);
   new_cols->lb = (double *) malloc(max_new_vars * DSIZE);
   new_cols->ub = (double *) malloc(max_new_vars * DSIZE);
   new_cols->matbeg = (int *) malloc((max_new_vars+1) * ISIZE);
   new_cols->matbeg[0] = 0;
   new_cols->matind = (int *) malloc(max_new_nzcnt * ISIZE);
   new_cols->matval = (double *) malloc(max_new_nzcnt * DSIZE);
   new_cols->userind = (int *) malloc(max_new_vars * ISIZE);

   userind_sort_extra(p);

   /* Collect the non-base lpcuts */
   cutnum = m - p->base.cutnum;
   cuts = (cut_data **) lp_data->tmp.p1; /* m (actually, cutnum < m) */
   rows = lp_data->rows + p->base.cutnum;
   for (i = cutnum - 1; i >= 0; i--)
      cuts[i] = rows[i].cut;

   colind = lp_data->tmp.i1; /* m */
   colval = lp_data->tmp.d; /* 2*m */

   must_add = FALSE;
   dual_feas = TDF_HAS_ALL;
   check_ub(p);
   gap = p->has_ub ? p->ub - p->par.granularity - lp_data->objval :
      SYM_INFINITY;

   /*========================================================================*\
    * Now we loop through every single variable, and get those to be
    * added into the new_cols structure.
    *
    * In the loop 'i' runs on the base vars, 'j' on the extra vars and 'k'
    * on the not_fixed ones.
    * In every iteration we compute nextind, the next variable we know of,
    * based on i and j.
    * -- If we have run out of non-fixables and
    *   - nf_status == NF_CHECK_UNTIL_LAST then, we know of all non-fixable
    *     non-base variable and we are past them, then only base vars can be
    *     left to be checked. And nextind is processed.
    *   - nf_status == NF_CHECK_AFTER_LAST then there are more non-fixables,
    *     but we don't know what they are, then the user has to tell about the
    *     next variable.
    * -- If nextind is smaller than the next non-fixable we know of then
    * nextind is processed.
    * -- If it is greater then the next non-fixable has to be processed.
    \*=======================================================================*/

   curind = prevind = -1;

   for (i = 0, j = 0, k = 0; TRUE; prevind = curind){
      switch ((i < bvarnum ? 1 : 0) + (j < extranum ? 2 : 0)){
       case 0: /* none left */
	 nextind = -1; basevar = FALSE; break;
       case 1: /* only base vars left */
	 nextind = vars[i]->userind; basevar = TRUE; break;
       case 2: /* only extra vars left */
	 nextind = extra[j]->userind; basevar = FALSE; break;
       case 3: /* both of them left */
	 if (vars[i]->userind < extra[j]->userind){
	    nextind = vars[i]->userind; basevar = TRUE;
	 }else{
	    nextind = extra[j]->userind; basevar = FALSE;
	 }
	 break;
      }

      /*=====================================================================*\
       * If we have a chance to prove TDF   or
       * If we proved NOT_TDF but still have enough space to add new cols,
       * then the user the user generates the next col (or says that the next
       * col is what we have next), otherwise the next col is what we have
       * next.
      \*=====================================================================*/

      if ((dual_feas != NOT_TDF) ||
	  (dual_feas == NOT_TDF && new_vars < max_ndf_vars)){
	 if (k < not_fixed_num){
	    /* If anything is left on the not_fixed list then compare it to
	     * the next in matrix (nextind) and get the smaller one */
	    next_not_fixed = not_fixed[k];
	    if (nextind == -1 || nextind > next_not_fixed){
	       /* FIXME: Catch the error message for this function */
	       curind = generate_column_u(p, cutnum, cuts,
					  prevind, next_not_fixed,
					  GENERATE_NEXTIND,
					  colval, colind, &collen, &obj,
					  &lb, &ub);
	       k++;
	    }else{
	       if (nextind == next_not_fixed)
		  k++;
	       curind = nextind;
	    }
	 }else{
	    /*If nothing is left in not_fixed then things depend on nf_status*/
	    if (nf_status == NF_CHECK_UNTIL_LAST){
	       curind = nextind;
	    }else{ /* NF_CHECK_AFTER_LAST */
	       /* FIXME: Catch the error message for this function */
	       curind = generate_column_u(p, cutnum, cuts,
					  prevind, nextind,
					  GENERATE_REAL_NEXTIND,
					  colval, colind, &collen, &obj,
					  &lb, &ub);
	    }
	 }
      }else{
	 curind = nextind;
      }

      /* Now curind is the one we work with. If it's negative then there are
       * no more cols. */
      if (curind < 0)
	 break;

      if (curind == nextind){
	 /* The next col is what we have next, i.e., it is in the matrix.
	  * We've got to check it only if it is fixed to LB or UB. */
	 ind = basevar ? i : j + bvarnum;
	 red_cost = dj[ind];
	 if (status[ind] & TEMP_FIXED_TO_LB || status[ind] & TEMP_FIXED_TO_UB){
	    if (status[ind] & TEMP_FIXED_TO_LB){
	       if (red_cost < -lpetol){
		  if (dual_feas != NOT_TDF){
		     dual_feas = NOT_TDF;
		     rel_lb = rel_ub = new_vars = new_nzcnt = 0;
		  }
		  new_cols->rel_lb_ind[rel_lb++] = ind;
	       }else if (dual_feas != NOT_TDF){
		  if (red_cost < gap){
		     new_cols->rel_lb_ind[rel_lb++] = ind;
		  }else{
		     new_cols->rel_lb_ind[rel_lb++] = - ind - 1;
		  }
	       }
	    }else{ /* TEMP_FIXED_TO_UB */
	       if (red_cost > lpetol){
		  if (dual_feas != NOT_TDF){
		     dual_feas = NOT_TDF;
		     rel_lb = rel_ub = new_vars = new_nzcnt = 0;
		  }
		  new_cols->rel_ub_ind[rel_ub++] = ind;
	       }else if (dual_feas != NOT_TDF){
		  if (red_cost > -gap){
		     new_cols->rel_ub_ind[rel_ub++] = ind;
		  }else{
		     new_cols->rel_ub_ind[rel_ub++] = - ind - 1;
		  }
	       }
	    }
	 }
	 if (basevar)
	    i++;
	 else
	    j++;

      }else{ /* the next col is not in the matrix */
	 red_cost = obj - dot_product(colval, colind, collen, dual);
	 if (red_cost < -lpetol){
	    if (dual_feas != NOT_TDF){
	       dual_feas = NOT_TDF;
	       rel_lb = rel_ub = new_vars = new_nzcnt = 0;
	    }
	    must_add = TRUE;
	 }else if (dual_feas != NOT_TDF && red_cost < gap){
	    if (new_vars == max_nfix_vars){
	       /* Run out of space!! */
	       dual_feas = TDF_NOT_ALL;
#ifdef STATISTICS
	       if (!nfix)
		  nfix = new_vars;
	       nfix++;
#endif
	    }else{
	       must_add = TRUE;
	    }
	 }
	 if (must_add){
	    new_cols->objx[new_vars] = obj;
	    new_cols->lb[new_vars] = lb;
	    new_cols->ub[new_vars] = ub;
	    new_cols->matbeg[new_vars + 1] = new_cols->matbeg[new_vars]+collen;
	    memcpy((char *)(new_cols->matind+new_cols->matbeg[new_vars]),
		   (char *)colind, collen * ISIZE);
	    memcpy((char *)(new_cols->matval+new_cols->matbeg[new_vars]),
		   (char *)colval, collen * DSIZE);
	    new_nzcnt += collen;
	    new_cols->userind[new_vars++] = curind;
	    must_add = FALSE;
	 }
	 basevar = FALSE;
      }
      /* Add this variable to the not_fixed list if it cannot be permanently
       fixed*/
      if (red_cost > -gap && red_cost < gap && !basevar){
	 if (tmp_not_fixed_num < p->par.not_fixed_storage_size){
	    tmp_not_fixed[tmp_not_fixed_num++] = curind;
	 }else{
	    new_nf_status = NF_CHECK_AFTER_LAST;
	 }
      }
   }

   new_cols->num_vars = new_vars;
   new_cols->nzcnt = new_nzcnt;
   new_cols->rel_lb = rel_lb;
   new_cols->rel_ub = rel_ub;

   switch (new_cols->dual_feas = dual_feas){
    case NOT_TDF:
      PRINT(p->par.verbosity, 5, ("price_all_vars() : NOT_TDF.  [ %i ]\n",
				  new_vars + rel_ub + rel_lb));
      p->vars_at_lb -= rel_lb;
      p->vars_at_ub -= rel_ub;
      for (i = rel_lb - 1; i >= 0; i--){
	 if (! (status[new_cols->rel_lb_ind[i]] & NOT_REMOVABLE))
	    p->vars_deletable--;
      }
      break;

    case TDF_NOT_ALL:
      PRINT(p->par.verbosity, 5, ("price_all_vars() : TDF_NOT_ALL.\n"));
#ifdef STATISTICS
      PRINT(p->par.verbosity, 5, ("     ( nonfix / maxnonfix : %i / %i )\n",
				  nfix, max_nfix_vars));
#endif

      /*=====================================================================*\
       * If total dual feasibility has been proved but there are too
       * many vars to add, then (although we don't want to release
       * those already fixed but not permanently fixable) we can
       * permanently fix those which are marked such way.
      \*=====================================================================*/

      for (k = 0, i = 0; i < rel_lb; i++){
	 if ((j = new_cols->rel_lb_ind[i]) < 0){
	    j = -j-1;
	    status[j] ^= TEMP_PERM_LB__SWITCH;
	 }else{
	    new_cols->rel_lb_ind[k++] = j;
	 }
	 if (! (status[j] & NOT_REMOVABLE))
	    p->vars_deletable--;
      }
      new_cols->rel_lb = k;
      for (k = 0, i = 0; i < rel_ub; i++){
	 if ((j = new_cols->rel_ub_ind[i]) < 0){
	    status[-j-1] ^= TEMP_PERM_UB__SWITCH;
	 }else{
	    new_cols->rel_ub_ind[k++] = j;
	 }
      }
      new_cols->rel_ub = k;

      /*=====================================================================*\
       * If we come here, there should only be variables to add in the case
       * the LP is infeasible. Since we can't add all the variables in, we may
       * as well just wait for the call to restore_feasibility and only add in
       * the one (if there is one) that destroys the proof of infeasibility.
      \*=====================================================================*/

      new_vars = rel_ub = rel_lb = 0;

      /* Update which variables have to be priced out. */
      lp_data->nf_status = new_nf_status;
      lp_data->not_fixed_num = tmp_not_fixed_num;
      itmp = lp_data->not_fixed;
      lp_data->not_fixed = tmp_not_fixed;
      tmp_not_fixed = itmp;

      break;

    case TDF_HAS_ALL:
      /* Now if total dual feasibility is proved and there aren't too
       * many extra variables to be added... */
      lp_data->nf_status = NF_CHECK_NOTHING;
      lp_data->not_fixed_num = 0;
      for (k = 0, i = 0; i < rel_lb; i++){
	 if ((j = new_cols->rel_lb_ind[i]) < 0){
	    j = -j-1;
	    status[j] ^= TEMP_PERM_LB__SWITCH;
	 }else{
	    new_cols->rel_lb_ind[k++] = j;
	 }
	 if (! (status[j] & NOT_REMOVABLE))
	    p->vars_deletable--;
      }
      new_cols->rel_lb = rel_lb = k;
      for (k = 0, i = 0; i < rel_ub; i++){
	 if ((j = new_cols->rel_ub_ind[i]) < 0){
	    status[-j-1] ^= TEMP_PERM_UB__SWITCH;
	 }else{
	    new_cols->rel_ub_ind[k++] = j;
	 }
      }
      new_cols->rel_ub = rel_ub = k;

      p->vars_at_lb = 0; /* they are either permanently fixed or released */
      p->vars_at_ub = 0; /* they are either permanently fixed or released */

      PRINT(p->par.verbosity, 5, ("price_all_vars() : TDF_HAS_ALL.  [ %i ]\n",
				  new_vars + rel_ub + rel_lb));
      break;
   }

   FREE(tmp_not_fixed);

   if (rel_lb || rel_ub)
      PRINT(p->par.verbosity, 1,
	    ("Released %i 0-variables and %i 1-variables.\n", rel_lb, rel_ub));

   if (new_vars || rel_lb || rel_ub)
      add_col_set(p, new_cols);

   return(new_cols);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function is called with an infeasible problem in lp_data, and after
 * we found that the lp is total dual feasible but there were too many
 * variables to add (i.e., TDF_NOT_ALL).
 * Therefore in new_cols there will be a bunch of vars wich are not permanently
 * fixable so we should start to check those first.
 *
 * The final goal here is to find a variable which destroys a proof of primal
 * infeasibility
\*===========================================================================*/

int restore_lp_feasibility(lp_prob *p, our_col_set *new_cols)
{
   LPdata *lp_data = p->lp_data;
   double lpetol = lp_data->lpetol;
   char *status = lp_data->status;
   double *dual = lp_data->dualsol;

   int bvarnum = p->base.varnum;
   int extranum = lp_data->n - bvarnum;
   var_desc **vars = lp_data->vars;
   var_desc **extra = vars + bvarnum;

   int next_not_fixed, *not_fixed = lp_data->not_fixed;
   int  nf_status = lp_data->nf_status;
   int not_fixed_num = lp_data->not_fixed_num;

   int cutnum;
   cut_data **cuts;

   int infind, violation;

   int collen, *colind;
   double obj, lb = 0, ub = 0, *colval, *binvrow;

   double gap, red_cost, prod;

   char basevar = TRUE; /* just to keep gcc quiet */
   int prevind, curind, nextind = -1; /* just to keep gcc quiet */

   int i, j, k;

   /* Get a proof of infeasibility and get the corresponding row of the basis
    * inverse */

   violation = get_proof_of_infeas(lp_data, &infind);

   /*========================================================================*\
    * Collect the non-base lpcuts. We would have to do the same as in
    * price_all_vars(), but this function is called right after that, and thus
    * lp_data->tmp.p1 still has the pointers to the cuts :-).
    * And, price_all_vars did NOT resize the tmp arrays as it has failed if
    * we came to this function.
    * Also itmpm and dtmpm were reallocated there. (for big enough)
   \*========================================================================*/

   cutnum = lp_data->m - p->base.cutnum;
   cuts = (cut_data **) lp_data->tmp.p1;
   colind = lp_data->tmp.i1;
   colval = lp_data->tmp.d;
   binvrow = lp_data->tmp.d + lp_data->m;

   get_binvrow(lp_data, infind, binvrow);

   check_ub(p);
   gap = p->has_ub ? p->ub - p->par.granularity - lp_data->objval :
      SYM_INFINITY;

   /* First check those released from their lower bound in price_all_vars(),
    * and see if they destroy the proof of infeas */
   for (i=new_cols->rel_lb-1; i>=0; i--){
      j = new_cols->rel_lb_ind[i];
      get_column(lp_data, j, colval, colind, &collen, &obj);
      prod = dot_product(colval, colind, collen, binvrow);
      if ((violation == LOWER_THAN_LB && prod < -lpetol) ||
	  (violation == HIGHER_THAN_UB && prod > lpetol)){
	 /* OK, just release this variable */
	 PRINT(p->par.verbosity, 2,
	       ("RELEASED_FROM_LB while restoring feasibility.\n"));
	 new_cols->num_vars = new_cols->rel_ub = new_cols->rel_lb = 0;
	 change_ub(lp_data, j, lp_data->vars[j]->ub);
	 status[j] ^= NOT_FIXED__TEMP_LB__SWITCH;
	 release_var(lp_data, j, MOVE_TO_LB);
	 return(TRUE);
      }
   }
   new_cols->rel_lb = 0; /*We don't need these anymore*/

   /* Now check those released from their upper bound */
   for (i=new_cols->rel_ub-1; i>=0; i--){
      j = new_cols->rel_ub_ind[i];
      get_column(lp_data, j, colval, colind, &collen, &obj);
      prod = dot_product(colval, colind, collen, binvrow);
      if ((violation == LOWER_THAN_LB && prod > lpetol) ||
	  (violation == HIGHER_THAN_UB && prod < -lpetol)){
	 /* OK, just release this variable */
	 PRINT(p->par.verbosity, 2,
	       ("RELEASED_FROM_UB while restoring feasibility.\n"));
	 new_cols->num_vars = new_cols->rel_ub = new_cols->rel_lb = 0;
	 change_lb(lp_data, j, lp_data->vars[j]->lb);
	 status[j] ^= NOT_FIXED__TEMP_UB__SWITCH;
	 release_var(lp_data, j, MOVE_TO_UB);
	 return(TRUE);
      }
   }
   new_cols->rel_ub = 0; /*We don't need these anymore*/

   /* Now check the ones described in the new_vars part of new_cols.
    * These are either already added, or we got that far */
   for (i=0; i<new_cols->num_vars; i++){
      prod = dot_product(new_cols->matval + new_cols->matbeg[i],
			 new_cols->matind + new_cols->matbeg[i],
			 new_cols->matbeg[i+1] - new_cols->matbeg[i],
			 binvrow);
      if ((violation == LOWER_THAN_LB && prod < -lpetol) ||
	  (violation == HIGHER_THAN_UB && prod > lpetol)){
	 PRINT(p->par.verbosity, 2,
	       ("1 variable added while restoring feasibility.\n"));
	 new_cols->rel_ub = new_cols->rel_lb = 0;
	 new_cols->num_vars = 1;
	 if (i > 0){
	    new_cols->userind[0] = new_cols->userind[i];
	    new_cols->objx[0] = new_cols->objx[i];
	    new_cols->lb[0] = lb;
	    new_cols->ub[0] = ub;
	    memmove(new_cols->matind, new_cols->matind + new_cols->matbeg[i],
		    new_cols->nzcnt * ISIZE);
	    memmove(new_cols->matval, new_cols->matval + new_cols->matbeg[i],
		    new_cols->nzcnt * DSIZE);
	    new_cols->matbeg[1] = new_cols->nzcnt;
	 }
	 new_cols->nzcnt = new_cols->matbeg[i+1] - new_cols->matbeg[i];
	 add_col_set(p, new_cols);
	 return(TRUE);
      }
   }

   /* OK, we are out of the if, so we still didn't get rid of the proof. */

   userind_sort_extra(p);

   /* Just to avoid copying */
   colind = new_cols->matind;
   colval = new_cols->matval;

   /*========================================================================*\
    * Go through all the columns not in the matrix starting from prevind,
    * i.e., we start where price_all_vars gave up collecting not fixables.
    * Those in the matrix are already tested; they were listed in
    * rel_{lb,ub}_ind.
    * To do this first get the right i,j,k for a loop awfully similar to the
    * one in price_all_vars.
   \*========================================================================*/

   prevind = new_cols->userind[new_cols->num_vars-1];
   i = bvarnum > 0 ? bfind(prevind, p->base.userind, bvarnum) + 1 : 0;
   for (j = extranum - 1; j >= 0 && extra[j]->userind > prevind; j--); j++;
   k = not_fixed_num > 0 ? bfind(prevind, not_fixed, not_fixed_num) + 1 : 0;

   for (; ; prevind = curind){
      if (k == not_fixed_num && nf_status != NF_CHECK_AFTER_LAST)
	 /* nothing can help now... */
	 break;
      switch ((i < bvarnum ? 1 : 0) + (j < extranum ? 2 : 0)){
       case 0: /* none left */
	 nextind = -1; break;
       case 1: /* only base vars left */
	 nextind = vars[i]->userind; basevar = TRUE; break;
       case 2: /* only extra vars left */
	 nextind = extra[j]->userind; basevar = FALSE; break;
       case 3: /* both of them left */
	 if (vars[i]->userind < extra[j]->userind){
	    nextind = vars[i]->userind; basevar = TRUE;
	 }else{
	    nextind = extra[j]->userind; basevar = FALSE;
	 }
	 break;
      }
      if (k < not_fixed_num){
	 next_not_fixed = not_fixed[k];
	 if (nextind == -1 || nextind > next_not_fixed){
	    /* FIXME: Catch the error message for this function */
	    curind = generate_column_u(p, cutnum, cuts,
				       prevind, next_not_fixed,
				       GENERATE_NEXTIND, colval, colind,
				       &collen, &obj, &lb, &ub);
	    k++;
	 }else{
	    if (nextind == next_not_fixed)
	       k++;
	    curind = nextind;
	 }
      }else{ /* no we know that NF_CHECK_AFTER_LAST */
	 /* FIXME: Catch the error message for this function */
	 curind = generate_column_u(p, cutnum, cuts,
				    prevind, nextind, GENERATE_REAL_NEXTIND,
				    colval, colind, &collen, &obj, &lb, &ub);
      }

      /* Now curind is the one we work with. If it's negative then there are
       * no more cols. */
      if (curind < 0)
	 break;

      if (curind == nextind){
	 /* no point in testing curind: it is either in the matrix or has
	  * already been tested */
	 if (basevar)
	    i++;
	 else
	    j++;
      }else{
	 prod = dot_product(colval, colind, collen, binvrow);
	 if ((violation == LOWER_THAN_LB && prod < -lpetol) ||
	     (violation == HIGHER_THAN_UB && prod > lpetol)){
	    red_cost = obj - dot_product(colval, colind, collen, dual);
	    if (red_cost < gap){ /* It is at 0 level anyway and dual feas*/
	       PRINT(p->par.verbosity, 2,
		     ("1 variable added while restoring feasibility.\n"));
	       new_cols->num_vars = 1;
	       new_cols->userind[0] = curind;
	       new_cols->objx[0] = obj;
	       new_cols->matbeg[1] = collen;
	       new_cols->nzcnt = collen;
	       add_col_set(p, new_cols);
	       return(TRUE);
	    }
	 }
      }
   }

   /* We came out of the loop ==> primal feasibility cannot be restored */
   return(FALSE);
}

/*===========================================================================*/

void userind_sort_extra(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   if (lp_data->n > p->base.varnum + 1){
      if (lp_data->ordering == COLIND_ORDERED){
	 qsort((char *)(lp_data->vars + p->base.varnum),
	       lp_data->n - p->base.varnum,
	       sizeof(var_desc *), var_uind_comp);
	 lp_data->ordering = USERIND_ORDERED;
      }
   }else{
      lp_data->ordering = COLIND_AND_USERIND_ORDERED;
   }
}

/*===========================================================================*/

void colind_sort_extra(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   if (lp_data->n > p->base.varnum + 1){
      if (lp_data->ordering == USERIND_ORDERED){
	 qsort((char *)(lp_data->vars + p->base.varnum),
	       lp_data->n - p->base.varnum,
	       sizeof(var_desc *), var_cind_comp);
	 lp_data->ordering = COLIND_ORDERED;
      }
   }else{
      lp_data->ordering = COLIND_AND_USERIND_ORDERED;
   }
}

/*===========================================================================*/

int var_uind_comp(const void *v0, const void *v1)
{
   return((*(var_desc **)v0)->userind - (*(var_desc **)v1)->userind);
}

/*===========================================================================*/

int var_cind_comp(const void *v0, const void *v1)
{
   return((*(var_desc **)v0)->colind - (*(var_desc **)v1)->colind);
}


/*===========================================================================*/
#ifdef COMPILE_IN_LP
int save_root_reduced_costs(lp_prob *p)
{
   int         *indices;
   double      *values, *lb, *ub;
   tm_prob     *tm      = p->tm;
   rc_desc     *rc = NULL;
   int          pos;
   int         *tind = p->lp_data->tmp.i1;
   int          cnt = 0, i, j;
   int          n = p->lp_data->n;
   var_desc   **vars = p->lp_data->vars;
   double       lpetol = p->lp_data->lpetol;
   double      *lp_lb, *lp_ub;
   double      *dj = p->lp_data->dj;

   get_bounds(p->lp_data);
   lp_lb = p->lp_data->lb;
   lp_ub = p->lp_data->ub;
   for (i = 0; i < n; i++){
      if (vars[i]->is_int && lp_ub[i]-lp_lb[i]>lpetol &&
            (dj[i] > lpetol || dj[i] < -lpetol)){
         tind[cnt] = i;
         cnt++;
      }
   }
   PRINT(p->par.verbosity, 5, ("there are %d non zero reduced costs for "
            "integer vars\n", cnt));

   if (cnt==0) {
      return 0;
   }
   indices = (int *)malloc(cnt*ISIZE);
   values = (double *)malloc(cnt*DSIZE);
   lb = (double *)malloc(cnt*DSIZE);
   ub = (double *)malloc(cnt*DSIZE);

   for (i = 0; i < cnt; i++){
      j = tind[i];
      indices[i] = vars[j]->userind;
      values[i] = dj[j];
      lb[i] = lp_lb[j];
      ub[i] = lp_ub[j];
   }
   /*
   for (int i=0;i<cnt;i++) {
      printf("var %d, %20.10f\n",indices[i],values[i]);
   }
   printf("\n\n");
   */

   if (!tm->reduced_costs) {
      tm->reduced_costs = (rc_desc *) malloc(sizeof(rc_desc));
      rc = tm->reduced_costs;
      rc->size    = 10;
      rc->num_rcs = 0;
      rc->indices = (int **)calloc(rc->size,sizeof(int *));
      rc->values  = (double **)calloc(rc->size,sizeof(double *));
      rc->lb      = (double **)calloc(rc->size,sizeof(double *));
      rc->ub      = (double **)calloc(rc->size,sizeof(double *));
      rc->obj     = (double *)malloc(rc->size*DSIZE);
      rc->cnt     = (int *)calloc(rc->size,ISIZE);
   } else {
      rc = tm->reduced_costs;
   }

   pos = rc->num_rcs%rc->size;
   if (rc->size==rc->num_rcs) {
      /* replace the oldest one with the new one */
      FREE(rc->indices[pos]);
      FREE(rc->values[pos]);
      FREE(rc->lb[pos]);
      FREE(rc->ub[pos]);
   }
   rc->indices[pos] = indices;
   rc->values[pos] = values;
   rc->lb[pos] = lb;
   rc->ub[pos] = ub;
   rc->cnt[pos] = cnt;
   rc->obj[pos] = p->lp_data->objval;
   if (rc->num_rcs < rc->size) {
      rc->num_rcs++;
   }
   return 0;
}

/*===========================================================================*/
int tighten_root_bounds(lp_prob *p)
{
   /*
    * using the reduced costs that are saved from the root node, try to
    * improve variable bounds.
    * should be called whenever ub is updated.
    * change only bounds for root. not for the current node. the bounds for
    * current node are updated in tighten_bounds()
    */
   int                  i, j, k, l;
   rc_desc             *rc = p->tm->reduced_costs;
   double               gap, max_change;
   double              *dj, *lb, *ub;
   int                 *saved_ind;
   int                 *ind = p->lp_data->tmp.i1;
   double              *bd = p->lp_data->tmp.d;
   char                *lu = p->lp_data->tmp.c;
   int                  cnt, total_changes = 0;
   double               lpetol = p->lp_data->lpetol;
   bounds_change_desc  *bnd_change;
   int                 *new_ind;
   int                  num_new_bounds;
   int                  verbosity = p->par.verbosity;
   int                 *oldindex;
   double              *oldvalue;
   char                *oldlu;

   if (!rc) {
      return 0;
   }

   if (!p->has_ub) {
      PRINT(verbosity, -1, ("tighten_root_bounds: cant tighten bounds if ub "
            "does not exist!\n"));
      return 0;
   }

   new_ind = (int *)malloc(p->mip->n*ISIZE);
   for (i=0; i<rc->num_rcs;i++) {
      gap = p->ub - rc->obj[i] - p->par.granularity;
      if (gap<=lpetol) {
         continue;
      }
      saved_ind = rc->indices[i];
      dj  = rc->values[i];
      lb = rc->lb[i];
      ub = rc->ub[i];
      cnt = 0;
      for (j=0; j<rc->cnt[i]; j++) {
         max_change = gap/dj[j];
         if (max_change > 0 && max_change < ub[j]-lb[j]){
            ind[cnt] = saved_ind[j];
            lu[cnt] = 'U';
            bd[cnt++] = floor(lb[j] + max_change);
         }else if (max_change < 0 && max_change > lb[j] - ub[j]){
            ind[cnt] = saved_ind[j];
            lu[cnt] = 'L';
            bd[cnt++] = ceil(ub[j] + max_change);
         }
      }
      PRINT(verbosity, 5, ("tighten_root_bounds: at node %d, tightening %d "
               "bounds in root\n", p->bc_index, cnt));
      if (cnt == 0) {
         continue;
      }
      /* add these changes to root node */
      if (p->tm->rootnode->desc.bnd_change) {
         bnd_change = p->tm->rootnode->desc.bnd_change;
      } else {
         p->tm->rootnode->desc.bnd_change = bnd_change =
            (bounds_change_desc *)calloc(1,sizeof(bounds_change_desc));
      }
      if (bnd_change->num_changes>0) {
         /*
          * update existing changes and store the new ones in a separate array
          */
         num_new_bounds=0;
         oldvalue = bnd_change->value;
         oldindex = bnd_change->index;
         oldlu    = bnd_change->lbub;
         for (k=0; k<cnt; k++) {
            for (j=0; j<bnd_change->num_changes; j++) {
               if (oldindex[j]==ind[k] && oldlu[j]==lu[k]){
                  if (lu[k]=='L' && oldvalue[j]<bd[k]) {
                     oldvalue[j]=bd[k];
                     total_changes++;
                  } else if (lu[k]=='U' && oldvalue[j]>bd[k]) {
                     oldvalue[j]=bd[k];
                     total_changes++;
                  }
                  break;
               }
            }
            if (j>=bnd_change->num_changes) {
               new_ind[num_new_bounds] = k;
               num_new_bounds++;
            }
         }
         /* those changes that dint already have an entry and stored now */
         if (num_new_bounds) {
            int new_cnt = num_new_bounds+bnd_change->num_changes;
            bnd_change->index = (int *)realloc(bnd_change->index,
                  ISIZE*new_cnt);
            bnd_change->lbub  = (char *)realloc(bnd_change->lbub,
                  CSIZE*new_cnt);
            bnd_change->value = (double *)realloc(bnd_change->value,
                  DSIZE*new_cnt);
            oldvalue = bnd_change->value;
            oldindex = bnd_change->index;
            oldlu    = bnd_change->lbub;
            l = bnd_change->num_changes;
            for (j=0; j<num_new_bounds; j++) {
               total_changes++;
               k = new_ind[j];
               oldindex[l] = ind[k];
               oldlu[l]    = lu[k];
               oldvalue[l] = bd[k];
               bnd_change->num_changes++;
               l++;
            }
         }
      } else {
         bnd_change->index = (int *)malloc(cnt*ISIZE);
         bnd_change->lbub  = (char *)malloc(cnt*CSIZE);
         bnd_change->value = (double *)malloc(cnt*DSIZE);
         bnd_change->index = (int *) memcpy(bnd_change->index, ind, ISIZE*cnt);
         bnd_change->lbub  = (char *) memcpy(bnd_change->lbub, lu, CSIZE*cnt);
         bnd_change->value = (double *) memcpy(bnd_change->value, bd,
               DSIZE*cnt);
         bnd_change->num_changes = cnt;
      }
   }
   if (verbosity>5 && p->tm->rootnode->desc.bnd_change!=NULL) {
      printf("tighten_root_bounds: root now has %d changes\n",
            p->tm->rootnode->desc.bnd_change->num_changes);
   }
   FREE(new_ind);
   return 0;
}
#endif
/*===========================================================================*/
/*===========================================================================*/
