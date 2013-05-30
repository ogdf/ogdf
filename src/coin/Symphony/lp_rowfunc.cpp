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

#include <memory.h>
#include <math.h>
#include <string.h>

#include "sym_lp.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_types.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP functions related to row operations.
\*===========================================================================*/

int check_row_effectiveness(lp_prob *p)
{
   int ineff_cnt_to_delete = p->par.ineff_cnt_to_delete;
   int orig_eff = p->par.base_constraints_always_effective;
   LPdata *lp_data = p->lp_data;
   //double *dualsol = lp_data->dualsol;
   double lpetol = lp_data->lpetol;
   double lpetol10 = 10*lpetol;
   row_data *row, *rows = lp_data->rows;
   int m = lp_data->m;

   int bcutnum = p->base.cutnum;

   double slack, *slacks = lp_data->slacks;
   int *free_rows;
   int ineffective, deletable, violated, i, j, k;

   char *stat;
   int  *now_ineff, *outrhsind, *inrhsind, *slackstat;

   stat = lp_data->tmp.c; /* m */

   /* slacks is already filled up. We got the slacks before calling
    * fix_variables.
    * Now based on their slack values, mark each row whether it's
    * violated, loose or tight */

   //int base_m = orig_eff ? bcutnum : 0;

   for (i = m - 1; i >= 0; i--){
      slack = slacks[i];
      switch (rows[i].cut->sense){
       case 'E':
	 if (slack < -lpetol10 || slack > lpetol10) stat[i] = VIOLATED_ROW;
	 else                                   stat[i] = TIGHT_ROW;
	 break;
       case 'L':
	 if (slack > lpetol10)       stat[i] = SLACK_ROW;
	 else if (slack > -lpetol10) stat[i] = TIGHT_ROW;
	 else                        stat[i] = VIOLATED_ROW;
	 break;
       case 'G':
	 if (slack < -lpetol10)     stat[i] = SLACK_ROW;
	 else if (slack < lpetol10) stat[i] = TIGHT_ROW;
	 else                       stat[i] = VIOLATED_ROW;
	 break;
       case 'R':
	 if (rows[i].cut->range < 0){
	    if (slack > rows[i].cut->range + lpetol10 || slack < -lpetol10)
	       stat[i] = SLACK_ROW;
	    else if (slack > rows[i].cut->range - lpetol10 || slack < lpetol10)
	       stat[i] = TIGHT_ROW;
	    else
	       stat[i] = VIOLATED_ROW;
	 }else{
	    if (slack < rows[i].cut->range - lpetol10 || slack > lpetol10)
	       stat[i] = SLACK_ROW;
	    else if (slack < rows[i].cut->range + lpetol10 || slack > - lpetol10)
	       stat[i] = TIGHT_ROW;
	    else
	       stat[i] = VIOLATED_ROW;
	 }
	 break;
      }
   }

   /* Now set the branch values appropriately */
   if (p->par.branch_on_cuts){
      for (i=m-1; i>=0; i--){
	 if (stat[i] == SLACK_ROW){
	    if ((rows[i].cut->branch & ALLOWED_TO_BRANCH_ON))
	       rows[i].cut->branch ^= SWITCH_CANDIDATE_ALLOWED;
	 }else{
	    if ((rows[i].cut->branch & CANDIDATE_FOR_BRANCH))
	       rows[i].cut->branch ^= SWITCH_CANDIDATE_ALLOWED;
	 }
      }
   }

   /*========================================================================*\
    * A few words of wisdom:
    *
    * If a row wasn't free then everything is nice.
    * If it was free then there are complications because its slack variable
    * will be in basis and the corresponding dual variable is 0.
    *
    * Keep in mind the objective: if violated then make it constraining,
    * if tight and free keep it free, if slack, make it free.
    *
    * Also, base constraints and branching cuts may be ineffective but they
    * are not deletable.
    *   So be careful!
   \*========================================================================*/

   violated = ineffective = 0;

   /* we'll first use slackstat then outrhsind. no conflict */
   slackstat = outrhsind = lp_data->tmp.i1;
   inrhsind = outrhsind + m;
   now_ineff = inrhsind + m;
   if (p->par.ineffective_constraints != NO_CONSTRAINT_IS_INEFFECTIVE){
      /* Deal with the violated rows */
      for (i = orig_eff ? bcutnum : 0; i < m; i++){
	 if (stat[i] == VIOLATED_ROW){ /* must have been free */
	    rows[i].free = FALSE;
	    rows[i].eff_cnt = 0;
	    rows[i].ineff_cnt = 0;
	    inrhsind[violated++] = i;
	 }
      }
      /* Collect the rows that are deemed ineffective now */
      switch (p->par.ineffective_constraints){
       case NONZERO_SLACKS_ARE_INEFFECTIVE:
	 for (i = orig_eff ? bcutnum : 0; i < m; i++){
	    if (stat[i] == SLACK_ROW ||
		(stat[i] == TIGHT_ROW && rows[i].free == TRUE)){
	       now_ineff[ineffective++] = i;
	    }else{
	       rows[i].eff_cnt++;
	    }
	 }
	 break;
       case BASIC_SLACKS_ARE_INEFFECTIVE:
	 /* for violated free rows the slack is in basis! */
	 get_basis(lp_data, NULL, slackstat);
	 for (i = orig_eff ? bcutnum : 0; i < m; i++){
	    if (slackstat[i] == SLACK_BASIC && stat[i] != VIOLATED_ROW){
	       now_ineff[ineffective++] = i;
	    }else{
	       rows[i].eff_cnt++;
	    }
	 }
	 break;
       case ZERO_DUAL_VALUES_ARE_INEFFECTIVE:
	 for (i = orig_eff ? bcutnum : 0; i < m; i++){
 	    if (fabs(lp_data->dualsol[i]) < lpetol && stat[i] != VIOLATED_ROW){
 	       now_ineff[ineffective++] = i;
 	    }else{
 	       rows[i].eff_cnt++;
	    }
	 }
	 break;
      }
      /* Now violated rows have eff_cnt = 1 (not that it matters...) */
   }

   deletable = k = 0;
   for (j = ineffective - 1; j >= 0; j--){

      row = rows + (i = now_ineff[j]);

      if(p->bc_level > 100 && !(row->deletable))row->deletable = TRUE;

      if (!row->free && row->deletable){
	 row->free = TRUE;
	 row->ineff_cnt = stat[i] == TIGHT_ROW ? 0 : ((MAXINT) >> 1);
	 outrhsind[k++] = i;
      }
      row->ineff_cnt++;
      if (i >= bcutnum && ! (row->cut->branch & CUT_BRANCHED_ON) &&
	  row->deletable && row->ineff_cnt >= ineff_cnt_to_delete )
        deletable++;
   }

   /* stat is not used any more so its location can be used in
      constrain_row_set and free_row_set, but for integer tmp array
      they have to use the area behind in/outrhsind. This area was used
      for now_ineff, but we don't need that anymore either. */

   if (violated > 0)
      constrain_row_set(lp_data, violated, inrhsind);

   if (k > 0)
      free_row_set(lp_data, k, outrhsind);

   PRINT(p->par.verbosity, 3,
	 ("Row effectiveness: rownum: %i ineff: %i deletable: %i\n",
	  m, ineffective, deletable));
   if (p->par.verbosity > 6 && ineffective){
      printf("   Ineffective row(s):");
      for (i=0; i<m; i++){
	 if (rows[i].free)
	    printf(" %i", i);
      }
      printf("\n");
   }

   /*------------------------------------------------------------------------*\
    * Finally, remove the deletable rows if there are enough to remove
   \*------------------------------------------------------------------------*/

   if (deletable > p->par.mat_row_compress_ratio * m &&
       deletable > p->par.mat_row_compress_num){
      PRINT(p->par.verbosity, 3, ("   Removing deletable rows ...\n"));
      if (p->par.branch_on_cuts)
	 p->slack_cuts = (cut_data **) realloc(p->slack_cuts,
			 (p->slack_cut_num + deletable) * sizeof(cut_data *));

      free_rows = lp_data->tmp.i1;
      if (bcutnum > 0)
	 memset(free_rows, FALSE, bcutnum * ISIZE);
      /* remember, by now every ineffective row is free and ineff_cnt is
	 positive only for free rows */
      for (k = i = bcutnum; i < m; i++){
	 row = rows + i;
	 if (row->free && ! (row->cut->branch & CUT_BRANCHED_ON) &&
	     row->ineff_cnt >= ineff_cnt_to_delete){
	    free_rows[i] = TRUE;
	    if (row->cut->branch & CANDIDATE_FOR_BRANCH){
#ifdef DO_TESTS
	       if (!p->par.branch_on_cuts)
		  printf("No branch_on_cuts but a CANDIDATE_FOR_BRANCH!\n\n");
#endif
	       p->slack_cuts[p->slack_cut_num++] = row->cut;
	       row->cut = NULL;
	    }else{
#ifdef COMPILE_IN_LP /*we don't want to free rows that have a name if we are
		       using shared memory because they are still being used*/
	       if (row->cut->name < 0)
#endif
		  free_cut(&(row->cut));
	    }
	 }else{
	    free_rows[i] = FALSE;
	    rows[k++] = rows[i];
	 }
      }
      delete_rows(lp_data, deletable, free_rows);
      p->lp_stat.cuts_deleted_from_lps += deletable;
      if (p->bc_level > 0) {
         p->lp_stat.num_cuts_slacked_out_in_path += deletable;
      }
   }
   PRINT(p->par.verbosity, 3, ("\n"));

   return(violated);
}

/*===========================================================================*/

void add_row_set(lp_prob *p, waiting_row **wrows, int length)
{
   int i;
   row_data *row;

   add_waiting_rows(p, wrows, length);

   row = p->lp_data->rows + (p->lp_data->m - length);

   for (i=0; i<length; i++, row++){
      row->free = FALSE;
      row->cut = wrows[i]->cut;
      row->eff_cnt = 1;
      row->deletable = wrows[i]->cut->deletable;
      wrows[i]->cut = NULL;
   }

   free_waiting_rows(wrows, length);
}

/*===========================================================================*/

void add_new_rows_to_waiting_rows(lp_prob *p, waiting_row **new_rows,
				  int new_row_num)
{
   new_row_num = compute_violations(p, new_row_num, new_rows);

   if (new_row_num > 0){
      /* check to be sure there is enough room in the row set data
       * structure for the new rows -- otherwise reallocate memory */
      REALLOC(p->waiting_rows, waiting_row *, p->waiting_rows_size,
	      p->waiting_row_num + new_row_num, BB_BUNCH);
      memcpy((p->waiting_rows + p->waiting_row_num), new_rows,
	     new_row_num * sizeof(waiting_row *));
      p->waiting_row_num += new_row_num;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * This function orders the waiting rows so that the cuts sent by the same
 * process are grouped together, but otherwise preserving the order of
 * arrival.
 * The ordering is done in ascending order wrt the source_pid field of each
 * waiting_row. Newly arriving cuts have the correct value here, cuts already
 * in the local pool get MAXINT, so they are considered last.

 * NOTE: This ensures that results are reproducible, even with with multiple
 * cut pools/generators, as long as we never time out.
\*===========================================================================*/

void order_waiting_rows_based_on_sender(lp_prob *p)
{
   waiting_row **wrows = p->waiting_rows;
   waiting_row *wtmp;
   const int wrownum = p->waiting_row_num;
   int i, j;
   /* Do a simple bubble sort */
   for (i = 1; i < wrownum; ++i) {
      wtmp = wrows[i];
      for (j = i - 1; j >= 0; --j) {
	 if (wtmp->source_pid >= wrows[j]->source_pid) {
	    break;
	 } else {
	    wrows[j+1] = wrows[j];
	 }
      }
      wrows[j+1] = wtmp;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Order the cuts in waiting_rows based on their violation then pick
 * the best k, and add those to the problem
\============================================================================*/

int add_best_waiting_rows(lp_prob *p)
{
   int i, added_rows;
   row_data *rows;
   int max_cut_num_per_iter = (p->bc_level<1)?p->par.max_cut_num_per_iter_root:
                                            p->par.max_cut_num_per_iter;

   added_rows = MIN(max_cut_num_per_iter, p->waiting_row_num);
   if (added_rows < p->waiting_row_num)
      qsort((char *)p->waiting_rows, p->waiting_row_num,
	    sizeof(waiting_row *), waiting_row_comp);
   if (added_rows){
      print_stat_on_cuts_added_u(p, added_rows);
      add_row_set(p, p->waiting_rows, added_rows);
      rows = p->lp_data->rows + (p->lp_data->m - added_rows);
      for (i=0; i<added_rows; i++){
	 rows[i].eff_cnt = 1;
      }
      if (added_rows < p->waiting_row_num)
	 memmove(p->waiting_rows, p->waiting_rows + added_rows,
	       (p->waiting_row_num - added_rows) * sizeof(waiting_row *));
      p->waiting_row_num -= added_rows;
   }
   return(added_rows);
}

/*===========================================================================*/

void add_waiting_rows(lp_prob *p, waiting_row **wrows, int add_row_num)
{
   LPdata *lp_data = p->lp_data;
   char *sense;
   double *rhs, *rmatval;
   int *rmatbeg, *rmatind;
   int i, nzcnt;
   waiting_row *wrow;

   for (nzcnt=0, i=add_row_num-1; i>=0; i--)
      nzcnt += wrows[i]->nzcnt;

   size_lp_arrays(lp_data, TRUE, FALSE, add_row_num, 0, nzcnt);

   sense = lp_data->tmp.c; /* m */
   rhs = lp_data->tmp.d; /* m */
   REMALLOC(lp_data->tmp.dv, double, lp_data->tmp.dv_size, nzcnt,
         5*(int)BB_BUNCH);
   rmatval = lp_data->tmp.dv; /* nzcnt */
   rmatbeg = lp_data->tmp.i1;
   REMALLOC(lp_data->tmp.iv, int, lp_data->tmp.iv_size, nzcnt, 5*(int)BB_BUNCH);
   rmatind = lp_data->tmp.iv;

   *rmatbeg = 0;
   for (i = 0; i < add_row_num; i++){
      wrow = wrows[i];
      rhs[i] = wrow->cut->rhs;
      sense[i] = wrow->cut->sense;
      memcpy(rmatind + rmatbeg[i], wrow->matind, wrow->nzcnt * ISIZE);
      memcpy(rmatval + rmatbeg[i], wrow->matval, wrow->nzcnt * DSIZE);
      rmatbeg[i+1] = rmatbeg[i] + wrow->nzcnt;
   }

   add_rows(lp_data, add_row_num, nzcnt, rhs, sense, rmatbeg, rmatind,rmatval);

   for (i = add_row_num - 1; i >= 0; i--){
      if (sense[i] == 'R')
	 change_range(lp_data, lp_data->m+i, wrows[i]->cut->range);
   }
}

/*===========================================================================*\
 * This function is compares two waiting rows. Needed for ordering the cuts by
 * degree of violation.
\*===========================================================================*/

int waiting_row_comp(const void *wr0, const void *wr1)
{
   double v0 = (*((waiting_row **)wr0))->violation;
   double v1 = (*((waiting_row **)wr1))->violation;
   return(v0 < v1 ? 1 : (v0 > v1 ?  -1 : 0));
}

/*===========================================================================*/

int compute_violations(lp_prob *p, int new_row_num, waiting_row **new_rows)
{
   waiting_row *wrow;
   int *matind, i, j;
   double *matval, lpetol = p->lp_data->lpetol, lhs, *x = p->lp_data->x;
   cut_data *cut;

   for (i = 0; i < new_row_num; ){
      wrow = new_rows[i];
      matind = wrow->matind;
      matval = wrow->matval;
      for (lhs=0, j = wrow->nzcnt-1; j>=0; j--)
	 lhs += matval[j] * x[matind[j]];
      cut = wrow->cut;
      switch (cut->sense){
       case 'L': wrow->violation = lhs - cut->rhs;       break;
       case 'G': wrow->violation = cut->rhs - lhs;       break;
       case 'E': wrow->violation = fabs(lhs - cut->rhs); break;
       case 'R':
	 wrow->violation =
	    lhs < cut->rhs ? cut->rhs - lhs : lhs - cut->rhs - cut->range;
	 break;
      }
      if  (wrow->violation < lpetol){
	 free_waiting_row(new_rows+i);
	 new_rows[i] = new_rows[--new_row_num];
      }else{
	 i++;
      }
   }
   return(new_row_num);
}

/*===========================================================================*/

void compress_slack_cuts(lp_prob *p)
{
   int i, snum = p->slack_cut_num;
   cut_data **slack_cuts = p->slack_cuts;

   for (i=0; i<snum; ){
      if (slack_cuts[i] == NULL){
	 slack_cuts[i] = slack_cuts[--snum];
      }else{
	 i++;
      }
   }
   p->slack_cut_num = snum;
}

