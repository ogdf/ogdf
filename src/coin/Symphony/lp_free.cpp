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


#include "sym_lp.h"
#include "sym_types.h"
#include "sym_macros.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains LP functions related to freeing data structures
\*===========================================================================*/

void free_cut(cut_data **cut)
{
   if (*cut){
      if ((*cut)->coef){
	 FREE((*cut)->coef);
      }
      FREE(*cut);
   }
}

/*===========================================================================*/

void free_cuts(cut_data **cuts, int cut_num)
{
   int i;
   if (cuts)
      for (i=cut_num-1; i>=0; i--)
	 if (cuts[i])
#ifdef COMPILE_IN_LP
	    if (cuts[i]->name < 0 || cuts[i]->branch & CUT_BRANCHED_ON)
#endif
	       free_cut(cuts+i);
}

/*===========================================================================*/

void free_col_set(our_col_set **colset)
{
   if (*colset){
      our_col_set *cols = *colset;
      FREE(cols->rel_lb_ind);
      FREE(cols->rel_ub_ind);
      FREE(cols->userind);
      FREE(cols->objx);
      FREE(cols->lb);
      FREE(cols->ub);
      FREE(cols->matbeg);
      FREE(cols->matind);
      FREE(cols->matval);
      FREE(*colset);
   }
}

/*===========================================================================*/

void free_candidate(branch_obj **cand)
{
   int i;

   if (*cand){
      branch_obj *can = *cand;
#ifdef COMPILE_FRAC_BRANCHING
      for (i = can->child_num-1; i >= 0; i--){
	 if (can->frac_num[i]){
	    FREE(can->frac_ind[i]);
	    FREE(can->frac_val[i]);
	 }
      }
#endif
      free_waiting_row(&(can->row));
#ifndef MAX_CHILDREN_NUM
      FREE(can->sense); FREE(can->rhs); FREE(can->range); FREE(can->branch);

      if (can->solutions){
	 for (i = can->child_num-1; i >= 0; i--){
#else
      if (can->solutions){
         for (i = MAX_CHILDREN_NUM - 1; i >= 0; i--){
#endif
	    FREE(can->sol_inds[i]);
	    FREE(can->solutions[i]);
	 }
      }

#ifdef SENSITIVITY_ANALYSIS
#ifndef MAX_CHILDREN_NUM
      if (can->duals){
	 for (i = can->child_num-1; i >= 0; i--){
#else
      if (can->duals){
	 for (i = MAX_CHILDREN_NUM - 1; i >= 0; i--){
#endif
	    FREE(can->duals[i]);
	 }
      }
#endif

      FREE(can->sol_sizes);
      FREE(can->sol_inds);
      FREE(can->solutions);
#ifdef SENSITIVITY_ANALYSIS
      FREE(can->duals);
#endif

      FREE(*cand);
   }
}

/*===========================================================================*/

void free_candidate_completely(branch_obj **cand)
{
   if (*cand){
#ifndef MAX_CHILDREN_NUM
      branch_obj *can = *cand;
#endif
#ifndef MAX_CHILDREN_NUM
      FREE(can->objval);
      FREE(can->termcode);
      FREE(can->feasible);
      FREE(can->iterd);
#  ifdef COMPILE_FRAC_BRANCHING
      FREE(can->frac_num); FREE(can->frac_ind); FREE(can->frac_val);
#  endif
#endif
      free_candidate(cand);
   }
}

/*===========================================================================*/

void free_waiting_row(waiting_row **wrow)
{
   waiting_row *wr = *wrow;
   if (wr){
      FREE(wr->matval);
      FREE(wr->matind);
      free_cut(&wr->cut);
      free(wr);
      *wrow = NULL;
   }
}

/*===========================================================================*/

void free_waiting_rows(waiting_row **rows, int row_num)
{
   int i;
   if (rows)
      for (i=row_num-1; i>=0; i--)
	 free_waiting_row(rows+i);
}

/*===========================================================================*/

void free_waiting_row_array(waiting_row ***rows, int row_num)
{
   free_waiting_rows(*rows, row_num);
   FREE(*rows);
}

/*===========================================================================*/

void free_node_desc(node_desc **desc)
{
   if (*desc){
      node_desc *n = *desc;
      FREE(n->cutind.list);
      FREE(n->uind.list);
      if (n->nf_status == NF_CHECK_AFTER_LAST ||
	  n->nf_status == NF_CHECK_UNTIL_LAST)
	 FREE(n->not_fixed.list);
      if (n->basis.basis_exists){
	 FREE(n->basis.basevars.list);
	 FREE(n->basis.basevars.stat);
	 FREE(n->basis.extravars.list);
	 FREE(n->basis.extravars.stat);
	 FREE(n->basis.baserows.list);
	 FREE(n->basis.baserows.stat);
	 FREE(n->basis.extrarows.list);
	 FREE(n->basis.extrarows.stat);
      }
      if (n->desc_size > 0)
	 FREE(n->desc);
      if (n->bnd_change) {
         FREE(n->bnd_change->index);
         FREE(n->bnd_change->lbub);
         FREE(n->bnd_change->value);
         FREE(n->bnd_change);
      }
      FREE(*desc);
   }
}

/*===========================================================================*/

void free_node_dependent(lp_prob *p)
{
   LPdata *lp_data = p->lp_data;
   int i;

   free_node_desc(&p->desc);
   for (i = p->base.cutnum; i < lp_data->m; i++){
#ifdef COMPILE_IN_LP
      if (lp_data->rows[i].cut->name < 0 ||
	  lp_data->rows[i].cut->branch & CUT_BRANCHED_ON)
#endif
	 free_cut(&lp_data->rows[i].cut);
#ifdef COMPILE_IN_LP
      else
	 lp_data->rows[i].cut = NULL;
#endif
   }
   if (p->par.branch_on_cuts && p->slack_cut_num > 0){
      free_cuts(p->slack_cuts, p->slack_cut_num);
      p->slack_cut_num = 0;
   }
   // necessary to purge waiting rows, otherwise these may get added to the
   // node that is solved next time.
   if (p->waiting_row_num>0) {
      free_waiting_rows(p->waiting_rows, p->waiting_row_num);
      p->waiting_row_num = 0;
      FREE(p->waiting_rows);
   }

   unload_lp_prob(lp_data);
}

/*===========================================================================*/

void free_lp(lp_prob *p)
{
   int i;

   free_prob_dependent_u(p);
   free_waiting_row_array(&p->waiting_rows, p->waiting_row_num);
   for (i = p->lp_data->maxn - 1; i >= 0; i--)
      FREE(p->lp_data->vars[i]);
   FREE(p->lp_data->vars);
#ifdef COMPILE_IN_LP
   for (i = p->base.cutnum - 1; i >= 0; i--)
      free_cut(&(p->lp_data->rows[i].cut));
   free_node_desc(&p->desc);
#else
   for (i = p->lp_data->m - 1; i >= 0; i--)
      free_cut(&(p->lp_data->rows[i].cut));
   FREE(p->bdesc);
#endif
   FREE(p->lp_data->rows);
   close_lp_solver(p->lp_data);
   free_lp_arrays(p->lp_data);
   if (p->par.lp_data_mip_is_copied == TRUE) {
      free_mip_desc(p->lp_data->mip);
   }
   FREE(p->lp_data->mip);
   FREE(p->lp_data);
#if !(defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP))
   free_mip_desc(p->mip);
   FREE(p->mip);
   FREE(p->base.userind);
#endif
   FREE(p->best_sol.xind);
   FREE(p->best_sol.xval);
   if (p->par.branch_on_cuts){
      FREE(p->slack_cuts);
   }
   FREE(p->obj_history);
   FREE(p);
}

