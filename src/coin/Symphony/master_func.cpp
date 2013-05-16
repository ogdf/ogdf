/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* This file was developed by Menal Guzelsoy for the SYMPHONY OSI interface. */
/*                                                                           */
/* (c) Copyright 2000-2011 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#define COMPILING_FOR_MASTER

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __PVM__
#include <pvmtev.h>
#endif

#include "sym_master.h"
#include "sym_master_u.h"
#include "sym_macros.h"
#include "sym_pack_cut.h"
#include "sym_pack_array.h"
#include "sym_lp_solver.h"
//#include "sym_lp.h"
#include "sym_tm.h"

/*===========================================================================*/
/*===========================================================================*/

/* not used now! to be used later!*/
int resolve_node(sym_environment *env, bc_node *node)
{
   node_desc * desc = &node->desc;
   LPdata *lp_data = (LPdata*)calloc(1, sizeof(LPdata));
   lp_data->mip = create_copy_mip_desc(env->mip); //FIXME!!!
   lp_sol *best_sol = &(env->warm_start->best_sol);
   branch_desc *bpath;
   branch_obj *bobj;
   bc_node **path, *n;
   int level = node->bc_level;

   int *matbeg = 0, *matind= 0, size, nzcnt = 0, return_value, iterd = 0;
   double *matval = 0, colsol;
   int i, j = 0, cnt = 0, *xind = 0;
   double *xval = 0, *rhs = 0, lpetol = 9.9999999999999995e-07;
   char *sense = 0;
   cut_data *cut;

   /*------------------------------------------------------------------------*\
    * Now go through the branching stuff
   \*----------------------------------------------------------------------- */

   path = (bc_node **) malloc((2*(level+1)+BB_BUNCH)*sizeof(bc_node *));

   for (i = level, n = node; i >= 0; n = n->parent, i--)
      path[i] = n;

   /*------------------------------------------------------------------------*\
    * Read in the basis.
    * This is cplex style. sorry about it... Still, it
    * might be ok if {VAR,SLACK}_{B,LB,UB} are properly defined
   \*----------------------------------------------------------------------- */

   node_desc *new_desc;

   int varexp_ind = 0, cutexp_ind = 0, nfexp_ind = 0;
   int bv_ind = 0, br_ind = 0, ev_ind = 0, er_ind = 0;
   array_desc extravar = { EXPLICIT_LIST, 0, 0, NULL };
   array_desc extrarow = { EXPLICIT_LIST, 0, 0, NULL };
   array_desc not_fixed = { EXPLICIT_LIST, 0, 0, NULL };
   basis_desc basis;

   int *list, *stat;

   char deal_with_nf = (desc->nf_status == NF_CHECK_AFTER_LAST ||
			desc->nf_status == NF_CHECK_UNTIL_LAST);

   memset((char *)(&basis), 0, sizeof(basis_desc));

   /*------------------------------------------------------------------------*\
    * First go up in the search tree to the root and record for every field
    * the node where the first explicit description occurs.
   \*------------------------------------------------------------------------*/

   if (desc->uind.type == NO_DATA_STORED){
      varexp_ind = -1;
   }else{
      for (i = level, n = node; !varexp_ind && i > 0; n = n->parent, i--)
	 if (n->desc.uind.type == EXPLICIT_LIST)
	    varexp_ind = i;
   }
   if (desc->cutind.type == NO_DATA_STORED){
      cutexp_ind = -1;
   }else{
      for (i = level, n = node; !cutexp_ind && i > 0; n = n->parent, i--)
	 if (n->desc.cutind.type == EXPLICIT_LIST)
	    cutexp_ind = i;
   }
   if (deal_with_nf){
      for (i = level, n = node; !nfexp_ind && i > 0; n = n->parent, i--)
	 if (n->desc.not_fixed.type == EXPLICIT_LIST)
	    nfexp_ind = i;
   }
   if ((basis.basis_exists = desc->basis.basis_exists) == TRUE){
      for (i = level, n = node; !bv_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.basevars.type == EXPLICIT_LIST)
	    bv_ind = i;
      for (i = level, n = node; !br_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.baserows.type == EXPLICIT_LIST)
	    br_ind = i;
      for (i = level, n = node; !ev_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.extravars.type == EXPLICIT_LIST)
	    ev_ind = i;
      for (i = level, n = node; !er_ind && i > 0; n = n->parent, i--)
	 if (n->desc.basis.extrarows.type == EXPLICIT_LIST)
	    er_ind = i;
   }else{
      ev_ind = er_ind = level;
   }

   /* An upper estimate on the total length of arrays */
   if (varexp_ind >= 0){
      extravar.size = (n = path[varexp_ind]) -> desc.uind.size;
      for (i = varexp_ind + 1; i <= level; i++)
	 extravar.size += path[i]->desc.uind.added;
   }
   if (cutexp_ind >= 0){
      extrarow.size = (n = path[cutexp_ind]) -> desc.cutind.size;
      for (i = cutexp_ind + 1; i <= level; i++)
	 extrarow.size += path[i]->desc.cutind.added;
   }
   if (deal_with_nf && nfexp_ind >= 0){
      not_fixed.size = (n = path[nfexp_ind]) -> desc.not_fixed.size;
      for (i = nfexp_ind + 1; i <= level; i++)
	 not_fixed.size += path[i]->desc.not_fixed.added;
   }else{
      not_fixed.size = 0;
   }

   /* If the LP function is compiled into the tree manager as a single
      executable, then we allocate new memory for these arrays since
      these arrays will be used directly instead of being passed
      through PVM. */
   if (extravar.size){
      extravar.list = (int *) malloc(extravar.size*ISIZE);
      if (basis.basis_exists)
	 basis.extravars.stat = (int *) malloc(extravar.size*ISIZE);
   }
   if (extrarow.size){
      extrarow.list = (int *) malloc(extrarow.size*ISIZE);
      if (basis.basis_exists)
	 basis.extrarows.stat = (int *) malloc(extrarow.size*ISIZE);
   }
   if (not_fixed.size)
      not_fixed.list = (int *) malloc(not_fixed.size*ISIZE);
   if (env->base->varnum && basis.basis_exists)
      basis.basevars.stat = (int *) malloc(env->base->varnum*ISIZE);
   if (env->base->cutnum && basis.basis_exists)
      basis.baserows.stat = (int *) malloc(env->base->cutnum*ISIZE);

   /* The extra variables (uind) and the corresponding basis part */
   if (varexp_ind >= 0){
      extravar.size = (n = path[varexp_ind]) -> desc.uind.size;
      if (extravar.size > 0)
	 memcpy(extravar.list, n->desc.uind.list, ISIZE * extravar.size);
      for (i = varexp_ind + 1; i <= ev_ind; i++)
	 modify_list(&extravar, &path[i]->desc.uind);
      if (basis.basis_exists){
	 /* at this point i == ev_ind */

	 if (path[ev_ind]->desc.basis.extravars.size > 0)
	    memcpy(basis.extravars.stat,
		   path[ev_ind]->desc.basis.extravars.stat,
		   path[ev_ind]->desc.basis.extravars.size * ISIZE);
	 for (i = ev_ind + 1; i <= level; i++){
	    modify_list_and_stat(&extravar, basis.extravars.stat,
				 &path[i]->desc.uind,
				 &path[i]->desc.basis.extravars);
	 }
	 /* Although we send an explicit list, the type is sent over to show
	    to the LP process how extravars are stored in TM */
	 basis.extravars.type = node->desc.basis.extravars.type;
	 basis.extravars.size = extravar.size;
	 basis.extravars.list = NULL;
	 /* Now extravar.list/extravar.size and basis.extravars are OK */
	 /* Fix basis.basevars */
	 basis.basevars.type = EXPLICIT_LIST;
	 basis.basevars.size = path[bv_ind]->desc.basis.basevars.size;
	 basis.basevars.list = NULL;
	 if (basis.basevars.size > 0){
	    memcpy(basis.basevars.stat,
		   path[bv_ind]->desc.basis.basevars.stat,
		   basis.basevars.size * ISIZE);
	    for (i = bv_ind + 1; i <= level; i++){
	       list = path[i]->desc.basis.basevars.list;
	       stat = path[i]->desc.basis.basevars.stat;
	       for (j = path[i]->desc.basis.basevars.size - 1; j >= 0; j--)
		  basis.basevars.stat[list[j]] = stat[j];
	    }
	 }
      }
   }


   /* Now take care of cutind and the corresponding basis part */
   if (cutexp_ind >= 0){
      extrarow.size = (n = path[cutexp_ind]) -> desc.cutind.size;
      if (extrarow.size > 0)
	 memcpy(extrarow.list, n->desc.cutind.list, ISIZE * extrarow.size);
      for (i = cutexp_ind + 1; i <= er_ind; i++)
	 modify_list(&extrarow, &path[i]->desc.cutind);
      if (basis.basis_exists){
	 /* at this point i == er_ind */
	 if (path[er_ind]->desc.basis.extrarows.size > 0)
	    memcpy(basis.extrarows.stat,
		   path[er_ind]->desc.basis.extrarows.stat,
		   path[er_ind]->desc.basis.extrarows.size * ISIZE);
	 for (i = er_ind + 1; i <= level; i++){
	    modify_list_and_stat(&extrarow, basis.extrarows.stat,
				 &path[i]->desc.cutind,
				 &path[i]->desc.basis.extrarows);
	 }
	 /* Same trick as above */
	 basis.extrarows.type = node->desc.basis.extrarows.type;
	 basis.extrarows.size = extrarow.size;
	 basis.extrarows.list = NULL;
	 /* Now extrarow.list/extrarow.size and basis.extrarows are OK */
	 /* Fix basis.baserows */
	 basis.baserows.type = EXPLICIT_LIST;
	 basis.baserows.size = path[br_ind]->desc.basis.baserows.size;
	 basis.baserows.list = NULL;
	 if (basis.baserows.size > 0){
	    memcpy(basis.baserows.stat,
		   path[br_ind]->desc.basis.baserows.stat,
		   basis.baserows.size * ISIZE);
	    for (i = br_ind + 1; i <= level; i++){
	       list = path[i]->desc.basis.baserows.list;
	       stat = path[i]->desc.basis.baserows.stat;
	       for (j = path[i]->desc.basis.baserows.size - 1; j >= 0; j--)
		  basis.baserows.stat[list[j]] = stat[j];
	    }
	 }
      }
   }

   /* Finally the not fixed ones */
   if (deal_with_nf){
      not_fixed.size = (n = path[nfexp_ind]) -> desc.not_fixed.size;
      if (not_fixed.size > 0)
	 memcpy(not_fixed.list, n->desc.not_fixed.list, ISIZE*not_fixed.size);
      for (i = nfexp_ind + 1; i <= level; i++)
	 modify_list(&not_fixed, &path[i]->desc.not_fixed);
   }

   new_desc = (node_desc *) calloc(1,sizeof(node_desc));

   new_desc->nf_status = desc->nf_status;
   new_desc->basis = basis;
   if (deal_with_nf)
      new_desc->not_fixed = not_fixed;
   new_desc->uind = extravar;
   new_desc->cutind = extrarow;

#if 0
   /* The cuts themselves */
   if (extrarow.size > 0){
      new_desc->cuts = (cut_data **)
	 malloc(extrarow.size*sizeof(cut_data *));
      for (i = 0; i < extrarow.size; i++){
	 new_desc->cuts[i] = tm->cuts[extrarow.list[i]];
      }
   }
#endif

   /* User defined description */
   new_desc->desc_size = desc->desc_size;
   if (new_desc->desc_size > 0)
      memcpy((char *)new_desc->desc, (char *)desc->desc, new_desc->desc_size);


   /*------------------------------------------------------------------------*\
    * Load the lp problem (load_lp is an lp solver dependent routine).
   \*----------------------------------------------------------------------- */

   //   lp_data->mip->m += new_desc->cutind.size;

   lp_data->m = lp_data->mip->m;
   lp_data->n = lp_data->mip->n;

   open_lp_solver(lp_data);
   load_lp_prob(lp_data, 0, 0);


   /*------------------------------------------------------------------------*\
    * Now go through the branching stuff
   \*----------------------------------------------------------------------- */

   bpath = (branch_desc *) malloc
      ((2*(level+1)+BB_BUNCH)*sizeof(branch_desc));

   for (i = 0; i < level; i++, bpath++){
      for (j = path[i]->bobj.child_num - 1; j >= 0; j--)
	 if (path[i]->children[j] == path[i+1])
	    break;
      bobj = &path[i]->bobj;
      bpath->type = bobj->type;
      bpath->name = bobj->name;
      bpath->sense = bobj->sense[j];
      bpath->rhs = bobj->rhs[j];
      bpath->range = bobj->range[j];
      bpath->branch = bobj->branch[j];
   }

   bpath = bpath - level;
   if (level){
      for (i = 0; i < level; i++, bpath++){
	 if (bpath->type == BRANCHING_VARIABLE){
	    j = bpath->name;
	    switch (bpath->sense){
	     case 'E':
		change_lbub(lp_data, j, bpath->rhs, bpath->rhs);
		break;
	     case 'L':
		change_ub(lp_data, j, bpath->rhs);
	       break;
	     case 'G':
		change_lb(lp_data, j, bpath->rhs);
	       break;
	     case 'R':
		change_lbub(lp_data, j, bpath->rhs, bpath->rhs + bpath->range);
		break;
	    }
	 }else{ /* BRANCHING_CUT */
	    j = bpath->name;
	    change_row(lp_data, j, bpath->sense, bpath->rhs,bpath->range);
	 }
      }
   }
   bpath=bpath-level;
   /*------------------------------------------------------------------------*\
    * Add cuts here
    * FIXME! ASSUMING ALL THE CUTS ARE EXPLICIT ROW!
   \*----------------------------------------------------------------------- */
   desc = new_desc;

   if (desc->cutind.size > 0){
      size = desc->cutind.size;
      sense  = (char*) malloc(size*CSIZE);
      rhs = (double*) malloc(size*DSIZE);
      matbeg = (int *) calloc(size + 1, ISIZE);
      matbeg[0] = 0;

      for (i = 0; i<env->warm_start->cut_num; i++){
	 for(j=0; j<desc->cutind.size; j++){
	    if (i == desc->cutind.list[j]){
	       cut = env->warm_start->cuts[i];
	       nzcnt = ((int *) (cut->coef))[0];
	       sense[j] = cut->sense;
	       rhs[j] = cut->rhs;
	       matbeg[j+1] = matbeg[j] + nzcnt;
	    }
	 }
      }

      matind = (int *) malloc(nzcnt*ISIZE);
      matval = (double *) malloc(nzcnt*DSIZE);

      for (i = 0; i<env->warm_start->cut_num; i++){
	 for(j=0; j<desc->cutind.size; j++){
	    if (i == desc->cutind.list[j]){
	       cut = env->warm_start->cuts[i];
	       nzcnt = matbeg[j+1] - matbeg[j];
	       memcpy(matind + matbeg[j], (int *) (cut->coef + ISIZE),
		      ISIZE * nzcnt);
	       memcpy(matval + matbeg[j],
		      (double *) (cut->coef + (1 + nzcnt) * ISIZE),
		      DSIZE * nzcnt);
	    }
	 }
      }
      nzcnt = matbeg[j];
      add_rows(lp_data, size, nzcnt, rhs, sense, matbeg, matind, matval);
   }

   /*----------------------------------------------------------------------- */
   /* Load The Basis */
   /*----------------------------------------------------------------------- */

   if (desc->basis.basis_exists == TRUE){
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

   return_value = dual_simplex(lp_data, &iterd);

   if(return_value == LP_D_UNBOUNDED || return_value == LP_ABANDONED ||
      return_value == LP_D_INFEASIBLE){
      node->feasibility_status = INFEASIBLE_PRUNED;
      node->node_status = NODE_STATUS__PRUNED;
   }

   if(return_value == LP_OPTIMAL || return_value == LP_D_OBJLIM ||
      return_value == LP_D_ITLIM){
      lp_data->x = (double *)malloc(DSIZE*lp_data->n);
      get_x(lp_data);
      for(i = lp_data->n; i>=0; i--){
	 colsol = lp_data->x[i];
	 if(lp_data->mip->is_int[i]){
	    if(colsol-floor(colsol) > env->par.lp_par.granularity &&
	       ceil(colsol)-colsol > env->par.lp_par.granularity){
	       break;
	    }
	 }
      }
      if(i<0){
	 node->node_status = NODE_STATUS__PRUNED;
	 node->feasibility_status = FEASIBLE_PRUNED;

	 if((env->warm_start->has_ub &&
	     lp_data->objval < env->warm_start->ub)||
	    !env->warm_start->has_ub){

	    if(!env->warm_start->has_ub){
	       env->warm_start->has_ub = TRUE;
	       best_sol->has_sol = TRUE;
	    }

	    env->warm_start->ub = best_sol->objval = lp_data->objval;

	    FREE(best_sol->xind);
	    FREE(best_sol->xval);

	    xind = (int *)malloc(ISIZE*lp_data->n);
	    xval = (double *)malloc(DSIZE*lp_data->n);

	    for (i = 0; i < lp_data->n; i++){
	       if (lp_data->x[i] > lpetol || lp_data->x[i] < -lpetol){
		  xind[cnt] = i;
		  xval[cnt++] = lp_data->x[i];
	       }
	    }
	    best_sol->xind = (int *)malloc(ISIZE*cnt);
	    best_sol->xval = (double *)malloc(DSIZE*cnt);
	    memcpy(best_sol->xind, xind, ISIZE*cnt);
	    memcpy(best_sol->xval, xval, DSIZE*cnt);
	    best_sol->xlength = cnt;
	    best_sol->xlevel = node->bc_level;
	    best_sol->xindex = node->bc_index;
	    best_sol->lpetol = lpetol;
	 }
      }
   }

   for(i = 0; i < 2*(level+1)+BB_BUNCH; i++){
      path[i]=0;
   }

   FREE(path);
   FREE(bpath);
   FREE(sense);
   FREE(rhs);
   FREE(matind);
   FREE(matbeg);
   FREE(matval);
   /* FIXME- for now just copy free_node_desc here. Decide where to carry
      resolve_node() */
   //   free_node_desc(&desc);
     if (desc){
      node_desc *n = desc;
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
      FREE(desc);
   }
   node->lower_bound = lp_data->objval;
   free_mip_desc(lp_data->mip);
   free_lp_arrays(lp_data);
   close_lp_solver(lp_data);
   FREE(xind);
   FREE(xval);
   FREE(lp_data);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int update_tree_bound(sym_environment *env, bc_node *root, int *cut_num, int *cuts_ind, char *cru_vars,
		      int change_type)
{

   int i, resolve = 0;
   char deletable = TRUE;

   if (root){

      check_trim_tree(env, root, cut_num, cuts_ind, change_type);

      if (root->node_status == NODE_STATUS__PRUNED ||
	  root->node_status == NODE_STATUS__INTERRUPTED ||
	  root->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
	  root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){
	 if(change_type == OBJ_COEFF_CHANGED || change_type == RHS_CHANGED ||
	    change_type == COL_BOUNDS_CHANGED ||  change_type == COLS_ADDED){
	   if (root->feasibility_status == FEASIBLE_PRUNED ||
	       root->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
	       root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){

	      check_better_solution(env, root, FALSE, change_type);
	   }

	   if (root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION &&
	       root->bobj.child_num > 0){
	      for(i = 0; i<root->bobj.child_num; i++){
		 if(!update_tree_bound(env, root->children[i], cut_num, cuts_ind, cru_vars, change_type)){
		    deletable = FALSE;
		 }
	      }
	      if(change_type == COL_BOUNDS_CHANGED && root->bobj.child_num > 0){
		 update_branching_decisions(env, root, change_type);
	      }

	   } else{
	      if(root->node_status == NODE_STATUS__WSPRUNED)
		 root->node_status = NODE_STATUS__PRUNED;
	      else
		 root->node_status = NODE_STATUS__WARM_STARTED;
	      if(resolve == 0){
		 root->lower_bound = -MAXDOUBLE;
	      }else{
		 resolve_node(env, root);
	      }
	   }
	   root->feasibility_status = 0; // or? ROOT_NODE;
	 }
      } else{
	 if(root->bobj.child_num > 0){
	    if (env->mip->var_type_modified){
	       if(!env->mip->is_int[root->children[0]->bobj.name]){
		  for(i = 0; i<root->bobj.child_num; i++){
		     ws_free_subtree(env, root->children[i], change_type, TRUE, TRUE);
		  }
		  root->bobj.child_num = 0;
		  root->node_status = NODE_STATUS__WARM_STARTED;
		  if(resolve == 0){
		     root->lower_bound = -MAXDOUBLE;
		  }else{
		     resolve_node(env, root);
		  }
	       }
	    } else {
	       if(change_type == COL_BOUNDS_CHANGED && root->bobj.child_num > 0){
		  update_branching_decisions(env, root, change_type);
	       }
	       for(i = 0; i<root->bobj.child_num; i++){
		  if(!update_tree_bound(env, root->children[i], cut_num, cuts_ind, cru_vars, change_type)){
		     deletable = FALSE;
		  }
	       }
	    }
	 }else{
	    if(root->node_status == NODE_STATUS__WSPRUNED)
	       root->node_status = NODE_STATUS__PRUNED;
	    else
	       root->node_status = NODE_STATUS__WARM_STARTED;
	    if(resolve == 0){
	       root->lower_bound = -MAXDOUBLE;
	    }else{
	       resolve_node(env, root);
	    }
	 }
      }

      /* should be before resolve!!!*/
      if(change_type == COLS_ADDED){
	 update_node_desc(env, root, change_type);
      }

      if(env->warm_start->trim_tree == ON_CRU_VARS){
	 if(deletable && root->bobj.child_num){
	    for(i = 0; i<root->bobj.child_num; i++){
	       ws_free_subtree(env, root->children[i], change_type, FALSE, TRUE);
	    }
	    root->node_status = NODE_STATUS__WARM_STARTED;
	    if(resolve == 0){
	       root->lower_bound = -MAXDOUBLE;
	    }else{
	       resolve_node(env, root);
	    }
	    root->bobj.child_num = 0;
	    if(root->bc_level){
	       if(cru_vars[root->parent->bobj.name]){
		  deletable = FALSE;
	       }
	    }
	 }
	 if(!deletable && root->bobj.child_num){
	    for(i = 0; i<root->bobj.child_num; i++){
	       register_cuts(root->children[i], cut_num, cuts_ind);
	    }
	 }
	 if(root->bc_level){
	    if(cru_vars[root->parent->bobj.name]){
	       deletable = FALSE;
	    }
	 }
      }
   }
   return deletable;
}

/*===========================================================================*/
/*===========================================================================*/


void register_cuts(bc_node *root, int *cut_num,  int *cuts_ind){

   int i, r_cnum = root->desc.cutind.size;
   int *c_list = root->desc.cutind.list;
   int c_ind;
   if(r_cnum > 0){
      for(i = 0; i < r_cnum; i++){
	 c_ind = c_list[i];
	 if(cuts_ind[c_ind] < 0){
	    cuts_ind[c_ind] = c_list[i] = *cut_num;
	    (*cut_num)++;
	 }else{
	    c_list[i] = cuts_ind[c_ind];
	 }
      }
   }
}

/*===========================================================================*/
/*===========================================================================*/

void update_node_desc(sym_environment *env, bc_node *root, int change_type){

   node_desc *desc = &(root->desc);
   int i, tmp_size, colnum = env->mip->n;
   int added_colnum = env->mip->new_col_num;
   int old_colnum = colnum - added_colnum;

   if(change_type == COLS_ADDED && added_colnum > 0){
      if(desc->uind.type == EXPLICIT_LIST){
	 tmp_size = desc->uind.size;
	 if(tmp_size < colnum){
	    desc->uind.list = (int *) realloc((int *)desc->uind.list,
					      (tmp_size + added_colnum) * ISIZE);
	    for(i = 0; i < added_colnum; i++){
	       desc->uind.list[tmp_size + i] = old_colnum + i;
	    }
	    desc->uind.size += added_colnum;
	 }
      }
      if(desc->basis.basis_exists &&
	 desc->basis.extravars.type == EXPLICIT_LIST){
	 tmp_size = desc->basis.extravars.size;
	 if(tmp_size < colnum){
	   /*  if(desc->basis.extravars.list){ */
/* 	       desc->basis.extravars.list =  */
/* 		  (int*)realloc((int*)desc->basis.extravars.list,  */
/* 				(tmp_size + added_colnum) *ISIZE); */
/* 	       for(i = 0; i < added_colnum; i++){ */
/* 		  desc->basis.extravars.list[tmp_size + i] = old_colnum + i; */
/* 	       } */
/* 	    } */
	    desc->basis.extravars.stat =
	       (int*)realloc((int*)desc->basis.extravars.stat,
			     (tmp_size + added_colnum) *ISIZE);
	    for(i = 0; i < added_colnum; i++){
	       desc->basis.extravars.stat[tmp_size + i] = VAR_AT_LB;
	    }
	    desc->basis.extravars.size += added_colnum;
	 }
      }
   }
}

/*===========================================================================*/
/*===========================================================================*/
void update_branching_decisions(sym_environment *env, bc_node *root, int change_type)
{
   int i, j;
   int type = 0;
   int resolve = 0;
   int deleted = 0;

   if(change_type == COL_BOUNDS_CHANGED && root->bobj.child_num > 0){
      for(i = 0; i<root->bobj.child_num; i++){
	 if(root->bobj.type == BRANCHING_VARIABLE){
	    switch (root->bobj.sense[i]){
	     case 'E':
		printf("error1-update_warm_start_tree\n");
		exit(0);
	     case 'L':
		if(type == 0){
		   if(root->bobj.rhs[i] < env->mip->lb[root->bobj.name]){
		      root->bobj.rhs[i] = ceil(env->mip->lb[root->bobj.name]);
		   }
		   else if(root->bobj.rhs[i] > env->mip->ub[root->bobj.name]){
		      root->bobj.rhs[i] = ceil(env->mip->ub[root->bobj.name]) - 1;
		   }
		}else {
		   if(root->bobj.rhs[i] < env->mip->lb[root->bobj.name]){
		      root->children[i]->node_status = NODE_STATUS__WSPRUNED;
		      root->children[i]->feasibility_status = INFEASIBLE_PRUNED;
		      for(j = 0; j<root->children[i]->bobj.child_num; j++){
			 ws_free_subtree(env, root->children[i]->children[j], change_type, FALSE, TRUE);
		      }
		      root->children[i]->bobj.child_num = 0;
		   }
		   else if(root->bobj.rhs[i] > env->mip->ub[root->bobj.name]){
		      root->bobj.rhs[i] = floor(env->mip->ub[root->bobj.name]);
		   }
		}
		break;
	     case 'G':

		if(type == 0){
		   if(root->bobj.rhs[i] > env->mip->ub[root->bobj.name]){
		      root->bobj.rhs[i] = floor(env->mip->ub[root->bobj.name]);
		   }
		   else if(root->bobj.rhs[i] < env->mip->lb[root->bobj.name]){
		      root->bobj.rhs[i] = floor(env->mip->lb[root->bobj.name]) + 1;
		   }
		}else {
		   if(root->bobj.rhs[i] > env->mip->ub[root->bobj.name]){
		      root->children[i]->node_status = NODE_STATUS__WSPRUNED;
		      root->children[i]->feasibility_status = INFEASIBLE_PRUNED;
		      for(j = 0; j<root->children[i]->bobj.child_num; j++){
			 ws_free_subtree(env, root->children[i]->children[j], change_type, FALSE, TRUE);
		      }
		      root->children[i]->bobj.child_num = 0;
		   }
		   else if(root->bobj.rhs[i] < env->mip->lb[root->bobj.name]){
		      root->bobj.rhs[i] = ceil(env->mip->lb[root->bobj.name]);
		   }
		}

		break;
	     case 'R':
		printf("error2-update_warm_start_tree\n");
		exit(0);
	    }
	 } else {
	    printf("error3-update_warm_start_tree\n");
	    exit(0);
	 }
      }
      root->bobj.child_num -= deleted;
      if(root->bobj.child_num <= 0){
	 if(root->node_status !=NODE_STATUS__WSPRUNED)
	    root->node_status = NODE_STATUS__WARM_STARTED;
	 if(resolve == 0){
	    root->lower_bound = MAXDOUBLE;
	 }else{
	    resolve_node(env, root);
	 }
      }
   }
}

/*===========================================================================*/
/*===========================================================================*/
void check_trim_tree(sym_environment *env, bc_node *root, int *cut_num,
		     int *cuts_ind, int change_type)
{
   int i;
   char trim_type = env->warm_start->trim_tree;
   int level = env->warm_start->trim_tree_level;
   int index = env->warm_start->trim_tree_index;
   problem_stat *stat = &(env->warm_start->stat);
   int trim_subtree = FALSE;

   if(trim_type){
      register_cuts(root, cut_num, cuts_ind);
   }

   if(trim_type == TRIM_INDEX){
      if(root->bobj.child_num > 0){
	 for (i = 0; i < root->bobj.child_num; i++){
	    if(root->children[i]->bc_index <= index){
	       break;
	    }
	 }
	 if(i >= root->bobj.child_num){
	    trim_subtree = TRUE;
	 }
      }
   }else if(trim_type == TRIM_LEVEL){
      if(root->bc_level >= level){
	 trim_subtree = TRUE;
      }
   }

   if(trim_subtree && root->bobj.child_num){
      for (i = 0; i < root->bobj.child_num; i++){
	 ws_free_subtree(env, root->children[i], change_type, TRUE, FALSE);
      }
      root->bobj.child_num = 0;
   }else{
      for (i = 0; i < root->bobj.child_num; i++){
	 root->children[i]->bc_index = stat->tree_size++;
	 stat->created++;
      }
   }

   if (root->node_status == NODE_BRANCHED_ON && root->bobj.child_num){
      stat->analyzed++;
   }
}
/*===========================================================================*/
/*===========================================================================*/

void cut_ws_tree_index(sym_environment *env, bc_node *root, int index,
		       problem_stat *stat, int change_type)
{

  int i;

  if (root){
     if (root->node_status != NODE_STATUS__CANDIDATE){
	stat->analyzed++;
     }
     if(root->bobj.child_num > 0){
	for (i = 0; i < root->bobj.child_num; i++){
	   if(root->children[i]->bc_index <= index){
	      break;
	   }
	}
	if(i < root->bobj.child_num){

	   for (i = 0; i < root->bobj.child_num; i++){
	      root->children[i]->bc_index = ++stat->tree_size;
	      stat->created++;
	   }
	   for (i = root->bobj.child_num - 1; i >= 0; i--){
	      cut_ws_tree_index(env, root->children[i], index, stat,
				change_type);
	   }
	} else{
	   for (i = root->bobj.child_num - 1; i >= 0; i--)
	      ws_free_subtree(env, root->children[i], change_type, TRUE, FALSE);
	   root->bobj.child_num = 0;

	   if (root->node_status == NODE_STATUS__BRANCHED_ON){
	      root->node_status = NODE_STATUS__WARM_STARTED;
	   }
	}
     }
  }
}

/*===========================================================================*/
/*===========================================================================*/

void cut_ws_tree_level(sym_environment *env, bc_node *root, int level,
		       problem_stat *stat, int change_type)
{

   int i;
   if(root){
      if (root->node_status != NODE_STATUS__CANDIDATE){
	 stat->analyzed++;
      }

      if(root->bc_level < level){
	 for (i = 0; i < root->bobj.child_num; i++){
	    root->children[i]->bc_index = ++stat->tree_size;
	    stat->created++;
	 }
	 for (i = root->bobj.child_num - 1; i >= 0; i--){
	    cut_ws_tree_level(env, root->children[i], level, stat, change_type);
	 }
      }
      if(root->bc_level == level){
	 for (i = root->bobj.child_num - 1; i >= 0; i--)
	    ws_free_subtree(env, root->children[i], change_type, TRUE, FALSE);
	 root->bobj.child_num = 0;

	 if (root->node_status == NODE_STATUS__BRANCHED_ON){
	    root->node_status = NODE_STATUS__WARM_STARTED;
	 }
      }
   }
}

/*===========================================================================*/
/*===========================================================================*/

void ws_free_subtree(sym_environment *env, bc_node *root, int change_type, int check_solution,
		     int update_stats)
{

   int i;

   if (root == NULL) return;

   if(check_solution){
      if (root->feasibility_status == FEASIBLE_PRUNED ||
	  root->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
	  root->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){
	 check_better_solution(env, root, TRUE, change_type);
      }
   }

   for (i = root->bobj.child_num - 1; i >= 0; i--){
      ws_free_subtree(env, root->children[i], change_type, check_solution, update_stats);
   }

   if(update_stats){
      env->warm_start->stat.analyzed--;
      env->warm_start->stat.created--;
      env->warm_start->stat.tree_size--;
   }

   free_tree_node(root);
}

/*===========================================================================*/
/*===========================================================================*/
void check_better_solution(sym_environment * env, bc_node *root, int delete_node,
			   int change_type){

   int i, j, new_solution = FALSE;
   double upper_bound = 0.0, lpetol = 9.9999999999999995e-07, valuesi=0;
   double obj[2] = {0.0, 0.0}, gamma, tau, objval;
   int *matbeg, *matind;
   double *rowact = NULL, *matval, *colsol = NULL;
   int feasible = TRUE;


#ifdef USE_SYM_APPLICATION
   /* for now, just assume it is infeasible */

#else

   MIPdesc *mip = env->mip;
   lp_sol *best_sol = &(env->warm_start->best_sol);


   if(env->mip->var_type_modified == TRUE) {
      for(i = root->sol_size; i >=0; i++){
	 if(mip->is_int[root->sol_ind[i]]){
	    valuesi = root->sol[i];
	    if (valuesi-floor(valuesi) > lpetol &&
		ceil(valuesi)-valuesi > lpetol){
	       break;
	    }
	 }
      }
      feasible = i < 0 ? TRUE : FALSE;
   }

   if(feasible){
      if(change_type == OBJ_COEFF_CHANGED || change_type == COLS_ADDED){
	 if(!env->par.mc_warm_start){
	    for(i = 0; i<root->sol_size; i++){
	       upper_bound += mip->obj[root->sol_ind[i]] * root->sol[i];
	    }
	    if((env->warm_start->has_ub &&
		upper_bound<env->warm_start->ub)||
	       !env->warm_start->has_ub){

	       if(!env->warm_start->has_ub){
		  env->warm_start->has_ub = TRUE;
		  best_sol->has_sol = TRUE;
	       }

	       env->warm_start->ub = upper_bound;
	       best_sol->objval = upper_bound;
	       new_solution = true;
	    }

	 } else {

	    gamma = env->par.lp_par.mc_gamma;
	    tau = env->par.lp_par.mc_tau;

	    for(i = 0; i<root->sol_size; i++){
	       obj[0] += mip->obj1[root->sol_ind[i]]*root->sol[i];
	       obj[1] += mip->obj2[root->sol_ind[i]]*root->sol[i];
	    }

	    if(gamma != -1.0){
	    objval = gamma*obj[0] + tau*obj[1];
	    } else{
	       objval = tau*obj[1];
	    }
	    if ((env->has_mc_ub &&  (objval < env->mc_ub - lpetol -
				     MAX(0, MIN(gamma, tau)))) ||
		!env->has_mc_ub){

	       if(!env->has_mc_ub){
	       env->has_mc_ub = TRUE;
	       env->warm_start->has_ub = TRUE;
	       best_sol->has_sol = TRUE;
	       }

	       env->mc_ub = env->warm_start->ub = best_sol->objval = objval;
	       env->obj[0] = obj[0];
	       env->obj[1] = obj[1];
	       new_solution = TRUE;
	    }
	 }
      }else if(change_type == RHS_CHANGED){

	 colsol = (double*)calloc(env->mip->n, DSIZE);
	 rowact = (double*) calloc(env->mip->m, DSIZE);
	 matbeg = env->mip->matbeg;
	 matval = env->mip->matval;
	 matind = env->mip->matind;

	 for(i = 0; i<root->sol_size; i++){
	    upper_bound += mip->obj[root->sol_ind[i]] * root->sol[i];
	    colsol[root->sol_ind[i]] = root->sol[i];
	 }

	 for(i = 0; i < env->mip->n; i++){
	    for(j = matbeg[i]; j<matbeg[i+1]; j++){
	       rowact[matind[j]] += matval[j] * colsol[i];
	    }
	 }

	 for(i = 0; i < env->mip->m; i++){
	    switch(env->mip->sense[i]){
	     case 'L':
		if (rowact[i] > env->mip->rhs[i] + lpetol)
		   feasible = FALSE;
		break;
	     case 'G':
		if (rowact[i] < env->mip->rhs[i] - lpetol)
		   feasible = FALSE;
		break;
	     case 'E':
		if (!((rowact[i] > env->mip->rhs[i] - lpetol) &&
		      (rowact[i] < env->mip->rhs[i] + lpetol)))
		   feasible = FALSE;
		break;
	     case 'R':
		if (rowact[i] > env->mip->rhs[i] + lpetol ||
		    rowact[i] < env->mip->rhs[i] - env->mip->rngval[i] - lpetol)
		   feasible = FALSE;
		break;
	     case 'N':
	     default:
		break;
	    }

	    if (!feasible)
	       break;
	 }

	 if(feasible) {
	    if((env->warm_start->has_ub &&
		upper_bound<env->warm_start->ub)||
	    !env->warm_start->has_ub){

	       if(!env->warm_start->has_ub){
		  env->warm_start->has_ub = TRUE;
		  best_sol->has_sol = TRUE;
	       }

	    env->warm_start->ub = upper_bound;
	    best_sol->objval = upper_bound;
	    new_solution = true;
	    }
	 }
      }else if(change_type == COL_BOUNDS_CHANGED){

	 colsol = (double*)calloc(env->mip->n, DSIZE);

	 for(i = 0; i<root->sol_size; i++){
	    upper_bound += mip->obj[root->sol_ind[i]] * root->sol[i];
	    colsol[root->sol_ind[i]] = root->sol[i];
	 }

	 for(i = 0; i < env->mip->n; i++){
	    if ((colsol[i] < env->mip->lb[i] - lpetol) ||
		(colsol[i] > env->mip->ub[i] + lpetol))
	       feasible = FALSE;
	 }

	 if(feasible) {
	    if((env->warm_start->has_ub &&
		upper_bound < env->warm_start->ub)||
	       !env->warm_start->has_ub){

	       if(!env->warm_start->has_ub){
		  env->warm_start->has_ub = TRUE;
		  best_sol->has_sol = TRUE;
	       }

	       env->warm_start->ub = upper_bound;
	       best_sol->objval = upper_bound;
	       new_solution = true;
	    }
	 }
      }
   }

   if(new_solution){

      best_sol->xlevel = root->bc_level;
      best_sol->xindex = root->bc_index;
      best_sol->xlength = root->sol_size;
      best_sol->lpetol = lpetol;

      FREE(best_sol->xind);
      FREE(best_sol->xval);

      //      if(delete_node){
      best_sol->xind = root->sol_ind;
      best_sol->xval = root->sol;
      root->sol_ind = 0;
      root->sol = 0;
	 //      } else {
	 //	 best_sol->xind = (int *)malloc(ISIZE*root->sol_size);
	 //	 best_sol->xval = (double *)malloc(DSIZE*root->sol_size);
	 //	 memcpy(best_sol->xind, root->sol_ind, ISIZE*root->sol_size);
	 //	 memcpy(best_sol->xval, root->sol, DSIZE* root->sol_size);
      //}
   }
   FREE(rowact);
   FREE(colsol);
   FREE(root->sol);
   FREE(root->sol_ind);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int copy_node(bc_node * n_to, bc_node *n_from)
{

   if (!n_to || !n_from){
      printf("copy_node(): Empty node_structure(s)!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);;
   }

   n_to->bc_index = n_from->bc_index;
   n_to->bc_level = n_from->bc_level;

   n_to->lp = n_from->lp;
   n_to->cg = n_from->cg;
   n_to->cp = n_from->cp;


   n_to->lower_bound = n_from->lower_bound;
   n_to->opt_estimate = n_from->opt_estimate;
   n_to->node_status = n_from->node_status;
   n_to->feasibility_status = n_from->feasibility_status;
   n_to->sol_size = n_from->sol_size;

   if (n_from->feasibility_status == FEASIBLE_PRUNED ||
       n_from->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
       n_from->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){
      if(n_from->sol){
	 n_to->sol = (double *)malloc(n_from->sol_size*DSIZE);
	 n_to->sol_ind = (int *)malloc(n_from->sol_size*ISIZE);
	 memcpy(n_to->sol, n_from->sol, n_from->sol_size * DSIZE);
	 memcpy(n_to->sol_ind, n_from->sol_ind, n_from->sol_size * ISIZE);
      }
   }

#ifdef TRACE_PATH
   n_to->optimal_path = n_from->optimal_path;
#endif

   n_to->bobj = n_from->bobj;

#if defined (COMPILING_FOR_LP) || defined(COMPILE_IN_LP)

   //FIXME, Do we need this while writing to file

#if 0
   if (n_from->bobj.row){
      n_to->bobj.row = (waiting_row*) malloc(sizeof(waiting_row));
      memcpy(n_to->bobj.row, n_from->bobj.row, sizeof(waiting_row));

      n_to->bobj.row->matind = (int*)malloc(sizeof(int)*n_to->bobj.row->nzcnt);
      n_to->bobj.row->matval =
	 (double*)malloc(sizeof(double)*n_to->bobj.row->nzcnt);

      memcpy(n_to->bobj.row->matind, n_from->bobj.row->matind,
	     ISIZE*n_to->bobj.row->nzcnt);
      memcpy(n_to->bobj.row->matval, n_from->bobj.row->matval,
	     DSIZE*n_to->bobj.row->nzcnt);

      n_to->bobj.row->cut = (cut_data*)malloc( sizeof(cut_data));
      memcpy(n_to->bobj.row->cut, n_from->bobj.row->cut, sizeof(cut_data));

      n_to->bobj.row->cut->coef =
	 (char*)malloc(sizeof(char)*n_to->bobj.row->cut->size);
      memcpy(n_to->bobj.row->cut->coef, n_from->bobj.row->cut->coef,
	     CSIZE*n_to->bobj.row->cut->size);
   }
#endif

#endif

#ifndef MAX_CHILDREN_NUM
   n_to->bobj.sense = (char*)malloc(n_to->bobj.child_num*CSIZE);
   n_to->bobj.rhs = (double *) malloc(n_to->bobj.child_num*DSIZE);
   n_to->bobj.range = (double *) malloc(n_to->bobj.child_num*DSIZE);
   n_to->bobj.branch = (int *) malloc(n_to->bobj.child_num*ISIZE);

   n_to->bobj.objval = (double*)malloc(n_to->bobj.child_num*DSIZE);
   n_to->bobj.termcode = (int *) malloc(n_to->bobj.child_num*ISIZE);
   n_to->bobj.iterd = (int *) malloc(n_to->bobj.child_num*ISIZE);
   n_to->bobj.feasible = (int *) malloc(n_to->bobj.child_num*ISIZE);

#endif

   memcpy(n_to->bobj.sense, n_from->bobj.sense,
	  n_to->bobj.child_num*CSIZE);
   memcpy(n_to->bobj.rhs, n_from->bobj.rhs,
	  n_to->bobj.child_num*DSIZE);
   memcpy(n_to->bobj.range, n_from->bobj.range,
	  n_to->bobj.child_num*DSIZE);
   memcpy(n_to->bobj.branch, n_from->bobj.branch,
	  n_to->bobj.child_num*ISIZE);

#if 0
   /* We shouldn't need to copy this. It's only used in the LP */
   memcpy(n_to->bobj.objval, n_from->bobj.objval,
	  n_to->bobj.child_num*DSIZE);
   memcpy(n_to->bobj.termcode, n_from->bobj.termcode,
	  n_to->bobj.child_num*ISIZE);
   memcpy(n_to->bobj.iterd, n_from->bobj.iterd,
	  n_to->bobj.child_num*ISIZE);
   memcpy(n_to->bobj.feasible, n_from->bobj.feasible,
	  n_to->bobj.child_num*ISIZE);
#endif

   n_to->desc = n_from->desc;

   if (n_to->desc.uind.size){
      n_to->desc.uind.list = (int *) malloc(n_to->desc.uind.size*ISIZE);
      memcpy( n_to->desc.uind.list,  n_from->desc.uind.list,
	      n_to->desc.uind.size*ISIZE);
   }
   if(n_to->desc.basis.basis_exists){
      if (n_to->desc.basis.basevars.size){
	 n_to->desc.basis.basevars.stat =
	    (int *) malloc(n_to->desc.basis.basevars.size*ISIZE);
	 memcpy( n_to->desc.basis.basevars.stat,
		 n_from->desc.basis.basevars.stat,
		 n_to->desc.basis.basevars.size*ISIZE);
	 if (n_to->desc.basis.basevars.type == WRT_PARENT){
	    n_to->desc.basis.basevars.list =
	       (int *) malloc(n_to->desc.basis.basevars.size*ISIZE);
	    memcpy( n_to->desc.basis.basevars.list,
		    n_from->desc.basis.basevars.list,
		    n_to->desc.basis.basevars.size*ISIZE);
	 }
      }

      if (n_to->desc.basis.extravars.size){
	 n_to->desc.basis.extravars.stat =
	    (int *) malloc(n_to->desc.basis.extravars.size*ISIZE);
	 memcpy( n_to->desc.basis.extravars.stat,
		 n_from->desc.basis.extravars.stat,
		 n_to->desc.basis.extravars.size*ISIZE);
	 if (n_to->desc.basis.extravars.type == WRT_PARENT){
	    n_to->desc.basis.extravars.list =
	       (int *) malloc(n_to->desc.basis.extravars.size*ISIZE);
	    memcpy( n_to->desc.basis.extravars.list,
		    n_from->desc.basis.extravars.list,
		    n_to->desc.basis.extravars.size*ISIZE);
	 }
      }

      if (n_to->desc.basis.baserows.size){
	 n_to->desc.basis.baserows.stat =
	    (int *) malloc(n_to->desc.basis.baserows.size*ISIZE);
	 memcpy( n_to->desc.basis.baserows.stat,
		 n_from->desc.basis.baserows.stat,
		 n_to->desc.basis.baserows.size*ISIZE);
	 if (n_to->desc.basis.baserows.type == WRT_PARENT){
	    n_to->desc.basis.baserows.list =
	       (int *) malloc(n_to->desc.basis.baserows.size*ISIZE);
	    memcpy( n_to->desc.basis.baserows.list,
		    n_from->desc.basis.baserows.list,
		    n_to->desc.basis.baserows.size*ISIZE);
	 }
      }

      if (n_to->desc.basis.extrarows.size){
	 n_to->desc.basis.extrarows.stat =
	    (int *) malloc(n_to->desc.basis.extrarows.size*ISIZE);
	 memcpy( n_to->desc.basis.extrarows.stat,
		 n_from->desc.basis.extrarows.stat,
		 n_to->desc.basis.extrarows.size*ISIZE);
	 if (n_to->desc.basis.extrarows.type == WRT_PARENT){
	    n_to->desc.basis.extrarows.list =
	       (int *) malloc(n_to->desc.basis.extrarows.size*ISIZE);
	    memcpy( n_to->desc.basis.extrarows.list,
		    n_from->desc.basis.extrarows.list,
		    n_to->desc.basis.extrarows.size*ISIZE);
	 }
      }
   }
   if (n_to->desc.not_fixed.size){
	 n_to->desc.not_fixed.list =
	 (int *) malloc(n_to->desc.not_fixed.size*ISIZE);
      memcpy( n_to->desc.not_fixed.list,  n_from->desc.not_fixed.list,
	      n_to->desc.not_fixed.size*ISIZE);
   }

   if (n_to->desc.cutind.size){
      n_to->desc.cutind.list = (int *) malloc(n_to->desc.cutind.size*ISIZE);
      memcpy( n_to->desc.cutind.list,  n_from->desc.cutind.list,
	      n_to->desc.cutind.size*ISIZE);
   }

   if (n_to->desc.desc_size){
      n_to->desc.desc = (char*) malloc(n_to->desc.desc_size*CSIZE);
      memcpy(n_to->desc.desc, n_from->desc.desc,
	     n_to->desc.desc_size*CSIZE);
   }

   if(n_to->desc.bnd_change){
      n_to->desc.bnd_change = (bounds_change_desc *)
	 calloc(sizeof(bounds_change_desc),1);
      if(n_from->desc.bnd_change->num_changes){
	 n_to->desc.bnd_change->index =(int *)malloc
	    (ISIZE *n_from->desc.bnd_change->num_changes);
	 n_to->desc.bnd_change->lbub =(char *)malloc
	    (CSIZE *n_from->desc.bnd_change->num_changes);
	 n_to->desc.bnd_change->value =(double *)malloc
	    (DSIZE *n_from->desc.bnd_change->num_changes);
	 memcpy( n_to->desc.bnd_change->index,
		 n_from->desc.bnd_change->index,
		 ISIZE *n_from->desc.bnd_change->num_changes);
	 memcpy( n_to->desc.bnd_change->lbub,
		 n_from->desc.bnd_change->lbub,
		 CSIZE *n_from->desc.bnd_change->num_changes);
	 memcpy( n_to->desc.bnd_change->value,
		 n_from->desc.bnd_change->value,
		 DSIZE *n_from->desc.bnd_change->num_changes);
      }
      n_to->desc.bnd_change->num_changes =
	 n_from->desc.bnd_change->num_changes;
   }


   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int copy_tree(bc_node *root_to, bc_node *root_from)
{
   int i, childNum;

   if (!root_to || !root_from){
      printf("copy_tree(): Empty root node(s)!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   if (root_from){
      copy_node(root_to, root_from);
      childNum = root_to->bobj.child_num;
      if (childNum) {
	 root_to->children = (bc_node **) calloc(sizeof(bc_node*), childNum);
	 for (i = 0; i < childNum; i++){
	    root_to->children[i] = (bc_node *) calloc(1, sizeof(bc_node));
	    root_to->children[i]->parent = root_to;
	    copy_tree(root_to->children[i], root_from->children[i]);
	 }
      }
   }
   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int write_node(bc_node *node, FILE*f)
{
   int i;

   if (!node){
      printf("write_node(): Empty node!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   fprintf(f,"\n\n");

   fprintf(f," NODE_INDEX      : %i\n",node->bc_index);
   fprintf(f," NODE_LEVEL      : %i\n",node->bc_level);
   fprintf(f," LOWER_BOUND     : %.4f\n",node->lower_bound);
   fprintf(f," NODE_STATUS     : %i\n",(int)node->node_status);
   fprintf(f," NODE_LP         : %i\n",node->lp);
   fprintf(f," NODE_CG         : %i\n",node->cg);
   fprintf(f," NODE_CP         : %i\n",node->cp);
   fprintf(f," OPT_ESTIMATE    : %.4f\n",node->opt_estimate);


#ifdef TRACE_PATH
   fprintf(f," OPTIMAL_PATH    : %c\n",node->optimal_path);
#endif
   if (node->parent){
      fprintf(f," PARENT_INDEX    : %i\n",node->parent->bc_index);
   }
   else{
      fprintf(f," PARENT_INDEX    : -1\n");
   }

   fprintf(f, "CHILDREN:        %i %i %i\n", (int)node->bobj.type,
           node->bobj.name, node->bobj.child_num);
   for (i = 0; i < node->bobj.child_num; i++){
      fprintf(f, "%i %c %f %f %i\n", node->children[i]->bc_index,
              node->bobj.sense[i], node->bobj.rhs[i],
              node->bobj.range[i], node->bobj.branch[i]);
   }

   fprintf(f," NODE_DESCRIPTION                 : %i\n",node->desc.nf_status);
   fprintf(f," USER_INDICES(Type,Size,Added)    : %i %i %i\n",
	   (int)node->desc.uind.type, node->desc.uind.size,
	   node->desc.uind.added);

   for (i = 0; i < node->desc.uind.size; i++){
      fprintf(f," %i", node->desc.uind.list[i]);
   }
   fprintf(f, "\n");

   fprintf(f," NOT_FIXED(Type,Size,Added)   : %i %i %i\n",
	   (int)node->desc.not_fixed.type, node->desc.not_fixed.size,
	   node->desc.not_fixed.added);

   for (i = 0; i < node->desc.not_fixed.size; i++){
      fprintf(f," %i", node->desc.not_fixed.list[i]);
   }

   fprintf(f, "\n");

   fprintf(f," CUT_INDICES(Type,Size,Added)   : %i %i %i\n",
	   (int)node->desc.cutind.type, node->desc.cutind.size,
	   node->desc.cutind.added);

   for (i = 0; i < node->desc.cutind.size; i++){
      fprintf(f," %i", node->desc.cutind.list[i]);
   }
   fprintf(f, "\n");

   fprintf(f," BASIS          : %i\n", (int)node->desc.basis.basis_exists);
   fprintf(f," BASE_VARIABLES : %i %i\n",(int)node->desc.basis.basevars.type,
      node->desc.basis.basevars.size);
   if (node->desc.basis.basevars.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.basevars.size; i++){
	 fprintf(f," %i %i", node->desc.basis.basevars.list[i],
		 node->desc.basis.basevars.stat[i]);
      }
   }
   else{
      for (i = 0; i < node->desc.basis.basevars.size; i++){
	 fprintf(f," %i", node->desc.basis.basevars.stat[i]);
      }
   }
   fprintf(f,"\n");

   fprintf(f," EXTRA_VARIABLES : %i %i\n",
	   (int)node->desc.basis.extravars.type,
	   node->desc.basis.extravars.size);
   if (node->desc.basis.extravars.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.extravars.size; i++){
	 fprintf(f," %i %i", node->desc.basis.extravars.list[i],
	    node->desc.basis.extravars.stat[i]);
      }
   }
   else{
      for (i = 0; i < node->desc.basis.extravars.size; i++){
	 fprintf(f," %i", node->desc.basis.extravars.stat[i]);
      }
   }

   fprintf(f,"\n");
   fprintf(f," BASE_ROWS      : %i %i\n",
	   (int)node->desc.basis.baserows.type,
	   node->desc.basis.baserows.size);
   if (node->desc.basis.baserows.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.baserows.size; i++){
	 fprintf(f," %i %i", node->desc.basis.baserows.list[i],
		 node->desc.basis.baserows.stat[i]);
      }
   }
   else{
      for (i = 0; i < node->desc.basis.baserows.size; i++){
	 fprintf(f," %i", node->desc.basis.baserows.stat[i]);
      }
   }
   fprintf(f,"\n");

   fprintf(f," EXTRA_ROWS       : %i %i\n",
	   (int)node->desc.basis.extrarows.type,
	   node->desc.basis.extrarows.size);
   if (node->desc.basis.extrarows.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.extrarows.size; i++){
	 fprintf(f," %i %i", node->desc.basis.extrarows.list[i],
	    node->desc.basis.extrarows.stat[i]);
      }
   }
   else{
      for (i = 0; i < node->desc.basis.extrarows.size; i++){
	 fprintf(f," %i", node->desc.basis.extrarows.stat[i]);
      }
   }

   fprintf(f,"\n");
   fprintf(f," USER_DESC_SIZE_&_ELEMENTS       : %i\n",
	   node->desc.desc_size);
   for(i = 0; i<node->desc.desc_size;i++){
      fprintf(f," %i", (int)node->desc.desc[i]);
   }

   fprintf(f,"\n");

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int write_tree(bc_node *root, FILE *f)
{
   int i;
   if (!root){
      printf("write_tree(): Empty root node!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   write_node(root, f);

   for(i=0; i<root->bobj.child_num; i++){
      write_tree(root->children[i], f);
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int read_node(bc_node * node, FILE * f)
{
   char str[80];
   int i=0, num=0, ch=0;

   if (!node || !f){
      printf("read_node(): Empty node or unable to read from file!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   fscanf(f,"%s %s %i", str, str, &node->bc_index);
   fscanf(f,"%s %s %i", str, str, &node->bc_level);
   fscanf(f,"%s %s %lf", str, str, &node->lower_bound);
   fscanf(f,"%s %s %i", str, str, &ch);
   node->node_status = (char)ch;
   fscanf(f,"%s %s %i", str, str, &node->lp);
   fscanf(f,"%s %s %i", str, str, &node->cg);
   fscanf(f,"%s %s %i", str, str, &node->cp);
   fscanf(f,"%s %s %lf", str, str, &node->opt_estimate);

#ifdef TRACE_PATH
   fscanf(f,"%s %s %c", str, str, &node->optimal_path);
#endif
   fscanf(f,"%s %s %i", str, str, &num);
   fscanf(f,"%s %s %i %i %i", str, str, &ch, &node->bobj.name,
	  &node->bobj.child_num);
   node->bobj.type = (char)ch;
   if (node->bobj.child_num){
#ifndef MAX_CHILDREN_NUM
      node->bobj.sense = (char*)malloc(node->bobj.child_num*CSIZE);
      node->bobj.rhs = (double *) malloc(node->bobj.child_num*DSIZE);
      node->bobj.range = (double *) malloc(node->bobj.child_num*DSIZE);
      node->bobj.branch = (int *) malloc(node->bobj.child_num*ISIZE);
#endif
      for(i = 0; i < node->bobj.child_num; i++){
	 fscanf(f,"%i %c %lf %lf %i", &num, &node->bobj.sense[i],
		&node->bobj.rhs[i], &node->bobj.range[i],
		&node->bobj.branch[i]);
      }
   }

   fscanf(f,"%s %s %i", str, str, &node->desc.nf_status);
   fscanf(f,"%s %s %i %i %i", str, str, &ch, &node->desc.uind.size,
	  &node->desc.uind.added);
   node->desc.uind.type =  (char)ch;

   if (node->desc.uind.size){
      node->desc.uind.list = (int *) malloc(node->desc.uind.size*ISIZE);
      for (i = 0; i < node->desc.uind.size; i++){
	 fscanf(f, "%i", &node->desc.uind.list[i]);
      }
   }

   fscanf(f,"%s %s %i %i %i", str, str, &ch, &node->desc.not_fixed.size,
	  &node->desc.not_fixed.added);
   node->desc.not_fixed.type = (char)ch;

   if (node->desc.not_fixed.size){
      node->desc.not_fixed.list =
	 (int *) malloc(node->desc.not_fixed.size*ISIZE);
      for (i = 0; i < node->desc.not_fixed.size; i++){
	 fscanf(f, "%i", &node->desc.not_fixed.list[i]);
      }
   }

   fscanf(f,"%s %s %i %i %i", str, str, &ch, &node->desc.cutind.size,
	  &node->desc.cutind.added);
   node->desc.cutind.type = (char) ch;
   if (node->desc.cutind.size){
      node->desc.cutind.list = (int *) malloc(node->desc.cutind.size*ISIZE);
      for (i = 0; i < node->desc.cutind.size; i++){
	 fscanf(f, "%i", &node->desc.cutind.list[i]);
      }   }

   fscanf(f,"%s %s %i", str, str, &ch);
   node->desc.basis.basis_exists = (char)ch;
   fscanf(f,"%s %s %i %i", str, str, &ch, &node->desc.basis.basevars.size);
   node->desc.basis.basevars.type = (char)ch;
   if (node->desc.basis.basevars.size){
      node->desc.basis.basevars.stat =
	 (int *) malloc(node->desc.basis.basevars.size*ISIZE);
      if (node->desc.basis.basevars.type == WRT_PARENT){
	 node->desc.basis.basevars.list =
	    (int *) malloc(node->desc.basis.basevars.size*ISIZE);
	 for (i = 0; i < node->desc.basis.basevars.size; i++){
	    fscanf(f, "%i %i", &node->desc.basis.basevars.list[i],
	       &node->desc.basis.basevars.stat[i]);
	 }
      }
      else{
	 for (i = 0; i < node->desc.basis.basevars.size; i++)
	    fscanf(f, "%i", &node->desc.basis.basevars.stat[i]);
      }
   }

   fscanf(f,"%s %s %i %i", str, str, &ch, &node->desc.basis.extravars.size);
   node->desc.basis.extravars.type = (char)ch;
   if (node->desc.basis.extravars.size){
      node->desc.basis.extravars.stat =
	 (int *) malloc(node->desc.basis.extravars.size*ISIZE);
      if (node->desc.basis.extravars.type == WRT_PARENT){
	 node->desc.basis.extravars.list =
	    (int *) malloc(node->desc.basis.extravars.size*ISIZE);
	 for (i = 0; i < node->desc.basis.extravars.size; i++){
	    fscanf(f, "%i %i", &node->desc.basis.extravars.list[i],
		   &node->desc.basis.extravars.stat[i]);
	 }
      }else{
	 for (i = 0; i < node->desc.basis.extravars.size; i++)
	    fscanf(f, "%i", &node->desc.basis.extravars.stat[i]);
      }
   }

   fscanf(f,"%s %s %i %i", str, str, &ch, &node->desc.basis.baserows.size);
   node->desc.basis.baserows.type = (char)ch;
   if (node->desc.basis.baserows.size){
      node->desc.basis.baserows.stat =
	 (int *) malloc(node->desc.basis.baserows.size*ISIZE);
      if (node->desc.basis.baserows.type == WRT_PARENT){
	 node->desc.basis.baserows.list =
	    (int *) malloc(node->desc.basis.baserows.size*ISIZE);
	 for (i = 0; i < node->desc.basis.baserows.size; i++)
	    fscanf(f, "%i %i", &node->desc.basis.baserows.list[i],
		   &node->desc.basis.baserows.stat[i]);
      }else{
	 for (i = 0; i < node->desc.basis.baserows.size; i++)
	    fscanf(f, "%i", &node->desc.basis.baserows.stat[i]);
      }
   }

   fscanf(f,"%s %s %i %i", str, str, &ch, &node->desc.basis.extrarows.size);
   node->desc.basis.extrarows.type = (char)ch;
   if (node->desc.basis.extrarows.size){
      node->desc.basis.extrarows.stat =
	 (int *) malloc(node->desc.basis.extrarows.size*ISIZE);
      if (node->desc.basis.extrarows.type == WRT_PARENT){
	 node->desc.basis.extrarows.list =
	    (int *) malloc(node->desc.basis.extrarows.size*ISIZE);
	 for (i = 0; i < node->desc.basis.extrarows.size; i++)
	    fscanf(f, "%i %i", &node->desc.basis.extrarows.list[i],
		   &node->desc.basis.extrarows.stat[i]);
      }else{
	 for (i = 0; i < node->desc.basis.extrarows.size; i++)
	    fscanf(f, "%i", &node->desc.basis.extrarows.stat[i]);
      }
   }

   fscanf(f,"%s %s %i", str, str, &node->desc.desc_size);
   if (node->desc.desc_size){
      node->desc.desc =
	 (char *) malloc(node->desc.desc_size*CSIZE);
      for(i = 0; i<node->desc.desc_size; i++){
	 fscanf(f, "%i", &ch);
	 node->desc.desc[i] = (char)ch;
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int read_tree(bc_node * root, FILE *f)
{
   if (!root || !f){
      printf("read_tree(): Empty node or unable to write!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   read_node(root, f);

   int i, childNum = root->bobj.child_num;

   if (childNum!=0){
      root->children = (bc_node **) malloc(sizeof(bc_node*)*childNum);
      for (i = 0; i < childNum; i++){
	 root->children[i] = (bc_node *) calloc(1,sizeof(bc_node));
	 root->children[i]->parent = root;
	 read_tree(root->children[i], f);
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int set_param(sym_environment *env, char *line)
{
   int i;
   char key[MAX_LINE_LENGTH +1], value[MAX_LINE_LENGTH +1];
   double timeout;
   str_int colgen_str[COLGEN_STR_SIZE] = COLGEN_STR_ARRAY;
   str_int compare_can_str[COMPARE_CAN_STR_SIZE] = COMPARE_CAN_STR_ARRAY;
   tm_params *tm_par = &env->par.tm_par;
   lp_params *lp_par = &env->par.lp_par;
   cg_params *cg_par = &env->par.cg_par;
   cp_params *cp_par = &env->par.cp_par;
   dg_params *dg_par = &env->par.dg_par;
   prep_params *prep_par = &env->par.prep_par;

   strcpy(key,"");
   sscanf(line,"%s%s", key, value);

   /***********************************************************************
    ***                    Global parameters                            ***
    ***********************************************************************/
   if (strcmp(key, "verbosity") == 0){
      READ_INT_PAR(env->par.verbosity);
      tm_par->verbosity = lp_par->verbosity = cg_par->verbosity =
	 cp_par->verbosity = env->par.verbosity;
      prep_par->verbosity = env->par.verbosity;
      return(0);
   }
   else if (strcmp(key, "test") == 0){
     READ_INT_PAR(env->par.test);
     return(0);
   }
   else if (strcmp(key, "random_seed") == 0){
      READ_INT_PAR(env->par.random_seed);
      tm_par->random_seed = env->par.random_seed;
      return(0);
   }
   else if (strcmp(key, "granularity") == 0){
      READ_DBL_PAR(tm_par->granularity);
      lp_par->granularity = tm_par->granularity;
      return(0);
   }

      /***********************************************************************
       ***                    Master parameters                            ***
       ***********************************************************************/
   else if (strcmp(key, "upper_bound") == 0 ||
	    strcmp(key, "M_upper_bound") == 0){
      READ_DBL_PAR(env->ub);
      env->has_ub = TRUE;
      return(0);
   }
   else if (strcmp(key, "upper_bound_estimate") == 0 ||
	    strcmp(key, "M_upper_bound_estimate") == 0){
      READ_DBL_PAR(env->ub_estimate);
      env->has_ub_estimate = TRUE;
      return(0);
   }
   else if (strcmp(key, "lower_bound") == 0 ||
	    strcmp(key, "M_lower_bound") == 0){
      READ_DBL_PAR(env->lb);
      return(0);
   }
   else if (strcmp(key, "obj_offset") == 0 ||
	    strcmp(key, "M_obj_offset") == 0){
      READ_DBL_PAR(env->obj_offset);
      return(0);
   }
   else if (strcmp(key, "M_verbosity") == 0){
      READ_INT_PAR(env->par.verbosity);
      return(0);
   }
   else if (strcmp(key, "M_random_seed") == 0){
      READ_INT_PAR(env->par.random_seed);
      return(0);
   }
   else if (strcmp(key, "param_file") == 0){
     read_string(env->par.param_file, line, MAX_FILE_NAME_LENGTH);
     return(0);
   }
   else if (strcmp(key, "tm_executable_name") == 0 ||
	    strcmp(key, "tm_exe") == 0 ||
	    strcmp(key, "M_tm_exe") == 0 ||
	    strcmp(key, "M_tm_executable_name") == 0){
      read_string(env->par.tm_exe, line, MAX_FILE_NAME_LENGTH);
      return(0);
   }
   else if (strcmp(key, "dg_executable_name") == 0 ||
	    strcmp(key, "dg_exe") == 0 ||
	    strcmp(key, "M_dg_exe") == 0 ||
	    strcmp(key, "M_dg_executable_name") == 0){
      read_string(env->par.dg_exe, line, MAX_FILE_NAME_LENGTH);
      return(0);
   }
   else if (strcmp(key, "tm_debug") == 0 ||
	    strcmp(key, "M_tm_debug") == 0){
      READ_INT_PAR(env->par.tm_debug);
      if (env->par.tm_debug) env->par.tm_debug = 4;
      return(0);
   }
   else if (strcmp(key, "dg_debug") == 0 ||
	    strcmp(key, "M_dg_debug") == 0){
      READ_INT_PAR(env->par.dg_debug);
      if (env->par.dg_debug) env->par.dg_debug = 4;
      return(0);
   }
   else if (strcmp(key, "tm_machine") == 0 ||
	    strcmp(key, "M_tm_machine") == 0){
	 read_string(env->par.tm_machine, line, MACH_NAME_LENGTH);
	 env->par.tm_machine_set = TRUE;
      return(0);
   }
   else if (strcmp(key, "dg_machine") == 0 ||
	    strcmp(key, "M_dg_machine") == 0){
      read_string(env->par.dg_machine, line, MACH_NAME_LENGTH);
      env->par.dg_machine_set = TRUE;
      return(0);
   }

   else if (strcmp(key, "pvm_trace") == 0 ||
	    strcmp(key, "M_pvm_trace") == 0){
      READ_INT_PAR(env->par.pvm_trace);
      return(0);
   }
   else if (strcmp(key, "do_branch_and_cut") == 0 ||
	    strcmp(key, "M_do_branch_and_cut") == 0){
      READ_INT_PAR(env->par.do_branch_and_cut);
      return(0);
   }
   else if (strcmp(key, "do_draw_graph") == 0 ||
	    strcmp(key, "M_do_draw_graph") == 0){
      READ_INT_PAR(env->par.do_draw_graph);
      return(0);
   }
   else if (strcmp(key, "use_permanent_cut_pools") == 0 ||
	    strcmp(key, "M_use_permanent_cut_pools") == 0){
      READ_INT_PAR(env->par.use_permanent_cut_pools);
      return(0);
   }
   else if (strcmp(key, "mc_compare_solution_tolerance") == 0 ||
	    strcmp(key, "M_mc_compare_solution_tolerance") == 0){
      READ_DBL_PAR(env->par.mc_compare_solution_tolerance);
      return(0);
   }
   else if (strcmp(key, "mc_binary_search_tolerance") == 0 ||
	    strcmp(key, "M_mc_binary_search_tolerance") == 0){
      READ_DBL_PAR(env->par.mc_binary_search_tolerance);
      return(0);
   }
   else if (strcmp(key, "mc_search_order") == 0 ||
	    strcmp(key, "M_mc_search_order") == 0){
      READ_INT_PAR(env->par.mc_search_order);
      return(0);
   }
   else if (strcmp(key, "mc_warm_start") == 0 ||
	     strcmp(key, "M_mc_warm_start") == 0){
      READ_INT_PAR(env->par.mc_warm_start);
      return(0);
   }
   else if (strcmp(key, "mc_warm_start_rule") == 0 ||
	     strcmp(key, "M_mc_warm_start_rule") == 0){
      READ_INT_PAR(env->par.mc_warm_start_rule);
      return(0);
   }
   else if (strcmp(key, "trim_warm_tree") == 0 ||
	     strcmp(key, "M_trim_warm_tree") == 0){
      READ_INT_PAR(env->par.trim_warm_tree);
      return(0);
   }

   /***********************************************************************
    ***                 DrawGraph parameters                            ***
    ***********************************************************************/

   else if (strcmp(key, "source_path") == 0 ||
	    strcmp(key, "DG_source_path") == 0){
      read_string(dg_par->source_path, line, MAX_FILE_NAME_LENGTH);
      return(0);
   }
   else if (strcmp(key, "echo_commands") == 0 ||
	    strcmp(key, "DG_echo_commands") == 0){
      READ_INT_PAR(dg_par->echo_commands);
      return(0);
   }
   else if (strcmp(key, "canvas_width") == 0 ||
	    strcmp(key, "DG_canvas_width") == 0){
      READ_INT_PAR(dg_par->canvas_width);
      return(0);
   }
   else if (strcmp(key, "canvas_height") == 0 ||
	    strcmp(key, "DG_canvas_height") == 0){
      READ_INT_PAR(dg_par->canvas_height);
      return(0);
   }
   else if (strcmp(key, "viewable_width") == 0 ||
	    strcmp(key, "DG_viewable_width") == 0){
      READ_INT_PAR(dg_par->viewable_width);
      return(0);
   }
   else if (strcmp(key, "viewable_height") == 0 ||
	    strcmp(key, "DG_viewable_height") == 0){
      READ_INT_PAR(dg_par->viewable_width);
      return(0);
   }
   else if (strcmp(key, "disp_nodelabels") == 0 ||
	    strcmp(key, "DG_disp_nodelabels") == 0){
      READ_INT_PAR(dg_par->disp_nodelabels);
      return(0);
   }
   else if (strcmp(key, "disp_nodeweights") == 0 ||
	    strcmp(key, "DG_disp_nodeweights") == 0){
      READ_INT_PAR(dg_par->disp_nodeweights);
      return(0);
   }
   else if (strcmp(key, "disp_edgeweights") == 0 ||
	    strcmp(key, "DG_disp_edgeweights") == 0){
      READ_INT_PAR(dg_par->disp_edgeweights);
      return(0);
   }
   else if (strcmp(key, "node_dash") == 0 ||
	    strcmp(key, "DG_node_dash") == 0){
      read_string(dg_par->node_dash, line, MAX_DASH_PATTERN_LENGTH);
      return(0);
   }
   else if (strcmp(key, "edge_dash") == 0 ||
	    strcmp(key, "DG_edge_dash") == 0){
      read_string(dg_par->edge_dash, line, MAX_DASH_PATTERN_LENGTH);
      return(0);
   }
   else if (strcmp(key, "node_radius") == 0 ||
	    strcmp(key, "DG_node_radius") == 0){
      READ_INT_PAR(dg_par->node_radius);
      return(0);
   }
   else if (strcmp(key, "interactive_mode") == 0 ||
	    strcmp(key, "DG_interactive_mode") == 0){
      READ_INT_PAR(dg_par->interactive_mode);
      return(0);
   }
   else if (strcmp(key, "mouse_tracking") == 0 ||
	    strcmp(key, "DG_mouse_tracking") == 0){
      READ_INT_PAR(dg_par->mouse_tracking);
      return(0);
   }
   else if (strcmp(key, "scale_factor") == 0 ||
	    strcmp(key, "DG_scale_factor") == 0){
      READ_DBL_PAR(dg_par->scale_factor);
      return(0);
   }
   else if (strcmp(key, "nodelabel_font") == 0 ||
	    strcmp(key, "DG_nodelabel_font") == 0){
      read_string(dg_par->nodelabel_font, line, MAX_FONT_LENGTH);
      return(0);
   }
   else if (strcmp(key, "nodeweight_font") == 0 ||
	       strcmp(key, "DG_nodeweight_font") == 0){
      read_string(dg_par->nodeweight_font, line, MAX_FONT_LENGTH);
      return(0);
   }
   else if (strcmp(key, "edgeweight_font") == 0 ||
	    strcmp(key, "DG_edgeweight_font") == 0){
      read_string(dg_par->edgeweight_font, line, MAX_FONT_LENGTH);
      return(0);
   }

   /***********************************************************************
    ***                  Treemanager parameters                         ***
    ***********************************************************************/
   if (strcmp(key, "TM_verbosity") == 0){
      READ_INT_PAR(tm_par->verbosity);
      return(0);
   }
   else if (strcmp(key, "TM_granularity") == 0){
      READ_DBL_PAR(tm_par->granularity);
      lp_par->granularity = tm_par->granularity;
      return(0);
   }
   else if (strcmp(key, "lp_executable_name") == 0 ||
	    strcmp(key, "lp_exe") == 0 ||
	    strcmp(key, "TM_lp_exe") == 0 ||
	    strcmp(key, "TM_lp_executable_name") == 0){
      read_string(tm_par->lp_exe, line, MAX_FILE_NAME_LENGTH);
      return(0);
   }
   else if (strcmp(key, "cg_executable_name") == 0 ||
	    strcmp(key, "cg_exe") == 0 ||
	    strcmp(key, "TM_cg_exe") == 0 ||
	    strcmp(key, "TM_cg_executable_name") == 0){
      read_string(tm_par->cg_exe, line, MAX_FILE_NAME_LENGTH);
      return(0);
   }
   else if (strcmp(key, "cp_executable_name") == 0 ||
	    strcmp(key, "cp_exe") == 0 ||
	    strcmp(key, "TM_cp_exe") == 0 ||
	    strcmp(key, "TM_cp_executable_name") == 0){
      read_string(tm_par->cp_exe, line, MAX_FILE_NAME_LENGTH);
      return(0);
   }
   else if (strcmp(key, "lp_debug") == 0 ||
	    strcmp(key, "TM_lp_debug") == 0){
      READ_INT_PAR(tm_par->lp_debug);
      if (tm_par->lp_debug) tm_par->lp_debug = 4;
      return(0);
   }
   else if (strcmp(key, "cg_debug") == 0 ||
	    strcmp(key, "TM_cg_debug") == 0){
      READ_INT_PAR(tm_par->cg_debug);
      if (tm_par->cg_debug) tm_par->cg_debug = 4;
      return(0);
   }
   else if (strcmp(key, "cp_debug") == 0 ||
	    strcmp(key, "TM_cp_debug") == 0){
      READ_INT_PAR(tm_par->cp_debug);
      if (tm_par->cp_debug) tm_par->cp_debug = 4;
      return(0);
   }
   else if (strcmp(key, "max_active_nodes") == 0 ||
	    strcmp(key, "TM_max_active_nodes") == 0){
#if !defined(COMPILE_IN_LP) || defined _OPENMP
      READ_INT_PAR(tm_par->max_active_nodes);
#else
      printf("\nWarning: Trying to use multiple processors with ");
      printf("sequential build...\n");
#endif
      return(0);
   }
   else if (strcmp(key, "max_cp_num") == 0 ||
	    strcmp(key, "TM_max_cp_num") == 0){
      READ_INT_PAR(tm_par->max_cp_num);
      return(0);
   }
   else if (strcmp(key, "lp_mach_num") == 0 ||
	    strcmp(key, "TM_lp_mach_num") == 0){
      READ_INT_PAR(tm_par->lp_mach_num);
      return(0);
   }
   else if (strcmp(key, "cg_mach_num") == 0 ||
	    strcmp(key, "TM_cg_mach_num") == 0){
      READ_INT_PAR(tm_par->cg_mach_num);
      return(0);
   }
   else if (strcmp(key, "cp_mach_num") == 0 ||
	    strcmp(key, "TM_cp_mach_num") == 0){
      READ_INT_PAR(tm_par->cp_mach_num);
      return(0);
   }
#ifndef COMPILE_IN_CG
   else if (strcmp(key, "use_cg") == 0 ||
	    strcmp(key, "TM_use_cg") == 0 ||
	    strcmp(key, "LP_use_cg") == 0){
      READ_INT_PAR(tm_par->use_cg);
      lp_par->use_cg = tm_par->use_cg;
      return(0);
   }
#endif
   else if (strcmp(key, "TM_random_seed") == 0){
      READ_INT_PAR(tm_par->random_seed);
      return(0);
   }
   else if (strcmp(key, "unconditional_dive_frac") == 0 ||
	    strcmp(key, "TM_unconditional_dive_frac") == 0){
      READ_DBL_PAR(tm_par->unconditional_dive_frac);
      return(0);
   }
   else if (strcmp(key, "diving_strategy") == 0 ||
	    strcmp(key, "TM_diving_strategy") == 0){
      READ_INT_PAR(tm_par->diving_strategy);
      return(0);
   }
   else if (strcmp(key, "diving_k") == 0 ||
	    strcmp(key, "TM_diving_k") == 0){
      READ_INT_PAR(tm_par->diving_k);
      return(0);
   }
   else if (strcmp(key, "diving_threshold") == 0 ||
	    strcmp(key, "TM_diving_threshold") == 0){
      READ_DBL_PAR(tm_par->diving_threshold);
      return(0);
   }
   else if (strcmp(key, "node_selection_rule") == 0 ||
	    strcmp(key, "TM_node_selection_rule") == 0){
      READ_INT_PAR(tm_par->node_selection_rule);
      return(0);
   }
   else if (strcmp(key, "keep_description_of_pruned") == 0 ||
	    strcmp(key, "TM_keep_description_of_pruned") == 0){
      READ_INT_PAR(tm_par->keep_description_of_pruned);
      return(0);
   }
   else if(strcmp(key, "keep_warm_start") == 0){
      if (value[0]){
	 tm_par->keep_description_of_pruned =
	   lp_par->keep_description_of_pruned = KEEP_IN_MEMORY;
         lp_par->should_reuse_lp = FALSE; //by asm4
         return(0);
   } else
	 tm_par->keep_description_of_pruned =
	   lp_par->keep_description_of_pruned = DISCARD;
      return(0);
   }
   else if (strcmp(key, "warm_start") == 0 ||
	    strcmp(key, "TM_warm_start") == 0){
      READ_INT_PAR(tm_par->warm_start);
      return(0);
   }
   else if (strcmp(key, "warm_start_node_limit") == 0 ||
	    strcmp(key, "TM_warm_start_node_limit") == 0){
      READ_INT_PAR(tm_par->warm_start_node_limit);
      return(0);
   }
   else if (strcmp(key, "warm_start_node_level") == 0 ||
	    strcmp(key, "TM_warm_start_node_level") == 0){
      READ_INT_PAR(tm_par->warm_start_node_level);
      return(0);
   }
   else if (strcmp(key, "warm_start_node_level_ratio") == 0 ||
	    strcmp(key, "TM_warm_start_node_level_ratio") == 0){
      READ_DBL_PAR(tm_par->warm_start_node_level_ratio);
      return(0);
   }
   else if (strcmp(key, "warm_start_node_ratio") == 0 ||
	    strcmp(key, "TM_warm_start_node_ratio") == 0){
      READ_DBL_PAR(tm_par->warm_start_node_ratio);
      return(0);
   }
   else if (strcmp(key, "vbc_emulation") == 0 ||
	    strcmp(key, "TM_vbc_emulation") == 0){
      READ_INT_PAR(tm_par->vbc_emulation);
      return(0);
   }
   else if (strcmp(key, "logging_interval") == 0 ||
	    strcmp(key, "TM_logging_interval") == 0){
      READ_INT_PAR(tm_par->logging_interval);
      return(0);
   }
   else if (strcmp(key, "logging") == 0 ||
	    strcmp(key, "TM_logging") == 0){
      READ_INT_PAR(tm_par->logging);
      return(0);
   }
   else if (strcmp(key, "price_in_root") == 0 ||
	    strcmp(key, "TM_price_in_root") == 0){
      READ_INT_PAR(tm_par->price_in_root);
      return(0);
   }
   else if (strcmp(key, "trim_search_tree") == 0 ||
	    strcmp(key, "TM_trim_search_tree") == 0){
      READ_INT_PAR(tm_par->trim_search_tree);
      return(0);
   }
   else if (strcmp(key, "colgen_in_first_phase") == 0 ||
	    strcmp(key, "TM_colgen_in_first_phase") == 0){
      READ_INT_PAR(tm_par->colgen_strat[0]);
      return(0);
   }
   else if (strcmp(key, "colgen_in_second_phase") == 0 ||
	    strcmp(key, "TM_colgen_in_second_phase") == 0){
      READ_INT_PAR(tm_par->colgen_strat[1]);
      return(0);
   }
   else if (strcmp(key, "colgen_in_first_phase_str") == 0 ||
	    strcmp(key, "TM_colgen_in_first_phase_str") == 0){
      READ_STRINT_PAR(tm_par->colgen_strat[0],
		      colgen_str, COLGEN_STR_SIZE, value);
      return(0);
   }
   else if (strcmp(key, "colgen_in_second_phase_str") == 0 ||
	    strcmp(key, "TM_colgen_in_second_phase_str") == 0){
      READ_STRINT_PAR(tm_par->colgen_strat[1],
		      colgen_str, COLGEN_STR_SIZE, value);
      return(0);
   }
   else if (strcmp(key, "time_limit") == 0 ||
	    strcmp(key, "TM_time_limit") == 0){
      READ_DBL_PAR(tm_par->time_limit);
      lp_par->time_limit = tm_par->time_limit;
      return(0);
   }
   else if (strcmp(key, "node_limit") == 0 ||
	    strcmp(key, "TM_node_limit") == 0){
      READ_INT_PAR(tm_par->node_limit);
      return(0);
   }
   else if (strcmp(key, "gap_limit") == 0 ||
	    strcmp(key, "TM_gap_limit") == 0){
      READ_DBL_PAR(tm_par->gap_limit);
      return(0);
   }
   else if (strcmp(key, "find_first_feasible") == 0 ||
	    strcmp(key, "TM_find_first_feasible") == 0){
      READ_INT_PAR(tm_par->find_first_feasible);
      return(0);
   }
   else if (strcmp(key, "sensitivity_analysis") == 0 ||
	    strcmp(key, "TM_sensitivity_analysis") == 0 ){
#ifdef SENSITIVITY_ANALYSIS
      READ_INT_PAR(tm_par->sensitivity_analysis);
      if(tm_par->sensitivity_analysis){
	 tm_par->keep_description_of_pruned =
	    lp_par->keep_description_of_pruned = KEEP_IN_MEMORY;
      }else{
	 tm_par->keep_description_of_pruned = DISCARD;
      }
      lp_par->sensitivity_analysis = tm_par->sensitivity_analysis;
#else
      printf("Warning: Sensitivity analysis features are not currently\n");
      printf("         enabled. You must rebuild SYMPHONY with these\n");
      printf("         features enabled in order to use them.\n\n");
#endif
      return(0);
   }

   /***********************************************************************
    ***                      LP parameters                              ***
    ***********************************************************************/
   if (strcmp(key, "LP_verbosity") == 0){
      READ_INT_PAR(lp_par->verbosity);
      return(0);
   }
   else if (strcmp(key, "LP_granularity") == 0){
      READ_DBL_PAR(lp_par->granularity);
      tm_par->granularity = lp_par->granularity;
      return(0);
   }
   else if (strcmp(key, "set_obj_upper_lim") == 0 ||
	    strcmp(key, "LP_set_obj_upper_lim") == 0){
      READ_INT_PAR(lp_par->set_obj_upper_lim);
      return(0);
   }
   else if (strcmp(key, "do_primal_heuristic") == 0 ||
	    strcmp(key, "LP_do_primal_heuristic") == 0){
      READ_INT_PAR(lp_par->do_primal_heuristic);
      return(0);
   }
   else if (strcmp(key, "scaling") == 0 ||
	    strcmp(key, "LP_scaling") == 0){
      READ_INT_PAR(lp_par->scaling);
      return(0);
   }
   else if (strcmp(key, "fastmip") == 0 ||
	    strcmp(key, "LP_fastmip") == 0){
      READ_INT_PAR(lp_par->fastmip);
      return(0);
   }
   else if (strcmp(key, "should_warmstart_chain") == 0 ||
	    strcmp(key, "LP_should_warmstart_chain") == 0){
      READ_INT_PAR(lp_par->should_warmstart_chain);
      return(0);
   }
   else if (strcmp(key, "try_to_recover_from_error") == 0 ||
	    strcmp(key, "LP_try_to_recover_from_error") == 0){
      READ_INT_PAR(lp_par->try_to_recover_from_error);
      return(0);
   }
   else if (strcmp(key, "problem_type") == 0 ||
	    strcmp(key, "LP_problem_type") == 0){
      READ_INT_PAR(lp_par->problem_type);
      return(0);
   }
   else if (strcmp(key, "not_fixed_storage_size") == 0 ||
	    strcmp(key, "LP_not_fixed_storage_size") == 0 ||
	    strcmp(key, "TM_not_fixed_storage_size") == 0 ){
      READ_INT_PAR(lp_par->not_fixed_storage_size);
      tm_par->not_fixed_storage_size = lp_par->not_fixed_storage_size;
      return(0);
   }
   else if (strcmp(key, "cut_pool_check_frequency") == 0 ||
	    strcmp(key, "LP_cut_pool_check_frequency") == 0){
      READ_INT_PAR(lp_par->cut_pool_check_freq);
      return(0);
   }
   else if (strcmp(key, "load_balance_level") == 0 ||
	    strcmp(key, "LP_load_balance_level") == 0){
      READ_INT_PAR(lp_par->load_balance_level);
      return(0);
   }
   else if (strcmp(key, "load_balance_iterations") == 0 ||
	    strcmp(key, "LP_load_balance_iterations") == 0){
      READ_INT_PAR(lp_par->load_balance_iterations);
      return(0);
   }
   else if (strcmp(key, "load_balance_compare_candidates") == 0 ||
	    strcmp(key, "LP_load_balance_compare_candidates") == 0){
      READ_INT_PAR(lp_par->load_balance_compare_candidates);
      return(0);
   }
   else if (strcmp(key, "fractional_diving_ratio") == 0 ||
	    strcmp(key, "LP_fractional_diving_ratio") == 0){
      READ_DBL_PAR(lp_par->fractional_diving_ratio);
      return(0);
   }
   else if (strcmp(key, "fractional_diving_num") == 0 ||
	    strcmp(key, "LP_fractional_diving_num") == 0){
      READ_INT_PAR(lp_par->fractional_diving_num);
      return(0);
   }
   else if (strcmp(key, "max_non_dual_feas_to_add_frac") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_frac") == 0){
      READ_DBL_PAR(lp_par->max_non_dual_feas_to_add_frac);
      return(0);
   }
   else if (strcmp(key, "max_cols_to_add_min") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_min") == 0){
      READ_INT_PAR(lp_par->max_non_dual_feas_to_add_min);
      return(0);
   }
   else if (strcmp(key, "max_non_dual_feas_to_add_max") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_max") == 0){
      READ_INT_PAR(lp_par->max_non_dual_feas_to_add_max);
      return(0);
   }
   else if (strcmp(key, "max_not_fixable_to_add_frac") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_frac") == 0){
      READ_DBL_PAR(lp_par->max_not_fixable_to_add_frac);
      return(0);
   }
   else if (strcmp(key, "max_not_fixable_to_add_min") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_min") == 0){
      READ_INT_PAR(lp_par->max_not_fixable_to_add_min);
      return(0);
   }
   else if (strcmp(key, "max_not_fixable_to_add_max") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_max") == 0){
      READ_INT_PAR(lp_par->max_not_fixable_to_add_max);
      return(0);
   }

   else if (strcmp(key, "mat_col_compress_num") == 0 ||
	    strcmp(key, "LP_mat_col_compress_num") == 0){
      READ_INT_PAR(lp_par->mat_col_compress_num);
      return(0);
   }
   else if (strcmp(key, "mat_col_compress_ratio") == 0 ||
	    strcmp(key, "LP_mat_col_compress_ratio") == 0){
      READ_DBL_PAR(lp_par->mat_col_compress_ratio);
      return(0);
   }
   else if (strcmp(key, "mat_row_compress_num") == 0 ||
	    strcmp(key, "LP_mat_row_compress_num") == 0){
      READ_INT_PAR(lp_par->mat_row_compress_num);
      return(0);
   }
   else if (strcmp(key, "mat_row_compress_ratio") == 0 ||
	    strcmp(key, "LP_mat_row_compress_ratio") == 0){
      READ_DBL_PAR(lp_par->mat_row_compress_ratio);
      return(0);
   }

   else if (strcmp(key, "tailoff_gap_backsteps") == 0 ||
	    strcmp(key, "LP_tailoff_gap_backsteps") == 0){
      READ_INT_PAR(lp_par->tailoff_gap_backsteps);
      return(0);
   }
   else if (strcmp(key, "tailoff_obj_backsteps") == 0 ||
	    strcmp(key, "LP_tailoff_obj_backsteps") == 0){
      READ_INT_PAR(lp_par->tailoff_obj_backsteps);
      return(0);
   }
   else if (strcmp(key, "tailoff_gap_frac") == 0 ||
	    strcmp(key, "LP_tailoff_gap_frac") == 0){
      READ_DBL_PAR(lp_par->tailoff_gap_frac);
      return(0);
   }
   else if (strcmp(key, "tailoff_obj_frac") == 0 ||
	    strcmp(key, "LP_tailoff_obj_frac") == 0){
      READ_DBL_PAR(lp_par->tailoff_obj_frac);
      return(0);
   }
   else if (strcmp(key, "tailoff_absolute") == 0 ||
	    strcmp(key, "LP_tailoff_absolute") == 0){
      READ_DBL_PAR(lp_par->tailoff_absolute);
      return(0);
   }
   else if (strcmp(key, "tailoff_max_no_iterative_impr_iters_root") == 0 ||
	    strcmp(key, "LP_tailoff_max_no_iterative_impr_iters_root") == 0){
      READ_INT_PAR(lp_par->tailoff_max_no_iterative_impr_iters_root);
      return(0);
   }

   else if (strcmp(key, "ineff_cnt_to_delete") == 0 ||
	    strcmp(key, "LP_ineff_cnt_to_delete") == 0){
      READ_INT_PAR(lp_par->ineff_cnt_to_delete);
      return(0);
   }
   else if (strcmp(key, "eff_cnt_before_cutpool") == 0 ||
	    strcmp(key, "LP_eff_cnt_before_cutpool") == 0){
      READ_INT_PAR(lp_par->eff_cnt_before_cutpool);
      return(0);
   }
   else if (strcmp(key, "ineffective_constraints") == 0 ||
	    strcmp(key, "LP_ineffective_constraints") == 0){
      READ_INT_PAR(lp_par->ineffective_constraints);
      return(0);
   }
   else if (strcmp(key, "base_constraints_always_effective") == 0 ||
	    strcmp(key, "LP_base_constraints_always_effective") == 0){
      READ_INT_PAR(lp_par->base_constraints_always_effective);
      return(0);
   }

   else if (strcmp(key, "branch_on_cuts") == 0 ||
	    strcmp(key, "LP_branch_on_cuts") == 0){
      READ_INT_PAR(lp_par->branch_on_cuts);
      return(0);
   }
   else if (strcmp(key, "discard_slack_cuts") == 0 ||
	    strcmp(key, "LP_discard_slack_cuts") == 0){
      READ_INT_PAR(lp_par->discard_slack_cuts);
      return(0);
   }

   /* timeouts on receiving cuts */
   else if (strcmp(key, "first_lp_first_cut_time_out") == 0 ||
	    strcmp(key, "LP_first_lp_first_cut_time_out") == 0){
      READ_DBL_PAR(timeout);
      if (timeout == -1){
	 lp_par->first_lp.first_cut_time_out = 0;
      }else{
	 lp_par->first_lp.first_cut_time_out = timeout;
      }
      return(0);
   }
   else if (strcmp(key, "first_lp_all_cuts_time_out") == 0 ||
	    strcmp(key, "LP_first_lp_all_cuts_time_out") == 0){
      READ_DBL_PAR(timeout);
      if (timeout == -1){
	 lp_par->first_lp.all_cuts_time_out = 0;
      }else{
	 lp_par->first_lp.all_cuts_time_out = timeout;
      }
      return(0);
   }
   else if (strcmp(key, "later_lp_first_cut_time_out") == 0 ||
	    strcmp(key, "LP_later_lp_first_cut_time_out") == 0){
      READ_DBL_PAR(timeout);
      if (timeout == -1){
	 lp_par->later_lp.first_cut_time_out = 0;
      }else{
	 lp_par->later_lp.first_cut_time_out = timeout;
      }
      return(0);
   }
   else if (strcmp(key, "later_lp_all_cuts_time_out") == 0 ||
	    strcmp(key, "LP_later_lp_all_cuts_time_out") == 0){
      READ_DBL_PAR(timeout);
      if (timeout == -1){
	 lp_par->later_lp.all_cuts_time_out = 0;
      }else{
	 lp_par->later_lp.all_cuts_time_out = timeout;
      }
      return(0);
   }

   else if (strcmp(key, "no_cut_timeout") == 0 ||
	    strcmp(key, "LP_no_cut_timeout") == 0){
      lp_par->first_lp.first_cut_time_out = 0;
      lp_par->first_lp.all_cuts_time_out = 0;
      lp_par->later_lp.first_cut_time_out = 0;
      lp_par->later_lp.all_cuts_time_out = 0;
      return(0);
   }
   else if (strcmp(key, "all_cut_timeout") == 0 ||
	    strcmp(key, "LP_all_cut_timeout") == 0){
      READ_DBL_PAR(timeout);
      lp_par->first_lp.first_cut_time_out = timeout;
      lp_par->first_lp.all_cuts_time_out = timeout;
      lp_par->later_lp.first_cut_time_out= timeout;
      lp_par->later_lp.all_cuts_time_out = timeout;
      return(0);
   }

   else if (strcmp(key, "max_cut_num_per_iter") == 0 ||
	    strcmp(key, "LP_max_cut_num_per_iter") == 0){
      READ_INT_PAR(lp_par->max_cut_num_per_iter);
      return(0);
   }

   else if (strcmp(key, "max_cut_num_per_iter_root") == 0 ||
	    strcmp(key, "LP_max_cut_num_per_iter_root") == 0){
      READ_INT_PAR(lp_par->max_cut_num_per_iter_root);
      return(0);
   }

   else if (strcmp(key, "min_root_cut_rounds") == 0 ||
	    strcmp(key, "LP_min_root_cut_rounds") == 0){
      READ_INT_PAR(lp_par->min_root_cut_rounds);
      return(0);
   }

   else if (strcmp(key, "max_cut_length") == 0 ||
	    strcmp(key, "LP_max_cut_length") == 0){
      READ_INT_PAR(lp_par->max_cut_length);
      return(0);
   }

   /* variable fixing parameters */
   else if (strcmp(key, "do_reduced_cost_fixing") == 0 ||
	    strcmp(key, "LP_do_reduced_cost_fixing") == 0){
      READ_INT_PAR(lp_par->do_reduced_cost_fixing);
      return(0);
   }
   else if (strcmp(key, "gap_as_ub_frac") == 0 ||
	    strcmp(key, "LP_gap_as_ub_frac") == 0){
      READ_DBL_PAR(lp_par->gap_as_ub_frac);
      return(0);
   }
   else if (strcmp(key, "gap_as_last_gap_frac") == 0 ||
	    strcmp(key, "LP_gap_as_last_gap_frac") == 0){
      READ_DBL_PAR(lp_par->gap_as_last_gap_frac);
      return(0);
   }
   else if (strcmp(key, "do_logical_fixing") == 0 ||
	    strcmp(key, "LP_do_logical_fixing") == 0){
      READ_INT_PAR(lp_par->do_logical_fixing);
      return(0);
   }
   else if (strcmp(key, "fixed_to_ub_before_logical_fixing") == 0 ||
	    strcmp(key, "LP_fixed_to_ub_before_logical_fixing") == 0){
      READ_INT_PAR(lp_par->fixed_to_ub_before_logical_fixing);
      return(0);
   }
   else if (strcmp(key, "fixed_to_ub_frac_before_logical_fixing")==0 ||
	    strcmp(key, "LP_fixed_to_ub_frac_before_logical_fixing")==0){
      READ_DBL_PAR(lp_par->fixed_to_ub_frac_before_logical_fixing);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_gomory_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_gomory_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_gomory_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_knapsack_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_knapsack_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_knapsack_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_oddhole_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_oddhole_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_oddhole_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_probing_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_probing_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_probing_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_clique_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_clique_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_clique_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_mir_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_mir_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_mir_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_twomir_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_twomir_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_twomir_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_flow_and_cover_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_flow_and_cover_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_flowcover_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_flowcover_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_flowcover_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_flowcover_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_rounding_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_rounding_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_rounding_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_lift_and_project_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_lift_and_project_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_lift_and_project_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_landp_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_landp_cuts") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_landp_cuts);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_gomory_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_gomory_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_gomory_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_knapsack_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_knapsack_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_knapsack_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_oddhole_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_oddhole_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_oddhole_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_probing_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_probing_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_probing_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_clique_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_clique_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_clique_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_mir_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_mir_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_mir_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_twomir_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_twomir_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_twomir_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_flowcover_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_flowcover_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_flowcover_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_rounding_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_rounding_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_rounding_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_lift_and_project_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_lift_and_project_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_lift_and_project_cuts_freq);
      return(0);
   }

   else if (strcmp(key, "generate_cgl_landp_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_landp_cuts_freq") == 0){
      READ_INT_PAR(lp_par->cgl.generate_cgl_landp_cuts_freq);
      return(0);
   }

    else if (strcmp(key, "gomory_max_depth") == 0 ||
	    strcmp(key, "LP_gomory_max_depth") == 0){
      READ_INT_PAR(lp_par->cgl.gomory_max_depth);
      return(0);
   }

   else if (strcmp(key, "knapsack_max_depth") == 0 ||
	    strcmp(key, "LP_knapsack_max_depth") == 0){
      READ_INT_PAR(lp_par->cgl.knapsack_max_depth);
      return(0);
   }

   else if (strcmp(key, "oddhole_max_depth") == 0 ||
	    strcmp(key, "LP_oddhole_max_depth") == 0){
      READ_INT_PAR(lp_par->cgl.oddhole_max_depth);
      return(0);
   }

   else if (strcmp(key, "probing_max_depth") == 0 ||
	    strcmp(key, "LP_probing_max_depth") == 0){
      READ_INT_PAR(lp_par->cgl.probing_max_depth);
      return(0);
   }

   else if (strcmp(key, "clique_max_depth") == 0 ||
	    strcmp(key, "LP_clique_max_depth") == 0){
      READ_INT_PAR(lp_par->cgl.clique_max_depth);
      return(0);
   }

   else if (strcmp(key, "twomir_max_depth") == 0 ||
	    strcmp(key, "LP_twomir_max_depth") == 0){
      READ_INT_PAR(lp_par->cgl.twomir_max_depth);
      return(0);
   }

   else if (strcmp(key, "flowcover_max_depth") == 0 ||
	    strcmp(key, "LP_flowcover_max_depth") == 0){
      READ_INT_PAR(lp_par->cgl.flowcover_max_depth);
      return(0);
   }

   else if (strcmp(key, "max_chain_backtrack") == 0 ||
	    strcmp(key, "LP_max_chain_backtrack") == 0){
      READ_INT_PAR(lp_par->cgl.max_chain_backtrack);
      return(0);
   }

   else if (strcmp(key, "max_chain_trial_num") == 0 ||
	    strcmp(key, "LP_max_chain_trial_num") == 0){
      READ_INT_PAR(lp_par->cgl.max_chain_trial_num);
      return(0);
   }

   else if (strcmp(key, "chain_trial_freq") == 0 ||
	    strcmp(key, "chain_trial_freq") == 0){
      READ_INT_PAR(lp_par->cgl.chain_trial_freq);
      return(0);
   }

   else if (strcmp(key, "chain_weighted_gap") == 0 ||
	    strcmp(key, "LP_chain_weighted_gap") == 0){
      READ_DBL_PAR(lp_par->cgl.chain_weighted_gap);
      return(0);
   }

   else if (strcmp(key, "max_presolve_iter") == 0 ||
	    strcmp(key, "LP_max_presolve_iter") == 0){
      READ_INT_PAR(lp_par->max_presolve_iter);
      lp_par->user_set_max_presolve_iter = TRUE;
      return(0);
   }

   /* user-defined function defaults */
   if (strcmp(key, "is_feasible_default") == 0 ||
	    strcmp(key, "LP_is_feasible_default") == 0){
      READ_INT_PAR(lp_par->is_feasible_default);
      return(0);
   }
   else if (strcmp(key, "send_feasible_solution_default") == 0 ||
	    strcmp(key, "LP_send_feasible_solution_default") == 0){
      READ_INT_PAR(lp_par->send_feasible_solution_default);
      return(0);
   }
   else if (strcmp(key, "display_solution_default") == 0 ||
	    strcmp(key, "LP_display_solution_default") == 0){
      READ_INT_PAR(lp_par->display_solution_default);
      return(0);
   }
   else if (strcmp(key, "shall_we_branch_default") == 0 ||
	    strcmp(key, "LP_shall_we_branch_default") == 0){
      READ_INT_PAR(lp_par->shall_we_branch_default);
      return(0);
   }
   else if (strcmp(key, "select_candidates_default") == 0 ||
	    strcmp(key, "LP_select_candidates_default") == 0){
      READ_INT_PAR(lp_par->select_candidates_default);
      return(0);
   }
   else if (strcmp(key, "strong_branching_cand_num") == 0){
      READ_INT_PAR(lp_par->strong_branching_cand_num_max);
      lp_par->strong_branching_cand_num_min =
	 lp_par->strong_branching_cand_num_max;
      lp_par->user_set_strong_branching_cand_num = TRUE;
      lp_par->strong_branching_red_ratio = 0;
      return(0);
   }
   else if (strcmp(key, "strong_branching_cand_num_max") == 0 ||
	    strcmp(key, "LP_strong_branching_cand_num_max") == 0){
      READ_INT_PAR(lp_par->strong_branching_cand_num_max);
      return(0);
   }
   else if (strcmp(key, "strong_branching_cand_num_min") == 0 ||
	    strcmp(key, "LP_strong_branching_cand_num_min") == 0){
      READ_INT_PAR(lp_par->strong_branching_cand_num_min);
      return(0);
   }
   else if (strcmp(key, "strong_br_min_level") == 0) {
      READ_INT_PAR(lp_par->strong_br_min_level);
      return(0);
   }
   else if (strcmp(key, "strong_br_all_candidates_level") == 0) {
      READ_INT_PAR(lp_par->strong_br_all_candidates_level);
      return(0);
   }
   else if (strcmp(key,"strong_branching_red_ratio") == 0 ||
	    strcmp(key,"LP_strong_branching_red_ratio") == 0){
      READ_DBL_PAR(lp_par->strong_branching_red_ratio);
      return(0);
   }
   else if (strcmp(key,"strong_branching_high_low_weight") == 0 ||
	    strcmp(key,"LP_strong_branching_high_low_weight") == 0){
      READ_DBL_PAR(lp_par->strong_branching_high_low_weight);
      return(0);
   }
   else if (strcmp(key,"use_hot_starts") == 0) {
      READ_INT_PAR(lp_par->use_hot_starts);
      return(0);
   }
   else if (strcmp(key,"should_use_rel_br") == 0) {
      READ_INT_PAR(lp_par->should_use_rel_br);
      return(0);
   }
   else if (strcmp(key,"rel_br_threshold") == 0) {
      READ_INT_PAR(lp_par->rel_br_threshold);
      return(0);
   }
   else if (strcmp(key,"rel_br_cand_threshold") == 0) {
      READ_INT_PAR(lp_par->rel_br_cand_threshold);
      return(0);
   }
   else if (strcmp(key, "compare_candidates_default") == 0 ||
	    strcmp(key, "LP_compare_candidates_default") == 0){
      READ_INT_PAR(lp_par->compare_candidates_default);
      return(0);
   }
   else if (strcmp(key, "compare_candidates_default_str") == 0 ||
	    strcmp(key, "LP_compare_candidates_default_str") == 0){
      READ_STRINT_PAR(lp_par->compare_candidates_default,
		      compare_can_str, COMPARE_CAN_STR_SIZE, value);
      return(0);
   }
   else if (strcmp(key, "select_child_default") == 0 ||
	    strcmp(key, "LP_select_child_default") == 0){
      READ_INT_PAR(lp_par->select_child_default);
      return(0);
   }
   else if (strcmp(key, "pack_lp_solution_default") == 0 ||
	    strcmp(key, "LP_pack_lp_solution_default") == 0){
      READ_INT_PAR(lp_par->pack_lp_solution_default);
      return(0);
   }
   else if (strcmp(key, "multi_criteria") == 0 ||
	    strcmp(key, "LP_multi_criteria") == 0 ){
      READ_INT_PAR(lp_par->multi_criteria);
      env->par.multi_criteria = lp_par->multi_criteria;
      return(0);
   }
   else if (strcmp(key, "mc_find_supported_solutions") == 0 ||
	    strcmp(key, "LP_mc_find_supported_solutions") == 0 ){
      READ_INT_PAR(lp_par->mc_find_supported_solutions);
      return(0);
   }
   else if (strcmp(key, "mc_add_optimality_cuts") == 0 ||
	    strcmp(key, "LP_mc_add_optimality_cuts") == 0 ){
      READ_INT_PAR(lp_par->mc_add_optimality_cuts);
      return(0);
   }
   else if (strcmp(key, "mc_gamma") == 0 ||
	    strcmp(key, "LP_mc_gamma") == 0 ){
      READ_DBL_PAR(lp_par->mc_gamma);
      return(0);
   }
   else if (strcmp(key, "mc_tau") == 0 ||
	    strcmp(key, "LP_mc_tau") == 0 ){
      READ_DBL_PAR(lp_par->mc_tau);
      return(0);
   }
   else if (strcmp(key, "mc_rho") == 0 ||
	    strcmp(key, "LP_mc_rho") == 0 ){
      READ_DBL_PAR(lp_par->mc_rho);
      return(0);
   }
   else if (strcmp(key, "fp_enabled") == 0) {
      READ_INT_PAR(lp_par->fp_enabled);
      return(0);
   }
   else if (strcmp(key, "fp_frequency") == 0) {
      READ_INT_PAR(lp_par->fp_frequency);
      return(0);
   }
   else if (strcmp(key, "fp_max_cycles") == 0) {
      READ_INT_PAR(lp_par->fp_max_cycles);
      return(0);
   }
   else if (strcmp(key, "fp_time_limit") == 0) {
      READ_DBL_PAR(lp_par->fp_time_limit);
      return(0);
   }
   else if (strcmp(key, "fp_display_interval") == 0) {
      READ_DBL_PAR(lp_par->fp_display_interval);
      return(0);
   }
   else if (strcmp(key, "fp_flip_fraction") == 0) {
      READ_DBL_PAR(lp_par->fp_flip_fraction);
      return(0);
   }
   else if (strcmp(key, "fp_max_initial_time") == 0) {
      READ_DBL_PAR(lp_par->fp_max_initial_time);
      return(0);
   }
   else if (strcmp(key, "fp_min_gap") == 0) {
      READ_DBL_PAR(lp_par->fp_min_gap);
      return(0);
   }

   /***********************************************************************
    ***                     cut_gen parameters                          ***
    ***********************************************************************/
   else if (strcmp(key, "CG_verbosity") == 0){
      READ_INT_PAR(cg_par->verbosity);
      return(0);
   }
   else if (strcmp(key, "do_findcuts") == 0 ||
	    strcmp(key, "CG_do_findcuts") == 0){
      READ_INT_PAR(cg_par->do_findcuts);
      return(0);
   }

   /***********************************************************************
    ***                      cutpool parameters                         ***
    ***********************************************************************/
   if (strcmp(key, "CP_verbosity") == 0){
      READ_INT_PAR(cp_par->verbosity);
      return(0);
   }
   else if (strcmp(key, "cp_warm_start") == 0 ||
	    strcmp(key, "CP_warm_start") == 0){
      READ_INT_PAR(cp_par->warm_start);
      return(0);
   }
   else if (strcmp(key, "cp_logging") == 0 ||
	    strcmp(key, "CP_logging") == 0){
      READ_INT_PAR(cp_par->logging);
      return(0);
   }
   else if (strcmp(key, "block_size") == 0 ||
	    strcmp(key, "CP_block_size") == 0){
      READ_INT_PAR(cp_par->block_size);
      return(0);
   }
   else if (strcmp(key, "max_size") == 0 ||
	    strcmp(key, "CP_max_size") == 0){
      READ_INT_PAR(cp_par->max_size);
      return(0);
   }
   else if (strcmp(key, "max_number_of_cuts") == 0 ||
	    strcmp(key, "CP_max_number_of_cuts") == 0){
      READ_INT_PAR(cp_par->max_number_of_cuts);
      return(0);
   }
   else if (strcmp(key, "cuts_to_check") == 0 ||
	    strcmp(key, "cuts_to_check") == 0){
      READ_INT_PAR(cp_par->cuts_to_check);
      return(0);
   }
   else if (strcmp(key, "delete_which") == 0 ||
	    strcmp(key, "CP_delete_which") == 0){
      READ_INT_PAR(cp_par->delete_which);
      return(0);
   }
   else if (strcmp(key, "touches_until_deletion") == 0 ||
	    strcmp(key, "CP_touches_until_deletion") == 0){
      READ_INT_PAR(cp_par->touches_until_deletion);
      return(0);
   }
   else if (strcmp(key, "min_to_delete") == 0 ||
	    strcmp(key, "CP_min_to_delete") == 0){
      READ_INT_PAR(cp_par->min_to_delete);
      return(0);
   }
   else if (strcmp(key, "check_which") == 0 ||
	       strcmp(key, "CP_check_which") == 0){
      READ_INT_PAR(cp_par->check_which);
      return(0);
   }

   /*************************************************************************
    ***                     preprocessing - parameters                    ***
    *************************************************************************/

   if (strcmp(key, "prep_level") == 0){
      READ_INT_PAR(prep_par->level);
      return(0);
   }
   else if (strcmp(key, "prep_dive_level") == 0){
      READ_INT_PAR(prep_par->dive_level);
      return(0);
   }
   else if (strcmp(key, "prep_impl_dive_level") == 0){
      READ_INT_PAR(prep_par->impl_dive_level);
      return(0);
   }
   else if (strcmp(key, "prep_impl_limit") == 0){
      READ_INT_PAR(prep_par->impl_limit);
      return(0);
   }
   else if (strcmp(key, "prep_iter_limit") == 0){
      READ_INT_PAR(prep_par->iteration_limit);
      return(0);
   }
   else if (strcmp(key, "prep_do_probing") == 0){
      READ_INT_PAR(prep_par->do_probe);
      return(0);
   }
   else if (strcmp(key, "prep_do_sr") == 0){
      READ_INT_PAR(prep_par->do_single_row_rlx);
      return(0);
   }
   else if (strcmp(key, "prep_verbosity") == 0){
      READ_INT_PAR(prep_par->verbosity);
      return(0);
   }
   else if (strcmp(key, "prep_reduce_mip") == 0){
      READ_INT_PAR(prep_par->reduce_mip);
      return(0);
   }
   else if (strcmp(key, "prep_probing_verbosity") == 0){
      READ_INT_PAR(prep_par->probe_verbosity);
      return(0);
   }
   else if (strcmp(key, "prep_probing_level") == 0){
      READ_INT_PAR(prep_par->probe_level);
      return(0);
   }
   else if (strcmp(key, "prep_display_stats") == 0){
      READ_INT_PAR(prep_par->display_stats);
      return(0);
   }
   else if (strcmp(key, "max_sr_cnt") == 0){
      READ_INT_PAR(prep_par->max_sr_cnt);
      return(0);
   }
   else if (strcmp(key, "max_aggr_row_cnt") == 0){
      READ_INT_PAR(prep_par->max_aggr_row_cnt);
      return(0);
   }
   else if (strcmp(key, "keep_row_ordered") == 0){
      READ_INT_PAR(prep_par->keep_row_ordered);
      return(0);
   }
   else if (strcmp(key, "write_mps") == 0){
      READ_INT_PAR(prep_par->write_mps);
      return(0);
   }
   else if (strcmp(key, "write_lp") == 0){
      READ_INT_PAR(prep_par->write_lp);
      return(0);
   }
   else if (strcmp(key, "prep_time_limit") == 0){
      READ_INT_PAR(prep_par->time_limit);
      return(0);
   }

   return(FUNCTION_TERMINATED_ABNORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc *create_copy_warm_start(warm_start_desc *ws)
{
   int i, num=0, allocated_cut_num = 0;
   warm_start_desc * ws_copy;

   if (!ws){
      printf("create_copy_warm_start():");
      printf("The warm start description is empty!\n");
      return(NULL);
   }

   ws_copy = (warm_start_desc*)calloc(1,sizeof(warm_start_desc));
   memcpy(ws_copy, ws, sizeof(warm_start_desc));
   num = ws_copy->cut_num;
   allocated_cut_num = ws_copy->allocated_cut_num;
   ws_copy->cuts = (cut_data**)
      calloc(allocated_cut_num, sizeof(cut_data*));
   for(i = 0; i<num; i++){
      ws_copy->cuts[i] = (cut_data*)calloc(1,sizeof(cut_data));
      memcpy(ws_copy->cuts[i], ws->cuts[i], sizeof(cut_data));
      ws_copy->cuts[i]->coef =
	 (char*)calloc(ws_copy->cuts[i]->size,CSIZE);
      memcpy(ws_copy->cuts[i]->coef, ws->cuts[i]->coef,
	     CSIZE* ws_copy->cuts[i]->size);
   }
   ws_copy->rootnode = (bc_node*)calloc(1,sizeof(bc_node));
   copy_tree(ws_copy->rootnode, ws->rootnode);


   if(ws->best_sol.xlength){
      ws_copy->best_sol.xind = (int*) malloc (ISIZE * ws->best_sol.xlength);
      ws_copy->best_sol.xval = (double*) malloc (DSIZE * ws->best_sol.xlength);

      memcpy(ws_copy->best_sol.xind, ws->best_sol.xind,
	     ISIZE * ws->best_sol.xlength);
      memcpy(ws_copy->best_sol.xval, ws->best_sol.xval,
	     DSIZE * ws->best_sol.xlength);
   }

   return(ws_copy);
}

/*===========================================================================*/
/*===========================================================================*/

MIPdesc *create_copy_mip_desc(MIPdesc * mip)
{
   MIPdesc * mip_copy;
   int i;

   if (mip){
      mip_copy = (MIPdesc*) calloc(1, sizeof(MIPdesc));
      memcpy(mip_copy, mip, sizeof(MIPdesc));

      if (mip->n){
	 mip_copy->obj       = (double *) malloc(DSIZE * mip_copy->n);

	 mip_copy->ub        = (double *) malloc(DSIZE * mip_copy->n);
	 mip_copy->lb        = (double *) malloc(DSIZE * mip_copy->n);
	 mip_copy->is_int    = (char *)   malloc(CSIZE * mip_copy->n);
	 mip_copy->matbeg    = (int *)    malloc(ISIZE * (mip_copy->n + 1));

	 memcpy(mip_copy->obj,    mip->obj,    DSIZE * mip_copy->n);
	 memcpy(mip_copy->ub,     mip->ub,     DSIZE * mip_copy->n);
	 memcpy(mip_copy->lb,     mip->lb,     DSIZE * mip_copy->n);
	 memcpy(mip_copy->is_int, mip->is_int, CSIZE * mip_copy->n);
	 memcpy(mip_copy->matbeg, mip->matbeg, ISIZE * (mip_copy->n + 1));

	 if(mip->obj1){
	    mip_copy->obj1      = (double *) malloc(DSIZE * mip_copy->n);
	    memcpy(mip_copy->obj1,   mip->obj1,   DSIZE * mip_copy->n);
	 }
	 if(mip->obj2){
	    mip_copy->obj2      = (double *) malloc(DSIZE * mip_copy->n);
	    memcpy(mip_copy->obj2,   mip->obj2,   DSIZE * mip_copy->n);
	 }
      }

      if (mip->m){
	 mip_copy->rhs    = (double *) malloc(DSIZE * mip_copy->m);
	 mip_copy->sense  = (char *)   malloc(CSIZE * mip_copy->m);
	 mip_copy->rngval = (double *) malloc(DSIZE * mip_copy->m);

	 memcpy(mip_copy->rhs, mip->rhs,       DSIZE * mip_copy->m);
	 memcpy(mip_copy->sense, mip->sense,   CSIZE * mip_copy->m);
	 memcpy(mip_copy->rngval, mip->rngval, DSIZE * mip_copy->m);
      }

      if (mip->nz){
	 mip_copy->matval     = (double *) malloc(DSIZE*mip_copy->nz);
	 mip_copy->matind     = (int *)    malloc(ISIZE*mip_copy->nz);
	 memcpy(mip_copy->matval, mip->matval, DSIZE * mip_copy->nz);
	 memcpy(mip_copy->matind, mip->matind, ISIZE * mip_copy->nz);
      }

      /* will not be used anywhere other than presolve where they are
	 always initialized */
      mip_copy->mip_inf = 0;
      mip->cru_vars = 0;
      mip->orig_sense = 0;
      mip->orig_ind = 0;
#if 0
      if (mip->row_matbeg){
	 mip_copy->row_matbeg  = (int *) malloc(ISIZE * (mip_copy->m + 1));
	 mip_copy->row_matind = (int *)    malloc(ISIZE*mip_copy->nz);
	 mip_copy->row_matval = (double *) malloc(DSIZE*mip_copy->nz);
	 mip_copy->row_lengths = (int *) malloc(ISIZE * mip_copy->m);
	 mip_copy->col_lengths    = (int *)    malloc(ISIZE * mip_copy->n);

	 memcpy(mip_copy->row_matbeg, mip->row_matbeg,
		ISIZE*(mip_copy->m + 1));
	 memcpy(mip_copy->row_matind, mip->row_matind, ISIZE*mip_copy->nz);
	 memcpy(mip_copy->row_matval, mip->row_matval, DSIZE*mip_copy->nz);
	 memcpy(mip_copy->row_lengths, mip->row_lengths, ISIZE*mip_copy->m);
	 memcpy(mip_copy->col_lengths, mip->col_lengths, ISIZE * mip_copy->n);
      }
#endif

      if (mip->colname){
	 mip_copy->colname = (char**)calloc(sizeof(char*), mip_copy->n);

	 for(i=0; i<mip_copy->n; i++){
	    /* FIXME! Resctricting col_name to 20 chars! */
	    if(mip->colname[i]){
	       mip_copy->colname[i] = (char*)malloc(CSIZE*20);
	       strncpy(mip_copy->colname[i], mip->colname[i], 20);
	       mip_copy->colname[i][19] = 0;
	    }
	 }
      }

      if(mip->fixed_n){
	 memcpy(mip_copy->fixed_ind, mip->fixed_ind, ISIZE*mip->fixed_n);
	 memcpy(mip_copy->fixed_val, mip->fixed_val, DSIZE*mip->fixed_n);
      }
   }
   else{
      printf("create_copy_mip_desc():");
      printf("Trying to copy an empty mip desc!\n");
      return(NULL);
   }

   return(mip_copy);
}

/*===========================================================================*/
/*===========================================================================*/

sym_environment *create_copy_environment (sym_environment *env)
{
   int i, j, num;
   sym_environment *env_copy;
   params *par;
   lp_sol *sol;
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_CP)
   cp_cut_data *cp_cut;
   cut_data *cut;
#endif

   if (!env){
      printf("create_copy_sym_environment(): The given problem is empty!\n");
      printf("Unable to copy.\n");
      return(NULL);
   }
   env_copy = (sym_environment*) calloc(1, sizeof(sym_environment));
   memcpy(env_copy, env, sizeof(sym_environment));

   //   initialize_u(env_copy);

   /* Note that, if some modifications on the user function have been done
      after initialization, it will not be reflected here, since SYMPHONY
      doesn't know anoything about the user structure! For a temporary
      solution, the user pointer will be directed to the original user
      structure! So, be careful from now on that, further modifications
      on the user structure of either the original or the clone env will
      affect the both! */

   env_copy->user = env->user;

   /*========================================================================*/
   /*   copy params */

   par = &(env_copy->par);

   if(par->tm_par.lp_mach_num)
      par->tm_par.lp_machs =
	 (char**)malloc(sizeof(char*)*par->tm_par.lp_mach_num);
   if(par->tm_par.cg_mach_num)
      par->tm_par.cg_machs =
	 (char**)malloc(sizeof(char*)*par->tm_par.cg_mach_num);
   if(par->tm_par.cp_mach_num)
      par->tm_par.cp_machs =
	 (char**)malloc(sizeof(char*)*par->tm_par.cp_mach_num);
   for(i = 0; i<par->tm_par.lp_mach_num; i++){
      par->tm_par.lp_machs[i] =
	 (char*)malloc(CSIZE*(MACH_NAME_LENGTH+1));
      memcpy(par->tm_par.lp_machs[i], env->par.tm_par.lp_machs[i],
	     CSIZE*(MACH_NAME_LENGTH+1));
   }

   for(i = 0; i<par->tm_par.cg_mach_num; i++){
      par->tm_par.cg_machs[i] =
	 (char*)malloc(CSIZE*(MACH_NAME_LENGTH+1));
      memcpy(par->tm_par.cg_machs[i], env->par.tm_par.cg_machs[i],
	     CSIZE*(MACH_NAME_LENGTH+1));
   }

   for(i = 0; i<par->tm_par.cp_mach_num; i++){
      par->tm_par.cp_machs[i] =
	 (char*)malloc(CSIZE*(MACH_NAME_LENGTH+1));
      memcpy(par->tm_par.cp_machs[i], env->par.tm_par.cp_machs[i],
	     CSIZE*(MACH_NAME_LENGTH+1));
   }


   /*========================================================================*/

   /* copy lp_sol */

   sol = &(env_copy->best_sol);
   if(sol->xlength){
      sol->xind = (int *)malloc(ISIZE * sol->xlength);
      sol->xval = (double *)malloc(DSIZE * sol->xlength);
      memcpy(sol->xind, env->best_sol.xind, ISIZE*sol->xlength);
      memcpy(sol->xval, env->best_sol.xval, DSIZE*sol->xlength);
   }

   /*========================================================================*/

   /* copy mip */
   if (env->mip){
      //free_mip_desc(env_copy->mip);
      if(env->prep_mip){
	 env_copy->prep_mip = create_copy_mip_desc(env->prep_mip);
	 env_copy->orig_mip = create_copy_mip_desc(env->orig_mip);
	 env_copy->mip = env_copy->orig_mip;
      }else{
	 env_copy->mip = create_copy_mip_desc(env->mip);
	 env_copy->prep_mip = env_copy->orig_mip = 0;
      }
   }

   /*========================================================================*/

   /* copy base_desc */

   if (env->base){
      env_copy->base = (base_desc*) calloc(1, sizeof(base_desc));
      memcpy(env_copy->base, env->base, sizeof(base_desc));

      if (env->base->varnum){
	 env_copy->base->userind = (int *) malloc(ISIZE*env->base->varnum);
	 memcpy(env_copy->base->userind, env->base->userind,
		ISIZE*env->base->varnum);
      }
   }

   /*========================================================================*/

   /* copy root_desc */

   if (env->rootdesc){
      node_desc *desc = env_copy->rootdesc =
	 (node_desc *) calloc(1, sizeof(node_desc));
      memcpy(desc, env->rootdesc, sizeof(node_desc));

      if (desc->uind.size){
	 desc->uind.list = (int *) malloc(desc->uind.size*ISIZE);
	 memcpy( desc->uind.list,  env->rootdesc->uind.list,
		 desc->uind.size*ISIZE);
      }

      if (desc->not_fixed.size){
	 desc->not_fixed.list =
	    (int *) malloc(desc->not_fixed.size*ISIZE);
	 memcpy( desc->not_fixed.list,  env->rootdesc->not_fixed.list,
		 desc->not_fixed.size*ISIZE);
      }

      if (desc->cutind.size){
	 desc->cutind.list = (int *) malloc(desc->cutind.size*ISIZE);
	 memcpy( desc->cutind.list,  env->rootdesc->cutind.list,
		 desc->cutind.size*ISIZE);
      }

      if (desc->desc_size){
	 desc->desc = (char*) malloc(desc->desc_size*CSIZE);
	 memcpy(desc->desc, env->rootdesc->desc,
		desc->desc_size*CSIZE);
      }

      if(desc->bnd_change){
	 desc->bnd_change = (bounds_change_desc *)
	    calloc(sizeof(bounds_change_desc),1);
	 if(env->rootdesc->bnd_change->num_changes){
	    desc->bnd_change->index =(int *)malloc
	       (ISIZE *env->rootdesc->bnd_change->num_changes);
	    desc->bnd_change->lbub =(char *)malloc
	       (CSIZE *env->rootdesc->bnd_change->num_changes);
	    desc->bnd_change->value =(double *)malloc
	       (DSIZE *env->rootdesc->bnd_change->num_changes);
	    memcpy(desc->bnd_change->index,
		   env->rootdesc->bnd_change->index,
		   ISIZE *env->rootdesc->bnd_change->num_changes);
	    memcpy(desc->bnd_change->lbub,
		   env->rootdesc->bnd_change->lbub,
		   CSIZE *env->rootdesc->bnd_change->num_changes);
	    memcpy(desc->bnd_change->value,
		   env->rootdesc->bnd_change->value,
		   DSIZE *env->rootdesc->bnd_change->num_changes);
	 }
	 desc->bnd_change->num_changes =
	    env->rootdesc->bnd_change->num_changes;
      }
   }

   /*========================================================================*/
   /* jump the tm */


   /*========================================================================*/
   /* copy the warm start */

   if (env->warm_start){
      env_copy->warm_start = create_copy_warm_start(env->warm_start);
   }

   /*========================================================================*/
   /*copy the cut pool */

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_CP)
   if (env_copy->par.tm_par.max_cp_num > 1){
      env_copy->cp =
	 (cut_pool **) malloc(env_copy->par.tm_par.max_cp_num*
			      sizeof(cut_pool *));
      for (i = 0; i < env_copy->par.tm_par.max_cp_num; i++){
	 env_copy->cp[i] = (cut_pool *) calloc(1, sizeof(cut_pool));
	 env_copy->cp[i]->par = env_copy->par.cp_par;
#ifdef USE_SYM_APPLICATION
	 user_send_cp_data(env_copy->user, &env_copy->cp[i]->user);
#else
	 env_copy->cp[i]->user = env_copy->user;
#endif
      }
      num = env_copy->par.tm_par.max_cp_num;
   }else{
      num = 0;
   }

   if (num){
      for (i = 0; i < num; i++){
	 memcpy(env_copy->cp[i], env->cp[i], sizeof(cut_pool));
	 env_copy->cp[i]->cuts =
	    (cp_cut_data**)malloc(env_copy->cp[i]->allocated_cut_num*
				  sizeof(cp_cut_data*));
	    for(j = 0; j<env_copy->cp[i]->cut_num; j++){
	       env_copy->cp[i]->cuts[j] =
	       (cp_cut_data*)calloc(1, sizeof(cp_cut_data));
	       cp_cut = env_copy->cp[i]->cuts[j];
	       memcpy(cp_cut, env->cp[i]->cuts[j], sizeof(cp_cut_data));
	       cp_cut->cut.coef = (char*)malloc(cp_cut->cut.size*CSIZE);
	       memcpy(cp_cut->cut.coef, env->cp[i]->cuts[j]->cut.coef,
		      cp_cut->cut.size*CSIZE);
	    }

	 sol = &(env_copy->cp[i]->cur_sol);

	 sol->xind = (int *)malloc(ISIZE * sol->max_sol_length);
	 sol->xval = (double *)malloc(DSIZE * sol->max_sol_length);
	 memcpy(sol->xind, env->cp[i]->cur_sol.xind, ISIZE*sol->max_sol_length);
	 memcpy(sol->xval, env->cp[i]->cur_sol.xval, DSIZE*sol->max_sol_length);

	 num = env_copy->cp[i]->cuts_to_add_num;
	 if (num){
	    env_copy->cp[i]->cuts_to_add =
	       (cut_data**)malloc(num * sizeof(cut_data*));
	    for(j = 0; j<num; j++){
	       env_copy->cp[i]->cuts_to_add[j] =
		  (cut_data*)calloc(1, sizeof(cut_data));
	       cut = env_copy->cp[i]->cuts_to_add[j];
	       memcpy(cut,env_copy->cp[i]->cuts_to_add[j], sizeof(cut_data));
	       cut->coef = (char*)malloc(cut->size*CSIZE);
	       memcpy(cut->coef, env->cp[i]->cuts_to_add[j]->coef,
		      cut->size*CSIZE);
	    }
	 }
      }
   }
#endif

   return(env_copy);
}

/*===========================================================================*/
/*===========================================================================*/
double get_lb_for_new_rhs(bc_node *root, MIPdesc *mip, int cnt, int *ind,
			  double *val)
{
#ifdef SENSITIVITY_ANALYSIS
   int i, j, retval;
   double min = SYM_INFINITY;
   bc_node * child;

   if(root){

      root->C_LP = root->lower_bound;
      for(i=0; i<cnt; i++){
	 root->C_LP += root->duals[ind[i]]*(val[i] - mip->rhs[ind[i]]);
      }

      for(i = 0; i < root->bobj.child_num; i++){

	 child = root->children[i];
	 child->C_LP = 0.0;

	 if(child->node_status == NODE_STATUS__PRUNED){
	    if(child->feasibility_status == FEASIBLE_PRUNED ||
	       child->feasibility_status == OVER_UB_PRUNED){

	       if (child->feasibility_status == FEASIBLE_PRUNED){
		  for(j = 0; j< child->sol_size; j++){
		     child->C_LP += mip->obj[child->sol_ind[j]] * child->sol[j];
		  }
	       } else{
		  printf("OVER_UB_PRUNED!\n");
		  child->C_LP = child->lower_bound;
	       }

	       for(j=0; j<cnt; j++){
		  child->C_LP += child->duals[ind[j]] *
		     (val[j]- mip->rhs[ind[j]]);
	       }
	       child->B_IP = child->C_LP;
	    }
	    else if (child->feasibility_status == INFEASIBLE_PRUNED){

	       child->B_IP = SYM_INFINITY;
	       retval = check_feasibility_new_rhs(child, mip, cnt, ind, val);
	       /* if(retval == LP_D_UNBOUNDED || retval == LP_ABANDONED ||
		   retval == LP_D_INFEASIBLE){
		   child->B_IP = SYM_INFINITY;
	       } */
	       if(retval == LP_OPTIMAL || retval == LP_D_OBJLIM ||
		  retval == LP_D_ITLIM){
		  child->B_IP = -SYM_INFINITY;
	       }
	    }
	    else {
	       printf("get_lb_for_new_rhs(): Unknown error!\n");
	       exit(1);
	    }
	 }
	 else {
	    child->B_IP = get_lb_for_new_rhs(child, mip, cnt, ind, val);
	 }

	 if(child->B_IP < min)
	    min = child->B_IP;
      }

      if(root->C_LP > min )
	 return (root->C_LP);
      else
	 return (min);
   }

   return (min);
#else
   printf("get_lb_for_new_rhs():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(-SYM_INFINITY);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

double get_ub_for_new_obj(bc_node *root, MIPdesc *mip, int cnt,
				int *ind, double *val)
{
#ifdef SENSITIVITY_ANALYSIS
   int i, j, n;
   double inf = SYM_INFINITY;
   double min = inf;
   bc_node * child;
   double valuesi = 0.0, lpetol =  9.9999999999999995e-07;
   double objval = inf;
   double objval_pr = inf;
   double * sol = NULL;
   if(root){

      for(n = 0; n < root->bobj.child_num; n++){

	 child = root->children[n];

	 if(child->node_status == NODE_STATUS__PRUNED ||
	    child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){

	    if (child->feasibility_status == FEASIBLE_PRUNED ||
		child->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
		child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){

	       /* too much memory to keep the OVER_UB_PRUNED solutions -
		  if sensitivity-analysis is set,
		  try the feasibility test and if succeeded then keep them*/

	       if (sol)
		  FREE(sol);
	       sol = (double*)calloc(mip->n, DSIZE);

	       objval = 0.0;
	       for(j = 0; j< child->sol_size; j++){
		  objval += mip->obj[child->sol_ind[j]] * child->sol[j];
		  sol[child->sol_ind[j]] = child->sol[j];
	       }

	       for(j=0; j<cnt; j++){
		  objval += sol[ind[j]] *
		     (val[j] - mip->obj[ind[j]]);
	       }

	       if(child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){
		  objval_pr = get_ub_for_new_obj(child, mip, cnt, ind, val);
		  objval = objval < objval_pr ? objval : objval_pr;
	       }
	    }

	 }else{
	    objval = get_ub_for_new_obj(child, mip, cnt, ind, val);
	 }

	 if(objval < min){
	    min = objval;
	 }
      }
   }

   if (sol)
      FREE(sol);

   return (min);
#else
   printf("get_ub_for_new_obj():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(SYM_INFINITY);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

/* send in a row oriented mip description! */
double get_ub_for_new_rhs(bc_node *root, MIPdesc *mip, int cnt,
			  int *ind, double *val)
{
#ifdef SENSITIVITY_ANALYSIS
   int i, j, k, n;
   double inf = SYM_INFINITY;
   double min = inf;
   bc_node * child;
   double valuesi = 0.0, lpetol =  9.9999999999999995e-07;
   double objval = inf, row_val = 0.0;
   double objval_pr;
   int feasible;
   int nonzeros;
   double * sol = NULL;

   int *matbeg = mip->matbeg, *matind = mip->matind;
   double *matval = mip->matval;

   if(root){

      for(n = 0; n < root->bobj.child_num; n++){

	 child = root->children[n];

	 if(child->node_status == NODE_STATUS__PRUNED ||
	    child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){

	    if (child->feasibility_status == FEASIBLE_PRUNED ||
		child->feasibility_status == PRUNED_HAS_CAN_SOLUTION ||
		child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){

	       if (sol)
		  FREE(sol);
	       sol = (double*)calloc(mip->n, DSIZE);
	       objval = 0.0;

	       for(j = 0; j< child->sol_size; j++){
		  objval += mip->obj[child->sol_ind[j]] * child->sol[j];
		  sol[child->sol_ind[j]] = child->sol[j];
	       }

	       /* see whether it is feasible for the new rhs! */
	       feasible = TRUE;

	       for(i=0; i<cnt; i++){
		  row_val = 0.0;
		  for(j=matbeg[ind[i]]; j<matbeg[ind[i]+1]; j++){
		     row_val += matval[j] * sol[matind[j]];
		  }
		  switch(mip->sense[ind[i]]){
		   case 'L':
		      if(row_val > val[i] + lpetol){
			 feasible = FALSE;
		      }
		      break;
		   case 'G':
		      if(row_val < val[i] - lpetol){
			 feasible = FALSE;
		      }
		      break;
		   case 'E':
		      if(!((row_val > val[i] - lpetol) &&
			   (row_val < val[i] + lpetol))){
			 feasible = FALSE;
		      }
		      break;
		   case 'R':
		      if(row_val > val[i] + lpetol
			 || row_val < val[i] - mip->rngval[i] - lpetol){
			 feasible = FALSE;
		      }
		      break;
		   case 'N':
		      break;
		  }
		  if(!feasible){
		     break;
		  }
	       }

	       if(!feasible){
		  objval = inf;
	       }

	       if(child->feasibility_status == NOT_PRUNED_HAS_CAN_SOLUTION){
		  objval_pr = get_ub_for_new_rhs(child, mip, cnt, ind, val);
		  objval = objval < objval_pr ? objval : objval_pr;
	       }
	    }
	 } else {
	    objval = get_ub_for_new_rhs(child, mip, cnt, ind, val);
	 }

	 if(objval < min){
	    min = objval;
	 }
      }
   }

   return (min);
#else
   printf("get_ub_for_new_rhs():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(SYM_INFINITY);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int check_feasibility_new_rhs(bc_node * node, MIPdesc * mip,
				 int cnt, int *ind, double *val)
{
#ifdef SENSITIVITY_ANALYSIS
   int i, j;
   int level = node->bc_level;
   bc_node **path, *n;
   branch_desc *bpath ;
   branch_obj * bobj;
   int retval, iterd;
   LPdata * lp_data;

   double * old_obj = mip->obj;
   //   mip->obj = (double*)calloc(mip->n, DSIZE);

   lp_data = (LPdata *) malloc (sizeof(LPdata));
   lp_data->mip = mip;
   lp_data->n = mip->n;
   lp_data->m = mip->m;

   path = (bc_node **) malloc((2*(level+1)+BB_BUNCH)*sizeof(bc_node *));
   bpath = (branch_desc *) malloc
      ((2*(level+1)+BB_BUNCH)*sizeof(branch_desc));

   for (i = level, n = node; i >= 0; n = n->parent, i--)
      path[i] = n;

   for (i = 0; i < level; i++, bpath++){
      for (j = path[i]->bobj.child_num - 1; j >= 0; j--)
	 if (path[i]->children[j] == path[i+1])
	    break;
      bobj = &path[i]->bobj;
      bpath->type = bobj->type;
      bpath->name = bobj->name;
      bpath->sense = bobj->sense[j];
      bpath->rhs = bobj->rhs[j];      bpath->range = bobj->range[j];
      bpath->branch = bobj->branch[j];
   }

   //load_the problem

   open_lp_solver(lp_data);
   load_lp_prob(lp_data, 0, 0);

   //change the rhs!

   //FIXME! cange_rhs needs lp_data->tmp.c and lp_data->tmp.d???

   lp_data->tmp.c = (char*) calloc(mip->m, CSIZE);
   lp_data->tmp.d = (double*) calloc(mip->m, DSIZE);

   change_rhs(lp_data, cnt, ind, val);

   //add the branching changes

   bpath = bpath - level;
   if(level){
      for (i = 0; i < level; i++, bpath++){
	 //	 bpath = bpath + i;
	 if (bpath->type == BRANCHING_VARIABLE){
	    j = bpath->name;   //assuming no extra vars!
	    switch (bpath->sense){
	    case 'E':
	       change_lbub(lp_data, j, bpath->rhs, bpath->rhs);
	       break;
	     case 'L':
	       change_ub(lp_data, j, bpath->rhs);
	       break;
	     case 'G':
	       change_lb(lp_data, j, bpath->rhs);
	       break;
	     case 'R':
	       change_lbub(lp_data, j, bpath->rhs, bpath->rhs + bpath->range);
	       break;
	    }
	 }else{ /* BRANCHING_CUT */
	    j = bpath->name;
	    change_row(lp_data, j, bpath->sense, bpath->rhs, bpath->range);
	 }
      }
   }

   //see whether it is feasible!

   retval = dual_simplex(lp_data, &iterd);
   close_lp_solver(lp_data);

   lp_data->mip = NULL;
   FREE(lp_data);
   bpath -= level;
   FREE(bpath);

   for(i=0; i<(level+1); i++){
      if(path[i])
	 path[i]=NULL;
   }

   FREE(path);

   return (retval);
#else
   printf("check_feasibility_new_rhs():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(0);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int trim_warm_tree(sym_environment *env, bc_node *n)
{
   int i, not_pruned = 0;

   /* There isn't anything to do if this is a leaf. */
   if (n->bobj.child_num == 0)
      return(0);

   /* There isn't anything to do if all children are pruned, and we are
      better off to go down if only one is not pruned. */
   for (i = n->bobj.child_num - 1; i >= 0; i--)
      if (n->children[i]->node_status != NODE_STATUS__PRUNED)
	 if (++not_pruned > 1)
	    break;
   if (not_pruned == 0)
      return(0);
   if (not_pruned == 1){
      for (i = n->bobj.child_num - 1; i >= 0; i--)
	 if (n->children[i]->node_status != NODE_STATUS__PRUNED){
	    trim_warm_tree(env, n->children[i]);
	    break;
	 }
      return(0);
   }

   /* So there are at least two not pruned. */
   for (i = n->bobj.child_num - 1; i >= 0; i--)
      if (n->children[i]->lower_bound + env->par.tm_par.granularity <
	  env->warm_start->ub)
	 break;

   /* if all children have high objval */
   if (i < 0){
      /* get rid of the children */
      for (i = n->bobj.child_num - 1; i >= 0; i--)
	 free_subtree(n->children[i]);
      /* free the children description */
      FREE(n->children);
      n->bobj.child_num = 0;
#ifndef MAX_CHILDREN_NUM
      FREE(n->bobj.sense);
      FREE(n->bobj.rhs);
      FREE(n->bobj.range);
      FREE(n->bobj.branch);
#endif
   }else{
      /* try to trim every child */
      for (i = n->bobj.child_num - 1; i >= 0; i--)
	 trim_warm_tree(env, n->children[i]);
   }
   return(0);
}

/*===========================================================================*/
/*===========================================================================*/

