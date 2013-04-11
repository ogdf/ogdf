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

#ifndef _CUT_GEN_H
#define _CUT_GEN_H

#include "symphony.h"
#include "sym_types.h"
#include "sym_cg_params.h"
#include "sym_cg_u.h"
#include "sym_lp_solver.h"

/*===========================================================================*/

/*stores the data needed by the cut_generator*/
typedef struct CG_PROB{
   int            proc_index;
   void          *user;
   int            msgtag;
   int            master;
   int            draw_graph;    /* the tid of DrawGraph */
   int            tree_manager;  /* the tid of the tree manager */
   cg_params      par;           /* the parameters for the cut generator */
   char           has_ub;        /* is there an upper bound */
   double         ub;            /* the current best upper bound if there
				    is one */
   double         tt;
   lp_sol         cur_sol;
#ifdef COMPILE_IN_CG
   int           cuts_to_add_num;
   cut_data    **cuts_to_add;
   int           cuts_to_add_size;
#endif
}cg_prob;

/*===========================================================================*/
/*==================== CG basic functions (cg_func.c) =======================*/
/*===========================================================================*/

cg_prob *get_cg_ptr PROTO((cg_prob **cg_list));
void cg_initialize PROTO((cg_prob *p, int master_tid));
void cg_close PROTO((cg_prob * p));
cut_data *create_explicit_cut PROTO((int nzcnt, int *indices, double *values,
				     double rhs, double range, char sense,
				     char send_to_cp));
int cg_add_explicit_cut PROTO((int nzcnt, int *indices, double *values,
			       double rhs, double range, char sense,
			       char send_to_cp, int *num_cuts, int *alloc_cuts,
			       cut_data ***cuts));
int cg_add_user_cut PROTO((cut_data *new_cut, int *num_cuts, int *alloc_cuts,
			   cut_data ***cuts));

/*===========================================================================*/
/*=============== CG communication functions (cg_proccomm.c) ================*/
/*===========================================================================*/

int cg_process_message PROTO((cg_prob *p, int r_bufid));
int cg_send_cut PROTO((cut_data *new_cut, int *num_cuts, int *alloc_cuts,
		       cut_data ***cuts));

/*===========================================================================*/
/*==================== LP wrapper functions (cg_wrapper.c) ==================*/
/*===========================================================================*/

int receive_cg_data_u PROTO((cg_prob *p));
int receive_lp_solution_cg_u PROTO((cg_prob *p));
int free_cg_u PROTO((cg_prob *p));
int find_cuts_u PROTO((cg_prob *p, LPdata *lp_data, int *num_cuts));
int check_validity_of_cut_u PROTO((cg_prob *p, cut_data *new_cut));

#endif
