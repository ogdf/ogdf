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

#ifndef _CG_U_H
#define _CG_U_H

#include "sym_types.h"

int cg_add_explicit_cut PROTO((int nzcnt, int *indices, double *values,
			       double rhs, double range, char sense,
			       char send_to_cp, int *num_cuts, int *alloc_cuts,
			       cut_data ***cuts));
int cg_add_user_cut PROTO((cut_data *new_cut, int *num_cuts, int *alloc_cuts,
			   cut_data ***cuts));
int cg_send_cut PROTO((cut_data *new_cut, int *num_cuts, int *alloc_cuts,
		       cut_data ***cuts));

/*===========================================================================*/
/*======================= User supplied functions ===========================*/
/*===========================================================================*/

int user_receive_cg_data PROTO((void **user, int dg_id));
int user_free_cg PROTO((void **user));
int user_find_cuts PROTO((void *user, int varnum, int iter_num, int level,
			  int index, double objval, int *indices,
			  double *values, double ub, double lpetol,
			  int *num_cuts, int *alloc_cuts, cut_data ***cuts));
int user_receive_lp_solution_cg PROTO((void *user));
#ifdef CHECK_CUT_VALIDITY
int user_check_validity_of_cut PROTO((void *user, cut_data *new_cut));
#endif

#endif
