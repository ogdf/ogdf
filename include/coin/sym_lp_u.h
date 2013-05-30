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

#ifndef LP_U_H
#define LP_U_H

#include "sym_proto.h"
#include "sym_types.h"
#include "sym_lp_solver.h"

/*----------- Generic selection rules to be used by the user ----------------*/

void branch_close_to_half PROTO((int max_cand_num, int *cand_num,
				 branch_obj ***candidates));
void branch_close_to_half_and_expensive PROTO((int max_cand_num, int *cand_num,
					       branch_obj ***candidates));
void branch_close_to_one_and_cheap PROTO((int max_cand_num, int *cand_num,
					  branch_obj ***candidates));

/*===========================================================================*/
/*========================= User supplied functions =========================*/
/*===========================================================================*/

int user_receive_lp_data PROTO((void **user));
int user_free_lp PROTO((void **user));
int user_create_subproblem PROTO((void *user, int *indices, MIPdesc *mip,
				  int *maxn, int *maxm, int *maxnz));
int user_is_feasible PROTO((void *user, double lpetol, int varnum,
			    int *indices, double *values, int *feasible,
			    double *true_objval, char branching,
			    double *heur_solution));
int user_send_feasible_solution PROTO((void *user, double lpetol, int varnum,
				       int *indices, double *values));
int user_display_lp_solution PROTO((void *user, int which_sol, int varnum,
				    int *indices, double *values));
int user_shall_we_branch PROTO((void *user, double lpetol, int cutnum,
				int slacks_in_matrix_num,
				cut_data **slacks_im_matrix, int slack_cut_num,
				cut_data **slack_cuts, int varnum,
				var_desc **vars, double *x, char *status,
				int *cand_num, branch_obj ***candidates,
				int *action));
int user_select_candidates PROTO((void *user, double lpetol, int cutnum,
				  int slacks_in_matrix_num,
				  cut_data **slacks_im_matrix,
				  int slack_cut_num, cut_data **slack_cuts,
				  int varnum, var_desc **vars, double *x,
				  char *status, int *cand_num,
				  branch_obj ***candidates, int *action,
				  int bc_level));
int user_compare_candidates PROTO((void *user, branch_obj *can1,
				   branch_obj *can2, double ub,
				   double granularity, int *which_is_better));
int user_select_child PROTO((void *user, double ub, branch_obj *can,
			     char *action));
int user_print_branch_stat PROTO((void *user, branch_obj *can, cut_data *cut,
				  int varnum, var_desc **vars, char *action));
int user_add_to_desc PROTO((void *user, int *desc_size, char **desc));
int user_same_cuts PROTO((void *user, cut_data *cut1, cut_data *cut2,
			  int *same_cuts));
int user_unpack_cuts PROTO((void *user, int from, int type, int varnum,
			    var_desc **vars, int cutnum, cut_data **cuts,
			    int *new_row_num, waiting_row ***new_rows));
int user_send_lp_solution PROTO((void *user, int varnum, var_desc **vars,
				 double *x, int where));
int user_logical_fixing PROTO((void *user, int varnum, var_desc **vars,
			       double *x, char *status, int *fixed_num));
int user_generate_column PROTO((void *user, int generate_what, int cutnum,
				cut_data **cuts, int prevind, int nextind,
				int *real_nextind, double *colval, int *colind,
				int *collen, double *obj, double *lb,
				double *ub));
int user_generate_cuts_in_lp PROTO((void *user, LPdata *lp_data, int varnum,
				    var_desc **vars,double *x,int *new_cut_num,
				    cut_data ***new_cuts));
int user_print_stat_on_cuts_added PROTO((void *user, int rownum,
					 waiting_row **rows));
int user_purge_waiting_rows PROTO((void *user, int rownum,
				   waiting_row **rows, char *deleten));

#endif
