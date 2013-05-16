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

#ifndef MASTER_U_H
#define MASTER_U_H

#include "sym_proto.h"

/*===========================================================================*/
/*======================= User supplied functions ===========================*/
/*===========================================================================*/

void user_usage PROTO((void));
int user_initialize PROTO((void **user));
int user_free_master PROTO((void **user));
int user_readparams PROTO((void *user, char *filename, int argc, char **argv));
int user_io PROTO((void *user));
int user_init_draw_graph PROTO((void *user, int dg_id));
int user_start_heurs PROTO((void *user, double *ub, double *ub_estimate));
int user_initialize_root_node PROTO((void *user, int *basevarnum, int **basevars,
				     int *basecutnum, int *extravarnum,
				     int **extravars, char *obj_sense,
				     double *obj_offset, char ***col_names,
				     int *colgen_strat));
int user_receive_feasible_solution PROTO((void *user, int msgtag, double cost,
					  int numvars, int *indices,
					  double *values));
int user_send_lp_data PROTO((void *user, void **user_lp));
int user_send_cg_data PROTO((void *user, void **user_cg));
int user_send_cp_data PROTO((void *user, void **user_cp));
int user_display_solution PROTO((void *user, double lpetol, int varnum,
				 int *indices, double *values, double objval));
int user_process_own_messages PROTO((void *user, int msgtag));
int user_send_feas_sol PROTO((void *user, int *feas_sol_size, int **feas_sol));
int user_ws_update_cuts PROTO((void *user, int *size, char **coef, double * rhs,
			       char *sense, char type, int new_col_num,
			       int change_type));


#endif
