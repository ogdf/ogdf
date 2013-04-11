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

#ifndef _CP_U_H
#define _CP_U_H

#include "sym_proto.h"
#include "sym_types.h"

/*===========================================================================*/
/*====================== User supplied functions ============================*/
/*===========================================================================*/

int user_receive_cp_data PROTO((void **user));
int user_free_cp PROTO((void **user));
int user_prepare_to_check_cuts PROTO((void *user, int varnum, int *indices,
				      double *values));
int user_check_cut PROTO((void *user, double lpetol, int varnum, int *indices,
			  double *values, cut_data *cut, int *is_violated,
			  double *quality));
int user_finished_checking_cuts PROTO((void *user));
int user_receive_lp_solution_cp PROTO((void *user));

#endif
