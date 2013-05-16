/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2005-2011 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _SYM_API_H
#define _SYM_API_H

#define COMPILING_FOR_MASTER

#ifdef PROTO
#undef PROTO
#endif
#define PROTO(x) x

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************                  Return Values                       **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*----------------------- Global return codes -------------------------------*/
#define FUNCTION_TERMINATED_NORMALLY      0
#define FUNCTION_TERMINATED_ABNORMALLY   -1
#define ERROR__USER                      -100

/*-------------- Return codes for sym_parse_comand_line() -------------------*/
#define ERROR__OPENING_PARAM_FILE        -110
#define ERROR__PARSING_PARAM_FILE        -111

/*----------------- Return codes for sym_load_problem() ---------------------*/
#define ERROR__READING_GMPL_FILE         -120
#define ERROR__READING_WARM_START_FILE   -121
#define ERROR__READING_MPS_FILE          -122
#define ERROR__READING_LP_FILE           -123

/*-------------------- Return codes for sym_solve() -------------------------*/
#define TM_NO_PROBLEM                     225
#define TM_NO_SOLUTION                    226
#define TM_OPTIMAL_SOLUTION_FOUND         227
#define TM_TIME_LIMIT_EXCEEDED            228
#define TM_NODE_LIMIT_EXCEEDED            229
#define TM_TARGET_GAP_ACHIEVED            230
#define TM_FOUND_FIRST_FEASIBLE           231
#define TM_FINISHED                       232
#define TM_UNFINISHED                     233
#define TM_FEASIBLE_SOLUTION_FOUND        234
#define TM_SIGNAL_CAUGHT                  235
#define TM_UNBOUNDED                      236
#define PREP_OPTIMAL_SOLUTION_FOUND       237
#define PREP_NO_SOLUTION                  238
#define TM_ERROR__NO_BRANCHING_CANDIDATE -250
#define TM_ERROR__ILLEGAL_RETURN_CODE    -251
#define TM_ERROR__NUMERICAL_INSTABILITY  -252
#define TM_ERROR__COMM_ERROR             -253
#define TM_ERROR__USER                   -275
#define PREP_ERROR                       -276

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************                  General Constants                   **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef ANYONE
#define ANYONE   -1
#endif
#ifndef ANYTHING
#define ANYTHING -1
#endif

#define DSIZE sizeof(double)
#define ISIZE sizeof(int)
#define CSIZE sizeof(char)

#ifndef BITSPERBYTE
#define	BITSPERBYTE 8
#endif
#ifndef BITS
#define	BITS(type) (BITSPERBYTE * (int)sizeof (type))
#endif

#ifdef  HIBITI
#undef  HIBITI
#endif
#define	HIBITI (1U << (BITS(int) - 1))
#ifdef  MAXINT
#undef  MAXINT
#endif
#define	MAXINT ((int)(~(HIBITI)))
#ifdef MAXDOUBLE
#undef MAXDOUBLE
#endif
#define	MAXDOUBLE 1.79769313486231570e+308

#define SYM_INFINITY                 1e20

#define BIG_DBL                      1e40

#define SYM_MINIMIZE                 0
#define SYM_MAXIMIZE                 1

/*--------------------- return values for user-written functions ------------*/
#define USER_ERROR              -5
#define USER_SUCCESS            -4
#define USER_NO_PP              -3
#define USER_AND_PP             -2
#define USER_DEFAULT            -1

/*------------ search order options for multi-criteria problems -------------*/
#define MC_FIFO                  0
#define MC_LIFO                  1

/*------------  warm_starting options for multi-criteria problems -------------*/
#define MC_WS_UTOPIA_FIRST               0
#define MC_WS_UTOPIA_BOTH_FIXED          1
#define MC_WS_UTOPIA_BOTH                2
#define MC_WS_BEST_CLOSE                 3

/*------------------------ compare_candidates -------------------------------*/
#define BIGGEST_DIFFERENCE_OBJ   0
#define LOWEST_LOW_OBJ           1
#define HIGHEST_LOW_OBJ          2
#define LOWEST_HIGH_OBJ          3
#define HIGHEST_HIGH_OBJ         4
#define HIGH_LOW_COMBINATION     9

/*--------------------------- select_child ----------------------------------*/
#define PREFER_LOWER_OBJ_VALUE   0
#define PREFER_HIGHER_OBJ_VALUE  1

/*-------------------- generate_cuts_in_lp defaults -------------------------*/
#define GENERATE_CGL_CUTS                  20
#define DO_NOT_GENERATE_CGL_CUTS           21

/*-------------------- xxx_cuts_generation_levels ---------------------------*/
#define DO_NOT_GENERATE        -1
#define GENERATE_DEFAULT        0
#define GENERATE_IF_IN_ROOT     1
#define GENERATE_ONLY_IN_ROOT   2
#define GENERATE_ALWAYS         3
#define GENERATE_PERIODICALLY   4

/*------------------------- node selection rules ----------------------------*/
#define LOWEST_LP_FIRST       0
#define HIGHEST_LP_FIRST      1
#define BREADTH_FIRST_SEARCH  2
#define DEPTH_FIRST_SEARCH    3
#define BEST_FIRST_SEARCH     4
#define DEPTH_FIRST_THEN_BEST_FIRST 5

/*-------------------------- diving_strategy --------------------------------*/
#define BEST_ESTIMATE         0
#define COMP_BEST_K           1
#define COMP_BEST_K_GAP       2

/*--------------- parameter values for feasibility pump heuristic -----------*/
#define SYM_FEAS_PUMP_DEFAULT    1       /* use fp using the default rules   */
#define SYM_FEAS_PUMP_REPEATED   2       /* use fp till the end of solve     */
#define SYM_FEAS_PUMP_TILL_SOL   3       /* use fp till a solution is found  */
#define SYM_FEAS_PUMP_DISABLE   -1       /* dont use fp */

typedef struct MIPDESC MIPdesc;
typedef struct WARM_START_DESC warm_start_desc;
typedef struct SYM_ENVIRONMENT sym_environment;

/*===========================================================================*/
/*===================== Interface functions (master.c) ======================*/
/*===========================================================================*/

void sym_version PROTO((void));
sym_environment *sym_open_environment PROTO((void));
int sym_set_defaults PROTO((sym_environment *env));
int sym_parse_command_line PROTO((sym_environment *env, int argc,
				  char **argv));
int sym_set_user_data PROTO((sym_environment *env, void *user));
int sym_get_user_data PROTO((sym_environment *env, void **user));
int sym_read_mps PROTO((sym_environment *env, char *infile));
int sym_read_lp PROTO((sym_environment *env, char *infile));
int sym_read_gmpl PROTO((sym_environment *env, char *modelfile,
			 char *datafile));
int sym_write_mps PROTO((sym_environment *env, char *infile));
int sym_write_lp PROTO((sym_environment *env, char *infile));

int sym_load_problem PROTO((sym_environment *env));
int sym_find_initial_bounds PROTO((sym_environment *env));

int sym_solve PROTO((sym_environment *env));
int sym_warm_solve PROTO((sym_environment *env));
int sym_mc_solve PROTO((sym_environment *env));

int sym_create_permanent_cut_pools PROTO((sym_environment *env, int *cp_num));
int sym_close_environment PROTO((sym_environment *env));
int sym_explicit_load_problem PROTO((sym_environment *env, int numcols,
				     int numrows, int *start, int *index,
				     double *value, double *collb,
				     double *colub, char *is_int, double *obj,
				     double *obj2, char *rowsen,
				     double *rowrhs, double *rowrng,
				     char make_copy));

int sym_is_abandoned PROTO((sym_environment *env));
int sym_is_proven_optimal PROTO((sym_environment *env));
int sym_is_proven_primal_infeasible PROTO((sym_environment *env));
int sym_is_iteration_limit_reached PROTO((sym_environment *env));
int sym_is_time_limit_reached PROTO((sym_environment *env));
int sym_is_target_gap_achieved PROTO((sym_environment *env));

int sym_get_status PROTO((sym_environment *env));
int sym_get_num_cols PROTO((sym_environment *env, int *numcols));
int sym_get_num_rows PROTO((sym_environment *env, int *numrows));
int sym_get_num_elements PROTO((sym_environment *env, int *numelems));
int sym_get_col_lower PROTO((sym_environment *env, double *collb));
int sym_get_col_upper PROTO((sym_environment *env, double *colub));
int sym_get_row_sense PROTO((sym_environment *env, char *rowsen));
int sym_get_rhs PROTO((sym_environment *env, double *rowrhs));
int sym_get_matrix PROTO((sym_environment *env, int *nz, int *matbeg,
			  int *matind, double *matval));
int sym_get_row_range PROTO((sym_environment *env, double *rowrng));
int sym_get_row_lower PROTO((sym_environment *env, double *rowlb));
int sym_get_row_upper PROTO((sym_environment *env, double *rowub));
int sym_get_obj_coeff PROTO((sym_environment *env, double *obj));
int sym_get_obj2_coeff PROTO((sym_environment *env, double *obj2));
int sym_get_obj_sense PROTO((sym_environment *env, int *sense));

int sym_is_continuous PROTO((sym_environment *env, int index, int *value));
int sym_is_binary PROTO((sym_environment *env, int index, int *value));
int sym_is_integer PROTO((sym_environment *env, int index, char *value));

double sym_get_infinity PROTO(());

int sym_get_col_solution PROTO((sym_environment *env, double *colsol));
int sym_get_row_activity PROTO((sym_environment *env, double *rowact));
int sym_get_obj_val PROTO((sym_environment *env, double *objval));
int sym_get_primal_bound PROTO((sym_environment *env, double *ub));
int sym_get_iteration_count PROTO((sym_environment *env, int *numnodes));

int sym_set_obj_coeff PROTO((sym_environment *env, int index, double value));
int sym_set_obj2_coeff PROTO((sym_environment *env, int index, double value));
int sym_set_col_lower PROTO((sym_environment *env, int index, double value));
int sym_set_col_upper PROTO((sym_environment *env, int index, double value));
int sym_set_row_lower PROTO((sym_environment *env, int index, double value));
int sym_set_row_upper PROTO((sym_environment *env, int index, double value));
int sym_set_row_type PROTO((sym_environment *env, int index, char rowsense,
			     double rowrhs, double rowrng));
int sym_set_obj_sense PROTO((sym_environment *env, int sense));
int sym_set_col_solution PROTO((sym_environment *env, double * colsol));
int sym_set_primal_bound PROTO((sym_environment *env, double bound));
int sym_set_continuous PROTO((sym_environment *env, int index));
int sym_set_integer PROTO((sym_environment *env, int index));
int sym_set_col_names PROTO((sym_environment *env, char **colname));
int sym_add_col PROTO((sym_environment *env, int numelems, int *indices,
		       double *elements, double collb, double colub,
		       double obj, char is_int, char *name));
int sym_add_row PROTO((sym_environment *env, int numelems, int *indices,
		       double *elements, char rowsen, double rowrhs,
		       double rowrng));
int sym_delete_cols PROTO((sym_environment *env, int num, int * indices));
int sym_delete_rows PROTO((sym_environment *env, int num, int * indices));

int sym_write_warm_start_desc PROTO((warm_start_desc *ws, char *file));
warm_start_desc *sym_read_warm_start PROTO((char *file));

void sym_delete_warm_start PROTO((warm_start_desc *ws));
warm_start_desc *sym_get_warm_start PROTO((sym_environment *env,
					   int copy_warm_start));

int sym_set_warm_start PROTO((sym_environment *env, warm_start_desc *ws));

int sym_set_int_param PROTO((sym_environment *env, const char *key, int value));
int sym_set_dbl_param PROTO((sym_environment *env, const char *key, double value));
int sym_set_str_param PROTO((sym_environment *env, const char *key, const char *value));

int sym_get_int_param PROTO((sym_environment *env, const char *key, int *value));
int sym_get_dbl_param PROTO((sym_environment *env, const char *key, double *value));
int sym_get_str_param PROTO((sym_environment *env, const char *key, char **value));

int sym_get_lb_for_new_rhs PROTO((sym_environment *env, int cnt,
				  int *new_rhs_ind, double *new_rhs_val,
				  double *lb_for_new_rhs));
int sym_get_ub_for_new_rhs PROTO((sym_environment *env, int cnt,
				  int *new_rhs_ind, double *new_rhs_val,
				  double *ub_for_new_rhs));
#if 0
int sym_get_lb_for_new_obj PROTO((sym_environment *env, int cnt,
				  int *new_obj_ind, double *new_obj_val,
				  double *lb_for_new_obj));
#endif
int sym_get_ub_for_new_obj PROTO((sym_environment *env, int cnt,
				  int *new_obj_ind, double *new_obj_val,
				  double *ub_for_new_obj));

warm_start_desc *sym_create_copy_warm_start PROTO((warm_start_desc * ws));
MIPdesc *sym_create_copy_mip_desc PROTO((sym_environment *env));
MIPdesc *sym_get_presolved_mip_desc PROTO((sym_environment *env));
sym_environment * sym_create_copy_environment PROTO((sym_environment *env));

int sym_test PROTO((sym_environment *env, int *test_status));

#endif
