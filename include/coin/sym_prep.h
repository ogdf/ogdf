/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Menal Guzelsoy                                 */
/*                                                                           */
/* (c) Copyright 2006-2011 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
/* last modified: January 09, menal*/

#ifndef SYM_PREP__H
#define SYM_PREP__H

/* return codes of the functions*/

#define PREP_FUNC_SUCCESS   0
#define PREP_FUNC_ERROR    -1

#include "symphony.h"
#include "sym_types.h"
#include "sym_constants.h"
#include "sym_prep_params.h"

#ifdef INF
#undef INF
#endif
#ifdef SYM_INFINITY
#define INF SYM_INFINITY
#else
#define INF 1e20
#endif

#define PREP_QUIT(f)          \
  ((f != PREP_UNMODIFIED && f != PREP_MODIFIED) ? TRUE : FALSE)

/*===========================================================================*/

/* for sr internal use */
/* Menal's single-row (sr) stuff */

#define SR_NO_UPDATES       0
#define SR_BOUNDS_UPDATED   1
#define SR_INFEAS           2

#define POS_VAL             0
#define ZERO_VAL            1
#define NEG_VAL             2

#define SR_MIN              0
#define SR_MAX              1

#define RND_FLOOR           0
#define RND_CEIL            1

#define VAR_NEW             0
#define VAR_IN_BOUND        1

#define UB_SIDE             0
#define LB_SIDE             1
#define BOTH_SIDE           2

/* status of a variable in sr problem */
#define SR_VAR_IN           0
#define SR_VAR_IN_FIXED_UB  1
#define SR_VAR_IN_FIXED_LB  2
#define SR_VAR_IN_FRAC      3
#define SR_VAR_FIXED_UB     4
#define SR_VAR_FIXED_LB     5

/*===========================================================================*/

/* modification types on a variable */

#define FIX_NO_BOUND        0
#define FIX_BINARY          1
#define FIX_OTHER           2
#define FIX_FIXABLE         3
#define IMPROVE_UB          4
#define IMPROVE_LB          5
#define IMPROVE_COEF        6
#define FIX_AGGREGATE       7

/* for a range of variables */
#define FIX_ROW_LB          8
#define FIX_ROW_UB          9

/*===========================================================================*/
/* data structure to keep the statistics */
/*===========================================================================*/
typedef struct PREP_STATS
{
   int rows_deleted;
   int vars_fixed;
   int coeffs_nulled;
   int bounds_integerized;
   int vars_aggregated;
   int vars_integerized;

   /* regarding coeffs changes and bounds tightening */
   int coeffs_changed;
   char *nz_coeff_changed;

   int bounds_tightened;

   /* regarding uboundedness and infeasiblity */
   int col_infeas_ind;
   int row_infeas_ind;
   int col_unbound_ind;
   int col_numeric_ind;
}prep_stats;

/*===========================================================================*/
/* Single Row Relaxation data structure */
/* Under development */
/*===========================================================================*/
typedef struct SRDESC{

   int prob_type;
   char sense;
   double rhs;

   int max_n;    /* all variables which are not fixed yet */
   double *obj_max;
   double *matval_max;
   double *ratio_max;
   int *matind_max;
   char *reversed_max;
   //  int *ratio_type_max;
   double ub_offset;
   double rhs_max;
   double sum_c_max;
   double sum_a_max;

   char ub_updated;
   double ub;

   int min_n;
   double *obj_min;
   double *matval_min;
   double *ratio_min;
   int *matind_min;
   char *reversed_min;
   //   int *ratio_type_min;
   double lb_offset;
   double rhs_min;
   double sum_c_min;
   double sum_a_min;

   char lb_updated;
   double lb;

   /* for sorting purposes */
   int * fixed_ind;
   int * tmp_ind;

   /* for variable fixing, bound tightening purpose*/

   int * var_stat_max;
   int * var_stat_min;

   double *var_obj_max;
   double *var_matval_max;

   double *var_obj_min;
   double *var_matval_min;

   double *var_min_opt; /* for solving the same problem for
			   each variable fixed
			*/
   double *var_max_opt;

}SRdesc;

/*===========================================================================*/
/* Preprocessing data structure  */
/*===========================================================================*/
typedef struct PREPDesc
{
   MIPdesc * mip;
   MIPdesc * orig_mip;
   prep_stats stats;
   prep_params params;

   /* for logical fixing */
   int impl_limit;
   //int impl_var_cnt; /* fixed ones */
   IMPlist *list; /* the list under inspection */
   int      impl_col_ind;
   prep_stats impl_stats;
   int      impl_row_cnt;
   int      impl_var_cnt;
   char     *impl_vars;

   ROWinfo *impl_rows;
   COLinfo *impl_cols;

   double *impl_ub;
   double *impl_lb;

   char *ulist_checked;
   char *llist_checked;

   /* trying single/aggr row relaxations to improve bounds*/
   int max_sr_cnt;
   int max_aggr_cnt;
   SRdesc *sr; /* for 'L', 'G' constraints */
   SRdesc *d_sr; /* additionally, for 'E' constraints */
   /* for subproblems checking purposes */
   char *rows_checked;
   double alloc_time;

   /* will need for sorting etc */
   int * user_col_ind;
   int * user_row_ind;
   double alloc2_time;
   double impl_array_time;
   double impl_cols_time;
   double impl_rows_time;
}PREPdesc;

/*===========================================================================*/
/* Data structure to keep relevant info of a column */

/* User accessible environment to manage preprocessing
   -added for user in another package
   -not used in symphony */
/*=========================================================================*/

typedef struct PREP_ENVIRONMENT{
   PREPdesc * P;
   prep_stats stats;
   prep_params params;
   int termcode;
}prep_environment;

/*===========================================================================*/

/* presolve the MIP formulation stored in the current environment */
int sym_presolve(sym_environment *env);

/* some data structures in root description are initialized before calling
 * preprocessor. update these after the preprocessor has changed the problem
 * size. */
int prep_update_rootdesc(sym_environment *env);

/*load a problem through MIP model description arrays*/
int prep_load_problem(prep_environment *prep, int numcols, int numrows,
		      int *start, int *index, double *value,
		      double *collb, double *colub, char *is_int,
		      double *obj, double obj_offset, char *rowsen,
		      double *rowrhs, double *rowrng, char make_copy);

/*==========================================================================*/
/*==========================================================================*/
/*==========================================================================*/

/* internal functions */

/*presolve the desc */
int prep_solve_desc(PREPdesc *P);

/* initialize the presolve description */
int prep_initialize_mipinfo(PREPdesc *P);

/* get the row oriented matrix description*/
int prep_fill_row_ordered(PREPdesc *P);

/* final touchup on the description*/
int prep_cleanup_desc(PREPdesc *P);

/* integerize the variable bounds */
int prep_integerize_bounds(PREPdesc *P);

/* integerize a continuous variable */
int prep_integerize_var(PREPdesc *P, int col_ind);

/* the main presolve loop*/
int prep_basic(PREPdesc *P);

/* try to improve the bounds of a variable*/
int prep_improve_variable(PREPdesc *P, int col_ind, int row_ind, int a_loc,
			  int dive_level, char check_improve, char impl_mode,
			  char use_sr_bounds,
			  double sr_ub, double sr_lb, int use_mip);

/* check if the given row is redundant */
int prep_check_redundancy(PREPdesc *P, int row_ind, char use_sr_bounds,
			  double sr_ub, double sr_lb, char impl_mode,
			  int dive_level);

/* if a column is modified, then update the model*/
int prep_modified_cols_update_info(PREPdesc *P, int col_cnt, int *col_start,
				   int row_ind, int dive_level,
				   double fixed_bound,  int fix_type,
				   char check_redundancy, char impl_mode);

/* for the unbounded variables, check if we can tighten their bounds*/

int prep_force_row_bounds(PREPdesc *P, int row_ind, int col_ind, int a_loc);

/* update the matrix when a row is proved to be redundant*/
int prep_deleted_row_update_info(MIPdesc *mip, int row_ind);

/* try to find duplicate rows and columns */
int prep_delete_duplicate_rows_cols(PREPdesc *P, char check_rows,
				    char check_cols);
/* utility functions */
void prep_sos_fill_var_cnt(PREPdesc *P);
void prep_sos_fill_row(ROWinfo *row, int alloc_size, int size,
		       int *ind);

double prep_rnd_integral(double val, double etol, char rnd_type);
int  prep_get_row_bounds(MIPdesc *mip, int r_ind, double etol);
char prep_is_equal(double lval, double rval, double etol);
char prep_is_integral(double val, double etol);

/* reporting functions */
int prep_declare_fixed_var(int col_ind, char *name, double fixed_bound);
int prep_declare_redundant_row(ROWinfo row, int row_ind, char sense,
			       double rhs);
int prep_declare_coef_change(int row_ind, int col_ind,
			     char *name, double a_val,
			     double rhs);
int prep_report(PREPdesc *P, int termcode);

/* implications - under development*/
int prep_add_to_impl_list(IMPlist *list, int ind, int fix_type,
			  double val);
int prep_initialize_impl_lists(PREPdesc *P);

/* experimental - under development */
int prep_solve_sr_rlx(PREPdesc *P, int row_cnt, int *row_indices);

/*==========================================================================*/
/*==========================================================================*/
/*==========================================================================*/

/* functions to solve single row relaxations of the model*/
/* ---- under development ----- not entirely tested*/


/* initialize/allocate SR description */
void sr_initialize(SRdesc **sr, int n);
void sr_allocate(SRdesc **sr, int n);

/* solve the single row (indexed by row_ind) relaxation*/
int sr_solve_bounded_prob(PREPdesc *P, SRdesc *sr, SRdesc *d_sr,
			  int obj_ind, int row_ind,
			  int *r_matbeg, int *r_matind, double *r_matval,
			  COLinfo *cols, double *ub, double *lb, double etol);
/* internal functions: */
/* add a new column to the problem: the description of column is passed in
   through function arguments*/
int sr_add_new_col(SRdesc *sr, SRdesc *d_sr, double c_val, double a_val,
		   int col_ind, char col_var_type, double col_ub,
		   double col_lb, char sense, int col_type,
		   int col_bound_type);
/* add a new column to the problem: the description of column is passed in
   through function arguments - here we know that it is bounded*/
int sr_add_new_bounded_col(SRdesc *sr, double c_val, double a_val,
			   int col_ind,
			   double rhs_ub_offset, double rhs_lb_offset,
			   double obj_ub_offset, double obj_lb_offset,
			   double col_ub, double col_lb, int obj_sense,
			   char var_type);
/* helper functions */
int sr_find_opt_bounded(PREPdesc *P, SRdesc *sr, int obj_ind,
			double *ub, double *lb);

int sr_solve_open_prob(PREPdesc *P, SRdesc *sr, int obj_ind,
		       int row_ind, int *r_matbeg,
		       int *r_matind, double *r_matval, COLinfo *cols,
		       double *ub, double *lb, double etol);

void free_prep_desc(PREPdesc *P);
void free_sr_desc(SRdesc *sr);
void free_imp_list(IMPlist **list);
#endif
