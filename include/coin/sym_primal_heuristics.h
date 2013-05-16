/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY MILP Solver Framework.                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Ashutosh Mahajan                               */
/*                                                                           */
/* (c) Copyright 2006-2011 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

/*
 * TODO:
 * change ifdef _FEASI.. to ifdef _HEURISTICS
 */

#ifndef _PRIMAL_HEURISTICS_H
#define _PRIMAL_HEURISTICS_H
#include "sym_lp_solver.h"
#include "sym_lp.h"
#include "sym_types.h"

/* feasibility pump */
typedef struct FP_VARS {
   char          is_bin;
   char          is_int;
   int           xplus;
   int           xminus;
}FPvars;

typedef struct FP_DATA {
   FPvars      **fp_vars;       /* an array of fp_vars */
   int           n0;            /* no. of vars in orignial lp */
   int           m0;
   int           n;             /* no. of vars in pumping lp */
   int           m;             /* no. of constraints in pumping lp */
   int           iter;
   int           numNonBinInts;
   int           numInts;
   int          *index_list;
   int         **x_bar_ind;     /* array containing previous x_bars */
   double      **x_bar_val;     /* array containing previous x_bars */
   int          *x_bar_len;     /* rounded x_lp */
   double       *alpha_p;       /* previous alphas */
   double       *x_lp;          /* solution of pumpling lp */
   double       *x_ip;          /* rounded x_lp */
   double       *mip_obj;       /* normalized original obj */
   double       *obj;           /* obj function for pumping lp */
   char          can_check_sos; /* whether we can check sos rows while fixing bin vars */
   char         *sos_row_filled;/*to keep track of the sos variables while flipping */
   char         *sos_var_fixed_zero;/*to keep track of the sos variables while flipping */
   double        norm_c;        /* norm of mip_obj */
   double        alpha;
   double        alpha_decr;
   int           verbosity;
   double        flip_fraction;
   double        norm;
   int           iterd;
}FPdata;

/*  solution pool */
int sp_add_solution PROTO((lp_prob *p, int cnt, int *indices, double *values, double obj_value, int bc_index));
int sp_delete_solution PROTO((sp_desc *sp, int position));
int sp_is_solution_in_sp PROTO((lp_prob *p, int cnt, int *indices, double *values, double obj_value));
#ifdef COMPILE_IN_LP
int sp_initialize(tm_prob *tm);
#endif
int sp_free_sp(sp_desc *sp);

/* feasibility pump */
int feasibility_pump (lp_prob *p, char *found_better_solution, double &solution_value, double *betterSolution);
int fp_round (lp_prob *p, FPdata *fp_data, LPdata *lp_data);
int fp_is_feasible (LPdata *lp_data, const CoinPackedMatrix *matrix, const double *r_low, const double *r_up, FPdata *fp_data, char *is_feasible );
int fp_initialize_lp_solver(lp_prob *p, LPdata *new_lp_data, FPdata *fp_data);
int fp_solve_lp(LPdata *lp_data, FPdata *fp_data, char *is_feasible) ;
int fp_should_call_fp(lp_prob *p, int branching, int *should_call,
      char is_last_iter);
int fp_add_obj_row(LPdata *new_lp_data, int n, const double *obj, double rhs);
int fp_can_sos_var_fix(lp_prob *p, FPdata *fp_data, int ind, int *filled_row_count);
int fp_fix_sos_var(lp_prob *p, FPdata *fp_data, int ind);
#endif
