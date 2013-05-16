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

#ifndef _LPSOLVER_H
#define _LPSOLVER_H

#include "sym_proto.h"
#include "sym_types.h"
#include "sym_lp_params.h"

#define LP_MAX_ITER 9999999

#ifdef __CPLEX__

/*****************************************************************************/
/*******              here are the definitions for CPLEX               *******/
/*****************************************************************************/

#include <cplex.h>

void CPX_check_error PROTO((const char *erring_func));

#elif defined(__OSL__)

/*****************************************************************************/
/*******              here are the definitions for OSL                 *******/
/*****************************************************************************/

#include <ekk_c_api.h>

void OSL_check_error PROTO((const char *erring_func));

#elif defined(__OSI_CPLEX__) || defined(__OSI_OSL__) || defined(__OSI_CLP__) \
|| defined(__OSI_XPRESS__) || defined(__OSI_SOPLEX__) || defined(__OSI_VOL__) \
|| defined(__OSI_DYLP__) || defined (__OSI_GLPK__)

/*****************************************************************************/
/*******              here are the definitions for OSI                 *******/
/*****************************************************************************/

#include "OsiSolverInterface.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinMpsIO.hpp"
#include "CoinLpIO.hpp"

#ifdef USE_CGL_CUTS
#include "OsiCuts.hpp"
#include "CglCutGenerator.hpp"
#include "CglLiftAndProject.hpp"
#include "CglSimpleRounding.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglClique.hpp"
#include "CglGomory.hpp"
#include "CglOddHole.hpp"
#include "CglKnapsackCover.hpp"
#include "CglProbing.hpp"
#include "CglFlowCover.hpp"
#include "CglTwomir.hpp"
#include "CglLandP.hpp"
#include "CglRedSplit.hpp"
#endif

#ifdef __OSI_CPLEX__
#include "OsiCpxSolverInterface.hpp"
typedef OsiCpxSolverInterface OsiXSolverInterface;
#endif

#ifdef __OSI_OSL__
#include "OsiOslSolverInterface.hpp"
typedef OsiOslSolverInterface OsiXSolverInterface;
#endif

#ifdef __OSI_CLP__
#include "OsiClpSolverInterface.hpp"
typedef OsiClpSolverInterface OsiXSolverInterface;
#endif

#ifdef __OSI_XPRESS__
#include "OsiXprSolverInterface.hpp"
typedef OsiXprSolverInterface OsiXSolverInterface;
#endif

#ifdef __OSI_SOPLEX__
#include "OsiSpxSolverInterface.hpp"
typedef OsiSpxSolverInterface OsiXSolverInterface;
#endif

#ifdef __OSI_VOL__
#include "OsiVolSolverInterface.hpp"
typedef OsiVolSolverInterface OsiXSolverInterface;
#endif

#ifdef __OSI_DYLP__
#include "OsiDylpSolverInterface.hpp"
typedef OsiDylpSolverInterface OsiXSolverInterface;
#endif

#ifdef __OSI_GLPK__
#include "OsiGlpkSolverInterface.hpp"
typedef OsiGlpkSolverInterface OsiXSolverInterface;
#endif

#else

#error ###################################
#error # Undefined or unknown LP solver.
#error # Please edit SYMPHONY/Makefile
#error # and define LP_SOLVER properly.
#error ###################################

#endif

/*****************************************************************************/
/*******                  end LP solver definitions                    *******/
/*****************************************************************************/

#ifdef USE_GLPMPL
extern "C"
{
   #include "glpk.h"
}
#endif

/* Temporary storage */

typedef struct TEMPORARY{
   char      *c;           /* max(2m,n) */
   int       *i1;          /* 3m+2n */
   int       *i2;          /* m */
   double    *d;           /* max(2m,2n) */
   void     **p1;          /* m */
   void     **p2;          /* m */

   char      *cv;          /* variable */
   int        cv_size;
   int       *iv;          /* variable (>= */
   int        iv_size;
   double    *dv;          /* variable */
   int        dv_size;
}temporary;

/* The LP solver data */

typedef struct LPDATA{
   /* First, the problem pointers */
#ifdef __CPLEX__
   CPXENVptr  cpxenv;
   CPXLPptr   lp;
#endif
#ifdef __OSL__
   EKKContext *env;
   EKKModel   *lp;
#endif
#if defined(__OSI_CPLEX__) || defined(__OSI_OSL__) || defined(__OSI_CLP__) \
|| defined(__OSI_XPRESS__) || defined(__OSI_SOPLEX__) || defined(__OSI_VOL__) \
|| defined(__OSI_DYLP__) || defined (__OSI_GLPK__)
   OsiXSolverInterface * si;
#endif
   double     lpetol;
   char       lp_is_modified;
   char       col_set_changed;
   double     objval;
   int        termcode;
   MIPdesc   *mip;
   int        n;           /* number of columns without slacks */
   int        maxn;
   int        m;           /* number of rows */
   int        maxm;
   int        nz;          /* number of nonzeros */
   int        maxnz;       /* space is allocated for this many nonzeros */
   double    *random_hash;
   double    *heur_solution; /* space for heur solution */

   char       ordering;    /* COLIND_AND_USERIND_ORDERED, COLIND_ORDERED or
			      USERIND_ORDERED */
   var_desc **vars;        /* maxn */ /* BB */

   int        not_fixed_num;
   int       *not_fixed;
   int        nf_status;

   char      *status;      /* maxn */ /* BB */
   double    *x;           /* maxn */ /* BB */
   double    *dj;          /* maxn */ /* BB */
   double    *dualsol;     /* maxm */ /* BB */
   double    *slacks;      /* maxm */
   double    *ub;
   double    *lb;

   row_data  *rows;      /* maxm */

   temporary  tmp;
#ifdef PSEUDO_COSTS
   double     *pseudo_costs_one;
   double     *pseudo_costs_zero;
#endif
   int         lp_count;
   cgl_params  cgl;

}LPdata;

/*****************************************************************************/
/*******                    common definitions                         *******/
/*****************************************************************************/

double dot_product PROTO((double *val, int *ind, int collen, double *col));
void free_lp_arrays PROTO((LPdata *lp_data));
void free_mip_desc PROTO((MIPdesc *mip));
void size_lp_arrays PROTO((LPdata *lp_data, char do_realloc, char set_max,
			     int row_num, int col_num, int nzcnt));
void open_lp_solver PROTO((LPdata *lp_data));
void close_lp_solver PROTO((LPdata *lp_data));
void load_lp_prob PROTO((LPdata *lp_data, int scaling, int fastmip));
int reset_lp_prob PROTO ((LPdata *lp_data, int scaling, int fastmip));
int save_lp PROTO((LPdata *lp_data));
void unload_lp_prob PROTO((LPdata *lp_data));
void load_basis PROTO((LPdata *lp_data, int *cstat, int *rstat));
void refactorize PROTO((LPdata *lp_data));
void add_rows PROTO((LPdata *lp_data, int rcnt, int nzcnt, double *rhs,
		     char *sense, int *rmatbeg, int *rmatind,double *rmatval));
void add_cols PROTO((LPdata *lp_data, int ccnt, int nzcnt, double *obj,
		     int *cmatbeg, int *cmatind, double *cmatval,
		     double *lb, double *ub, char *where_to_move));
void change_row PROTO((LPdata *lp_data, int row_ind,
		       char sense, double rhs, double range));
void change_col PROTO((LPdata *lp_data, int col_ind,
		       char sense, double lb, double ub));
int initial_lp_solve PROTO((LPdata *lp_data, int *iterd));
int dual_simplex PROTO((LPdata *lp_data, int *iterd));
int solve_hotstart PROTO((LPdata *lp_data, int *iterd));
int mark_hotstart PROTO((LPdata *lp_data));
int unmark_hotstart PROTO((LPdata *lp_data));
void btran PROTO((LPdata *lp_data, double *col));
void get_binvcol PROTO((LPdata *lp_data, int j, double *col));
void get_binvrow PROTO((LPdata *lp_data, int i, double *row));
void get_basis PROTO((LPdata *lp_data, int *cstat, int *rstat));
void set_obj_upper_lim PROTO((LPdata *lp_data, double lim));
void set_itlim PROTO((LPdata *lp_data, int itlim));
void set_itlim_hotstart PROTO((LPdata *lp_data, int itlim));
void get_column PROTO((LPdata *lp_data, int j,
		       double *colval, int *colind, int *collen, double *cj));
void get_row PROTO((LPdata *lp_data, int i,
		    double *rowval, int *rowind, int *rowlen,
		    double *rowub, double *rowlb));
int get_proof_of_infeas PROTO((LPdata *lp_data, int *infind));
void get_x PROTO((LPdata *lp_data));
void get_dj_pi PROTO((LPdata *lp_data));
void get_slacks PROTO((LPdata *lp_data));
void change_range PROTO((LPdata *lp_data, int rowind, double value));
void change_rhs PROTO((LPdata *lp_data,
		       int rownum, int *rhsind, double *rhsval));
void change_sense PROTO((LPdata *lp_data, int cnt, int *index, char *sense));
void change_bounds PROTO((LPdata *lp_data,
			  int cnt, int *index, char *lu, double *bd));
void change_lbub PROTO((LPdata *lp_data, int j, double lb, double ub));
void change_ub PROTO((LPdata *lp_data, int j, double ub));
void change_lb PROTO((LPdata *lp_data, int j, double lb));
void get_ub PROTO((LPdata *lp_data, int j, double *ub));
void get_lb PROTO((LPdata *lp_data, int j, double *lb));
void get_bounds PROTO((LPdata *lp_data));
void get_objcoef PROTO((LPdata *lp_data, int j, double *objcoef));
void get_objcoeffs(LPdata *lp_data);
void change_objcoeff(LPdata *lp_data, const int* indexFirst,
      const int* indexLast, double *coeffs);
void get_rhs_rng_sense(LPdata *lp_data);
int copy_lp_data(LPdata *lp_data, LPdata *new_data);
void delete_rows PROTO((LPdata *lp_data, int deletable, int *free_rows));
void delete_rows_with_ind PROTO((LPdata *lp_data, int deletable, int *rowind));
int delete_cols PROTO((LPdata *lp_data, int delnum, int *delstat));
void release_var PROTO((LPdata *lp_data, int j, int where_to_move));
void free_row_set PROTO((LPdata *lp_data, int length, int *index));
void constrain_row_set PROTO((LPdata *lp_data, int length, int *index));
int read_mps PROTO((MIPdesc *mip, char *infile, char *probname));
int read_lp PROTO((MIPdesc *mip, char *infile, char *probname));
void write_mps PROTO((LPdata *lp_data, char *fname));
void write_mip_desc_mps PROTO((MIPdesc *mip, char *fname));
void write_mip_desc_lp PROTO((MIPdesc *mip, char *fname));
void write_sav PROTO((LPdata *lp_data, char *fname));
#ifdef USE_CGL_CUTS
void generate_cgl_cuts(LPdata *lp_data, int *num_cuts, cut_data ***cuts,
		       char send_to_pool, int bc_index, int bc_level,
                       int node_iter_limit, int max_cuts_before_resolve,
                       double ub, int *bnd_changes,
                       lp_stat_desc *lp_stat, node_times *comp_times,
                       int verbosity);
int check_cuts(OsiCuts &cutlist, LPdata *lp_data, int bc_level, int
      *num_cuts, cut_data ***cuts, char send_to_pool, int *bnd_changes,
      lp_stat_desc *lp_stat, node_times *compe_times, int verbosity);
int should_generate_this_cgl_cut(int cut_num, int max_cuts_before_resolve,
      int generation_flag, int freq, int bc_level, int bc_index,
      int cuts_in_root, int *should_generate);
/*
void generate_cgl_cuts PROTO((LPdata * lp_data, int *num_cuts,
			      cut_data ***cuts, char send_to_pool,
			      int is_rootnode, lp_stat_desc *lp_stat,
                              node_times *comp_times, int verbosity));
*/
#endif
#ifdef USE_GLPMPL
int read_gmpl PROTO((MIPdesc *mip, char *modelfile, char *datafile,
		     char *probname));
#endif
#endif
