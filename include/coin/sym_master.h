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

#ifndef _MASTER_H
#define _MASTER_H

#define COMPILING_FOR_MASTER

#include <stdio.h>

#include "symphony.h"
#include "sym_types.h"
#include "sym_macros.h"
#include "sym_master_params.h"
#ifdef COMPILE_IN_TM
#include "sym_tm.h"
#endif

/*===========================================================================*\
 * The problem data structure contains the data for a problem instance, as
 * well as some of the tours that have been generated.
\*===========================================================================*/

typedef struct SYM_ENVIRONMENT{
   void            *user;
   int              my_tid;
   int              tm_tid;
   int              dg_tid;
   params           par;         /* problem parameters */
   prob_times       comp_times;  /* keeps track of the computation times for
				    the problem */
   int              has_ub;
   double           ub;
   lp_sol           best_sol;
   int              has_mc_ub;
   double           mc_ub;
   double           obj[2];
   double           utopia[2];
   char             has_ub_estimate;
   double           ub_estimate;
   double           lb;
   double           obj_offset;

   MIPdesc         *mip; /*For holding the description when read in from MPS
			   - also the working copy */

   MIPdesc         *orig_mip; /*For holding the original description if
				presolve is used*/
   MIPdesc         *prep_mip; /* For holding the presolved description if
				 presolve is used*/

   char             probname[81];

   base_desc       *base;
   node_desc       *rootdesc;

   int              termcode;

   warm_start_desc *warm_start;

   double          mc_time;

#ifdef COMPILE_IN_TM
   tm_prob         *tm;
#ifdef COMPILE_IN_CP
   cut_pool       **cp;
#endif
#endif
}sym_environment;

/*===========================================================================*/
/*=================== Master I/O functions (readparams.c) ===================*/
/*===========================================================================*/

void usage PROTO((void));
int parse_command_line PROTO((sym_environment *env, int argc, char **argv));
void read_string PROTO((char *target, char *line, int maxlen));
void print_statistics PROTO((node_times *tim, problem_stat *stat,
                            lp_stat_desc *lp_stat, double ub,
			     double lb, double initial_time,
			     double start_time, double finish_time,
			     double obj_offset, char obj_sense, int has_ub,
                             sp_desc *solpool));

/*===========================================================================*/
/*=============== Master wrapper functions (master_wrapper.c) ===============*/
/*===========================================================================*/

int initialize_u PROTO((sym_environment *env));
int free_master_u PROTO((sym_environment *env));
int readparams_u PROTO((sym_environment *env, int argc, char **argv));
int io_u PROTO((sym_environment *env));
int init_draw_graph_u PROTO((sym_environment *env));
int start_heurs_u PROTO((sym_environment *env));
int display_solution_u PROTO((sym_environment *env, int thread_num));
int initialize_root_node_u PROTO((sym_environment *env));
int receive_feasible_solution_u PROTO((sym_environment *env, int msgtag));
int send_lp_data_u PROTO((sym_environment *env, int sender));
int send_cg_data_u PROTO((sym_environment *env, int sender));
int send_cp_data_u PROTO((sym_environment *env, int sender));
int send_sp_data_u PROTO((sym_environment *env, int sender));
int process_own_messages_u PROTO((sym_environment *env, int msgtag));

/*===========================================================================*/
/*=================== Master helper functions (master_func.c) ===============*/
/*===========================================================================*/

int resolve_node PROTO((sym_environment *env, bc_node * node));
int update_tree_bound PROTO((sym_environment *env, bc_node *root, int *cut_num,
			      int *cut_ind, char *cru_vars, int change_type));
void register_cuts PROTO((bc_node *root, int *cut_num,  int *cuts_ind));
void update_node_desc PROTO((sym_environment *env, bc_node *root,
			     int change_type));
void update_branching_decisions PROTO((sym_environment *env, bc_node *root,
			      int change_type));
void check_trim_tree PROTO((sym_environment *env, bc_node *root, int *cut_num,
			    int *cuts_ind, int change_type));
void cut_ws_tree_index PROTO((sym_environment *env, bc_node *root, int index,
			      problem_stat * stat, int change_type));
void cut_ws_tree_level PROTO((sym_environment *env, bc_node *root, int level,
			      problem_stat * stat, int change_type));
void ws_free_subtree PROTO((sym_environment *env, bc_node *root,
			    int change_type, int check_solution, int update_stats));
void check_better_solution PROTO((sym_environment * env, bc_node *root,
				  int delete_node, int change_type));
int copy_node PROTO((bc_node * n_to, bc_node *n_from));
int copy_tree PROTO((bc_node *root_to, bc_node *root_from));
int read_node PROTO((bc_node * node, FILE *f));
int read_tree PROTO((bc_node * root, FILE *f));
int write_node PROTO((bc_node *node, FILE *f));
int write_tree PROTO((bc_node *root, FILE *f));

int set_param PROTO((sym_environment *env,  char *line));

warm_start_desc *create_copy_warm_start PROTO((warm_start_desc * ws));
MIPdesc *create_copy_mip_desc PROTO((MIPdesc *mip));
sym_environment *create_copy_environment PROTO((sym_environment *env));

double get_lb_for_new_rhs PROTO((bc_node *root, MIPdesc *mip, int cnt,
				 int *ind, double *val));
double get_ub_for_new_rhs PROTO((bc_node *root, MIPdesc *mip, int cnt,
				 int *ind, double *val));
#if 0
double get_lb_for_new_obj PROTO((bc_node *root, MIPdesc *mip, int cnt,
				 int *ind, double *val));
#endif
double get_ub_for_new_obj PROTO((bc_node *root, MIPdesc *mip, int cnt,
				 int *ind, double *val));
int check_feasibility_new_rhs PROTO((bc_node * node, MIPdesc * mip,
					int cnt, int *ind, double *val));
int trim_warm_tree PROTO((sym_environment *env, bc_node *n));
#endif
