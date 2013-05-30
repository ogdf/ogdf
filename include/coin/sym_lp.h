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

#ifndef _LP_H
#define _LP_H

#define COMPILING_FOR_LP

#include "symphony.h"
#include "sym_timemeas.h"
#include "sym_lp_params.h"
#include "sym_types.h"
#include "sym_lp_solver.h"
#include "sym_lp_u.h"

#ifdef COMPILE_IN_CG
#include "sym_cg.h"
#endif

#ifdef COMPILE_IN_CP
#include "sym_cp.h"
#endif

#ifdef COMPILE_IN_LP
#include "sym_tm.h"
#endif

/*===========================================================================*/
/*===========================================================================*/

typedef struct OUR_COL_SET{
   int           dual_feas;

   int           rel_lb;
   int          *rel_lb_ind;
   int           rel_ub;
   int          *rel_ub_ind;

   int           num_vars;
   int          *userind;

   double       *objx;
   double       *lb;
   double       *ub;
   int          *matbeg;
   int          *matind;
   double       *matval;
   int           nzcnt;
}our_col_set;

/*===========================================================================*/
typedef struct LP_PROB{
   double       str_time;

   int           proc_index;
   void         *user;

   lp_params     par;

   int           has_ub;
   double        ub;

   double        root_objval;

   int           phase;

   double        start_time;

   base_desc     base;

   branch_desc  *bdesc;        /* there are p->bc_level branch_desc's */
   bounds_change_desc *bnd_changes; /* there are p->bc_level bnd_change's */
   int           bdesc_size;

   int           master;
   int           draw_graph;
   int           tree_manager;
   int           cut_pool;
   int           cut_gen;

#ifdef COMPILE_IN_CG
   cg_prob      *cgp;
#endif
#ifdef COMPILE_IN_LP
   tm_prob      *tm;
#endif
   lp_sol        best_sol;
   double        obj[2];
   double        utopia[2];
   int           has_mc_ub;
   double        mc_ub;

   double        tt;
   node_times    comp_times;
   lp_stat_desc  lp_stat;

   node_desc    *desc;
   int           bc_index;
   int           bc_level;
   int           vars_at_ub;
   int           vars_at_lb;
   int           vars_deletable; /* a subset of vars at LB */

   int           dive;
   int           colgen_strategy;
   char          colgen_happened;
   char          colset_changed;

   int           iter_num;
   int           node_iter_num;
   int           bound_changes_in_iter;
   int           vars_recently_fixed_to_ub;
   LPdata       *lp_data;
   MIPdesc      *mip; /* Holds the MIP description when read in from MPS */

   double        last_gap;
   double       *obj_history;
   int           has_tailoff;
   int           obj_no_impr_iters;

   //double        str_br_impr_count; /* if we dont have 3 consequent
   //				       improvements, return to initials */

   /*========================================================================*\
    * The following fields refer to the cuts/rows arrived to the LP,
    * but not yet added
   \*========================================================================*/

   int           waiting_row_num;
   waiting_row **waiting_rows;
   int           waiting_rows_size;

   /*========================================================================*\
    * The following fields refer to the rows that once have been added to
    * the problem, but became slack
   \*========================================================================*/

   int         slack_cut_num;
   cut_data  **slack_cuts;
   int         slack_cuts_size;

   /* pseudo costs and reliability measures */
   double         *pcost_down;
   double         *pcost_up;
   int            *br_rel_down;
   int            *br_rel_up;
   int            *br_rel_cand_list;
   char            str_br_check;
   int            *br_rel_down_min_level;
   int            *br_rel_up_min_level;
   double          str_check_obj;
   int             str_check_trial;
   int             str_check_freq;
   int             str_check_cnt;
}lp_prob;

/*===========================================================================*/
/*============= LP general purpose functions (lp_genfunc.c) =================*/
/*===========================================================================*/

lp_prob *get_lp_ptr PROTO((lp_prob **lp_list));
int lp_initialize PROTO((lp_prob *p, int master_tid));
int process_chain PROTO((lp_prob *p));
int fathom_branch PROTO((lp_prob *p));
int check_bounds PROTO((lp_prob *p, int *termcode));
int fathom PROTO((lp_prob *p, int primal_feasible));
int repricing PROTO((lp_prob *p));
int bfind PROTO((int key, int *table, int size));
int collect_nonzeros PROTO((lp_prob *p, double *x, int *tind, double *tx));
int collect_fractions PROTO((lp_prob *p, double *x, int *tind, double *tx));
node_desc *create_explicit_node_desc PROTO((lp_prob *p));
int check_tailoff PROTO((lp_prob *p));
int round_solution PROTO((lp_prob *p, double *solution_value,
			  double *betterSolution));
int local_search PROTO((lp_prob *p, double *solution_value,
			double *col_solution, double *better_solution));
void lp_exit PROTO((lp_prob *p));
void lp_close PROTO((lp_prob *p));
int generate_cgl_cuts_new PROTO((lp_prob *p, int *num_cuts, cut_data ***cuts,
      int send_to_pool, int *bound_changes));
int should_use_cgl_generator PROTO ((lp_prob *p, int *should_generate,
      int which_generator, void *generator));
int generate_cgl_cut_of_type PROTO((lp_prob *p, int i, OsiCuts *cutlist_p,
         int *was_tried));
int check_and_add_cgl_cuts PROTO((lp_prob *p, int i, cut_data ***cuts, int *num_cuts, int *bound_changes, OsiCuts *cutlist, int send_to_pool));
int should_stop_adding_cgl_cuts PROTO((lp_prob *p, int i, int *should_stop));
int add_col_cuts PROTO((lp_prob *p, OsiCuts *cutlist, int *bound_changes));
int update_pcost PROTO ((lp_prob *p));
int str_br_bound_changes PROTO((lp_prob *p, int num_bnd_changes,
         double *bnd_val, int *bnd_ind, char *bnd_sense));

/*===========================================================================*/
/*======== LP functions related to variable management (lp_varfunc.c) =======*/
/*===========================================================================*/

void add_col_set PROTO((lp_prob *p, our_col_set *new_cols));
void colind_sort_extra PROTO((lp_prob *p));
void userind_sort_extra PROTO((lp_prob *p));
void tighten_bounds PROTO((lp_prob *p));
int save_root_reduced_costs(lp_prob *p);
int tighten_root_bounds(lp_prob *p);
our_col_set *price_all_vars PROTO((lp_prob *p));
int restore_lp_feasibility PROTO((lp_prob *p, our_col_set *new_cols));
void userind_sort_extra PROTO((lp_prob *p));
void colind_sort_extra PROTO((lp_prob *p));
int var_uind_comp PROTO((const void *v0, const void *v1));
int var_cind_comp PROTO((const void *v0, const void *v1));

/*===========================================================================*/
/*========== LP functions related to row management (lp_rowfunc.c) ==========*/
/*===========================================================================*/

int check_row_effectiveness PROTO((lp_prob *p));
void add_row_set PROTO((lp_prob *p, waiting_row **wrows, int length));
void add_new_rows_to_waiting_rows PROTO((lp_prob *p, waiting_row **new_rows,
					 int new_row_num));
void order_waiting_rows_based_on_sender PROTO((lp_prob *p));
int add_best_waiting_rows PROTO((lp_prob *p));
void add_waiting_rows PROTO((lp_prob *p, waiting_row **wrows,int add_row_num));
int waiting_row_comp PROTO((const void *wr0, const void *wr1));
int compute_violations PROTO((lp_prob *p,
			      int new_row_num, waiting_row **new_rows));
void compress_slack_cuts PROTO((lp_prob *p));

/*===========================================================================*/
/*================= LP branching functions (lp_branch.c) ====================*/
/*===========================================================================*/

void add_slacks_to_matrix PROTO((lp_prob *p, int cand_num,
				 branch_obj **candidates));
int add_violated_slacks PROTO((lp_prob *p, int cand_num,
			       branch_obj **candidates));
int select_branching_object PROTO((lp_prob *p, int *cuts,
				   branch_obj **can));
int should_continue_strong_branching PROTO((lp_prob *p, int i, int cand_num,
                                     double st_time, int total_iters,
                                     int *should_continue));
int strong_branch(lp_prob *p, int branch_var, double lb, double ub,
		  double new_lb, double new_ub, double *obj, int should_use_hot_starts,
                  int *termstatus, int *iterd);
int branch PROTO((lp_prob *p, int cuts));
int col_gen_before_branch PROTO((lp_prob *p, int *new_vars));

/*----------- Generic selection rules to be used by the user ----------------*/

void branch_close_to_half PROTO((lp_prob *p, int max_cand_num, int *cand_num,
				 branch_obj ***candidates));
void branch_close_to_half_and_expensive PROTO((lp_prob *p, int max_cand_num,
					       int *cand_num,
					       branch_obj ***candidates));
void branch_close_to_one_and_cheap PROTO((lp_prob *p, int max_cand_num,
					  int *cand_num,
					  branch_obj ***candidates));
/*===========================================================================*/
/*================ LP communication functions (lp_proccomm.c) ===============*/
/*===========================================================================*/

void check_ub PROTO((lp_prob *p));
int process_message PROTO((lp_prob *p, int r_bufid, int *pindex, int *pitnum));
void lp_process_ub_message PROTO((lp_prob *p));
int receive_active_node PROTO((lp_prob *p));
int receive_cuts PROTO((lp_prob *p, int first_lp, int no_more_cuts_count));
void send_node_desc PROTO((lp_prob *p, int node_type));
array_desc pack_array_desc_diff PROTO((array_desc *ad, array_desc *new_ad,
				       int *itmp));
basis_desc pack_basis_diff PROTO((node_desc *oldnode, node_desc *newnode,
				  char uind_type, char cutind_type,
				  int *itmp));
char pack_base_diff PROTO((int *size, int *oldstat, int *newstat, int *itmp));
char pack_extra_diff PROTO((array_desc *olddesc, int *oldstat,
			    array_desc *newdesc, int *newstat,
			    char oldbasis_type_in_tm, char newdesc_type_in_tm,
			    int *itmp, int *size));
void send_branching_info PROTO((lp_prob *p, branch_obj *can, char *action,
				int *keep));
void send_lp_is_free PROTO((lp_prob *p));
void send_cuts_to_pool PROTO((lp_prob *p, int eff_cnt_limit));
int add_bound_changes_to_desc PROTO((node_desc *new_tm_desc, lp_prob *p));
int update_cut_parameters(lp_prob *p);

/*===========================================================================*/
/*======================= Freeing things (lp_free.c) ========================*/
/*===========================================================================*/

void free_cut PROTO((cut_data **lpcut));
void free_waiting_row PROTO((waiting_row **wrow));
void free_waiting_rows PROTO((waiting_row **rows, int row_num));
void free_waiting_row_array PROTO((waiting_row ***rows, int row_num));
void free_cuts PROTO((cut_data **lpcuts, int cut_num));
void free_col_set PROTO((our_col_set **colset));
void free_candidate PROTO((branch_obj **cand));
void free_candidate_completely PROTO((branch_obj **cand));
void free_node_dependent PROTO((lp_prob *p));
void free_node_desc PROTO((node_desc **desc));
void free_lp PROTO((lp_prob *p));

/*===========================================================================*/
/*==================== LP wrapper functions (lp_wrapper.c) ==================*/
/*===========================================================================*/

int receive_lp_data_u PROTO((lp_prob *p));
void free_prob_dependent_u PROTO((lp_prob *p));
int comp_cut_name PROTO((const void *c0, const void *c1));
int create_subproblem_u PROTO((lp_prob *p));
int is_feasible_u PROTO((lp_prob *p, char branching, char is_last_iter));
void send_feasible_solution_u PROTO((lp_prob *p, int xlevel, int xindex,
				     int xiter_num, double lpetol,
				     double new_ub, int cnt, int *xind,
				     double *xval));
void display_lp_solution_u PROTO((lp_prob *p, int which_sol));
int select_candidates_u PROTO((lp_prob *p, int *cuts, int *new_vars,
			       int *cand_num, branch_obj ***candidates));
int compare_candidates_u PROTO((lp_prob *p, double oldobjval,
				branch_obj *best,branch_obj *can));
int select_child_u PROTO((lp_prob *p, branch_obj *can, char *action));
void print_branch_stat_u PROTO((lp_prob *p, branch_obj *can, char *action));
void add_to_desc_u PROTO((lp_prob *p, node_desc *desc));
int same_cuts_u PROTO((lp_prob *p, waiting_row *wrow1, waiting_row *wrow2));
void unpack_cuts_u PROTO((lp_prob *p, int from, int type,
			  int cut_num, cut_data **cuts,
			  int *new_row_num, waiting_row ***new_rows));
int send_lp_solution_u PROTO((lp_prob *p, int tid));
void logical_fixing_u PROTO((lp_prob *p));
int generate_column_u PROTO((lp_prob *p, int lpcutnum, cut_data **cuts,
			     int prevind, int nextind, int generate_what,
			     double *colval, int *colind, int *collen,
			     double *obj, double *ub, double *lb));
void print_stat_on_cuts_added_u PROTO((lp_prob *p, int added_rows));
void purge_waiting_rows_u PROTO((lp_prob *p));
int generate_cuts_in_lp_u PROTO((lp_prob *p));
int analyze_multicriteria_solution PROTO((lp_prob *p, int *indices,
					   double *values, int length,
					   double *true_objval, double etol,
					   char branching));
#endif
