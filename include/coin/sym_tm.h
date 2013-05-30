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

#ifndef _TREE_MANAGER
#define _TREE_MANAGER

#define COMPILING_FOR_TM

#include "symphony.h"
#include "sym_tm_params.h"
#include "sym_types.h"
#ifdef COMPILE_IN_CG
#include "sym_cg.h"
#endif
#ifdef COMPILE_IN_CP
#include "sym_cp.h"
#endif

/*===========================================================================*/

typedef struct PROCESS_SET{
   int        procnum;       /* the number of processes spawned */
   int       *procs;         /* the tid's of the processes */
   int        free_num;      /* how many of them are free */
   int       *free_ind;      /* the indices of the free ones */
}process_set;

/*===========================================================================*/

typedef struct TM_TEMP{
   int       *i;
   int        i_size;
   char      *c;
   int        c_size;
   double    *d;
   int        d_size;
}tm_temp;

/*===========================================================================*\
 * Problem data structure
\*===========================================================================*/

struct LP_PROB;

typedef struct TM_PROB{
   tm_params       par;
   int             master;
   int             has_ub;
   char            has_ub_estimate;
   double          start_time;
   double          ub;       /* the best global upper bound found */
   double          lb;       /* the best global lower bound known */
   lp_sol          best_sol;
   double          obj_offset; /* constant to be added to the objective value*/
   char            obj_sense;  /* objective sense*/
   double          ub_estimate;
   process_set     lp;
   process_set     cg;
   process_set     cp;

#ifdef COMPILE_IN_LP
   struct LP_PROB **lpp;
   int            opt_thread_num;
#ifdef COMPILE_IN_CG
   cg_prob      **cgp;
#endif
#endif
#ifdef COMPILE_IN_CP
   cut_pool     **cpp;
#endif

   int            *nodes_per_cp;        /* for each cut_pool it contains how
					   many nodes are assigned to it */
   int            *active_nodes_per_cp; /* same for active_nodes */

   bc_node        *rootnode;

   int             bvarnum;             /* number of base variables */
   int             bcutnum;             /* number of base constraints */

   int             phase;               /* the current phase */

   int             active_node_num;
   bc_node       **active_nodes;        /* contains a list of the nodes
					   currently being processed */

   int             samephase_candnum;   /* nodes still to be processed in */
   bc_node       **samephase_cand;      /* the current phase */
   int             samephase_cand_size;

   int             nextphase_candnum;   /* nodes to be processed next phase*/
   bc_node       **nextphase_cand;
   int             nextphase_cand_size;

   /* The list of cuts that were active in at least one search tree node */
   int             cut_num;
   int             allocated_cut_num;
   cut_data      **cuts;

   problem_stat    stat;

   node_times      comp_times;         /* keeps track of the computation times
			                  for the problem */
   lp_stat_desc    lp_stat;
   rc_desc        *reduced_costs;

   /* pseudo costs and reliability measures */
   double         *pcost_down;
   double         *pcost_up;
   int            *br_rel_down;
   int            *br_rel_up;
   int            *br_rel_cand_list;
   int            *br_rel_down_min_level;
   int            *br_rel_up_min_level;

   /* some temporary stuff */
   bc_node      ***rpath;
   int            *rpath_size;
   branch_desc   **bpath;
   int            *bpath_size;

#ifdef TRACE_PATH
   int             feas_sol_size;
   int            *feas_sol;
#endif

   tm_temp         tmp;
   /* solution pool */
   sp_desc  *sp;
}tm_prob;

/*===========================================================================*/
/*==================== TM basic functions (tm_func.c) =======================*/
/*===========================================================================*/

int tm_initialize PROTO((tm_prob *tm, base_desc *base,
			 node_desc *root_desc));
int solve PROTO((tm_prob *tm));
void print_tree_status PROTO((tm_prob *tm));
void calculate_widths PROTO((bc_node *node, int *widths));
int start_node PROTO((tm_prob *tm, int thread_num));
bc_node *del_best_node PROTO((tm_prob *tm));
void insert_new_node PROTO((tm_prob *tm, bc_node *new_node));
int node_compar PROTO((int rule, bc_node *node0, bc_node *node1));
int assign_pool PROTO((tm_prob *tm, int oldpool, process_set *pools,
		       int *active_nodes_per_pool, int *nodes_per_pool));
int generate_children PROTO((tm_prob *tm, bc_node *node, branch_obj *bobj,
			     double *objval, int *feasible, char *action,
			     int olddive, int *keep, int new_branching_cut));
char shall_we_dive PROTO((tm_prob *tm, double objval));
int purge_pruned_nodes PROTO((tm_prob *tm, bc_node *node, int category));
int find_process_index PROTO((process_set *pset, int tid));
void mark_lp_process_free PROTO((tm_prob *tm, int lp, int cp));
int add_cut_to_list PROTO((tm_prob *tm, cut_data *cut));
void install_new_ub PROTO((tm_prob *tm, double new_ub, int opt_thread_num,
			   int bc_index, char branching, int feasible));
int find_tree_lb PROTO((tm_prob *tm));

/*--------------- Function related to merging descriptions ------------------*/

void merge_descriptions PROTO((node_desc *old_node, node_desc *new_node));
void merge_base_stat PROTO((double_array_desc *dad,
			    double_array_desc *moddad));
void merge_extra_array_and_stat PROTO((array_desc *ad, double_array_desc *dad,
				       array_desc *modad,
				       double_array_desc *moddad));
void merge_double_array_descs PROTO((double_array_desc *dad,
				     double_array_desc *newdad));
void merge_arrays PROTO((array_desc *array, array_desc *adesc));
void modify_list PROTO((array_desc *origad, array_desc *modad));
void modify_list_and_stat PROTO((array_desc *origad, int *origstat,
				 array_desc *modad,
				 double_array_desc *moddad));

/*--------------- Functions related to two-phase algorithm ------------------*/

int tasks_before_phase_two PROTO((tm_prob *tm));
int trim_subtree PROTO((tm_prob *tm, bc_node *n));
int mark_subtree PROTO((tm_prob *tm, bc_node *n));
void propagate_nf_status PROTO((bc_node *n, int nf_status));

/*------------------ Functions related to logging ---------------------------*/

void write_log_files PROTO((tm_prob *tm));
int write_pruned_nodes PROTO((tm_prob *tm, bc_node *node));
int write_node PROTO((bc_node *node, char *file, FILE *f, char append));
int write_subtree PROTO((bc_node *node, char *file, FILE* f, char append,
			 int logging));
int write_tm_cut_list PROTO((tm_prob *tm, char *file, char append));
int write_tm_info PROTO((tm_prob *tm, char *file, FILE* f, char append));
int write_base PROTO((base_desc *base, char *file, FILE* f, char append));
int read_node PROTO((tm_prob *tm, bc_node *node, FILE *f, int **children));
int read_subtree PROTO((tm_prob *tm, bc_node *node, FILE *f));
int read_tm_cut_list PROTO((tm_prob *tm, char * file));
int read_tm_info PROTO((tm_prob *tm, FILE *f));
int read_base PROTO((base_desc *base, FILE *f));

/*-------------- Functions related to freeing and closing -------------------*/

void free_tm PROTO((tm_prob *tm));
void free_subtree PROTO((bc_node *n));
void free_tree_node PROTO((bc_node *n));
void free_basis PROTO((basis_desc *basis));
int tm_close PROTO((tm_prob *tm, int termcode));

/*===========================================================================*/
/*============== TM communication functions (tm_proccomm.c) =================*/
/*===========================================================================*/

process_set start_processes PROTO((tm_prob *tm,
				   int procnum, char *procname, int procdebug,
				   int machnum, char **mach));
void stop_processes PROTO((process_set *pset));
char processes_alive PROTO((tm_prob *tm));
void send_active_node PROTO((tm_prob *tm, bc_node *node, int colgen_strat,
			     int thread_num));
void receive_node_desc PROTO((tm_prob *tm, bc_node *n));
void process_branching_info PROTO((tm_prob *tm, bc_node *node));
char process_messages PROTO((tm_prob *tm, int r_bufid));
void process_ub_message PROTO((tm_prob *tm));
void unpack_cut_set PROTO((tm_prob *tm, int sender, int cutnum,
			   row_data *rows));
int receive_lp_timing PROTO((tm_prob *tm));

void sym_catch_c PROTO((int num));
int merge_bound_changes PROTO((bounds_change_desc **bnd_change_ptr,
                               bounds_change_desc  *p_bnd_change));
#endif
