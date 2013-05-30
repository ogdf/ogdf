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

#define COMPILING_FOR_MASTER

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#ifdef __PVM__
#include <pvmtev.h>
#endif

#include "symphony.h"
#include "SymConfig.h"
#include "sym_proccomm.h"
#include "sym_timemeas.h"
#include "sym_messages.h"
#include "sym_macros.h"
#include "sym_qsort.h"
#include "sym_pack_cut.h"
#include "sym_pack_array.h"
#include "sym_master.h"
#include "sym_master_u.h"
#include "sym_lp_solver.h"
#include "sym_primal_heuristics.h"
#include "sym_prep.h"
#ifdef COMPILE_IN_TM
#include "sym_tm.h"
#ifdef COMPILE_IN_LP
#include "sym_lp.h"
#endif
#endif

#ifndef TEV_INIT_MASK
/* We must have pvm3.4 where it is called TEV_MASK_INIT */
#  define TEV_INIT_MASK(m)  TEV_MASK_INIT(m)
#  define TEV_SET_MASK(m,k)  TEV_MASK_SET(m,k)
#  define TEV_MCAST0  TEV_MCAST
#  define TEV_RECV0   TEV_RECV
#  define TEV_SEND0   TEV_SEND
#  define TEV_NRECV0  TEV_NRECV
#endif

/*===========================================================================*\
 * This file implements the SYMPHONY callable API
\*===========================================================================*/

/*===========================================================================*/
/*===========================================================================*/

void sym_version(void)
{
   printf("\n");
   printf("==  Welcome to the SYMPHONY MILP Solver \n");
   printf("==  Copyright 2000-2011 Ted Ralphs and others \n");
   printf("==  All Rights Reserved. \n");
   printf("==  Distributed under the Eclipse Public License 1.0 \n");
   if (strcmp(SYMPHONY_VERSION, "trunk")){
      printf("==  Version: %s \n", SYMPHONY_VERSION);
   }else{
      printf("==  Version: Trunk (unstable) \n");
   }
   printf("==  Build Date: %s \n", __DATE__);
#ifdef SYMPHONY_SVN_REV
   printf("==  Revision Number: %d \n", SYMPHONY_SVN_REV);
#endif
   printf("\n");
}

/*===========================================================================*/
/*===========================================================================*/

sym_environment *sym_open_environment()
{
   sym_environment *env;
#if (!defined(COMPILE_IN_LP) || !defined(COMPILE_IN_CG) || \
   !defined(COMPILE_IN_CP)) && defined(__PVM__)
   int xpvm_tid;
   Pvmtmask trace_mask;
#endif

   setvbuf(stdout, (char *)NULL, _IOLBF, 2);

   env = (sym_environment *) calloc(1, sizeof(sym_environment));

#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)

   env->my_tid = register_process();   /* Enroll this process */

#ifdef __PVM__
   pvm_catchout(stdout); /* Tells PVM to treat all output from the children of
			    this process as output from this process itself*/
#endif
#endif

#if 0
   sym_version();
#endif

   if (initialize_u(env) == FUNCTION_TERMINATED_NORMALLY){
      return(env);
   }else{
      FREE(env);
      return(NULL);
   }

   /* This next set of commands has to be executed if we want to create a PVM
      trace file for viewing in xpvm (this is a very slow process) */

#if (!defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)) && defined(__PVM__)
   if (env->par.pvm_trace){
      if ((xpvm_tid = pvm_gettid((char *)"xpvm", 0)) > 0){
	 pvm_setopt(PvmSelfTraceTid, xpvm_tid);
	 pvm_setopt(PvmSelfTraceCode, 666);
	 pvm_setopt(PvmSelfOutputTid, xpvm_tid);
	 pvm_setopt(PvmSelfOutputCode, 667);
	 pvm_setopt(PvmTraceTid, xpvm_tid);
	 pvm_setopt(PvmTraceCode, 666);
	 pvm_setopt(PvmOutputTid, xpvm_tid);
	 pvm_setopt(PvmOutputCode, 667);
	 TEV_INIT_MASK(trace_mask);
	 TEV_SET_MASK(trace_mask, TEV_MCAST0);
	 TEV_SET_MASK(trace_mask, TEV_RECV0);
	 TEV_SET_MASK(trace_mask, TEV_SEND0);
	 TEV_SET_MASK(trace_mask, TEV_NRECV0);
	 pvm_settmask(PvmTaskSelf, trace_mask);
	 pvm_settmask(PvmTaskChild, trace_mask);
      }else{
	 PVM_ERROR(xpvm_tid);
      }
   }
#endif

}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_defaults(sym_environment *env)
{
   int termcode = 0;

   tm_params *tm_par = &env->par.tm_par;
   lp_params *lp_par = &env->par.lp_par;
   cg_params *cg_par = &env->par.cg_par;
   cp_params *cp_par = &env->par.cp_par;
   dg_params *dg_par = &env->par.dg_par;
   prep_params *prep_par = &env->par.prep_par;


   /************************* Global defaults ********************************/
   env->ub = 0;
   env->has_ub = FALSE;
   env->lb = -MAXDOUBLE;
   env->termcode = TM_NO_PROBLEM;
   env->par.verbosity = 0;
   env->par.random_seed = 17;
   env->par.tm_machine_set = FALSE;
   env->par.dg_machine_set = FALSE;
   strcpy(env->par.tm_exe, "symphony_tm");
#ifdef COMPILE_IN_LP
   strcat(env->par.tm_exe, "_lp");
#ifdef COMPILE_IN_CG
   strcat(env->par.tm_exe, "_cg");
#endif
#endif
#ifdef COMPILE_IN_CP
   strcat(env->par.tm_exe, "_cp");
#endif
   strcpy(env->par.dg_exe, "symphony_dg");
   env->par.tm_debug = 0;
   env->par.dg_debug = 0;
   env->par.pvm_trace = 0;
   env->par.do_branch_and_cut = 1;
   env->par.do_draw_graph = FALSE;
   env->par.use_permanent_cut_pools = FALSE;
   env->par.multi_criteria = FALSE;
   env->par.mc_binary_search_tolerance = 0;
   env->par.mc_compare_solution_tolerance = .001;
   env->par.mc_search_order = MC_FIFO;
   env->par.mc_warm_start = FALSE;
   env->par.mc_warm_start_rule = MC_WS_UTOPIA_FIRST;
   env->par.trim_warm_tree = FALSE;
   env->par.test = FALSE;
   /************************** treemanager defaults **************************/
   tm_par->verbosity = 0;
   tm_par->granularity = 0.000001;
   strcpy(tm_par->lp_exe, "symphony_lp");
#ifdef COMPILE_IN_CG
   strcat(tm_par->lp_exe, "_cg");
#endif
   strcpy(tm_par->cg_exe, "symphony_cg");
   strcpy(tm_par->cp_exe, "symphony_cp");
   tm_par->lp_debug = 0;
   tm_par->cg_debug = 0;
   tm_par->cp_debug = 0;
   tm_par->max_active_nodes = 1;
   tm_par->max_cp_num = 1;
   tm_par->lp_mach_num = 0;
   tm_par->lp_machs = NULL;
   tm_par->cg_mach_num = 0;
   tm_par->cg_machs = NULL;
   tm_par->cp_mach_num = 0;
   tm_par->cp_machs = NULL;

   tm_par->use_cg = FALSE;
   tm_par->random_seed = 17;
   tm_par->unconditional_dive_frac = 0;
   tm_par->diving_strategy = BEST_ESTIMATE;
   tm_par->diving_k = 1;
   tm_par->diving_threshold = 0.05;
   tm_par->node_selection_rule = LOWEST_LP_FIRST;
   tm_par->keep_description_of_pruned = DISCARD;

   tm_par->warm_start = FALSE;
   tm_par->warm_start_node_ratio = 0.0;
   tm_par->warm_start_node_limit = INT_MAX; //(int)SYM_INFINITY;
   tm_par->warm_start_node_level_ratio = 0.0;
   tm_par->warm_start_node_level = INT_MAX; //(int)SYM_INFINITY;

   tm_par->logging = NO_LOGGING;
   tm_par->logging_interval = 1800;
   tm_par->vbc_emulation = NO_VBC_EMULATION;
   tm_par->price_in_root = FALSE;
   tm_par->trim_search_tree = FALSE;
   tm_par->colgen_strat[0] = (FATHOM__DO_NOT_GENERATE_COLS__DISCARD  |
			      BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
   tm_par->colgen_strat[1] = (FATHOM__DO_NOT_GENERATE_COLS__DISCARD  |
			      BEFORE_BRANCH__DO_NOT_GENERATE_COLS);
   tm_par->not_fixed_storage_size = 2048;
   tm_par->time_limit = lp_par->time_limit = -1.0;
   tm_par->node_limit = -1;
   tm_par->gap_limit = -1.0;
   // tm_par->gap_limit = 0.0;
   tm_par->find_first_feasible = FALSE;
   tm_par->sensitivity_analysis = FALSE;

   /************************** lp defaults ***********************************/
   lp_par->verbosity = 0;
   lp_par->granularity = tm_par->granularity;
   lp_par->use_cg = tm_par->use_cg;
   lp_par->set_obj_upper_lim = TRUE;
   lp_par->do_primal_heuristic = TRUE;
   lp_par->scaling = -1; /* CPLEX'ism ... don't scale */
   lp_par->fastmip = 1; /* CPLEX'ism ... set it to 1 */
   lp_par->should_warmstart_chain = TRUE; /* see header file for description */
   lp_par->should_reuse_lp = FALSE; /* see header file for description */
#ifdef COMPILE_IN_LP
   lp_par->should_reuse_lp = TRUE; /* see header file for description */
#endif
#ifdef _OPENMP
   lp_par->should_reuse_lp = FALSE; /* see header file for description */
#endif
#ifdef USE_SYM_APPLICATION
   lp_par->should_reuse_lp = FALSE; /* see header file for description */
#endif
   lp_par->try_to_recover_from_error = TRUE;
   lp_par->problem_type = ZERO_ONE_PROBLEM;
   lp_par->keep_description_of_pruned = tm_par->keep_description_of_pruned;
   lp_par->not_fixed_storage_size = tm_par->not_fixed_storage_size;
   lp_par->cut_pool_check_freq = 10;
   lp_par->load_balance_level = -1;
   lp_par->load_balance_iterations = -1;
   lp_par->load_balance_compare_candidates = HIGHEST_LOW_OBJ;
   lp_par->fractional_diving_ratio = 0.02;
   lp_par->fractional_diving_num = 0;
   lp_par->max_non_dual_feas_to_add_frac = 0.05;
   lp_par->max_non_dual_feas_to_add_min = 20;
   lp_par->max_non_dual_feas_to_add_max = 200;
   lp_par->max_not_fixable_to_add_frac = 0.1;
   lp_par->max_not_fixable_to_add_min = 100;
   lp_par->max_not_fixable_to_add_max = 500;
   lp_par->mat_col_compress_num = 50;
   lp_par->mat_col_compress_ratio = .05;
   /*
    * changed row compression so that poor cuts are fully purged and are not
    * merely lying around with very high rhs values -- asm4
    */
   lp_par->mat_row_compress_num = 0;
   lp_par->mat_row_compress_ratio = .00001;
   lp_par->tailoff_gap_backsteps = 2;
   lp_par->tailoff_gap_frac = .99;
   lp_par->tailoff_obj_backsteps = 3;
   lp_par->tailoff_obj_frac = .75;
   lp_par->tailoff_absolute = 0.0001;
   lp_par->tailoff_max_no_iterative_impr_iters_root = 3;
   lp_par->ineff_cnt_to_delete = 0;
   lp_par->eff_cnt_before_cutpool = 3;
   lp_par->ineffective_constraints = BASIC_SLACKS_ARE_INEFFECTIVE;
   lp_par->base_constraints_always_effective = TRUE;
   lp_par->branch_on_cuts = FALSE;
   lp_par->discard_slack_cuts = DISCARD_SLACKS_BEFORE_NEW_ITERATION;
   lp_par->first_lp.first_cut_time_out = 0;
   lp_par->first_lp.all_cuts_time_out = 0;
   lp_par->later_lp.first_cut_time_out = 5;
   lp_par->later_lp.first_cut_time_out = 0;
   lp_par->later_lp.all_cuts_time_out = 1;
   lp_par->later_lp.all_cuts_time_out = 0;
   lp_par->max_cut_num_per_iter = 50;
   lp_par->max_cut_num_per_iter_root = 500;
   lp_par->min_root_cut_rounds = 100;
   lp_par->tried_long_cuts = FALSE;
   lp_par->max_cut_length = 100;
   lp_par->do_reduced_cost_fixing = TRUE;
   lp_par->gap_as_ub_frac = .1;
   lp_par->gap_as_last_gap_frac = .7;
   lp_par->do_logical_fixing = 1;
   lp_par->fixed_to_ub_before_logical_fixing = 1;
   lp_par->fixed_to_ub_frac_before_logical_fixing = .01;

   lp_par->cgl.generate_cgl_cuts = TRUE;
   lp_par->cgl.generate_cgl_gomory_cuts = GENERATE_DEFAULT;
   lp_par->cgl.generate_cgl_redsplit_cuts = DO_NOT_GENERATE;
   lp_par->cgl.generate_cgl_knapsack_cuts = GENERATE_DEFAULT;
   lp_par->cgl.generate_cgl_oddhole_cuts = DO_NOT_GENERATE;
   lp_par->cgl.generate_cgl_clique_cuts = GENERATE_DEFAULT;
   lp_par->cgl.generate_cgl_probing_cuts = GENERATE_DEFAULT;
   lp_par->cgl.generate_cgl_mir_cuts = DO_NOT_GENERATE;
   lp_par->cgl.generate_cgl_twomir_cuts = GENERATE_ONLY_IN_ROOT;
   lp_par->cgl.generate_cgl_flowcover_cuts = GENERATE_DEFAULT;
   lp_par->cgl.generate_cgl_rounding_cuts = DO_NOT_GENERATE;
   lp_par->cgl.generate_cgl_lift_and_project_cuts = DO_NOT_GENERATE;
   lp_par->cgl.generate_cgl_landp_cuts = DO_NOT_GENERATE;

   lp_par->cgl.probing_is_expensive = FALSE;
   lp_par->cgl.probing_root_max_look = 100;

   lp_par->cgl.gomory_max_depth = 500;
   lp_par->cgl.probing_max_depth = 100;
   lp_par->cgl.flowcover_max_depth = 50;
   lp_par->cgl.twomir_max_depth = 50;
   lp_par->cgl.clique_max_depth = 50;
   lp_par->cgl.oddhole_max_depth = 50;
   lp_par->cgl.knapsack_max_depth = 50;

   lp_par->cgl.generate_cgl_gomory_cuts_freq =
      lp_par->cgl.generate_cgl_redsplit_cuts_freq =
      lp_par->cgl.generate_cgl_knapsack_cuts_freq =
      lp_par->cgl.generate_cgl_oddhole_cuts_freq =
      lp_par->cgl.generate_cgl_clique_cuts_freq =
      lp_par->cgl.generate_cgl_probing_cuts_freq =
      lp_par->cgl.generate_cgl_mir_cuts_freq =
      lp_par->cgl.generate_cgl_twomir_cuts_freq =
      lp_par->cgl.generate_cgl_flowcover_cuts_freq =
      lp_par->cgl.generate_cgl_rounding_cuts_freq =
      lp_par->cgl.generate_cgl_lift_and_project_cuts_freq =
      lp_par->cgl.generate_cgl_landp_cuts_freq = 5;

   lp_par->cgl.gomory_generated_in_root = FALSE;
   lp_par->cgl.redsplit_generated_in_root = FALSE;
   lp_par->cgl.knapsack_generated_in_root = FALSE;
   lp_par->cgl.oddhole_generated_in_root = FALSE;
   lp_par->cgl.probing_generated_in_root = FALSE;
   lp_par->cgl.mir_generated_in_root = FALSE;
   lp_par->cgl.twomir_generated_in_root = FALSE;
   lp_par->cgl.clique_generated_in_root = FALSE;
   lp_par->cgl.flowcover_generated_in_root = FALSE;
   lp_par->cgl.rounding_generated_in_root = FALSE;
   lp_par->cgl.lift_and_project_generated_in_root = FALSE;
   lp_par->cgl.landp_generated_in_root = FALSE;

   lp_par->cgl.use_chain_strategy = TRUE;
   lp_par->cgl.chain_status = CGL_CHAIN_START;
   lp_par->cgl.max_chain_backtrack = 1;
   lp_par->cgl.max_chain_trial_num = 10;
   lp_par->cgl.chain_trial_freq = 10;
   lp_par->cgl.chain_weighted_gap = 9.333e-6;

   lp_par->multi_criteria = FALSE;
   lp_par->mc_find_supported_solutions = FALSE;
   lp_par->mc_add_optimality_cuts = TRUE;
   lp_par->mc_gamma = 1;       /* Determines the weight on objective 1 */
   lp_par->mc_tau   = 0;       /* Determines the weight on objective 2 */
   lp_par->mc_rho   = 0.00001; /* For augmented Chebyshev norm */

#ifdef __OSI_GLPK__
   lp_par->max_presolve_iter = -1;
#else
   lp_par->max_presolve_iter = 40;
#endif

   lp_par->is_feasible_default = TEST_INTEGRALITY;
   lp_par->send_feasible_solution_default = SEND_NONZEROS;
   lp_par->display_solution_default = DISP_NOTHING;
   lp_par->shall_we_branch_default = USER__BRANCH_IF_TAILOFF;
   lp_par->select_candidates_default = USER__CLOSE_TO_HALF;
   lp_par->strong_branching_cand_num_max = 20;
   lp_par->strong_branching_cand_num_min = 5;
   lp_par->strong_branching_red_ratio = 1;
   lp_par->strong_branching_high_low_weight = 0.8; // alpha*min + (1-alpha)*max
   lp_par->user_set_strong_branching_cand_num = FALSE;
   lp_par->user_set_max_presolve_iter = FALSE;
   /*
    * strong branching is carried out for candidate variables even when pseudo
    * costs are reliably known when the depth is less than this number
   */

   lp_par->strong_br_min_level = 4;
   lp_par->strong_br_all_candidates_level = 6;
   lp_par->use_hot_starts = TRUE;
   lp_par->should_use_rel_br = FALSE;
#ifdef COMPILE_IN_LP
   lp_par->should_use_rel_br = TRUE;
#endif
#ifdef _OPENMP
   lp_par->should_use_rel_br = FALSE;
#endif
   lp_par->rel_br_override_default = TRUE;
   lp_par->rel_br_override_max_solves = 200;
   lp_par->rel_br_chain_backtrack = 5;
   lp_par->rel_br_min_imp = 0.0133;
   lp_par->rel_br_max_imp = 0.30;

   lp_par->rel_br_threshold = 8;
   lp_par->rel_br_max_solves = 20;      /* stop after these many LP-solve calls
                                           regardless of improvement */
   lp_par->rel_br_cand_threshold = 10;  /* stop doing LP-solve if last
                                           10 LP-solves didnt help. */
   lp_par->compare_candidates_default = HIGH_LOW_COMBINATION;
   lp_par->select_child_default = PREFER_LOWER_OBJ_VALUE;
   lp_par->pack_lp_solution_default = SEND_NONZEROS;
   lp_par->sensitivity_analysis = FALSE;

   /* feasibility pump */
   lp_par->fp_enabled          = SYM_FEAS_PUMP_DEFAULT;
   lp_par->fp_max_cycles       = 100;
   lp_par->fp_time_limit       = 50;
   lp_par->fp_display_interval = 10;
   lp_par->fp_poor_sol_lim_fac = 10;
   lp_par->fp_flip_fraction    = 0.1;
   lp_par->fp_frequency        = 10;
   lp_par->fp_max_initial_time = 100;
   lp_par->fp_min_gap          = 0.5;                   /* 1% gap */

   /************************** cut_gen defaults *****************************/
   cg_par->verbosity = 0;
   cg_par->do_findcuts = TRUE;

   /************************** cutpool defaults ******************************/
   cp_par->verbosity = 0;
   cp_par->warm_start = FALSE;
   cp_par->logging = FALSE;
   cp_par->block_size = 5000;
   cp_par->max_size = 2000000;
   cp_par->max_number_of_cuts = 10000;
   cp_par->cuts_to_check = 1000;
   cp_par->delete_which = DELETE_BY_QUALITY;
   cp_par->touches_until_deletion = 10;
   cp_par->min_to_delete = 1000;
   cp_par->check_which = CHECK_ALL_CUTS;

   /********************** draw_graph defaults  ******************************/
   strcpy(dg_par->source_path, ".");
   dg_par->echo_commands = FALSE;
   dg_par->canvas_width = 1000;
   dg_par->canvas_height = 700;
   dg_par->viewable_width = 600;
   dg_par->viewable_height = 400;
   dg_par->disp_nodelabels = 1;
   dg_par->disp_nodeweights = 1;
   dg_par->disp_edgeweights = 1;
   dg_par->node_dash[0] = 0;
   dg_par->edge_dash[0] = 0;
   dg_par->node_radius = 8;
   dg_par->interactive_mode = 1;
   dg_par->mouse_tracking = 1;
   dg_par->scale_factor = 1;
   strcpy(dg_par->nodelabel_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");
   strcpy(dg_par->nodeweight_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");
   strcpy(dg_par->edgeweight_font,
	  "-adobe-helvetica-bold-r-normal--11-80-*-*-*-*-*-*");

   /********************* preprocessor defaults ******************************/
   prep_par->level = 5;
   prep_par->dive_level = 5;
   prep_par->impl_dive_level = 0;
   prep_par->impl_limit = 50;
   prep_par->do_probe = 1;
   prep_par->verbosity = 1;
   prep_par->reduce_mip = 1;
   prep_par->probe_verbosity = 0;
   prep_par->probe_level = 1;
   prep_par->display_stats = 0;
   prep_par->iteration_limit = 10;
   prep_par->etol = tm_par->granularity;
   prep_par->do_single_row_rlx = 0;
   prep_par->single_row_rlx_ratio = 0.1;
   prep_par->max_sr_cnt = 5;
   prep_par->do_aggregate_row_rlx = 0;
   prep_par->max_aggr_row_cnt = 0;
   prep_par->max_aggr_row_ratio = 0.1;
   prep_par->keep_row_ordered = 1;
   prep_par->keep_track = 0;
   prep_par->time_limit = 50;
   prep_par->write_mps = 0;
   prep_par->write_lp = 0;

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_parse_command_line(sym_environment *env, int argc, char **argv)
{
   int termcode = 0;

   CALL_WRAPPER_FUNCTION( readparams_u(env, argc, argv) );

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_user_data(sym_environment *env, void *user)
{
   if (user == NULL){
      return(ERROR__USER);
   }

   env->user = user;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_user_data(sym_environment *env, void **user)
{
   if (env->user == NULL){
      return(ERROR__USER);
   }

   *user = env->user;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_read_mps(sym_environment *env, char *infile)
{

  strncpy(env->par.infile, infile, MAX_FILE_NAME_LENGTH);
  strcpy(env->par.datafile, "");
  env->par.file_type = MPS_FORMAT;
  return(sym_load_problem(env));
}

/*===========================================================================*/
/*===========================================================================*/

int sym_read_lp(sym_environment *env, char *infile)
{

  strncpy(env->par.infile, infile, MAX_FILE_NAME_LENGTH);
  strcpy(env->par.datafile, "");
  env->par.file_type = LP_FORMAT;
  return(sym_load_problem(env));
}

/*===========================================================================*/
/*===========================================================================*/

int sym_read_gmpl(sym_environment *env, char *modelfile, char *datafile)
{
  strncpy(env->par.infile, modelfile, MAX_FILE_NAME_LENGTH);
  strncpy(env->par.datafile, datafile, MAX_FILE_NAME_LENGTH);
  env->par.file_type = GMPL_FORMAT;
  return(sym_load_problem(env));
}

/*===========================================================================*/
/*===========================================================================*/

int sym_write_mps(sym_environment *env, char *infile)
{
   write_mip_desc_mps(env->mip, infile);
   return 0;
}


/*===========================================================================*/
/*===========================================================================*/

int sym_write_lp(sym_environment *env, char *infile)
{
   write_mip_desc_lp(env->mip, infile);
   return 0;

}

/*===========================================================================*/
/*===========================================================================*/

int sym_load_problem(sym_environment *env)
{
   double t = 0;
   int termcode = 0;

   /*------------------------------------------------------------------------*\
    *                         start reading in problem
   \*------------------------------------------------------------------------*/

   (void) used_time(&t);

   /* Get the problem data */
   CALL_WRAPPER_FUNCTION( io_u(env) );

   /* Start up the graphics window*/
#if !defined(_MSC_VER) && !defined (__MNO_CYGWIN)
   CALL_WRAPPER_FUNCTION( init_draw_graph_u(env) );
#endif

   /*------------------------------------------------------------------------*\
    * Have the user generate the base and root description
   \*------------------------------------------------------------------------*/

   CALL_WRAPPER_FUNCTION( initialize_root_node_u(env) );

   env->comp_times.readtime = used_time(&t);

   env->termcode = TM_NO_SOLUTION;
   env->mip->is_modified = TRUE;

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_find_initial_bounds(sym_environment *env)
{
   double total_time = 0;
   int termcode = 0;

   /* Finds the upper and lower bounds for the problem */
   CALL_WRAPPER_FUNCTION( start_heurs_u(env) );

   if (!env->par.do_branch_and_cut){
      printf("\n****************************************************\n");
      printf(  "* Heuristics Finished!!!!!!!                       *\n");
      printf(  "* Now displaying stats and best solution....       *\n");
      printf(  "****************************************************\n\n");
      total_time += env->comp_times.ub_overhead + env->comp_times.ub_heurtime;
      total_time += env->comp_times.lb_overhead + env->comp_times.lb_heurtime;
#if !defined(_MSC_VER) && !defined (__MNO_CYGWIN) /* FIXME: CPU timing doesn't work in Windows */
      printf( "  Problem IO     %.3f\n", env->comp_times.readtime);
      printf( "  Overhead: UB   %.3f\n", env->comp_times.ub_overhead);
      printf( "            LB   %.3f\n", env->comp_times.lb_overhead);
      printf( "  Runtime:  UB   %.3f\n", env->comp_times.ub_heurtime);
      printf( "            LB   %.3f\n", env->comp_times.lb_heurtime);
      printf( "  Total User Time    %.3f\n", total_time);
#endif
      if (env->has_ub){
	 if (env->mip->obj_sense == SYM_MAXIMIZE){
	    printf( "Lower Bound: %.3f\n", -env->ub + env->mip->obj_offset);
	 }else{
	    printf( "Upper Bound: %.3f\n", env->ub + env->mip->obj_offset);
	 }
      }
      CALL_WRAPPER_FUNCTION( display_solution_u(env, 0) );
      if (env->par.tm_par.lp_machs)
	 FREE(env->par.tm_par.lp_machs[0]);
      FREE(env->par.tm_par.lp_machs);
   }

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_solve(sym_environment *env)
{

   int s_bufid, r_bufid, bytes, msgtag = 0, sender, termcode = 0, temp, i;
   int lp_data_sent = 0, cg_data_sent = 0, cp_data_sent = 0;
#ifndef COMPILE_IN_TM
   char repricing, node_type;
#else
   tm_prob *tm;
#endif
   double start_time, lb;
   struct timeval timeout = {10, 0};
   double total_time = 0;

   int granularity = 0;

   node_desc *rootdesc = env->rootdesc;
   base_desc *base = env->base;

   start_time = wall_clock(NULL);

   double *tmp_sol;
   lp_sol *best_sol = &(env->best_sol);

   if(best_sol->has_sol && env->mip->is_modified){
      FREE(best_sol->xind);
      FREE(best_sol->xval);
      best_sol->has_sol = FALSE;
   }

   env->mip->is_modified = FALSE;

#ifndef USE_SYM_APPLICATION

   /* we send environment in just because we may need to
      update rootdesc and so...*/

   if (!env->par.multi_criteria){
       termcode = sym_presolve(env);
   }else{
       env->par.prep_par.level = 0;
   }

   if(termcode == PREP_INFEAS || termcode == PREP_UNBOUNDED ||
      termcode == PREP_SOLVED || termcode == PREP_NUMERIC_ERROR ||
      termcode == PREP_OTHER_ERROR){

      env->mip = env->orig_mip;
      env->orig_mip = 0;

      if(termcode == PREP_INFEAS){
	 return(env->termcode = PREP_NO_SOLUTION);
      }else if(termcode == PREP_UNBOUNDED){
	 return(env->termcode = PREP_UNBOUNDED);
      }else if(termcode == PREP_SOLVED){
	 best_sol->has_sol = TRUE;
	 best_sol->xind = (int *) malloc(ISIZE *
					     env->prep_mip->fixed_n);
	 best_sol->xval = (double *) malloc(DSIZE *
						env->prep_mip->fixed_n);

	 best_sol->xlength = env->prep_mip->fixed_n;
	 memcpy(best_sol->xind, env->prep_mip->fixed_ind, ISIZE *
		env->prep_mip->fixed_n);
	 memcpy(best_sol->xval, env->prep_mip->fixed_val, ISIZE *
		env->prep_mip->fixed_n);

	 return(env->termcode = PREP_OPTIMAL_SOLUTION_FOUND);
      }else if(termcode == PREP_NUMERIC_ERROR){
	 return(env->termcode = PREP_ERROR);
      }
   }

   if(termcode == PREP_OTHER_ERROR || env->par.prep_par.level <= 0){
      if(env->prep_mip){
	 free_mip_desc(env->prep_mip);
	 FREE(env->prep_mip);
	 env->prep_mip = 0;
      }
   }

#endif

   if (env->par.verbosity >= -1){
      printf("Solving...\n\n");
   }

#ifndef COMPILE_IN_TM
   /*------------------------------------------------------------------------*\
    * Start the tree manager and send the parameters
   \*------------------------------------------------------------------------*/

   if (env->par.tm_machine_set){
      spawn(env->par.tm_exe, (char **)NULL, env->par.tm_debug | TaskHost,
	    env->par.tm_machine, 1, &env->tm_tid);
   }else{
      spawn(env->par.tm_exe, (char **)NULL, env->par.tm_debug, (char *)NULL, 1,
	    &env->tm_tid);
   }
   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&env->par.tm_par), sizeof(tm_params));
   send_char_array(&env->has_ub, 1);
   if (env->has_ub)
      send_dbl_array(&env->ub, 1);
   send_char_array(&env->has_ub_estimate, 1);
   if (env->has_ub_estimate)
      send_dbl_array(&env->ub_estimate, 1);
   if (env->par.tm_par.lp_mach_num)
      send_char_array(env->par.tm_par.lp_machs[0],
		      env->par.tm_par.lp_mach_num*MACH_NAME_LENGTH);
   if (env->par.tm_par.cg_mach_num)
      send_char_array(env->par.tm_par.cg_machs[0],
		      env->par.tm_par.cg_mach_num*MACH_NAME_LENGTH);
   if (env->par.tm_par.cp_mach_num)
      send_char_array(env->par.tm_par.cp_machs[0],
		      env->par.tm_par.cp_mach_num*MACH_NAME_LENGTH);
   send_int_array(&base->varnum, 1);
   send_int_array(&base->cutnum, 1);
#ifdef TRACE_PATH
   {
      int feas_sol_size;
      int *feas_sol;

#ifdef USE_SYM_APPLICATION
      if (user_send_feas_sol(env->user, &feas_sol_size, &feas_sol)==USER_NO_PP){
	 send_int_array(&feas_sol_size, 1);
	 if (feas_sol_size){
	    send_int_array(feas_sol, feas_sol_size);
	 }
      }
#endif
   }
#endif
   send_msg(env->tm_tid, TM_DATA);

   /*------------------------------------------------------------------------*\
    * Send out the root node
   \*------------------------------------------------------------------------*/

   if (!env->par.warm_start){
      repricing = FALSE;
      node_type = ROOT_NODE;

      s_bufid = init_send(DataInPlace);
      send_char_array(&repricing, 1);
      send_char_array(&node_type, 1);
      send_dbl_array(&env->lb, 1);
      send_int_array(&rootdesc->nf_status, 1);
      pack_array_desc(&rootdesc->uind);
      if (rootdesc->nf_status == NF_CHECK_AFTER_LAST ||
	  rootdesc->nf_status == NF_CHECK_UNTIL_LAST)
	 pack_array_desc(&rootdesc->not_fixed);
      pack_array_desc(&rootdesc->cutind);
      pack_basis(&rootdesc->basis, TRUE);
      send_int_array(&rootdesc->desc_size, 1);
      if (rootdesc->desc_size)
	 send_char_array(rootdesc->desc, rootdesc->desc_size);
      if (rootdesc->cutind.size > 0){ /* Hey, we have cuts! Pack them, too. */
	 /* Pack their number again, so we can call unpack_cut_set in TM */
	 int i;
	 send_int_array(&rootdesc->cutind.size, 1);
	 for (i = 0; i < rootdesc->cutind.size; i++)
	    pack_cut(rootdesc->cuts[i]);
      }
      send_msg(env->tm_tid, TM_ROOT_DESCRIPTION);
      freebuf(s_bufid);
   }
#else

   /*------------------------------------------------------------------------*\
    * Create the treemanager and copy the problem data
   \*------------------------------------------------------------------------*/

   /*
    * set granularity.
    * TODO: move this to preprocessor when it becomes available
    * TODO: put a new flag that checks if user wants to override this
    * TODO: find if granularity could be 0.1 or 0.2 or ... instead of just
    *       1.0, 2.0, ...
    */
   if (env->mip && env->mip->obj && env->par.tm_par.granularity<=0.000001
       && !env->par.multi_criteria) {
      for (int i=0;i<env->mip->n;i++) {
         double coeff = env->mip->obj[i];
         if (fabs(coeff)>0.000001) {
            if (env->mip->is_int[i]) {
               if (fabs(floor(coeff+0.5)-coeff)<0.000001) {
                  granularity = sym_gcd(granularity,(int)floor(coeff+0.5));
               } else {
                  granularity = 0;
                  break;
               }
            } else if (env->mip->ub[i]-env->mip->lb[i]>0.000001) {
               granularity=0;
               break;
            }
         } // else do nothing
      }
      /*
       * if granularity >= 1, set it at granularity - epsilon, otherwise set at
       * epsilon
       */
      env->par.tm_par.granularity = env->par.lp_par.granularity =
         fabs((double)granularity-0.000001);
   }
   PRINT(env->par.verbosity, 0, ("granularity set at %f\n",
            env->par.tm_par.granularity));

   if (env->par.tm_par.node_selection_rule == BEST_FIRST_SEARCH){
      env->par.tm_par.node_selection_rule = LOWEST_LP_FIRST;
   }

   env->tm = tm = (tm_prob *) calloc(1, sizeof(tm_prob));

   tm->par = env->par.tm_par;

   if ((tm->has_ub = env->has_ub))
	tm->ub = env->ub;
   if ((tm->has_ub_estimate = env->has_ub_estimate))
      tm->ub_estimate = env->ub_estimate;
   tm->lb = env->lb;

   if(env->obj_offset){
      env->mip->obj_offset += env->obj_offset;
   }

   tm->obj_offset = env->mip->obj_offset;
   tm->obj_sense = env->mip->obj_sense;
   tm->master = env->my_tid;

#ifdef COMPILE_IN_LP
   CALL_WRAPPER_FUNCTION( send_lp_data_u(env, 0) );
#ifdef _OPENMP
   lp_data_sent = env->par.tm_par.max_active_nodes;
#else
   lp_data_sent = 1;
#endif
#ifdef COMPILE_IN_CG
   CALL_WRAPPER_FUNCTION( send_cg_data_u(env, 0) );
#ifdef _OPENMP
   cg_data_sent = env->par.tm_par.max_active_nodes;
#else
   cg_data_sent = 1;
#endif
#endif
#endif
#ifdef COMPILE_IN_CP
   if (env->cp && env->par.use_permanent_cut_pools){
      tm->cpp = env->cp;
   }else{
      CALL_WRAPPER_FUNCTION( send_cp_data_u(env, 0) );
   }
#ifdef _OPENMP
   cp_data_sent = env->par.tm_par.max_cp_num;
#else
   cp_data_sent = 1;
#endif
#endif

   // Check stored solution to see if it is still feasible

   if (best_sol->has_sol){
      tmp_sol = (double *) calloc(env->mip->n, DSIZE);
      for (i = 0; i < best_sol->xlength; i++){
	 tmp_sol[best_sol->xind[i]] = best_sol->xval[i];
      }
      sym_set_col_solution(env, tmp_sol);
   }

   //   memset(&(env->best_sol), 0, sizeof(lp_sol));

   if (env->warm_start && env->par.tm_par.warm_start){
      /* Load warm start info */
      tm->rootnode = env->warm_start->rootnode;
      tm->cuts = env->warm_start->cuts;
      tm->cut_num = env->warm_start->cut_num;
      tm->allocated_cut_num = env->warm_start->allocated_cut_num;
      tm->stat = env->warm_start->stat;
      tm->comp_times = env->warm_start->comp_times;
      tm->lb = env->warm_start->lb;
      if (env->warm_start->has_ub){
	 if (env->warm_start->ub < tm->ub || !tm->has_ub){
	    tm->ub = env->warm_start->ub;
	 }
	 tm->has_ub = TRUE;
      }
      if (best_sol->objval > env->warm_start->best_sol.objval){
	 FREE(best_sol->xind);
	 FREE(best_sol->xval);
	 env->best_sol = env->warm_start->best_sol;
      }
      tm->phase = env->warm_start->phase;
   }else if (env->warm_start){
      /* Otherwise, free what was saved */
      free_subtree(env->warm_start->rootnode);
      if(env->warm_start->best_sol.xlength){
	 FREE(env->warm_start->best_sol.xind);
	 FREE(env->warm_start->best_sol.xval);
      }
      if (env->warm_start->cuts){
	 for (i = env->warm_start->cut_num - 1; i >= 0; i--)
	    if (env->warm_start->cuts[i]){
	       FREE(env->warm_start->cuts[i]->coef);
	       FREE(env->warm_start->cuts[i]);
	    }
	 FREE(env->warm_start->cuts);
      }
   }
   /* Now the tree manager owns everything */
   FREE(env->warm_start);

   if ((termcode = tm_initialize(tm , base, rootdesc)) < 0){
      tm_close(tm, termcode);

      if (env->par.do_draw_graph){
	 s_bufid = init_send(DataInPlace);
	 send_msg(env->dg_tid, CTOI_YOU_CAN_DIE);
	 freebuf(s_bufid);
      }

      if (env->par.tm_par.lp_machs)
	 FREE(env->par.tm_par.lp_machs[0]);
      FREE(env->par.tm_par.lp_machs);
      if (env->par.tm_par.cg_machs)
	 FREE(env->par.tm_par.cg_machs[0]);
      FREE(env->par.tm_par.cg_machs);
      if (env->par.tm_par.cp_machs)
	 FREE(env->par.tm_par.cp_machs[0]);
      FREE(env->par.tm_par.cp_machs);

      free_tm(tm);

      env->termcode = termcode;

      return(termcode);
   }


#ifdef TRACE_PATH
   {
      int feas_sol_size;
      int *feas_sol;
#ifdef USE_SYM_APPLICATION
      if (user_send_feas_sol(env->user,&feas_sol_size,&feas_sol)==USER_NO_PP){
	 tm->feas_sol_size = feas_sol_size;
	 tm->feas_sol = (int *) calloc (tm->feas_sol_size, sizeof(int));
	 memcpy((char *)tm->feas_sol, (char *)feas_sol, feas_sol_size * ISIZE);
      }
#endif
   }
#endif
#endif

   /*------------------------------------------------------------------------*\
    * Wait for messages
   \*------------------------------------------------------------------------*/

#ifdef COMPILE_IN_TM
   while (!(lp_data_sent == env->par.tm_par.max_active_nodes) ||
	  !(cg_data_sent == env->par.tm_par.max_active_nodes) ||
	  !(cp_data_sent == env->par.tm_par.max_cp_num)){
#else
   do{
            /* } unconfuse vi */
#endif
      r_bufid = treceive_msg(ANYONE, ANYTHING, &timeout);
      if (r_bufid == 0){
#ifndef COMPILE_IN_TM
	 if (pstat(env->tm_tid) != PROCESS_OK){
	    printf("\nThe treemanager has died :-(\n\n");
#else
	 if (!processes_alive(env->tm)){
            /* } unconfuse vi */
#endif
	    termcode = msgtag = SOMETHING_DIED;
	    break;
	 }else{
	    continue;
	 }
      }
      bufinfo(r_bufid, &bytes, &msgtag, &sender);

      switch (msgtag){
       case FEASIBLE_SOLUTION_NONZEROS:
       case FEASIBLE_SOLUTION_USER:
	 CALL_WRAPPER_FUNCTION( receive_feasible_solution_u(env, msgtag) );
	 if (env->par.verbosity >= -1){
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
	    CALL_WRAPPER_FUNCTION( display_solution_u(env,
						env->tm->opt_thread_num) );
#else
	    CALL_WRAPPER_FUNCTION( display_solution_u(env, 0) );
#endif
	 }
	 break;

       case REQUEST_FOR_LP_DATA:
	 /* An LP process has been started and asks for all necessary data */
	 CALL_WRAPPER_FUNCTION( send_lp_data_u(env, sender) );
	 lp_data_sent++;
	 break;

       case REQUEST_FOR_CG_DATA:
	 /* A CG process has been started and asks for all necessary data */
	 CALL_WRAPPER_FUNCTION( send_cg_data_u(env, sender) );
	 cg_data_sent++;
	 break;

       case REQUEST_FOR_CP_DATA:
	 /* A CP process has been started and asks for all necessary data */
	 CALL_WRAPPER_FUNCTION( send_cp_data_u(env, sender) );
	 cp_data_sent++;
	 break;

       case TM_FIRST_PHASE_FINISHED:
	 receive_char_array((char *)(&env->comp_times.bc_time),
			     sizeof(node_times));
	 receive_dbl_array(&lb, 1);
	 if (lb > env->lb) env->lb = lb;
	 receive_char_array((char *)&env->warm_start->stat,
			    sizeof(problem_stat));
	 printf( "\n");
	 printf( "****************************************************\n");
	 printf( "* Branch and Cut First Phase Finished!!!!          *\n");
	 printf( "* Now displaying stats and best solution...        *\n");
	 printf( "****************************************************\n\n");

	 print_statistics(&(env->comp_times.bc_time), &(env->warm_start->stat),
                          NULL,
			  env->ub, env->lb, 0, start_time, wall_clock(NULL),
			  env->mip->obj_offset, env->mip->obj_sense,
			  env->has_ub, NULL);
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
	 CALL_WRAPPER_FUNCTION( display_solution_u(env,
						   env->tm->opt_thread_num) );
#else
	 CALL_WRAPPER_FUNCTION( display_solution_u(env, 0) );
#endif
	 break;

       case SOMETHING_DIED:
       case TM_TIME_LIMIT_EXCEEDED:
       case TM_SIGNAL_CAUGHT:
       case TM_NODE_LIMIT_EXCEEDED:
       case TM_TARGET_GAP_ACHIEVED:
       case TM_FOUND_FIRST_FEASIBLE:
       case TM_OPTIMAL_SOLUTION_FOUND:
       case TM_ERROR__NO_BRANCHING_CANDIDATE:
       case TM_ERROR__ILLEGAL_RETURN_CODE:
       case TM_ERROR__NUMERICAL_INSTABILITY:
       case TM_ERROR__COMM_ERROR:
       case TM_ERROR__USER:
	 receive_char_array((char *)(&env->comp_times.bc_time),
			    sizeof(node_times));
	 receive_dbl_array(&lb, 1);
	 if (lb > env->lb) env->lb = lb;
	 receive_char_array((char *)&env->warm_start->stat,
			    sizeof(problem_stat));
	 break;

       default:
	 CALL_WRAPPER_FUNCTION( process_own_messages_u(env, msgtag) );
	 break;
      }
      freebuf(r_bufid);

#ifndef COMPILE_IN_TM
   }while (msgtag != TM_OPTIMAL_SOLUTION_FOUND && msgtag != SOMETHING_DIED &&
	   msgtag != TM_TIME_LIMIT_EXCEEDED &&
	   msgtag != TM_SIGNAL_CAUGHT &&
	   msgtag != TM_NODE_LIMIT_EXCEEDED &&
	   msgtag != TM_TARGET_GAP_ACHIEVED &&
	   msgtag != TM_FOUND_FIRST_FEASIBLE &&
	   msgtag != TM_ERROR__NO_BRANCHING_CANDIDATE &&
	   msgtag != TM_ERROR__ILLEGAL_RETURN_CODE &&
	   msgtag != TM_ERROR__NUMERICAL_INSTABLITY &&
	   msgtag != TM_ERROR__COMM_ERROR &&
	   msgtag != TM_ERROR__USER);

   termcode = msgtag;
#else
      /* unconfuse vi { */
   }

   /*------------------------------------------------------------------------*\
    * Solve the problem and receive solutions
   \*------------------------------------------------------------------------*/
#ifdef COMPILE_IN_LP
   sp_initialize(tm);
#endif

   tm->start_time += start_time;

   termcode = solve(tm);
   tm_close(tm, termcode);

   /* Save the warm start info */
   env->warm_start = (warm_start_desc *) calloc (1, sizeof(warm_start_desc));
   env->warm_start->rootnode = tm->rootnode;
   env->warm_start->cuts = env->tm->cuts;
   env->warm_start->cut_num = env->tm->cut_num;
   env->warm_start->allocated_cut_num = env->tm->allocated_cut_num;
   env->warm_start->stat = tm->stat;
   env->warm_start->phase = tm->phase;
   env->warm_start->lb = tm->lb;
   if ((env->warm_start->has_ub = tm->has_ub)!=0){
      env->warm_start->ub = tm->ub;
   }
   env->par.tm_par.warm_start = FALSE;

#ifdef COMPILE_IN_LP
   int thread_num;
   thread_num = env->tm->opt_thread_num;
   if (env->tm->lpp[thread_num]){
      env->par.lp_par.cgl = env->tm->lpp[thread_num]->par.cgl;
      if (env->tm->lpp[thread_num]->best_sol.has_sol){
	 FREE(best_sol->xind);
	 FREE(best_sol->xval);
	 env->best_sol =
	    env->tm->lpp[thread_num]->best_sol;
      }else {
	 env->tm->lpp[thread_num]->best_sol = env->best_sol;
      }
   }
#else
   if (env->tm->best_sol.has_sol){
     FREE(best_sol->xind);
     FREE(best_sol->xval);
     env->best_sol = env->tm->best_sol;
   }
#endif

   if (best_sol->has_sol) {
      memcpy(&env->warm_start->best_sol, &env->best_sol, sizeof(lp_sol) *1);
      env->warm_start->best_sol.xind = 0;
      env->warm_start->best_sol.xval = 0;
      if(best_sol->xlength){
	 env->warm_start->best_sol.xind = (int *) malloc(ISIZE * best_sol->xlength);
	 env->warm_start->best_sol.xval = (double *) malloc(DSIZE *
							    best_sol->xlength);
	 memcpy(env->warm_start->best_sol.xind,
		best_sol->xind, ISIZE * best_sol->xlength);
	 memcpy(env->warm_start->best_sol.xval,
		best_sol->xval, DSIZE * best_sol->xlength);
      }
   }

   tm->rootnode = NULL;
   tm->cuts = NULL;
   tm->cut_num = tm->allocated_cut_num = 0;
#ifdef COMPILE_IN_CP
   if (env->cp && env->par.use_permanent_cut_pools){
      tm->cpp = NULL;
   }
#endif

#if !defined(COMPILE_IN_LP) && 0
   /* This is not needed anymore */
   if (termcode != SOMETHING_DIED){
      int old_termcode = termcode;
      do{
	 r_bufid = receive_msg(ANYONE, ANYTHING);
	 if (r_bufid == 0){
	    printf("\nError receiving solution ...\n");
	    break;
	 }
	 bufinfo(r_bufid, &bytes, &msgtag, &sender);
	 if (msgtag == FEASIBLE_SOLUTION_NONZEROS ||
	     msgtag == FEASIBLE_SOLUTION_USER){
	    CALL_WRAPPER_FUNCTION( receive_feasible_solution_u(env, msgtag) );
	 }
      }while (msgtag != FEASIBLE_SOLUTION_NONZEROS &&
	      msgtag != FEASIBLE_SOLUTION_USER);
      termcode = old_termcode;
   }
#endif

   /* FIXME: Set the correct termcode. This can't be done in the treemanager
      because it doesn't know whether a solution was found. This should be
      changed. */
   if (termcode == TM_FINISHED){
      if (tm->par.find_first_feasible && best_sol->has_sol){
	 termcode = TM_FOUND_FIRST_FEASIBLE;
      }else if (best_sol->has_sol){
	 termcode = TM_OPTIMAL_SOLUTION_FOUND;
      }else{
	 termcode = TM_NO_SOLUTION;
      }
   }
#if 0
   /* Not sure of the reason for this */
   else if((termcode == TM_ERROR__NUMERICAL_INSTABILITY ||
	    termcode == SOMETHING_DIED) &&
	   best_sol->xlength ){
     termcode = TM_FEASIBLE_SOLUTION_FOUND;
   }
#endif

#endif

   /*------------------------------------------------------------------------*\
    * Display the the results and solution data
   \*------------------------------------------------------------------------*/

   if (env->par.verbosity >= -1 ){
      printf("\n****************************************************\n");
      if (termcode == TM_OPTIMAL_SOLUTION_FOUND){
	 printf(  "* Optimal Solution Found                           *\n");
      }else if (termcode == TM_NO_SOLUTION || termcode == TM_UNBOUNDED){
	 printf(  "* Branch and Cut Finished                          *\n");
      }else if (termcode == TM_TIME_LIMIT_EXCEEDED){
	 printf(  "* Time Limit Reached                               *\n");
      }else if (termcode == TM_NODE_LIMIT_EXCEEDED){
	 printf(  "* Node Limit Reached                               *\n");
      }else if (termcode == TM_SIGNAL_CAUGHT){
	 printf(  "* Abort Requested                                  *\n");
      }else if (termcode == TM_TARGET_GAP_ACHIEVED){
	 printf(  "* Target Gap Achieved                              *\n");
      }else if (termcode == TM_FOUND_FIRST_FEASIBLE){
	 printf(  "* Stopping After Finding First Feasible Solution   *\n");
      }else if (termcode == TM_ERROR__NO_BRANCHING_CANDIDATE ||
		termcode == TM_ERROR__ILLEGAL_RETURN_CODE ||
		termcode == TM_ERROR__NUMERICAL_INSTABILITY ||
		termcode == TM_ERROR__COMM_ERROR ||
		termcode == TM_ERROR__USER){
	 printf(  "* Terminated abnormally with error message %i      *\n",
		  termcode);
      }else{
	 printf(  "* A process has died abnormally -- halting         *\n");
      }
      if (env->par.verbosity >=0 ){
	 printf(  "* Now displaying stats and best solution found...  *\n");
      }
      printf(  "****************************************************\n\n");
      if (env->par.verbosity >=0 ){
	 total_time  = env->comp_times.readtime;
	 total_time += env->comp_times.ub_overhead + env->comp_times.ub_heurtime;
	 total_time += env->comp_times.lb_overhead + env->comp_times.lb_heurtime;

#if !defined(_MSC_VER) && defined (__MNO_CYGWIN) /* FIXME: CPU timing doesn't work in Windows */
	 printf( "====================== Misc Timing =========================\n");
	 printf( "  Problem IO        %.3f\n", env->comp_times.readtime);
#if 0
	 printf( "  UB overhead:      %.3f\n", env->comp_times.ub_overhead);
	 printf( "  UB runtime:       %.3f\n", env->comp_times.ub_heurtime);
	 printf( "  LB overhead:      %.3f\n", env->comp_times.lb_overhead);
	 printf( "  LB runtime:       %.3f\n", env->comp_times.lb_heurtime);
#endif
#endif
      }
   }

   env->termcode = termcode;

#ifdef COMPILE_IN_TM
      if (tm->lb > env->lb) env->lb = tm->lb;
      if(env->par.verbosity >=0 ) {
	 print_statistics(&(tm->comp_times), &(tm->stat),
                          &(tm->lp_stat),
                          tm->ub, env->lb,
			  total_time, start_time, wall_clock(NULL),
			  env->mip->obj_offset, env->mip->obj_sense,
			  tm->has_ub, tm->sp);
      }
      temp = termcode;
#ifdef COMPILE_IN_LP
      sp_free_sp(tm->sp);
      FREE(tm->sp);
#endif

      if(env->par.verbosity >=-1 ) {
#ifdef COMPILE_IN_LP
	 CALL_WRAPPER_FUNCTION( display_solution_u(env, env->tm->opt_thread_num) );
#else
	 CALL_WRAPPER_FUNCTION( display_solution_u(env, 0) );
#endif
      }
#else
      if(env->par.verbosity >=0 ) {
	 print_statistics(&(env->comp_times.bc_time), &(env->warm_start->stat),
                          NULL,
			  env->ub, env->lb, 0, start_time, wall_clock(NULL),
			  env->mip->obj_offset, env->mip->obj_sense,
                          env->has_ub, NULL);
	 CALL_WRAPPER_FUNCTION( display_solution_u(env, 0) );
      }
#endif
   termcode = temp;
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   if (env->tm && env->tm->lpp[env->tm->opt_thread_num]){
      env->tm->lpp[env->tm->opt_thread_num]->best_sol.xlength = 0;
      env->tm->lpp[env->tm->opt_thread_num]->best_sol.xind = NULL;
      env->tm->lpp[env->tm->opt_thread_num]->best_sol.xval = NULL;
   }
#endif


#ifndef USE_SYM_APPLICATION

   if(env->par.prep_par.level > 0){
      if(env->orig_mip){
	 env->mip = env->orig_mip;
	 env->orig_mip = 0;
      }
   }

#endif

   env->has_ub = FALSE;
   env->ub = 0.0;
   env->lb = -MAXDOUBLE;

   if (env->par.do_draw_graph){
      s_bufid = init_send(DataInPlace);
      send_msg(env->dg_tid, CTOI_YOU_CAN_DIE);
      freebuf(s_bufid);
   }

   if (env->par.tm_par.lp_machs)
      FREE(env->par.tm_par.lp_machs[0]);
   FREE(env->par.tm_par.lp_machs);
   if (env->par.tm_par.cg_machs)
      FREE(env->par.tm_par.cg_machs[0]);
   FREE(env->par.tm_par.cg_machs);
   if (env->par.tm_par.cp_machs)
      FREE(env->par.tm_par.cp_machs[0]);
   FREE(env->par.tm_par.cp_machs);

   free_tm(tm);

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/
int sym_warm_solve(sym_environment *env)
{

   int i, change_type;
   int node_limit, analyzed, depth, index, rated, level, level_rated;
   cut_data **upd_cuts;
   int loc, ws_cnum, cut_num = 0, *cut_ind = NULL, *tmp_ind = NULL;
   /* first check for the updates! */
   char *cru_vars = NULL;
   double etol = 1e-04;
   if(env->par.tm_par.keep_description_of_pruned != KEEP_IN_MEMORY){

      return(sym_solve(env));

   }else{

      if (env->warm_start){
	 env->par.tm_par.warm_start = TRUE;
      } else {
	return(sym_solve(env));
      }

      if(env->mip->change_num){
	 env->has_ub = FALSE;
	 env->ub = 0.0;
	 env->lb = -MAXDOUBLE;

	 env->warm_start->has_ub = env->best_sol.has_sol =
	    env->warm_start->best_sol.has_sol = FALSE;
	 env->warm_start->ub = env->warm_start->best_sol.objval = 0.0;
	 env->warm_start->lb = -MAXDOUBLE;
	 env->warm_start->best_sol.xlength = 0;
	 FREE(env->warm_start->best_sol.xind);
	 FREE(env->warm_start->best_sol.xval);
      }else {
	 env->has_ub = env->warm_start->has_ub;
	 env->ub = env->warm_start->ub;
	 env->lb = env->warm_start->lb;
      }

      if(env->par.multi_criteria){
	 env->has_ub = env->has_mc_ub;
	 env->ub = env->mc_ub;
      }

      for(i = 0; i < env->mip->change_num; i++){
	 change_type = env->mip->change_type[i];
	 if(change_type == RHS_CHANGED || change_type == COL_BOUNDS_CHANGED ||
	    change_type == OBJ_COEFF_CHANGED || change_type == COLS_ADDED){
	    if(change_type == OBJ_COEFF_CHANGED){
	       if(env->par.lp_par.do_reduced_cost_fixing && !env->par.multi_criteria){
		  printf("sym_warm_solve(): SYMPHONY can not resolve for the\n");
		  printf("obj coeff change when reduced cost fixing is on,");
		  printf("for now!\n");
		  return(FUNCTION_TERMINATED_ABNORMALLY);
	       }
	    } else{
	       if(env->par.lp_par.cgl.generate_cgl_cuts){
		  printf("sym_warm_solve(): SYMPHONY can not resolve for the\n");
		  printf("rhs or column bounds change when cuts exist, for now!\n");
		  return(FUNCTION_TERMINATED_ABNORMALLY);
	       }
	    }
	    if(!env->mip->cru_vars_num){
	       analyzed = env->warm_start->stat.analyzed;
	       depth = env->warm_start->stat.max_depth;
	       rated = (int)(env->par.tm_par.warm_start_node_ratio * analyzed);
	       level_rated = (int)(env->par.tm_par.warm_start_node_level_ratio * depth);
	       node_limit = env->par.tm_par.warm_start_node_limit;
	       level = env->par.tm_par.warm_start_node_level;
	       index = node_limit <= rated ? node_limit : rated ;
	       level = level <= level_rated ? level : level_rated;

	       if ((level > 0 && level < depth) || index > 0) {
		  if ( level > 0 && level < depth) {
		     env->warm_start->trim_tree = TRIM_LEVEL;
		     env->warm_start->trim_tree_level = level;
		     //cut_ws_tree_level(env, env->warm_start->rootnode, level,
		     //	    &(env->warm_start->stat), change_type);
		     env->warm_start->stat.max_depth = level;
		  } else {
		     if (index < analyzed) {
			if (!index) index = 1;
			env->warm_start->trim_tree = TRIM_INDEX;
			env->warm_start->trim_tree_index = index;
			//   cut_ws_tree_index(env, env->warm_start->rootnode, index,
			//	       &(env->warm_start->stat), change_type);
		     }
		  }
	       }
	    }else{
	       env->warm_start->trim_tree = ON_CRU_VARS;
	       cru_vars = (char *)calloc(CSIZE,env->mip->n);
	       for(i = 0; i < env->mip->cru_vars_num; i++){
		  cru_vars[env->mip->cru_vars[i]] = TRUE;
	       }
	    }

	    ws_cnum = env->warm_start->cut_num;
	    if(env->warm_start->trim_tree && ws_cnum){
	       cut_ind = (int *)malloc(ISIZE*ws_cnum);
	       memset(cut_ind, -1, ISIZE*ws_cnum);
	    }
	    env->warm_start->stat.analyzed =
	       env->warm_start->stat.created =
	       env->warm_start->stat.tree_size = 1; //for root node */

	    update_tree_bound(env, env->warm_start->rootnode, &cut_num, cut_ind, cru_vars, change_type);

	    /* FIXME!!!! feasible solutions are getting lost in a sequence of warm-solve---
	       for a temporary fix, increase ub a litte... */
	    if(env->warm_start->has_ub){
	       env->warm_start->ub += etol;
	    }

	    if(cut_num > 0){
	       upd_cuts = (cut_data **)malloc(sizeof(cut_data *)*env->warm_start->allocated_cut_num);
	       tmp_ind = (int *)malloc(ISIZE*ws_cnum);
	       for(i = 0; i < ws_cnum; i++){
		  tmp_ind[i] = i;
	       }
	       qsort_ii(cut_ind, tmp_ind, ws_cnum);

	       for(i = 0; i < cut_num; i++){
		  loc = tmp_ind[ws_cnum - cut_num + i];
		  upd_cuts[i] = env->warm_start->cuts[loc];
		  upd_cuts[i]->name = i;
		  env->warm_start->cuts[loc] = 0;
	       }
	       for (i = env->warm_start->cut_num - 1; i >= 0; i--){
		  if (env->warm_start->cuts[i]){
		     FREE(env->warm_start->cuts[i]->coef);
		  }
		  FREE(env->warm_start->cuts[i]);
	       }
	       FREE(env->warm_start->cuts);
	       env->warm_start->cuts = upd_cuts;
	       env->warm_start->cut_num = cut_num;
	    } else{
	       if(env->warm_start->trim_tree && env->warm_start->cut_num){
		  for (i = env->warm_start->cut_num - 1; i >= 0; i--){
		     if (env->warm_start->cuts[i]){
			FREE(env->warm_start->cuts[i]->coef);
		     }
		     FREE(env->warm_start->cuts[i]);
		  }
		  //	  FREE(env->warm_start->cuts);
		  //env->warm_start->cuts = 0;
		  env->warm_start->cut_num = 0;
	       }
	    }

#ifdef USE_SYM_APPLICATION
	    cut_data * cut;
	    if(change_type == COLS_ADDED || change_type == RHS_CHANGED){
	       for(i = 0; i < env->warm_start->cut_num; i++){
		  cut = env->warm_start->cuts[i];
		  user_ws_update_cuts(env->user, &(cut->size), &(cut->coef),
				      &(cut->rhs), &(cut->sense), cut->type,
				      env->mip->new_col_num,
				      change_type);
	       }
	    }
#endif
	    env->warm_start->trim_tree = DO_NOT_TRIM;
	    env->mip->change_num = 0;
	    env->mip->var_type_modified = FALSE;
	    env->mip->new_col_num = 0;
	    if(env->mip->cru_vars_num){
	       FREE(env->mip->cru_vars);
	       env->mip->cru_vars_num = 0;
	    }
	 } else{
	    printf("sym_warm_solve():");
	    printf("Unable to re-solve this type of modification,for now!\n");
	    return(FUNCTION_TERMINATED_ABNORMALLY);
	 }
      }
   }

   /* Uncommented for now! */
#if 0
   if (env->par.trim_warm_tree) {
      trim_warm_tree(env, env->warm_start->rootnode);
   }
#endif

   FREE(cru_vars);
   FREE(cut_ind);
   FREE(tmp_ind);

   return(sym_solve(env));
}

/*===========================================================================*/
/* These data types are for multi-criteria problems and are only used here   */
/*===========================================================================*/

typedef struct SOLUTION_DATA{
   double  obj[2];
   double  gamma;
   double  tau;
   int     length;
   int    *indices;
   double *values;
}solution_data;

/*===========================================================================*/

typedef struct SOLUTION_PAIRS{
   int solution1;
   int solution2;
   double gamma1;
   double gamma2;
}solution_pairs;

/*===========================================================================*/

 typedef struct WS_ITEM{
   warm_start_desc * ws;
   struct WS_ITEM * next;
   double gamma;
}ws_item;

/*===========================================================================*/

#define MAX_NUM_PAIRS 10000
#define MAX_NUM_SOLUTIONS 10000
#define MAX_NUM_INFEASIBLE 10000

/*===========================================================================*/
/*===========================================================================*/

int sym_mc_solve(sym_environment *env)
{
   int i, cp_num;
   double gamma, gamma0, gamma1, tau, slope;
   double start_time;
   warm_start_desc *ws = NULL, *ws1 = NULL, *ws2 = NULL;
   ws_item *head = NULL, *tail = NULL, *item = NULL, *temp = NULL;
   solution_data solutions[MAX_NUM_PAIRS];
   int numsolutions = 0, numprobs = 0, numinfeasible = 0;
   solution_pairs pairs[MAX_NUM_PAIRS];
   int numpairs = 0, cur_position = 0, first = 0, last = 0, previous = 0;
   int *indices;
   double *values;
   int length, termcode;
   int solution1, solution2;
   double utopia[2];
   double compare_sol_tol, ub = 0.0;
   int binary_search = FALSE;

   for (i = 0; i < env->mip->n; i++){
      if (env->mip->obj2[i] != 0){
	 break;
      }
   }
   if (i == env->mip->n){
      printf("Second objective function is identically zero.\n");
      printf("Switching to standard branch and bound.\n\n");
      return(sym_solve(env));
   }

   sym_set_int_param(env, "multi_criteria", TRUE);
   memcpy((char *)env->mip->obj1, (char *)env->mip->obj, DSIZE*env->mip->n);
   if (!env->par.lp_par.mc_find_supported_solutions){
      env->base->cutnum += 2;
      env->rootdesc->uind.size++;
      env->rootdesc->uind.list = (int *) realloc(env->rootdesc->uind.list,
					 env->rootdesc->uind.size*ISIZE);
      env->rootdesc->uind.list[env->rootdesc->uind.size-1] = env->mip->n;
   }

   start_time = wall_clock(NULL);

   /* Set some parameters */
   compare_sol_tol = env->par.mc_compare_solution_tolerance;
   if (env->par.lp_par.mc_find_supported_solutions){
      env->par.lp_par.mc_rho = 0;
   }
   env->par.tm_par.granularity = env->par.lp_par.granularity =
      -MAX(env->par.lp_par.mc_rho, compare_sol_tol);

   if (env->par.mc_binary_search_tolerance > 0){
      binary_search = TRUE;
   }
   if (env->par.verbosity >= 0){
      if (env->par.mc_binary_search_tolerance > 0){
	 printf("Using binary search with tolerance = %f...\n",
		env->par.mc_binary_search_tolerance);
      }
      if (env->par.mc_search_order == MC_LIFO){
	 printf("Using LIFO search order...\n");
      }else{
	 printf("Using FIFO search order...\n");
      }
      if (env->par.lp_par.mc_rho > 0){
	 printf("Using augmented Chebyshev weight %.8f\n",
		env->par.lp_par.mc_rho);
      }
      if (env->par.use_permanent_cut_pools){
	 printf("Saving the global cut pool between iterations...\n");
	 sym_create_permanent_cut_pools(env, &cp_num);
      }
      printf("\n");
   }

   /* First, calculate the utopia point */
   env->par.lp_par.mc_gamma = 1.0;
   env->par.lp_par.mc_tau = -1.0;

   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = 1.0 tau = 0.0 \n");
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* For now, close reduced cost fixing and keep_warm_start param if
      warm starting is to be used...
   */
   if (env->par.lp_par.mc_find_supported_solutions &&
       env->par.mc_warm_start){
      sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
      sym_set_int_param(env, "keep_warm_start", TRUE);
   }

   /* Solve */
   env->utopia[0] = 0;
   env->utopia[1] = -MAXINT;
   if ((termcode = sym_solve(env)) < 0){
      env->base->cutnum -=2;
      env->rootdesc->uind.size--;
      return(termcode);
   }


   if(env->par.lp_par.mc_find_supported_solutions &&
      env->par.mc_warm_start){

      ws = sym_get_warm_start(env, FALSE);

      switch (env->par.mc_warm_start_rule){
       case 0:
       case 1:
	  sym_set_int_param(env, "do_reduced_cost_fixing", TRUE);
	  sym_set_int_param(env, "keep_warm_start", FALSE);
       case 2:
	  ws1=ws2=ws;
	  break;
       case 3:
	  head = tail = (ws_item *) calloc (1, sizeof(ws_item));
	  head->ws = tail->ws = ws;
	  head->gamma = tail->gamma = 1.0;
	  break;
       default:
	  break;
      }
   }
   numprobs++;

   /* Store the solution */
   length = solutions[numsolutions].length = env->best_sol.xlength;
   indices = solutions[numsolutions].indices = (int *) calloc(length, ISIZE);
   values = solutions[numsolutions].values = (double *) calloc(length, DSIZE);
   memcpy((char *) indices, env->best_sol.xind, length * ISIZE);
   memcpy((char *) values, env->best_sol.xval, length * DSIZE);
   solutions[numsolutions].gamma = 1.0;
   solutions[numsolutions].tau = 0.0;
   solutions[numsolutions].obj[0] = env->obj[0];
   solutions[numsolutions++].obj[1] = env->obj[1];
   utopia[0] = env->obj[0];

   env->par.lp_par.mc_gamma = -1.0;
   env->par.lp_par.mc_tau = 1.0;

   printf("***************************************************\n");
   printf("***************************************************\n");
   printf("Now solving with gamma = 0.0 tau = 1.0 \n");
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Resolve */
   env->utopia[0] = -MAXINT;
   env->utopia[1] = 0;
   if (env->par.lp_par.mc_find_supported_solutions){

      for (i = 0; i < env->mip->n; i++){
	 sym_set_obj_coeff(env, i, env->mip->obj2[i] +
			   env->par.lp_par.mc_rho*(env->mip->obj1[i] +
						   env->mip->obj2[i]));
      }

      if (env->par.mc_warm_start && env->par.mc_warm_start_rule == 0){
	 sym_set_warm_start(env, ws);
	 if ((termcode = sym_warm_solve(env)) < 0){
	    switch(env->par.mc_warm_start_rule){
	     case 0:
		sym_delete_warm_start(ws);
		break;
	     case 3:
		sym_delete_warm_start(head->ws);
		FREE(head);
		break;
	    }
	    return(termcode);
	 }
      }else{
	 if ((termcode = sym_solve(env)) < 0){
	    return(termcode);
	 }
      }
   }else{
      if ((termcode = sym_solve(env)) < 0){
	 env->base->cutnum -=2;
	 env->rootdesc->uind.size--;
	 return(termcode);
      }
   }

   if(env->par.lp_par.mc_find_supported_solutions &&
      env->par.mc_warm_start){

      switch (env->par.mc_warm_start_rule){
       case 0:
	  break;
       case 1:
       case 2:
	  ws2=sym_get_warm_start(env, FALSE);
	  break;
       case 3:
	  item = (ws_item *) calloc (1, sizeof(ws_item));
	  item->ws = sym_get_warm_start(env, FALSE);
	  item->gamma = 0.0;
	  head->next = tail = item;
	  break;
       default:
	  break;
      }
   }
   numprobs++;

   /* Store the solution */
   length = solutions[numsolutions].length = env->best_sol.xlength;
   indices = solutions[numsolutions].indices = (int *) calloc(length, ISIZE);
   values = solutions[numsolutions].values = (double *) calloc(length, DSIZE);
   memcpy((char *) indices, env->best_sol.xind, length * ISIZE);
   memcpy((char *) values, env->best_sol.xval, length * DSIZE);
   solutions[numsolutions].gamma = 0.0;
   solutions[numsolutions].tau = 1.0;
   solutions[numsolutions].obj[0] = env->obj[0];
   solutions[numsolutions++].obj[1] = env->obj[1];
   utopia[1] = env->obj[1];

   env->utopia[1] = utopia[1];
   env->utopia[0] = utopia[0];

   printf("***************************************************\n");
   printf("***************************************************\n");
   if(env->mip->obj_sense == SYM_MAXIMIZE){
      printf("Utopia point has first  objective value %.3f\n", -utopia[0]);
      printf("                 second objective value %.3f\n", -utopia[1]);
   }else{
      printf("Utopia point has first  objective value %.3f\n", utopia[0]);
      printf("                 second objective value %.3f\n", utopia[1]);
   }
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   /* Add the first pair to the list */
   if (solutions[0].obj[0] != solutions[1].obj[0]){
      if (binary_search){
	 pairs[first].gamma1 = 1.0;
	 pairs[first].gamma2 = 0.0;
      }
      pairs[first].solution1 = 0;
      pairs[first].solution2 = 1;
      first = last = 0;
      numpairs = 1;
   }else{
      numpairs = 0;
   }

   /* Keep taking pairs off the list and processing them until there are none
      left */
   while (numpairs > 0 && numpairs < MAX_NUM_PAIRS &&
	  numsolutions < MAX_NUM_SOLUTIONS &&
	  numinfeasible < MAX_NUM_INFEASIBLE){

      if (env->par.mc_search_order == MC_LIFO){
	 solution1 = pairs[last].solution1;
	 solution2 = pairs[last].solution2;
	 cur_position = last;
	 if (--last < 0){
	    last = MAX_NUM_PAIRS - 1;
	 }
	 numpairs--;
      }else{
	 solution1 = pairs[first].solution1;
	 solution2 = pairs[first].solution2;
	 cur_position = first;
	 if (++first > MAX_NUM_PAIRS-1)
	    first = 0;
	 numpairs--;
      }

      if (binary_search){
	 gamma = (pairs[cur_position].gamma1 + pairs[cur_position].gamma2)/2;
      }else if (!env->par.lp_par.mc_find_supported_solutions){
	 gamma = (utopia[1] - solutions[solution1].obj[1])/
	    (utopia[0] - solutions[solution2].obj[0] +
	     utopia[1] - solutions[solution1].obj[1]);
      }else{
	 slope = (solutions[solution1].obj[1] -
		  solutions[solution2].obj[1])/
	    (solutions[solution2].obj[0] -
	     solutions[solution1].obj[0]);
	 gamma = slope/(1+slope);
      }
      tau = 1 - gamma;

      env->par.lp_par.mc_gamma = gamma;
      env->par.lp_par.mc_tau = tau;

      /* Find upper bound */

      env->has_mc_ub = env->has_ub = FALSE;
      env->mc_ub = env->ub = MAXDOUBLE;
      if (!binary_search){
	 for (i = 0; i < numsolutions; i++){
	    if (!env->par.lp_par.mc_find_supported_solutions){
	       ub = MAX(gamma*(solutions[i].obj[0] - utopia[0]),
			tau*(solutions[i].obj[1] - utopia[1]));
	    }else{
	       ub = gamma*solutions[i].obj[0] + tau*solutions[i].obj[1] +
		  env->par.lp_par.mc_rho * (solutions[i].obj[0] +
					  solutions[i].obj[1]);
	    }
	    if (ub + env->par.lp_par.mc_rho * (solutions[i].obj[0] +
					     solutions[i].obj[1]) < env->ub){
	       env->has_mc_ub = env->has_ub = TRUE;
	       env->ub = ub + env->par.lp_par.mc_rho *
		  (solutions[i].obj[0] + solutions[i].obj[1]) - compare_sol_tol;
	       env->obj[0] = solutions[i].obj[0];
	       env->obj[1] = solutions[i].obj[1];
	       env->mc_ub = ub;
	    }
	 }
      }

      printf("***************************************************\n");
      printf("***************************************************\n");
      printf("Now solving with gamma = %.6f tau = %.6f \n", gamma, tau);
      printf("***************************************************\n");
      printf("***************************************************\n\n");

      env->obj[0] = env->obj[1] = 0.0;

      if (env->par.lp_par.mc_find_supported_solutions){

	 for (i = 0; i < env->mip->n; i++){
	    sym_set_obj_coeff(env, i, gamma*env->mip->obj1[i]
			      + tau*env->mip->obj2[i]
			      + env->par.lp_par.mc_rho*(env->mip->obj1[i]
							+ env->mip->obj2[i]));
	 }

	 if(env->par.mc_warm_start){

	    switch(env->par.mc_warm_start_rule){
	     case 0:
		sym_set_warm_start(env, ws);
		break;
	     case 1:
		if (gamma > 0.5 ){
		   sym_set_warm_start(env, ws1);
		} else {
		   sym_set_warm_start(env, ws2);
		}
		break;
	     case 2:
		sym_delete_warm_start(env->warm_start);
		if (gamma > 0.5 ){
		   env->warm_start = ws1;
		} else {
		   env->warm_start = ws2;
		}
		break;
	     case 3:

		if (gamma > 0.0){
		   for(item = head;;item = item->next){
		      if (gamma >= item->next->gamma){
			 break;
		      }
		   }
		   if ((item->gamma - gamma) <= (gamma - item->next->gamma)){
		      sym_set_warm_start(env, item->ws);
		   } else {
		      sym_set_warm_start(env, item->next->ws);
		   }
		} else {
		   sym_set_warm_start(env, tail->ws);
		}
		break;
	     default:
		break;
	    }

	    if ((termcode = sym_warm_solve(env)) < 0){

	       /* FIXME! copy best_sol.xind and .xval from env to warm_start*/
	       memset(&(env->best_sol), 0, sizeof(lp_sol));

	       switch(env->par.mc_warm_start_rule){
		case 0:
		   sym_delete_warm_start(ws);
		   break;
		case 1:
		case 2:
		   sym_delete_warm_start(ws1);
		   sym_delete_warm_start(ws2);
		   break;
		case 3:
		   while(TRUE){
		      sym_delete_warm_start(head->ws);
		      if(head != tail){
			 temp = head->next;
			 head->next = 0;
			 FREE(head);
			 head = temp;
		      } else {
			 FREE(head);
			 break;
		      }
		   }
		   break;
		default:
		   break;
	       }
	       return(termcode);
	    }
	 }else{
	    if ((termcode = sym_solve(env)) < 0){
	       return(termcode);
	    }
	 }
      } else{
	 if ((termcode = sym_solve(env)) < 0){
	    env->base->cutnum -=2;
	    env->rootdesc->uind.size--;
	    return(termcode);
	 }
      }

      switch(env->par.mc_warm_start_rule){
       case 0:
       case 1:
	  break;
       case 2:
	  if (gamma > 0.5 ){
	     ws1 = sym_get_warm_start(env,FALSE);
	  } else {
	     ws2 = sym_get_warm_start(env,FALSE);
	  }
	  break;
       case 3:
	  temp = (ws_item *) calloc (1, sizeof(ws_item));
	  temp->ws = sym_get_warm_start(env, FALSE);
	  temp->gamma = gamma;
	  temp->next = item->next;
	  item->next = temp;
	  break;
       default:
	  break;
      }

      numprobs++;

      if (binary_search){
	 if (env->obj[0] - solutions[solution1].obj[0] <
	     compare_sol_tol &&
	     solutions[solution1].obj[1] - env->obj[1] <
	     compare_sol_tol){
	    if (pairs[cur_position].gamma1 - gamma >
		env->par.mc_binary_search_tolerance){
	       if (++last > MAX_NUM_PAIRS - 1)
		  last = 0;
	       pairs[last].solution1 = solution1;
	       pairs[last].solution2 = solution2;
	       pairs[last].gamma1 = gamma;
	       pairs[last].gamma2 = pairs[cur_position].gamma2;
	       numpairs++;
	    }
	    continue;
	 }
	 if (solutions[solution2].obj[0] - env->obj[0] < compare_sol_tol
	     && env->obj[1] - solutions[solution2].obj[1] <
	     compare_sol_tol){
	    if (gamma - pairs[cur_position].gamma2 >
		env->par.mc_binary_search_tolerance){
	       if (++last > MAX_NUM_PAIRS - 1)
		  last = 0;
	       pairs[last].solution1 = solution1;
	       pairs[last].solution2 = solution2;
	       pairs[last].gamma1 = pairs[cur_position].gamma1;
	       pairs[last].gamma2 = gamma;
	       numpairs++;
	    }
	    continue;
	 }
      }else{
	 if (env->obj[0] == 0.0 && env->obj[1] == 0.0){
	    numinfeasible++;
	    continue;
	 }else if (env->obj[0] - solutions[solution1].obj[0] <
		   compare_sol_tol &&
		   solutions[solution1].obj[1] - env->obj[1] <
		   compare_sol_tol){
	    numinfeasible++;
	    continue;
	 }else if (solutions[solution2].obj[0] - env->obj[0] <
		   compare_sol_tol &&
		   env->obj[1] - solutions[solution2].obj[1] <
		   compare_sol_tol){
	    numinfeasible++;
	    continue;
	 }
      }

      /* Insert new solution */
      numinfeasible = 0;
      if (last + 2 == MAX_NUM_PAIRS){
	 last = 0;
	 previous = MAX_NUM_PAIRS - 1;
      }else if (last + 2 == MAX_NUM_PAIRS + 1){
	 last = 1;
	 previous = 0;
      }else{
	 last += 2;
	 previous = last - 1;
      }
      if (binary_search){
	 pairs[previous].gamma1 = pairs[cur_position].gamma1;
	 pairs[previous].gamma2 = gamma;
	 pairs[last].gamma1 = gamma;
	 pairs[last].gamma2 = pairs[cur_position].gamma2;
      }
      pairs[previous].solution1 = solution1;
      pairs[previous].solution2 = solution2;
      pairs[last].solution1 = solution2;
      pairs[last].solution2 = solution2+1;
      numpairs += 2;
      for (i = numsolutions; i > solution2; i--){
	 solutions[i] = solutions[i-1];
      }
      numsolutions++;
      if (env->par.mc_search_order == MC_FIFO){
	 if (first < last){
	    for (i = first; i < last - 1; i++){
	       if (pairs[i].solution1 >= solution2){
		  pairs[i].solution1++;
	       }
	       if (pairs[i].solution2 >= solution2){
		  pairs[i].solution2++;
	       }
	    }
	 }else{
	    for (i = first; i < MAX_NUM_PAIRS - (last == 0 ? 1 : 0); i++){
	       if (pairs[i].solution1 >= solution2){
		  pairs[i].solution1++;
	       }
	       if (pairs[i].solution2 >= solution2){
		  pairs[i].solution2++;
	       }
	    }
	    for (i = 0; i < last - 1; i++){
	       if (pairs[i].solution1 >= solution2){
		  pairs[i].solution1++;
	       }
	       if (pairs[i].solution2 >= solution2){
		  pairs[i].solution2++;
	       }
	    }
	 }
      }

      length = solutions[solution2].length = env->best_sol.xlength;
      indices = solutions[solution2].indices = (int *) calloc(length, ISIZE);
      values = solutions[solution2].values = (double *) calloc(length, DSIZE);
      memcpy((char *) indices, env->best_sol.xind, length * ISIZE);
      memcpy((char *) values, env->best_sol.xval, length * DSIZE);
      solutions[solution2].gamma = gamma;
      solutions[solution2].tau = tau;
      solutions[solution2].obj[0] = env->obj[0];
      solutions[solution2].obj[1] = env->obj[1];
   }

   printf("\n********************************************************\n");

   if (numsolutions >= MAX_NUM_SOLUTIONS){
      printf("Maximum number of solutions (%i) reached\n\n",
	     MAX_NUM_SOLUTIONS);
   }

   if (numinfeasible >= MAX_NUM_INFEASIBLE){
      printf("Maximum number of infeasible subproblems (%i) reached\n\n",
	     MAX_NUM_INFEASIBLE);
   }

   if (numpairs >= MAX_NUM_PAIRS){
      printf("Maximum number of solution pairs (%i) reached\n\n",
	     MAX_NUM_PAIRS);
      printf("\n********************************************************\n");
      if (!env->par.lp_par.mc_find_supported_solutions){
	 printf(  "* Found set of non-dominated solutions!!!!!!! *\n");
      }else{
	 printf(  "* Found set of supported solutions!!!!!!!     *\n");
      }
   }else{
      printf("\n********************************************************\n");
      if (!env->par.lp_par.mc_find_supported_solutions){
	 printf(  "* Found complete set of non-dominated solutions!!!!!!! *\n");
      }else{
	 printf(  "* Found complete set of supported solutions!!!!!!!     *\n");
      }
   }
   printf(  "* Now displaying stats...                              *\n");
   printf(  "********************************************************\n\n");

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_CP)
   if (env->par.use_permanent_cut_pools){
      for (i = 0; i < env->par.tm_par.max_cp_num; i++){
	 env->comp_times.bc_time.cut_pool += env->cp[i]->cut_pool_time;
	 env->warm_start->stat.cuts_in_pool += env->cp[i]->cut_num;
      }
   }
#endif

   if (!env->par.multi_criteria){
      print_statistics(&(env->comp_times.bc_time), &(env->warm_start->stat),
                       NULL, 0.0, 0.0, 0, start_time, wall_clock(NULL),
                       env->mip->obj_offset, env->mip->obj_sense, env->has_ub,
                       NULL);
   } else{
      printf("Total WallClock Time         %.3f\n", wall_clock(NULL) -
	     start_time);
   }

   printf("\nNumber of subproblems solved: %i\n", numprobs);
   printf("Number of solutions found: %i\n\n", numsolutions);

   printf("***************************************************\n");
   printf("***************************************************\n");
   if (!env->par.lp_par.mc_find_supported_solutions){
      printf("Displaying non-dominated solution values and breakpoints\n");
   }else{
      printf("Displaying supported solution values and breakpoints\n");
   }
   printf("***************************************************\n");
   printf("***************************************************\n\n");

   gamma0 = 1.0;
   for (i = 0; i < numsolutions - 1; i++){
      if (!env->par.lp_par.mc_find_supported_solutions){
	 gamma1 = (utopia[1] - solutions[i].obj[1])/
	    (utopia[0] - solutions[i+1].obj[0] +
	     utopia[1] - solutions[i].obj[1]);
      }else{
	 slope = (solutions[i].obj[1] -
		  solutions[i+1].obj[1])/
	    (solutions[i+1].obj[0] -
	     solutions[i].obj[0]);
	 gamma1 = slope/(1+slope);
      }
      if(env->mip->obj_sense == SYM_MAXIMIZE){
	 printf("First Objective: %.3f Second Objective: %.3f ",
		-solutions[i].obj[0], -solutions[i].obj[1]);
      }else{
	 printf("First Objective: %.3f Second Objective: %.3f ",
		solutions[i].obj[0], solutions[i].obj[1]);
      }

      printf("Range: %.6f - %.6f\n", gamma1, gamma0);
      gamma0 = gamma1;
   }
   if(env->mip->obj_sense == SYM_MAXIMIZE){
      printf("First Objective: %.3f Second Objective: %.3f ",
	     -solutions[i].obj[0], -solutions[i].obj[1]);
   }else{
      printf("First Objective: %.3f Second Objective: %.3f ",
	     solutions[i].obj[0], solutions[i].obj[1]);
   }

   printf("Range: %.6f - %.6f\n", 0.0, gamma0);

   for (i = 0 ; i < numsolutions; i++){
      FREE(solutions[i].values);
      FREE(solutions[i].indices);
   }
   if (env->par.lp_par.mc_find_supported_solutions){
       if(env->par.mc_warm_start){

	  /* FIXME! copy best_sol.xind and .xval from env to warm_start*/
	  memset(&(env->best_sol), 0, sizeof(lp_sol));

	  switch(env->par.mc_warm_start_rule){
	   case 0:
	      sym_delete_warm_start(ws);
	      break;
	   case 1:
	   case 2:
	      sym_delete_warm_start(ws1);
	      sym_delete_warm_start(ws2);
	      break;
	   case 3:
	      while(TRUE){
		 sym_delete_warm_start(head->ws);
		 if(head != tail){
		    temp = head->next;
		    head->next = 0;
		    FREE(head);
		    head = temp;
		 } else {
		    FREE(head);
		    break;
		 }
	      }
	      break;
	   default:
	      break;
	  }
       }
   }else {
      env->base->cutnum -=2;
      env->rootdesc->uind.size--;
   }

   return(TM_OPTIMAL_SOLUTION_FOUND);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_create_permanent_cut_pools(sym_environment *env, int * cp_num)
{

   *cp_num = 0;

#if !(defined(COMPILE_IN_TM) && defined(COMPILE_IN_CP))
   return(FUNCTION_TERMINATED_ABNORMALLY);
#else
   int i;

   if (env->par.tm_par.max_cp_num){
      env->cp =
	 (cut_pool **) malloc(env->par.tm_par.max_cp_num*sizeof(cut_pool *));
      for (i = 0; i < env->par.tm_par.max_cp_num; i++){
	 env->cp[i] = (cut_pool *) calloc(1, sizeof(cut_pool));
	 env->cp[i]->par = env->par.cp_par;
#ifdef USE_SYM_APPLICATION
	 CALL_USER_FUNCTION( user_send_cp_data(env->user, &env->cp[i]->user) );
#else
	 env->cp[i]->user = env->user;
#endif
      }
      *cp_num = env->par.tm_par.max_cp_num;
      return(FUNCTION_TERMINATED_NORMALLY);
   }else{
      printf("sym_create_permanent_cut_pools(): \"max_cp_num\" param was not set!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int sym_close_environment(sym_environment *env)
{
   int termcode = 0;

   CALL_WRAPPER_FUNCTION( free_master_u(env) );

   FREE(env);

#if (!defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP) ||                   \
    !defined(COMPILE_IN_CG) || !defined(COMPILE_IN_CP)) && defined(__PVM__)
   pvm_catchout(0);
   comm_exit();
#endif

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_explicit_load_problem(sym_environment *env, int numcols, int numrows,
			      int *start, int *index, double *value,
			      double *collb, double *colub, char *is_int,
			      double *obj, double *obj2, char *rowsen,
			      double *rowrhs, double *rowrng, char make_copy)
{
   int termcode = 0;
   double t = 0, inf = SYM_INFINITY;
   int i = 0;

   if ((!numcols && !numrows) || numcols < 0 || numrows <0){
      printf("sym_explicit_load_problem():The given problem is empty or incorrect ");
      printf("problem description!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   (void)used_time(&t);

   env->mip->m  = numrows;
   env->mip->n  = numcols;

   if (make_copy){

      if(numcols){
	 env->mip->obj    = (double *) calloc(numcols, DSIZE);
	 env->mip->obj1   = (double *) calloc(numcols, DSIZE);
	 env->mip->obj2   = (double *) calloc(numcols, DSIZE);
	 env->mip->ub     = (double *) calloc(numcols, DSIZE);
	 env->mip->lb     = (double *) calloc(numcols, DSIZE);
	 env->mip->is_int = (char *)   calloc(CSIZE, numcols);

	 if (obj){
	    memcpy(env->mip->obj,  obj,  DSIZE * numcols);
	 }

	 if (obj2){
	    memcpy(env->mip->obj2, obj2, DSIZE * numcols);
	 }

	 if (colub){
	    memcpy(env->mip->ub, colub, DSIZE * numcols);
	 }else{
	    for(i = 0; i<env->mip->n; i++){
	       env->mip->ub[i] = inf;
	    }
	 }

	 if(collb){
	    memcpy(env->mip->lb, collb, DSIZE * numcols);
	 }

	 if (is_int){
	    memcpy(env->mip->is_int, is_int, CSIZE * numcols);
	 }
      }

      if(numrows){

	 env->mip->rhs    = (double *) calloc(numrows, DSIZE);
	 env->mip->sense  = (char *)   malloc(CSIZE * numrows);
	 env->mip->rngval = (double *) calloc(numrows, DSIZE);

	 if (rowsen){
	    memcpy(env->mip->sense, rowsen, CSIZE * numrows);
	 }else{
	    memset(env->mip->sense, 'N', CSIZE *numrows);
	 }

	 if(rowrhs){
	    memcpy(env->mip->rhs, rowrhs, DSIZE * numrows);
	 }

	 if (rowrng){
	    memcpy(env->mip->rngval, rowrng, DSIZE * numrows);
	 }
      }

      //user defined matind, matval, matbeg--fill as column ordered

      if(start){

	 env->mip->nz = start[numcols];
	 env->mip->matbeg = (int *) calloc(ISIZE, (numcols + 1));
	 env->mip->matval = (double *) calloc(DSIZE,start[numcols]);
	 env->mip->matind = (int *)    calloc(ISIZE,start[numcols]);

	 memcpy(env->mip->matbeg, start, ISIZE *(numcols + 1));
	 memcpy(env->mip->matval, value, DSIZE *start[numcols]);
	 memcpy(env->mip->matind, index, ISIZE *start[numcols]);
      }

   }else{

      if (obj){
	 env->mip->obj = obj;
      }else{
	 env->mip->obj    = (double *) calloc(numcols, DSIZE);
      }

      env->mip->obj1   = (double *) calloc(numcols, DSIZE);

      if (obj2){
	 env->mip->obj2 = obj2;
      }else{
	 env->mip->obj2   = (double *) calloc(numcols, DSIZE);
      }

      if (rowsen){
	 env->mip->sense = rowsen;
      }else{
	 env->mip->sense  = (char *) malloc(CSIZE * numrows);
	 memset(env->mip->sense, 'N', CSIZE *numrows);
      }

      if(rowrhs){
	 env->mip->rhs = rowrhs;
      }else{
	 env->mip->rhs = (double *) calloc(numrows, DSIZE);
      }

      if (rowrng){
	 env->mip->rngval = rowrng;
      }else{
	 env->mip->rngval = (double *) calloc(numrows, DSIZE);
      }

      if (colub){
	 env->mip->ub = colub;
      }else{
	 env->mip->ub = (double *) calloc(numcols, DSIZE);
	 for(i = 0; i<env->mip->n; i++){
	    env->mip->ub[i] = inf;
	 }
      }

      if (collb){
	 env->mip->lb = collb;
      }else{
	 env->mip->lb = (double *) calloc(numcols, DSIZE);
      }

      if (is_int){
	 env->mip->is_int = is_int;
      }else{
	 env->mip->is_int = (char *)   calloc(CSIZE, numcols);
      }

      if(start){
	 env->mip->nz = start[numcols];
	 env->mip->matbeg = start;
	 env->mip->matval = value;
	 env->mip->matind = index;
      }
   }

   /* Start up the graphics window*/
#if !defined(_MSC_VER) && !defined (__MNO_CYGWIN)
   CALL_WRAPPER_FUNCTION( init_draw_graph_u(env) );
#endif

   if(env->mip->obj_sense == SYM_MAXIMIZE){
      for (i = 0; i < numcols; i++){
	 env->mip->obj[i] *= -1.0;
	 env->mip->obj2[i] *= -1.0;
      }
   }

   /*------------------------------------------------------------------------*\
    * Have the user generate the base and root description
   \*------------------------------------------------------------------------*/

   CALL_WRAPPER_FUNCTION(initialize_root_node_u(env) );

   env->comp_times.readtime = used_time(&t);

   env->termcode = TM_NO_SOLUTION;
   env->mip->is_modified = TRUE;

   return termcode;
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_abandoned(sym_environment *env)
{

   switch(env->termcode){
    case PREP_ERROR:
    case SOMETHING_DIED:
    case TM_ERROR__NUMERICAL_INSTABILITY:
       return(TRUE);
    default:
       break;
   }

   return(FALSE);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_proven_optimal(sym_environment *env)
{

   switch(env->termcode){
    case TM_OPTIMAL_SOLUTION_FOUND:
    case PREP_OPTIMAL_SOLUTION_FOUND:
      return(TRUE);
    default:
       break;
   }

   return(FALSE);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_proven_primal_infeasible(sym_environment *env)
{

   switch(env->termcode){
    case TM_NO_SOLUTION:
    case PREP_NO_SOLUTION:
      return(TRUE);
    default:
       break;
   }

   return(FALSE);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_iteration_limit_reached(sym_environment *env)
{

   switch(env->termcode){
    case TM_NODE_LIMIT_EXCEEDED:
    case TM_FOUND_FIRST_FEASIBLE:
       return(TRUE);
    default:
       break;
   }

   return(FALSE);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_time_limit_reached(sym_environment *env)
{

   switch(env->termcode){
    case TM_TIME_LIMIT_EXCEEDED:
       return(TRUE);
    default:
       break;
   }

   return(FALSE);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_target_gap_achieved(sym_environment *env)
{

   switch(env->termcode){
    case TM_TARGET_GAP_ACHIEVED:
    case TM_OPTIMAL_SOLUTION_FOUND:
       return(TRUE);
    default:
       break;
   }

   return(FALSE);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_status(sym_environment *env)
{
   return (env->termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_num_cols(sym_environment *env, int *numcols)
{

   if (!env->mip){
      if(env->par.verbosity >= 1){
	 printf("sym_get_num_cols():There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *numcols = env->mip->n;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_num_rows(sym_environment *env, int *numrows)
{
   if (!env->mip){
      if(env->par.verbosity >= 1){
	 printf("sym_get_num_rows():There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *numrows = env->mip->m;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_num_elements(sym_environment *env, int *numelems)
{
   if (!env->mip){
      if(env->par.verbosity >= 1){
	 printf("sym_get_num_elements():There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *numelems = env->mip->nz;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_col_lower(sym_environment *env, double *collb)
{
   if (!env->mip || !env->mip->n || !env->mip->lb){
      if(env->par.verbosity >= 1){
	 printf("sym_get_col_lower():There is no loaded mip description or\n");
	 printf("there is no loaded column description!\n");
      }
     return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   memcpy(collb, env->mip->lb, DSIZE * env->mip->n);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_col_upper(sym_environment *env, double *colub)
{
   if (!env->mip || !env->mip->n || !env->mip->ub){
      if(env->par.verbosity >= 1){
	 printf("sym_get_col_upper():There is no loaded mip description or\n");
	 printf("there is no loaded column description!\n");
      }
    return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   memcpy(colub, env->mip->ub, DSIZE * env->mip->n);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_row_sense(sym_environment *env, char *rowsen)
{
   if (!env->mip || !env->mip->m || !env->mip->sense){
      if(env->par.verbosity >= 1){
	 printf("sym_get_row_sense():There is no loaded mip description or\n");
	 printf("there is no loaded row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   memcpy(rowsen, env->mip->sense, CSIZE * env->mip->m);

   return(FUNCTION_TERMINATED_NORMALLY);

}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_rhs(sym_environment *env, double *rowrhs)
{
   if (!env->mip || !env->mip->m || !env->mip->rhs){
      if(env->par.verbosity >= 1){
	 printf("sym_get_rhs():There is no loaded mip description or\n");
	 printf("there is no loaded row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   memcpy(rowrhs, env->mip->rhs, DSIZE * env->mip->m);

   return(FUNCTION_TERMINATED_NORMALLY);

}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_matrix(sym_environment *env, int *nz, int *matbeg, int *matind,
		   double *matval)
{
   if (!env->mip || !env->mip->m || !env->mip->n || !env->mip->matbeg){
      if(env->par.verbosity >= 1){
	 printf("sym_get_rhs():There is no loaded mip description or\n");
	 printf("there is no loaded matrix description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *nz = env->mip->nz;

   memcpy(matbeg, env->mip->matbeg, ISIZE * (env->mip->n + 1));
   memcpy(matind, env->mip->matind, ISIZE * (*nz));
   memcpy(matval, env->mip->matval, DSIZE * (*nz));

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_row_range(sym_environment *env, double *rowrng)
{
   if (!env->mip || !env->mip->m){
      if(env->par.verbosity >= 1){
	 printf("sym_get_row_range():There is no loaded mip description or\n");
	 printf("there is no loaded row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   memcpy(rowrng, env->mip->rngval, DSIZE * env->mip->m);

   return(FUNCTION_TERMINATED_NORMALLY);

}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_row_lower(sym_environment *env, double *rowlb)
{
   if (!env->mip || !env->mip->m || !env->mip->rhs){
      if(env->par.verbosity >= 1){
	 printf("sym_get_row_lower():There is no loaded mip description or\n");
	 printf("there is no loaded row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   double upper;
   int i;
   double rhs, range, inf = SYM_INFINITY;
   char sense;

   for ( i = env->mip->m - 1; i >= 0; --i ){

	 rhs   = env->mip->rhs[i];
	 range = env->mip->rngval[i];
	 sense = env->mip->sense[i];

	 switch (sense) {
	  case 'E':
	     rowlb[i] = upper = rhs;
	     break;
	  case 'L':
	     rowlb[i] = -inf;
	     upper = rhs;
	     break;
	  case 'G':
	     rowlb[i] = rhs;
	     upper = inf;
	     break;
	  case 'R':
	     rowlb[i] = rhs - range;
	     upper = rhs;
	     break;
	  case 'N':
	     rowlb[i] = -inf;
	     upper = inf;
	     break;
	 }
      }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_row_upper(sym_environment *env, double *rowub)
{
   if (!env->mip || !env->mip->m || !env->mip->rhs){
      if(env->par.verbosity >= 1){
	 printf("sym_get_row_upper():There is no loaded mip description or\n");
	 printf("there is no loaded row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   double lower;
   int i;
   double rhs, range, inf = SYM_INFINITY;
   char sense;

   for ( i = env->mip->m - 1; i >= 0; --i )
      {

	 rhs   = env->mip->rhs[i];
	 range = env->mip->rngval[i];
	 sense = env->mip->sense[i];

	 switch (sense) {
	  case 'E':
	     lower = rowub[i] = rhs;
	     break;
	  case 'L':
	     lower = -inf;
	     rowub[i] = rhs;
	     break;
	  case 'G':
	     lower = rhs;
	     rowub[i] = inf;
	     break;
	  case 'R':
	     lower = rhs - range;
	     rowub[i] = rhs;
	     break;
	  case 'N':
	     lower = -inf;
	     rowub[i] = inf;
	     break;
	 }
      }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_obj_coeff(sym_environment *env, double *obj)
{
   if (!env->mip || !env->mip->n || !env->mip->obj){
      if(env->par.verbosity >= 1){
	 printf("sym_get_obj_coeff():There is no loaded mip description or\n");
	 printf("there is no loaded obj vector description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   memcpy(obj, env->mip->obj, DSIZE*env->mip->n);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_obj2_coeff(sym_environment *env, double *obj2)
{
   if (!env->mip || !env->mip->n || !env->mip->obj2){
      if(env->par.verbosity >= 1){
	 printf("sym_get_obj2_coeff():There is no loaded mip description or\n");
	 printf("or there is no loaded second obj vector description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   memcpy(obj2, env->mip->obj2, DSIZE*env->mip->n);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_obj_sense(sym_environment *env, int *sense)
{
   if (!env->mip){
      if(env->par.verbosity >= 1){
	 printf("sym_get_obj_sense():There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *sense = 1;

   if (env->mip->obj_sense == SYM_MAXIMIZE) {
      *sense = -1;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_continuous(sym_environment *env, int index, int *value)
{
   if (!env->mip || index < 0 || index > env->mip->n || !env->mip->n ||
       !env->mip->is_int){
      if(env->par.verbosity >= 1){
	 printf("sym_is_continuous():There is no loaded mip description or\n");
	 printf("index is out of range or no column description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *value = FALSE;

   if (!env->mip->is_int[index]){
      *value = TRUE;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_binary(sym_environment *env, int index, int *value)
{
   if (!env->mip || index < 0 || index >= env->mip->n){
      if(env->par.verbosity >= 1){
	 printf("sym_is_binary(): Index out of range\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   if (!env->mip->n || !env->mip->is_int || !env->mip->ub || !env->mip->lb){
      if(env->par.verbosity >= 1){
	 printf("sym_is_binary(): There is no loaded mip description\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *value = FALSE;

   if (env->mip->is_int[index] && env->mip->lb[index] == 0.0 &&
       env->mip->ub[index] == 1.0) {
      *value = TRUE;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_is_integer(sym_environment *env, int index, char *value)
{
   if (!env->mip || index < 0 || index >= env->mip->n){
      if(env->par.verbosity >= 1){
	 printf("sym_is_binary(): Index out of range\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   if (!env->mip->n || !env->mip->is_int){
      if(env->par.verbosity >= 1){
	 printf("sym_is_binary(): There is no loaded mip description\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *value = env->mip->is_int[index];

   return(FUNCTION_TERMINATED_NORMALLY);

}

/*===========================================================================*/
/*===========================================================================*/

double sym_get_infinity()
{
   return(SYM_INFINITY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_col_solution(sym_environment *env, double *colsol)
{
   int i;
   lp_sol sol;

   sol = env->best_sol;

   if (!sol.xlength || sol.xlength && (!sol.xind || !sol.xval)){
      if(env->par.verbosity >= 1){
	 printf("sym_get_col_solution(): There is no solution!\n");
      }
      if(env->mip->n){
	 memcpy(colsol, env->mip->lb, DSIZE*env->mip->n);
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   else{
      if (!sol.has_sol){
	 printf("sym_get_col_solution(): Stored solution may not be feasible!\n");
      }
      memset(colsol, 0, DSIZE*env->mip->n);
      if(sol.xlength){
	 if(!env->prep_mip){
	    for( i = 0; i<sol.xlength; i++){
	       colsol[sol.xind[i]] = sol.xval[i];
	    }
	 }else{
	    for( i = 0; i<sol.xlength; i++){
	       colsol[env->prep_mip->orig_ind[sol.xind[i]]] = sol.xval[i];
	    }
	    for(i = 0; i < env->prep_mip->fixed_n; i++){
	       colsol[env->prep_mip->fixed_ind[i]] =
		  env->prep_mip->fixed_val[i];
	    }
	 }
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_row_activity(sym_environment *env, double *rowact)
{
   double * colsol;
   int i, j;

   int * matbeg;
   double * matval;
   int * matind;


   if (!env->mip || !env->mip->n){
      if(env->par.verbosity >= 1){
	 printf("sym_get_row_activity():There is no loaded mip description or\n");
	 printf("no column description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   colsol = (double *)malloc(DSIZE*env->mip->n);

   sym_get_col_solution(env, colsol);

   matbeg = env->mip->matbeg;
   matval = env->mip->matval;
   matind = env->mip->matind;

   memset(rowact, 0, DSIZE*env->mip->m);

   for(i = 0; i<env->mip->n; i++){
      for(j = matbeg[i]; j<matbeg[i+1]; j++){
	 rowact[matind[j]] += matval[j] * colsol[i];
      }
   }

   return (FUNCTION_TERMINATED_NORMALLY);

}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_obj_val(sym_environment *env, double *objval)
{
   int i;

   if (env->best_sol.has_sol){
      *objval = (env->mip->obj_sense == SYM_MINIMIZE ? env->best_sol.objval :
		 -env->best_sol.objval) +
	 (env->prep_mip ? env->prep_mip->obj_offset : env->mip->obj_offset);
   }else{
      if(env->par.verbosity >= 1){
	 printf("sym_get_obj_val(): There is no solution!\n");
      }
      /* return collb * objcoeff! */
      *objval = 0;
      for(i = 0; i<env->mip->n; i++){
	 *objval += env->mip->obj[i] * env->mip->lb[i];
      }
      *objval = env->mip->obj_sense ==
	 SYM_MINIMIZE ? *objval : -(*objval) ;
      return (FUNCTION_TERMINATED_ABNORMALLY);
   }

   return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_primal_bound(sym_environment *env, double *ub)
{

   if (!env->mip){
      if(env->par.verbosity >= 1){
	 printf("sym_get_primal_bound():There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   if (env->has_ub){
      *ub = env->mip->obj_sense == SYM_MINIMIZE ? env->ub : -env->ub ;
   }else{
      *ub = env->mip->obj_sense == SYM_MINIMIZE ? SYM_INFINITY : -SYM_INFINITY;
   }

   return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_iteration_count(sym_environment *env, int *numnodes)
{
   if (!env->warm_start){
      if(env->par.verbosity >= 1){
	 printf("sym_get_iteration_count():");
	 printf("There is no post-solution information available!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   *numnodes = env->warm_start->stat.analyzed;

   return (FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_obj_coeff(sym_environment *env, int index, double value)
{

   int i;

   if (!env->mip || !env->mip->n || index > env->mip->n || index < 0 ||
       !env->mip->obj){
      if(env->par.verbosity >= 1){
	 printf("sym_set_obj_coeff():There is no loaded mip description or\n");
	 printf("index is out of range or no column description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   if (env->mip->obj_sense == SYM_MAXIMIZE){
      env->mip->obj[index] = -value;
   }else{
      env->mip->obj[index] = value;
   }

   if (env->mip->change_num){
      for(i = env->mip->change_num - 1 ; i >=0 ; i--){
	 if (env->mip->change_type[i] == OBJ_COEFF_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 env->mip->change_type[env->mip->change_num++] = OBJ_COEFF_CHANGED;
      }
   }
   else{
      env->mip->change_type[env->mip->change_num++] = OBJ_COEFF_CHANGED;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_obj2_coeff(sym_environment *env, int index, double value)
{

   if (!env->mip || !env->mip->n || index > env->mip->n || index < 0 ||
       !env->mip->obj2){
      if(env->par.verbosity >= 1){
	 printf("sym_set_obj_coeff():There is no loaded mip description or\n");
	 printf("index is out of range or no column description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   if (env->mip->obj_sense == SYM_MAXIMIZE){
      env->mip->obj2[index] = -value;
   }else{
      env->mip->obj2[index] = value;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_lower(sym_environment *env, int index, double value)
{
   int i;

   if (!env->mip || !env->mip->n || index > env->mip->n || index < 0 ||
       !env->mip->lb){
      if(env->par.verbosity >= 1){
	 printf("sym_set_col_lower():There is no loaded mip description or\n");
	 printf("index is out of range or no column description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   env->mip->lb[index] = value;

   if (env->mip->change_num){
      for(i = env->mip->change_num - 1 ; i >=0 ; i--){
	 if (env->mip->change_type[i] == COL_BOUNDS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 env->mip->change_type[env->mip->change_num++] = COL_BOUNDS_CHANGED;
      }
   }
   else{
      env->mip->change_type[env->mip->change_num++] = COL_BOUNDS_CHANGED;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_upper(sym_environment *env, int index, double value)
{
   int i;

   if (!env->mip || !env->mip->n || index > env->mip->n || index < 0 ||
       !env->mip->ub){
      if(env->par.verbosity >= 1){
	 printf("sym_set_col_upper():There is no loaded mip description!\n");
	 printf("index is out of range or no column description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   env->mip->ub[index] = value;

   if (env->mip->change_num){
      for(i = env->mip->change_num - 1 ; i >=0 ; i--){
	 if (env->mip->change_type[i] == COL_BOUNDS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 env->mip->change_type[env->mip->change_num++] = COL_BOUNDS_CHANGED;
      }
   }
   else{
      env->mip->change_type[env->mip->change_num++] = COL_BOUNDS_CHANGED;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_row_lower(sym_environment *env, int index, double value)
{
   double rhs, range, lower = 0, upper = 0, inf = SYM_INFINITY;
   char   sense;
   int i;

   if (!env->mip || !env->mip->m || index > env->mip->m || index < 0 ||
       !env->mip->rhs){
      if(env->par.verbosity >= 1){
	 printf("sym_set_row_lower():There is no loaded mip description or\n");
	 printf("index is out of range or no row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   rhs   = env->mip->rhs[index];
   range = env->mip->rngval[index];
   sense = env->mip->sense[index];

   switch (sense) {
    case 'E':
       lower = upper = rhs;
       break;
    case 'L':
       lower = -inf;
       upper = rhs;
       break;
    case 'G':
       lower = rhs;
       upper = inf;
       break;
    case 'R':
       lower = rhs - range;
       upper = rhs;
       break;
    case 'N':
       lower = -inf;
       upper = inf;
       break;
   }

   if ( lower != value ) {
      lower = value;
      range = 0.0;
      if (lower > -inf) {
	 if (upper < inf) {
	    rhs = upper;
	    if (upper==lower) {
	       sense = 'E';
	    }else{
	       sense = 'R';
	       range = upper - lower;
	    }
	 }else{
	    sense = 'G';
	    rhs = lower;
	 }
      }else{
	 if (upper < inf) {
	    sense = 'L';
	    rhs = upper;
	 }else{
	    sense = 'N';
	    rhs = 0.0;
	 }
      }

      env->mip->sense[index] = sense;
      env->mip->rhs[index] = rhs;
      env->mip->rngval[index] = range;
   }

   if (env->mip->change_num){
      for(i = env->mip->change_num - 1 ; i >=0 ; i--){
	 if (env->mip->change_type[i] == RHS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 env->mip->change_type[env->mip->change_num] = RHS_CHANGED;
	 env->mip->change_num++;
      }
   }
   else{
      env->mip->change_type[env->mip->change_num] = RHS_CHANGED;
      env->mip->change_num++;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_row_upper(sym_environment *env, int index, double value)
{
   double rhs, range, lower = 0, upper = 0, inf = SYM_INFINITY;
   char   sense;
   int i;

   if (!env->mip || !env->mip->m || index > env->mip->m || index < 0 ||
       !env->mip->rhs){
      if(env->par.verbosity >= 1){
	 printf("sym_set_row_upper():There is no loaded mip description or\n");
	 printf("index is out of range or no row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   rhs   = env->mip->rhs[index];
   range = env->mip->rngval[index];
   sense = env->mip->sense[index];

   switch (sense) {
    case 'E':
       lower = upper = rhs;
       break;
    case 'L':
       lower = -inf;
       upper = rhs;
       break;
    case 'G':
       lower = rhs;
       upper = inf;
       break;
    case 'R':
       lower = rhs - range;
       upper = rhs;
       break;
    case 'N':
       lower = -inf;
       upper = inf;
       break;
   }

   /*   convertSenseToBound( sense, rhs, range, lower, upper );*/

   if ( upper != value ) {
      upper = value;
      /* convertBountToSense */
      range = 0.0;
      if (lower > -inf) {
	 if (upper < inf) {
	    rhs = upper;
	    if (upper==lower) {
	       sense = 'E';
	    }else{
	       sense = 'R';
	       range = upper - lower;
	    }
	 }else{
	    sense = 'G';
	    rhs = lower;
	 }
      }else{
	 if (upper < inf) {
	    sense = 'L';
	    rhs = upper;
	 }else{
	    sense = 'N';
	    rhs = 0.0;
	 }
      }

      env->mip->sense[index] = sense;
      env->mip->rhs[index] = rhs;
      env->mip->rngval[index] = range;
   }

   if (env->mip->change_num){
      for(i = env->mip->change_num - 1 ; i >=0 ; i--){
	 if (env->mip->change_type[i] == RHS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 env->mip->change_type[env->mip->change_num] = RHS_CHANGED;
	 env->mip->change_num++;

      }
   }
   else{
      env->mip->change_type[env->mip->change_num] = RHS_CHANGED;
      env->mip->change_num++;

   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_row_type(sym_environment *env, int index, char rowsense,
		     double rowrhs, double rowrng)

{

   int i;

   if (!env->mip || !env->mip->m || index > env->mip->m || index < 0 ||
       !env->mip->rhs){
      if(env->par.verbosity >= 1){
	 printf("sym_set_row_type():There is no loaded mip description or\n");
	 printf("index is out of range or no row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   env->mip->sense[index] = rowsense;
   env->mip->rhs[index] = rowrhs;
   env->mip->rngval[index] = rowrng;


   if (env->mip->change_num){
      for(i = env->mip->change_num - 1 ; i >=0 ; i--){
	 if (env->mip->change_type[i] == RHS_CHANGED){
	    break;
	 }
      }
      if (i < 0 ){
	 env->mip->change_type[env->mip->change_num] = RHS_CHANGED;
	 env->mip->change_num++;

      }
   }
   else{
      env->mip->change_type[env->mip->change_num] = RHS_CHANGED;
      env->mip->change_num++;

   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_obj_sense(sym_environment *env, int sense)
{

   int i;

   if (!env->mip){
      if(env->par.verbosity >= 1){
	 printf("sym_set_obj_type():There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   if (sense==-1 && env->mip->obj_sense != SYM_MAXIMIZE ){
      for (i = 0; i < env->mip->n; i++){
	 env->mip->obj[i] *= -1.0;
	 env->mip->obj2[i] *= -1.0;
      }
      env->mip->obj_sense = SYM_MAXIMIZE;
   }
   else if (sense != -1 && env->mip->obj_sense != SYM_MINIMIZE ){
      /* assume it to be min problem */
      for (i = 0; i < env->mip->n; i++){
	 env->mip->obj[i] *= -1.0;
	 env->mip->obj2[i] *= -1.0;
      }
      env->mip->obj_sense = SYM_MINIMIZE;
   }

   //else - do nothing!

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_solution(sym_environment *env, double * colsol)
{
   int i, j, nz = 0,*matBeg, *matInd;
   double value, *rowAct = NULL, *matVal;
   char feasible;
   double lpetol =  9.9999999999999995e-07;
   lp_sol * sol;
   int * tmp_ind;

   if (!env->mip || !env->mip->n){
      if(env->par.verbosity >= 1){
	 printf("sym_set_col_solution(): There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   /* Feasibility Check*/

   /* step 1. check for bounds and integrality */
   for (i = env->mip->n - 1; i >= 0; i--){
      if (colsol[i] < env->mip->lb[i] - lpetol ||
	  colsol[i] > env->mip->ub[i] + lpetol)
	 break;
      if (!env->mip->is_int[i])
	 continue; /* Not an integer variable */
      value = colsol[i];
      if (colsol[i] > env->mip->lb[i] && colsol[i] < env->mip->ub[i]
	  && colsol[i]-floor(colsol[i]) > lpetol &&
	  ceil(colsol[i])-colsol[i] > lpetol){
	 break;
      }
   }

   feasible = i < 0 ? true : false;

   /* step 2. check for the constraint matrix */

   if (feasible){
      rowAct = (double*) calloc(env->mip->m, DSIZE);
      matBeg = env->mip->matbeg;
      matVal = env->mip->matval;
      matInd = env->mip->matind;

      for(i = 0; i < env->mip->n; i++){
	 for(j = matBeg[i]; j<matBeg[i+1]; j++){
	    rowAct[matInd[j]] += matVal[j] * colsol[i];
	 }
      }

      for(i = 0; i < env->mip->m; i++){
	 switch(env->mip->sense[i]){
	  case 'L':
	     if (rowAct[i] > env->mip->rhs[i] + lpetol)
		feasible = FALSE;
	     break;
	  case 'G':
	     if (rowAct[i] < env->mip->rhs[i] - lpetol)
		feasible = FALSE;
	     break;
	  case 'E':
	     if (!((rowAct[i] > env->mip->rhs[i] - lpetol) &&
		   (rowAct[i] < env->mip->rhs[i] + lpetol)))
		feasible = FALSE;
	     break;
	  case 'R':
	     if (rowAct[i] > env->mip->rhs[i] + lpetol ||
		 rowAct[i] < env->mip->rhs[i] - env->mip->rngval[i] - lpetol)
		feasible = FALSE;
	     break;
	  case 'N':
	  default:
	     break;
	 }

	 if (!feasible)
	    break;
      }
   }


   tmp_ind = (int*)malloc(ISIZE*env->mip->n);

   for (i = 0; i < env->mip->n; i++){
      if (colsol[i] > lpetol || colsol[i] < - lpetol){
	 tmp_ind[nz] = i;
	 nz++;
      }
   }

   sol = &(env->best_sol);
   if(sol->xlength){
      FREE(sol->xind);
      FREE(sol->xval);
   }

   sol->xlength = nz;
   sol->objval = 0.0;
   sol->has_sol = FALSE;

   if(nz){
      sol->xval = (double*)calloc(nz,DSIZE);
      sol->xind = (int*)malloc(ISIZE*nz);
      memcpy(sol->xind, tmp_ind, ISIZE*nz);
      for (i = 0; i < nz; i++){
	 sol->xval[i] = colsol[tmp_ind[i]];
	 sol->objval += sol->xval[i] * env->mip->obj[tmp_ind[i]];
      }
   }

   if (feasible){
      /* now, it is feasible, set the best_sol to colsol */
      //FIXME

      if (env->has_ub_estimate){
	 if (env->ub_estimate > sol->objval)
	    env->ub_estimate = sol->objval; //no need for this, I guess.
      }
      else{
	 env->has_ub_estimate = TRUE;
	 env->ub_estimate = sol->objval; //no need for this, I guess.
      }

      if (env->has_ub){
	 if (env->ub > sol->objval)
	    env->ub = sol->objval;
      }
      else{
	 env->has_ub = TRUE;
	 env->ub = sol->objval;
      }
      sol->has_sol = TRUE;
   }else{
      //      env->best_sol.objval = SYM_INFINITY;
      env->best_sol.objval = 0.0;
   }

   if (rowAct){
      FREE(rowAct);
   }

   FREE(tmp_ind);
   env->mip->is_modified = FALSE;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_primal_bound(sym_environment *env, double bound)
{

   if (!env->mip){
      if(env->par.verbosity >= 1){
	 printf("sym_set_primal_bound():There is no loaded mip description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   bound = env->mip->obj_sense == SYM_MINIMIZE ? bound : -bound;

   if (!env->has_ub || bound < env->ub){
      env->has_ub = TRUE;
      env->ub = bound;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_continuous(sym_environment *env, int index)
{
   if (!env->mip || !env->mip->n || index > env->mip->n || index < 0 ||
       !env->mip->is_int){
      if(env->par.verbosity >= 1){
	 printf("sym_set_continuous():There is no loaded mip description or\n");
	 printf("index is out of range or no row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   env->mip->is_int[index] = FALSE;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_integer(sym_environment *env, int index)
{

   if (!env->mip || !env->mip->n || index > env->mip->n || index < 0 ||
       !env->mip->is_int){
      if(env->par.verbosity >= 1){
	 printf("sym_set_integer():There is no loaded mip description or\n");
	 printf("index is out of range or no row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   env->mip->is_int[index] = TRUE;
   env->mip->var_type_modified = TRUE;

   return(FUNCTION_TERMINATED_NORMALLY);
}
/*===========================================================================*/
/*===========================================================================*/

int sym_set_col_names(sym_environment * env, char **colname)
{
   int j;

   if (!env->mip || !env->mip->n || !colname){
      if(env->par.verbosity >= 1){
	 printf("sym_set_col_names():There is no loaded mip description or");
	 printf("an empty name array given!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   if(env->mip->colname){
      for (j = 0; j < env->mip->n; j++){
	 if(env->mip->colname[j]){
	    FREE(env->mip->colname[j]);
	 }
      }
      FREE(env->mip->colname);
   }

   env->mip->colname = (char **)  calloc(sizeof(char *), env->mip->n);

   for (j = 0; j < env->mip->n; j++){
      if(colname[j]){
	 env->mip->colname[j] = (char *) malloc(CSIZE * 21);
	 strncpy(env->mip->colname[j], colname[j], 20);
	 env->mip->colname[j][20] = 0;
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}
/*===========================================================================*/
/*===========================================================================*/

int sym_add_col(sym_environment *env, int numelems, int *indices,
			double *elements, double collb, double colub,
			double obj, char is_int, char *name)
{
   int i, k, n, m, nz, *matBeg, *matInd;
   double *matVal, *colLb, *colUb, *objN, *obj1N, *obj2N;
   char *isInt, **colName;
   int * user_indices, *user_size;

   if ((numelems && !indices) || numelems < 0){
      if(env->par.verbosity >= 1){
	 printf("sym_add_col(): Incorrect column description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   /* order the indices! */
   if(numelems){
      qsort_id(indices, elements, numelems);
   }

   if(!env->mip->n && !env->mip->m){
      n = 1;
      if(numelems){
	 m = indices[numelems - 1];
	 matBeg = (int*) calloc(ISIZE, n+1);
	 matBeg[n] = numelems;
      }
      else{
	 m = 0;
	 matBeg = NULL;
      }
      return(sym_explicit_load_problem(env, n, m, matBeg, indices, elements,
				       &collb, &colub, &is_int, &obj, NULL,
				       NULL, NULL, NULL, TRUE));

   } else{

      n = env->mip->n;
      nz = env->mip->nz;

      user_indices = env->rootdesc->uind.list;
      user_size = &env->rootdesc->uind.size;
      (*user_size) += 1;
      env->rootdesc->uind.list = (int *) malloc(ISIZE*(*user_size));
      memcpy(env->rootdesc->uind.list, user_indices, ISIZE*(*user_size - 1));
      env->rootdesc->uind.list[*user_size - 1] = n;

      colLb = (double*) malloc(DSIZE*(n+1));
      colUb = (double*) malloc(DSIZE*(n+1));
      objN = (double*) malloc(DSIZE*(n+1));
      obj1N = (double*) calloc(DSIZE, (n+1));
      obj2N = (double*) calloc(DSIZE, (n+1));
      isInt = (char*) calloc(CSIZE, (n+1));

      if(n){
	 memcpy(colLb, env->mip->lb, DSIZE*n);
	 memcpy(colUb, env->mip->ub, DSIZE*n);
	 memcpy(objN, env->mip->obj, DSIZE*n);
	 memcpy(obj1N, env->mip->obj1, DSIZE*n);
	 memcpy(obj2N, env->mip->obj2, DSIZE*n);
	 memcpy(isInt, env->mip->is_int, CSIZE*n);
      }

      matBeg = (int*) calloc(ISIZE,(n+2));

      if(numelems){

	 /* if it is out of row size? need additional rows!*/
	 k = indices[numelems-1] + 1 - env->mip->m;
	 if( k > 0){
	    for(i=0; i<k; i++){
	       sym_add_row(env, 0, NULL, NULL, 'N', 0.0, 0.0);
	    }
	 }

	 matInd = (int*) malloc(ISIZE*(nz+numelems));
	 matVal = (double*) malloc(DSIZE*(nz+numelems));

	 if(nz){
	    memcpy(matInd, env->mip->matind, ISIZE*nz);
	    memcpy(matVal, env->mip->matval, DSIZE*nz);
	 }

	 memcpy(matInd + nz, indices, ISIZE*numelems);
	 memcpy(matVal + nz, elements, DSIZE*numelems);

	 FREE(env->mip->matind);
	 FREE(env->mip->matval);
	 env->mip->matind = matInd;
	 env->mip->matval = matVal;
      }

      if(nz){
	 memcpy(matBeg, env->mip->matbeg, ISIZE*(n+1));
      }

      matBeg[n+1] = matBeg[n] + numelems;
      colLb[n] = collb;
      colUb[n] = colub;
      objN[n] = obj;
      isInt[n] = is_int;

      if(n){
	 FREE(env->mip->matbeg);
	 FREE(env->mip->lb);
	 FREE(env->mip->ub);
	 FREE(env->mip->obj);
	 FREE(env->mip->obj1);
	 FREE(env->mip->obj2);
	 FREE(env->mip->is_int);
	 FREE(user_indices);
      }

      env->mip->n = n+1;
      env->mip->nz = nz + numelems;
      env->mip->matbeg = matBeg;
      env->mip->lb =  colLb;
      env->mip->ub = colUb;
      env->mip->obj = objN;
      env->mip->obj1 = obj1N;
      env->mip->obj2 = obj2N;
      env->mip->is_int = isInt;

      /* take care of the name */

      if(env->mip->colname || name){
	 colName = (char**) calloc(sizeof(char*),(n+1));
	 if(env->mip->colname){
	    for (i = 0; i < n; i++){
	       if(env->mip->colname[i]){
		  colName[i] = (char *) malloc(CSIZE * 21);
		  strncpy(colName[i], env->mip->colname[i],21);
		  colName[i][20] = 0;
		  FREE(env->mip->colname[i]);
	       }
	    }
	 }

	 if(name){
	    colName[n] = (char *) malloc(CSIZE * 21);
	    strncpy(colName[n], name, 21);
	    colName[n][20] = 0;
	 }

	 FREE(env->mip->colname);
	 env->mip->colname = colName;
      }
   }

   if (env->mip->change_num){
      if(env->mip->change_type[0] == OBJ_COEFF_CHANGED){ /* will be treated
							    same for now */
	 env->mip->change_type[0] = COLS_ADDED;
      }
      for(i = env->mip->change_num - 1 ; i >=0 ; i--){
	 if (env->mip->change_type[i] == COLS_ADDED){
	    break;
	 }
      }
      if (i < 0 ){
	 env->mip->change_type[env->mip->change_num++] = COLS_ADDED;
      }
   }
   else{
      env->mip->change_type[env->mip->change_num++] = COLS_ADDED;
   }

   env->mip->new_col_num++;
   env->mip->is_modified = TRUE;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_add_row(sym_environment *env, int numelems, int *indices,
		double *elements, char rowsen, double rowrhs,
		double rowrng)
{
   int i, j, k, m, n, nz, *matBeg, *matInd, *lengths;
   double *matVal, *rhs, *range;
   char *sense;

   if ((numelems && !indices) || numelems < 0){
      if(env->par.verbosity >= 1){
	 printf("sym_add_row():Incorrect row description!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   /* order the indices */
   if(numelems){
      qsort_id(indices, elements, numelems);
   }


   if(!env->mip->n && !env->mip->m){
      m = 1;
      if(numelems){
	 n = indices[numelems - 1];
	 matBeg = (int *) calloc(ISIZE, n+1);
	 matInd = (int *)calloc(ISIZE, numelems);
	 for(i = 0, j = 0; i<n; i++){
	    if(j < numelems && indices[j] == i){
	       matBeg[i+1] = matBeg[i] + 1;
	       j++;
	    }else{
	       matBeg[i+1] = matBeg[i];
	    }
	 }
	 if(j!=numelems){
	    printf("sym_add_row(): Unknown Problem!\n");
	    return(FUNCTION_TERMINATED_ABNORMALLY);
	 }
      }else{
	 n = 0;
	 matBeg = NULL;
	 matInd = NULL;
      }
      return(sym_explicit_load_problem(env, n, m, matBeg, matInd, elements,
				       NULL, NULL, NULL, NULL, NULL,
				       &rowsen, &rowrhs, &rowrng, TRUE));
   }else{

      m = env->mip->m;
      nz = env->mip->nz;

      env->base->cutnum +=1;

      if(numelems){

	 /* if it is out of row size? need additional rows!*/
	 k = indices[numelems-1] + 1 - env->mip->n;
	 if(k > 0){
	    for(i = 0; i < k; i++){
	       sym_add_col(env, 0, NULL, NULL, 0.0, SYM_INFINITY, 0.0, FALSE,
			   NULL);
	    }
	    env->mip->is_modified = TRUE;
	 }

	 n = env->mip->n;

	 matBeg = (int*) calloc (n+1, ISIZE);
	 matInd = (int*) malloc(ISIZE*(nz+numelems));
	 matVal = (double*) malloc(DSIZE*(nz+numelems));
	 lengths = (int*) calloc (ISIZE,n);

	 if(env->mip->matbeg){
	    for(i = 0; i<n; i++){
	       lengths[i] = env->mip->matbeg[i+1] - env->mip->matbeg[i];
	    }
	 }

	 for(i = 0; i<numelems; i++){
	    lengths[indices[i]]++;
	 }

	 for(i = 0, j = 0; i<n; i++){
	    matBeg[i+1] = matBeg[i] + lengths[i];
	    if(env->mip->matbeg && env->mip->matind && env->mip->matval){
	       memcpy(matInd + matBeg[i], env->mip->matind +
		      env->mip->matbeg[i],
		      ISIZE * (env->mip->matbeg[i+1]-env->mip->matbeg[i]));
	       memcpy(matVal + matBeg[i], env->mip->matval +
		      env->mip->matbeg[i],
		      DSIZE * (env->mip->matbeg[i+1]-env->mip->matbeg[i]));

	    }
	    if (j < numelems && indices[j] == i){
	       matInd[matBeg[i+1]-1] = m;
	       matVal[matBeg[i+1]-1] = elements[j];
	       j++;
	    }
	 }

	 if(j!=numelems){
	    printf("sym_add_row(): Unknown Problem!\n");
	    return(FUNCTION_TERMINATED_ABNORMALLY);
	 }

	 /*can use FREE_mip_desc???*/
	 FREE(env->mip->matbeg);
	 FREE(env->mip->matind);
	 FREE(env->mip->matval);
	 FREE(lengths);

	 env->mip->nz = nz + numelems;
	 env->mip->matbeg = matBeg;
	 env->mip->matind = matInd;
	 env->mip->matval = matVal;
      }

      sense = (char*) malloc(CSIZE*(m+1));
      rhs = (double*) malloc(DSIZE*(m+1));
      range = (double*) malloc(DSIZE*(m+1));

      if(m){
	 memcpy(sense, env->mip->sense, CSIZE*m);
	 memcpy(range, env->mip->rngval, DSIZE*m);
	 memcpy(rhs, env->mip->rhs, DSIZE*m);
      }

      env->mip->m = m+1;
      sense[m] = rowsen;
      rhs[m] = rowrhs;
      range[m] = rowrng;

      FREE(env->mip->sense);
      FREE(env->mip->rhs);
      FREE(env->mip->rngval);

      env->mip->sense =  sense;
      env->mip->rhs = rhs;
      env->mip->rngval = range;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}
/*===========================================================================*/
/*===========================================================================*/

/* Important: The indices given here are with respect to the current
   not the original user indices! */

int sym_delete_cols(sym_environment *env, int num, int * indices)
{

   int i, j, k, n, nz, num_to_delete = 0, *matBeg, *matInd, *lengths;
   //FIXME! how about base varnum? If they are to be deleted???
   double *matVal, *colLb, *colUb, *objN, *obj1N, *obj2N;
   char *isInt, **colName;

   if (num <= 0){
      return(FUNCTION_TERMINATED_NORMALLY);
   }

   if (!env->mip || !env->mip->n || !env->base || !env->rootdesc ||
       num > env->mip->n || !env->mip->matbeg){
      if(env->par.verbosity >= 1){
	 printf("sym_delete_cols(): No mip description has been loaded\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   int bvarnum = env->base->varnum, bind = 0;
   int user_size = env->rootdesc->uind.size, uind = 0;
   int *bvar_ind = env->base->userind;
   int *user_ind = env->rootdesc->uind.list;

   /* sort the indices in case they are not given sorted! */

   qsort_i(indices, num);

   /* First, adjust the index lists */
   /* Warning: This resets the user indices to be equal to the real indices.
      This is usually fine for generic MIPS, but may not work for applications.
      Names stay the same, however */

   n = env->mip->n;

   for (i = 0, j = 0; i < bvarnum && j < num; i++){
      if (indices[j] == i){
	 j++;
      }else{
	 bvar_ind[bind] = bind;
	 bind++;
      }
   }

   if (j == num){
      for (; i < bvarnum; i++){
	 bvar_ind[bind] = bind;
	 bind++;
      }
      uind = user_size;
   }else{
      for (; i < n && j < num; i++){
	 if (indices[j] == i){
	    j++;
	 }else{
	    user_ind[uind] = uind+bind;
	    uind++;
	 }
      }
      for (; i < n; i++){
	 user_ind[uind] = uind+bind;
	 uind++;
      }
   }

   if (j < num){
      printf("sym_delete_cols() Error: Column index may be out of range.\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

#if 0
   if(i + j != n){
      printf("sym_delete_cols() Error: Unknown problem!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
#endif

   if (bind == bvarnum && uind == user_size){
      printf("sym_delete_cols() Warning: No columns deleted.\n");
      return (FUNCTION_TERMINATED_NORMALLY);
   }

   if (bind < bvarnum){
      env->base->userind = (int *) realloc (bvar_ind, ISIZE * bind);
      env->base->varnum = bind;
   }
   if (uind < user_size){
      env->rootdesc->uind.list = (int *) realloc (user_ind, ISIZE * uind);
      env->rootdesc->uind.size = uind;
   }

   /* Now adjust the MIP description */

   lengths = (int*) malloc (ISIZE*n);

   for (i = 0; i < n; i++){
      lengths[i] = env->mip->matbeg[i+1] - env->mip->matbeg[i];
   }

   nz = env->mip->nz;

   for (i = 0; i < num; i++){
      if (indices[i] < n){
	 num_to_delete += lengths[indices[i]];
      }else{
	 /*FIXME*/
	 printf("sym_delete_cols(): Error. Column index is out of range!\n");
	 return(FUNCTION_TERMINATED_ABNORMALLY);
      }
   }

   matBeg =  env->mip->matbeg;
   matInd =  env->mip->matind;
   matVal =  env->mip->matval;
   colLb =   env->mip->lb;
   colUb =   env->mip->ub;
   objN =    env->mip->obj;
   obj1N =   env->mip->obj1;
   obj2N =   env->mip->obj2;
   isInt =   env->mip->is_int;
   colName = env->mip->colname;

   matBeg[0] = 0;

   for(i = 0, j = 0, k = 0; j < num; i++){
      if (indices[j] == i){
	 j++;
	 continue;
      }
      matBeg[k+1] = matBeg[k] + lengths[i];
      memmove(matInd + matBeg[k], matInd + matBeg[i], ISIZE * lengths[i]);
      memmove(matVal + matBeg[k], matVal + matBeg[i], DSIZE * lengths[i]);
      colLb[k] = colLb[i];
      colUb[k] = colUb[i];
      objN[k] = objN[i];
      isInt[k] = isInt[i];
      colName[k] = colName[i];
      k++;
   }

   for(; i < n; i++, k++){
      matBeg[k+1] = matBeg[k] + lengths[i];
      memmove(matInd + matBeg[k], matInd + matBeg[i], ISIZE * lengths[i]);
      memmove(matVal + matBeg[k], matVal + matBeg[i], DSIZE * lengths[i]);
      colLb[k] = colLb[i];
      colUb[k] = colUb[i];
      objN[k] = objN[i];
      isInt[k] = isInt[i];
      colName[k] = colName[i];
   }

   if (obj1N){
      for(i = 0, j = 0, k = 0; j < num; i++){
	 if (indices[j] == i){
	    j++;
	    continue;
	 }
	 obj1N[k] = obj1N[i];
	 k++;
      }

      for(; i < n; i++, k++){
	 obj1N[k] = obj1N[i];
      }
   }

   if (obj2N){
      for(i = 0, j = 0, k = 0; j < num; i++){
	 if (indices[j] == i){
	    j++;
	    continue;
	 }
	 obj2N[k] = obj2N[i];
	 k++;
      }

      for(; i < n; i++, k++){
	 obj2N[k] = obj2N[i];
      }
   }

   n = env->mip->n = n - num;
   nz = env->mip->nz = nz - num_to_delete;
   env->mip->matbeg = (int *) realloc(matBeg, (n+1)*ISIZE);
   env->mip->matind = (int *) realloc(matInd, nz*ISIZE);
   env->mip->matval = (double *) realloc(matVal, nz*DSIZE);
   env->mip->lb = (double *) realloc(colLb, n*DSIZE);
   env->mip->ub = (double *) realloc(colUb, n*DSIZE);
   env->mip->obj = (double *) realloc(objN, n*DSIZE);
   env->mip->is_int = (char *) realloc(isInt, n*CSIZE);
   env->mip->colname = (char **) realloc(colName, n*sizeof(char *));

   free(lengths);

   env->mip->is_modified = TRUE;

   return(FUNCTION_TERMINATED_NORMALLY);

}

/*===========================================================================*/
/*===========================================================================*/
int sym_delete_rows(sym_environment *env, int num, int * indices)
{

   int i, j, k, n, m, nz, new_num_elements = 0, new_num_rows = 0;
   int *matBeg, *matInd, *new_rows = 0;
   double *matVal, *rhs, *range;
   char *sense;

   if (num <= 0){
      return(FUNCTION_TERMINATED_NORMALLY);
   }

   if (!env->mip || !env->mip->m || !env->base || num > env->mip->m){
      if(env->par.verbosity >= 1){
	 printf("sym_delete_rows():There is no loaded mip or base description \n");
	 printf("or the number of rows or num exceeds the real row number!\n");
      }
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   //FIXME!
   env->base->cutnum -= num;

   if (!env->mip->matbeg){
      /* We don't have a generic MIP description */
      return (FUNCTION_TERMINATED_NORMALLY);
   }

   n = env->mip->n;
   m = env->mip->m;
   nz = env->mip->nz;

   matBeg = env->mip->matbeg;
   matInd = env->mip->matind;
   matVal = env->mip->matval;
   sense = env->mip->sense;
   rhs = env->mip->rhs;
   range = env->mip->rngval;

   /* sort the indices in case they are not given sorted! */

   qsort_i(indices, num);

   new_rows = (int *) malloc(m*ISIZE);
   for (new_num_rows = 0, i = 0, k = 0; i < m && k < num; i++){
      if (indices[k] == i){
	 new_rows[i] = -1;
	 k++;
      }else{
	 new_rows[i] = new_num_rows++;
      }
   }

   for (; i < m; i++){
      new_rows[i] = new_num_rows++;
   }

   if (k < num){
      printf("sym_delete_rows() Error: Row index may be out of range.\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   for (new_num_elements = 0, i = 0, j = 0; i < n; i++){
      for (; j < matBeg[i+1]; j++){
	 if (new_rows[matInd[j]] < 0){
	    continue;
	 }else{
	    matInd[new_num_elements] = new_rows[matInd[j]];
	    matVal[new_num_elements++] = matVal[j];
	 }
      }
      j = matBeg[i+1];
      matBeg[i+1] = new_num_elements;
   }
   //   matBeg[n] = new_num_elements;

   for (i = 0; i < m; i++){
      if (new_rows[i] >= 0){
	 sense[new_rows[i]] = sense[i];
	 rhs[new_rows[i]] = rhs[i];
	 range[new_rows[i]] = range[i];
      }
   }

   if (new_num_rows != m - num){
      printf("sym_delete_rows(): Unknown error!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   env->mip->m  = new_num_rows;
   env->mip->nz = new_num_elements;

   env->mip->rhs    = (double *) realloc(rhs, DSIZE * new_num_rows);
   env->mip->sense  = (char *)   realloc(sense, CSIZE * new_num_rows);
   env->mip->rngval = (double *) realloc(range, DSIZE * new_num_rows);

   env->mip->matval = (double *) realloc(matVal, DSIZE*new_num_elements);
   env->mip->matind = (int *)    realloc(matInd, ISIZE*new_num_elements);

   FREE(new_rows);

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_write_warm_start_desc(warm_start_desc *ws, char *file)
{

   FILE * f = NULL;
   int i, j;
   cut_data ** cuts;
   problem_stat stat;
   node_times compT;

   f = fopen(file, "w");

   if (!ws){
      printf("There is no loaded warmStart to write!\n");
      fclose(f);
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   else{
      fprintf(f, "########################################################\n");
      fprintf(f, " BOUND INFO \n");
      fprintf(f, "########################################################\n");
      fprintf(f, " PHASE      : %i\n", ws->phase);
      fprintf(f, " LB         : %.4f\n", ws->lb);
      fprintf(f, " HAS_UB     : %i\n", (int)ws->has_ub);
      fprintf(f, " UB         : %.4f\n\n", ws->ub);

      fprintf(f, "########################################################\n");
      fprintf(f, " CUT INFO \n");
      fprintf(f, "########################################################\n");
      fprintf(f, " CUT_NUM             : %i\n", ws->cut_num);
      fprintf(f, " ALLOCATED_CUT_NUM   : %i\n\n",
	      ws->allocated_cut_num);

      //FIXME! WHAT TYPE A CUT CAN BE OTHER THAN EXPLICIT_ROW

      cuts = ws->cuts;

      for(i=0; i<ws->cut_num; i++){
	 fprintf(f, " CUT %i : \n",i);
	 fprintf(f, " SIZE        : %i \n",cuts[i]->size);
	 fprintf(f, " ELEMENTS    : ");
	 for(j=0; j<cuts[i]->size; j++){
	    fprintf(f," %i",(int)cuts[i]->coef[j]);
	 }
	 fprintf(f, "\n");
	 fprintf(f, " RHS         : %.4f \n",cuts[i]->rhs);
	 fprintf(f, " RANGE       : %.4f \n",cuts[i]->range);
	 fprintf(f, " TYPE        : %i \n",(int)cuts[i]->type);
	 fprintf(f, " SENSE       : %c \n",cuts[i]->sense);
	 fprintf(f, " DELETABLE   : %i \n",(int)cuts[i]->deletable);
	 fprintf(f, " BRANCH      : %i \n",(int)cuts[i]->branch);
	 fprintf(f, " NAME        : %i \n\n",cuts[i]->name);
      }

      fprintf(f, "########################################################\n");
      fprintf(f, " PROBLEM STATISTICS \n");
      fprintf(f, "########################################################\n");

      stat= ws->stat;


      fprintf(f," ROOT_LB                : %.4f\n", stat.root_lb);
      fprintf(f," CUTS_IN_POOL           : %i\n", stat.cuts_in_pool);
      fprintf(f," MAXIMIM_DEPTH          : %i\n", stat.max_depth);
      fprintf(f," DIVING_CHAINS          : %i\n", stat.chains);
      fprintf(f," DIVING_STOPS           : %i\n", stat.diving_halts);
      fprintf(f," TREE_SIZE              : %i\n", stat.tree_size);
      fprintf(f," CREATED_NODES          : %i\n", stat.created);
      fprintf(f," ANALYZED_NODES         : %i\n", stat.analyzed);
      fprintf(f," LEAVES_BEFORE_TRIMMING : %i\n", stat.leaves_before_trimming);
      fprintf(f," LEAVES_BEFORE_TRIMMING : %i\n", stat.leaves_after_trimming);
      fprintf(f," NOT_FIXED_VARIABLE_NUM : %i\n", stat.vars_not_priced);
      fprintf(f," NF_STATUS_OF_ROOT      : %i\n\n", (int)stat.nf_status);

      fprintf(f, "########################################################\n");
      fprintf(f, " COMPUTATION TIMES \n");
      fprintf(f, "########################################################\n");

      compT = ws->comp_times;

      fprintf(f," COMMUNICATION       : %.4f\n",compT.communication);
      fprintf(f," LP                  : %.4f\n",compT.lp);
      fprintf(f," SEPARATION          : %.4f\n",compT.separation);
      fprintf(f," FIXING              : %.4f\n",compT.fixing);
      fprintf(f," PRICING             : %.4f\n",compT.pricing);
      fprintf(f," STRONG_BRANCHING    : %.4f\n",compT.strong_branching);
      fprintf(f," WALL_CLOCK_LP       : %.4f\n",compT.wall_clock_lp);
      fprintf(f," RAMP_UP_TM          : %.4f\n",compT.ramp_up_tm);
      fprintf(f," RAMP_UP_LP          : %.4f\n",compT.ramp_up_lp);
      fprintf(f," RAMP_DOWN_TIME      : %.4f\n",compT.ramp_down_time);
      fprintf(f," IDLE_DIVING         : %.4f\n",compT.idle_diving);
      fprintf(f," IDLE_NODE           : %.4f\n",compT.idle_node);
      fprintf(f," IDLE_NAMES          : %.4f\n",compT.idle_names);
      fprintf(f," IDLE_CUTS           : %.4f\n",compT.idle_cuts);
      fprintf(f," START_NODE          : %.4f\n",compT.start_node);
      fprintf(f," CUT_POOL            : %.4f\n\n",compT.cut_pool);

      fprintf(f, "########################################################\n");
      fprintf(f, " TREE DESCRIPTION \n");
      fprintf(f, "########################################################\n");

      write_tree(ws->rootnode, f);
      fclose(f);
      return(FUNCTION_TERMINATED_NORMALLY);
   }
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc * sym_read_warm_start(char *file)
{
   FILE * f;
   char str[80];
   int i=0, j=0, num=0, ch=0;
   int temp =0;
   cut_data *cut;
   problem_stat stat;
   node_times compT;
   warm_start_desc *ws;

   if (!(f = fopen(file, "r"))){
      printf("sym_read_warm_start():");
      printf("Can not open the warm start file to read!\n");
      return(NULL);
   }
   else{

      ws = (warm_start_desc *) calloc(1, sizeof(warm_start_desc));

      /* bound info */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %i", str, str, &ws->phase);
      fscanf(f,"%s %s %lf", str, str, &ws->lb);
      fscanf(f,"%s %s %i", str, str, &ch);
      ws->has_ub = (char)ch;
      fscanf(f,"%s %s %lf", str, str, &ws->ub);

      /* cut info */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %i", str, str, &ws->cut_num);
      fscanf(f,"%s %s %i", str, str, &temp);
      ws->allocated_cut_num = temp;

      if (temp){
	 ws->cuts = (cut_data **)malloc(temp *sizeof(cut_data *));
	 for(i = 0; i < ws->cut_num; i++){
	    cut = (cut_data*)malloc(sizeof(cut_data));
	    fscanf(f,"%s %i %s", str, &num, str);
	    fscanf(f,"%s %s %i", str, str, &cut->size);
	    cut->coef = (char*)malloc(CSIZE*cut->size);
	    fscanf(f,"%s %s", str, str);

	    for(j=0; j<cut->size; j++){
	       fscanf(f,"%i", &ch);
	       cut->coef[j] = (char)ch;
	    }
	    fscanf(f,"%s %s %lf", str, str, &cut->rhs);
	    fscanf(f,"%s %s %lf", str, str, &cut->range);
	    fscanf(f,"%s %s %i", str, str, &ch);
	    cut->type = (char)ch;
	    fscanf(f,"%s %s %c", str, str, &cut->sense);
	    fscanf(f,"%s %s %i", str, str, &ch);
	    cut->deletable=(char)ch;
	    fscanf(f,"%s %s %i", str, str, &ch);
	    cut->branch = (char)ch;
	    fscanf(f,"%s %s %i", str, str, &cut->name);

	    ws->cuts[i] = cut;
	 }
      }

      /* problem stats */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %lf", str, str, &stat.root_lb);
      fscanf(f,"%s %s %i", str, str, &stat.cuts_in_pool);
      fscanf(f,"%s %s %i", str, str, &stat.max_depth);
      fscanf(f,"%s %s %i", str, str, &stat.chains);
      fscanf(f,"%s %s %i", str, str, &stat.diving_halts);
      fscanf(f,"%s %s %i", str, str, &stat.tree_size);
      fscanf(f,"%s %s %i", str, str, &stat.created);
      fscanf(f,"%s %s %i", str, str, &stat.analyzed);
      fscanf(f,"%s %s %i", str, str, &stat.leaves_before_trimming);
      fscanf(f,"%s %s %i", str, str, &stat.leaves_after_trimming);
      fscanf(f,"%s %s %i", str, str, &stat.vars_not_priced);
      fscanf(f,"%s %s %i", str, str, &ch);
      stat.nf_status = (char)ch;

      ws->stat = stat;

      /* computation times */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      fscanf(f,"%s %s %lf", str, str, &compT.communication);
      fscanf(f,"%s %s %lf", str, str, &compT.lp);
      fscanf(f,"%s %s %lf", str, str, &compT.separation);
      fscanf(f,"%s %s %lf", str, str, &compT.fixing);
      fscanf(f,"%s %s %lf", str, str, &compT.pricing);
      fscanf(f,"%s %s %lf", str, str, &compT.strong_branching);
      fscanf(f,"%s %s %lf", str, str, &compT.wall_clock_lp);
      fscanf(f,"%s %s %lf", str, str, &compT.ramp_up_tm);
      fscanf(f,"%s %s %lf", str, str, &compT.ramp_up_lp);
      fscanf(f,"%s %s %lf", str, str, &compT.ramp_down_time);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_diving);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_node);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_names);
      fscanf(f,"%s %s %lf", str, str, &compT.idle_cuts);
      fscanf(f,"%s %s %lf", str, str, &compT.start_node);
      fscanf(f,"%s %s %lf", str, str, &compT.cut_pool);

      ws->comp_times = compT;

      /* tree description */
      fscanf(f,"%s %s %s %s", str, str, str, str);
      ws->rootnode = (bc_node*)calloc(1,sizeof(bc_node));
      read_tree(ws->rootnode, f);
   }

   fclose(f);
   return(ws);
}

/*===========================================================================*/
/*===========================================================================*/

void sym_delete_warm_start(warm_start_desc *ws)
{
   int i, temp;
   if (ws) {
      if (ws->rootnode) {
	 free_subtree(ws->rootnode);
      }
      if (ws->cuts){
	 temp = ws->cut_num;
	 for(i = 0; i < ws->cut_num; i++){
	    if (ws->cuts[i]){
	       if (ws->cuts[i]->coef){
		  FREE(ws->cuts[i]->coef);
	       }
	    }
	    FREE(ws->cuts[i]);
	 }
	 FREE(ws->cuts);
      }

      if(ws->best_sol.xlength){
	 FREE(ws->best_sol.xind);
	 FREE(ws->best_sol.xval);
      }
      FREE(ws);
   }

   ws = 0;
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc * sym_get_warm_start(sym_environment *env, int copy_warm_start)
{

   warm_start_desc * ws;

   if (!env->warm_start){
      printf("sym_get_warm_start_desc():");
      printf("The env. warm start description is empty!\n");
      return(NULL);
   }

   if (copy_warm_start){
      ws = create_copy_warm_start(env->warm_start);
   }
   else{
      ws = env->warm_start;
      env->warm_start = 0;
   }

   return(ws);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_warm_start (sym_environment *env, warm_start_desc *ws)
{

   if (!ws){
      printf("sym_set_warm_start():The warm_start desc. is empty!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }

   warm_start_desc * ws_copy = create_copy_warm_start(ws);
   sym_delete_warm_start(env->warm_start);
   env->warm_start = ws_copy;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_int_param(sym_environment *env, const char *key, int value)
{
   int termcode;
   char *line = (char*)malloc(CSIZE*(MAX_LINE_LENGTH+1));
   sprintf(line, "%s %d", key, value);
   termcode = set_param(env, line);
   FREE(line);

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_dbl_param(sym_environment *env, const char *key, double value)
{
   int termcode;
   char *line = (char*)malloc(CSIZE*(MAX_LINE_LENGTH+1));
   sprintf(line, "%s %.30f", key, value);
   termcode = set_param(env, line);
   FREE(line);

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_set_str_param(sym_environment *env, const char *key, const char *value)
{
   int termcode;
   char *line = (char*)malloc(CSIZE*(MAX_LINE_LENGTH+1));
   sprintf(line, "%s %s", key, value);
   termcode = set_param(env, line);
   FREE(line);

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_int_param(sym_environment *env, const char *key, int *value)
{

   tm_params *tm_par = &env->par.tm_par;
   lp_params *lp_par = &env->par.lp_par;
   cg_params *cg_par = &env->par.cg_par;
   cp_params *cp_par = &env->par.cp_par;

   dg_params *dg_par = &env->par.dg_par;

   if (strcmp(key, "verbosity") == 0){
      *value = env->par.verbosity;
      return(0);
   }
   else if (strcmp(key, "random_seed") == 0){
      *value = env->par.random_seed;
      return(0);
   }

   /***********************************************************************
    ***                    Master params                            ***
    ***********************************************************************/
   else if (strcmp(key, "M_verbosity") == 0){
      *value = env->par.verbosity;
      return(0);
   }
   else if (strcmp(key, "M_random_seed") == 0){
      *value = env->par.random_seed;
      return(0);
   }
   else if (strcmp(key, "tm_debug") == 0 ||
	    strcmp(key, "M_tm_debug") == 0){
      *value = env->par.tm_debug;
      return(0);
   }
   else if (strcmp(key, "dg_debug") == 0 ||
	    strcmp(key, "M_dg_debug") == 0){
      *value = env->par.dg_debug;
      return(0);
   }
   else if (strcmp(key, "pvm_trace") == 0 ||
	    strcmp(key, "M_pvm_trace") == 0){
      *value = env->par.pvm_trace;
      return(0);
   }
   else if (strcmp(key, "do_branch_and_cut") == 0 ||
	    strcmp(key, "M_do_branch_and_cut") == 0){
      *value = env->par.do_branch_and_cut;
      return(0);
   }
   else if (strcmp(key, "do_draw_graph") == 0 ||
	    strcmp(key, "M_do_draw_graph") == 0){
      *value = env->par.do_draw_graph;
      return(0);
   }
   else if (strcmp(key, "use_permanent_cut_pools") == 0 ||
	    strcmp(key, "M_use_permanent_cut_pools") == 0){
      *value = env->par.use_permanent_cut_pools;
      return(0);
   }
   else if (strcmp(key, "mc_search_order") == 0 ||
	    strcmp(key, "M_mc_search_order") == 0){
      *value = env->par.mc_search_order;
      return(0);
   }
   else if (strcmp(key, "mc_warm_start") == 0 ||
	     strcmp(key, "M_mc_warm_start") == 0){
      *value = env->par.mc_warm_start;
      return(0);
   }
   else if (strcmp(key, "mc_warm_start_rule") == 0 ||
	    strcmp(key, "M_mc_warm_start_rule") == 0){
      *value = env->par.mc_warm_start_rule;
      return(0);
   }
   else if (strcmp(key, "trim_warm_tree") == 0 ||
	     strcmp(key, "M_trim_warm_tree") == 0){
      *value = env->par.trim_warm_tree;
      return(0);
   }

   /***********************************************************************
    ***                 DrawGraph params                            ***
    ***********************************************************************/
   else if (strcmp(key, "echo_commands") == 0 ||
	    strcmp(key, "DG_echo_commands") == 0){
      *value = dg_par->echo_commands;
      return(0);
   }
   else if (strcmp(key, "canvas_width") == 0 ||
	    strcmp(key, "DG_canvas_width") == 0){
      *value = dg_par->canvas_width;
      return(0);
   }
   else if (strcmp(key, "canvas_height") == 0 ||
	    strcmp(key, "DG_canvas_height") == 0){
      *value = dg_par->canvas_height;
      return(0);
   }
   else if (strcmp(key, "viewable_width") == 0 ||
	    strcmp(key, "DG_viewable_width") == 0){
      *value = dg_par->viewable_width;
      return(0);
   }
   else if (strcmp(key, "viewable_height") == 0 ||
	    strcmp(key, "DG_viewable_height") == 0){
      *value = dg_par->viewable_width;
      return(0);
   }
   else if (strcmp(key, "disp_nodelabels") == 0 ||
	    strcmp(key, "DG_disp_nodelabels") == 0){
      *value = dg_par->disp_nodelabels;
      return(0);
   }
   else if (strcmp(key, "disp_nodeweights") == 0 ||
	    strcmp(key, "DG_disp_nodeweights") == 0){
      *value = dg_par->disp_nodeweights;
      return(0);
   }
   else if (strcmp(key, "disp_edgeweights") == 0 ||
	    strcmp(key, "DG_disp_edgeweights") == 0){
      *value = dg_par->disp_edgeweights;
      return(0);
   }
   else if (strcmp(key, "node_radius") == 0 ||
	    strcmp(key, "DG_node_radius") == 0){
      *value = dg_par->node_radius;
      return(0);
   }
   else if (strcmp(key, "interactive_mode") == 0 ||
	    strcmp(key, "DG_interactive_mode") == 0){
      *value = dg_par->interactive_mode;
      return(0);
   }
   else if (strcmp(key, "mouse_tracking") == 0 ||
	    strcmp(key, "DG_mouse_tracking") == 0){
      *value = dg_par->mouse_tracking;
      return(0);
   }

   /***********************************************************************
    ***                  Treemanager params                         ***
    ***********************************************************************/

   if (strcmp(key, "TM_verbosity") == 0){
      *value = tm_par->verbosity;
      return(0);
   }
   else if (strcmp(key, "lp_debug") == 0 ||
	    strcmp(key, "TM_lp_debug") == 0){
      *value = tm_par->lp_debug;
      return(0);
   }
   else if (strcmp(key, "cg_debug") == 0 ||
	    strcmp(key, "TM_cg_debug") == 0){
      *value = tm_par->cg_debug;
      return(0);
   }
   else if (strcmp(key, "cp_debug") == 0 ||
	    strcmp(key, "TM_cp_debug") == 0){
      *value = tm_par->cp_debug;
      return(0);
   }
   else if (strcmp(key, "max_active_nodes") == 0 ||
	    strcmp(key, "TM_max_active_nodes") == 0){
      *value = tm_par->max_active_nodes;
      return(0);
   }
   else if (strcmp(key, "max_cp_num") == 0 ||
	    strcmp(key, "TM_max_cp_num") == 0){
      *value = tm_par->max_cp_num;
      return(0);
   }
   else if (strcmp(key, "lp_mach_num") == 0 ||
	    strcmp(key, "TM_lp_mach_num") == 0){
      *value = tm_par->lp_mach_num;
      return(0);
   }
   else if (strcmp(key, "cg_mach_num") == 0 ||
	    strcmp(key, "TM_cg_mach_num") == 0){
      *value = tm_par->cg_mach_num;
      return(0);
   }
   else if (strcmp(key, "cp_mach_num") == 0 ||
	    strcmp(key, "TM_cp_mach_num") == 0){
      *value = tm_par->cp_mach_num;
      return(0);
   }
#ifndef COMPILE_IN_CG
   else if (strcmp(key, "use_cg") == 0 ||
	    strcmp(key, "TM_use_cg") == 0 ||
	    strcmp(key, "LP_use_cg") == 0){
      *value = tm_par->use_cg;
      return(0);
   }
#endif
   else if (strcmp(key, "TM_random_seed") == 0){
      *value = tm_par->random_seed;
      return(0);
   }
   else if (strcmp(key, "diving_strategy") == 0 ||
	    strcmp(key, "TM_diving_strategy") == 0){
      *value = tm_par->diving_strategy;
      return(0);
   }
   else if (strcmp(key, "diving_k") == 0 ||
	    strcmp(key, "TM_diving_k") == 0){
      *value = tm_par->diving_k;
      return(0);
   }
   else if (strcmp(key, "node_selection_rule") == 0 ||
	    strcmp(key, "TM_node_selection_rule") == 0){
      *value = tm_par->node_selection_rule;
      return(0);
   }
   else if (strcmp(key, "keep_description_of_pruned") == 0 ||
	    strcmp(key, "TM_keep_description_of_pruned") == 0){
      *value = tm_par->keep_description_of_pruned;
      return(0);
   }
   else if (strcmp(key, "keep_warm_start") == 0){
      if (tm_par->keep_description_of_pruned == KEEP_IN_MEMORY){
	 *value = TRUE;
      }else{
	 *value = FALSE;
      }
      return(0);
   }
   else if (strcmp(key, "warm_start") == 0 ||
	    strcmp(key, "TM_warm_start") == 0){
      *value = tm_par->warm_start;
      return(0);
   }
   else if (strcmp(key, "vbc_emulation") == 0 ||
	    strcmp(key, "TM_vbc_emulation") == 0){
      *value = tm_par->vbc_emulation;
      return(0);
   }
   else if (strcmp(key, "logging_interval") == 0 ||
	    strcmp(key, "TM_logging_interval") == 0){
      *value = tm_par->logging_interval;
      return(0);
   }
   else if (strcmp(key, "logging") == 0 ||
	    strcmp(key, "TM_logging") == 0){
      *value = tm_par->logging;
      return(0);
   }
   else if (strcmp(key, "price_in_root") == 0 ||
	    strcmp(key, "TM_price_in_root") == 0){
      *value = tm_par->price_in_root;
      return(0);
   }
   else if (strcmp(key, "trim_search_tree") == 0 ||
	    strcmp(key, "TM_trim_search_tree") == 0){
      *value = tm_par->trim_search_tree;
      return(0);
   }
   else if (strcmp(key, "colgen_in_first_phase") == 0 ||
	    strcmp(key, "TM_colgen_in_first_phase") == 0){
      *value = tm_par->colgen_strat[0];
      return(0);
   }

   else if (strcmp(key, "colgen_in_second_phase") == 0 ||
	    strcmp(key, "TM_colgen_in_second_phase") == 0){
      *value = tm_par->colgen_strat[1];
      return(0);
   }
   else if (strcmp(key, "node_limit") == 0 ||
	    strcmp(key, "TM_node_limit") == 0){
      *value = tm_par->node_limit;
      return(0);
   }
   else if (strcmp(key, "find_first_feasible") == 0 ||
	    strcmp(key, "TM_find_first_feasible") == 0){
      *value = tm_par->find_first_feasible;
      return(0);
   }
   else if (strcmp(key, "sensitivity_analysis") == 0 ||
	    strcmp(key, "TM_sensitivity_analysis") == 0 ){
      *value = tm_par->sensitivity_analysis;
      return(0);
   }

   /***********************************************************************
    ***                      LP params                              ***
    ***********************************************************************/
   if (strcmp(key, "LP_verbosity") == 0){
      *value = lp_par->verbosity;
   }
   else if (strcmp(key, "set_obj_upper_lim") == 0 ||
	    strcmp(key, "LP_set_obj_upper_lim") == 0){
      *value = lp_par->set_obj_upper_lim;
      return(0);
   }
   else if (strcmp(key, "do_primal_heuristic") == 0 ||
	    strcmp(key, "LP_do_primal_heuristic") == 0){
      *value = lp_par->do_primal_heuristic;
      return(0);
   }
   else if (strcmp(key, "scaling") == 0 ||
	    strcmp(key, "LP_scaling") == 0){
      *value = lp_par->scaling;
      return(0);
   }
   else if (strcmp(key, "fastmip") == 0 ||
	    strcmp(key, "LP_fastmip") == 0){
      *value = lp_par->fastmip;
      return(0);
   }
   else if (strcmp(key, "should_warmstart_chain") == 0 ||
	    strcmp(key, "LP_should_warmstart_chain") == 0){
      *value = lp_par->should_warmstart_chain;
      return(0);
   }
   else if (strcmp(key, "should_reuse_lp") == 0 ||
	    strcmp(key, "LP_should_reuse_lp") == 0){
      *value = lp_par->should_reuse_lp;
      return(0);
   }
   else if (strcmp(key, "try_to_recover_from_error") == 0 ||
	    strcmp(key, "LP_try_to_recover_from_error") == 0){
      *value = lp_par->try_to_recover_from_error;
      return(0);
   }
   else if (strcmp(key, "problem_type") == 0 ||
	    strcmp(key, "LP_problem_type") == 0){
      *value = lp_par->problem_type;
      return(0);
   }
   else if (strcmp(key, "not_fixed_storage_size") == 0 ||
	    strcmp(key, "LP_not_fixed_storage_size") == 0 ||
	    strcmp(key, "TM_not_fixed_storage_size") == 0 ){
      *value = lp_par->not_fixed_storage_size;
      return(0);
   }
   else if (strcmp(key, "cut_pool_check_frequency") == 0 ||
	    strcmp(key, "LP_cut_pool_check_frequency") == 0){
      *value = lp_par->cut_pool_check_freq;
      return(0);
   }
   else if (strcmp(key, "load_balance_level") == 0 ||
	    strcmp(key, "LP_load_balance_level") == 0){
      *value = lp_par->load_balance_level;
      return(0);
   }
   else if (strcmp(key, "load_balance_iterations") == 0 ||
	    strcmp(key, "LP_load_balance_iterations") == 0){
      *value = lp_par->load_balance_iterations;
      return(0);
   }
   else if (strcmp(key, "load_balance_compare_candidates") == 0 ||
	    strcmp(key, "LP_load_balance_compare_candidates") == 0){
      *value = lp_par->load_balance_compare_candidates;
      return(0);
   }
   else if (strcmp(key, "fractional_diving_num") == 0 ||
	    strcmp(key, "LP_fractional_diving_num") == 0){
      *value = lp_par->fractional_diving_num;
      return(0);
   }
   else if (strcmp(key, "max_cols_to_add_min") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_min") == 0){
      *value = lp_par->max_non_dual_feas_to_add_min;
      return(0);
   }
   else if (strcmp(key, "max_non_dual_feas_to_add_max") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_max") == 0){
      *value = lp_par->max_non_dual_feas_to_add_max;
      return(0);
   }
   else if (strcmp(key, "max_not_fixable_to_add_min") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_min") == 0){
      *value = lp_par->max_not_fixable_to_add_min;
      return(0);
   }
   else if (strcmp(key, "max_not_fixable_to_add_max") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_max") == 0){
      *value = lp_par->max_not_fixable_to_add_max;
      return(0);
   }

   else if (strcmp(key, "mat_col_compress_num") == 0 ||
	    strcmp(key, "LP_mat_col_compress_num") == 0){
      *value = lp_par->mat_col_compress_num;
      return(0);
   }
   else if (strcmp(key, "mat_row_compress_num") == 0 ||
	    strcmp(key, "LP_mat_row_compress_num") == 0){
      *value = lp_par->mat_row_compress_num;
      return(0);
   }
   else if (strcmp(key, "tailoff_gap_backsteps") == 0 ||
	    strcmp(key, "LP_tailoff_gap_backsteps") == 0){
      *value = lp_par->tailoff_gap_backsteps;
      return(0);
   }
   else if (strcmp(key, "tailoff_obj_backsteps") == 0 ||
	    strcmp(key, "LP_tailoff_obj_backsteps") == 0){
      *value = lp_par->tailoff_obj_backsteps;
      return(0);
   }
   else if (strcmp(key, "ineff_cnt_to_delete") == 0 ||
	    strcmp(key, "LP_ineff_cnt_to_delete") == 0){
      *value = lp_par->ineff_cnt_to_delete;
      return(0);
   }
   else if (strcmp(key, "eff_cnt_before_cutpool") == 0 ||
	    strcmp(key, "LP_eff_cnt_before_cutpool") == 0){
      *value = lp_par->eff_cnt_before_cutpool;
      return(0);
   }
   else if (strcmp(key, "ineffective_constraints") == 0 ||
	    strcmp(key, "LP_ineffective_constraints") == 0){
      *value = lp_par->ineffective_constraints;
      return(0);
   }
   else if (strcmp(key, "base_constraints_always_effective") == 0 ||
	    strcmp(key, "LP_base_constraints_always_effective") == 0){
      *value = lp_par->base_constraints_always_effective;
      return(0);
   }
   else if (strcmp(key, "branch_on_cuts") == 0 ||
	    strcmp(key, "LP_branch_on_cuts") == 0){
      *value = lp_par->branch_on_cuts;
      return(0);
   }
   else if (strcmp(key, "discard_slack_cuts") == 0 ||
	    strcmp(key, "LP_discard_slack_cuts") == 0){
      *value = lp_par->discard_slack_cuts;
      return(0);
   }
   else if (strcmp(key, "max_cut_num_per_iter") == 0 ||
	    strcmp(key, "LP_max_cut_num_per_iter") == 0){
      *value = lp_par->max_cut_num_per_iter;
      return(0);
   }
   else if (strcmp(key, "max_cut_num_per_iter_root") == 0 ||
	    strcmp(key, "LP_max_cut_num_per_iter_root") == 0){
      *value = lp_par->max_cut_num_per_iter_root;
      return(0);
   }
   else if (strcmp(key, "min_root_cut_rounds") == 0) {
      *value = lp_par->min_root_cut_rounds;
      return(0);
   }
   else if (strcmp(key, "max_cut_length") == 0) {
      *value = lp_par->max_cut_length;
      return(0);
   }

   /* variable fixing params */
   else if (strcmp(key, "do_reduced_cost_fixing") == 0 ||
	    strcmp(key, "LP_do_reduced_cost_fixing") == 0){
      *value = lp_par->do_reduced_cost_fixing;
      return(0);
   }
   else if (strcmp(key, "do_logical_fixing") == 0 ||
	    strcmp(key, "LP_do_logical_fixing") == 0){
      *value = lp_par->do_logical_fixing;
      return(0);
   }
   else if (strcmp(key, "fixed_to_ub_before_logical_fixing") == 0 ||
	    strcmp(key, "LP_fixed_to_ub_before_logical_fixing") == 0){
      *value = lp_par->fixed_to_ub_before_logical_fixing;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_cuts") == 0){
      *value = cg_par->do_findcuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_gomory_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_gomory_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_gomory_cuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_redsplit_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_redsplit_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_redsplit_cuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_knapsack_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_knapsack_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_knapsack_cuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_oddhole_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_oddhole_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_oddhole_cuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_probing_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_probing_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_probing_cuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_clique_cuts") == 0 ||
            strcmp(key, "LP_generate_cgl_clique_cuts") == 0){
     *value = lp_par->cgl.generate_cgl_clique_cuts;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_mir_cuts") == 0 ||
            strcmp(key, "LP_generate_cgl_mir_cuts") == 0){
     *value = lp_par->cgl.generate_cgl_mir_cuts;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_twomir_cuts") == 0 ||
            strcmp(key, "LP_generate_cgl_twomir_cuts") == 0){
     *value = lp_par->cgl.generate_cgl_twomir_cuts;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_flowcover_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_flowcvber_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_flowcover_cuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_rounding_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_rounding_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_rounding_cuts;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_lift_and_project_cuts") == 0 ||
	    strcmp(key, "LP_generate_cgl_lift_and_project_cuts") == 0){
      *value = lp_par->cgl.generate_cgl_lift_and_project_cuts;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_landp_cuts") == 0 ||
            strcmp(key, "LP_generate_cgl_landp_cuts") == 0){
     *value = lp_par->cgl.generate_cgl_landp_cuts;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_gomory_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_gomory_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_gomory_cuts_freq;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_redsplit_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_redsplit_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_gomory_cuts_freq;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_knapsack_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_knapsack_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_knapsack_cuts_freq;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_oddhole_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_oddhole_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_oddhole_cuts_freq;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_probing_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_probing_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_probing_cuts_freq;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_clique_cuts_freq") == 0 ||
            strcmp(key, "LP_generate_cgl_clique_cuts_freq") == 0){
     *value = lp_par->cgl.generate_cgl_clique_cuts_freq;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_mir_cuts_freq") == 0 ||
            strcmp(key, "LP_generate_cgl_mir_cuts_freq") == 0){
     *value = lp_par->cgl.generate_cgl_mir_cuts_freq;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_twomir_cuts_freq") == 0 ||
            strcmp(key, "LP_generate_cgl_twomir_cuts_freq") == 0){
     *value = lp_par->cgl.generate_cgl_twomir_cuts_freq;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_flowcover_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_flowcvber_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_flowcover_cuts_freq;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_rounding_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_rounding_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_rounding_cuts_freq;
      return(0);
   }
   else if (strcmp(key, "generate_cgl_lift_and_project_cuts_freq") == 0 ||
	    strcmp(key, "LP_generate_cgl_lift_and_project_cuts_freq") == 0){
      *value = lp_par->cgl.generate_cgl_lift_and_project_cuts_freq;
     return(0);
   }
   else if (strcmp(key, "generate_cgl_landp_cuts_freq") == 0 ||
            strcmp(key, "LP_generate_cgl_landp_cuts_freq") == 0){
     *value = lp_par->cgl.generate_cgl_landp_cuts_freq;
     return(0);
   }
   else if (strcmp(key, "max_presolve_iter") == 0 ||
	    strcmp(key, "LP_max_presolve_iter") == 0){
      *value = lp_par->max_presolve_iter;
      return(0);
   }

   /* user-defined function defaults */
   else if (strcmp(key, "is_feasible_default") == 0 ||
	    strcmp(key, "LP_is_feasible_default") == 0){
      *value = lp_par->is_feasible_default;
      return(0);
   }
   else if (strcmp(key, "send_feasible_solution_default") == 0 ||
	    strcmp(key, "LP_send_feasible_solution_default") == 0){
      *value = lp_par->send_feasible_solution_default;
      return(0);
   }
   else if (strcmp(key, "display_solution_default") == 0 ||
	    strcmp(key, "LP_display_solution_default") == 0){
      *value = lp_par->display_solution_default;
      return(0);
   }
   else if (strcmp(key, "shall_we_branch_default") == 0 ||
	    strcmp(key, "LP_shall_we_branch_default") == 0){
      *value = lp_par->shall_we_branch_default;
      return(0);
   }
   else if (strcmp(key, "select_candidates_default") == 0 ||
	    strcmp(key, "LP_select_candidates_default") == 0){
      *value = lp_par->select_candidates_default;
      return(0);
   }
   else if (strcmp(key, "strong_branching_cand_num") == 0){
      *value = lp_par->strong_branching_cand_num_max;
      return(0);
   }
   else if (strcmp(key, "strong_branching_cand_num_max") == 0 ||
	    strcmp(key, "LP_strong_branching_cand_num_max") == 0){
      *value = lp_par->strong_branching_cand_num_max;
      return(0);
   }
   else if (strcmp(key, "strong_branching_cand_num_min") == 0 ||
	    strcmp(key, "LP_strong_branching_cand_num_min") == 0){
      *value = lp_par->strong_branching_cand_num_min;
      return(0);
   }
   else if (strcmp(key, "user_set_strong_branching_cand_num") == 0) {
      *value = lp_par->user_set_strong_branching_cand_num;
      return(0);
   }
   else if (strcmp(key, "user_set_max_presolve_iter") == 0) {
      *value = lp_par->user_set_max_presolve_iter;
      return(0);
   }
   else if (strcmp(key, "strong_br_min_level") == 0) {
      *value = lp_par->strong_br_min_level;
      return(0);
   }
   else if (strcmp(key, "strong_br_all_candidates_level") == 0) {
      *value = lp_par->strong_br_all_candidates_level;
      return(0);
   }
   else if (strcmp(key,"use_hot_starts") == 0) {
      *value = lp_par->use_hot_starts;
      return(0);
   }
   else if (strcmp(key,"rel_br_threshold") == 0) {
      *value = lp_par->rel_br_threshold;
      return(0);
   }
   else if (strcmp(key,"rel_br_cand_threshold") == 0) {
      *value = lp_par->rel_br_cand_threshold;
      return(0);
   }
   else if (strcmp(key,"should_use_rel_br") == 0) {
      *value = lp_par->should_use_rel_br;
      return(0);
   }
   else if (strcmp(key, "compare_candidates_default") == 0 ||
	    strcmp(key, "LP_compare_candidates_default") == 0){
      *value = lp_par->compare_candidates_default;
      return(0);
   }

   else if (strcmp(key, "select_child_default") == 0 ||
	    strcmp(key, "LP_select_child_default") == 0){
      *value = lp_par->select_child_default;
      return(0);
   }
   else if (strcmp(key, "pack_lp_solution_default") == 0 ||
	       strcmp(key, "LP_pack_lp_solution_default") == 0){
      *value = lp_par->pack_lp_solution_default;
      return(0);
   }
   else if (strcmp(key, "multi_criteria") == 0 ||
	    strcmp(key, "LP_multi_criteria") == 0 ){
      *value = lp_par->multi_criteria;
      return(0);
   }
   else if (strcmp(key, "mc_find_supported_solutions") == 0 ||
	    strcmp(key, "LP_mc_find_supported_solutions") == 0 ){
      *value = lp_par->mc_find_supported_solutions;
      return(0);
   }
   else if (strcmp(key, "mc_add_optimality_cuts") == 0 ||
	    strcmp(key, "LP_mc_add_optimality_cuts") == 0 ){
      *value = lp_par->mc_add_optimality_cuts;
      return(0);
   }
   else if (strcmp(key, "fp_enabled") == 0) {
      *value = lp_par->fp_enabled;
      return(0);
   }
   else if (strcmp(key, "fp_frequency") == 0) {
      *value = lp_par->fp_frequency;
      return(0);
   }
   else if (strcmp(key, "fp_max_cycles") == 0) {
      *value = lp_par->fp_max_cycles;
      return(0);
   }

   /***********************************************************************
    ***                     cut_gen params                          ***
    ***********************************************************************/
   if (strcmp(key, "CG_verbosity") == 0){
      *value = cg_par->verbosity;
      return(0);
   }
   else if (strcmp(key, "do_findcuts") == 0 ||
	    strcmp(key, "CG_do_findcuts") == 0){
      *value = cg_par->do_findcuts;
      return(0);
   }

   /***********************************************************************
    ***                      cutpool params                         ***
    ***********************************************************************/
   else if (strcmp(key, "CP_verbosity") == 0){
      *value = cp_par->verbosity;
      return(0);
   }
   else if (strcmp(key, "cp_warm_start") == 0 ||
	    strcmp(key, "CP_warm_start") == 0){
      *value = cp_par->warm_start;
      return(0);
   }
   else if (strcmp(key, "cp_logging") == 0 ||
	    strcmp(key, "CP_logging") == 0){
      *value = cp_par->logging;
      return(0);
   }
   else if (strcmp(key, "block_size") == 0 ||
	    strcmp(key, "CP_block_size") == 0){
      *value = cp_par->block_size;
      return(0);
   }
   else if (strcmp(key, "max_size") == 0 ||
	    strcmp(key, "CP_max_size") == 0){
      *value = cp_par->max_size;
      return(0);
   }
   else if (strcmp(key, "max_number_of_cuts") == 0 ||
	    strcmp(key, "CP_max_number_of_cuts") == 0){
      *value = cp_par->max_number_of_cuts;
      return(0);
   }
   else if (strcmp(key, "cuts_to_check") == 0 ||
	    strcmp(key, "cuts_to_check") == 0){
      *value = cp_par->cuts_to_check;
      return(0);
   }
   else if (strcmp(key, "delete_which") == 0 ||
	    strcmp(key, "CP_delete_which") == 0){
      *value = cp_par->delete_which;
      return(0);
   }
   else if (strcmp(key, "touches_until_deletion") == 0 ||
	    strcmp(key, "CP_touches_until_deletion") == 0){
      *value = cp_par->touches_until_deletion;
      return(0);
   }
   else if (strcmp(key, "min_to_delete") == 0 ||
	    strcmp(key, "CP_min_to_delete") == 0){
      *value = cp_par->min_to_delete;
      return(0);
   }
   else if (strcmp(key, "check_which") == 0 ||
         strcmp(key, "CP_check_which") == 0){
      *value = cp_par->check_which;
   }

   return (FUNCTION_TERMINATED_ABNORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_dbl_param(sym_environment *env, const char *key, double *value)
{

   tm_params *tm_par = &env->par.tm_par;
   lp_params *lp_par = &env->par.lp_par;
   //cg_params *cg_par = &env->par.cg_par;
   //cp_params *cp_par = &env->par.cp_par;

   dg_params *dg_par = &env->par.dg_par;

   if (strcmp(key, "granularity") == 0){
      *value = tm_par->granularity;
      return(0);
   }
   else if (strcmp(key, "upper_bound") == 0 ||
	    strcmp(key, "M_upper_bound") == 0){
      *value = env->ub;
      return(0);
   }
   else if (strcmp(key, "upper_bound_estimate") == 0 ||
	    strcmp(key, "M_upper_bound_estimate") == 0){
      *value = env->ub_estimate;
      return(0);
   }
   else if (strcmp(key, "lower_bound") == 0 ||
	    strcmp(key, "M_lower_bound") == 0){
      *value = env->lb;
      return(0);
   }
   else if (strcmp(key, "obj_offset") == 0 ||
	    strcmp(key, "M_obj_offset") == 0){
      *value = env->obj_offset;
      return(0);
   }
   else if (strcmp(key, "scale_factor") == 0 ||
	    strcmp(key, "DG_scale_factor") == 0){
      *value = dg_par->scale_factor;
      return(0);
   }
   else if (strcmp(key, "mc_compare_solution_tolerance") == 0 ||
	    strcmp(key, "M_mc_compare_solution_tolerance") == 0){
      *value = env->par.mc_compare_solution_tolerance;
      return(0);
   }
   else if (strcmp(key, "mc_binary_search_tolerance") == 0 ||
	    strcmp(key, "M_mc_binary_search_tolerance") == 0){
      *value = env->par.mc_binary_search_tolerance;
      return(0);
   }

   /***********************************************************************
    ***                  Treemanager params                         ***
    ***********************************************************************/
   else if (strcmp(key, "TM_granularity") == 0){
      *value = tm_par->granularity;
      return(0);
   }
   else if (strcmp(key, "unconditional_dive_frac") == 0 ||
	    strcmp(key, "TM_unconditional_dive_frac") == 0){
      *value = tm_par->unconditional_dive_frac;
      return(0);
   }
   else if (strcmp(key, "diving_threshold") == 0 ||
	    strcmp(key, "TM_diving_threshold") == 0){
     *value =  tm_par->diving_threshold;
      return(0);
   }
   else if (strcmp(key, "time_limit") == 0 ||
	    strcmp(key, "TM_time_limit") == 0){
     *value =  tm_par->time_limit;
      return(0);
   }
   else if (strcmp(key, "gap_limit") == 0 ||
	    strcmp(key, "TM_gap_limit") == 0){
      *value = tm_par->gap_limit;
      return(0);
   }

   /***********************************************************************
    ***                      LP params                              ***
    ***********************************************************************/
   else if (strcmp(key, "LP_granularity") == 0){
      *value = lp_par->granularity;
      return(0);
   }
   else if (strcmp(key, "fractional_diving_ratio") == 0 ||
	    strcmp(key, "LP_fractional_diving_ratio") == 0){
      *value = lp_par->fractional_diving_ratio;
      return(0);
   }
   else if (strcmp(key, "max_non_dual_feas_to_add_frac") == 0 ||
	    strcmp(key, "LP_max_non_dual_feas_to_add_frac") == 0){
      *value = lp_par->max_non_dual_feas_to_add_frac;
      return(0);
   }
   else if (strcmp(key, "max_not_fixable_to_add_frac") == 0 ||
	    strcmp(key, "LP_max_not_fixable_to_add_frac") == 0){
      *value = lp_par->max_not_fixable_to_add_frac;
      return(0);
   }
   else if (strcmp(key, "mat_col_compress_ratio") == 0 ||
	    strcmp(key, "LP_mat_col_compress_ratio") == 0){
      *value = lp_par->mat_col_compress_ratio;
      return(0);
   }
   else if (strcmp(key, "mat_row_compress_ratio") == 0 ||
	    strcmp(key, "LP_mat_row_compress_ratio") == 0){
      *value = lp_par->mat_row_compress_ratio;
      return(0);
   }
   else if (strcmp(key, "tailoff_gap_frac") == 0 ||
	    strcmp(key, "LP_tailoff_gap_frac") == 0){
      *value = lp_par->tailoff_gap_frac;
      return(0);
   }
   else if (strcmp(key, "tailoff_obj_frac") == 0 ||
	    strcmp(key, "LP_tailoff_obj_frac") == 0){
      *value = lp_par->tailoff_obj_frac;
      return(0);
   }
   else if (strcmp(key, "tailoff_absolute") == 0 ||
	    strcmp(key, "LP_tailoff_absolute") == 0){
      *value = lp_par->tailoff_absolute;
      return(0);
   }
   else if (strcmp(key, "tailoff_max_no_iterative_impr_iters_root") == 0 ||
	    strcmp(key, "LP_tailoff_max_no_iterative_impr_iters_root") == 0){
      *value = lp_par->tailoff_max_no_iterative_impr_iters_root;
      return(0);
   }

   /* timeouts on receiving cuts */
   else if (strcmp(key, "first_lp_first_cut_time_out") == 0 ||
	    strcmp(key, "LP_first_lp_first_cut_time_out") == 0){
      *value = lp_par->first_lp.first_cut_time_out;
      return(0);
   }
   else if (strcmp(key, "first_lp_all_cuts_time_out") == 0 ||
	    strcmp(key, "LP_first_lp_all_cuts_time_out") == 0){
      *value = lp_par->first_lp.all_cuts_time_out;
      return(0);
   }
   else if (strcmp(key, "later_lp_first_cut_time_out") == 0 ||
	    strcmp(key, "LP_later_lp_first_cut_time_out") == 0){
      *value = lp_par->later_lp.first_cut_time_out;
      return(0);
   }
   else if (strcmp(key, "later_lp_all_cuts_time_out") == 0 ||
	    strcmp(key, "LP_later_lp_all_cuts_time_out") == 0){
      *value = lp_par->later_lp.all_cuts_time_out;
      return(0);
   }

   else if (strcmp(key, "gap_as_ub_frac") == 0 ||
	    strcmp(key, "LP_gap_as_ub_frac") == 0){
      *value = lp_par->gap_as_ub_frac;
      return(0);
   }
   else if (strcmp(key, "gap_as_last_gap_frac") == 0 ||
	    strcmp(key, "LP_gap_as_last_gap_frac") == 0){
      *value = lp_par->gap_as_last_gap_frac;
      return(0);
   }
   else if (strcmp(key, "fixed_to_ub_frac_before_logical_fixing")==0 ||
	    strcmp(key, "LP_fixed_to_ub_frac_before_logical_fixing")==0){
      *value = lp_par->fixed_to_ub_frac_before_logical_fixing;
      return(0);
   }
   else if (strcmp(key,"strong_branching_red_ratio") == 0 ||
	    strcmp(key,"LP_strong_branching_red_ratio") == 0){
      *value = lp_par->strong_branching_red_ratio;
      return(0);
   }
   else if (strcmp(key,"strong_branching_high_low_weight") == 0 ||
	    strcmp(key,"LP_strong_branching_high_low_weight") == 0){
      *value = lp_par->strong_branching_high_low_weight;
      return(0);
   }
  else if (strcmp(key, "mc_gamma") == 0 ||
	    strcmp(key, "LP_mc_gamma") == 0 ){
      *value = lp_par->mc_gamma;
      return(0);
   }
   else if (strcmp(key, "mc_tau") == 0 ||
	    strcmp(key, "LP_mc_tau") == 0 ){
      *value = lp_par->mc_tau;
      return(0);
   }
   else if (strcmp(key, "mc_rho") == 0 ||
	    strcmp(key, "LP_mc_rho") == 0 ){
      *value = lp_par->mc_rho;
      return(0);
   }
   else if (strcmp(key, "fp_time_limit") == 0) {
      *value = lp_par->fp_time_limit;
      return(0);
   }
   else if (strcmp(key, "fp_flip_fraction") == 0) {
      *value = lp_par->fp_flip_fraction;
      return(0);
   }
   else if (strcmp(key, "fp_max_initial_time") == 0) {
      *value = lp_par->fp_max_initial_time;
      return(0);
   }
   else if (strcmp(key, "fp_display_time") == 0) {
      *value = lp_par->fp_display_interval;
      return(0);
   }
   else if (strcmp(key, "fp_min_gap") == 0) {
      *value = lp_par->fp_min_gap;
      return(0);
   }


   return (FUNCTION_TERMINATED_ABNORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_str_param(sym_environment *env, const char *key, char **value)
{

   tm_params *tm_par = &env->par.tm_par;
   //lp_params *lp_par = &env->par.lp_par;
   //cg_params *cg_par = &env->par.cg_par;
   //cp_params *cp_par = &env->par.cp_par;

   dg_params *dg_par = &env->par.dg_par;

   if (strcmp(key, "problem_name") == 0){
      *value = env->probname;
      return(0);
   }
   else if(strcmp(key, "infile_name") == 0){
      *value = env->par.infile;
      return(0);
   }
   else if (strcmp(key, "tm_executable_name") == 0 ||
	    strcmp(key, "tm_exe") == 0 ||
	    strcmp(key, "M_tm_exe") == 0 ||
	    strcmp(key, "M_tm_executable_name") == 0){
      *value = env->par.tm_exe;
      return(0);
   }
   else if (strcmp(key, "dg_executable_name") == 0 ||
	    strcmp(key, "dg_exe") == 0 ||
	    strcmp(key, "M_dg_exe") == 0 ||
	    strcmp(key, "M_dg_executable_name") == 0){
      *value = env->par.dg_exe;
      return(0);
   }
   else if (strcmp(key, "tm_machine") == 0 ||
	    strcmp(key, "M_tm_machine") == 0){
      *value = env->par.tm_machine;
      return(0);
   }
   else if (strcmp(key, "dg_machine") == 0 ||
	    strcmp(key, "M_dg_machine") == 0){
      *value = env->par.dg_machine;
      return(0);
   }
   else if (strcmp(key, "param_file") == 0 ||
	    strcmp(key, "M_param_file") == 0){
      *value = env->par.param_file;
      return(0);
   }

   /***********************************************************************
    ***                 DrawGraph params                            ***
    ***********************************************************************/

   else if (strcmp(key, "source_path") == 0 ||
	    strcmp(key, "DG_source_path") == 0){
      *value = dg_par->source_path;
      return(0);
   }
   else if (strcmp(key, "node_dash") == 0 ||
	    strcmp(key, "DG_node_dash") == 0){
      *value = dg_par->node_dash;
      return(0);
   }
   else if (strcmp(key, "edge_dash") == 0 ||
	    strcmp(key, "DG_edge_dash") == 0){
      *value = dg_par->edge_dash;
      return(0);
   }
   else if (strcmp(key, "nodelabel_font") == 0 ||
	    strcmp(key, "DG_nodelabel_font") == 0){
      *value = dg_par->nodelabel_font;
      return(0);
   }
   else if (strcmp(key, "nodeweight_font") == 0 ||
	    strcmp(key, "DG_nodeweight_font") == 0){
      *value = dg_par->nodeweight_font;
      return(0);
   }
   else if (strcmp(key, "edgeweight_font") == 0 ||
	    strcmp(key, "DG_edgeweight_font") == 0){
      *value = dg_par->edgeweight_font;
      return(0);
   }

   /***********************************************************************
    ***                  Treemanager params                         ***
    ***********************************************************************/
   else if (strcmp(key, "lp_executable_name") == 0 ||
	    strcmp(key, "lp_exe") == 0 ||
	    strcmp(key, "TM_lp_exe") == 0 ||
	    strcmp(key, "TM_lp_executable_name") == 0){
      *value = tm_par->lp_exe;
      return(0);
   }
   else if (strcmp(key, "cg_executable_name") == 0 ||
	    strcmp(key, "cg_exe") == 0 ||
	    strcmp(key, "TM_cg_exe") == 0 ||
	    strcmp(key, "TM_cg_executable_name") == 0){
      *value = tm_par->cg_exe;
      return(0);
   }
   else if (strcmp(key, "cp_executable_name") == 0 ||
	    strcmp(key, "cp_exe") == 0 ||
	    strcmp(key, "TM_cp_exe") == 0 ||
	    strcmp(key, "TM_cp_executable_name") == 0){
      *value = tm_par->cp_exe;
      return(0);
   }
   return (FUNCTION_TERMINATED_ABNORMALLY);
}

/*===========================================================================*/
/*===========================================================================*/

warm_start_desc *sym_create_copy_warm_start(warm_start_desc *ws)
{
   return(create_copy_warm_start(ws));
}

/*===========================================================================*/
/*===========================================================================*/

MIPdesc *sym_create_copy_mip_desc(sym_environment *env)
{
   return(create_copy_mip_desc(env->mip));
}

/*===========================================================================*/
/*===========================================================================*/

MIPdesc *sym_get_presolved_mip_desc(sym_environment *env)
{
   return(env->prep_mip);
}

/*===========================================================================*/
/*===========================================================================*/

sym_environment * sym_create_copy_environment (sym_environment *env)
{
   return(create_copy_environment(env));
}

/*===========================================================================*/
/*===========================================================================*/

int sym_get_lb_for_new_rhs(sym_environment *env, int cnt, int *new_rhs_ind,
			      double *new_rhs_val, double *lb_for_new_rhs)
{
#ifdef SENSITIVITY_ANALYSIS
#ifdef USE_CGL_CUTS
   printf("sym_get_lb_for_new_rhs():\n");
   printf("SYMPHONY can not do sensitivity analysis when cuts are present, for now!\n");
   return(FUNCTION_TERMINATED_ABNORMALLY);
#else
   if (!env || !env->mip ||
      !env->par.tm_par.sensitivity_analysis){
      printf("sym_get_lb_for_new_rhs():\n");
      printf("Trying to read an empty problem, an empty problem description");
      printf(" or tree nodes were not kept in memory!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   else{
      if (!env->warm_start){
	 printf("sym_get_lb_for_new_rhs():\n");
	 printf("No available warm start data to do sens. analysis. \n");
	 return(FUNCTION_TERMINATED_ABNORMALLY);
      }
      else{
	 /* check if we only have the root node, then no need to call
	    recursive algorithm */
	 int i;
	 if(env->warm_start->stat.analyzed == 1) {
	    *lb_for_new_rhs =  env->warm_start->rootnode->lower_bound;
	    for(i=0; i<cnt; i++){
	       *lb_for_new_rhs +=
		  env->warm_start->rootnode->duals[new_rhs_ind[i]]*
		  (new_rhs_val[i] - env->mip->rhs[new_rhs_ind[i]]);
	    }
	 } else {
	    *lb_for_new_rhs =
	       get_lb_for_new_rhs(env->warm_start->rootnode, env->mip, cnt,
				  new_rhs_ind, new_rhs_val);
	 }
      }
      return(FUNCTION_TERMINATED_NORMALLY);
   }
#endif
#else
   printf("sym_get_lb_for_new_rhs():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(FUNCTION_TERMINATED_ABNORMALLY);
#endif
 }

/*===========================================================================*/
/*===========================================================================*/

int sym_get_ub_for_new_rhs(sym_environment *env, int cnt, int *new_rhs_ind,
			   double *new_rhs_val, double *ub_for_new_rhs)
{
#ifdef SENSITIVITY_ANALYSIS
   int *matbeg = NULL, *matind = NULL, nz, i, j, k;
   double *matval = NULL;

   if (!env || !env->mip ||
       !env->par.tm_par.sensitivity_analysis){
      printf("sym_get_ub_for_new_rhs():\n");
      printf("Trying to read an empty problem, an empty problem description");
      printf(" or tree nodes were not kept in memory!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   else{
      if (!env->warm_start){
	 printf("sym_get_ub_for_new_rhs():\n");
	 printf("No available warm start data to do sens. analysis. \n");
	 return(FUNCTION_TERMINATED_ABNORMALLY);
      }
      else{

	 /* prepare to send a row oriented mip description in to
	    get_ub_for_new_rhs() */

	 matbeg = env->mip->matbeg;
	 matind = env->mip->matind;
	 matval = env->mip->matval;

	 env->mip->matbeg = (int *) calloc(ISIZE, (env->mip->m+1));
	 env->mip->matind = (int *) malloc (ISIZE* env->mip->nz);
	 env->mip->matval = (double *) malloc (DSIZE* env->mip->nz);

	 nz = 0;
	 for(j = 0; j < env->mip->n; j++){
	    for(i = 0; i < env->mip->m; i++){
	       for(k = matbeg[i]; k < matbeg[i+1]; k++){
		  if(matind[k] == j){
		     env->mip->matind[nz] = i;
		     env->mip->matval[nz] = matval[k];
		     nz++;
		     break;
		  }
	       }
	    }
	    env->mip->matbeg[j+1] = nz;
	 }

	 *ub_for_new_rhs =
	    get_ub_for_new_rhs(env->warm_start->rootnode, env->mip, cnt,
			       new_rhs_ind, new_rhs_val);

	 FREE(env->mip->matbeg);
	 FREE(env->mip->matind);
	 FREE(env->mip->matval);
	 env->mip->matbeg = matbeg;
	 env->mip->matind = matind;
	 env->mip->matval = matval;
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
#else
   printf("sym_get_ub_for_new_rhs():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(FUNCTION_TERMINATED_ABNORMALLY);
#endif
}

/*===========================================================================*/
/*===========================================================================*/
#if 0
int sym_get_lb_for_new_obj(sym_environment *env, int cnt,
				 int *new_obj_ind,
				 double *new_obj_val,
				 double *lb_for_new_obj)
{
#ifdef SENSITIVITY_ANALYSIS
   double ub;

   if (!env || !env->mip ||
       !env->par.tm_par.sensitivity_analysis){
      printf("sym_get_lb_for_new_obj():\n");
      printf("Trying to read an empty problem, an empty problem description");
      printf(" or tree nodes were not kept in memory!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   else{
      if (!env->warm_start){
	 printf("sym_get_lb_for_new_obj():\n");
	 printf("No available warm start data to do sens. analysis. \n");
	 return(FUNCTION_TERMINATED_ABNORMALLY);
      }
      else{
	 ub = get_ub_for_new_obj(env->warm_start->rootnode, env->mip, cnt,
				 new_obj_ind, new_obj_val);
	 *lb_for_new_obj =
	    get_lb_for_new_obj(env->warm_start->rootnode, env->mip, cnt,
			       new_obj_ind, new_obj_val);
	 if(*lb_for_new_obj > ub){
	    *lb_for_new_obj = ub;
	 }
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
#else
   printf("sym_get_lb_for_new_obj():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(FUNCTION_TERMINATED_ABNORMALLY);
#endif
}
#endif
/*===========================================================================*/
/*===========================================================================*/

int sym_get_ub_for_new_obj(sym_environment *env, int cnt,
				 int *new_obj_ind,
				 double *new_obj_val,
				 double *ub_for_new_obj)
{
#ifdef SENSITIVITY_ANALYSIS
   if (!env || !env->mip ||
       !env->par.tm_par.sensitivity_analysis){
      printf("sym_get_ub_for_new_obj():\n");
      printf("Trying to read an empty problem, an empty problem description");
      printf(" or tree nodes were not kept in memory!\n");
      return(FUNCTION_TERMINATED_ABNORMALLY);
   }
   else{
      if (!env->warm_start){
	 printf("sym_get_ub_for_new_obj():\n");
	 printf("No available warm start data to do sens. analysis. \n");
	 return(FUNCTION_TERMINATED_ABNORMALLY);
      }
      else{
	 *ub_for_new_obj =
	    get_ub_for_new_obj(env->warm_start->rootnode, env->mip, cnt,
			       new_obj_ind, new_obj_val);
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
#else
   printf("sym_get_ub_for_new_obj():\n");
   printf("Sensitivity analysis features are not enabled.\n");
   printf("Please rebuild SYMPHONY with these features enabled\n");
   return(FUNCTION_TERMINATED_ABNORMALLY);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int sym_test(sym_environment *env, int *test_status)
{

  int termcode = 0, verbosity;
  int i, file_num = 12;
  char mps_files[12][MAX_FILE_NAME_LENGTH +1] = {
    "air03", "dcmulti", "egout", "flugpl", "khb05250", "l152lav",
    "lseu", "mod010", "p0033", "p0201", "stein27", "vpm1" };

  double sol[12] = {340160, 188182, 568.101, 1201500,
			  106940226, 4722, 1120, 6548,
			  3089, 7615, 18, 20};

  char *mps_dir = (char*)malloc(CSIZE*(MAX_FILE_NAME_LENGTH+1));
  char *infile = (char*)malloc(CSIZE*(MAX_FILE_NAME_LENGTH+1));
  double *obj_val = (double *)calloc(DSIZE,file_num);
  double tol = 1e-03;

  //size_t size = 1000;
  int size = 1000;
  char* buf = 0;

  *test_status = 0;
  verbosity = sym_get_int_param(env, "verbosity", &verbosity);

  while (true) {
     buf = (char*)malloc(CSIZE*size);
     if (getcwd(buf, size))
	break;
     FREE(buf);
     buf = 0;
     size = 2*size;
  }
  char dirsep = buf[0] == '/' ? '/' : '\\';
  FREE(buf);

  if (strcmp(env->par.test_dir, "") == 0){
     if (dirsep == '/')
	strcpy(mps_dir, "../../Data/miplib3");
     else
	strcpy(mps_dir, "..\\..\\Data\\miplib3");
  } else{
    strcpy(mps_dir, env->par.test_dir);
  }

  for(i = 0; i<file_num; i++){

    if(env->mip->n){
      free_master_u(env);
      strcpy(env->par.infile, "");
      env->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));
    }
    sym_set_defaults (env);
    sym_set_int_param(env, "verbosity", -10);

    strcpy(infile, "");
    if (dirsep == '/')
       sprintf(infile, "%s%s%s", mps_dir, "/", mps_files[i]);
    else
       sprintf(infile, "%s%s%s", mps_dir, "\\", mps_files[i]);
    if((termcode = sym_read_mps(env, infile)) < 0)
      return(termcode);

    printf("\nSolving %s...\n\n", mps_files[i]);

    if((termcode = sym_solve(env)) < 0)
      return(termcode);

    sym_get_obj_val(env, &obj_val[i]);

    if((obj_val[i] < sol[i] + tol) &&
       (obj_val[i] > sol[i] - tol)){
      printf("\nSuccess! %s solved correctly...\n\n", mps_files[i]);
    } else {
       printf("\nFailure! Solver returned solution value: %f", obj_val[i]);
       printf("\n         True solution value is:         %f\n\n", sol[i]);
       *test_status = 1;
    }
  }

  FREE(mps_dir);
  FREE(infile);
  FREE(obj_val);

  sym_set_int_param(env, "verbosity", verbosity);

  return(termcode);

}

