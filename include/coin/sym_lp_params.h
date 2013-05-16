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

#ifndef _LP_PARAMS_H
#define _LP_PARAMS_H

#include "sym_constants.h"
#include "sym_timemeas.h"

/*---------------------------------------------------------------------------*\
| The list of parameters associated with processing a node in the branch and  |
| cut tree. See the README file for an explanation of the parameters          |
\*---------------------------------------------------------------------------*/

typedef struct CUT_TIME_OUT{
   double            first_cut_time_out;
   double            all_cuts_time_out;
}cut_time_out;

typedef struct CGL_PARAMS{
   /* Cut generation in LP */
   int               generate_cgl_cuts;
   int               max_depth_for_cgl_cuts;
   int               generate_cgl_gomory_cuts;
   int               generate_cgl_redsplit_cuts;
   int               generate_cgl_knapsack_cuts;
   int               generate_cgl_oddhole_cuts;
   int               generate_cgl_probing_cuts;
   int               generate_cgl_mir_cuts;
   int               generate_cgl_twomir_cuts;
   int               generate_cgl_clique_cuts;
   int               generate_cgl_flowcover_cuts;
   int               generate_cgl_rounding_cuts;
   int               generate_cgl_lift_and_project_cuts;
   int               generate_cgl_landp_cuts;

   int               generate_cgl_gomory_cuts_freq;
   int               generate_cgl_redsplit_cuts_freq;
   int               generate_cgl_knapsack_cuts_freq;
   int               generate_cgl_oddhole_cuts_freq;
   int               generate_cgl_probing_cuts_freq;
   int               generate_cgl_mir_cuts_freq;
   int               generate_cgl_twomir_cuts_freq;
   int               generate_cgl_clique_cuts_freq;
   int               generate_cgl_flowcover_cuts_freq;
   int               generate_cgl_rounding_cuts_freq;
   int               generate_cgl_lift_and_project_cuts_freq;
   int               generate_cgl_landp_cuts_freq;

   int               gomory_generated_in_root;
   int               redsplit_generated_in_root;
   int               knapsack_generated_in_root;
   int               oddhole_generated_in_root;
   int               probing_generated_in_root;
   int               mir_generated_in_root;
   int               twomir_generated_in_root;
   int               clique_generated_in_root;
   int               flowcover_generated_in_root;
   int               rounding_generated_in_root;
   int               lift_and_project_generated_in_root;
   int               landp_generated_in_root;

   int               probing_is_expensive;
   int               probing_root_max_look;

   int               gomory_max_depth;
   int               probing_max_depth;
   int               flowcover_max_depth;
   int               twomir_max_depth;
   int               clique_max_depth;
   int               oddhole_max_depth;
   int               knapsack_max_depth;

   int               use_chain_strategy;
   int               chain_status;
   int               max_chain_backtrack;
   int               max_chain_trial_num;
   int               chain_trial_freq;
   int               chain_check_index;
   double            chain_weighted_gap;
}cgl_params;

typedef struct LP_PARAMS{
   int               verbosity;
   double            granularity;
   int               use_cg;
   int               set_obj_upper_lim;
   int               do_primal_heuristic;
   double            time_limit;

   int               lp_data_mip_is_copied;
   /* TRUE: save the base model after root solve and then load it each time we
    * start a new chain (dive). FALSE: load the model from scratch. Basis
    * information is loaded separately in both cases for a warm start. Cant be
    * set by user.
    */
   int               should_reuse_lp;

   /* these two are passed directly to the lp solver */
   int               scaling;
   int               fastmip;

   /*
    * should we do initial_solve() or dual_simplex() when we start a new
    * chain. both have pros and cons and asm4 is not sure what to do.
    */
   int               should_warmstart_chain;


   int               try_to_recover_from_error;
   /* ZERO_ONE_PROBLEM / INTEGER_PROBLEM / MIXED_INTEGER_PROBLEM */
   int               problem_type;
   int               keep_description_of_pruned;

   int               not_fixed_storage_size;

   int               cut_pool_check_freq;

   int               load_balance_level;
   int               load_balance_iterations;
   int               load_balance_compare_candidates;

   double            fractional_diving_ratio;
   int               fractional_diving_num;

   /* parameters constraining the growth of the matrix */
   double            max_non_dual_feas_to_add_frac;
   int               max_non_dual_feas_to_add_min;
   int               max_non_dual_feas_to_add_max;
   double            max_not_fixable_to_add_frac;
   int               max_not_fixable_to_add_min;
   int               max_not_fixable_to_add_max;

   int               mat_col_compress_num;
   double            mat_col_compress_ratio;
   int               mat_row_compress_num;
   double            mat_row_compress_ratio;

   /* parameters governing tailing off checking */
   int               tailoff_gap_backsteps;
   double            tailoff_gap_frac;
   int               tailoff_obj_backsteps;
   double            tailoff_obj_frac;
   double            tailoff_absolute;
   int               tailoff_max_no_iterative_impr_iters_root;

   int               ineff_cnt_to_delete;
   int               eff_cnt_before_cutpool;
   int               ineffective_constraints;
   int               base_constraints_always_effective;

   int               branch_on_cuts;   /* TRUE / FALSE */
   int               discard_slack_cuts;

   cut_time_out      first_lp;
   cut_time_out      later_lp;

   int               max_cut_num_per_iter;
   int               max_cut_num_per_iter_root;
   int               min_root_cut_rounds;
   int               max_cut_length;
   int               tried_long_cuts;

   /* Reduced cost and logical fixing parameters */
   int               do_reduced_cost_fixing;
   double            gap_as_ub_frac;
   double            gap_as_last_gap_frac;
   int               do_logical_fixing;
   int               fixed_to_ub_before_logical_fixing; /* OK */
   double            fixed_to_ub_frac_before_logical_fixing; /* OK */

   /* CGL parameters */
   cgl_params        cgl;

   /* Parameters affecting branching */
   int               max_presolve_iter;

   /*Defaults for the user supplied routines*/
   int               is_feasible_default;
   int               send_feasible_solution_default;
   int               display_solution_default;
   int               shall_we_branch_default;
   int               select_candidates_default;
   int               strong_branching_cand_num_min;
   int               strong_branching_cand_num_max;
   double            strong_branching_red_ratio;
   double            strong_branching_high_low_weight;
   int               use_hot_starts;
   int               strong_br_all_candidates_level;
   int               strong_br_min_level;
   int               user_set_strong_branching_cand_num;
   int               user_set_max_presolve_iter;
   int               should_use_rel_br;
   int               rel_br_override_default;
   int               rel_br_override_max_solves;
   int               rel_br_chain_backtrack;
   double            rel_br_min_imp;
   double            rel_br_max_imp;

   int               rel_br_threshold; /* how many times to do strong branching
                                          on each variable before using pseudo
                                          cost estimates */
   int               rel_br_cand_threshold; /* how many candidates to solve
                                               using strong branching without
                                               any improvement in score before
                                               stopping */
   int               rel_br_max_solves; /* stop after these many LP-solve calls
                                           regardless of improvement */

   int               compare_candidates_default;
   int               select_child_default;
   int               pack_lp_solution_default;

   /* Multi-criteria parameters */
   int               multi_criteria;
   int               mc_find_supported_solutions;
   int               mc_add_optimality_cuts;
   double            mc_rho;   /* For augmented Chebyshev norm */
   double            mc_gamma; /* Weight on first objective */
   double            mc_tau;   /* Weight on second objective */

   int               sensitivity_analysis;

   /* feasibility pump parameters */
   int               fp_enabled;
   int               fp_frequency;
   int               fp_max_cycles;
   int               fp_poor_sol_lim_fac;
   double            fp_time_limit;
   double            fp_display_interval;
   double            fp_flip_fraction;
   double            fp_max_initial_time;
   double            fp_min_gap;

   /* to avoid nested for loops, check if userind's are in order */
   int               is_userind_in_order;
}lp_params;

#endif
