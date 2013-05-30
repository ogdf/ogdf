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

#ifndef _MASTER_PARAMS_H
#define _MASTER_PARAMS_H

#include "sym_tm_params.h"
#include "sym_cp_params.h"
#include "sym_cg_params.h"
#include "sym_lp_params.h"
#include "sym_dg_params.h"
#include "sym_prep_params.h"

/*===========================================================================*\
 * The params structure contains all of the user-specified parameters
 * to be read in from the parameter file. See the README file for an
 * explanation of the parameters
\*===========================================================================*/

typedef struct PARAMS{
   int        warm_start;
   int        verbosity;
   char       param_file[MAX_FILE_NAME_LENGTH +1];
   int        random_seed;
   cp_params  cp_par;
   cg_params  cg_par;
   lp_params  lp_par;
   tm_params  tm_par;
   dg_params  dg_par;
   prep_params prep_par;

   char       tm_exe[MAX_FILE_NAME_LENGTH +1];
   char       dg_exe[MAX_FILE_NAME_LENGTH +1];
   int        tm_debug;
   int        dg_debug;
   int        tm_machine_set;
   char       tm_machine[MACH_NAME_LENGTH +1];
   int        dg_machine_set;
   char       dg_machine[MACH_NAME_LENGTH +1];
   int        pvm_trace;
   int        do_branch_and_cut;
   int        do_draw_graph;
   int        use_permanent_cut_pools;
   char       infile[MAX_FILE_NAME_LENGTH +1]; /* For MPS file name
						  or GNUMP modelfile */
   int        file_type;
   char       datafile[MAX_FILE_NAME_LENGTH +1]; /* GNUMP datafile */
   char       test_dir[MAX_FILE_NAME_LENGTH +1];  /* Test files directory */
   double     obj_offset;

   /* Multi-criteria parameters */
   int        multi_criteria;
   int        mc_search_order;
   double     mc_compare_solution_tolerance;
   double     mc_binary_search_tolerance;
   int        mc_warm_start;
   int        mc_warm_start_rule;
   int        trim_warm_tree;
   int        test;
}params;

#endif
