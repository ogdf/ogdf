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

#define COMPILING_FOR_TM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if !defined(_MSC_VER) && !defined(__MNO_CYGWIN) && defined(SIGHANDLER)
#include <signal.h>
#if !defined(HAS_SRANDOM)
extern int srandom PROTO((unsigned seed));
#endif
#if !defined(HAS_RANDOM)
extern long random PROTO((void));
#endif
#endif
#ifdef _OPENMP
#include "omp.h"
#endif
#if !defined (_MSC_VER)
#include <unistd.h>            /* this defines sleep() */
#endif

#include "sym_tm.h"
#include "sym_constants.h"
#include "sym_types.h"
#include "sym_macros.h"
#include "sym_messages.h"
#include "sym_proccomm.h"
#include "sym_timemeas.h"
#include "sym_pack_cut.h"
#include "sym_pack_array.h"
#ifdef COMPILE_IN_LP
#include "sym_lp.h"
#endif
#ifdef COMPILE_IN_TM
#include "sym_master.h"
#else
#include "sym_cp.h"
#endif

int c_count = 0;

/*===========================================================================*/

/*===========================================================================*\
 * This file contains basic functions associated with the treemanager process
\*===========================================================================*/

/*===========================================================================*/

/*===========================================================================*\
 * This function receives the intitial parameters and data and sets up
 * the tree manager data structures, etc.
\*===========================================================================*/

int tm_initialize(tm_prob *tm, base_desc *base, node_desc *rootdesc)
{
#ifndef COMPILE_IN_TM
   int r_bufid, bytes, msgtag, i;
#endif
   FILE *f = NULL;
   tm_params *par;
   bc_node *root = (bc_node *) calloc(1, sizeof(bc_node));
#ifdef COMPILE_IN_LP
   int i;
#else
#ifdef COMPILE_IN_TM
   int i;
#endif
   int s_bufid;
#endif
   int *termcodes = NULL;
#if !defined(_MSC_VER) && !defined(__MNO_CYGWIN) && defined(SIGHANDLER)
   signal(SIGINT, sym_catch_c);
#endif
   par = &tm->par;

#ifdef _OPENMP
   tm->rpath =
      (bc_node ***) calloc(par->max_active_nodes, sizeof(bc_node **));
   tm->rpath_size = (int *) calloc(par->max_active_nodes, sizeof(int));
   tm->bpath =
      (branch_desc **) calloc(par->max_active_nodes, sizeof(branch_desc *));
   tm->bpath_size = (int *) calloc(par->max_active_nodes, sizeof(int));
   termcodes = (int *) calloc(par->max_active_nodes, sizeof(int));
#else
   tm->rpath = (bc_node ***) calloc(1, sizeof(bc_node **));
   tm->rpath_size = (int *) calloc(1, sizeof(int));
   tm->bpath = (branch_desc **) calloc(1, sizeof(branch_desc *));
   tm->bpath_size = (int *) calloc(1, sizeof(int));
   termcodes = (int *) calloc(1, sizeof(int));
#endif

   /*------------------------------------------------------------------------*\
    * Receives data from the master
   \*------------------------------------------------------------------------*/

#ifdef COMPILE_IN_TM
   tm->bvarnum = base->varnum;
   tm->bcutnum = base->cutnum;
#else
   r_bufid = receive_msg(ANYONE, TM_DATA);
   bufinfo(r_bufid, &bytes, &msgtag, &tm->master);
   receive_char_array((char *)par, sizeof(tm_params));
   receive_char_array(&tm->has_ub, 1);
   if (tm->has_ub)
      receive_dbl_array(&tm->ub, 1);
   receive_char_array(&tm->has_ub_estimate, 1);
   if (tm->has_ub_estimate)
      receive_dbl_array(&tm->ub_estimate, 1);
   READ_STR_LIST(par->lp_mach_num, MACH_NAME_LENGTH,
		     par->lp_machs[0], par->lp_machs);
   READ_STR_LIST(par->cg_mach_num, MACH_NAME_LENGTH,
		     par->cg_machs[0], par->cg_machs);
   READ_STR_LIST(par->cp_mach_num, MACH_NAME_LENGTH,
		     par->cp_machs[0], par->cp_machs);
   receive_int_array(&tm->bvarnum, 1);
   receive_int_array(&tm->bcutnum, 1);
#ifdef TRACE_PATH
   receive_int_array(&tm->feas_sol_size, 1);
   if (tm->feas_sol_size){
      tm->feas_sol = (int *) calloc (tm->feas_sol_size, sizeof(int));
      receive_int_array(tm->feas_sol, tm->feas_sol_size);
   }
#endif
   freebuf(r_bufid);
#endif

   SRANDOM(par->random_seed);

#ifdef COMPILE_IN_LP
#ifdef _OPENMP
   omp_set_dynamic(FALSE);
   omp_set_num_threads(par->max_active_nodes);
#else
   par->max_active_nodes = 1;
#endif
   tm->active_nodes = (bc_node **) calloc(par->max_active_nodes, sizeof(bc_node *));
#ifndef COMPILE_IN_TM
   tm->lpp = (lp_prob **)
      malloc(par->max_active_nodes * sizeof(lp_prob *));
   for (i = 0; i < par->max_active_nodes; i++){
      tm->lpp[i] = (lp_prob *) calloc(1, sizeof(lp_prob));
      tm->lpp[i]->proc_index = i;
   }
#ifdef COMPILE_IN_CG
   tm->cgp = (cg_prob **) malloc(par->max_active_nodes * sizeof(cg_prob *));
   for (i = 0; i < par->max_active_nodes; i++)
      tm->lpp[i]->cgp = tm->cgp[i] = (cg_prob *) calloc(1, sizeof(cg_prob));
   par->use_cg = FALSE;
#endif
#endif
#pragma omp parallel for shared(tm)
   for (i = 0; i < par->max_active_nodes; i++){
      if ((termcodes[i] = lp_initialize(tm->lpp[i], 0)) < 0){
	 printf("LP initialization failed with error code %i in thread %i\n\n",
		termcodes[i], i);
      }
      tm->lpp[i]->tm = tm;
   }
   tm->lp.free_num = par->max_active_nodes;
   for (i = 0; i < par->max_active_nodes; i++){
      if (termcodes[i] < 0){
	 int tmp = termcodes[i];
	 FREE(termcodes);
	 return(tmp);
      }
   }
#else
   tm->active_nodes =
      (bc_node **) malloc(par->max_active_nodes * sizeof(bc_node *));

   /*------------------------------------------------------------------------*\
    * Start the lp, cg processes and send cg tid's to the lp's.
    * Also, start the cp, sp processes.
   \*------------------------------------------------------------------------*/

   tm->lp = start_processes(tm, par->max_active_nodes, par->lp_exe,
			    par->lp_debug, par->lp_mach_num, par->lp_machs);
#endif

#pragma omp critical (cut_pool)
   if (!tm->cuts){
      tm->cuts = (cut_data **) malloc(BB_BUNCH * sizeof(cut_data *));
   }

   if (par->use_cg){
#ifndef COMPILE_IN_CG
      tm->cg = start_processes(tm, par->max_active_nodes, par->cg_exe,
			       par->cg_debug, par->cg_mach_num, par->cg_machs);

#ifdef COMPILE_IN_LP
      for (i = 0; i < par->max_active_nodes; i++)
	 tm->lpp[i]->cut_gen = tm->cg.procs[i];
#else
      for (i = 0; i < tm->lp.procnum; i++){
	 s_bufid = init_send(DataInPlace);
	 send_int_array(tm->cg.procs + i, 1);
	 send_msg(tm->lp.procs[i], LP__CG_TID_INFO);
      }
#endif
#endif
   }

   if (par->max_cp_num){
#ifdef COMPILE_IN_CP
#ifndef COMPILE_IN_TM
      tm->cpp = (cut_pool **) malloc(par->max_cp_num * sizeof(cut_pool *));
#endif
      for (i = 0; i < par->max_cp_num; i++){
#ifndef COMPILE_IN_TM
	 tm->cpp[i] = (cut_pool *) calloc(1, sizeof(cut_pool));
#endif
	 cp_initialize(tm->cpp[i], tm->master);
      }
      tm->cp.free_num = par->max_cp_num;
      tm->cp.procnum = par->max_cp_num;
      tm->cp.free_ind = (int *) malloc(par->max_cp_num * ISIZE);
      for (i = par->max_cp_num - 1; i >= 0; i--)
	 tm->cp.free_ind[i] = i;
#else
      tm->cp = start_processes(tm, par->max_cp_num, par->cp_exe,
			      par->cp_debug, par->cp_mach_num, par->cp_machs);
#endif
      tm->nodes_per_cp = (int *) calloc(tm->par.max_cp_num, ISIZE);
      tm->active_nodes_per_cp = (int *) calloc(tm->par.max_cp_num, ISIZE);
   }else{
#ifdef COMPILE_IN_CP
      tm->cpp = (cut_pool **) calloc(1, sizeof(cut_pool *));
#endif
   }

   /*------------------------------------------------------------------------*\
    * Receive the root node and send out initial data to the LP processes
   \*------------------------------------------------------------------------*/

   FREE(termcodes);
   if (tm->par.warm_start){
      if (!tm->rootnode){
	 if (!(f = fopen(tm->par.warm_start_tree_file_name, "r"))){
	    printf("Error reading warm start file %s\n\n",
		   tm->par.warm_start_tree_file_name);
	    return(ERROR__READING_WARM_START_FILE);
	 }
	 read_tm_info(tm, f);
      }else{
	 free(root);
	 root = tm->rootnode;
      }
      read_subtree(tm, root, f);
      if (f)
	 fclose(f);
      if (!tm->rootnode){
	 if (!read_tm_cut_list(tm, tm->par.warm_start_cut_file_name)){
	    printf("Error reading warm start file %s\n\n",
		   tm->par.warm_start_cut_file_name);
	    return(ERROR__READING_WARM_START_FILE);
	 }
      }
      tm->rootnode = root;
      if(root->node_status != NODE_STATUS__WARM_STARTED)
	root->node_status = NODE_STATUS__ROOT;
   }else{
#ifdef COMPILE_IN_TM
      (tm->rootnode = root)->desc = *rootdesc;
      /* Copy the root description in case it is still needed */
      root->desc.uind.list = (int *) malloc(rootdesc->uind.size*ISIZE);
      memcpy((char *)root->desc.uind.list, (char *)rootdesc->uind.list,
	     rootdesc->uind.size*ISIZE);
      root->bc_index = tm->stat.created++;
      root->lower_bound = -MAXDOUBLE;
      tm->stat.tree_size++;
      insert_new_node(tm, root);
      tm->phase = 0;
      tm->lb = 0;
#else
      r_bufid = receive_msg(tm->master, TM_ROOT_DESCRIPTION);
      receive_node_desc(tm, root);
      if (root->desc.cutind.size > 0){ /* Hey we got cuts, too! Unpack them. */
	 unpack_cut_set(tm, 0, 0, NULL);
      }
      freebuf(r_bufid);
#endif
#ifdef TRACE_PATH
      root->optimal_path = TRUE;
#endif
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the main loop that solves the problem
\*===========================================================================*/

int solve(tm_prob *tm)
{
#ifndef COMPILE_IN_LP
   int r_bufid;
#endif
   int termcode = 0;
   double start_time = tm->start_time;
   double no_work_start, ramp_up_tm = 0, ramp_down_time = 0;
   char ramp_down = FALSE, ramp_up = TRUE;
   double then, then2, then3, now;
   double timeout2 = 30, timeout3 = tm->par.logging_interval, timeout4 = 10;

   /*------------------------------------------------------------------------*\
    * The Main Loop
   \*------------------------------------------------------------------------*/

   no_work_start = wall_clock(NULL);

   termcode = TM_UNFINISHED;
   for (; tm->phase <= 1; tm->phase++){
      if (tm->phase == 1 && !tm->par.warm_start){
	 if ((termcode = tasks_before_phase_two(tm)) ==
	     FUNCTION_TERMINATED_NORMALLY){
	    termcode = TM_FINISHED; /* Continue normally */
	 }
      }
      then  = wall_clock(NULL);
      then2 = wall_clock(NULL);
      then3 = wall_clock(NULL);
#pragma omp parallel default(shared)
{
#ifdef _OPENMP
      int i, thread_num = omp_get_thread_num();
#else
      int i, thread_num = 0;
#endif
      while (tm->active_node_num > 0 || tm->samephase_candnum > 0){
	 /*------------------------------------------------------------------*\
	  * while there are nodes being processed or while there are nodes
	  * waiting to be processed, continue to execute this loop
	 \*------------------------------------------------------------------*/
	 i = NEW_NODE__STARTED;
	 while (tm->lp.free_num > 0 && (tm->par.time_limit >= 0.0 ?
		(wall_clock(NULL) - start_time < tm->par.time_limit) : TRUE) &&
		(tm->par.node_limit >= 0 ?
		tm->stat.analyzed < tm->par.node_limit : TRUE) &&
		((tm->has_ub && (tm->par.gap_limit >= 0.0)) ?
		 fabs(100*(tm->ub-tm->lb)/tm->ub) > tm->par.gap_limit : TRUE)
		&& !(tm->par.find_first_feasible && tm->has_ub) && c_count <= 0){
	    if (tm->samephase_candnum > 0){
#pragma omp critical (tree_update)
	       i = start_node(tm, thread_num);
	    }else{
	       i = NEW_NODE__NONE;
	    }

	    if (i != NEW_NODE__STARTED)
	       break;

	    if (ramp_up){
	       ramp_up_tm += (wall_clock(NULL) -
				no_work_start) * (tm->lp.free_num + 1);
	    }
	    if (ramp_down){
	       ramp_down_time += (wall_clock(NULL) -
				  no_work_start) * (tm->lp.free_num + 1);
	    }

	    if (!tm->lp.free_num){
	       ramp_down = FALSE;
	       ramp_up = FALSE;
	    }else if (ramp_up){
	       no_work_start = wall_clock(NULL);
	    }else{
	       ramp_down = TRUE;
	       no_work_start = wall_clock(NULL);
	    }
#ifdef COMPILE_IN_LP
#ifdef _OPENMP
	    if (tm->par.verbosity > 0)
	       printf("Thread %i now processing node %i\n", thread_num,
		      tm->lpp[thread_num]->bc_index);
#endif

	    if(tm->par.node_selection_rule == DEPTH_FIRST_THEN_BEST_FIRST &&
	       tm->has_ub){
	      tm->par.node_selection_rule = LOWEST_LP_FIRST;
	    }

	    switch(process_chain(tm->lpp[thread_num])){

	     case FUNCTION_TERMINATED_NORMALLY:
	       break;

	     case ERROR__NO_BRANCHING_CANDIDATE:
	       termcode = TM_ERROR__NO_BRANCHING_CANDIDATE;
	       break;

	     case ERROR__ILLEGAL_RETURN_CODE:
	       termcode = TM_ERROR__ILLEGAL_RETURN_CODE;
	       break;

	     case ERROR__NUMERICAL_INSTABILITY:
	       termcode = TM_ERROR__NUMERICAL_INSTABILITY;
	       break;

	     case ERROR__COMM_ERROR:
	       termcode = TM_ERROR__COMM_ERROR;

	     case ERROR__USER:
	       termcode = TM_ERROR__USER;
	       break;

	     case ERROR__DUAL_INFEASIBLE:
	       if(tm->lpp[thread_num]->bc_index < 1 ) {
		  termcode = TM_UNBOUNDED;
	       }else{
		  termcode = TM_ERROR__NUMERICAL_INSTABILITY;
	       }
	       break;
	    }
#endif
#pragma omp master
{
            now = wall_clock(NULL);
	    if (now - then2 > timeout2){
	       if(tm->par.verbosity >= -1 ){
		  print_tree_status(tm);
	       }
	       then2 = now;
	    }
	    if (now - then3 > timeout3){
	       write_log_files(tm);
	       then3 = now;
	    }
}
	 }

	 if (c_count > 0){
	    termcode = TM_SIGNAL_CAUGHT;
	    c_count = 0;
	    break;
	 }

	 if (tm->par.time_limit >= 0.0 &&
	     wall_clock(NULL) - start_time > tm->par.time_limit &&
	     termcode != TM_FINISHED){
	    termcode = TM_TIME_LIMIT_EXCEEDED;
	    break;
	 }

	 if (tm->par.node_limit >= 0 && tm->stat.analyzed >=
	     tm->par.node_limit && termcode != TM_FINISHED){
	    if (tm->active_node_num + tm->samephase_candnum > 0){
	       termcode = TM_NODE_LIMIT_EXCEEDED;
	    }else{
	       termcode = TM_FINISHED;
	    }
	    break;
	 }

	 if (tm->par.find_first_feasible && tm->has_ub){
	    termcode = TM_FINISHED;
	    break;
	 }

	 if (i == NEW_NODE__ERROR){
	    termcode = SOMETHING_DIED;
	    break;
	 }

	 if (tm->has_ub && (tm->par.gap_limit >= 0.0)){
            find_tree_lb(tm);
	    if (fabs(100*(tm->ub-tm->lb)/tm->ub) <= tm->par.gap_limit){
	       if (tm->lb < tm->ub){
		  termcode = TM_TARGET_GAP_ACHIEVED;
	       }else{
		  termcode = TM_FINISHED;
	       }
	       break;
	    }
	 }
	 if (i == NEW_NODE__NONE && tm->active_node_num == 0)
	    break;

#ifndef COMPILE_IN_LP

	 struct timeval timeout = {5, 0};
	 r_bufid = treceive_msg(ANYONE, ANYTHING, &timeout);
	 if (r_bufid && !process_messages(tm, r_bufid)){
            find_tree_lb(tm);
	    termcode = SOMETHING_DIED;
	    break;
	 }
#endif
	 now = wall_clock(NULL);
	 if (now - then > timeout4){
	    if (!processes_alive(tm)){
               find_tree_lb(tm);
               termcode = SOMETHING_DIED;
	       break;
	    }
	    then = now;
	 }
#pragma omp master
{
         for (i = 0; i < tm->par.max_active_nodes; i++){
	    if (tm->active_nodes[i]){
	       break;
	    }
	 }
	 if (i == tm->par.max_active_nodes){
	    tm->active_node_num = 0;
	 }
	 if (now - then2 > timeout2){
	    if(tm->par.verbosity >=0 ){
	       print_tree_status(tm);
	    }
	    then2 = now;
	 }
	 if (now - then3 > timeout3){
	    write_log_files(tm);
	    then3 = now;
	 }
}
      }
}

      if(termcode == TM_UNBOUNDED) break;

      if (tm->samephase_candnum + tm->active_node_num == 0){
	 termcode = TM_FINISHED;
      }
      if (tm->nextphase_candnum == 0)
	 break;
      if (termcode != TM_UNFINISHED)
	 break;
   }
   find_tree_lb(tm);
   tm->comp_times.ramp_up_tm = ramp_up_tm;
   tm->comp_times.ramp_down_time = ramp_down_time;
   write_log_files(tm);

   return(termcode);
}

/*===========================================================================*/

/*==========================================================================*\
 * Write out the log files
\*==========================================================================*/

void write_log_files(tm_prob *tm)
{
#if !defined(COMPILE_IN_LP) || !defined(COMPILE_IN_CP)
   int s_bufid;
#endif

   if (tm->par.logging){
      write_tm_info(tm, tm->par.tree_log_file_name, NULL, FALSE);
      write_subtree(tm->rootnode, tm->par.tree_log_file_name, NULL, TRUE,
		    tm->par.logging);
      if (tm->par.logging != VBC_TOOL)
	 write_tm_cut_list(tm, tm->par.cut_log_file_name, FALSE);
   }

   if (tm->par.max_cp_num > 0 && tm->par.cp_logging){
#if defined(COMPILE_IN_LP) && defined(COMPILE_IN_CP)
      write_cp_cut_list(tm->cpp[0], tm->cpp[0]->par.log_file_name,
			FALSE);
#else
      s_bufid = init_send(DataInPlace);
      send_msg(tm->cp.procs[0], WRITE_LOG_FILE);
#endif
   }
}

/*===========================================================================*/

/*==========================================================================*\
 * Prints out the current size of the tree and the gap                      *
\*==========================================================================*/

void print_tree_status(tm_prob *tm)
{
   double elapsed_time;
   double obj_ub = SYM_INFINITY, obj_lb = -SYM_INFINITY;

#ifdef SHOULD_SHOW_MEMORY_USAGE
   int i;
   int pid;
   int tmp_int;
   long unsigned vsize;
   char tmp_str[100], proc_filename[100];
   FILE *proc_file;
   double vsize_in_mb;
#endif

#if 0
   int *widths;
   double *gamma;
   int last_full_level = 0, max_width = 0, num_nodes_estimate = 1;
   int first_waist_level = 0, last_waist_level = 0, waist_level = 0;
   double average_node_time, estimated_time_remaining, user_time = 0.0;

   widths = (int *) calloc (tm->stat.max_depth + 1, ISIZE);
   gamma = (double *) calloc (tm->stat.max_depth + 1, DSIZE);

   calculate_widths(tm->rootnode, widths);

   last_full_level = tm->stat.max_depth;
   for (i = tm->stat.max_depth - 1; i > 0; i--){
      if ((double)(widths[i])/(double)(widths[i - 1]) < 2){
	 last_full_level = i - 1;
      }
      if (widths[i] > max_width){
	 max_width = widths[i];
	 last_waist_level = i;
	 first_waist_level = i;
      }
      if (widths[i] == max_width){
	 first_waist_level = i;
      }
   }

   waist_level = (first_waist_level + last_waist_level)/2;

   for (i = 0; i < tm->stat.max_depth; i++){
      if (i < last_full_level){
	 gamma[i] = 2.0;
      }else if (i < waist_level){
	 gamma[i] = 2.0 - (double)((i - last_full_level + 1))/
	    (double)((waist_level - last_full_level + 1));
      }else{
	 gamma[i] = 1.0 - (double)(i - waist_level + 1)/
	    (double)(tm->stat.max_depth - waist_level + 1);
      }
   }

   for (i = 1; i < tm->stat.max_depth; i++){
      gamma[i] *= gamma[i - 1];
      num_nodes_estimate += (int)(gamma[i] + 0.5);
   }

   elapsed_time = wall_clock(NULL) - tm->start_time;
   average_node_time = elapsed_time/tm->stat.analyzed;

   estimated_time_remaining =
      MAX(average_node_time*(num_nodes_estimate - tm->stat.analyzed), 0);
#else

   elapsed_time = wall_clock(NULL) - tm->start_time;

#endif

#ifdef SHOULD_SHOW_MEMORY_USAGE
   pid = getpid();
   //printf("process id = %d\n",pid);
   sprintf(proc_filename,"/proc/%d/stat",pid);
   proc_file = fopen (proc_filename, "r");
   fscanf (proc_file, "%d %s %s", &tmp_int, tmp_str, tmp_str);
   for (i=0; i<19;i++) {
      fscanf (proc_file, "%d", &tmp_int);
   }
   fscanf (proc_file, "%lu", &vsize);
   fclose(proc_file);
   //printf("vsize = %lu\n",vsize);
   vsize_in_mb = vsize/1024.0/1024.0;
   if (tm->stat.max_vsize<vsize_in_mb) {
      tm->stat.max_vsize = vsize_in_mb;
   }
   printf("memory: %.2f MB ", vsize_in_mb);
#endif

   printf("done: %i ", tm->stat.analyzed-tm->active_node_num);
   printf("left: %i ", tm->samephase_candnum+tm->active_node_num);
   if (tm->has_ub) {
      if (tm->obj_sense == SYM_MAXIMIZE){
         obj_lb = -tm->ub + tm->obj_offset;
	 printf("lb: %.2f ", obj_lb);
      }else{
         obj_ub = tm->ub + tm->obj_offset;
	 printf("ub: %.2f ", obj_ub);
      }
   } else {
      if (tm->obj_sense == SYM_MAXIMIZE){
	 printf("lb: ?? ");
      }else{
	 printf("ub: ?? ");
      }
   }
   find_tree_lb(tm);
   if(tm->lb > -SYM_INFINITY){
      if (tm->obj_sense == SYM_MAXIMIZE){
	 obj_ub = -tm->lb + tm->obj_offset;
	 printf("ub: %.2f ", obj_ub);
      }else{
	 obj_lb = tm->lb + tm->obj_offset;
	 printf("lb: %.2f ", obj_lb);
      }
   }else{
      if (tm->obj_sense == SYM_MAXIMIZE){
	 printf("ub: ?? ");
      }else{
	 printf("lb: ?? ");
      }
   }
   if (tm->has_ub && tm->ub && tm->lb > -SYM_INFINITY){
      printf("gap: %.2f ", fabs(100*(obj_ub-obj_lb)/obj_ub));
   }
   printf("time: %i\n", (int)(elapsed_time));
#if 0
   printf("Estimated nodes remaining:         %i\n", num_nodes_estimate);
   printf("Estimated time remaining:          %i\n",
	  (int)(estimated_time_remaining));
#endif

   if (tm->par.vbc_emulation == VBC_EMULATION_FILE){
      FILE *f;
#pragma omp critical(write_vbc_emulation_file)
      if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 PRINT_TIME(tm, f);
	 fprintf(f, "L %.2f \n", tm->lb);
	 fclose(f);
      }
   }else if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$L %.2f\n", tm->lb);
   }

#if 0
   FREE(widths);
   FREE(gamma);
#endif
}

/*===========================================================================*/

void calculate_widths(bc_node *node, int* widths)
{
   int i;

   widths[node->bc_level] += 1;
   for (i = 0; i < node->bobj.child_num; i ++){
      calculate_widths(node->children[i], widths);
   }
}
/*===========================================================================*/

/*===========================================================================*\
 * This function picks the "best" node off the active node list
\*===========================================================================*/

int start_node(tm_prob *tm, int thread_num)
{
   int lp_ind, get_next, ind;
   bc_node *best_node = NULL;
   double time;

   time = wall_clock(NULL);

   /*------------------------------------------------------------------------*\
    * First choose the "best" node from the list of candidate nodes.
    * If the list for the current phase is empty then we return NEW_NODE__NONE.
    * Also, if the lower bound on the "best" node is above the current UB then
    * we just move that node the list of next phase candidates.
   \*------------------------------------------------------------------------*/

   get_next = TRUE;
   while (get_next){
      if ((best_node = del_best_node(tm)) == NULL)
	 return(NEW_NODE__NONE);

      if (best_node->node_status == NODE_STATUS__WARM_STARTED){
	 if(best_node->lower_bound >= MAXDOUBLE)
	    break;
      }

      /* if no UB yet or lb is lower than UB then go ahead */
      if (!tm->has_ub ||
	  (tm->has_ub && best_node->lower_bound < tm->ub-tm->par.granularity))
	 break;
      /* ok, so we do have an UB and lb is higher than the UB. */
      /* in this switch we assume that there are only two phases! */
      switch (((best_node->desc.nf_status) << 8) + tm->phase){
       case (NF_CHECK_NOTHING << 8) + 0: /* prune these */
       case (NF_CHECK_NOTHING << 8) + 1:
	  if(!tm->par.sensitivity_analysis){
	     if (tm->par.max_cp_num > 0 && best_node->cp){
#ifdef COMPILE_IN_CP
		ind = best_node->cp;
#else
		ind = find_process_index(&tm->cp, best_node->cp);
#endif
		tm->nodes_per_cp[ind]--;
		if (tm->nodes_per_cp[ind] + tm->active_nodes_per_cp[ind] == 0)
		   tm->cp.free_ind[tm->cp.free_num++] = ind;
	     }
	     best_node->node_status = NODE_STATUS__PRUNED;
	     best_node->feasibility_status = OVER_UB_PRUNED;

	     if (tm->par.verbosity > 0){
		printf("++++++++++++++++++++++++++++++++++++++++++++++++++\n");
		printf("+ TM: Pruning NODE %i LEVEL %i instead of sending it.\n",
		       best_node->bc_index, best_node->bc_level);
		printf("++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	     }
	     if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL ||
		 tm->par.keep_description_of_pruned == KEEP_ON_DISK_FULL ||
		 tm->par.keep_description_of_pruned == DISCARD){
		if (tm->par.keep_description_of_pruned ==
		    KEEP_ON_DISK_VBC_TOOL ||
		    tm->par.keep_description_of_pruned == KEEP_ON_DISK_FULL){
#pragma omp critical (write_pruned_node_file)
		   write_pruned_nodes(tm, best_node);
		}
#if 0
		if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW) {
		   purge_pruned_nodes(tm, best_node, VBC_PRUNED_FATHOMED);
		} else {
		   purge_pruned_nodes(tm, best_node, VBC_PRUNED);
		}
#else
		   purge_pruned_nodes(tm, best_node, VBC_PRUNED);
#endif
	     }
	     break;
	  }
       case (NF_CHECK_ALL            << 8) + 1: /* work on these */
       case (NF_CHECK_UNTIL_LAST     << 8) + 1:
       case (NF_CHECK_AFTER_LAST     << 8) + 1:
	 get_next = FALSE;
	 break;

       default:
	 /* i.e., phase == 0 and nf_status != NF_CHECK_NOTHING */
	  if (!(tm->par.colgen_strat[0] & FATHOM__GENERATE_COLS__RESOLVE)){
	     REALLOC(tm->nextphase_cand, bc_node *, tm->nextphase_cand_size,
		     tm->nextphase_candnum+1, BB_BUNCH);
	     tm->nextphase_cand[tm->nextphase_candnum++] = best_node;
	  }else{
	     get_next = FALSE;
	  }
	  break;
      }
   }

   /* Assign a free lp process */
#ifdef COMPILE_IN_LP
   lp_ind = thread_num;
#else
   lp_ind = tm->lp.free_ind[--tm->lp.free_num];
   best_node->lp = tm->lp.procs[lp_ind];
   best_node->cg = tm->par.use_cg ? tm->cg.procs[lp_ind] : 0;
#endif
   /* assign pools, too */

   best_node->cp = assign_pool(tm, best_node->cp, &tm->cp,
			       tm->active_nodes_per_cp, tm->nodes_per_cp);
   if (best_node->cp < 0) return(NEW_NODE__ERROR);


   /* It's time to put together the node and send it out */
   tm->active_nodes[lp_ind] = best_node;
   tm->active_node_num++;
   tm->stat.analyzed++;

   send_active_node(tm,best_node,tm->par.colgen_strat[tm->phase],thread_num);

   tm->comp_times.start_node += wall_clock(NULL) - time;

   return(NEW_NODE__STARTED);
}

/*===========================================================================*/

/*===========================================================================*\
 * Returns the "best" active node and deletes it from the list
\*===========================================================================*/

bc_node *del_best_node(tm_prob *tm)
{
   bc_node **list = tm->samephase_cand;
   int size = tm->samephase_candnum;
   bc_node *temp = NULL, *best_node;
   int pos, ch;
   int rule = tm->par.node_selection_rule;

   if (size == 0)
      return(NULL);

   best_node = list[1];

   temp = list[1] = list[size];

   tm->samephase_candnum = --size;

   if (tm->par.verbosity > 10)
      if (tm->samephase_candnum % 10 == 0)
	 printf("\nTM: tree size: %i , %i\n\n",
		tm->samephase_candnum, tm->nextphase_candnum);

   pos = 1;
   while ((ch=2*pos) < size){
      if (node_compar(rule, list[ch], list[ch+1]))
	 ch++;
      if (node_compar(rule, list[ch], temp)){
	 list[pos] = temp;
	 return(best_node);
      }
      list[pos] = list[ch];
      pos = ch;
   }
   if (ch == size){
      if (node_compar(rule, temp, list[ch])){
	 list[pos] = list[ch];
	 pos = ch;
      }
   }
   list[pos] = temp;
   return(best_node);
}

/*===========================================================================*/

/*===========================================================================*\
 * Insert a new active node into the active node list (kept as a binary tree)
\*===========================================================================*/

void insert_new_node(tm_prob *tm, bc_node *node)
{
   int pos, ch, size = tm->samephase_candnum;
   bc_node **list;
   int rule = tm->par.node_selection_rule;

   tm->samephase_candnum = pos = ++size;

   if (tm->par.verbosity > 10)
      if (tm->samephase_candnum % 10 == 0)
	 printf("\nTM: tree size: %i , %i\n\n",
		tm->samephase_candnum, tm->nextphase_candnum);

   REALLOC(tm->samephase_cand, bc_node *,
	   tm->samephase_cand_size, size + 1, BB_BUNCH);
   list = tm->samephase_cand;

   while ((ch=pos>>1) != 0){
      if (node_compar(rule, list[ch], node)){
	 list[pos] = list[ch];
	 pos = ch;
      }else{
	 break;
      }
   }
   list[pos] = node;
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the node comparison function used to order the list of active
 * Nodes are ordered differently depending on what the comparison rule is
\*===========================================================================*/

int node_compar(int rule, bc_node *node0, bc_node *node1)
{
   switch(rule){
    case LOWEST_LP_FIRST:
      return(node1->lower_bound < node0->lower_bound ? 1:0);
    case HIGHEST_LP_FIRST:
      return(node1->lower_bound > node0->lower_bound ? 1:0);
    case BREADTH_FIRST_SEARCH:
      return(node1->bc_level < node0->bc_level ? 1:0);
    case DEPTH_FIRST_SEARCH:
    case DEPTH_FIRST_THEN_BEST_FIRST:
      return(node1->bc_level > node0->bc_level ? 1:0);
   }
   return(0); /* fake return */
}

/*===========================================================================*/

/*===========================================================================*\
 * Nodes by default inherit their parent's pools. However if there is a free
 * pool then the node is moved over to the free pool.
\*===========================================================================*/

int assign_pool(tm_prob *tm, int oldpool, process_set *pools,
		int *active_nodes_per_pool, int *nodes_per_pool)
{
   int oldind = -1, ind, pool;
#ifndef COMPILE_IN_CP
   int s_bufid, r_bufid;
   struct timeval timeout = {5, 0};
#endif

   if (pools->free_num == 0){
      /* No change in the pool assigned to this node */
      return(oldpool);
   }

   if (oldpool > 0){
#ifdef COMPILE_IN_CP
      oldind = oldpool;
#else
      oldind = find_process_index(pools, oldpool);
#endif
      if (nodes_per_pool[oldind] == 1){
	 nodes_per_pool[oldind]--;
	 active_nodes_per_pool[oldind]++;
	 return(oldpool);
      }
   }

   ind = pools->free_ind[--pools->free_num];
#ifdef COMPILE_IN_CP
   pool = ind;
#else
   pool = pools->procs[ind];
#endif

   if (! oldpool){
      /* If no pool is assigned yet then just assign the free one */
      active_nodes_per_pool[ind] = 1;
      return(pool);
   }

   /* finally when we really move the node from one pool to another */
   nodes_per_pool[oldind]--;
   active_nodes_per_pool[ind] = 1;
#ifdef COMPILE_IN_CP
   /*FIXME: Multiple Pools won't work in shared memory mode until I fill this
     in.*/
#else
   s_bufid = init_send(DataInPlace);
   send_int_array(&oldpool, 1);
   send_msg(pool, POOL_YOU_ARE_USELESS);
   s_bufid = init_send(DataInPlace);
   send_int_array(&pool, 1);
   send_msg(oldpool, POOL_COPY_YOURSELF);
   freebuf(s_bufid);

   do{
      r_bufid = treceive_msg(pool, POOL_USELESSNESS_ACKNOWLEDGED, &timeout);
      if (r_bufid == 0)
	 if (pstat(pool) != PROCESS_OK) return(NEW_NODE__ERROR);
   }while (r_bufid == 0);
   freebuf(r_bufid);
#endif

   return(pool);
}

/*===========================================================================*/

/*===========================================================================*\
 * Takes the branching object description and sets up data structures
 * for the resulting children and adds them to the list of candidates.
\*===========================================================================*/

int generate_children(tm_prob *tm, bc_node *node, branch_obj *bobj,
		      double *objval, int *feasible, char *action,
		      int olddive, int *keep, int new_branching_cut)
{
   node_desc *desc;
   int np_cp = 0, np_sp = 0;
   int dive = DO_NOT_DIVE, i;
   bc_node *child;
   int child_num;
#ifdef TRACE_PATH
   int optimal_path = -1;
#endif

   /* before we start to generate the children we must figure out if we'll
    * dive so that we can put the kept child into the right location */
   if (*keep >= 0 && (olddive == CHECK_BEFORE_DIVE || olddive == DO_DIVE))
      dive = olddive == DO_DIVE ? DO_DIVE : shall_we_dive(tm, objval[*keep]);

   node->children = (bc_node **) calloc(bobj->child_num, sizeof(bc_node *));
   if (node->bc_level == tm->stat.max_depth)
      tm->stat.max_depth++;
   child_num = bobj->child_num;
#ifdef TRACE_PATH
   if (node->optimal_path && tm->feas_sol_size){
      for (i = 0; i < tm->feas_sol_size; i++)
	 if (tm->feas_sol[i] == bobj->name)
	    break;
      if (i < tm->feas_sol_size)
	 optimal_path = 1;
      else
	 optimal_path = 0;
      printf("\n\nNode %i is on the optimal path\n\n",
	     tm->stat.tree_size + optimal_path);
   }
#endif
   for (i = 0; i < child_num; i++){
      child = node->children[i] = (bc_node *) calloc(1, sizeof(bc_node));
      child->bc_index = tm->stat.tree_size++;
      child->bc_level = node->bc_level + 1;
      child->lower_bound = objval[i];
#ifdef COMPILE_IN_LP
      child->update_pc = bobj->is_est[i] ? TRUE : FALSE;
#endif
      child->parent = node;
      if (tm->par.verbosity > 10){
	 printf("Generating node %i from %i...\n", child->bc_index,
		node->bc_index);
      }
      if (tm->par.vbc_emulation == VBC_EMULATION_FILE){
	 FILE *f;
#pragma omp critical(write_vbc_emulation_file)
	 if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	    printf("\nError opening vbc emulation file\n\n");
	 }else{
	    PRINT_TIME(tm, f);
	    fprintf(f, "N %i %i %i\n", node->bc_index+1, child->bc_index+1,
		  feasible[i] ? VBC_FEAS_SOL_FOUND :
		  ((dive != DO_NOT_DIVE && *keep == i) ?
		   VBC_ACTIVE_NODE : VBC_CAND_NODE));
            fclose(f);
	 }
      } else if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW) {
	 FILE *f;
#pragma omp critical(write_vbc_emulation_file)
	 if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	    printf("\nError opening vbc emulation file\n\n");
	 }else{
	    PRINT_TIME2(tm, f);
	    char reason[50];
	    char branch_dir = 'M';
	    sprintf (reason, "%s %i %i", "candidate", child->bc_index+1,
		     node->bc_index+1);
	    if (child->bc_index>0){
	       if (node->children[0]==child) {
		  branch_dir = node->bobj.sense[0];
		  /*branch_dir = 'L';*/
	       } else {
		  branch_dir = node->bobj.sense[1];
		  /*branch_dir = 'R';*/
	       }
	       if (branch_dir == 'G') {
		  branch_dir = 'R';
	       }
	    }
	    if (action[i] == PRUNE_THIS_CHILD_FATHOMABLE ||
		action[i] == PRUNE_THIS_CHILD_INFEASIBLE){
	       sprintf(reason,"%s %c", reason, branch_dir);
	    }else{
	       sprintf(reason,"%s %c %f", reason, branch_dir,
		       child->lower_bound);
	    }
	    fprintf(f,"%s\n",reason);
	    fclose(f);
	 }
      }else if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
	 printf("$N %i %i %i\n", node->bc_index+1, child->bc_index+1,
		feasible[i] ? VBC_FEAS_SOL_FOUND :
		((dive != DO_NOT_DIVE && *keep == i) ?
		 VBC_ACTIVE_NODE: VBC_CAND_NODE));
      }
#ifdef TRACE_PATH
      if (optimal_path == i)
	 child->optimal_path = TRUE;
#endif
      tm->stat.created++;
#ifndef ROOT_NODE_ONLY
      if (action[i] == PRUNE_THIS_CHILD ||
	  action[i] == PRUNE_THIS_CHILD_FATHOMABLE ||
	  action[i] == PRUNE_THIS_CHILD_INFEASIBLE ||
	  (tm->has_ub && tm->ub - tm->par.granularity < objval[i] &&
	   node->desc.nf_status == NF_CHECK_NOTHING)){
	 /* this last can happen if the TM got the new bound but it hasn't
	   * been propagated to the LP yet */
#else /*We only want to process the root node in this case - discard others*/
      if (TRUE){
#endif
	 if (tm->par.verbosity > 0){
	    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	    printf("+ TM: Pruning NODE %i LEVEL %i while generating it.\n",
		   child->bc_index, child->bc_level);
	    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	 }

	 child->node_status = NODE_STATUS__PRUNED;
#ifdef TRACE_PATH
	 if (child->optimal_path){
	    printf("\n\nAttempting to prune the optimal path!!!!!!!!!\n\n");
	    sleep(600);
	    if (tm->par.logging){
	       write_tm_info(tm, tm->par.tree_log_file_name, NULL, FALSE);
	       write_subtree(tm->rootnode, tm->par.tree_log_file_name, NULL,
			     TRUE, tm->par.logging);
	       write_tm_cut_list(tm, tm->par.cut_log_file_name, FALSE);
	    }
	    exit(1);
	 }
#endif
	 if (tm->par.keep_description_of_pruned == DISCARD ||
	     tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
	    child->parent = node;
	    if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL)
#pragma omp critical (write_pruned_node_file)
	       write_pruned_nodes(tm, child);
	    if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW) {
	       int vbc_node_pr_reason;
	       switch (action[i]) {
		case PRUNE_THIS_CHILD_INFEASIBLE:
		  vbc_node_pr_reason = VBC_PRUNED_INFEASIBLE;
		  break;
		case PRUNE_THIS_CHILD_FATHOMABLE:
		  vbc_node_pr_reason = VBC_PRUNED_FATHOMED;
		  break;
		default:
		  vbc_node_pr_reason = VBC_PRUNED;
	       }
	       /* following is no longer needed because this care is taken
		* care of in install_new_ub
		*/
	       /*
	       if (feasible[i]) {
		  vbc_node_pr_reason = VBC_FEAS_SOL_FOUND;
	       }
	       */
#pragma omp critical (tree_update)
	       purge_pruned_nodes(tm, child, vbc_node_pr_reason);
	    } else {
#pragma omp critical (tree_update)
	       purge_pruned_nodes(tm, child, feasible[i] ? VBC_FEAS_SOL_FOUND :
		     VBC_PRUNED);
	    }

	    if (--child_num == 0){
	       *keep = -1;
	       return(DO_NOT_DIVE);
	    }
	    if (*keep == child_num) *keep = i;
#ifdef TRACE_PATH
	    if (optimal_path == child_num) optimal_path = i;
#endif
	    action[i] = action[child_num];
	    objval[i] = objval[child_num];
	    feasible[i--] = feasible[child_num];
	    continue;
	 }
      }else{
	 child->node_status = NODE_STATUS__CANDIDATE;
	 /* child->lp = child->cg = 0;   zeroed out by calloc */
	 child->cp = node->cp;
      }
#ifdef DO_TESTS
      if (child->lower_bound < child->parent->lower_bound - .01){
	 printf("#######Error: Child's lower bound (%.3f) is less than ",
		child->lower_bound);
	 printf("parent's (%.3f)\n", child->parent->lower_bound);
      }
      if (child->lower_bound < tm->rootnode->lower_bound - .01){
	 printf("#######Error: Node's lower bound (%.3f) is less than ",
		child->lower_bound);
	 printf("root's (%.3f)\n", tm->rootnode->lower_bound);
      }
#endif

      /* child->children = NULL;   zeroed out by calloc */
      /* child->child_num = 0;   zeroed out by calloc */
      /* child->died = 0;   zeroed out by calloc */
      desc = &child->desc;
      /* all this is set by calloc
       * desc->uind.type = 0;            WRT_PARENT and no change
       * desc->uind.size = 0;
       * desc->uind.added = 0;
       * desc->uind.list = NULL;

       * desc->not_fixed.type = 0;       WRT_PARENT and no change
       * desc->not_fixed.size = 0;
       * desc->not_fixed.added = 0;
       * desc->not_fixed.list = NULL;

       * desc->cutind.type = 0;          WRT_PARENT and no change
       * desc->cutind.size = 0;
       * desc->cutind.added = 0;
       * desc->cutind.list = NULL;

       * desc->basis.basis_exists = FALSE;    This has to be validated!!!
       * desc->basis.{[base,extra][rows,vars]}
                    .type = 0;           WRT_PARENT and no change
                    .size = 0;
		    .list = NULL;
		    .stat = NULL;
       */

      if (node->desc.basis.basis_exists){
	 desc->basis.basis_exists = TRUE;
      }

      /* If we have a non-base, new branching cut then few more things
         might have to be fixed */
      if (new_branching_cut && bobj->name >= 0){
	 /* Fix cutind and the basis description */
	 desc->cutind.size = 1;
	 desc->cutind.added = 1;
	 desc->cutind.list = (int *) malloc(ISIZE);
	 desc->cutind.list[0] = bobj->name;
	 if (desc->basis.basis_exists){
	    desc->basis.extrarows.size = 1;
	    desc->basis.extrarows.list = (int *) malloc(ISIZE);
	    desc->basis.extrarows.list[0] = bobj->name;
	    desc->basis.extrarows.stat = (int *) malloc(ISIZE);
	    desc->basis.extrarows.stat[0] = SLACK_BASIC;
	 }
      }

      desc->desc_size = node->desc.desc_size;
      desc->desc = node->desc.desc;
      desc->nf_status = node->desc.nf_status;


#ifdef SENSITIVITY_ANALYSIS
      if (tm->par.sensitivity_analysis &&
	  action[i] != PRUNE_THIS_CHILD_INFEASIBLE){
	 child->duals = bobj->duals[i];
	 bobj->duals[i] = 0;
      }
#endif

      if (child->node_status != NODE_STATUS__PRUNED && feasible[i]){
	 if(tm->par.keep_description_of_pruned == KEEP_IN_MEMORY){
	    child->sol_size = bobj->sol_sizes[i];
	    child->sol_ind = bobj->sol_inds[i];
	    bobj->sol_inds[i]=0;
	    child->sol = bobj->solutions[i];
	    bobj->solutions[i] = 0;
	    child->feasibility_status = NOT_PRUNED_HAS_CAN_SOLUTION;
	 }
      }

      if (child->node_status == NODE_STATUS__PRUNED){

	 if(tm->par.keep_description_of_pruned == KEEP_IN_MEMORY){

	    child->feasibility_status = OVER_UB_PRUNED;

	    if (feasible[i]){
	       child->sol_size = bobj->sol_sizes[i];
	       child->sol_ind = bobj->sol_inds[i];
	       bobj->sol_inds[i] = 0;
	       child->sol = bobj->solutions[i];
	       bobj->solutions[i] = 0;
	       child->feasibility_status = FEASIBLE_PRUNED;
	    }

	    if (action[i] == PRUNE_THIS_CHILD_INFEASIBLE){
	       child->feasibility_status = INFEASIBLE_PRUNED;
	    }
	 }

#ifdef TRACE_PATH
	 if (child->optimal_path){
	    printf("\n\nAttempting to prune the optimal path!!!!!!!!!\n\n");
	    sleep(600);
	    if (tm->par.logging){
	       write_tm_info(tm, tm->par.tree_log_file_name, NULL, FALSE);
	       write_subtree(tm->rootnode, tm->par.tree_log_file_name, NULL,
			     TRUE, tm->par.logging);
	       write_tm_cut_list(tm, tm->par.cut_log_file_name, FALSE);
	    }
	    exit(1);
	 }
#endif
	 if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_FULL ||
	     tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
#pragma omp critical (write_pruned_node_file)
	    write_pruned_nodes(tm, child);
#pragma omp critical (tree_update)
	    if (tm->par.vbc_emulation== VBC_EMULATION_FILE_NEW) {
	       int vbc_node_pr_reason;
	       switch (action[i]) {
		case PRUNE_THIS_CHILD_INFEASIBLE:
		  vbc_node_pr_reason = VBC_PRUNED_INFEASIBLE;
		  break;
		case PRUNE_THIS_CHILD_FATHOMABLE:
		  vbc_node_pr_reason = VBC_PRUNED_FATHOMED;
		  break;
		default:
		  vbc_node_pr_reason = VBC_PRUNED;
	       }
	       /* following is no longer needed because this care is taken
		* care of in install_new_ub
		*/
	       /*
		  if (feasible[i]) {
		  vbc_node_pr_reason = VBC_FEAS_SOL_FOUND;
	       }
	       */
	       purge_pruned_nodes(tm, child, vbc_node_pr_reason);
	    } else {
	       purge_pruned_nodes(tm, child, feasible[i] ? VBC_FEAS_SOL_FOUND :
		     VBC_PRUNED);
	    }

	    if (--child_num == 0){
	       *keep = -1;
	       return(DO_NOT_DIVE);
	    }
	    if (*keep == child_num) *keep = i;
#ifdef TRACE_PATH
	    if (optimal_path == child_num) optimal_path = i;
#endif
	    action[i] = action[child_num];
	    objval[i] = objval[child_num];
	    feasible[i--] = feasible[child_num];
	 }
	 continue;
      }

      if (tm->phase == 0 &&
	  !(tm->par.colgen_strat[0] & FATHOM__GENERATE_COLS__RESOLVE) &&
	  (feasible[i] == LP_D_UNBOUNDED ||
	   (tm->has_ub && tm->ub - tm->par.granularity < child->lower_bound))){
	 /* it is kept for the next phase  (==> do not dive) */
	 if (*keep == i)
	    dive = DO_NOT_DIVE;
	 REALLOC(tm->nextphase_cand, bc_node *,
                 tm->nextphase_cand_size, tm->nextphase_candnum+1, BB_BUNCH);
	 tm->nextphase_cand[tm->nextphase_candnum++] = child;
	 np_cp++;
	 np_sp++;
      }else{
	 /* it will be processed in this phase (==> insert it if not kept) */
	 if (*keep != i || dive == DO_NOT_DIVE){
#pragma omp critical (tree_update)
	    insert_new_node(tm, child);
	    np_cp++;
	    np_sp++;
	 }
      }
   }
   if (node->cp)
#ifdef COMPILE_IN_CP
      tm->nodes_per_cp[node->cp] += np_cp;
#else
      tm->nodes_per_cp[find_process_index(&tm->cp, node->cp)] += np_cp;
#endif

   return(dive);
}

/*===========================================================================*/

/*===========================================================================*\
 * Determines whether or not the LP process should keep one of the
 * children resulting from branching or whether it should get a new node
 * from the candidate list.
\*===========================================================================*/

char shall_we_dive(tm_prob *tm, double objval)
{
   char dive;
   int i, k;
   double rand_num, average_lb;
   double cutoff = 0;
   double etol = 1e-3;

   if (tm->par.time_limit >= 0.0 &&
	wall_clock(NULL) - tm->start_time >= tm->par.time_limit){
      return(FALSE);
   }

   if (tm->par.node_limit >= 0 && tm->stat.analyzed >= tm->par.node_limit){
      return(FALSE);
   }

   if (tm->has_ub && (tm->par.gap_limit >= 0.0)){
      find_tree_lb(tm);
      if (100*(tm->ub-tm->lb)/(fabs(tm->ub)+etol) <= tm->par.gap_limit){
	 return(FALSE);
      }
   }

   rand_num = ((double)(RANDOM()))/((double)(MAXINT));
   if (tm->par.unconditional_dive_frac > 1 - rand_num){
      dive = CHECK_BEFORE_DIVE;
   }else{
      switch(tm->par.diving_strategy){
       case BEST_ESTIMATE:
	 if (tm->has_ub_estimate){
	    if (objval > tm->ub_estimate){
	       dive = DO_NOT_DIVE;
	       tm->stat.diving_halts++;
	    }else{
	       dive = CHECK_BEFORE_DIVE;
	    }
	    break;
	 }
       case COMP_BEST_K:
	 average_lb = 0;
#pragma omp critical (tree_update)
	 for (k = 0, i = MIN(tm->samephase_candnum, tm->par.diving_k);
	      i > 0; i--)
	    if (tm->samephase_cand[i]->lower_bound < MAXDOUBLE/2){
	       average_lb += tm->samephase_cand[i]->lower_bound;
	       k++;
	    }
	 if (k){
	    average_lb /= k;
	 }else{
	    dive = CHECK_BEFORE_DIVE;
	    break;
	 }
         if (fabs(average_lb) < etol) {
            average_lb = (average_lb > 0) ? etol : -etol;
            if (fabs(objval) < etol) {
               objval = (objval > 0) ? etol : -etol;
            }
         }
	 if (fabs((objval/average_lb)-1) > tm->par.diving_threshold){
	    dive = DO_NOT_DIVE;
	    tm->stat.diving_halts++;
	 }else{
	    dive = CHECK_BEFORE_DIVE;
	 }
	 break;
       case COMP_BEST_K_GAP:
	 average_lb = 0;
	 for (k = 0, i = MIN(tm->samephase_candnum, tm->par.diving_k);
	      i > 0; i--)
	    if (tm->samephase_cand[i]->lower_bound < MAXDOUBLE/2){
	       average_lb += tm->samephase_cand[i]->lower_bound;
	       k++;
	    }
	 if (k){
	    average_lb /= k;
	 }else{
	    dive = CHECK_BEFORE_DIVE;
	    break;
	 }
	 if (tm->has_ub)
	    cutoff = tm->par.diving_threshold*(tm->ub - average_lb);
	 else
	    cutoff = (1 + tm->par.diving_threshold)*average_lb;
	 if (objval > average_lb + cutoff){
	    dive = DO_NOT_DIVE;
	    tm->stat.diving_halts++;
	 }else{
	    dive = CHECK_BEFORE_DIVE;
	 }
	 break;
       default:
	 printf("Unknown diving strategy -- diving by default\n");
	 dive = DO_DIVE;
	 break;
      }
   }
   return(dive);
}

/*===========================================================================*/

/*===========================================================================*\
 * This routine is entirely for saving memory. If there is no need to
 * keep the description of the pruned nodes in memory, they are freed as
 * soon as they are no longer needed. This can set off a chain reaction
 * of other nodes that are no longer needed.
\*===========================================================================*/

int purge_pruned_nodes(tm_prob *tm, bc_node *node, int category)
{
   int i, new_child_num;
   branch_obj *bobj = &node->parent->bobj;
   char reason[30];
   char branch_dir = 'M';

   if (tm->par.vbc_emulation != VBC_EMULATION_FILE_NEW &&
	 (category == VBC_PRUNED_INFEASIBLE || category == VBC_PRUNED_FATHOMED
	  || category == VBC_IGNORE)) {
      printf("Error in purge_pruned_nodes.");
      printf("category refers to VBC_EMULATION_FILE_NEW");
      printf("when it is not used.\n");
      exit(456);
   }

   if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW) {
      switch (category) {
       case VBC_PRUNED_INFEASIBLE:
	 sprintf(reason,"%s","infeasible");
	 sprintf(reason,"%s %i %i",reason, node->bc_index+1,
	       node->parent->bc_index+1);
	 if (node->bc_index>0) {
	    if (node->parent->children[0]==node) {
	       branch_dir = node->parent->bobj.sense[0];
	       /*branch_dir = 'L';*/
	    } else {
	       branch_dir = node->parent->bobj.sense[1];
	       /*branch_dir = 'R';*/
	    }
	    if (branch_dir == 'G') {
	       branch_dir = 'R';
	    }
	 }
	 sprintf(reason,"%s %c %s", reason, branch_dir, "\n");
	 break;
       case VBC_PRUNED_FATHOMED:
	 sprintf(reason,"%s","fathomed");
	 sprintf(reason,"%s %i %i",reason, node->bc_index+1,
		 node->parent->bc_index+1);
	 if (node->bc_index>0) {
	    if (node->parent->children[0]==node) {
	       branch_dir = node->parent->bobj.sense[0];
	       /*branch_dir = 'L';*/
	    } else {
	       branch_dir = node->parent->bobj.sense[1];
	       /*branch_dir = 'R';*/
	    }
	    if (branch_dir == 'G') {
	       branch_dir = 'R';
	    }
	 }
	 sprintf(reason,"%s %c %s", reason, branch_dir, "\n");
	 break;
       case VBC_FEAS_SOL_FOUND:
	 /* This case has already been dealt in install_new_ub(), hence
	  * commented out
	  */
/*	 sprintf(reason,"%s","integer");
	 sprintf(reason,"%s %i %i",reason, node->bc_index+1,
		 node->parent->bc_index+1);
	 if (node->parent->children[0]==node) {
	    branch_dir = 'L';
	 } else {
	    branch_dir = 'R';
	 }
	 sprintf(reason,"%s %c %f\n", reason, branch_dir, tm->ub);
	 break;
 */
       default:
	 category = VBC_IGNORE;
	 break;
      }
   }

   if (node->parent == NULL){
      return(1);
   }
   if (category == VBC_IGNORE) {
#if 0
      PRINT(tm->par.verbosity, 1,
	    ("ignoring vbc update in purge_pruned_nodes"));
#endif
   } else if (tm->par.vbc_emulation == VBC_EMULATION_FILE){
      FILE *f;
#pragma omp critical(write_vbc_emulation_file)
      if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 PRINT_TIME(tm, f);
	 fprintf(f, "P %i %i\n", node->bc_index+1, category);
	 fclose(f);
      }
   } else if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$P %i %i\n", node->bc_index+1, category);
   } else if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW) {
      FILE *f;
#pragma omp critical(write_vbc_emulation_file)
      if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 PRINT_TIME2(tm, f);
	 fprintf(f, "%s", reason);
	 fclose(f);
      }
   }

   if ((new_child_num = --bobj->child_num) == 0){
      if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW) {
	 purge_pruned_nodes(tm, node->parent, VBC_IGNORE);
      } else {
	 purge_pruned_nodes(tm, node->parent, category);
      }
   }else{
      for (i = 0; i <= bobj->child_num; i++){
	 if (node->parent->children[i] == node){
	    if (i == new_child_num){
	       node->parent->children[i] = NULL;
	    }else{
	       node->parent->children[i]=node->parent->children[new_child_num];
	       bobj->sense[i] = bobj->sense[new_child_num];
	       bobj->rhs[i] = bobj->rhs[new_child_num];
	       bobj->range[i] = bobj->range[new_child_num];
	       bobj->branch[i] = bobj->branch[new_child_num];
	    }
	 }
      }
   }

   free_tree_node(node);
   return(1);
}

/*===========================================================================*\
 * This routine is for writing the pruned nodes to disk before deleting them
 * from memory.
\*===========================================================================*/

int write_pruned_nodes(tm_prob *tm, bc_node *node)
{
   FILE *f = NULL;
   branch_obj *bobj = &node->parent->bobj;

   if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_FULL ||
       tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
      if (!(f = fopen(tm->par.pruned_node_file_name, "a"))){
	 printf("\nError opening pruned node file\n\n");
	 return(0);
      }
   }

   if (node->parent == NULL){
      return(1);
   }

   if (bobj->child_num == 1){
      write_pruned_nodes(tm, node->parent);
   }

   if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_VBC_TOOL){
      if (node->parent)
	 fprintf(f, "%i %i\n", node->parent->bc_index + 1, node->bc_index + 1);
      fclose(f);
   }else if (tm->par.keep_description_of_pruned == KEEP_ON_DISK_FULL){
      write_node(node, tm->par.pruned_node_file_name, f, TRUE);
      fclose(f);
   }

   return(1);
}

/*===========================================================================*/

/*===========================================================================*\
 * Find the index of a particular process in a list
\*===========================================================================*/

int find_process_index(process_set *pset, int tid)
{
   int i = pset->procnum-1, *procs = pset->procs;

   for ( ; i >= 0 && procs[i] != tid; i--);
#ifdef DO_TESTS
   if (i == -1){
      printf("TM: process index not found !!!\n\n");
      exit(-5);
   }
#endif
   return(i);
}

/*===========================================================================*/

void mark_lp_process_free(tm_prob *tm, int lp, int cp)
{
   int ind;
   if (tm->cp.procnum > 0){
#ifdef COMPILE_IN_CP
      ind = cp;
#else
      ind = find_process_index(&tm->cp, cp);
#endif
      tm->active_nodes_per_cp[ind]--;
      if (tm->nodes_per_cp[ind] + tm->active_nodes_per_cp[ind] == 0)
	 tm->cp.free_ind[tm->cp.free_num++] = ind;
   }
   tm->active_nodes[lp] = NULL;
   tm->lp.free_ind[tm->lp.free_num++] = lp;
   tm->active_node_num--;
}

/*===========================================================================*/

int add_cut_to_list(tm_prob *tm, cut_data *cut)
{
#pragma omp critical (cut_pool)
   {
      REALLOC(tm->cuts, cut_data *, tm->allocated_cut_num, tm->cut_num + 1,
	      (tm->cut_num / tm->stat.created + 5) * BB_BUNCH);
      cut->name = tm->cut_num;
      tm->cuts[tm->cut_num++] = cut;
   }
   return(cut->name);
}

/*===========================================================================*/

/*===========================================================================*\
 * Installs a new upper bound and cleans up the candidate list
\*===========================================================================*/

void install_new_ub(tm_prob *tm, double new_ub, int opt_thread_num,
		    int bc_index, char branching, int feasible){
   bc_node *node, *temp, **list;
   int rule, pos, prev_pos, last, i;

   tm->has_ub = TRUE;
   tm->ub = new_ub;
#ifdef COMPILE_IN_LP
   tm->opt_thread_num = opt_thread_num;
#endif
   if (tm->par.vbc_emulation == VBC_EMULATION_FILE){
      FILE *f;
#pragma omp critical(write_vbc_emulation_file)
      if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 PRINT_TIME(tm, f);
	 fprintf(f, "U %.2f\n", new_ub);
	 fclose(f);
      }
   }else if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$U %.2f\n", new_ub);
   }else if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW &&
	     (feasible == IP_FEASIBLE || feasible == IP_HEUR_FEASIBLE)){
      FILE *f;
      //      char reason[30];
      char branch_dir = 'M';
      if (!(f = fopen(tm->par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else if ((feasible == IP_FEASIBLE && branching) ||
		(feasible == IP_HEUR_FEASIBLE)) {
#pragma omp critical(write_vbc_emulation_file)
	 PRINT_TIME2(tm, f);
	 fprintf(f, "%s %f %i\n", "heuristic", new_ub, bc_index+1);
      }else if (feasible == IP_FEASIBLE && !branching){
	 node = tm->active_nodes[opt_thread_num];
	 if (node->bc_index>0) {
	    if (node->parent->children[0]==node) {
	       branch_dir = node->parent->bobj.sense[0];
	       /*branch_dir = 'L';*/
	    } else {
	       branch_dir = node->parent->bobj.sense[1];
	       /*branch_dir = 'R';*/
	    }
	    if (branch_dir == 'G') {
	       branch_dir = 'R';
	    }
	 }
	 PRINT_TIME2(tm, f);
	 if (node->bc_index){
	    fprintf (f, "%s %i %i %c %f\n", "integer", node->bc_index+1,
		     node->parent->bc_index+1, branch_dir, new_ub);
	 }else{
	    fprintf (f, "%s %i %i %c %f\n", "integer", 1, 0, 'M', new_ub);
	 }

      }
      if (f){
	 fclose(f);
      }
   }

   /* Remove nodes that can now be fathomed from the list */
#pragma omp critical (tree_update)
   {
      rule = tm->par.node_selection_rule;
      list = tm->samephase_cand;
      char has_exchanged = FALSE;
      for (last = i = tm->samephase_candnum; i > 0; i--){
	 has_exchanged = FALSE;
	 node = list[i];
	 if (tm->has_ub &&
	     node->lower_bound >= tm->ub-tm->par.granularity){
	    if (i != last){
	       list[i] = list[last];
	       for (prev_pos = i, pos = i/2; pos >= 1;
		    prev_pos = pos, pos /= 2){
		  if (node_compar(rule, list[pos], list[prev_pos])){
		     temp = list[prev_pos];
		     list[prev_pos] = list[pos];
		     list[pos] = temp;
		     has_exchanged = TRUE;
		  }else{
		     break;
		  }
	       }
	    }
	    tm->samephase_cand[last] = NULL;
	    last--;
	    if (tm->par.verbosity > 0){
	       printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	       printf("+ TM: Pruning NODE %i LEVEL %i after new incumbent.\n",
		   node->bc_index, node->bc_level);
	       printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	    }
	    if (tm->par.keep_description_of_pruned == DISCARD ||
		tm->par.keep_description_of_pruned ==
		KEEP_ON_DISK_VBC_TOOL){
	       if (tm->par.keep_description_of_pruned ==
		   KEEP_ON_DISK_VBC_TOOL)
#pragma omp critical (write_pruned_node_file)
		  write_pruned_nodes(tm, node);
	       if (tm->par.vbc_emulation == VBC_EMULATION_FILE_NEW) {
		  purge_pruned_nodes(tm, node,
				     VBC_PRUNED_FATHOMED);
	       } else {
		  purge_pruned_nodes(tm, node, VBC_PRUNED);
	       }
	    }
	 }
	 if (has_exchanged) {
	    /*
	     * if exchanges have taken place, node[i] should be
	     * checked again for pruning
	     */
	    i++;
	 }
      }
      tm->samephase_candnum = last;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * This routine takes the description of the node after processing and
 * compares it to the description before processing. Then it updates
 * data structures appropriately. The other functions below are all
 * related to this one.
\*===========================================================================*/

void merge_descriptions(node_desc *old_node, node_desc *new_node)
{
   if (old_node->basis.basis_exists && new_node->basis.basis_exists){
      merge_base_stat(&old_node->basis.basevars, &new_node->basis.basevars);
      merge_extra_array_and_stat(&old_node->uind, &old_node->basis.extravars,
				 &new_node->uind, &new_node->basis.extravars);
      merge_base_stat(&old_node->basis.baserows, &new_node->basis.baserows);
      merge_extra_array_and_stat(&old_node->cutind,&old_node->basis.extrarows,
				 &new_node->cutind,&new_node->basis.extrarows);
   }else{
      old_node->basis = new_node->basis;
      merge_arrays(&old_node->uind, &new_node->uind);
      merge_arrays(&old_node->cutind, &new_node->cutind);
#ifdef COMPILE_IN_LP
      memset((char *)&(new_node->basis), 0, sizeof(basis_desc));
#endif
   }
   old_node->nf_status = new_node->nf_status;
   if (new_node->nf_status == NF_CHECK_AFTER_LAST ||
       new_node->nf_status == NF_CHECK_UNTIL_LAST){
      merge_arrays(&old_node->not_fixed, &new_node->not_fixed);
   }else{
      FREE(old_node->not_fixed.list);
   }
}

/*===========================================================================*/

void merge_base_stat(double_array_desc *dad, double_array_desc *moddad)
{
   if (moddad->type == EXPLICIT_LIST){
      FREE(dad->list);
      FREE(dad->stat);
      *dad = *moddad;
#ifdef COMPILE_IN_LP
      moddad->stat = NULL;
#endif
      return;
   }

   /* we've got a diff list against the old list */
   if (moddad->size > 0) { /* if there is a change */
      if (dad->type == EXPLICIT_LIST){
	 /* just overwrite the changes */
	 int i;
	 int *oldstat = dad->stat;
	 int  newsize = moddad->size;
	 int *newlist = moddad->list;
	 int *newstat = moddad->stat;
	 for (i = newsize - 1; i >= 0; i--)
	    oldstat[newlist[i]] = newstat[i];
      }else{
	 merge_double_array_descs(dad, moddad);
      }
   }
}

/*===========================================================================*/

void merge_extra_array_and_stat(array_desc *ad, double_array_desc *dad,
				array_desc *modad, double_array_desc *moddad)
{
   if (moddad->type != WRT_PARENT){
      /* If moddad is explicit then just use it */
      FREE(dad->list);
      FREE(dad->stat);
      *dad = *moddad;
#ifdef COMPILE_IN_LP
      moddad->stat = NULL;
#endif
   }else{
      /* So moddad is WRT. Then dad must be WRT, too!! */
      /* First throw out from dad everything that's just been deleted */
      int  newdeled = modad->size - modad->added;
      int *newtodel = modad->list + modad->added;
      int i, j, k, nextdel;
      if (newdeled > 0 && dad->size > 0){
	 int dsize = dad->size;
	 int *dlist = dad->list;
	 int *dstat = dad->stat;
	 k = j = 0;
	 for (i = 0; i < newdeled; ){
	    nextdel = newtodel[i];
	    for (; j < dsize && dlist[j] < nextdel; ){
	       dlist[k] = dlist[j];
	       dstat[k++] = dstat[j++];
	    }
	    if (j == dsize){
	       break;
	    }else{ /* in this case dlist[j] >= nextdel */
	       i++;
	       if (dlist[j] == nextdel)
		  j++;
	    }
	 }
	 while (j < dsize){
	    dlist[k] = dlist[j];
	    dstat[k++] = dstat[j++];
	 }
	 dad->size = k;
      }

      /* Now merge the remaining dad and moddad together */
      merge_double_array_descs(dad, moddad);
   }

   /* dad is fine now. Update ad */

   merge_arrays(ad, modad);
}

/*===========================================================================*/

/*===========================================================================*\
 * Merge the old and new changes together, in case of identical
 *  userindices the newer value overrides the older
\*===========================================================================*/

void merge_double_array_descs(double_array_desc *dad,
			      double_array_desc *moddad)
{
   if (moddad->size != 0){
      if (dad->size == 0){
	 *dad = *moddad;
#ifdef COMPILE_IN_LP
      moddad->stat = moddad->list = NULL;
#endif
      }else{
	 int i, j, k;
	 int  oldsize = dad->size;
	 int *oldlist = dad->list;
	 int *oldstat = dad->stat;
	 int newsize = moddad->size;
	 int *newlist = moddad->list;
	 int *newstat = moddad->stat;
	 int *dlist = dad->list = (int *) malloc((oldsize+newsize) * ISIZE);
	 int *dstat = dad->stat = (int *) malloc((oldsize+newsize) * ISIZE);

	 for (k = 0, i = j = 0; i < oldsize && j < newsize; ){
	    if (oldlist[i] < newlist[j]){
	       dlist[k] = oldlist[i];
	       dstat[k++] = oldstat[i++];
	    }else{
	       if (oldlist[i] == newlist[j])
		  i++;
	       dlist[k] = newlist[j];
	       dstat[k++] = newstat[j++];
	    }
	 }
	 while (i < oldsize){
	    dlist[k] = oldlist[i];
	    dstat[k++] = oldstat[i++];
	 }
	 while (j < newsize){
	    dlist[k] = newlist[j];
	    dstat[k++] = newstat[j++];
	 }
	 dad->size = k;
	 FREE(oldlist);
	 FREE(oldstat);
	 FREE(moddad->list);
	 FREE(moddad->stat);
      }
   }
}

/*===========================================================================*/

void merge_arrays(array_desc *array, array_desc *adesc)
{
   if (adesc->type != WRT_PARENT){
      /* Replace the old with the new one */
      FREE(array->list);
      *array = *adesc;
#ifdef COMPILE_IN_LP
      adesc->list = NULL;
#endif
      return;
   }

   if (adesc->size == 0){
      /* No change, leave the old description alone. */
      return;
   }

   if (array->size == 0){
      /* The old size is 0 (the new type is still WRT_PARENT).
	 If the old type was WRT_PARENT then we can simply replace it with
	 the new WRT_PARENT data.
	 If it was EXPLICIT_LIST with nothing in it, then... contradiction!
	 it would be shorter to explicitly list the new stuff then wrt an
	 empty explicit list. */
      *array = *adesc;
#ifdef COMPILE_IN_LP
      adesc->list = NULL;
#endif
   }else{
      /* OK, array is either WRT or EXP.
	 But!!! If array is EXP then array->added is set to array->size, so
	 we can handle it exactly as if it were WRT!
	 */
      /* Now comes the ugly part... we had an old WRT list and got a new
	 one wrt to the old (none of them is empty)... Create a correct one
	 now.
	 For extra vars/rows the basis status further complicates things.
	 If we had the basis stati stored as EXP, then we don't have to worry
	 about it, since we are going to receive an EXP. But if it was WRT
	 then we must delete the from the basis stati list the extras to be
	 deleted according to the new description! */
      int i, j, k, *list;

      int  newsize = adesc->size;
      int  newadded = adesc->added;
      int *newtoadd = adesc->list;
      int  newdeled = newsize - newadded;
      int *newtodel = newtoadd + newadded;
      int  oldadded = array->added;
      int *oldtoadd = array->list;
      int  olddeled = array->size - oldadded;
      int *oldtodel = oldtoadd + oldadded;

      /* cancel those both in oldtoadd and newdeled; also those both in
	 newtoadd and olddeled.
	 Then the unions of oldtoadd-newtodel and oldtodel-newtoadd will be
	 the correct list. */
      /* First count the number of collisions and mark the colliding ones */
      k = 0;
      for (i = j = 0; i < oldadded && j < newdeled; ){
	 if      (oldtoadd[i] < newtodel[j])  i++;
	 else if (oldtoadd[i] > newtodel[j])  j++;
	 else{
	    oldtoadd[i] = newtodel[j] = -1;
	    i++;
	    j++;
	    k++;
	 }
      }
      for (i = j = 0; i < newadded && j < olddeled; ){
	 if      (newtoadd[i] < oldtodel[j])  i++;
	 else if (newtoadd[i] > oldtodel[j])  j++;
	 else{
	    newtoadd[i] = oldtodel[j] = -1;
	    i++;
	    j++;
	    k++;
	 }
      }
      array->size = array->size + newsize - 2 * k;
      if (array->size == 0){
	 /* Nothing is left, great! */
	 array->added = 0;
	 FREE(adesc->list);
	 FREE(array->list);
      }else{
	 array->list = list = (int *) malloc(array->size * ISIZE);
	 for (i = j = k = 0; i < oldadded && j < newadded; ){
	    if      (oldtoadd[i] == -1)  i++;
	    else if (newtoadd[j] == -1)  j++;
	    else if (oldtoadd[i] < newtoadd[j])  list[k++] = oldtoadd[i++];
	    else /* can't be == */               list[k++] = newtoadd[j++];
	 }
	 for ( ; i < oldadded; i++)
	    if (oldtoadd[i] != -1) list[k++] = oldtoadd[i];
	 for ( ; j < newadded; j++)
	    if (newtoadd[j] != -1) list[k++] = newtoadd[j];
	 array->added = k;

	 for (i = j = 0; i < olddeled && j < newdeled; ){
	    if      (oldtodel[i] == -1)  i++;
	    else if (newtodel[j] == -1)  j++;
	    else if (oldtodel[i] < newtodel[j])  list[k++] = oldtodel[i++];
	    else /* can't be == */               list[k++] = newtodel[j++];
	 }
	 for ( ; i < olddeled; i++)
	    if (oldtodel[i] != -1) list[k++] = oldtodel[i];
	 for ( ; j < newdeled; j++)
	    if (newtodel[j] != -1) list[k++] = newtodel[j];

	 FREE(adesc->list); /* adesc->list */
	 FREE(oldtoadd); /* the old array->list */
      }
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * This routine modifies a list of integers, by deleting those listed in
 * todel and adding those listed in toadd. todel is a subset of iarray,
 * toadd has no common element with iarray.
\*===========================================================================*/

void modify_list(array_desc *origad, array_desc *modad)
{
   int j, k, l, nextdel;
   int added = modad->added;
   int *toadd = modad->list;
   int deled = modad->size - modad->added;
   int *todel = toadd + added;
   int size = origad->size;
   int *origlist = origad->list;

   if (deled){
      /* l is the location where we copy to; and k is where we copy from */
      l = k = 0;
      for (j = 0; j < deled; k++, j++){
	 nextdel = todel[j];
	 for (; origlist[k] != nextdel; origlist[l++] = origlist[k++]);
      }
      for (; k < size; origlist[l++] = origlist[k++]);
      size = l;
   }
   if (added){
      /* add toadd to origlist */
      for (l = size+added-1, k = added-1, j = size-1; k >= 0 && j >= 0; ){
	 if (origlist[j] > toadd[k])
	    origlist[l--] = origlist[j--];
	 else
	    origlist[l--] = toadd[k--];
      }
      if (k >= 0)
	 memcpy(origlist, toadd, (k+1) * ISIZE);
      size += added;
   }
   origad->size = size;
}

/*===========================================================================*/

void modify_list_and_stat(array_desc *origad, int *origstat,
			  array_desc *modad, double_array_desc *moddad)
{
   int i, j, k, l, nextdel;
   int added = modad->added;
   int *toadd = modad->list;
   int deled = modad->size - modad->added;
   int *todel = toadd + added;
   int size = origad->size;
   int *origlist = origad->list;

   /* First modify origad, and at the same time delete the appropriate entries
      from origstat as well as insert phony ones where needed */

   if (deled){
      /* l is the location where we copy to; and k is where we copy from */
      l = k = 0;
      for (j = 0; j < deled; k++, j++){
	 nextdel = todel[j];
	 for (; origlist[k] != nextdel; ){
	    origstat[l] = origstat[k];
	    origlist[l++] = origlist[k++];
	 }
      }
      for (; k < size; ){
	 origstat[l] = origstat[k];
	 origlist[l++] = origlist[k++];
      }
      size = l;
   }
   if (added){
      /* add toadd to origlist */
      for (l = size+added-1, k = added-1, j = size-1; k >= 0 && j >= 0; )
	 if (origlist[j] > toadd[k]){
	    origstat[l] = origstat[j];
	    origlist[l--] = origlist[j--];
	 }else{
	    origstat[l] = INVALID_BASIS_STATUS;
	    origlist[l--] = toadd[k--];
	 }
      if (k >= 0){
	 for ( ; k >= 0; ){
	    origstat[l] = INVALID_BASIS_STATUS;
	    origlist[l--] = toadd[k--];
	 }
      }
      size += added;
   }
   origad->size = size;

#ifdef DO_TM_BASIS_TESTS
   if (origad->size == 0 && moddad->size > 0){
      printf("TM: Problem with storing the basis!!\n\n");
      exit(990);
   }
#endif

   /* Now adjust the basis stati */
   if (origad->size > 0 && moddad->size > 0){
      int *modlist = moddad->list;
      int *modstat = moddad->stat;

      for (i = moddad->size - 1, j = origad->size - 1; i >= 0 && j >= 0; j--){
	 if (origlist[j] == modlist[i])
	    origstat[j] = modstat[i--];
#ifdef DO_TM_BASIS_TESTS
	 else if (origlist[j] < modlist[i]){
	    printf("TM: Problem with storing the basis!!\n\n");
	    exit(990);
	 }
#endif
      }
#ifdef DO_TM_BASIS_TESTS
      if (i >= 0){
	 printf("TM: Problem with storing the basis!!\n\n");
	 exit(990);
      }
#endif
   }
#ifdef DO_TM_BASIS_TESTS
   for (j = origad->size - 1; j >= 0; j--){
      if (origstat[j] == INVALID_BASIS_STATUS){
	 printf("TM: Problem with storing the basis!!\n\n");
	 exit(990);
      }
   }
#endif
}

/*===========================================================================*/

/*===========================================================================*\
 * Do some tasks before phase 2. Thes are:
 * - price out variables in the root if requested
 * - build up the samephase_cand binary tree
 * - inform everyone about it
\*===========================================================================*/

int tasks_before_phase_two(tm_prob *tm)
{
#if !defined(COMPILE_IN_TM)||(defined(COMPILE_IN_TM)&&!defined(COMPILE_IN_LP))
   int s_bufid;
#endif
#ifdef COMPILE_IN_LP
   int num_threads = 1;
#else
   int r_bufid, msgtag, bytes, sender, not_fixed_size;
#endif
   int i;
   bc_node *n;
   int termcode = 0;

#ifdef COMPILE_IN_LP
#ifdef _OPENMP
   num_threads =  omp_get_num_threads();
#endif
   for (i = 0; i < num_threads; i++){
      free_node_desc(&tm->lpp[i]->desc);
      tm->lpp[i]->phase = 1;
   }
#else
   /* Notify the LPs about the start of the second phase and get back the
      timing data for the first phase */
   s_bufid = init_send(DataInPlace);
   msend_msg(tm->lp.procs, tm->lp.procnum, LP__SECOND_PHASE_STARTS);
#endif

   if (tm->par.price_in_root && tm->has_ub){
      /* put together a message and send it out to an LP process. At this
	 point tm->rootnode->{lp,cg,cp,sp} should still be set to wherever
	 it was first processed */
#ifdef DO_TESTS
      if ((!tm->rootnode->lp) ||
	  (!tm->rootnode->cg && tm->par.use_cg) ||
	 (!tm->rootnode->cp && tm->cp.procnum > 0)){
	 printf("When trying to send root for repricing, the root doesn't\n");
	 printf("   have some process id correctly set!\n\n");
	 exit(-100);
      }
#endif

      send_active_node(tm, tm->rootnode, COLGEN_REPRICING, 0);
   }

   /* trim the tree */
   tm->stat.leaves_before_trimming = tm->nextphase_candnum;
   if (tm->par.trim_search_tree && tm->has_ub)
      tm->stat.tree_size -= trim_subtree(tm, tm->rootnode);

   /* while the LP is working, build up the samephase_cand binary tree */
   REALLOC(tm->samephase_cand, bc_node *,
	   tm->samephase_cand_size, tm->nextphase_candnum + 1, BB_BUNCH);
   for (i = 0; i < tm->nextphase_candnum; i++){
      if ((n = tm->nextphase_cand[i])){
	 if (n->bc_index >= 0){
	    insert_new_node(tm, n);
	 }else{
	    free_tree_node(n);
	 }
      }
   }
   tm->stat.leaves_after_trimming = tm->samephase_candnum;

   if ((termcode = receive_lp_timing(tm)) < 0)
      return(SOMETHING_DIED);

   if (tm->par.price_in_root && tm->has_ub){
      /* receive what the LP has to say, what is the new not_fixed list.
       * also, incorporate that list into the not_fixed field of everything */
#ifdef COMPILE_IN_LP
      switch(process_chain(tm->lpp[0])){

      case FUNCTION_TERMINATED_NORMALLY:
	 break;

      case ERROR__NO_BRANCHING_CANDIDATE:
	 return(TM_ERROR__NO_BRANCHING_CANDIDATE);

      case ERROR__ILLEGAL_RETURN_CODE:
	 return(TM_ERROR__ILLEGAL_RETURN_CODE);

      case ERROR__NUMERICAL_INSTABILITY:
	 return(TM_ERROR__NUMERICAL_INSTABILITY);

      case ERROR__USER:
	 return(TM_ERROR__USER);

      }
#else
      char go_on;
      int nsize, nf_status;

      do{
	 r_bufid = receive_msg(tm->rootnode->lp, ANYTHING);
	 bufinfo(r_bufid, &bytes, &msgtag, &sender);
	 switch (msgtag){
	  case LP__NODE_DESCRIPTION:
	    n = (bc_node *) calloc(1, sizeof(bc_node));
	    receive_node_desc(tm, n);
	    tm->stat.root_lb = n->lower_bound;
	    if (n->node_status == NODE_STATUS__PRUNED){
	       /* Field day! Proved optimality! */
	       free_subtree(tm->rootnode);
	       tm->rootnode = n;
	       tm->samephase_candnum = tm->nextphase_candnum = 0;
	       return (FUNCTION_TERMINATED_NORMALLY);
	    }
	    /* Otherwise in 'n' we have the new description of the root node.
	       We don't care about the cuts, just the not fixed variables.
	       We got to pay attention to changes in uind and not_fixed.
	       We won't change the uind in root but put every change into
	       not_fixed.
	       The new not_fixed list will comprise of the vars added to uind
	       and whatever is in the new not_fixed. */
	    if (n->desc.uind.size > 0){
	       array_desc *uind = &n->desc.uind;
	       array_desc *ruind = &tm->rootnode->desc.uind;
	       int usize = uind->size;
	       int rusize = ruind->size;
	       int *ulist = uind->list;
	       int *rulist = ruind->list;
	       int j, k;
	       /* Kick out from uind those in root's uind */
	       for (i = 0, j = 0, k = 0; i < usize && j < rusize; ){
		  if (ulist[i] < rulist[j]){
		     /* a new element in uind */
		     ulist[k++] = ulist[i++];
		  }else if (ulist[i] < rulist[j]){
		     /* something got kicked out of ruind */
		     j++;
		  }else{ /* ulist[i] == rulist[j] */
		     /* It just stayed there peacefully */
		     i++; j++;
		  }
	       }
	       if (i < usize){
		  /* The rest are new */
		  for ( ; i < usize; i++, k++)
		     ulist[k] = ulist[i];
	       }

	       if ((usize = k) > 0){
		  if ((nsize = n->desc.not_fixed.size) == 0){
		     /* All we got is from uind */
		     n->desc.not_fixed.size = usize;
		     n->desc.not_fixed.list = ulist;
		     uind->list = NULL;
		  }else{
		     /* Now merge whatever is left in ulist with not_fixed.
			Note that the two lists are disjoint. */
		     int *not_fixed = (int *) malloc((usize + nsize) * ISIZE);
		     int *nlist = n->desc.not_fixed.list;
		     for (not_fixed_size=i=j=k=0; i < usize && j < nsize;
			  not_fixed_size++){
			if (ulist[i] < nlist[j]){
			   not_fixed[k++] = ulist[i++];
			}else if (ulist[i] > nlist[j]){
			   not_fixed[k++] = nlist[j++];
			}else{
			   not_fixed[k++] = nlist[j++];
			   i++;
			}
		     }
		     if (i < usize)
			memcpy(not_fixed+k, ulist+i, (usize-i)*ISIZE);
		     if (j < nsize)
			memcpy(not_fixed+k, nlist+j, (nsize-j)*ISIZE);
		     FREE(nlist);
		     n->desc.not_fixed.size = not_fixed_size;
		     n->desc.not_fixed.list = not_fixed;
		  }
	       }
	    }

	    /* OK, now every new thingy is in n->desc.not_fixed */
	    nsize = n->desc.not_fixed.size;
	    if (nsize == 0){
	       /* Field day! Proved optimality!
		  Caveats:
		  This proves optimality, but the current tree may not contain
		  this proof, since the cuts used in pricing out might be
		  different from those originally in the root.
		  For now just accept this sad fact and report optimality.
		  Later, when the tree could be written out on disk, take care
		  of writing out BOTH root descriptions to prove optimality.
		  FIXME */
	       if (tm->par.keep_description_of_pruned){
		  /* We got to write it out here. */
	       }
	       free_tree_node(n);
	       tm->samephase_candnum = tm->nextphase_candnum = 0;
	       return(FUNCTION_TERMINATED_NORMALLY);
	    }else{
	       tm->rootnode->desc.not_fixed.list = n->desc.not_fixed.list;
	       n->desc.not_fixed.list = NULL;
	       if (nsize > tm->par.not_fixed_storage_size){
		  tm->rootnode->desc.not_fixed.size =
		     tm->par.not_fixed_storage_size;
		  nf_status = NF_CHECK_AFTER_LAST;
	       }else{
		  tm->rootnode->desc.not_fixed.size = nsize;
		  nf_status = NF_CHECK_UNTIL_LAST;
	       }
	    }
	    propagate_nf_status(tm->rootnode, nf_status);
	    tm->stat.nf_status = nf_status;
	    tm->stat.vars_not_priced = tm->rootnode->desc.not_fixed.size;
	    free_tree_node(n);
	    go_on = FALSE;
	    break;

	  case UPPER_BOUND:
	    process_ub_message(tm);
	    go_on = TRUE;
	    break;

	  case LP__CUT_NAMES_REQUESTED:
	    unpack_cut_set(tm, sender, 0, NULL);
	    go_on = TRUE;
	    break;

	  default: /* We shouldn't get anything else */
	    printf("Unexpected message at repricing! (%i)\n\n", msgtag);
	    return(ERROR__COMM_ERROR);
	 }
      }while (go_on);
#endif
   }

#ifdef COMPILE_IN_TM
   if (tm->samephase_candnum > 0){
      printf( "\n");
      printf( "**********************************************\n");
      printf( "* Branch and Cut First Phase Finished!!!!    *\n");
      printf( "* Now displaying stats and best solution...  *\n");
      printf( "**********************************************\n\n");

      print_statistics(&(tm->comp_times), &(tm->stat), &(tm->lp_stat),
                       tm->ub, tm->lb, 0,
		       tm->start_time, wall_clock(NULL),
		       tm->obj_offset,
		       tm->obj_sense, tm->has_ub,NULL);
   }
#else
   /* Report to the master all kind of statistics */
   s_bufid = init_send(DataInPlace);
   send_char_array((char *)&tm->comp_times, sizeof(node_times));
   send_char_array((char *)&tm->lp_stat, sizeof(lp_stat_desc));
   send_dbl_array(&tm->lb, 1);
   send_char_array((char *)&tm->stat, sizeof(tm_stat));
   send_msg(tm->master, TM_FIRST_PHASE_FINISHED);
#endif

   tm->nextphase_candnum = 0;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

/*===========================================================================*\
 * For now we'll trim the tree only between phases when nothing is processed.
 * Reason: It's kinda ugly to cut off a subtree when several of its nodes
 * might be processed. Nevertheless, this should be done sooner or later.
 *
 * FIXME!
 *
\*===========================================================================*/

int trim_subtree(tm_prob *tm, bc_node *n)
{
   int i, deleted = 0, not_pruned = 0;

   /* Theer isn't anything to do if this is a leaf. */
   if (n->bobj.child_num == 0)
      return(0);

   /* There isn't anything to do if all children are pruned, and we are
      better off to go down if only one is not pruned. */
   for (i = n->bobj.child_num - 1; i >= 0; i--)
      if (n->children[i]->node_status != NODE_STATUS__PRUNED)
	 if (++not_pruned > 1)
	    break;
   if (not_pruned == 0)
      return(0);
   if (not_pruned == 1){
      for (i = n->bobj.child_num - 1; i >= 0; i--)
	 if (n->children[i]->node_status != NODE_STATUS__PRUNED){
	    deleted = trim_subtree(tm, n->children[i]);
	    break;
	 }
      return(deleted);
   }

   /* So there are at least two not pruned. */
   for (i = n->bobj.child_num - 1; i >= 0; i--)
      if (n->children[i]->lower_bound + tm->par.granularity < tm->ub)
	 break;

   /* if all children have high objval */
   if (i < 0){
      /* put back this node into the nodes_per_... stuff */
      if (tm->par.max_cp_num > 0 && n->cp)
#ifdef COMPILE_IN_CP
	 tm->nodes_per_cp[n->cp]++;
#else
	 tm->nodes_per_cp[find_process_index(&tm->cp, n->cp)]++;
#endif
      /* also put the node on the nextphase list */
      REALLOC(tm->nextphase_cand, bc_node *,
	      tm->nextphase_cand_size, tm->nextphase_candnum+1, BB_BUNCH);
      tm->nextphase_cand[tm->nextphase_candnum++] = n;
      /* get rid of the children */
      for (i = n->bobj.child_num - 1; i >= 0; i--)
	 deleted += mark_subtree(tm, n->children[i]);
      /* free the children description */
      FREE(n->children);
      n->bobj.child_num = 0;
#ifndef MAX_CHILDREN_NUM
      FREE(n->bobj.sense);
      FREE(n->bobj.rhs);
      FREE(n->bobj.range);
      FREE(n->bobj.branch);
#endif
      FREE(n->bobj.solutions); //added by asm4
   }else{
      /* try to trim every child */
      for (i = n->bobj.child_num - 1; i >= 0; i--)
	 deleted += trim_subtree(tm, n->children[i]);
   }
   return(deleted);
}

/*===========================================================================*/

int mark_subtree(tm_prob *tm, bc_node *n)
{
   int i, deleted = 0;
   /* mark the node trimmed */
   if (n->bobj.child_num == 0){
      /* if this was a leaf and it is not pruned,
	 delete it from the nodes_per_... stuff */
      if (n->node_status != NODE_STATUS__PRUNED){
	 if (tm->par.max_cp_num > 0 && n->cp){
#ifdef COMPILE_IN_CP
	    i = n->cp;
#else
	    i = find_process_index(&tm->cp, n->cp);
#endif
	    tm->nodes_per_cp[i]--;
	    if (tm->nodes_per_cp[i] + tm->active_nodes_per_cp[i] == 0)
	       tm->cp.free_ind[tm->cp.free_num++] = i;
	 }
	 n->bc_index = -1;
      }else{
	 /* if it was pruned already the free it now */
	 free_tree_node(n);
      }
   }else{
      /* if non-leaf then process the children recursively and then free
         the node */
      for (i = n->bobj.child_num - 1; i >= 0; i--)
	 deleted += mark_subtree(tm, n->children[i]);
      free_tree_node(n);
   }
   return(++deleted);
}

/*===========================================================================*/

void propagate_nf_status(bc_node *n, int nf_status)
{
   int i;

   for (i = n->bobj.child_num - 1; i >= 0; i--)
      propagate_nf_status(n->children[i], nf_status);
   n->desc.nf_status = nf_status;
}


/*===========================================================================*/

/*===========================================================================*\
 * The functions below are for logging. They write out all necessary
 * data to a file so that a warmstart can be made.
\*===========================================================================*/

int write_node(bc_node *node, char *file, FILE* f, char append)
{
   int i;
   char close = FALSE;

   if (!f){
      if (!(f = fopen(file, append ? "a" : "w"))){
	 printf("\nError opening node file\n\n");
	 return(0);
      }
      close = TRUE;
   }

   if (append)
      fprintf(f, "\n");
   fprintf(f, "NODE INDEX:      %i\n", node->bc_index);
   fprintf(f, "NODE LEVEL:      %i\n", node->bc_level);
   fprintf(f, "LOWER BOUND:     %f\n", node->lower_bound);
   fprintf(f, "NODE STATUS:     %i\n", (int)node->node_status);
#ifdef TRACE_PATH
   fprintf(f, "OPTIMAL PATH:    %i\n", (int)node->optimal_path);
#endif
   if (node->parent)
      fprintf(f, "PARENT INDEX:    %i\n", node->parent->bc_index);
   else
      fprintf(f, "PARENT INDEX:    -1\n");
   fprintf(f, "CHILDREN:        %i %i %i\n", (int)node->bobj.type,
	   node->bobj.name, node->bobj.child_num);
   for (i = 0; i < node->bobj.child_num; i++)
      fprintf(f, "%i %c %f %f %i\n", node->children[i]->bc_index,
	      node->bobj.sense[i], node->bobj.rhs[i],
	      node->bobj.range[i], node->bobj.branch[i]);
   fprintf(f, "NODE DESCRIPTION: %i\n", node->desc.nf_status);
   fprintf(f, "USER INDICES:    %i %i %i\n", (int)node->desc.uind.type,
	   node->desc.uind.size, node->desc.uind.added);
   for (i = 0; i < node->desc.uind.size; i++)
      fprintf(f, "%i\n", node->desc.uind.list[i]);
   fprintf(f, "NOT FIXED:       %i %i %i\n", (int)node->desc.not_fixed.type,
	   node->desc.not_fixed.size, node->desc.not_fixed.added);
   for (i = 0; i < node->desc.not_fixed.size; i++)
      fprintf(f, "%i\n", node->desc.not_fixed.list[i]);
   fprintf(f, "CUT INDICES:     %i %i %i\n", (int)node->desc.cutind.type,
	   node->desc.cutind.size, node->desc.cutind.added);
   for (i = 0; i < node->desc.cutind.size; i++)
      fprintf(f, "%i\n", node->desc.cutind.list[i]);
   fprintf(f, "BASIS: %i\n", (int)node->desc.basis.basis_exists);
   fprintf(f, "BASE VARIABLES:  %i %i\n", (int)node->desc.basis.basevars.type,
	   node->desc.basis.basevars.size);
   if (node->desc.basis.basevars.type == WRT_PARENT)
      for (i = 0; i < node->desc.basis.basevars.size; i++)
	 fprintf(f, "%i %i\n", node->desc.basis.basevars.list[i],
		 node->desc.basis.basevars.stat[i]);
   else
      for (i = 0; i < node->desc.basis.basevars.size; i++)
	 fprintf(f, "%i\n", node->desc.basis.basevars.stat[i]);
   fprintf(f, "EXTRA VARIABLES: %i %i\n", (int)node->desc.basis.extravars.type,
	   node->desc.basis.extravars.size);
   if (node->desc.basis.extravars.type == WRT_PARENT)
      for (i = 0; i < node->desc.basis.extravars.size; i++)
	 fprintf(f, "%i %i\n", node->desc.basis.extravars.list[i],
		 node->desc.basis.extravars.stat[i]);
   else
      for (i = 0; i < node->desc.basis.extravars.size; i++)
	 fprintf(f, "%i\n", node->desc.basis.extravars.stat[i]);
   fprintf(f, "BASE ROWS:       %i %i\n", (int)node->desc.basis.baserows.type,
	   node->desc.basis.baserows.size);
   if (node->desc.basis.baserows.type == WRT_PARENT)
      for (i = 0; i < node->desc.basis.baserows.size; i++)
	 fprintf(f, "%i %i\n", node->desc.basis.baserows.list[i],
		 node->desc.basis.baserows.stat[i]);
   else
      for (i = 0; i < node->desc.basis.baserows.size; i++)
	 fprintf(f, "%i\n", node->desc.basis.baserows.stat[i]);
   fprintf(f, "EXTRA ROWS:      %i %i\n", (int)node->desc.basis.extrarows.type,
	   node->desc.basis.extrarows.size);
   if (node->desc.basis.extrarows.type == WRT_PARENT)
      for (i = 0; i < node->desc.basis.extrarows.size; i++)
	 fprintf(f, "%i %i\n", node->desc.basis.extrarows.list[i],
		 node->desc.basis.extrarows.stat[i]);
   else
      for (i = 0; i < node->desc.basis.extrarows.size; i++)
	 fprintf(f, "%i\n", node->desc.basis.extrarows.stat[i]);

   if (close)
      fclose(f);

   return(1);
}

/*===========================================================================*/

int read_node(tm_prob *tm, bc_node *node, FILE *f, int **children)
{
   int i, parent = 0, tmp = 0;
   char str1[10], str2[10];

   if (f){
      fscanf(f, "%s %s %i", str1, str2, &node->bc_index);
      fscanf(f, "%s %s %i", str1, str2, &node->bc_level);
      fscanf(f, "%s %s %lf", str1, str2, &node->lower_bound);
      fscanf(f, "%s %s %i", str1, str2, &tmp);
      node->node_status = (char)tmp;
#ifdef TRACE_PATH
      fscanf(f, "%s %s %i\n", str1, str2, &tmp);
      node->optimal_path = (char)tmp;
#endif
      fscanf(f, "%s %s %i", str1, str2, &parent);
      fscanf(f, "%s %i %i %i", str1, &tmp,
	     &node->bobj.name, &node->bobj.child_num);
      node->bobj.type = (char)tmp;
      if (node->bobj.child_num){
#ifndef MAX_CHILDREN_NUM
	 node->bobj.sense = malloc(node->bobj.child_num*sizeof(char));
	 node->bobj.rhs = (double *) malloc(node->bobj.child_num*DSIZE);
	 node->bobj.range = (double *) malloc(node->bobj.child_num*DSIZE);
	 node->bobj.branch = (int *) malloc(node->bobj.child_num*ISIZE);
#endif
	 *children = (int *) malloc(node->bobj.child_num*ISIZE);
	 for (i = 0; i < node->bobj.child_num; i++)
	    fscanf(f, "%i %c %lf %lf %i", *children+i, node->bobj.sense+i,
		   node->bobj.rhs+i, node->bobj.range+i, node->bobj.branch+i);
      }
      fscanf(f, "%s %s %i", str1, str2, &node->desc.nf_status);
      fscanf(f, "%s %s %i %i %i", str1, str2, &tmp, &node->desc.uind.size,
	     &node->desc.uind.added);
      node->desc.uind.type = (char)tmp;
      if (node->desc.uind.size){
	 node->desc.uind.list = (int *) malloc(node->desc.uind.size*ISIZE);
	 for (i = 0; i < node->desc.uind.size; i++)
	    fscanf(f, "%i", node->desc.uind.list+i);
      }
      fscanf(f, "%s %s %i %i %i", str1, str2, &tmp,
	     &node->desc.not_fixed.size, &node->desc.not_fixed.added);
      node->desc.not_fixed.type = (char)tmp;
      if (node->desc.not_fixed.size){
	 node->desc.not_fixed.list =
	    (int *) malloc(node->desc.not_fixed.size*ISIZE);
	 for (i = 0; i < node->desc.not_fixed.size; i++)
	    fscanf(f, "%i", node->desc.not_fixed.list+i);
      }
      fscanf(f, "%s %s %i %i %i", str1, str2, &tmp,
	     &node->desc.cutind.size, &node->desc.cutind.added);
      node->desc.cutind.type = (char)tmp;
      if (node->desc.cutind.size){
	 node->desc.cutind.list = (int *) malloc(node->desc.cutind.size*ISIZE);
	 for (i = 0; i < node->desc.cutind.size; i++)
	    fscanf(f, "%i", node->desc.cutind.list+i);
      }
      fscanf(f, "%s %i", str1, &tmp);
      node->desc.basis.basis_exists = (char)tmp;
      fscanf(f, "%s %s %i %i", str1, str2, &tmp,
	     &node->desc.basis.basevars.size);
      node->desc.basis.basevars.type = (char)tmp;
      if (node->desc.basis.basevars.size){
	 node->desc.basis.basevars.stat =
	    (int *) malloc(node->desc.basis.basevars.size*ISIZE);
	 if (node->desc.basis.basevars.type == WRT_PARENT){
	    node->desc.basis.basevars.list =
	       (int *) malloc(node->desc.basis.basevars.size*ISIZE);
	    for (i = 0; i < node->desc.basis.basevars.size; i++)
	       fscanf(f, "%i %i", node->desc.basis.basevars.list+i,
		      node->desc.basis.basevars.stat+i);
	 }else{
	    for (i = 0; i < node->desc.basis.basevars.size; i++)
	       fscanf(f, "%i", node->desc.basis.basevars.stat+i);
	 }
      }
      fscanf(f, "%s %s %i %i", str1, str2, &tmp,
	     &node->desc.basis.extravars.size);
      node->desc.basis.extravars.type = (char)tmp;
      if (node->desc.basis.extravars.size){
	 node->desc.basis.extravars.stat =
	    (int *) malloc(node->desc.basis.extravars.size*ISIZE);
	 if (node->desc.basis.extravars.type == WRT_PARENT){
	    node->desc.basis.extravars.list =
	       (int *) malloc(node->desc.basis.extravars.size*ISIZE);
	    for (i = 0; i < node->desc.basis.extravars.size; i++)
	       fscanf(f, "%i %i", node->desc.basis.extravars.list+i,
		      node->desc.basis.extravars.stat+i);
	 }else{
	    for (i = 0; i < node->desc.basis.extravars.size; i++)
	       fscanf(f, "%i", node->desc.basis.extravars.stat+i);
	 }
      }
      fscanf(f, "%s %s %i %i", str1, str2, &tmp,
	     &node->desc.basis.baserows.size);
      node->desc.basis.baserows.type = (char)tmp;
      if (node->desc.basis.baserows.size){
	 node->desc.basis.baserows.stat =
	    (int *) malloc(node->desc.basis.baserows.size*ISIZE);
	 if (node->desc.basis.baserows.type == WRT_PARENT){
	    node->desc.basis.baserows.list =
	       (int *) malloc(node->desc.basis.baserows.size*ISIZE);
	    for (i = 0; i < node->desc.basis.baserows.size; i++)
	       fscanf(f, "%i %i", node->desc.basis.baserows.list+i,
		      node->desc.basis.baserows.stat+i);
	 }else{
	    for (i = 0; i < node->desc.basis.baserows.size; i++)
	       fscanf(f, "%i", node->desc.basis.baserows.stat+i);
	 }
      }
      fscanf(f, "%s %s %i %i", str1, str2, &tmp,
	     &node->desc.basis.extrarows.size);
      node->desc.basis.extrarows.type = (char)tmp;
      if (node->desc.basis.extrarows.size){
	 node->desc.basis.extrarows.stat =
	    (int *) malloc(node->desc.basis.extrarows.size*ISIZE);
	 if (node->desc.basis.extrarows.type == WRT_PARENT){
	    node->desc.basis.extrarows.list =
	       (int *) malloc(node->desc.basis.extrarows.size*ISIZE);
	    for (i = 0; i < node->desc.basis.extrarows.size; i++)
	       fscanf(f, "%i %i", node->desc.basis.extrarows.list+i,
		      node->desc.basis.extrarows.stat+i);
	 }else{
	    for (i = 0; i < node->desc.basis.extrarows.size; i++)
	       fscanf(f, "%i", node->desc.basis.extrarows.stat+i);
	 }
      }
   }

   switch (node->node_status){
    case NODE_STATUS__HELD:
      REALLOC(tm->nextphase_cand, bc_node *,
	      tm->nextphase_cand_size, tm->nextphase_candnum+1, BB_BUNCH);
      tm->nextphase_cand[tm->nextphase_candnum++] = node;
      /* update the nodes_per_... stuff */
      /* the active_nodes_per_... will be updated when the LP__IS_FREE
	 message comes */
      if (node->cp)
#ifdef COMPILE_IN_CP
	 tm->nodes_per_cp[node->cp]++;
#else
	 tm->nodes_per_cp[find_process_index(&tm->cp, node->cp)]++;
#endif
      break;
    case NODE_STATUS__ROOT:
      tm->rootnode = node;
      break;
    case NODE_STATUS__WARM_STARTED:
    case NODE_STATUS__CANDIDATE:
#pragma omp critical (tree_update)
      insert_new_node(tm, node);
      break;
   }

   return(parent);
}

/*===========================================================================*/

int write_subtree(bc_node *root, char *file, FILE *f, char append, int logging)
{
   int i;
   char close = FALSE;

   if (!f){
      if (!(f = fopen(file, append ? "a" : "w"))){
	 printf("\nError opening subtree file\n\n");
	 return(0);
      }
      close = TRUE;
   }

   if (logging == VBC_TOOL){
      if (root->parent)
	 fprintf(f, "%i %i\n", root->parent->bc_index + 1, root->bc_index + 1);
   }else{
      write_node(root, file, f, append);
   }
   for (i = 0; i < root->bobj.child_num; i++)
      write_subtree(root->children[i], file, f, TRUE, logging);

   if (close)
      fclose(f);

   return(1);
}

/*===========================================================================*/

int read_subtree(tm_prob *tm, bc_node *root, FILE *f)
{
   int parent, i;
   int *children;

   parent = read_node(tm, root, f, &children);
   if (f && root->bobj.child_num){
      root->children = (bc_node **)
	 malloc(root->bobj.child_num*sizeof(bc_node *));
      for (i = 0; i < root->bobj.child_num; i++){
	 root->children[i] = (bc_node *) calloc(1, sizeof(bc_node));
	 root->children[i]->parent = root;
      }
   }
   for (i = 0; i < root->bobj.child_num; i++){
      read_subtree(tm, root->children[i], f);
   }

   return(parent);
}

/*===========================================================================*/

int write_tm_cut_list(tm_prob *tm, char *file, char append)
{
   FILE *f;
   int i, j;

   if (!(f = fopen(file, append ? "a" : "w"))){
      printf("\nError opening cut file\n\n");
      return(0);
   }

   fprintf(f, "CUTNUM: %i %i\n", tm->cut_num, tm->allocated_cut_num);
   for (i = 0; i < tm->cut_num; i++){
      fprintf(f, "%i %i %i %c %i %f %f\n", tm->cuts[i]->name,
	      tm->cuts[i]->size, (int)tm->cuts[i]->type, tm->cuts[i]->sense,
	      (int)tm->cuts[i]->branch, tm->cuts[i]->rhs, tm->cuts[i]->range);
      for (j = 0; j < tm->cuts[i]->size; j++)
	 fprintf(f, "%i ", (int)tm->cuts[i]->coef[j]);
      fprintf(f, "\n");
   }

   fclose(f);

   return(1);
}

/*===========================================================================*/

int read_tm_cut_list(tm_prob *tm, char *file)
{
   FILE *f;
   int i, j, tmp1 = 0, tmp2 = 0;
   char str[20];

   if (!(f = fopen(file, "r"))){
      printf("\nError opening cut file\n\n");
      return(0);
   }

   fscanf(f, "%s %i %i", str, &tm->cut_num, &tm->allocated_cut_num);
   tm->cuts = (cut_data **) malloc(tm->allocated_cut_num*sizeof(cut_data *));
   for (i = 0; i < tm->cut_num; i++){
      tm->cuts[i] = (cut_data *) malloc(sizeof(cut_data));
      fscanf(f, "%i %i %i %c %i %lf %lf", &tm->cuts[i]->name,
	     &tm->cuts[i]->size, &tmp1, &tm->cuts[i]->sense,
	     &tmp2, &tm->cuts[i]->rhs, &tm->cuts[i]->range);
      tm->cuts[i]->type = (char)tmp1;
      tm->cuts[i]->branch = (char)tmp2;
      tm->cuts[i]->coef = (char *) malloc(tm->cuts[i]->size*sizeof(char));
      for (j = 0; j < tm->cuts[i]->size; j++){
	 fscanf(f, "%i ", &tmp1);
	 tm->cuts[i]->coef[j] = (char)tmp1;
      }
   }

   fclose(f);

   return(1);
}

/*===========================================================================*/

int write_tm_info(tm_prob *tm, char *file, FILE* f, char append)
{
   char close = FALSE;

   if (!f){
      if (!(f = fopen(file, append ? "a" : "w"))){
	 printf("\nError opening TM info file\n\n");
	 return(0);
      }
      close = TRUE;
   }

   if (tm->par.logging == VBC_TOOL){
      fprintf(f, "#TYPE: COMPLETE TREE\n");
      fprintf(f, "#TIME: NOT\n");
      fprintf(f, "#BOUNDS: NONE\n");
      fprintf(f, "#INFORMATION: EXCEPTION\n");
      fprintf(f, "#NODE_NUMBER: NONE\n");
      if (close)
	 fclose(f);

      return(1);
   }

   fprintf(f, "UPPER BOUND: ");
   if (tm->has_ub)
      fprintf(f, "   %f\n", tm->ub);
   else
      fprintf(f, "\n");
   fprintf(f, "LOWER BOUND:    %f\n", tm->lb);
   fprintf(f, "PHASE:          %i\n", tm->phase);
   fprintf(f, "ROOT LB:        %f\n", tm->stat.root_lb);
   fprintf(f, "MAX DEPTH:      %i\n", tm->stat.max_depth);
   fprintf(f, "CHAINS:         %i\n", tm->stat.chains);
   fprintf(f, "DIVING HALTS:   %i\n", tm->stat.diving_halts);
   fprintf(f, "TREE SIZE:      %i\n", tm->stat.tree_size);
   fprintf(f, "NODES CREATED:  %i\n", tm->stat.created);
   fprintf(f, "NODES ANALYZED: %i\n", tm->stat.analyzed);
   fprintf(f, "LEAVES BEFORE:  %i\n", tm->stat.leaves_before_trimming);
   fprintf(f, "LEAVES AFTER:   %i\n", tm->stat.leaves_after_trimming);
   fprintf(f, "NF STATUS:      %i\n", (int)tm->stat.nf_status);
   fprintf(f, "TIMING:\n");
   fprintf(f, " COMM:          %f\n", tm->comp_times.communication);
   fprintf(f, " LP:            %f\n", tm->comp_times.lp);
   fprintf(f, " SEPARATION:    %f\n", tm->comp_times.separation);
   fprintf(f, " FIXING:        %f\n", tm->comp_times.fixing);
   fprintf(f, " PRICING:       %f\n", tm->comp_times.pricing);
   fprintf(f, " BRANCHING:     %f\n", tm->comp_times.strong_branching);
   fprintf(f, " CUT POOL:      %f\n", tm->comp_times.cut_pool);
   fprintf(f, " REAL TIME:     %f\n", wall_clock(NULL) - tm->start_time);

   if (close)
      fclose(f);

   return(1);
}

/*===========================================================================*/

int read_tm_info(tm_prob *tm, FILE *f)
{
   char str1[20], str2[20];
   int tmp = 0;
   double previous_elapsed_time = 0;

   if (!f)
      return(0);

   fscanf(f, "%s %s", str1, str2);
   if (fscanf(f, "%lf", &tm->ub) != 0)
      tm->has_ub = TRUE;
   fscanf(f, "%s %s %lf", str1, str2, &tm->lb);
   fscanf(f, "%s %i", str1, &tm->phase);
   fscanf(f, "%s %s %lf", str1, str2, &tm->stat.root_lb);
   fscanf(f, "%s %s %i", str1, str2, &tm->stat.max_depth);
   fscanf(f, "%s %i", str1, &tm->stat.chains);
   fscanf(f, "%s %s %i", str1, str2, &tm->stat.diving_halts);
   fscanf(f, "%s %s %i", str1, str2, &tm->stat.tree_size);
   fscanf(f, "%s %s %i", str1, str2, &tm->stat.created);
   fscanf(f, "%s %s %i", str1, str2, &tm->stat.analyzed);
   fscanf(f, "%s %s %i", str1, str2, &tm->stat.leaves_before_trimming);
   fscanf(f, "%s %s %i", str1, str2, &tm->stat.leaves_after_trimming);
   fscanf(f, "%s %s %i", str1, str2, &tmp);
   tm->stat.nf_status = (char)tmp;
   fscanf(f, "%s", str1);
   fscanf(f, "%s %lf", str1, &tm->comp_times.communication);
   fscanf(f, "%s %lf", str1, &tm->comp_times.lp);
   fscanf(f, "%s %lf", str1, &tm->comp_times.separation);
   fscanf(f, "%s %lf", str1, &tm->comp_times.fixing);
   fscanf(f, "%s %lf", str1, &tm->comp_times.pricing);
   fscanf(f, "%s %lf", str1, &tm->comp_times.strong_branching);
   fscanf(f, "%s %s %lf", str1, str2, &tm->comp_times.cut_pool);
   fscanf(f, "%s %s %lf\n", str1, str2, &previous_elapsed_time);
   tm->start_time -= previous_elapsed_time;

   return(1);
}

/*===========================================================================*/

int write_base(base_desc *base, char *file, FILE *f, char append)
{
   int i;
   char close = FALSE;

   if (!f){
      if (!(f = fopen(file, append ? "a" : "w"))){
	 printf("\nError opening base file\n\n");
	 return(0);
      }
      close = TRUE;
   }

   fprintf(f, "BASE DESCRIPTION: %i %i\n", base->varnum, base->cutnum);
   for (i = 0; i < base->varnum; i++)
      fprintf(f, "%i\n", base->userind[i]);

   if (close)
      fclose(f);

   return(1);
}

/*===========================================================================*/

int read_base(base_desc *base, FILE *f)
{
   char str1[20], str2[20];
   int i;

   fscanf(f, "%s %s %i %i", str1, str2, &base->varnum, &base->cutnum);
   base->userind = (int *) malloc(base->varnum*ISIZE);
   for (i = 0; i < base->varnum; i++)
      fscanf(f, "%i", base->userind+i);

   return(1);
}

/*===========================================================================*/

/*===========================================================================*\
 * Cleanup. Free the datastructure.
\*===========================================================================*/

void free_tm(tm_prob *tm)
{
   int i;
   cut_data **cuts = tm->cuts;
#ifdef _OPENMP
   int num_threads = tm->par.max_active_nodes;
#else
   int num_threads = 1;
#endif

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   for (i = 0; i < num_threads; i++)
      free_lp(tm->lpp[i]);
   FREE(tm->lpp);
#ifdef COMPILE_IN_CG
   FREE(tm->cgp);
#endif
#endif

   if (tm->par.lp_machs){
      FREE(tm->par.lp_machs[0]);
      FREE(tm->par.lp_machs);
   }
   if (tm->par.cg_machs){
      FREE(tm->par.cg_machs[0]);
      FREE(tm->par.cg_machs);
   }
   if (tm->par.cp_machs){
      FREE(tm->par.cp_machs[0]);
      FREE(tm->par.cp_machs);
   }
   FREE(tm->lp.procs);
   FREE(tm->lp.free_ind);
   FREE(tm->cg.procs);
   FREE(tm->cg.free_ind);
   FREE(tm->cp.procs);
   FREE(tm->cp.free_ind);
   FREE(tm->nodes_per_cp);
   FREE(tm->active_nodes_per_cp);

   FREE(tm->active_nodes);
   FREE(tm->samephase_cand);
   FREE(tm->nextphase_cand);

#ifndef COMPILE_IN_TM
   /* Go over the tree and free the nodes */
   free_subtree(tm->rootnode);
#endif

   /* Go over the cuts stored and free them all */
#pragma omp critical (cut_pool)
   if (cuts){
      for (i = tm->cut_num - 1; i >= 0; i--)
	 if (cuts[i]){
	    FREE(cuts[i]->coef);
	    FREE(cuts[i]);
	 }
      FREE(tm->cuts);
   }

#pragma omp critical (tmp_memory)
   {
      FREE(tm->tmp.i);
      FREE(tm->tmp.c);
      FREE(tm->tmp.d);
   }

   /*get rid of the added pointers for sens.analysis*/

   for (i = 0; i < num_threads; i++){
      if(tm->rpath[i])
	 if(tm->rpath[i][0])
	    tm->rpath[i][0] = NULL;
      FREE(tm->bpath[i]);
      FREE(tm->rpath[i]);
   }

   FREE(tm->rpath);
   FREE(tm->rpath_size);
   FREE(tm->bpath);
   FREE(tm->bpath_size);

   if (tm->reduced_costs) {
      for (i=0; i<tm->reduced_costs->num_rcs; i++) {
         FREE(tm->reduced_costs->indices[i]);
         FREE(tm->reduced_costs->values[i]);
         FREE(tm->reduced_costs->lb[i]);
         FREE(tm->reduced_costs->ub[i]);
      }
      FREE(tm->reduced_costs->indices);
      FREE(tm->reduced_costs->values);
      FREE(tm->reduced_costs->lb);
      FREE(tm->reduced_costs->ub);
      FREE(tm->reduced_costs->cnt);
      FREE(tm->reduced_costs->obj);
      FREE(tm->reduced_costs);
   }

   if (tm->pcost_down) {
      FREE(tm->pcost_down);
      FREE(tm->pcost_up);
      FREE(tm->br_rel_down);
      FREE(tm->br_rel_up);
      FREE(tm->br_rel_cand_list);
      FREE(tm->br_rel_down_min_level);
      FREE(tm->br_rel_up_min_level);
   }

   FREE(tm);
}

/*===========================================================================*/

void free_subtree(bc_node *n)
{
   int i;

   if (n == NULL) return;

   for (i = n->bobj.child_num - 1; i >= 0; i--)
      free_subtree(n->children[i]);
   free_tree_node(n);
}

/*===========================================================================*/

void free_tree_node(bc_node *n)
{

   FREE(n->sol);
   FREE(n->sol_ind);
#ifdef SENSITIVITY_ANALYSIS
   FREE(n->duals);
#endif
   FREE(n->children);
#ifndef MAX_CHILDREN_NUM
   FREE(n->bobj.sense);
   FREE(n->bobj.rhs);
   FREE(n->bobj.range);
   FREE(n->bobj.branch);
#endif
   FREE(n->bobj.solutions); //added by asm4

   FREE(n->desc.uind.list);
   free_basis(&n->desc.basis);
   FREE(n->desc.not_fixed.list);
   FREE(n->desc.cutind.list);
   FREE(n->desc.desc);
   if (n->desc.bnd_change) {
      FREE(n->desc.bnd_change->index);
      FREE(n->desc.bnd_change->lbub);
      FREE(n->desc.bnd_change->value);
      FREE(n->desc.bnd_change);
   }
   FREE(n);
}

/*===========================================================================*/

void free_basis(basis_desc *basis)
{
   FREE(basis->basevars.list);
   FREE(basis->basevars.stat);
   FREE(basis->extravars.list);
   FREE(basis->extravars.stat);
   FREE(basis->baserows.list);
   FREE(basis->baserows.stat);
   FREE(basis->extrarows.list);
   FREE(basis->extrarows.stat);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function shuts down the treemanager
\*===========================================================================*/

int tm_close(tm_prob *tm, int termcode)
{
#ifndef COMPILE_IN_TM
   int s_bufid;
#endif
#ifdef COMPILE_IN_LP
   lp_prob **lp = tm->lpp;
#endif
#ifndef COMPILE_IN_CP
   int r_bufid = 0, new_cuts;
   struct timeval timeout = {5, 0};
   double new_time;
#endif
   int i;

#if defined(DO_TESTS) && 0
   if (tm->cp.free_num != tm->cp.procnum)
      printf(" Something is fishy! tm->cp.freenum != tm->cp.procnum\n");
#endif

   if (tm->par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$#END_OF_OUTPUT");
   }

   /*------------------------------------------------------------------------*\
    * Kill the processes. Some of them will send back statistics.
   \*------------------------------------------------------------------------*/
#ifndef COMPILE_IN_LP
   stop_processes(&tm->lp);
#endif
#ifndef COMPILE_IN_CG
   stop_processes(&tm->cg);
#endif
#ifndef COMPILE_IN_CP
   stop_processes(&tm->cp);
#endif

   /*------------------------------------------------------------------------*\
    * Receive statistics from the cutpools
   \*------------------------------------------------------------------------*/
#ifdef COMPILE_IN_CP
   if (tm->cpp){
      for (i = 0; i < tm->par.max_cp_num; i++){
	 tm->comp_times.cut_pool += tm->cpp[i]->cut_pool_time;
	 tm->stat.cuts_in_pool += tm->cpp[i]->cut_num;
	 tm->cpp[i]->msgtag = YOU_CAN_DIE;
	 cp_close(tm->cpp[i]);
      }
      FREE(tm->cpp);
   }
#else
   for (i = 0; i < tm->par.max_cp_num;){
      r_bufid = treceive_msg(tm->cp.procs[i], POOL_TIME, &timeout);
      if (r_bufid > 0){
	 receive_dbl_array(&new_time, 1);
	 receive_int_array(&new_cuts, 1);
	 tm->comp_times.cut_pool += new_time;
	 tm->stat.cuts_in_pool += new_cuts;
	 i++;
      }else{
	 if (pstat(tm->cp.procs[i]) != PROCESS_OK)
	    i++;
      }
   }
   freebuf(r_bufid);
#endif
   /* Receive timing from the LPs */

   if (receive_lp_timing(tm) < 0){
      printf("\nWarning: problem receiving LP timing. LP process is dead\n\n");
   }

#ifdef COMPILE_IN_LP
   for (i = 0; i < tm->par.max_active_nodes; i ++){
      lp_close(lp[i]);
   }
#endif

   tm->stat.root_lb = tm->rootnode->lower_bound;
   find_tree_lb(tm);
   return(termcode);

#ifndef COMPILE_IN_TM
   /*------------------------------------------------------------------------*\
    * Send back the statistics to the master
   \*------------------------------------------------------------------------*/

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)&tm->comp_times, sizeof(node_times));
   send_char_array((char *)&tm->lp_stat, sizeof(lp_stat_desc));
   send_dbl_array(&tm->lb, 1);
   send_char_array((char *)&tm->stat, sizeof(tm_stat));
   send_msg(tm->master, termcode);
   freebuf(s_bufid);

   free_tm(tm);

#endif

   return(termcode);
}

/*===========================================================================*/
/*===========================================================================*/
#if !defined(_MSC_VER) && !defined(__MNO_CYGWIN) && defined(SIGHANDLER)
void sym_catch_c(int num)
{

   sigset_t mask_set;
   sigset_t old_set;

   /* SIGTSTP ? */
   signal(SIGINT, sym_catch_c);

   char temp [MAX_LINE_LENGTH + 1];

   sigfillset(&mask_set);
   sigprocmask(SIG_SETMASK, &mask_set, &old_set);

   strcpy(temp, "");
   printf("\nDo you want to abort immediately, exit gracefully (from the current solve call only), or continue? [a/e/c]: ");
   fflush(stdout);
   fgets(temp, MAX_LINE_LENGTH, stdin);
   if(temp[1] == '\n' && (temp[0] == 'a' || temp[0] == 'A')){
      printf("\nTerminating...\n");
      fflush(stdout);
      exit(0);
   }else if(temp[1] == '\n' && (temp[0] == 'e' || temp[0] == 'E')){
      c_count++;
   } else{
      printf("\nContinuing...\n");
      fflush(stdout);
      c_count = 0;
   }

}
#endif
/*===========================================================================*/
/*
 * Find the lowerbound of the current branch-and-cut tree and save it in
 * tm->lb
 */
int find_tree_lb(tm_prob *tm)
{
   double lb = MAXDOUBLE;
   bc_node **samephase_cand;

   if (tm->samephase_candnum > 0 || tm->active_node_num > 0) {
      if (tm->par.node_selection_rule == LOWEST_LP_FIRST) {
         lb = tm->samephase_cand[1]->lower_bound; /* [0] is a dummy */
      } else {
         samephase_cand = tm->samephase_cand;
         for (int i = tm->samephase_candnum; i >= 1; i--){
            if (samephase_cand[i]->lower_bound < lb) {
               lb = samephase_cand[i]->lower_bound;
            }
         }
      }
   } else {
      /* there are no more nodes left. */
      lb = tm->ub;
   }
   /*
   if (lb >= MAXDOUBLE / 2){
      lb = tm->ub;
   }
   */
   tm->lb = lb;
   return 0;
}
/*===========================================================================*/
