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

#include <stdlib.h>          /* malloc() is defined here in AIX ... */
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include "omp.h"
#endif

#include "sym_qsort.h"
#include "sym_messages.h"
#include "sym_proccomm.h"
#include "symphony.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_master.h"
#include "sym_master_u.h"
#include "sym_lp_solver.h"
#ifdef COMPILE_IN_TM
#include "sym_lp.h"
#endif

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the wrapper functions for the master process.
\*===========================================================================*/

int initialize_u(sym_environment *env)
{
   sym_set_defaults(env);

#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_initialize(&env->user) );
#else
   env->user = NULL;
#endif

   env->mip = (MIPdesc *) calloc(1, sizeof(MIPdesc));

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int readparams_u(sym_environment *env, int argc, char **argv)
{
   int i, user_res;
   char tmp, c, foundF, foundD;

   parse_command_line(env, argc, argv);

#ifdef USE_SYM_APPLICATION
   user_res = user_readparams(env->user, env->par.param_file, argc, argv);
#else
   user_res = USER_DEFAULT;
#endif

   switch(user_res){

    case USER_DEFAULT:

      foundF = foundD = FALSE;
      for (i = 1; i < argc; i++){
	 sscanf(argv[i], "%c %c", &tmp, &c);
	 if (tmp != '-')
	    continue;
	 switch (c) {
	  case 'L':
	    env->par.file_type = LP_FORMAT;
	  case 'F':
	    if (i < argc - 1){
	      sscanf(argv[i+1], "%c", &tmp);
	      if (tmp == '-'){
		printf("Warning: Missing argument to command-line switch -%c\n",
		       c);
	      }else{
		strncpy(env->par.infile, argv[++i],MAX_FILE_NAME_LENGTH);
		foundF = TRUE;
	      }
	    }else{
	      printf("Warning: Missing argument to command-line switch -%c\n",c);
	    }
	    break;
	 case 'D':
	   if (i < argc - 1){
	     sscanf(argv[i+1], "%c", &tmp);
	     if (tmp == '-'){
	       printf("Warning: Missing argument to command-line switch -%c\n",
		      c);
	     }else{
	       strncpy(env->par.datafile, argv[++i],MAX_FILE_NAME_LENGTH);
	       foundD = TRUE;
	     }
	   }else{
	     printf("Warning: Missing argument to command-line switch -%c\n",c);
	   }
	   break;
	 case 'T':
	   env->par.test = TRUE;
	   if(i+1 < argc){
	     sscanf(argv[i+1], "%c", &tmp);
	     if(tmp != '-'){
	       strncpy(env->par.test_dir, argv[++i],MAX_FILE_NAME_LENGTH);
	     }
	   }else{
	     printf("Warning: Missing argument to command-line switch -%c\n",c);
	   }
	   break;
	 default:
	   break;
	 }
	 if (foundF && foundD){
	   break;
	 }
      }
    case USER_SUCCESS:
    case USER_NO_PP:
    case USER_AND_PP:
      break;

    case USER_ERROR:

      printf("\n\n*********User error detected -- aborting***********\n\n");
      return(ERROR__USER);

    default:
      break;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int io_u(sym_environment *env)
{
   int err, user_res;

#ifdef USE_SYM_APPLICATION
   user_res = user_io(env->user);
#else
   user_res = USER_DEFAULT;
#endif

   switch( user_res ){

    case USER_DEFAULT:

      if (strcmp(env->par.infile, "") == 0){
	 printf("\nNo input file specified\n");
	 return (ERROR__READING_MPS_FILE);
      }

      if (env->par.verbosity >= 0){
	 printf("Reading input file...\n\n");
      }
      if (strcmp(env->par.datafile, "") == 0){
	 if (env->par.file_type == LP_FORMAT){
	    err = read_lp(env->mip, env->par.infile, env->probname);
	    env->par.file_type = MPS_FORMAT;
	    if (err != 0){
	       printf("\nErrors in reading LP file\n");
	       return (ERROR__READING_LP_FILE);
	    }
	 }else {
	    err = read_mps(env->mip, env->par.infile, env->probname);
	    if (err != 0){
	       printf("\nErrors in reading mps file\n");
	       return (ERROR__READING_MPS_FILE);
	    }
	 }
      }else{
#ifdef USE_GLPMPL
	 err = read_gmpl(env->mip, env->par.infile,
			 env->par.datafile, env->probname);
	 env->par.file_type = MPS_FORMAT;
	 if(!err){
	    printf("\nErrors in reading gmpl file\n");
	    return (ERROR__READING_GMPL_FILE);
	 }
#else
	 printf("ERROR: SYMPHONY can only read GMPL/AMPL files if GLPK is \n");
	 printf("installed and the USE_GLPMPL compiler define is set. \n");
	 printf("Exiting.\n\n");
	 return (ERROR__READING_GMPL_FILE);
#endif
      }

      break;

    case USER_ERROR:

      printf("\n\n*********User error detected -- aborting***********\n\n");
      return(ERROR__USER);
      break;

    case USER_SUCCESS:
    case USER_NO_PP:
    case USER_AND_PP:
    default:

      break;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int init_draw_graph_u(sym_environment *env)
{
   if (env->par.do_draw_graph){ /*start up the graphics window*/
      int s_bufid;
      if (env->par.dg_machine_set){
	 spawn(env->par.dg_exe, (char **)NULL, env->par.dg_debug | TaskHost,
	       env->par.dg_machine, 1, &env->dg_tid);
      }else{
	 spawn(env->par.dg_exe, (char **)NULL, env->par.dg_debug, (char *)NULL, 1,
	       &env->dg_tid);
      }
      s_bufid = init_send(DataInPlace);
      send_char_array((char *)&env->par.dg_par, sizeof(dg_params));
      send_msg(env->dg_tid, DG_DATA);
      freebuf(s_bufid);

#ifdef USE_SYM_APPLICATION
      if (env->dg_tid)
	 CALL_USER_FUNCTION( user_init_draw_graph(env->user, env->dg_tid) );
#endif
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int start_heurs_u(sym_environment *env)
{
   double ub = env->has_ub ? env->ub : -MAXDOUBLE;
   double ub_estimate = env->has_ub_estimate ? env->ub_estimate : -MAXDOUBLE;

#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_start_heurs(env->user, &ub, &ub_estimate) );
#endif

   if (!env->has_ub){
      if (ub > -MAXDOUBLE){
	 env->has_ub = TRUE;
	 env->ub = ub;
      }else{
	 env->ub = MAXDOUBLE;
      }
   }else if (ub < env->ub){
      env->ub = ub;
   }
   if (!env->has_ub_estimate){
      if (ub_estimate > -MAXDOUBLE){
	 env->has_ub_estimate = TRUE;
	 env->ub_estimate = ub_estimate;
      }
   }else if (ub_estimate < env->ub_estimate){
      env->ub_estimate = ub_estimate;
   }
   if (env->par.tm_par.vbc_emulation == VBC_EMULATION_FILE){
      FILE *f = NULL;
      if (!(f = fopen(env->par.tm_par.vbc_emulation_file_name, "a"))){
	 printf("\nError opening vbc emulation file\n\n");
      }else{
	 fprintf(f, "00:00:00.00 U %.2f \n", env->ub);
	 fclose(f);
      }
   }else if (env->par.tm_par.vbc_emulation == VBC_EMULATION_LIVE){
      printf("$U %.2f\n", env->ub);
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int initialize_root_node_u(sym_environment *env)
{
   int i, user_res;

   base_desc *base = env->base = (base_desc *) calloc(1, sizeof(base_desc));
   node_desc *root = env->rootdesc = (node_desc *) calloc(1,sizeof(node_desc));


#ifdef USE_SYM_APPLICATION
   user_res = user_initialize_root_node(env->user, &base->varnum,
					&base->userind,
					&base->cutnum, &root->uind.size,
					&root->uind.list, &env->mip->obj_sense,
					&env->mip->obj_offset,
					&env->mip->colname,
					env->par.tm_par.colgen_strat);
#else
   user_res = USER_DEFAULT;
#endif

   switch (user_res){

    case USER_ERROR:

      printf("\n\n*********User error detected -- aborting***********\n\n");
      return(ERROR__USER);

    case USER_SUCCESS:
    case USER_NO_PP:
    case USER_AND_PP:
      if (base->varnum)
	 qsort_i(base->userind, base->varnum);
      if (root->uind.size && !env->par.warm_start)
	 qsort_i(root->uind.list, root->uind.size);
      env->mip->n = base->varnum + root->uind.size;
      env->mip->m = base->cutnum;
#if 0
      /* FIXME: We'd like to create the user's problem description here,
	 but we can't because we don't have access to the user's LP
	 data structures. Right now, this means that we cannot answer
	 queries through the OSI interface. */

      userind = (int *) malloc((base->varnum + root->uind.size) * ISIZE);

      memcpy((char *)userind, (char *)base->userind, base->varnum * ISIZE);
      memcpy((char *)(userind + base->varnum), (char *)root->uind.list,
	     root->uind.size * ISIZE);

      user_create_subproblem(env->user, userind, env->mip, &maxn, &maxm,
			    &maxnz);
#endif
      break;

    case USER_DEFAULT:

       root->uind.size = env->mip->n;
       base->cutnum = env->mip->m;

#if 0
      if (env->mip->n && env->mip->m){
	 root->uind.size = env->mip->n;
	 base->cutnum = env->mip->m;
      }else if (!root->uind.size){
	 printf("Error setting up the root node.\n");
	 printf("User did not specify number of variables. Exiting.\n\n");
	 return(FUNCTION_TERMINATED_ABNORMALLY);
      }else if (!base->varnum){
	 printf("Error setting up the root node.\n");
	 printf("User did not specify number of base constraints.",
		"Exiting.\n\n");
	 return(FUNCTION_TERMINATED_ABNORMALLY);
      }
#endif
      if(root->uind.size){
	 root->uind.list = (int *) malloc(root->uind.size * ISIZE);
	 for (i = 0; i < root->uind.size; i++){
	    root->uind.list[i] = i;
	 }
      }

      base->varnum = 0;
      base->userind = NULL;

#if 0
      if (base->varnum){
	 qsort_i(base->userind, base->varnum);
      }
      if (root->uind.size && !env->par.warm_start){
	 qsort_i(root->uind.list, root->uind.size);
      }
#endif

      break;

    default:
      break;
   }

   if (env->par.warm_start){
      root->uind.size = 0;
      FREE(root->uind.list);
      return(FUNCTION_TERMINATED_NORMALLY);
   }

   root->uind.type = EXPLICIT_LIST;
   root->cutind.type = EXPLICIT_LIST;
   root->not_fixed.type = EXPLICIT_LIST;
   root->basis.basis_exists = FALSE;
   root->nf_status = NF_CHECK_NOTHING;
   root->nf_status = (env->par.tm_par.colgen_strat[0] & COLGEN__FATHOM) ?
                      NF_CHECK_ALL : NF_CHECK_NOTHING;

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int receive_feasible_solution_u(sym_environment *env, int msgtag)
{
   receive_int_array(&(env->best_sol.xlevel), 1);
   receive_int_array(&(env->best_sol.xindex), 1);
   receive_int_array(&(env->best_sol.xiter_num), 1);
   receive_dbl_array(&(env->best_sol.lpetol), 1);
   receive_dbl_array(&(env->best_sol.objval), 1);
   receive_int_array(&(env->best_sol.xlength), 1);
   if (env->best_sol.xlength > 0){
      FREE(env->best_sol.xind);
      FREE(env->best_sol.xval);
      env->best_sol.xind = (int *) malloc(env->best_sol.xlength * ISIZE);
      env->best_sol.xval = (double *) malloc(env->best_sol.xlength * DSIZE);
      receive_int_array(env->best_sol.xind, env->best_sol.xlength);
      receive_dbl_array(env->best_sol.xval, env->best_sol.xlength);
   }
   if (!env->has_ub || env->best_sol.objval < env->ub){
      env->has_ub = TRUE;
      env->ub = env->best_sol.objval;
   }
   env->best_sol.has_sol = TRUE;

   switch (msgtag){
    case FEASIBLE_SOLUTION_NONZEROS:
      break;

    case FEASIBLE_SOLUTION_USER:
      /* A feasible solution has been found in the LP process, and
       * it was packed by the user */

#ifdef USE_SYM_APPLICATION
      CALL_USER_FUNCTION( user_receive_feasible_solution(env->user, msgtag,
							 env->best_sol.objval,
							 env->best_sol.xlength,
							 env->best_sol.xind,
							 env->best_sol.xval) );
#endif
      break;
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int send_lp_data_u(sym_environment *env, int sender)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   int i;
   tm_prob *tm = env->tm;
   tm->par.max_active_nodes = env->par.tm_par.max_active_nodes;
#ifdef _OPENMP
   omp_set_dynamic(FALSE);
   omp_set_num_threads(tm->par.max_active_nodes);
#else
   tm->par.max_active_nodes = 1;
#endif

   tm->lpp = (lp_prob **) malloc(tm->par.max_active_nodes * sizeof(lp_prob *));

#pragma omp parallel for
   for (i = 0; i < tm->par.max_active_nodes; i ++){
      tm->lpp[i] = (lp_prob *) calloc(1, sizeof(lp_prob));
      tm->lpp[i]->proc_index = i;
      tm->lpp[i]->par = env->par.lp_par;

      if ((tm->lpp[i]->has_ub = env->has_ub)){
	 tm->lpp[i]->ub = env->ub;
      }else{
	 env->ub = - (MAXDOUBLE / 2);
      }
      if (env->par.multi_criteria){
	 if ((tm->lpp[i]->has_mc_ub = env->has_mc_ub)){
	    tm->lpp[i]->mc_ub = env->mc_ub;
	    tm->lpp[i]->obj[0] = env->obj[0];
	    tm->lpp[i]->obj[1] = env->obj[1];
	 }else{
	    env->mc_ub = - (MAXDOUBLE / 2);
	 }
	 tm->lpp[i]->utopia[0] = env->utopia[0];
	 tm->lpp[i]->utopia[1] = env->utopia[1];
      }
      tm->lpp[i]->draw_graph = env->dg_tid;
      tm->lpp[i]->base = *(env->base);
      tm->lpp[i]->mip = env->mip;
   }

#ifdef USE_SYM_APPLICATION
   for (i = 0; i < tm->par.max_active_nodes; i ++){
      CALL_USER_FUNCTION( user_send_lp_data(env->user, &(tm->lpp[i]->user)) );
   }
#endif

#else
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&env->par.lp_par), sizeof(lp_params));
   send_int_array(&env->has_ub, 1);
   if (env->has_ub)
      send_dbl_array(&env->ub, 1);
   if (env->par.multi_criteria){
      send_int_array(&env->has_mc_ub, 1);
      if (env->has_mc_ub){
	 send_dbl_array(&env->mc_ub, 1);
	 send_dbl_array(env->obj, 2);
      }
      send_dbl_array(env->utopia, 2);
   }
   send_int_array(&env->dg_tid, 1);
   send_int_array(&env->base->varnum, 1);
   if (env->base->varnum){
      send_int_array(env->base->userind, env->base->varnum);
   }
   send_int_array(&env->base->cutnum, 1);
   MIPdesc *mip = env->mip;
   send_int_array(&(mip->m), 1);
   send_int_array(&(mip->n), 1);
   send_int_array(&(mip->nz), 1);
   send_char_array(&(mip->obj_sense), 1);
   send_dbl_array(&(mip->obj_offset), 1);
   if (env->mip->matbeg){
      char has_desc = TRUE;
      char has_colnames = FALSE;
      send_char_array(&has_desc, 1);
      send_int_array(mip->matbeg, mip->n+1);
      send_int_array(mip->matind, mip->nz);
      send_dbl_array(mip->matval, mip->nz);
      send_dbl_array(mip->obj, mip->n);
      if (env->par.multi_criteria){
	 send_dbl_array(mip->obj, mip->n);
	 send_dbl_array(mip->obj2, mip->n);
      }
      send_dbl_array(mip->rhs, mip->m);
      send_char_array(mip->sense, mip->m);
      send_dbl_array(mip->rngval, mip->m);
      send_dbl_array(mip->ub, mip->n);
      send_dbl_array(mip->lb, mip->n);
      send_char_array(mip->is_int, mip->n);
      if (mip->colname){
	 int i;
	 has_colnames = TRUE;
	 send_char_array(&has_colnames, 1);
	 for (i = 0; i < mip->n; i++){
	    send_char_array(mip->colname[i], 8);
	 }
      }else{
	 send_char_array(&has_colnames, 1);
      }
   }else{
      char has_desc = FALSE;
      send_char_array(&has_desc, 1);
   }
#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_send_lp_data(env->user, NULL) );
#endif
   send_msg(sender, LP_DATA);
   freebuf(s_bufid);
#endif

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int send_cg_data_u(sym_environment *env, int sender)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined(COMPILE_IN_CG)
   int i;
   tm_prob *tm = env->tm;
   tm->cgp = (cg_prob **) malloc(tm->par.max_active_nodes*sizeof(cg_prob *));
#pragma omp parallel for
   for (i = 0; i < tm->par.max_active_nodes; i++){
      tm->lpp[i]->cgp = tm->cgp[i] = (cg_prob *) calloc(1, sizeof(cg_prob));

      tm->cgp[i]->par = env->par.cg_par;

      tm->cgp[i]->draw_graph = env->dg_tid;
   }
#ifdef USE_SYM_APPLICATION
   for (i = 0; i < tm->par.max_active_nodes; i++){
      CALL_USER_FUNCTION( user_send_cg_data(env->user,
					    &(tm->lpp[i]->cgp->user)) );
   }
#endif
#else
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&env->par.cg_par), sizeof(cg_params));
   send_int_array(&env->dg_tid, 1);
#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_send_cg_data(env->user, NULL) );
#endif
   send_msg(sender, CG_DATA);
   freebuf(s_bufid);
#endif

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int send_cp_data_u(sym_environment *env, int sender)
{
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_CP)
   int i;
   tm_prob *tm = env->tm;

   tm->cpp = (cut_pool **) malloc(env->par.tm_par.max_cp_num*sizeof(cut_pool *));
   for (i = 0; i < env->par.tm_par.max_cp_num; i++){
      tm->cpp[i] = (cut_pool *) calloc(1, sizeof(cut_pool));
      tm->cpp[i]->par = env->par.cp_par;
#ifdef USE_SYM_APPLICATION
      CALL_USER_FUNCTION( user_send_cp_data(env->user, &env->tm->cpp[i]->user) );
#endif
   }
#else
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_char_array((char *)(&env->par.cp_par), sizeof(cp_params));
#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_send_cp_data(env->user, NULL) );
#endif
   send_msg(sender, CP_DATA);
   freebuf(s_bufid);
#endif

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int display_solution_u(sym_environment *env, int thread_num)
{
   int user_res, i;
   lp_sol sol;

   memset(&sol, 0, sizeof(lp_sol));

   sol.xlength = 0;

   if (env->par.verbosity < -1){
       return(FUNCTION_TERMINATED_NORMALLY);
   }

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   if (env->tm && env->tm->lpp[thread_num]){
      sol = env->tm->lpp[thread_num]->best_sol;
      if (env->par.multi_criteria){
	 env->obj[0] = env->tm->lpp[thread_num]->obj[0];
	 env->obj[1] = env->tm->lpp[thread_num]->obj[1];
      }
   }
#else
   sol = env->best_sol;
#endif

   if (!sol.has_sol){
      switch(env->termcode){
       case TM_UNBOUNDED:
	 printf("\nThe problem is unbounded!\n\n");
	 return(FUNCTION_TERMINATED_NORMALLY);
       case TM_NO_SOLUTION:
	  printf("\nThe problem is infeasible!");
	  break;
       default:
	  break;
      }
      printf("\nNo Solution Found\n\n");
      return(FUNCTION_TERMINATED_NORMALLY);
   }

   printf("\nSolution Found: Node %i, Level %i\n", sol.xindex, sol.xlevel);
   if (env->par.multi_criteria){
      printf("First Objective: %.10f\n", env->obj[0]);
      printf("Second Objective: %.10f\n", env->obj[1]);
   }else{
      printf("Solution Cost: %.10f\n", env->mip->obj_sense == SYM_MINIMIZE ?
	     sol.objval + env->mip->obj_offset :
	     -sol.objval + env->mip->obj_offset);
   }
   qsort_id(sol.xind, sol.xval, sol.xlength);

#ifdef USE_SYM_APPLICATION
   user_res = user_display_solution(env->user, sol.lpetol, sol.xlength,
				    sol.xind, sol.xval, sol.objval);

#else
   user_res = USER_DEFAULT;
#endif

   if (env->par.verbosity > -1){
      switch(user_res){
       case USER_SUCCESS:
	  return(FUNCTION_TERMINATED_NORMALLY);
       case USER_DEFAULT:
	  if (sol.xlength){
	     if (env->mip->colname){
		printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
		printf("Column names and values of nonzeros in the solution\n");
		printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
		for (i = 0; i < sol.xlength; i++){
		   if (sol.xind[i] == env->mip->n){
		      continue;
		   }
		   printf("%8s %10.10f\n", env->mip->colname[sol.xind[i]],
			  sol.xval[i]);
		}
		for (i = 0; i < env->mip->fixed_n; i++){
		   printf("%8s %10.10f\n",
			  env->orig_mip->colname[env->mip->fixed_ind[i]],
			  env->mip->fixed_val[i]);
		}

		printf("\n");
	     }else{
		printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
		printf("User indices and values of nonzeros in the solution\n");
		printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
		for (i = 0; i < sol.xlength; i++){
		   if (sol.xind[i] == env->mip->n){
		      continue;
		   }
		   if(!env->prep_mip){
		      printf("%7d %10.10f\n", sol.xind[i], sol.xval[i]);
		   }else{
		      printf("%7d %10.10f\n",
			     env->prep_mip->orig_ind[sol.xind[i]],
			     sol.xval[i]);
		   }
		}
		for (i = 0; i < env->mip->fixed_n; i++){
		   printf("%7d %10.10f\n", env->mip->fixed_ind[i],
			  env->mip->fixed_val[i]);
		}
		printf("\n");
	     }

	     return(FUNCTION_TERMINATED_NORMALLY);
	  }else{
	     printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	     printf("All columns are zero in the solution!\n");
	     printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	     return(FUNCTION_TERMINATED_NORMALLY);
	  }
       case USER_ERROR:
	  return(FUNCTION_TERMINATED_NORMALLY);

       default:
	  return(FUNCTION_TERMINATED_NORMALLY);
      }
   }

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int process_own_messages_u(sym_environment *env, int msgtag)
{
#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_process_own_messages(env->user, msgtag) );
#endif
   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

int free_master_u(sym_environment *env)
{
   int i;

#ifdef USE_SYM_APPLICATION
   CALL_USER_FUNCTION( user_free_master(&env->user) );
#endif
   FREE(env->best_sol.xind);
   FREE(env->best_sol.xval);

   if (env->mip){
      free_mip_desc(env->mip);
      FREE(env->mip);
   }

   if(env->prep_mip){
      free_mip_desc(env->prep_mip);
      FREE(env->prep_mip);
   }

   if (env->rootdesc){
      FREE(env->rootdesc->desc);
      FREE(env->rootdesc->uind.list);
      FREE(env->rootdesc->not_fixed.list);
      FREE(env->rootdesc->cutind.list);
      FREE(env->rootdesc);
   }

   if (env->base){
      FREE(env->base->userind);
      FREE(env->base);
   }

#ifdef COMPILE_IN_TM
   if (env->warm_start){
      free_subtree(env->warm_start->rootnode);
      if(env->warm_start->best_sol.has_sol){
	 FREE(env->warm_start->best_sol.xind);
	 FREE(env->warm_start->best_sol.xval);
      }
      if (env->warm_start->cuts){
	 for (i = env->warm_start->cut_num - 1; i >= 0; i--){
	    if (env->warm_start->cuts[i]){
	       FREE(env->warm_start->cuts[i]->coef);
	    }
	    FREE(env->warm_start->cuts[i]);
	 }
      }

      FREE(env->warm_start->cuts);
      FREE(env->warm_start);
   }
#ifdef COMPILE_IN_CP
   if (env->cp){
      for (i = 0; i < env->par.tm_par.max_cp_num; i++){
	 env->cp[i]->msgtag = YOU_CAN_DIE;
	 cp_close(env->cp[i]);
      }
      FREE(env->cp);
   }
#endif
#endif

   return(FUNCTION_TERMINATED_NORMALLY);
}

/*===========================================================================*/

