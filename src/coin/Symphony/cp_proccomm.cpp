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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "sym_macros.h"
#include "sym_constants.h"
#include "sym_proccomm.h"
#include "sym_messages.h"
#include "sym_cp.h"
#include "sym_pack_cut.h"
#include "sym_timemeas.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains functions related to CP process communication.
\*===========================================================================*/

void cp_process_message(cut_pool *cp, int r_bufid)
{
   int s_bufid, bytes, new_tid;
   int size, i;
   char *buf, *bufc;
   cp_cut_data *cp_cut;
   double tt= 0;
   static struct timeval tout = {10, 0};
   int termcode = 0;

   bufinfo(r_bufid, &bytes, &cp->msgtag, &cp->cur_sol.lp);

   switch (cp->msgtag){
    case LP_SOLUTION_USER:
      cp->cut_pool_time += used_time(&tt);
      receive_int_array(&cp->cur_sol.xlevel, 1);
      receive_int_array(&cp->cur_sol.xindex, 1);
      receive_int_array(&cp->cur_sol.xiter_num, 1);
      receive_dbl_array(&cp->cur_sol.lpetol, 1);
      if ((termcode = receive_lp_solution_cp_u(cp)) < 0){
	 printf("Warning: User error detected in cut pool\n\n");
      }
      break;

    case LP_SOLUTION_NONZEROS:
    case LP_SOLUTION_FRACTIONS:
      /* receive an lp solution */
      cp->cut_pool_time += used_time(&tt);
      receive_int_array(&cp->cur_sol.xlevel, 1);
      receive_int_array(&cp->cur_sol.xindex, 1);
      receive_int_array(&cp->cur_sol.xiter_num, 1);
      receive_dbl_array(&cp->cur_sol.lpetol, 1);

      receive_int_array(&cp->cur_sol.xlength, 1);
      cp->cur_sol.xind = (int *) malloc(cp->cur_sol.xlength * ISIZE);
      cp->cur_sol.xval = (double *) malloc(cp->cur_sol.xlength * DSIZE);
      receive_int_array(cp->cur_sol.xind, cp->cur_sol.xlength);
      receive_dbl_array(cp->cur_sol.xval, cp->cur_sol.xlength);
      break;

    case PACKED_CUTS_TO_CP: /* receive a set of new (packed) cut */
      cut_pool_receive_cuts(cp, 0);
      freebuf(r_bufid);
      break;

   case WRITE_LOG_FILE:
      freebuf(r_bufid);
      if (cp->par.logging)
	 write_cp_cut_list(cp, cp->par.log_file_name, FALSE);
      break;

    case POOL_COPY_YOURSELF:
      /* This means that this cut pool will be split into two new cut
	 pools, each servicing a different subtree */
      receive_int_array(&new_tid, 1);
      freebuf(r_bufid);
      size = cp->cut_num * (int)sizeof(cp_cut_data);
      for (i=0; i<cp->cut_num; i++)
	 size += cp->cuts[i]->cut.size;
      buf = (char *) calloc(size, sizeof(char));

      bufc = buf;
      for (i=0; i<cp->cut_num; i++){
	 cp_cut = cp->cuts[i];
	 memcpy(bufc, (char *)cp_cut, sizeof(cp_cut_data));
	 bufc += sizeof(cp_cut_data);
	 memcpy(bufc, cp_cut->cut.coef, cp_cut->cut.size);
	 bufc += cp_cut->cut.size;
      }
      s_bufid = init_send(DataInPlace);
      send_int_array(&cp->cut_num, 1);
      send_int_array(&size, 1);
      send_char_array(buf, size);
      send_msg(new_tid, CUTPOOL_COPY);
      freebuf(s_bufid);
      FREE(buf);
      break;

    case YOU_CANNOT_DIE:
    case YOU_CAN_DIE:
      /* The program is exiting -- send back info */
      cp->cut_pool_time += used_time(&tt);
      cp->total_cut_num += cp->cut_num;
      cp_close(cp);
      if (cp->msgtag == YOU_CANNOT_DIE)
	 break;
      comm_exit();
      exit(1);

    case POOL_YOU_ARE_USELESS:
      receive_int_array(&new_tid, 1);
      freebuf(r_bufid);
      s_bufid = init_send(DataInPlace);
      send_msg(cp->tree_manager, POOL_USELESSNESS_ACKNOWLEDGED);

      cp->cut_pool_time += used_time(&tt);
      cp->total_cut_num += cp->cut_num;
      for (i = cp->cut_num - 1; i >= 0; i--){
	 FREE(cp->cuts[i]->cut.coef);
	 FREE(cp->cuts[i]);
      }

      do{
	 treceive_msg(new_tid, CUTPOOL_COPY, &tout);
	 if (! r_bufid){
	    if (pstat(new_tid) != PROCESS_OK){
	       printf("Other CP has died -- CP exiting\n\n");
	       exit(-602);
	    }
	 }
      }while (! r_bufid);
      receive_int_array(&cp->cut_num, 1);
      receive_int_array(&cp->size, 1);
      buf = (char *) calloc(cp->size, sizeof(char));
      receive_char_array(buf, cp->size);
      freebuf(r_bufid);

      if (cp->allocated_cut_num < cp->cut_num){
	 cp->allocated_cut_num = cp->cut_num + cp->par.block_size;
	 FREE(cp->cuts);
	 cp->cuts = (cp_cut_data **) malloc(cp->allocated_cut_num *
						   sizeof(cp_cut_data*));
      }
      bufc = buf;
      for (i=0; i<cp->cut_num; i++){
	 cp_cut = cp->cuts[i] =
	    (cp_cut_data *) malloc( sizeof(cp_cut_data));
	 memcpy(cp_cut, bufc, sizeof(cp_cut_data));
	 bufc += sizeof(cp_cut_data);
	 cp_cut->cut.coef = (char *) malloc(cp_cut->cut.size * CSIZE);
	 memcpy(cp_cut->cut.coef, bufc, cp_cut->cut.size);
	 bufc += cp_cut->cut.size;
      }
      FREE(buf);
      break;

    default:
      printf("Unrecognized message type!!! \n\n");
      break;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Send a packed cut to another process
\*===========================================================================*/

void cut_pool_send_cut(cut_pool *cp, cut_data *new_cut, int tid)
{
#ifdef COMPILE_IN_CP

   cut_data *tmp_cut;

   tmp_cut = (cut_data *) malloc (sizeof(cut_data));
   memcpy((char *)tmp_cut, (char *)new_cut, sizeof(cut_data));
   tmp_cut->coef = (char *) malloc (new_cut->size * sizeof(char));
   memcpy(tmp_cut->coef, new_cut->coef, new_cut->size);
   REALLOC(cp->cuts_to_add, cut_data *, cp->cuts_to_add_size,
	   cp->cuts_to_add_num + 1, BB_BUNCH);
   cp->cuts_to_add[cp->cuts_to_add_num++] = tmp_cut;

#else

   int s_bufid;

   s_bufid = init_send(DataInPlace);
   pack_cut(new_cut);
   send_msg(tid, PACKED_CUT);
   freebuf(s_bufid);

#endif
}

/*===========================================================================*/

/*==========================================================================*\
 *This function receives and adds a cut to the list after checking to
 *see of memory reallocation is necessary.
\*==========================================================================*/

void cut_pool_receive_cuts(cut_pool *cp, int bc_level)
{
   int i, cnt, level;
   int del_cuts = 0, deleted_duplicates = FALSE;
   cp_cut_data *cp_cut;

#ifdef COMPILE_IN_CP
   cnt = cp->cuts_to_add_num;
#else
   receive_int_array(&cnt, 1);
#endif

   if (cnt + cp->cut_num > cp->allocated_cut_num &&
       (cnt > cp->par.block_size ||
	cnt > cp->par.max_number_of_cuts - cp->par.cuts_to_check)){
      printf("Too many cuts have arrived to CP. Forget it...\n");
      printf("  [ cnt: %i   bl_size: %i   max: %i ]\n\n",
	     cnt, cp->par.block_size, cp->par.max_number_of_cuts);
#ifdef COMPILE_IN_CP
      for (i = cnt - 1; i >= 0; i--)
	 FREE(cp->cuts_to_add[i]);
      cp->cuts_to_add_num = 0;
#endif
      return;
   }

   while (TRUE){
      /* This loop (despite its appearance) is not endless. In fact, it is
	 entered at most four times. Check it out :-). */
      if (cp->cut_num + cnt <= cp->allocated_cut_num)
	 break;

      if (cp->allocated_cut_num + cnt + cp->par.block_size <=
	  cp->par.max_number_of_cuts){ /* be greedy */
	 cp->allocated_cut_num += cnt + cp->par.block_size;
	 cp->cuts = (cp_cut_data **) realloc
	    (cp->cuts, cp->allocated_cut_num * sizeof(cp_cut_data *));
	 break;
      }else if (cp->cut_num + cnt + cp->par.block_size <=
		cp->par.max_number_of_cuts){ /* be less greedy ... */
	 cp->allocated_cut_num = cp->cut_num + cnt + cp->par.block_size;
	 cp->cuts = (cp_cut_data **) realloc
	    (cp->cuts, cp->allocated_cut_num * sizeof(cp_cut_data *));
	 break;
      }else if (cnt < cp->par.block_size &&
		cp->cut_num + cp->par.block_size<=cp->par.max_number_of_cuts){
	 cp->allocated_cut_num = cp->cut_num + cp->par.block_size;
	 cp->cuts = (cp_cut_data **) realloc
	    (cp->cuts, cp->allocated_cut_num * sizeof(cp_cut_data *));
	 break;
      }else{
	 /* If the maximum number of cuts allowed in the pool is exceeded,
	    then the pool is purged to make room */
	 if (!deleted_duplicates){
	    del_cuts += delete_duplicate_cuts(cp);
	    deleted_duplicates = TRUE;
	 }else{
	    del_cuts += delete_ineffective_cuts(cp);
	 }
	 printf("Max num of cuts in CP pool exceeded -- deleted %i cuts\n",
		del_cuts);
      }
   }

#ifdef COMPILE_IN_CP
   level = bc_level;
#else
   receive_int_array(&level, 1);
#endif
   for (i = cnt - 1; i >= 0; i--, del_cuts = 0){
      cp_cut = (cp_cut_data *) malloc( sizeof(cp_cut_data));
#ifdef COMPILE_IN_CP
      memcpy((char *)(&cp_cut->cut), (char *)cp->cuts_to_add[i],
	     sizeof(cut_data));
      if (cp_cut->cut.size >0){
	 cp_cut->cut.coef = (char *) ((int *)malloc (cp_cut->cut.size+ISIZE));
	 memcpy(cp_cut->cut.coef, cp->cuts_to_add[i]->coef,
		cp->cuts_to_add[i]->size*sizeof(char));
      }
      FREE(cp->cuts_to_add[i]->coef);
      FREE(cp->cuts_to_add[i]);
#else
      cp_cut->cut.coef = NULL;
      unpack_cut(&cp_cut->cut);
#endif
      cp_cut->level = level;
      cp_cut->touches = cp_cut->check_num = 0;
      cp_cut->quality = 0.0;
#if 0
      if (cp_cut->cut.rhs <= cp->lpetol){
	 printf("cut_pool: cut arrived with 0 rhs... \n\n");
      }
#endif
      if (cp->size + cp_cut->cut.size + sizeof(cp_cut_data) >
	  cp->par.max_size){
	 /* If the maximum size of the cut pool is exceeded, then we attempt
	    to delete some cuts to make room */
	 if (!deleted_duplicates){
	    del_cuts += delete_duplicate_cuts(cp);
	    deleted_duplicates = TRUE;
	 }
	 while (cp->size + cp_cut->cut.size + sizeof(cp_cut_data) >
		cp->par.max_size){
	    del_cuts += delete_ineffective_cuts(cp);
	 }

	 if (cp->par.verbosity > 4)
	    printf("Maximum CP size exceeded -- deleted %i cuts, leaving %i\n",
		   del_cuts, cp->cut_num);
      }

      cp->cuts[cp->cut_num++] = cp_cut;
      cp->size += cp_cut->cut.size + (int)sizeof(cp_cut_data);
   }
}

