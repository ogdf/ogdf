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

#include <stdio.h>
#include <stdlib.h>

#include "sym_macros.h"
#include "sym_timemeas.h"
#include "sym_proccomm.h"
#include "sym_qsort.h"
#include "sym_messages.h"
#include "sym_pack_cut.h"
#include "sym_cp.h"
#include "sym_constants.h"

#ifdef CHAR_IS_SIGNED
#define MEMCMP(c0, c1, s) unsigned_memcmp(c0, c1, s)
#else
#include <memory.h>
#define MEMCMP(c0, c1, s) memcmp(c0, c1, s)
#endif

/*===========================================================================*/

/*===========================================================================*\
 * This file contains general functions used by the cut pool process.
\*===========================================================================*/

/*===========================================================================*/

void cp_initialize(cut_pool *cp, int master_tid)
{
#ifndef COMPILE_IN_CP
   int bytes;
#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   int s_bufid, r_bufid;
#endif
#endif
#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_LP)
   int s_bufid, r_bufid;
#endif

   /*------------------------------------------------------------------------*\
    * Receive tid info; request and receive problem specific data
   \*------------------------------------------------------------------------*/

#ifdef COMPILE_IN_CP

   cp->master = master_tid;

#else /*We only need to do the rest of this if the CP is running as a separate
	process*/

   /* set stdout to be line buffered */
   setvbuf(stdout, (char *)NULL, _IOLBF, 0);

   register_process();

   r_bufid = receive_msg(ANYONE, MASTER_TID_INFO);
   bufinfo(r_bufid, &bytes, &cp->msgtag, &cp->tree_manager);
   receive_int_array(&cp->master, 1);
   freebuf(r_bufid);

#endif

#if !defined(COMPILE_IN_TM) || !defined(COMPILE_IN_CP)

   s_bufid = init_send(DataInPlace);
   send_msg(cp->master, REQUEST_FOR_CP_DATA);
   freebuf(s_bufid);

   receive_cp_data_u(cp);

#endif

   if (cp->par.warm_start == READ_CP_LIST){
      read_cp_cut_list(cp, cp->par.warm_start_file_name);
   }else if (cp->par.warm_start == READ_TM_LIST){
      cp_read_tm_cut_list(cp, cp->par.warm_start_file_name);
   }else if (!cp->cuts){
      cp->cuts = (cp_cut_data **) calloc (cp->par.block_size,
					  sizeof(cp_cut_data *));
      cp->allocated_cut_num = cp->par.block_size;
   }
}

/*===========================================================================*/

int unsigned_memcmp(char *coef0, char *coef1, int size)
{
   register char *end0 = coef0 + size;

   for ( ; coef0 != end0; coef0++, coef1++)
      if (*coef0 != *coef1)
	 break;
   if (coef0 == end0)
      return(0);
   return ( (unsigned char)(*coef0) < (unsigned char)(*coef1) ? -1 : 1);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function compares the quality of two cuts.
\*===========================================================================*/

int cut_quality_cmp(const void *cut0ptr, const void *cut1ptr)
{
   cp_cut_data *cut0 = *((cp_cut_data **)cut0ptr);
   cp_cut_data *cut1 = *((cp_cut_data **)cut1ptr);

   return((int)(1000*(cut1->quality - cut0->quality)));
}

/*===========================================================================*/

/*===========================================================================*\
 * This function orders the cuts in the pool by the numerical measure of
 * "quality".
\*===========================================================================*/

void order_cuts_by_quality(cut_pool *cp)
{
   cp_cut_data **cuts = cp->cuts;

   /* order the cuts according to the function "cutcmp" */
   qsort(cuts, cp->cut_num, sizeof(cp_cut_data *), cut_quality_cmp);
}

/*===========================================================================*/

/*===========================================================================*\
 * For purposes of deleting duplicate cuts, this function compares two cuts
 * to see if they are the same or not. If the cuts are the same then the
 * return value is 0, otherwise nonzero.
\*===========================================================================*/

int cutcmp(const void *cut0ptr, const void *cut1ptr)
{
   cut_data *cut0 = *((cut_data **)cut0ptr);
   cut_data *cut1 = *((cut_data **)cut1ptr);
   int typediff, sizediff;

   if ((typediff = cut0->type - cut1->type))  return(typediff);
   if ((sizediff = cut0->size - cut1->size))  return(sizediff);
   return( MEMCMP(cut0->coef, cut1->coef, cut0->size) );
}

/*===========================================================================*/

int delete_ineffective_cuts(cut_pool *cp)
{
   cp_cut_data **cuts = cp->cuts;
   int num;
   int del_cuts = 0, tmp_del_cuts = 0, cuts_to_leave = 0;
   cp_cut_data **cp_cut1, **cp_cut2;
   int touches_until_deletion = cp->par.touches_until_deletion;
   int min_to_delete = cp->par.min_to_delete;

   if (min_to_delete > cp->cut_num)
      min_to_delete = (int)(0.2*cp->cut_num);

   switch (cp->par.delete_which){

    case DELETE_BY_QUALITY:

      order_cuts_by_quality(cp);

      cuts_to_leave = MIN(cp->par.cuts_to_check, cp->cut_num-min_to_delete);

      for (cp_cut1 = cuts + cuts_to_leave, num = cuts_to_leave;
	   num < cp->cut_num; cp_cut1++, num++){
	 del_cuts++;
	 cp->size -= (*cp_cut1)->cut.size;
	 FREE((*cp_cut1)->cut.coef);
	 FREE(*cp_cut1);
      }
      cp->cut_num -= del_cuts;
      cp->size -= del_cuts * (int) sizeof(cp_cut_data);
      break;

    case DELETE_BY_TOUCHES:
    default:

      while (del_cuts < min_to_delete){
	 for (tmp_del_cuts = 0, cp_cut1 = cp_cut2 = cuts,
		 num = cp->cut_num; num > 0; cp_cut2++, num--){
	    if ((*cp_cut2)->touches >= touches_until_deletion){
	       tmp_del_cuts++;
	       cp->size -= (*cp_cut2)->cut.size;
	       FREE((*cp_cut2)->cut.coef);
	       FREE(*cp_cut2);
	    }else{
	       *cp_cut1 = *cp_cut2;
	       cp_cut1++;
	    }
	 }
	 cp->cut_num -= tmp_del_cuts;
	 cp->size -= tmp_del_cuts * (int)sizeof(cp_cut_data);
	 del_cuts += tmp_del_cuts;
	 touches_until_deletion--;
      }
      break;

   }

   if (cp->par.verbosity > 5)
      printf("******* CUT_POOL : Deleted %i ineffective cuts leaving %i\n",
	     del_cuts, cp->cut_num);

   return(del_cuts);
}

/*===========================================================================*/

int delete_duplicate_cuts(cut_pool *cp)
{
   cp_cut_data **cuts = cp->cuts;
   int num;
   int del_cuts = 0;
   cp_cut_data **cp_cut1, **cp_cut2;
   int touches, level;

   /* order the cuts according to the function "cutcmp" */
   qsort(cuts, cp->cut_num, sizeof(cp_cut_data *), cutcmp);
   /* go through and remove duplicates */
   for(num = cp->cut_num-1, cp_cut1 = cuts, cp_cut2 = cp_cut1 + 1;
       num > 0; cp_cut2++, num--){
      switch (which_cut_to_delete(&(*cp_cut1)->cut, &(*cp_cut2)->cut)){
       case 0:
	 cp_cut1++;
	 *cp_cut1 = *cp_cut2;
	 break;
       case 1:
	 del_cuts++;
	 cp->size -= (*cp_cut1)->cut.size;
	 touches = MIN((*cp_cut1)->touches, (*cp_cut2)->touches);
	 level = MIN((*cp_cut1)->level, (*cp_cut2)->level);
	 FREE((*cp_cut1)->cut.coef);
	 FREE(*cp_cut1);
	 *cp_cut1 = *cp_cut2;
	 (*cp_cut1)->touches = touches;
	 (*cp_cut1)->level = level;
	 break;
       case 2:
	 del_cuts++;
	 cp->size -= (*cp_cut2)->cut.size;
	 touches = MIN((*cp_cut1)->touches, (*cp_cut2)->touches);
	 level = MIN((*cp_cut1)->level, (*cp_cut2)->level);
	 FREE((*cp_cut2)->cut.coef);
	 FREE(*cp_cut2);
	 (*cp_cut1)->touches = touches;
	 (*cp_cut1)->level = level;
	 break;
      }
   }

   cp->cut_num -= del_cuts;
   cp->size -= del_cuts * (int)sizeof(cp_cut_data);

   if (cp->par.verbosity > 5)
      printf("******* CUT_POOL : Deleted %i duplicate cuts leaving %i\n",
	     del_cuts, cp->cut_num);
   return(del_cuts);
}

/*===========================================================================*/

int which_cut_to_delete(cut_data *cut1, cut_data *cut2)
{
   if (cutcmp(&cut1, &cut2))
      return(0);

   return(cut1->sense == 'E' ? 2 :
	  cut2->sense == 'E' ? 1 :
	  cut1->sense != cut2->sense ? 0 :
	  cut1->sense == 'R' ? 0 :
	  cut1->sense == 'L' ? (cut1->rhs<=cut2->rhs ? 2:1) :
			       (cut1->rhs>=cut2->rhs ? 2:1));
}

/*===========================================================================*/

int write_cp_cut_list(cut_pool *cp, char *file, char append)
{
   FILE *f;
   int i, j;

   if (!(f = fopen(file, append ? "a" : "w"))){
      printf("\nError opening cut file\n\n");
      return(0);
   }

   fprintf(f, "CUTNUM: %i %i %i\n", cp->allocated_cut_num, cp->cut_num,
	   cp->size);
   for (i = 0; i < cp->cut_num; i++){
      fprintf(f, "%i %i %i %i %i %c %i %f %f\n", cp->cuts[i]->touches,
	      cp->cuts[i]->level, cp->cuts[i]->cut.name, cp->cuts[i]->cut.size,
	      (int)cp->cuts[i]->cut.type, cp->cuts[i]->cut.sense,
	      (int)cp->cuts[i]->cut.branch, cp->cuts[i]->cut.rhs,
	      cp->cuts[i]->cut.range);
      for (j = 0; j < cp->cuts[i]->cut.size; j++)
	 fprintf(f, "%i ", (int)cp->cuts[i]->cut.coef[j]);
      fprintf(f, "\n");
   }

   fclose(f);

   return(1);
}

/*===========================================================================*/

int read_cp_cut_list(cut_pool *cp, char *file)
{
   FILE *f;
   int i, j, tmp1 = 0, tmp2 = 0;
   char str[20];

   if (!(f = fopen(file, "r"))){
      printf("\nError opening cut file\n\n");
      return(0);
   }

   fscanf(f, "%s %i %i %i", str, &cp->allocated_cut_num, &cp->cut_num,
	  &cp->size);
   cp->cuts = (cp_cut_data **)
      malloc(cp->allocated_cut_num*sizeof(cp_cut_data *));
   for (i = 0; i < cp->cut_num; i++){
      cp->cuts[i] = (cp_cut_data *) malloc(sizeof(cp_cut_data));
      fscanf(f, "%i %i %i %i %i %c %i %lf %lf", &cp->cuts[i]->touches,
	     &cp->cuts[i]->level, &cp->cuts[i]->cut.name,
	     &cp->cuts[i]->cut.size, &tmp1, &cp->cuts[i]->cut.sense, &tmp2,
	     &cp->cuts[i]->cut.rhs, &cp->cuts[i]->cut.range);
      cp->cuts[i]->cut.type = (char) tmp1;
      cp->cuts[i]->cut.branch = (char) tmp2;
      cp->cuts[i]->cut.coef =
	 (char *) malloc(cp->cuts[i]->cut.size*sizeof(char));
      for (j = 0; j < cp->cuts[i]->cut.size; j++){
	 fscanf(f, "%i ", &tmp1);
	 cp->cuts[i]->cut.coef[j] = (char) tmp1;
      }
   }

   fclose(f);

   return(1);
}

/*===========================================================================*/

int cp_read_tm_cut_list(cut_pool *cp, char *file)
{
   FILE *f;
   int i, j, tmp1 = 0, tmp2 = 0;
   char str[20];

   if (!(f = fopen(file, "r"))){
      printf("\nError opening cut file\n\n");
      return(0);
   }

   cp->size = 0;
   fscanf(f, "%s %i %i", str, &cp->cut_num, &cp->allocated_cut_num);
   cp->cuts = (cp_cut_data **)
      malloc(cp->allocated_cut_num*sizeof(cp_cut_data *));
   for (i = 0; i < cp->cut_num; i++){
      cp->cuts[i] = (cp_cut_data *) calloc(1, sizeof(cp_cut_data));
      fscanf(f, "%i %i %i %c %i %lf %lf", &cp->cuts[i]->cut.name,
	     &cp->cuts[i]->cut.size, &tmp1, &cp->cuts[i]->cut.sense,
	     &tmp2, &cp->cuts[i]->cut.rhs, &cp->cuts[i]->cut.range);
      cp->cuts[i]->cut.type = (char)tmp1;
      cp->cuts[i]->cut.branch = (char)tmp2;
      cp->cuts[i]->cut.coef =
	(char *) malloc(cp->cuts[i]->cut.size*sizeof(char));
      cp->size += cp->cuts[i]->cut.size + (int) sizeof(cp_cut_data);
      for (j = 0; j < cp->cuts[i]->cut.size; j++){
	 fscanf(f, "%i ", &tmp1);
	 cp->cuts[i]->cut.coef[j] = (char)tmp1;
      }
   }

   fclose(f);

   return(1);
}

/*===========================================================================*/

void cp_close(cut_pool *cp)
{
#ifndef COMPILE_IN_CP
   int s_bufid;

   s_bufid = init_send(DataInPlace);
   send_dbl_array(&cp->cut_pool_time, 1);
   send_int_array(&cp->total_cut_num, 1);
   send_msg(cp->tree_manager, POOL_TIME);
   freebuf(s_bufid);
   if(cp->msgtag == YOU_CAN_DIE)
      free_cut_pool_u(cp);
#else
   FREE(cp->cuts_to_add);
   free_cut_pool_u(cp);
#endif
}
