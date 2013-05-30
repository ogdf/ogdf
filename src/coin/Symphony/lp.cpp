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

#ifndef COMPILE_IN_LP

#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "sym_lp.h"
#include "sym_proccomm.h"
#include "sym_messages.h"
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_types.h"
#include "sym_lp_solver.h"

/*===========================================================================*/

/*===========================================================================*\
 * This is the main() that is used if the LP is running as a separate
 * process. This file is only used in that case.
\*===========================================================================*/

int main(void)
{
   lp_prob *p;
   int r_bufid;
   double time, diff;
   struct timeval timeout = {10, 0};
   char first_node_rec = FALSE;
   int termcode;

   p = (lp_prob *) calloc(1, sizeof(lp_prob));

   p->start_time = wall_clock(NULL);

   if ((termcode = lp_initialize(p, 0)) < 0){
      printf("LP initialization failed with error code %i\n\n", termcode);
      lp_exit(p);
   }

   /*------------------------------------------------------------------------*\
    * Continue receiving node data and fathoming branches until this
    * process is killed
   \*------------------------------------------------------------------------*/

   p->phase = 0;
   while (TRUE){
      p->lp_data->col_set_changed = TRUE;
      /*---------------------------------------------------------------------*\
       * waits for an active node message but if there's anything left after
       * receiving that, those messages are processed, before going to
       * process_chain().
      \*---------------------------------------------------------------------*/
      time = wall_clock(NULL);
      do{
	 r_bufid = treceive_msg(ANYONE, ANYTHING, &timeout);
      }while (! process_message(p, r_bufid, NULL, NULL) );
      diff = wall_clock(NULL) - time;
      if (first_node_rec){
	 p->comp_times.idle_node += diff;
      }else{
	 first_node_rec = TRUE;
	 p->comp_times.ramp_up_lp += diff;
      }
      do{
	 r_bufid = nreceive_msg(ANYONE, ANYTHING);
	 if (r_bufid)
	    process_message(p, r_bufid, NULL, NULL);
      }while (r_bufid);

      p->comp_times.communication += used_time(&p->tt);

      if (process_chain(p) < 0){
	 printf("\nThere was an error in the LP process. Exiting now.\n\n");
	 /* There was an error in the LP. Abandon node. */
	 lp_exit(p);
      }
   }

   p->comp_times.wall_clock_lp = wall_clock(NULL) - p->start_time;

   lp_exit(p);

   return(0);
}

#endif
