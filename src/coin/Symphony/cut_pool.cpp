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

#ifndef COMPILE_IN_CP

#include <stdio.h>
#include <stdlib.h>

#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_timemeas.h"
#include "sym_proccomm.h"
#include "sym_messages.h"
#include "sym_cp.h"

/*===========================================================================*/

/*===========================================================================*\
 * This is the main() that is used if the CP is running as a separate
 * process. This file is only used in that case.
\*===========================================================================*/

int main(void)
{
   cut_pool *cp;
   int s_bufid, r_bufid;
   int num_cuts = 0;
   double tt = 0, get_cuts_time;
   struct timeval tout = {10, 0};

   cp = (cut_pool *) calloc(1, sizeof(cut_pool));

   cp_initialize(cp, 0);

   (void) used_time(&tt);

   /*------------------------------------------------------------------------*\
    * The main loop -- this keeps executing until the process dies
   \*------------------------------------------------------------------------*/

   while(TRUE){
      do{
	 r_bufid = treceive_msg(ANYONE, ANYTHING, &tout);
	 if (!r_bufid){
	    if (pstat(cp->tree_manager) != PROCESS_OK){
	       printf("TM has died -- CP exiting\n\n");
	       exit(-601);
	    }
	 }
      }while (! r_bufid);
      cp_process_message(cp, r_bufid);
      if (cp->msgtag==LP_SOLUTION_NONZEROS || cp->msgtag==LP_SOLUTION_USER ||
	  cp->msgtag==LP_SOLUTION_FRACTIONS){

	 num_cuts = check_cuts_u(cp, &cp->cur_sol);

	 if (cp->par.check_which == CHECK_ALL_CUTS ||
	     cp->par.check_which == CHECK_LEVEL ||
	     cp->par.check_which == CHECK_TOUCHES ||
	     cp->par.check_which == CHECK_LEVEL_AND_TOUCHES){
	    get_cuts_time = used_time(&tt);
	    s_bufid = init_send(DataInPlace);
	    send_int_array(&num_cuts, 1);
	    send_dbl_array(&get_cuts_time, 1);
	    send_int_array(&cp->cur_sol.xindex, 1);
	    send_int_array(&cp->cur_sol.xiter_num, 1);
	    send_msg(cp->cur_sol.lp, NO_MORE_CUTS);
	    freebuf(s_bufid);
	 }

	 if (++cp->reorder_count % 10 == 0){
	    delete_duplicate_cuts(cp);
	    order_cuts_by_quality(cp);
	    cp->reorder_count = 0;
	 }

	 FREE(cp->cur_sol.xind);
	 FREE(cp->cur_sol.xval);
      }
   }
   FREE(cp);

   return(0);
}

#endif
