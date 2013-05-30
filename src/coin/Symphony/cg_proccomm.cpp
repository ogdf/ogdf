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
#include <stdio.h>

#include "sym_macros.h"
#include "sym_constants.h"
#include "sym_pack_cut.h"
#include "sym_messages.h"
#include "sym_proccomm.h"
#include "sym_cg.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains function related to CG process communication.
\*===========================================================================*/

/*===========================================================================*\
 * Process a message from the queue -- this function is only used when the
 * cut generator is run as a separate process.
\*===========================================================================*/

int cg_process_message(cg_prob *p, int r_bufid)
{
   int bytes;

   bufinfo(r_bufid, &bytes, &p->msgtag, &p->cur_sol.lp);

   switch(p->msgtag){

    case YOU_CAN_DIE:
      cg_close(p);
      freebuf(r_bufid);
      comm_exit();
      exit(1);

    case LP_SOLUTION_NONZEROS:
    case LP_SOLUTION_FRACTIONS:
      /* receive a new LP solution for which cuts are to be generated */
      receive_int_array(&p->cur_sol.xlevel, 1);
      receive_int_array(&p->cur_sol.xindex, 1);
      receive_int_array(&p->cur_sol.xiter_num, 1);
      receive_dbl_array(&p->cur_sol.lpetol, 1);
      receive_dbl_array(&p->cur_sol.objval, 1);
      receive_char_array(&p->has_ub, 1);
      if (p->has_ub)
	 receive_dbl_array(&p->ub, 1);
      receive_int_array(&p->cur_sol.xlength, 1);
      REMALLOC(p->cur_sol.xind, int,
	       p->cur_sol.max_sol_length, p->cur_sol.xlength, BB_BUNCH);
      REMALLOC(p->cur_sol.xval, double,
	       p->cur_sol.max_sol_length, p->cur_sol.xlength, BB_BUNCH);
      receive_int_array(p->cur_sol.xind, p->cur_sol.xlength);
      receive_dbl_array(p->cur_sol.xval, p->cur_sol.xlength);
      freebuf(r_bufid);
      break;

    case LP_SOLUTION_USER:
      receive_int_array(&p->cur_sol.xlevel, 1);
      receive_int_array(&p->cur_sol.xindex, 1);
      receive_int_array(&p->cur_sol.xiter_num, 1);
      receive_dbl_array(&p->cur_sol.lpetol, 1);
      receive_dbl_array(&p->cur_sol.objval, 1);
      receive_char_array(&p->has_ub, 1);
      if (p->has_ub)
	 receive_dbl_array(&p->ub, 1);
      if (receive_lp_solution_cg_u(p) == USER_ERROR)
	 return(USER_ERROR);
      break;

    default:
      printf("Unrecognized message type %i from %i!!!\n",
	     p->msgtag, p->cur_sol.lp);
      break;
   }
   return(0);
}

