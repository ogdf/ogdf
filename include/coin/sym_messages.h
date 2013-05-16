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

#ifndef _MESSAGES_H
#define _MESSAGES_H

#define EMPTY_MSG_SEND(tid, msgtag) \
{ \
   int s_bufid; \
   if ((s_bufid = pvm_initsend (PvmDataRaw)) < 0) PVM_ERROR(); \
   if ((info = pvm_send((tid), (msgtag))) < 0) PVM_ERROR(); \
   if ((info = pvm_freebuf(s_bufid)) < 0) PVM_ERROR(); \
}

#define EMPTY_MSG_MCAST(tid, numtask, msgtag) \
{ \
   int s_bufid; \
   if ((s_bufid = pvm_initsend (PvmDataRaw)) < 0) PVM_ERROR(); \
   if ((info = pvm_mcast((tid), (numtask), (msgtag))) < 0) PVM_ERROR(); \
   if ((info = pvm_freebuf(s_bufid)) < 0) PVM_ERROR(); \
}

/*===========================================================================*
 * Message types
 * 1xx : general messages
 * 2xx : messages to/from master
 * 3xx : messages to/from tree manager
 * 4xx : rest of the messages
 *===========================================================================*/

/* we allow a process to die */
#define YOU_CAN_DIE 100

#define I_AM_DEAD 101

/*pretend to die*/
#define YOU_CANNOT_DIE 102

/* this is a new upper bound */
#define UPPER_BOUND 103

/* the tid of the master process is going to be sent */
#define MASTER_TID_INFO 104

/* tells the process to write a log file for a warm re-start in case of crash*/
#define WRITE_LOG_FILE 105


/* messages to set up an lp: lp-->master; master-->lp */
#define REQUEST_FOR_LP_DATA 200
#define LP_DATA             201
/* similar for cg */
#define REQUEST_FOR_CG_DATA 202
#define CG_DATA             203
/* similar for cp */
#define REQUEST_FOR_CP_DATA 204
#define CP_DATA             205
/* similar for sp */
#define REQUEST_FOR_SP_DATA 206
#define SP_DATA             207
/* similar for dg */
#define REQUEST_FOR_DG_DATA 208
#define DG_DATA             209

/* startup data for the TM */
#define TM_DATA             210

/* miscellaneous tm messages */
#define TM_ROOT_DESCRIPTION              211
#define TM_FIRST_PHASE_FINISHED          212

/*===========================================================================*
 * treemanager-->lp messages
 *===========================================================================*/
/* lp-->tm; Describes a particular search tree node */
#define LP__NODE_DESCRIPTION 300
/* lp-->tm; Describes the branching at this node */
#define LP__BRANCHING_INFO 301
/* lp-->tm; The LP is free to process a new node */
#define LP__IS_FREE 302
/* tm-->lp; 2nd phase started, from now on price*/
#define LP__SECOND_PHASE_STARTS 303

#define LP__CUT_NAMES_REQUESTED 304
#define LP__CUT_NAMES_SERVED  305

/* tm --> lp;  tm-->lp: this is your new active node, process it. */
#define LP__ACTIVE_NODE_DATA 306
/* tm-->lp: Instruction to the LP process whether to dive or not */
#define LP__DIVING_INFO 307

/* tm-->lp: The tid of the corresponding cut generator */
#define LP__CG_TID_INFO 308

/* lp-->tm: the newly sent active node is too expensive, hold it for the
            next phase */
#define LP__NODE_RESHELVED 309
/* lp-->tm: the newly sent active node is too expensive, discarded */
#define LP__NODE_DISCARDED 310
/* lp-->tm: timing data */
#define LP__TIMING 311


/*===========================================================================*
 * lp-->... messages
 *===========================================================================*/
/* lp-->master; different msgtypes for a feasible solution */
#define FEASIBLE_SOLUTION_NONZEROS  410
#define FEASIBLE_SOLUTION_FRACTIONS 411
#define FEASIBLE_SOLUTION_USER      412

/* lp-->cutgen,cutpool; a solution to be checked to find violated cuts
   and the solution is in a default format (userinds and values). Maybe
   nonzeros or fractions at this time */
#define LP_SOLUTION_NONZEROS  420
#define LP_SOLUTION_FRACTIONS 421
#define LP_SOLUTION_USER      422

/*===========================================================================*
 * treemanager-->cut/sol_pool messages
 *===========================================================================*/
/* notifies cut_pool to wait for a new set of cuts */
#define POOL_YOU_ARE_USELESS 501
#define POOL_USELESSNESS_ACKNOWLEDGED 502
/* tm asks a pool to copy itself into another pool */
#define POOL_COPY_YOURSELF 503
/* a pool reports back that it had finished and how much time it had used.
   a pool dies only if every node that belongs to it is fathomed */
#define POOL_TIME 504

/*****************************************************************************/

#define PACKED_CUT        600
#define PACKED_CUTS_TO_CP 601
#define CUTPOOL_COPY      602

#define NO_MORE_CUTS 605

/*****************************************************************************/

#define PACKED_COL 700
#define NO_MORE_COLS 701
#define CG_LP_SOLUTION 703

/*****************************************************************************/

#define SOMETHING_DIED 1000
#define NODE_DIED 1001


#endif
