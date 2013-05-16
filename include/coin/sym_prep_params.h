/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* The author of this file is Menal Guzelsoy                                 */
/*                                                                           */
/* (c) Copyright 2006-2011 Lehigh University. All Rights Reserved.           */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _PREP_PARAMS_H
#define _PREP_PARAMS_H

/*---------------------------------------------------------------------------*\
| The list of parameters associated with pre-processing                       |
|                                                                             |
\*---------------------------------------------------------------------------*/
typedef struct PREP_PARAMS{
   int               level;
   int               dive_level;
   int               impl_dive_level;
   int               impl_limit;
   int               do_probe;
   int               verbosity;
   int               reduce_mip;
   int               probe_verbosity;
   int               probe_level;
   int               display_stats;
   double            etol;
   int               keep_row_ordered;
   int               do_single_row_rlx;
   double            single_row_rlx_ratio;
   int               max_sr_cnt;
   char              do_aggregate_row_rlx;
   double            max_aggr_row_ratio;
   int               max_aggr_row_cnt;
   int               iteration_limit;
   int               keep_track;
   int               write_mps;
   int               write_lp;
   int               time_limit;
}prep_params;

#endif
