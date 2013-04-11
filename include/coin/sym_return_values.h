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

#ifndef _RETURN_VALUES_H
#define _RETURN_VALUES_H

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************                  Return Values                       **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*----------------------- Global return codes -------------------------------*/
#define FUNCTION_TERMINATED_NORMALLY      0
#define FUNCTION_TERMINATED_ABNORMALLY   -1
#define FUNCTION_OUT_OF_MEMORY           -2
#define ERROR__USER                      -100

/*-------------- Return codes for sym_parse_comand_line() -------------------*/
#define ERROR__OPENING_PARAM_FILE        -110
#define ERROR__PARSING_PARAM_FILE        -111

/*----------------- Return codes for sym_load_problem() ---------------------*/
#define ERROR__READING_GMPL_FILE         -120
#define ERROR__READING_WARM_START_FILE   -121
#define ERROR__READING_MPS_FILE          -122
#define ERROR__READING_LP_FILE           -123

/*-------------------- Return codes for sym_solve() -------------------------*/
#define TM_NO_PROBLEM                     225
#define TM_NO_SOLUTION                    226
#define TM_OPTIMAL_SOLUTION_FOUND         227
#define TM_TIME_LIMIT_EXCEEDED            228
#define TM_NODE_LIMIT_EXCEEDED            229
#define TM_TARGET_GAP_ACHIEVED            230
#define TM_FOUND_FIRST_FEASIBLE           231
#define TM_FINISHED                       232
#define TM_UNFINISHED                     233
#define TM_FEASIBLE_SOLUTION_FOUND        234
#define TM_SIGNAL_CAUGHT                  235
#define PREP_OPTIMAL_SOLUTION_FOUND       236
#define PREP_NO_SOLUTION                  237
#define TM_ERROR__NO_BRANCHING_CANDIDATE -250
#define TM_ERROR__ILLEGAL_RETURN_CODE    -251
#define TM_ERROR__NUMERICAL_INSTABILITY  -252
#define TM_ERROR__COMM_ERROR             -253
#define TM_ERROR__USER                   -275
#define PREP_ERROR                       -276
#endif
