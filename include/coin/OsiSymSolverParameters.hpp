/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Callable         */
/* Library.                                                                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2004-2006 Ted Ralphs and Lehigh University.                 */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* The authors of this file are Menal Guzelsoy and Ted Ralphs                */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef OsiSymSolverParameters_hpp
#define OsiSymSolverParameters_hpp

enum OsiSymIntParam {
   /** This controls the level of output */
   OsiSymVerbosity,
   OsiSymNodeLimit,
   OsiSymFindFirstFeasible,
   OsiSymSearchStrategy,
   OsiSymUsePermanentCutPools,
   OsiSymKeepWarmStart,
   OsiSymDoReducedCostFixing,
   OsiSymMCFindSupportedSolutions,
   OsiSymSensitivityAnalysis,
   OsiSymRandomSeed,
   OsiSymDivingStrategy,
   OsiSymDivingK,
   OsiSymDivingThreshold,
   OsiSymTrimWarmTree,
   OsiSymGenerateCglGomoryCuts,
   OsiSymGenerateCglKnapsackCuts,
   OsiSymGenerateCglOddHoleCuts,
   OsiSymGenerateCglProbingCuts,
   OsiSymGenerateCglFlowAndCoverCuts,
   OsiSymGenerateCglRoundingCuts,
   OsiSymGenerateCglLiftAndProjectCuts,
   OsiSymGenerateCglCliqueCuts
};

enum OsiSymDblParam {
   /** The granularity is the actual minimum difference in objective function
       value for two solutions that actually have do different objective
       function values. For integer programs with integral objective function
       coefficients, this would be 1, for instance. */ 
   OsiSymGranularity,
   OsiSymTimeLimit,
   OsiSymGapLimit,
   OsiSymUpperBound,
   OsiSymLowerBound
};

enum OsiSymStrParam {
};

#endif
