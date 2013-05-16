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
/* The authors of this file are Menal Guzelsoy and Ted Ralphs.               */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#include "OsiSymSolverInterface.hpp"

#include <iostream>
using namespace std;

#include "CoinMpsIO.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

#include "symphony.h"

//#############################################################################
// A couple of helper functions
// Taken from OsiCpxSolverInterface.cpp.
//#############################################################################

inline void freeCacheDouble( double*& ptr )
{
  if( ptr != NULL )
    {
      delete [] ptr;
      ptr = NULL;
    }
}

inline void freeCacheChar( char*& ptr )
{
  if( ptr != NULL )
    {
      delete [] ptr;
      ptr = NULL;
    }
}

inline void freeCacheMatrix( CoinPackedMatrix*& ptr )
{
  if( ptr != NULL )
    {
      delete ptr;
      ptr = NULL;
    }
}

/* Default constructor */
/*===========================================================================*/
/*===========================================================================*/

OsiSymSolverInterface::OsiSymSolverInterface()
{

   env_ = sym_open_environment();

   gutsOfConstructor();
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem()
{
   void *user;

   sym_load_problem(env_);
   sym_get_user_data(env_, &user);
   setApplicationData((void *) (user));
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::branchAndBound()
{
   freeCachedResults();
   sym_solve(env_);

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::resolve()
{
   freeCachedResults();
   sym_warm_solve(env_);
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getLbForNewRhs(int cnt, int *index,
						  double * value)
{
   double newBound;
   if (!sym_get_lb_for_new_rhs(env_, cnt, index, value, &newBound)){
      return (newBound);
   } else {
      return (-sym_get_infinity());
   }
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getUbForNewRhs(int cnt, int *index,
					     double * value)
{
   double newBound;
   if (!sym_get_ub_for_new_rhs(env_, cnt, index, value, &newBound)){
      return (newBound);
   } else {
      return (sym_get_infinity());
   }
}

/*===========================================================================*/
/*===========================================================================*/

#if 0
double OsiSymSolverInterface::getLbForNewObj(int cnt, int *index,
					     double * value)
{
   double newBound;
   if (!sym_get_lb_for_new_obj(env_, cnt, index, value, &newBound)){
      return (newBound);
   } else {
      return (-sym_get_infinity());
   }
}
#endif

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getUbForNewObj(int cnt, int *index,
					     double * value)
{
   double newBound;
   if (!sym_get_ub_for_new_obj(env_, cnt, index, value, &newBound)){
      return (newBound);
   } else {
      return (sym_get_infinity());
   }
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::multiCriteriaBranchAndBound()
{

   sym_mc_solve(env_);

}

/*===========================================================================*/
/*===========================================================================*/

OsiSymSolverInterface::~OsiSymSolverInterface()
{

   sym_close_environment(env_);

   gutsOfDestructor();

   env_ = 0;

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::reset()
{

   sym_close_environment(env_);

   env_ = sym_open_environment();

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setIntParam(OsiIntParam key, int value)
{

   const char * keyVal;

   switch(key) {

    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
       keyVal = "node_limit";
       break;
    case OsiNameDiscipline:
       return false ;
    case OsiLastIntParam:
       return false;

    default:
       return false;
   }

   return(!sym_set_int_param(env_, keyVal, value)? true :false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(OsiSymIntParam key, int value)
{

   const char *keyVal;

   switch(key){

    case OsiSymVerbosity:
       keyVal = "verbosity";
       break;
    case OsiSymNodeLimit:
       keyVal = "node_limit";
       break;
    case OsiSymFindFirstFeasible:
       keyVal = "find_first_feasible";
       break;
    case OsiSymSearchStrategy:
       keyVal = "node_selection_rule";
       break;
    case OsiSymUsePermanentCutPools:
       keyVal = "use_permanent_cut_pools";
       break;
    case OsiSymKeepWarmStart:
       keyVal = "keep_warm_start";
       break;
    case OsiSymDoReducedCostFixing:
       keyVal = "do_reduced_cost_fixing";
       break;
    case OsiSymMCFindSupportedSolutions:
       keyVal = "mc_find_supported_solutions";
       break;
    case OsiSymSensitivityAnalysis:
       keyVal = "sensitivity_analysis";
       break;
    case OsiSymRandomSeed:
       keyVal = "random_seed";
       break;
    case OsiSymDivingStrategy:
       keyVal = "diving_strategy";
       break;
    case OsiSymDivingK:
       keyVal = "diving_k";
       break;
    case OsiSymDivingThreshold:
       keyVal = "diving_threshold";
       break;
    case OsiSymTrimWarmTree:
       keyVal = "trim_warm_tree";
       break;
    case OsiSymGenerateCglGomoryCuts:
       keyVal = "generate_cgl_gomory_cuts";
       break;
   case OsiSymGenerateCglKnapsackCuts:
       keyVal = "generate_cgl_knapsack_cuts";
       break;
   case OsiSymGenerateCglOddHoleCuts:
       keyVal = "generate_cgl_oddhole_cuts";
       break;
   case OsiSymGenerateCglProbingCuts:
       keyVal = "generate_cgl_probing_cuts";
       break;
   case OsiSymGenerateCglCliqueCuts:
       keyVal = "generate_cgl_clique_cuts";
       break;
   case OsiSymGenerateCglFlowAndCoverCuts:
       keyVal = "generate_cgl_flow_and_cover_cuts";
       break;
   case OsiSymGenerateCglRoundingCuts:
       keyVal = "generate_cgl_rounding_cuts";
       break;
   case OsiSymGenerateCglLiftAndProjectCuts:
       keyVal = "generate_cgl_lift_and_project_cuts";
       break;
    default:
       return false;
   }

   return(!sym_set_int_param(env_, keyVal, value)? true:false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(const std::string key, int value)
{

   return(!sym_set_int_param(env_, key.c_str(), value)?
	  true :false );

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setDblParam(OsiDblParam key, double value)
{

   switch(key){

    case OsiDualTolerance:
    case OsiPrimalTolerance:
       sym_set_dbl_param(env_, "granularity", value);
       sym_set_dbl_param(env_, "LP_granularity", value);
       return true;
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    case OsiLastDblParam:
       return false;
    case OsiObjOffset:
       sym_set_dbl_param(env_, "obj_offset", -value);
       return true;

    default:
       return false;
   }

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(OsiSymDblParam key, double value)
{

   const char *keyVal;

   switch(key){

    case OsiSymGranularity:
       keyVal = "granularity";
       break;
    case OsiSymTimeLimit:
       keyVal = "time_limit";
       break;
    case OsiSymGapLimit:
       keyVal = "gap_limit";
       break;
    case OsiSymUpperBound:
       keyVal = "upper_bound";
       break;
    case OsiSymLowerBound:
       keyVal = "lower_bound";
       break;
    default:
       return false;
   }

   return(!sym_set_dbl_param(env_, keyVal, value)? true: false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(const std::string key, double value)
{

   return(!sym_set_dbl_param(env_, key.c_str(), value)
	  ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setStrParam(OsiStrParam key,
					const std::string & value)
{

   const char * keyVal;

   switch(key) {

    case OsiProbName:
       keyVal = "problem_name";
       break;
    case OsiSolverName:
    case OsiLastStrParam:
       return false;
    default:
       return false;
   }

   return(!sym_set_str_param(env_, keyVal, const_cast<char *>(value.c_str()))
	  ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(OsiSymStrParam key,
					   const std::string & value)
{

   //switch(key){
   //// case ' ':
   //default:
      return false;
   //}

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setSymParam(const std::string key,
					const std::string value)
{

   return(!sym_set_str_param(env_, key.c_str(),
			     value.c_str()) ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getIntParam(OsiIntParam key, int& value) const
{

   const char *keyVal;

   switch(key) {

    case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
       keyVal = "node_limit";
       break;
    case OsiNameDiscipline:
       return false ;
    case OsiLastIntParam:
       return false;
    default:
       break;
   }

   return(!sym_get_int_param(env_, keyVal, &value) ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(OsiSymIntParam key, int& value) const
{

   const char *keyVal;

   switch(key){

    case OsiSymVerbosity:
       keyVal = "verbosity";
       break;
    case OsiSymNodeLimit:
       keyVal = "node_limit";
       break;
    case OsiSymFindFirstFeasible:
       keyVal = "find_first_feasible";
       break;
    case OsiSymSearchStrategy:
       keyVal = "node_selection_rule";
       break;
    case OsiSymUsePermanentCutPools:
       keyVal = "use_permanent_cut_pools";
       break;
    case OsiSymKeepWarmStart:
       keyVal = "keep_warm_start";
       break;
    case OsiSymDoReducedCostFixing:
       keyVal = "do_reduced_cost_fixing";
       break;
    case OsiSymMCFindSupportedSolutions:
       keyVal = "mc_find_supported_solutions";
       break;
    case OsiSymSensitivityAnalysis:
       keyVal = "sensitivity_analysis";
       break;
    case OsiSymRandomSeed:
       keyVal = "random_seed";
       break;
    case OsiSymDivingStrategy:
       keyVal = "diving_strategy";
       break;
    case OsiSymDivingK:
       keyVal = "diving_k";
       break;
    case OsiSymDivingThreshold:
       keyVal = "diving_threshold";
       break;
    case OsiSymTrimWarmTree:
       keyVal = "trim_warm_tree";
    default:
       return false;
   }

   return(!sym_get_int_param(env_, keyVal, &value) ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(const std::string key,
					int& value) const
{

   return (!sym_get_int_param(env_, key.c_str(), &value)
	   ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getDblParam(OsiDblParam key, double& value) const
{

   switch(key){

    case OsiDualTolerance:
    case OsiPrimalTolerance:
       sym_get_dbl_param(env_, "LP_granularity", &value);
       return true;
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    case OsiLastDblParam:
       return false;

    case OsiObjOffset:
       sym_get_dbl_param(env_, "obj_offset", &value);
       value = -value;
       return true;
    default:
      return false;
   }

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(OsiSymDblParam key,
					double& value) const
{

   const char * keyVal;

   switch(key){
    case OsiSymGranularity:
       keyVal = "granularity";
       break;
    case OsiSymTimeLimit:
       keyVal = "time_limit";
       break;
    case OsiSymGapLimit:
       keyVal = "gap_limit";
       break;
    case OsiSymUpperBound:
       keyVal = "upper_bound";
       break;
    case OsiSymLowerBound:
       keyVal = "lower_bound";
       break;
    default:
       return false;
   }

   return(!sym_get_dbl_param(env_, keyVal, &value) ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(const std::string key,
					double& value) const
{

   return (!sym_get_dbl_param(env_, key.c_str(), &value)
	   ? true : false);

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getStrParam(OsiStrParam key,
					std::string& value) const
{

   const char * keyVal;
   char * val;

   switch(key) {
    case OsiProbName:
       keyVal = "problem_name";
       break;
    case OsiSolverName:
       value = "sym";
       return true;
    case OsiLastStrParam:
       return false;
    default:
       return false;
   }

   if(!sym_get_str_param(env_, keyVal, &val)){
      value = val;
      return true;
   }

   return false;

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(OsiSymStrParam key,
					   std::string& value) const
{

   //switch(key){
   // // case ' ':
   // default:
       return false;
   //}

}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::getSymParam(const std::string key,
					std::string& value) const
{

   char * val;
   if (!sym_get_str_param(env_, key.c_str(), &val)){
      value = val;
      return true;
   }

   return false;

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setInitialData()
{

   sym_set_defaults(env_);

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::parseCommandLine(int argc, char **argv)
{

   sym_parse_command_line(env_, argc, argv);

}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::readMps(const char * infile, const char * extension)
{

   return(OsiSolverInterface::readMps(infile, extension));

#if 0
   int termcode = 0;
   char *fn = new char [MAX_FILE_NAME_LENGTH+1];

   sprintf(fn, "%s%s%s", infile, ".", extension);

   termcode = sym_read_mps(env_, fn);

   delete [] fn;

   return (termcode);
#endif
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::readGMPL(const char * modelFile,
				    const char * dataFile)
{
   int termcode = 0;
   termcode = sym_read_gmpl(env_, const_cast<char*>(modelFile),
			    const_cast<char*>(dataFile));

   return (termcode);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::findInitialBounds()
{

   sym_find_initial_bounds(env_);

}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::createPermanentCutPools()
{

   int cpNum;
   if (!sym_create_permanent_cut_pools(env_, &cpNum)){
      return (cpNum);
   }else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::initialSolve()
{

   freeCachedResults();
   sym_solve(env_);

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
					const double* collb,
					const double* colub, const double* obj,
					const double* rowlb,
					const double* rowub)
{
   const double inf = getInfinity();

   int nrows = matrix.getNumRows();
//   int ncols = matrix.getNumCols();
   char   * rowSense = new char  [nrows];
   double * rowRhs   = new double[nrows];
   double * rowRange = new double[nrows];

   int i;
   for ( i = nrows - 1; i >= 0; --i )
      {
	 const double lower = rowlb ? rowlb[i] : -inf;
	 const double upper = rowub ? rowub[i] : inf;
	 convertBoundToSense( lower, upper, rowSense[i], rowRhs[i],
			      rowRange[i] );
      }

   loadProblem( matrix, collb, colub, obj, rowSense, rowRhs, rowRange );

   delete [] rowSense;
   delete [] rowRhs;
   delete [] rowRange;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
					  double*& collb, double*& colub,
					  double*& obj, double*& rowlb,
					  double*& rowub)

{
   loadProblem( *matrix, collb, colub, obj, rowlb, rowub );
   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowlb;  rowlb = 0;
   delete[] rowub;  rowub = 0;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
					const double* collb,
					const double* colub,const double* obj,
					const char* rowsen,
					const double* rowrhs,
					const double* rowrng)
{

   CoinPackedMatrix * symMatrix;
   bool isColOrdered = true;
   int i, numelem = 0, * matbeg = NULL, * matind = NULL;
   double * matval = NULL;

   if (!matrix.isColOrdered()){
      symMatrix = new CoinPackedMatrix();
      symMatrix->copyOf(matrix);      symMatrix->reverseOrdering();
      isColOrdered = false;
   }else{
      symMatrix = const_cast<CoinPackedMatrix *>(&matrix);
   }

   int numcols = symMatrix->getNumCols();
   int numrows = symMatrix->getNumRows();


   if(numcols == 0 || numrows == 0){
      cout<<"loadProblem():The given matrix is empty!"<<endl;
      return;
   }

   const int * lengths = symMatrix->getVectorLengths();
   const int * matbegS = symMatrix->getVectorStarts();
   const int * matindS = symMatrix->getIndices();
   const double * matvalS = symMatrix->getElements();

   for (i = 0; i<numcols; i++){
      numelem += lengths[i];
   }

   if (numelem){

      matbeg = new int [numcols + 1];
      matind = new int[numelem];
      matval = new double [numelem];

      matbeg[0] = 0;
      for (i = 0; i<numcols; i++){
	 matbeg[i+1] = matbeg[i] + lengths[i];
	 if (lengths[i]){
	    memcpy(matind + matbeg[i], matindS + matbegS[i] ,
		   sizeof(int) * lengths[i]);
	    memcpy(matval + matbeg[i], matvalS + matbegS[i] ,
		   sizeof(double) * lengths[i]);
	 }
      }
   }

   bool rowsen_faked = false;
   char *fake_rowsen = 0;

   if (!rowsen){
      fake_rowsen  = new char[numrows];
      memset(fake_rowsen, 'G', CSIZE*numrows);
      rowsen_faked = true;
   }

   loadProblem(numcols,numrows, matbeg, matind, matval, collb, colub, obj,
	       rowsen ? rowsen:fake_rowsen, rowrhs, rowrng);

   if (rowsen_faked)
      delete[] fake_rowsen;

   if(!isColOrdered)
      delete symMatrix;

   if(numelem){
      delete [] matbeg;
      delete [] matind;
      delete [] matval;
   }
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
					  double*& collb, double*& colub,
					  double*& obj, char*& rowsen,
					  double*& rowrhs, double*& rowrng)
{
   loadProblem( *matrix, collb, colub, obj, rowsen, rowrhs, rowrng );

   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowsen; rowsen = 0;
   delete[] rowrhs; rowrhs = 0;
   delete[] rowrng; rowrng = 0;
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const int numcols, const int numrows,
					const CoinBigIndex * start,
					const int* index, const double* value,
					const double* collb,
					const double* colub, const double* obj,
					const double* rowlb,
					const double* rowub)
{

   if(numcols == 0 || numrows == 0){
      cout<<"loadProblem():The given problem is empty!"<<endl;
      return;
   }

   const double inf = getInfinity();

   char   * sense = new char  [numrows];
   double * rhs   = new double[numrows];
   double * range = new double[numrows];

   int i;
   for ( i = numrows - 1; i >= 0; --i ){
      const double lower = rowlb ? rowlb[i] : -inf;
      const double upper = rowub ? rowub[i] : inf;

      /* convertBountToSense */
      range[i] = 0.0;
      if (lower > -inf) {
	 if (upper < inf) {
	    rhs[i] = upper;
	    if (upper==lower) {
	       sense[i] = 'E';
	    } else {
	       sense[i] = 'R';
	       range[i] = upper - lower;
	    }
	 } else {
	    sense[i] = 'G';
	    rhs[i] = lower;
	 }
      } else {
	 if (upper < inf) {
	    sense[i] = 'L';
	    rhs[i] = upper;
	 } else {
	    sense[i] = 'N';
	    rhs[i] = 0.0;
	 }
      }
      /*	 convertBoundToSense( lower, upper, rowSense[i], rowRhs[i],
		 rowRange[i] ); */
   }

   loadProblem(numcols,numrows, start, index, value, collb, colub, obj,
	       sense, rhs, range);

   delete [] sense;
   delete [] rhs;
   delete [] range;

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::loadProblem(const int numcols, const int numrows,
					const CoinBigIndex * start,
					const int* index, const double* value,
					const double* collb,
					const double* colub, const double* obj,
					const char* rowsen,
					const double* rowrhs,
					const double* rowrng)
{
   void *user;
   freeAllMemory();

   sym_explicit_load_problem(env_, numcols, numrows, const_cast<int*>(start),
			     const_cast<int*>(index),
			     const_cast<double*>(value),
			     const_cast<double*>(collb),
			     const_cast<double*>(colub), NULL,
			     const_cast<double*>(obj), NULL,
			     const_cast<char*>(rowsen),
			     const_cast<double*>(rowrhs),
			     const_cast<double*>(rowrng), true);


   sym_get_user_data(env_, &user);
   setApplicationData((void *) (user));
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isAbandoned() const
{

   if(sym_is_abandoned(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isProvenOptimal() const
{
   if(sym_is_proven_optimal(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isProvenPrimalInfeasible() const
{
   if(sym_is_proven_primal_infeasible(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

//bool OsiSymSolverInterface::isPrimalObjectiveLimitReached() const
//{
//   if(sym_is_target_gap_achieved(env_)){
//      return true;
//   }
//   else{
//      return false;
//   }
//}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isIterationLimitReached() const
{
   if(sym_is_iteration_limit_reached(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isTimeLimitReached() const
{
   if(sym_is_time_limit_reached(env_)){
      return true;
   }
   else{
      return false;
   }
}
/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isTargetGapReached() const
{
   if(sym_is_target_gap_achieved(env_)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface:: getNumCols() const
{
   int numCols;
   if (!sym_get_num_cols(env_, &numCols)){
      return (numCols);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getNumRows() const
{
   int numRows;
   if(!sym_get_num_rows(env_, &numRows)){
      return(numRows);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getNumElements() const
{
   int numElems;
   if(!sym_get_num_elements(env_, &numElems)){
      return(numElems);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColLower() const
{

   if(!collower_){
      collower_ = new double[getNumCols()];
   }

   if(!sym_get_col_lower(env_, collower_)){
      return (collower_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColUpper() const
{
   if(!colupper_){
      colupper_ = new double[getNumCols()];
   }

   if(!sym_get_col_upper(env_, colupper_)){
      return (colupper_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const char * OsiSymSolverInterface::getRowSense() const
{
   if(!rowsense_){
      rowsense_ = new char[getNumRows()];
   }

   if(!sym_get_row_sense(env_, rowsense_)){
      return (rowsense_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRightHandSide() const
{
   if(!rhs_){
      rhs_ = new double[getNumRows()];
   }

   if(!sym_get_rhs(env_, rhs_)){
      return (rhs_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowRange() const
{
   if(!rowrange_){
      rowrange_ = new double[getNumRows()];
   }

   if(!sym_get_row_range(env_, rowrange_)){
      return (rowrange_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowLower() const
{
   if(!rowlower_){
      rowlower_ = new double[getNumRows()];
   }

   if(!sym_get_row_lower(env_, rowlower_)){
      return (rowlower_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowUpper() const
{
   if(!rowupper_){
      rowupper_ = new double[getNumRows()];
   }

   if(!sym_get_row_upper(env_, rowupper_)){
      return (rowupper_);
   } else {
      return (0);
   }
}


/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowPrice() const
{
   if(!rowprice_){
      rowprice_ = new double[getNumRows()];
      memset(rowprice_, 0, getNumRows() * sizeof(double));
   }

   return rowprice_;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getReducedCost() const
{
   if(!colredcost_){
      colredcost_ = new double[getNumCols()];
      memset(colredcost_, 0, getNumCols() * sizeof(double));
   }

   return colredcost_;
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getObjCoefficients() const
{

   if(!obj_){
      obj_ = new double[getNumCols()];
   }

   if(!sym_get_obj_coeff(env_, obj_)){
      return (obj_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getObj2Coefficients() const
{

   if(!obj2_){
      obj2_ = new double[getNumCols()];
   }

   if(!sym_get_obj2_coeff(env_, obj2_)){
      return (obj2_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getObjSense() const
{
   int objSen;
   if(!sym_get_obj_sense(env_, &objSen)){
      return ((double)objSen);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isContinuous(int colIndex) const
{
   int value;
   if(!sym_is_continuous(env_, colIndex, &value)){
      return (value != 0);
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isBinary(int colIndex) const
{
   int value;
   if(!sym_is_binary(env_, colIndex, &value)){
      return (value != 0);
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isInteger(int colIndex) const
{
   char value;
   if(!sym_is_integer(env_, colIndex, &value)){
      return (value != 0);
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isIntegerNonBinary(int colIndex) const
{
   if(!isBinary(colIndex) && isInteger(colIndex)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::isFreeBinary(int colIndex) const
{

   if(isBinary(colIndex)){
      return true;
   }
   else{
      return false;
   }
}

/*===========================================================================*/
/*===========================================================================*/

const CoinPackedMatrix * OsiSymSolverInterface::getMatrixByRow() const
{
   if(!matrixByRow_){
      matrixByRow_ = new CoinPackedMatrix(*getMatrixByCol());
   }else{
      matrixByRow_->copyOf(*getMatrixByCol());
   }

   matrixByRow_->reverseOrdering();

   return matrixByRow_;

}
/*===========================================================================*/
/*===========================================================================*/

const CoinPackedMatrix * OsiSymSolverInterface::getMatrixByCol() const
{

   int m, n, nz, *matind, *matbeg;
   double *matval;


   m = getNumRows();
   n = getNumCols();
   nz = getNumElements();

   matbeg = new int[n + 1];
   matind = new int[nz];
   matval = new double[nz];

   sym_get_matrix(env_, &nz, matbeg, matind, matval);

   if(!matrixByCol_){
      matrixByCol_ =
	 new CoinPackedMatrix(true, m, n, nz,
			      matval, matind,
			      matbeg, 0);
   }else{
      matrixByCol_->copyOf(true, m, n, nz,
			   matval, matind,
			   matbeg, 0);
   }

   delete [] matbeg;
   delete [] matind;
   delete [] matval;

   return matrixByCol_;
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getInfinity() const
{
   return sym_get_infinity();
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getColSolution() const
{

   int n = getNumCols();

   if(!colsol_){
      colsol_ = new double[n];
   }
   if (sym_get_col_solution(env_, colsol_)){
      if (!getNumCols()){
	 return (0);
      }
   }

   return (colsol_);
}

/*===========================================================================*/
/*===========================================================================*/

const double * OsiSymSolverInterface::getRowActivity() const
{

   if(!rowact_){
      rowact_ = new double[getNumRows()];
   }

   if(!sym_get_row_activity(env_, rowact_)){
      return (rowact_);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getObjValue() const
{
   double objVal;
   if(sym_get_obj_val(env_, &objVal)){
      if(!getNumCols()){
	 return (0);
      }
   }

   return(objVal);
}

/*===========================================================================*/
/*===========================================================================*/

double OsiSymSolverInterface::getPrimalBound() const
{
   double ubPri;
   if(!sym_get_primal_bound(env_, &ubPri)){
      return (ubPri);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

int OsiSymSolverInterface::getIterationCount() const
{
   int numNodes;
   if(!sym_get_iteration_count(env_, &numNodes)){
      return (numNodes);
   } else {
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObjCoeff( int elementIndex,
					 double elementValue )
{
   freeCachedData(KEEPCACHED_ROW);
   sym_set_obj_coeff(env_, elementIndex, elementValue);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObj2Coeff( int elementIndex,
					 double elementValue )
{
   sym_set_obj2_coeff(env_, elementIndex, elementValue);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColLower( int elementIndex,
					 double elementValue )
{
   freeCachedData(KEEPCACHED_ROW);
   sym_set_col_lower(env_, elementIndex, elementValue);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColUpper( int elementIndex,
					 double elementValue )
{
   freeCachedData(KEEPCACHED_ROW);
   sym_set_col_upper(env_, elementIndex, elementValue);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowLower( int elementIndex,
					 double elementValue )
{
   freeCachedData(KEEPCACHED_COLUMN);
   sym_set_row_lower(env_, elementIndex, elementValue);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowUpper( int elementIndex,
					 double elementValue )
{
   freeCachedData(KEEPCACHED_COLUMN);
   sym_set_row_upper(env_, elementIndex, elementValue);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowType(int index, char sense,
				       double rightHandSide, double range)
{
   freeCachedData(KEEPCACHED_COLUMN);
   sym_set_row_type(env_, index, sense, rightHandSide, range);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setObjSense(double s)
{
   sym_set_obj_sense(env_, (int)s);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColSolution(const double *colsol)
{
   freeCachedResults();
   sym_set_col_solution(env_, const_cast<double *>(colsol));
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setPrimalBound(const double bound)
{
   sym_set_primal_bound(env_, bound);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setRowPrice(const double * rowprice)
{
   if( rowprice_ == NULL )
   {
      rowprice_ = new double[getNumRows()];
   }
   memcpy(rowprice_, rowprice, getNumRows() * sizeof(double));
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setContinuous(int index)
{
   sym_set_continuous(env_, index);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setInteger(int index)
{
   sym_set_integer(env_, index);
}
/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::setColName(char **colname)
{
   sym_set_col_names (env_, colname);
}

/*===========================================================================*/
/*===========================================================================*/
void OsiSymSolverInterface::addCol(const CoinPackedVectorBase& vec,
				   const double collb, const double colub,
				   const double obj)
{

   int numElements, *indices = 0;
   double *elements = 0;

   freeCachedData(KEEPCACHED_ROW);

   if((numElements = vec.getNumElements())){
      indices = const_cast<int*>(vec.getIndices());
      elements = const_cast<double*>(vec.getElements());
   }

   sym_add_col(env_, numElements, indices, elements, collb, colub, obj, false,
	       NULL);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::addRow(const CoinPackedVectorBase& vec,
				   const double rowlb, const double rowub)
{
   char    rowSen;
   double  rowRhs;
   double  rowRange;

   convertBoundToSense( rowlb, rowub, rowSen, rowRhs, rowRange);
   addRow(vec, rowSen, rowRhs, rowRange);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::addRow(const CoinPackedVectorBase& vec,
				   const char rowsen, const double rowrhs,
				   const double rowrng)
{
   int numElements, *indices = 0;
   double *elements = 0;

   freeCachedData(KEEPCACHED_COLUMN);

   if((numElements = vec.getNumElements())){
      indices = const_cast<int*>(vec.getIndices());
      elements = const_cast<double*>(vec.getElements());
   }

   sym_add_row(env_, numElements, indices, elements, rowsen, rowrhs, rowrng);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::deleteCols(const int num, const int * colIndices)
{
   freeCachedData(KEEPCACHED_ROW);

   sym_delete_cols(env_, num, const_cast<int*>(colIndices));
}

/*===========================================================================*/
/*===========================================================================*/
void OsiSymSolverInterface::deleteRows(const int num, const int * rowIndices)
{
   freeCachedData(KEEPCACHED_COLUMN);

   sym_delete_rows(env_, num, const_cast<int*>(rowIndices));
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::writeMps(const char *filename,
				     const char *extension,
				     double objSense) const
{

   int i, n;
   char ** colnames = 0;
   char ** rownames = 0;
   const CoinPackedMatrix * colMat = getMatrixByCol();
   char * isInteger;

   n = getNumCols();
   isInteger = new char[n];

   for( i = 0; i < n; i++){
      sym_is_integer(env_, i, &isInteger[i]);
   }

   CoinMpsIO mps;
   mps.setMpsData(*colMat, getInfinity(), getColLower(),
		  getColUpper(), getObjCoefficients(), isInteger,
		  getRowSense(), getRightHandSide(),
		  getRowRange(), colnames, rownames);

   string f(filename);
   string e(extension);
   string fullname = f + "." + e;

   mps.writeMps(fullname.c_str());

   delete [] isInteger;

}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::applyRowCut( const OsiRowCut & rc)
{
   double lb, ub;
   CoinPackedVector rowVec;

   freeCachedData(KEEPCACHED_COLUMN);

   rowVec = rc.row();
   lb = rc.lb();
   ub = rc.ub();

   addRow(rowVec, lb, ub);
}

/*===========================================================================*/
/*===========================================================================*/

void OsiSymSolverInterface::applyColCut( const OsiColCut & cc)
{

   int i;

   const CoinPackedVector & lbs = cc.lbs();
   const CoinPackedVector & ubs = cc.ubs();

   const int * colInd = lbs.getIndices();
   const double * colVal = lbs.getElements();

   freeCachedData(KEEPCACHED_ROW);

   for(i = 0; i<lbs.getNumElements(); i++){
      sym_set_col_lower(env_, colInd[i], colVal[i]);
   }

   colInd = ubs.getIndices();
   colVal = ubs.getElements();

   for(i = 0; i<ubs.getNumElements(); i++){
      sym_set_col_upper(env_, colInd[i], colVal[i]);
   }
}

/*===========================================================================*/
/*===========================================================================*/

CoinWarmStart * OsiSymSolverInterface::getWarmStart() const
{

   /* FIXME! Add a SymIntParam to determine whether keep the tree
      in case a getWarmStart() is called or not!
   */

   warm_start_desc * ws = sym_get_warm_start(env_, true);
   if (ws){
      SymWarmStart * symWS = new SymWarmStart(ws);
      sym_delete_warm_start(ws);
      return symWS;
   } else {
      sym_delete_warm_start(ws);
      return (0);
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool OsiSymSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{

   freeCachedResults();
   CoinWarmStart * wsC = const_cast<CoinWarmStart*> (warmstart);
   SymWarmStart *symWS = dynamic_cast<SymWarmStart*>(wsC);

   if (symWS==NULL) {
   	 cout << "setWarmStart(): No SymWarmStart was given!" << endl;
   	 return false;
   }

   warm_start_desc * ws  = symWS->getCopyOfWarmStartDesc();

   if(!ws){
      cout<<"setWarmStart(): An empty warmstart was given!"<<endl;
      return false;
   }

   sym_set_warm_start(env_, ws);
   sym_delete_warm_start(ws);

   return true;
}

/*===========================================================================*/
/* copy constructor, clone and assignment operator                           */
/*===========================================================================*/

OsiSymSolverInterface::OsiSymSolverInterface(const OsiSymSolverInterface & si)
{
   env_= sym_create_copy_environment(si.getSymphonyEnvironment());

   gutsOfConstructor();

   /* Note that, if the user structure was set by
      OsiSolverInterface::setApplicationData(), since SYMPHONY
      will know nothing about the user structure, it will not be possible
      to copy that!  For a temporary solution, the OsiSolverInterface::appData_
      will be directed to the original user structure! So,
      be careful from now on that, further modifications
      on the user structure of either the original or the clone OsiSym
      objects will affect the both! */
   setApplicationData(si.getApplicationData());
}

/*===========================================================================*/
/*===========================================================================*/

OsiSolverInterface * OsiSymSolverInterface::clone(bool copyData) const
{
   return (new OsiSymSolverInterface(*this));
}

/*===========================================================================*/
/*===========================================================================*/

OsiSymSolverInterface & OsiSymSolverInterface::operator=(const OsiSymSolverInterface& rhs)
{
   //   OsiSolverInterface * si_copy = const_cast<OsiSolverInterface *>(&rhs);
   //   OsiSymSolverInterface * sym = dynamic_cast<OsiSymSolverInterface*>(si_copy);

   if(this != &rhs){

      sym_close_environment(env_);
      gutsOfDestructor();
      env_= sym_create_copy_environment(rhs.getSymphonyEnvironment());
      gutsOfConstructor();

      /* Note that, if the user structure was set by
	 OsiSolverInterface::setApplicationData(), since SYMPHONY
	 will know nothing about the user structure, it will not be possible
	 to copy that!  For a temporary solution, the
	 OsiSolverInterface::appData_ will be directed to the original user
	 structure! So, be careful from now on that, further modifications
	 on the user structure of either the original or the clone OsiSym
	 objects will affect the both! */

      setApplicationData(rhs.getApplicationData());
   }

   return *this;
}

void OsiSymSolverInterface::gutsOfConstructor()
{
	obj_ = NULL;
	collower_ = NULL;
	colupper_ = NULL;
	colredcost_ = NULL;
	rowsense_ = NULL;
	rhs_ = NULL;
	rowrange_ = NULL;
	rowlower_ = NULL;
	rowupper_ = NULL;
	rowprice_ = NULL;
	colsol_ = NULL;
	rowact_ = NULL;
	matrixByRow_ = NULL;
	matrixByCol_ = NULL;
}

void OsiSymSolverInterface::gutsOfDestructor()
{
        freeAllMemory();

	assert( obj_ == NULL );
	assert( collower_ == NULL );
	assert( colupper_ == NULL );
  assert( colredcost_ == NULL );
	assert( rowsense_ == NULL );
	assert( rhs_ == NULL );
	assert( rowrange_ == NULL );
	assert( rowlower_ == NULL );
	assert( rowupper_ == NULL );
	assert( rowprice_ == NULL );
	assert( colsol_ == NULL );
	assert( rowact_ == NULL );
	assert( matrixByRow_ == NULL );
	assert( matrixByCol_ == NULL );
}

//-----------------------------------------------------------------------------
// free cached vectors
//-----------------------------------------------------------------------------

void OsiSymSolverInterface::freeCachedColRim()
{
   freeCacheDouble(obj_);
   freeCacheDouble(collower_);
   freeCacheDouble(colupper_);
   freeCacheDouble(colredcost_);
}

//-----------------------------------------------------------------------------

void OsiSymSolverInterface::freeCachedRowRim()
{
   freeCacheChar(rowsense_);
   freeCacheDouble(rhs_);
   freeCacheDouble(rowrange_);
   freeCacheDouble(rowlower_);
   freeCacheDouble(rowupper_);
   freeCacheDouble(rowprice_);
}

//-----------------------------------------------------------------------------

void OsiSymSolverInterface::freeCachedMatrix()
{
   freeCacheMatrix(matrixByRow_);
   freeCacheMatrix(matrixByCol_);
}

//-----------------------------------------------------------------------------

void OsiSymSolverInterface::freeCachedResults()
{
   freeCacheDouble(colsol_);
   freeCacheDouble(rowact_);
}

//-----------------------------------------------------------------------------

void OsiSymSolverInterface::freeCachedData( int keepCached )
{
   if( !(keepCached & OsiSymSolverInterface::KEEPCACHED_COLUMN) )
      freeCachedColRim();
   if( !(keepCached & OsiSymSolverInterface::KEEPCACHED_ROW) )
      freeCachedRowRim();
   if( !(keepCached & OsiSymSolverInterface::KEEPCACHED_MATRIX) )
      freeCachedMatrix();
   if( !(keepCached & OsiSymSolverInterface::KEEPCACHED_RESULTS) )
      freeCachedResults();
}

//-----------------------------------------------------------------------------

void OsiSymSolverInterface::freeAllMemory()
{
   freeCachedData();
}

/*===========================================================================*/
/*===========================================================================*/
