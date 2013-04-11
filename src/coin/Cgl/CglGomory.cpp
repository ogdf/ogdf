// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
//#define CGL_DEBUG 1
//#ifdef NDEBUG
//#undef NDEBUG
//#endif
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#define COIN_HAS_CLP_GOMORY
#ifdef COIN_HAS_CLP_GOMORY
#include "OsiClpSolverInterface.hpp"
#endif
#include "CoinFactorization.hpp"
#undef CLP_OSL
#if 1
#define CLP_OSL 1
#if CLP_OSL!=1&&CLP_OSL!=3
#undef CLP_OSL
#else
#include "CoinOslFactorization.hpp"
#endif
#endif
#include "CoinWarmStartBasis.hpp"
#include "CglGomory.hpp"
#include "CoinFinite.hpp"
#ifdef CGL_DEBUG_GOMORY
int gomory_try=CGL_DEBUG_GOMORY;
#endif
//-------------------------------------------------------------------
// Generate Gomory cuts
//------------------------------------------------------------------- 
void CglGomory::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info) const
{
#ifdef CGL_DEBUG_GOMORY
  gomory_try++;
#endif
  // Get basic problem information
  int numberColumns=si.getNumCols(); 
  
  // get integer variables and basis
  char * intVar = new char[numberColumns];
  int i;
  CoinWarmStart * warmstart = si.getWarmStart();
  CoinWarmStartBasis* warm =
    dynamic_cast<CoinWarmStartBasis*>(warmstart);
  const double * colUpper = si.getColUpper();
  const double * colLower = si.getColLower();
  //#define CLP_INVESTIGATE2
#ifndef CLP_INVESTIGATE2
  if ((info.options&16)!=0)
#endif
    printf("%d %d %d\n",info.inTree,info.options,info.pass);
  for (i=0;i<numberColumns;i++) {
    if (si.isInteger(i)) {
      if (colUpper[i]>colLower[i]+0.5) {
	if (fabs(colUpper[i]-1.0)<1.0e-12&&fabs(colLower[i])<1.0e-12) {
	  intVar[i]=1; //0-1
	} else  if (colLower[i]>=0.0) {
	  intVar[i] = 2; // other
	} else {
	  // negative bounds - I am not sure works
	  intVar[i] = 0;
	}
      } else {
	intVar[i] = 0;
      }
    } else {
      intVar[i]=0;
    }
  }
  const OsiSolverInterface * useSolver=&si;
#ifdef COIN_HAS_CLP_GOMORY
  double * objective = NULL;
  OsiClpSolverInterface * clpSolver = dynamic_cast<OsiClpSolverInterface *>(originalSolver_);
  int numberOriginalRows;
  if (clpSolver) {
    useSolver = originalSolver_;
    assert (gomoryType_);
    // check simplex is plausible
    if (!clpSolver->getNumRows()||numberColumns!=clpSolver->getNumCols()) {
      delete originalSolver_;
      originalSolver_=si.clone();
      clpSolver = dynamic_cast<OsiClpSolverInterface *>(originalSolver_);
      assert (clpSolver);
      useSolver = originalSolver_;
    }
    ClpSimplex * simplex = clpSolver->getModelPtr();
    numberOriginalRows = simplex->numberRows();
    int numberRows = si.getNumRows();
    assert (numberOriginalRows<=numberRows);
    // only do if different (unless type 2x)
    int gomoryType = gomoryType_%10;
    int whenToDo = gomoryType_/10;
    if (whenToDo==2 ||(numberRows>numberOriginalRows && whenToDo==1
		       && (info.options&512)==0) ||
	((info.options&1024)!=0 && (info.options&512)==0
	 && numberTimesStalled_<3)) {
      // bounds
      memcpy(simplex->columnLower(),colLower,numberColumns*sizeof(double));
      memcpy(simplex->columnUpper(),colUpper,numberColumns*sizeof(double));
      double * obj = simplex->objective();
      objective = CoinCopyOfArray(obj,numberColumns);
      const double * pi = si.getRowPrice();
      const CoinPackedMatrix * rowCopy = si.getMatrixByRow();
      const int * column = rowCopy->getIndices();
      const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
      const int * rowLength = rowCopy->getVectorLengths(); 
      const double * rowElements = rowCopy->getElements();
      const double * rowLower = si.getRowLower();
      const double * rowUpper = si.getRowUpper();
      int numberCopy;
      int numberAdd;
      double * rowLower2 = NULL;
      double * rowUpper2 = NULL;
      int * column2 = NULL;
      CoinBigIndex * rowStart2 = NULL;
      double * rowElements2 = NULL;
      char * copy = new char [numberRows-numberOriginalRows];
      memset(copy,0,numberRows-numberOriginalRows);
      if (gomoryType==2) {
	numberCopy=0;
	numberAdd=0;
	for (int iRow=numberOriginalRows;iRow<numberRows;iRow++) {
	  bool simple = true;
	  for (int k=rowStart[iRow];
	       k<rowStart[iRow]+rowLength[iRow];k++) {
	    double value = rowElements[k];
	    if (value!=floor(value+0.5)) {
	      simple=false;
	      break;
	    }
	  }
	  if (simple) {
	    numberCopy++;
	    numberAdd+=rowLength[iRow];
	    copy[iRow-numberOriginalRows]=1;
	  }
	}
	if (numberCopy) {
	  //printf("Using %d rows out of %d\n",numberCopy,numberRows-numberOriginalRows);
	  rowLower2 = new double [numberCopy];
	  rowUpper2 = new double [numberCopy];
	  rowStart2 = new CoinBigIndex [numberCopy+1];
	  rowStart2[0]=0;
	  column2 = new int [numberAdd];
	  rowElements2 = new double [numberAdd];
	}
      }
      numberCopy=0;
      numberAdd=0;
      for (int iRow=numberOriginalRows;iRow<numberRows;iRow++) {
	if (!copy[iRow-numberOriginalRows]) {
	  double value = pi[iRow];
	  for (int k=rowStart[iRow];
	       k<rowStart[iRow]+rowLength[iRow];k++) {
	    int iColumn=column[k];
	    obj[iColumn] -= value*rowElements[k];
	  }
	} else {
	  rowLower2[numberCopy]=rowLower[iRow];
	  rowUpper2[numberCopy]=rowUpper[iRow];
	  for (int k=rowStart[iRow];
	       k<rowStart[iRow]+rowLength[iRow];k++) {
	    column2[numberAdd]=column[k];
	    rowElements2[numberAdd++]=rowElements[k];
	  }
	  numberCopy++;
	  rowStart2[numberCopy]=numberAdd;
	}
      }
#if 0
      CoinThreadRandom randomNumberGenerator;
      const double * solution = si.getColSolution();
      for (int i=0;i<numberColumns;i++) {
	if (intVar[i]==1) {
	  double randomNumber = randomNumberGenerator.randomDouble();
	  //randomNumber = 0.001*floor(randomNumber*1000.0);
	  if (solution[i]>0.5)
	    obj[i] -= randomNumber*0.001*fabs(obj[i]);
	  else
	    obj[i] += randomNumber*0.001*fabs(obj[i]);
	}
      }
#endif
      if (numberCopy) {
	clpSolver->addRows(numberCopy,
			   rowStart2,column2,rowElements2,
			   rowLower2,rowUpper2);
	delete [] rowLower2 ;
	delete [] rowUpper2 ;
	delete [] column2 ;
	delete [] rowStart2 ;
	delete [] rowElements2 ;
      }
      delete [] copy;
      memcpy(simplex->primalColumnSolution(),si.getColSolution(),
	     numberColumns*sizeof(double));
      warm->resize(numberOriginalRows,numberColumns);
      clpSolver->setBasis(*warm);
      delete warm;
      simplex->setDualObjectiveLimit(COIN_DBL_MAX);
      simplex->setLogLevel(0);
      simplex->primal(1);
      //printf("Trying - %d its\n",simplex->numberIterations());
      //simplex->setLogLevel(0);
      warm=simplex->getBasis();
      warmstart=warm;
      assert (!simplex->status());
    } else {
      // don't do
      delete warmstart;
      warmstart=NULL;
      if ((info.options&1024)==0)
	numberTimesStalled_=0;
    }
  }
#endif
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&!debugger->onOptimalPath(si))
    debugger = NULL;
#else
  const OsiRowCutDebugger * debugger = NULL;
#endif
  int numberRowCutsBefore = cs.sizeRowCuts();

  if (warmstart)
    generateCuts(debugger, cs, *useSolver->getMatrixByCol(), 
		 *useSolver->getMatrixByRow(),
		 useSolver->getColSolution(),
		 useSolver->getColLower(), useSolver->getColUpper(), 
		 useSolver->getRowLower(), useSolver->getRowUpper(),
		 intVar,warm,info);
#ifdef COIN_HAS_CLP_GOMORY
  if (objective) {
    ClpSimplex * simplex = clpSolver->getModelPtr();
    memcpy(simplex->objective(),objective,numberColumns*sizeof(double));
    delete [] objective;
    // take out locally useless cuts
    const double * solution = si.getColSolution();
    double primalTolerance = 1.0e-7;
    int numberRowCutsAfter = cs.sizeRowCuts();
    for (int k = numberRowCutsAfter - 1; k >= numberRowCutsBefore; k--) {
      const OsiRowCut * thisCut = cs.rowCutPtr(k) ;
      double sum = 0.0;
      int n = thisCut->row().getNumElements();
      const int * column = thisCut->row().getIndices();
      const double * element = thisCut->row().getElements();
      assert (n);
      for (int i = 0; i < n; i++) {
	double value = element[i];
	sum += value * solution[column[i]];
      }
      if (sum > thisCut->ub() + primalTolerance) {
	sum = sum - thisCut->ub();
      } else if (sum < thisCut->lb() - primalTolerance) {
	sum = thisCut->lb() - sum;
      } else {
	sum = 0.0;
      }
      if (!sum) {
	// take out
	cs.eraseRowCut(k);
      }
    }
#ifdef CLP_INVESTIGATE2
    printf("OR %p pass %d inTree %c - %d cuts (but %d deleted)\n",
       originalSolver_,info.pass,info.inTree?'Y':'N',
       numberRowCutsAfter-numberRowCutsBefore,
       numberRowCutsAfter-cs.sizeRowCuts());
#endif
  }
#endif

  delete warmstart;
  delete [] intVar;
  if ((!info.inTree&&((info.options&4)==4||((info.options&8)&&!info.pass)))
      ||(info.options&16)!=0) {
    int limit = maximumLengthOfCutInTree();
    int numberRowCutsAfter = cs.sizeRowCuts();
    for (int i=numberRowCutsBefore;i<numberRowCutsAfter;i++) {
      int length = cs.rowCutPtr(i)->row().getNumElements();
      if (length<=limit)
	cs.rowCutPtr(i)->setGloballyValid();
    }
  }
  if ((gomoryType_%10)==2) {
    // back to original
    int numberRows = clpSolver->getNumRows();
    if (numberRows>numberOriginalRows) {
      int numberDelete = numberRows-numberOriginalRows;
      int * delRow = new int [numberDelete];
      for (int i=0;i<numberDelete;i++)
	delRow[i]=i+numberOriginalRows;
      clpSolver->deleteRows(numberDelete,delRow);
      delete [] delRow;
    }
  }
}

// Returns value - floor but allowing for small errors
inline double above_integer(double value) {
  double value2=floor(value);
  double value3=floor(value+0.5);
  if (fabs(value3-value)<1.0e-7*(fabs(value3)+1.0))
    return 0.0;
  return value-value2;
}
//-------------------------------------------------------------------
// Returns the greatest common denominator of two 
// positive integers, a and b, found using Euclid's algorithm 
//-------------------------------------------------------------------
static int gcd(int a, int b) 
{
  int remainder = -1;
#if CGL_DEBUG>1
  printf("gcd of %d and %d\n",a,b);
  int nLoop=0;
#endif
  // make sure a<=b (will always remain so)
  if(a > b) {
    // Swap a and b
    int temp = a;
    a = b;
    b = temp;
  }
  // if zero then gcd is nonzero (zero may occur in rhs of packed)
  if (!a) {
    if (b) {
      return b;
    } else {
      printf("**** gcd given two zeros!!\n");
      abort();
    }
  }
  while (remainder) {

#if CGL_DEBUG>1
    nLoop++;
    if (nLoop>50) {
      abort();
      return -1;
    }
#endif
    remainder = b % a;
    b = a;
    a = remainder;
  }
#if CGL_DEBUG>1
  printf("=> %d\n",b);
#endif
  return b;
}

//-------------------------------------------------------------------
// Returns the nearest rational with denominator < maxDenominator
//-------------------------------------------------------------------
typedef struct {
  int numerator;
  int denominator;
} Rational;
inline Rational nearestRational(double value, int maxDenominator) 
{
  Rational tryThis;
  Rational tryA;
  Rational tryB;
  double integerPart;

#if CGL_DEBUG>1
  printf("Rational of %g is ",value);
#endif
  int nLoop=0;

  tryA.numerator=0;
  tryA.denominator=1;
  tryB.numerator=1;
  tryB.denominator=0;

  if (fabs(value)<1.0e-10)
    return tryA;
  integerPart = floor(value);
  value -= integerPart;
  tryThis.numerator = tryB.numerator* static_cast<int> (integerPart) + tryA.numerator;
  tryThis.denominator = tryB.denominator* static_cast<int> (integerPart) + tryA.denominator;
  tryA = tryB;
  tryB = tryThis;

  while (value>1.0e-10 && tryB.denominator <=maxDenominator) {
    nLoop++;
    if (nLoop>50) {
      Rational bad;
      bad.numerator=-1;
      bad.denominator=-1;
#if CGL_DEBUG>1
      printf(" *** bad rational\n");
#endif
      return bad;
    }
    value = 1.0/value;
    integerPart = floor(value+1.0e-10);
    value -= integerPart;
    tryThis.numerator = tryB.numerator* static_cast<int> (integerPart) + tryA.numerator;
    tryThis.denominator = tryB.denominator* static_cast<int>(integerPart) + tryA.denominator;
    tryA = tryB;
    tryB = tryThis;
  }
  if (tryB.denominator <= maxDenominator) {
#if CGL_DEBUG>1
    printf("%d/%d\n",tryB.numerator,tryB.denominator);
#endif
    return tryB;
  } else {
#if CGL_DEBUG>1
    printf("%d/%d\n",tryA.numerator,tryA.denominator);
#endif
    return tryA;
  }
}
// Does actual work - returns number of cuts
int
CglGomory::generateCuts( 
#ifdef CGL_DEBUG
			const OsiRowCutDebugger * debugger, 
#else
			const OsiRowCutDebugger * , 
#endif
                         OsiCuts & cs,
                         const CoinPackedMatrix & columnCopy,
                         const CoinPackedMatrix & rowCopy,
                         const double * colsol,
                         const double * colLower, const double * colUpper,
                         const double * rowLower, const double * rowUpper,
			 const char * intVar,
                         const CoinWarmStartBasis* warm,
                         const CglTreeInfo info) const
{
  int infoOptions=info.options;
  bool globalCuts = (infoOptions&16)!=0;
  double testFixed = (!globalCuts) ? 1.0e-8 : -1.0;
  // get what to look at
  double away = info.inTree ? away_ : CoinMin(away_,awayAtRoot_);
  int numberRows=columnCopy.getNumRows();
  int numberColumns=columnCopy.getNumCols(); 
  int numberElements=columnCopy.getNumElements();
  // Allow bigger length on initial matrix (if special setting)
  //if (limit==512&&!info.inTree&&!info.pass)
  //limit=1024;
  // Start of code to create a factorization from warm start (A) ====
  // check factorization is okay
  CoinFactorization factorization;
#ifdef CLP_OSL
  CoinOslFactorization * factorization2=NULL;
  if (alternateFactorization_) {
    factorization2 = new CoinOslFactorization();
  }
#endif
  // We can either set increasing rows so ...IsBasic gives pivot row
  // or we can just increment iBasic one by one
  // for now let ...iBasic give pivot row
  int status=-100;
  // probably could use pivotVariables from OsiSimplexModel
  int * rowIsBasic = new int[numberRows];
  int * columnIsBasic = new int[numberColumns];
  int i;
  int numberBasic=0;
  for (i=0;i<numberRows;i++) {
    if (warm->getArtifStatus(i) == CoinWarmStartBasis::basic) {
      rowIsBasic[i]=1;
      numberBasic++;
    } else {
      rowIsBasic[i]=-1;
    }
  }
  for (i=0;i<numberColumns;i++) {
    if (warm->getStructStatus(i) == CoinWarmStartBasis::basic) {
      columnIsBasic[i]=1;
      numberBasic++;
    } else {
      columnIsBasic[i]=-1;
    }
  }
  //returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
  while (status<-98) {
#ifdef CLP_OSL
    if (!alternateFactorization_) {
#endif
      status=factorization.factorize(columnCopy,
				     rowIsBasic, columnIsBasic);
      if (status==-99) factorization.areaFactor(factorization.areaFactor() * 2.0);
#ifdef CLP_OSL
    } else {
      double areaFactor=1.0;
      status=factorization2->factorize(columnCopy,
				      rowIsBasic, columnIsBasic,areaFactor);
      if (status==-99) areaFactor *= 2.0;
    }
#endif
  } 
  if (status) {
#ifdef COIN_DEVELOP
    std::cout<<"Bad factorization of basis - status "<<status<<std::endl;
#endif
    delete [] rowIsBasic;
    delete [] columnIsBasic;
    return -1;
  }
  // End of creation of factorization (A) ====
  
#ifdef CLP_OSL
  double relaxation = !alternateFactorization_ ? factorization.conditionNumber() :
    factorization2->conditionNumber();
#else
  double relaxation = factorization.conditionNumber();
#endif
#ifdef COIN_DEVELOP_z
  if (relaxation>1.0e49)
    printf("condition %g\n",relaxation);
#endif
  relaxation *= conditionNumberMultiplier_;
  double bounds[2]={-COIN_DBL_MAX,0.0};
  int iColumn,iRow;

  const int * column = rowCopy.getIndices();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  const int * rowLength = rowCopy.getVectorLengths(); 
  const double * rowElements = rowCopy.getElements();
  const int * row = columnCopy.getIndices();
  const CoinBigIndex * columnStart = columnCopy.getVectorStarts();
  const int * columnLength = columnCopy.getVectorLengths(); 
  const double * columnElements = columnCopy.getElements();

  // we need to do book-keeping for variables at ub
  double tolerance = 1.0e-7;
  bool * swap= new bool [numberColumns];
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (columnIsBasic[iColumn]<0&&
	colUpper[iColumn]-colsol[iColumn]<=tolerance) {
      swap[iColumn]=true;
    } else {
      swap[iColumn]=false;
    }
  }

  // get row activities (could use solver but lets do here )
  double * rowActivity = new double [numberRows];
  CoinFillN(rowActivity,numberRows,0.0);
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double value = colsol[iColumn];
    int k;
    for (k=columnStart[iColumn];k<columnStart[iColumn]+columnLength[iColumn];k++) {
      iRow = row[k];
      rowActivity[iRow] += columnElements[k]*value;
    }
  }
  /* we need to mark rows:
     0) equality
     1) slack at lb (activity at ub)
     2) slack at ub (activity at lb)
     and 4 bit is set if slack must be integer

  */
  int * rowType = new int[numberRows];
  for (iRow=0;iRow<numberRows;iRow++) {
    if (rowIsBasic[iRow]<0&&rowUpper[iRow]>rowLower[iRow]+1.0e-7) {
      int type=0;
      double rhs=0.0;
      if (rowActivity[iRow]>=rowUpper[iRow]-1.0e-7) {
	type=1;
	rhs=rowUpper[iRow];
      } else if (rowActivity[iRow]<=rowLower[iRow]+1.0e-7) {
	type=2;
	rhs=rowLower[iRow];
      } else {
	// probably large rhs
	if (rowActivity[iRow]-rowLower[iRow]<
	    rowUpper[iRow]-rowActivity[iRow])
	  rowType[iRow]=2;
	else
	  rowType[iRow]=1;
#ifdef CGL_DEBUG
	assert (CoinMin(rowUpper[iRow]-rowActivity[iRow],
		    rowActivity[iRow]-rowUpper[iRow])<1.0e-5);
	//abort();
        continue;
#else
	continue;
#endif
      }
      if (above_integer(rhs)<1.0e-10) {
	// could be integer slack
	bool allInteger=true;
	int k;
	for (k=rowStart[iRow];
	     k<rowStart[iRow]+rowLength[iRow];k++) {
	  int iColumn=column[k];
	  if (!intVar[iColumn]||above_integer(rowElements[k])>1.0e-10) {
	    // not integer slacks
	    allInteger=false;
	    break;
	  }
	}
	if (allInteger) {
	  type |= 4;
	}
      }
      rowType[iRow]=type;
    } else {
      // row is equality or basic
      rowType[iRow]=0;
    }
  }

  // Start of code to create work arrays for factorization (B) ====
  // two vectors for updating (one is work)
  CoinIndexedVector work;
  CoinIndexedVector array;
  // make sure large enough
  work.reserve(numberRows);
  array.reserve(numberRows);
  int * arrayRows = array.getIndices();
  double * arrayElements = array.denseVector();
  // End of code to create work arrays (B) ====

  int numberAdded=0;
  // we also need somewhere to accumulate cut
  CoinIndexedVector cutVector;
  cutVector.reserve(numberColumns+1);
  int * cutIndex = cutVector.getIndices();
  double * cutElement = cutVector.denseVector(); 
  // and for packed form (as not necessarily in order)
  // also space for sort
  bool doSorted = (infoOptions&256)!=0;
  int lengthArray = static_cast<int>(numberColumns+1+((numberColumns+1)*sizeof(int))/sizeof(double));
  if (doSorted)
    lengthArray+=numberColumns;
  double * packed = new double[lengthArray]; 
  double * sort = packed+numberColumns+1;
  int * which = reinterpret_cast<int *>(doSorted ? (sort+numberColumns): (sort));
  double tolerance1=1.0e-6;
  double tolerance2=0.9;
  double tolerance3=1.0e-4;
  double tolerance6=1.0e-6;
  double tolerance9=1.0e-4;
#define MORE_GOMORY_CUTS 1
#ifdef CLP_INVESTIGATE2
  int saveLimit = info.inTree ? 50 : 1000;
#else
#if MORE_GOMORY_CUTS==2||MORE_GOMORY_CUTS==3
  int saveLimit;
#endif  
#endif  
  // get limit on length of cut
  int limit = 0;
  if (!limit_)
    dynamicLimitInTree_ = CoinMax(50,numberColumns/40);
  if (!info.inTree) {
    limit = limitAtRoot_;
    if (!info.pass) {
      tolerance1=1.0;
      tolerance2=1.0e-2;
      tolerance3=1.0e-6;
      tolerance6=1.0e-7;
      tolerance9=1.0e-5;
      if (!limit||limit>=500)
	limit=numberColumns;
    } else {
      if((infoOptions&32)==0/*&&numberTimesStalled_<3*/) {
	if (!limit) {
	  if(numberElements>8*numberColumns)
	    limit=numberColumns;
	  else
	    limit = CoinMax(1000,numberRows/4);
	}
      } else {
	limit=numberColumns;
	numberTimesStalled_++;
      } 
    }
  } else {
    limit = limit_;
    if (!limit) {
      if (!info.pass) 
	limit = dynamicLimitInTree_;
      else
	limit=50;
    }
  }
  // If big - allow for rows
  if (limit>=numberColumns)
    limit += numberRows;
#ifdef CLP_INVESTIGATE2
  if (limit>saveLimit&&!info.inTree&&(infoOptions&512)==0) 
    printf("Gomory limit changed from %d to %d, inTree %c, pass %d, r %d,c %d,e %d\n",
	   saveLimit,limit,info.inTree ? 'Y' : 'N',info.pass,
	   numberRows,numberColumns,numberElements);
#endif
  int nCandidates=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // This returns pivot row for columns or -1 if not basic (C) ====
    int iBasic=columnIsBasic[iColumn];
    if (iBasic>=0&&intVar[iColumn]) {
      double reducedValue=above_integer(colsol[iColumn]);
      //printf("col %d bas %d val %.18g\n",iColumn,iBasic,colsol[iColumn]);
      if(intVar[iColumn]&&reducedValue<1.0-away&&reducedValue>away) {
	if (doSorted)
	  sort[nCandidates]=fabs(0.5-reducedValue);
	which[nCandidates++]=iColumn;
      }
    }
  }
  int nTotalEls=COIN_INT_MAX;
  if (doSorted) {
    CoinSort_2(sort,sort+nCandidates,which);
    int nElsNow = columnCopy.getNumElements();
    int nAdd;
    int nAdd2;
    int nReasonable;
    int depth=info.level;
    if (depth<2) {
      nAdd=10000;
      if (info.pass>0)
	nAdd = CoinMin(nAdd,nElsNow+2*numberRows);
      nAdd2 = 5*numberColumns;
      nReasonable = CoinMax(nAdd2,nElsNow/8+nAdd);
      if (!depth&&!info.pass) {
	// allow more
	nAdd += nElsNow/2;
	nAdd2 += nElsNow/2;
	nReasonable += nElsNow/2;
	limit=numberRows+numberColumns;
      }
    } else {
      nAdd = 200;
      nAdd2 = 2*numberColumns;
      nReasonable = CoinMax(nAdd2,nElsNow/8+nAdd);
    }
    nTotalEls=nReasonable;
  }
#ifdef MORE_GOMORY_CUTS
  int saveTotalEls=nTotalEls;
#endif
#if MORE_GOMORY_CUTS==2||MORE_GOMORY_CUTS==3
  saveLimit=limit;
  if (doSorted)
    limit=numberRows+numberColumns;
  OsiCuts longCuts;
#endif
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
  OsiCuts secondaryCuts;
#endif
  for (int kColumn=0;kColumn<nCandidates;kColumn++) {
    if (nTotalEls<=0)
      break;  // Got enough
    iColumn=which[kColumn];
    double reducedValue=above_integer(colsol[iColumn]);;
    // This returns pivot row for columns or -1 if not basic (C) ====
    int iBasic=columnIsBasic[iColumn];
    double ratio=reducedValue/(1.0-reducedValue);
    if (iBasic>=0) {
  // Debug code below computes tableau column of basic ====
      int j;
#ifdef CGL_DEBUG
      {
	// put column into array
	array.setVector(columnLength[iColumn],row+columnStart[iColumn],
			columnElements+columnStart[iColumn]);
	// get column in tableau
#ifdef CLP_OSL
	if (!alternateFactorization_)
#endif
	  factorization.updateColumn ( &work, &array );
#ifdef CLP_OSL
	else
	  factorization2->updateColumn ( &work, &array );
#endif
	int nn=0;
	int numberInArray=array.getNumElements();
	for (j=0;j<numberInArray;j++) {
	  int indexValue=arrayRows[j];
	  double value=arrayElements[indexValue];
	  if (fabs(value)>1.0e-5) {
	    assert (fabs(value-1.0)<1.0e-7);
	    assert (indexValue==iBasic);
	    nn++;
	  }
	}
	assert (nn==1);
	array.clear();
	work.checkClear();
      }
#endif
      array.clear();
      assert(intVar[iColumn]&&reducedValue<1.0-away&&reducedValue>away);
      {
#ifdef CGL_DEBUG
	cutVector.checkClear();
#endif
	// get row of tableau
	double one =1.0;
	array.setVector(1,&iBasic,&one);
	int numberNonInteger=0;
	//Code below computes tableau row ====
	// get pi
#ifdef CLP_OSL
	if (!alternateFactorization_)
#endif
	  factorization.updateColumnTranspose ( &work, &array );
#ifdef CLP_OSL
	else
	  factorization2->updateColumnTranspose ( &work, &array );
#endif
	int numberInArray=array.getNumElements();
#ifdef CGL_DEBUG
	// check pivot on iColumn
	{
	  double value=0.0;
	  int k;
	  // add in row of tableau
	  for (k=columnStart[iColumn];
	       k<columnStart[iColumn]+columnLength[iColumn];k++) {
	    iRow = row[k];
	    value += columnElements[k]*arrayElements[iRow];
	  }
	  // should be 1
	  assert (fabs(value-1.0) < 1.0e-7);
	}
#endif
	double largestFactor=0.0;
	for (j=0;j<numberInArray;j++) {
	  int indexValue=arrayRows[j];
	  double value=arrayElements[indexValue];
	  largestFactor = CoinMax(largestFactor,fabs(value));
	}
	//reducedValue=colsol[iColumn];
	// coding from pg 130 of Wolsey 
	// adjustment to rhs
	double rhs=0.0;
	int number=0;
#ifdef CGL_DEBUG_GOMORY
	    if (!gomory_try)
	      printf("start for basic column %d\n",iColumn);
#endif
	// columns
	for (j=0;j<numberColumns;j++) {
	  if (columnIsBasic[j]<0&&colUpper[j]>colLower[j]+testFixed) {
	    double value=0.0;
	    int k;
	    // add in row of tableau
	    for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
	      iRow = row[k];
	      double value2 = columnElements[k]*arrayElements[iRow];
	      largestFactor = CoinMax(largestFactor,fabs(value2));
	      value += value2;
	    }
#ifdef CGL_DEBUG_GOMORY
	    if (!gomory_try&&value)
	      printf("col %d alpha %g colsol %g swap %c bounds %g %g\n",
		     j,value,colsol[j],swap[j] ? 'Y' : 'N',
		     colLower[j],colUpper[j]);
#endif
	    // value is entry in tableau row end (C) ====
	    if (fabs(value)<1.0e-16) {
	      // small value
	      continue;
	    } else {
	      // left in to stop over compilation?
	      //if (iColumn==-52) printf("for basic %d, column %d has alpha %g, colsol %g\n",
	      //		      iColumn,j,value,colsol[j]);
#if CGL_DEBUG>1
	      if (iColumn==52) printf("for basic %d, column %d has alpha %g, colsol %g\n",
				      iColumn,j,value,colsol[j]);
#endif
	      // deal with bounds
	      if (swap[j]) {
		//reducedValue -= value*colUpper[j];
		// negate
		value = - value;
	      } else {
		//reducedValue -= value*colLower[j];
	      }
#if CGL_DEBUG>1
	      if (iColumn==52) printf("%d value %g reduced %g int %d rhs %g swap %d\n",
				      j,value,reducedValue,intVar[j],rhs,swap[j]);
#endif
	      double coefficient;
	      if (intVar[j]) {
		// integer
		coefficient = above_integer(value);
		if (coefficient > reducedValue) {
		  coefficient = ratio * (1.0-coefficient);
		} 
	      } else {
		// continuous
		numberNonInteger++;
		if (value > 0.0) {
		  coefficient = value;
		} else {
		  //??? sign wrong in book
		  coefficient = -ratio*value;
		}
	      }
	      if (swap[j]) {
		// negate
		coefficient = - coefficient;
		rhs += colUpper[j]*coefficient;
	      } else {
		rhs += colLower[j]*coefficient;
	      }
	      if (fabs(coefficient)>= COIN_INDEXED_TINY_ELEMENT) {
		cutElement[j] = coefficient;
		cutIndex[number++]=j;
		// If too many - break from loop
		if (number>limit) 
		  break;
	      }
	    }
	  } else {
	    // basic
	    continue;
	  }
	}
	cutVector.setNumElements(number);
	// If too many - just clear vector and skip
	if (number>limit) {
	  cutVector.clear();
	  continue;
	}
	//check will be cut
	//reducedValue=above_integer(reducedValue);
	rhs += reducedValue;
	double violation = reducedValue;
#ifdef CGL_DEBUG
	std::cout<<"cut has violation of "<<violation
		 <<" value "<<colsol[iColumn]<<std::endl;
#endif
	// now do slacks part
	for (j=0;j<numberInArray;j++) {
	  iRow=arrayRows[j];
	  double value = arrayElements[iRow];
	  int type=rowType[iRow];
	  if (type&&fabs(value)>=1.0e-16) {
	    if ((type&1)==0) {
	      // negate to get correct coefficient
	      value = - value;
	    }
	    double coefficient;
	    if ((type&4)!=0) {
	      // integer
	      coefficient = above_integer(value);
	      if (coefficient > reducedValue) {
		coefficient = ratio * (1.0-coefficient);
	      } 
	    } else {
	      numberNonInteger++;
	      // continuous
	      if (value > 0.0) {
		coefficient = value;
	      } else {
		coefficient = -ratio*value;
	      }
	    }
	    if ((type&1)!=0) {
	      // slack at ub - treat as +1.0
	      rhs -= coefficient*rowUpper[iRow];
	    } else {
	      // negate yet again ?
	      coefficient = - coefficient;
	      rhs -= coefficient*rowLower[iRow];
	    }
	    int k;
	    for (k=rowStart[iRow];
		 k<rowStart[iRow]+rowLength[iRow];k++) {
	      int jColumn=column[k];
	      double value=rowElements[k];
	      cutVector.quickAdd(jColumn,-coefficient*value);
	    }
	  }
	}
	//check again and pack down
	// also change signs
	// also zero cutElement
	double sum=0.0;
	rhs = - rhs;
	int n = cutVector.getNumElements();
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
	double violation2=violation;
#endif
	number=0;
	for (j=0;j<n;j++) {
	  int jColumn =cutIndex[j];
	  double value=-cutElement[jColumn];
	  cutElement[jColumn]=0.0;
	  if (fabs(value)>1.0e-8) {
	    sum+=value*colsol[jColumn];
	    packed[number]=value;
	    cutIndex[number++]=jColumn;
          } else {
#define LARGE_BOUND 1.0e20
            // small - adjust rhs if rhs reasonable
            if (value>0.0&&colLower[jColumn]>-LARGE_BOUND) {
              rhs -= value*colLower[jColumn];
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
	      // weaken violation
	      violation2 -= fabs(value*(colsol[jColumn]-colLower[jColumn]));
#endif
            } else if (value<0.0&&colUpper[jColumn]<LARGE_BOUND) {
              rhs -= value*colUpper[jColumn];
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
	      // weaken violation
	      violation2 -= fabs(value*(colsol[jColumn]-colUpper[jColumn]));
#endif
            } else if (fabs(value)>1.0e-13) {
              // take anyway
              sum+=value*colsol[jColumn];
              packed[number]=value;
              cutIndex[number++]=jColumn;
            } 
          }
	}
	// Final test on number
	//if (number>limit)
	//continue;
	// say zeroed out
	cutVector.setNumElements(0);
	bool accurate2=false;
	double difference=fabs((sum-rhs)-violation);
	double useTolerance;
	if (tolerance1>0.99) {
	  // use absolute
	  useTolerance = tolerance;
	} else {
	  double rhs2=CoinMax(fabs(rhs),10.0);
	  useTolerance=rhs2*0.1*tolerance1;
	}
	bool accurate = (difference<useTolerance);
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
	double difference2=fabs((sum-rhs)-violation2);
#if MORE_GOMORY_CUTS==1
	if (difference2<useTolerance&&doSorted) 
#else
	if (difference2<useTolerance&&doSorted&&number<saveLimit) 
#endif
	  accurate2=true;
#endif
	if (sum >rhs+tolerance2*away&&
	    (accurate||accurate2)) {
	  //#ifdef CGL_DEBUG
#ifdef CGL_DEBUG
#if CGL_DEBUG<=1
	  if (number<=-10) {
#endif
	    for (j=0;j<number;j++) {
	      std::cout<<" ["<<cutIndex[j]<<","<<packed[j]<<"]";
	    }
	    std::cout<<" <= "<<rhs<<std::endl;
#if CGL_DEBUG<=1
	  }
#endif
#endif
	  if (!numberNonInteger&&number) {
#ifdef CGL_DEBUG
	    assert (sizeof(Rational)==sizeof(double));
#endif
	    Rational * cleaned = reinterpret_cast<Rational *> (cutElement);
	    int * xInt = reinterpret_cast<int *> (cutElement);
	    // cut should have an integer slack so try and simplify
	    // add in rhs and put in cutElements (remember to zero later)
	    cutIndex[number]=numberColumns+1;
	    packed[number]=rhs;
	    int numberNonSmall=0;
	    int lcm = 1;
	    
	    for (j=0;j<number+1;j++) {
	      double value=above_integer(fabs(packed[j]));
	      if (fabs(value)<tolerance3) {
		// too small
		continue;
	      } else {
		numberNonSmall++;
	      }
	      
	      cleaned[j]=nearestRational(value,100000);
	      if (cleaned[j].denominator<0) {
		// bad rational
		lcm=-1;
		break;
	      }
	      int thisGcd = gcd(lcm,cleaned[j].denominator);
	      // may need to check for overflow
	      lcm /= thisGcd;
	      lcm *= cleaned[j].denominator;
	    }
	    if (lcm>0&&numberNonSmall) {
	      double multiplier = lcm;
	      int nOverflow = 0; 
	      for (j=0; j<number+1;j++) {
		double value = fabs(packed[j]);
		double dxInt = value*multiplier;
		xInt[j]= static_cast<int> (dxInt+0.5); 
#if CGL_DEBUG>1
		printf("%g => %g   \n",value,dxInt);
#endif
		if (dxInt>1.0e9||fabs(dxInt-xInt[j])> 1.0e-8) {
		  nOverflow++;
		  break;
		}
	      }
	      
	      if (nOverflow){
#ifdef CGL_DEBUG
		printf("Gomory Scaling: Warning: Overflow detected \n");
#endif
		numberNonInteger=-1;
	      } else {
		
		// find greatest common divisor of the elements
		j=0;
		while (!xInt[j])
		  j++; // skip zeros
		int thisGcd = gcd(xInt[j],xInt[j+1]);
		j++;
		for (;j<number+1;j++) {
		  if (xInt[j])
		    thisGcd = gcd(thisGcd,xInt[j]);
		}
#if 0
		// Check nothing too illegal - FIX this
		for (j=0;j<number+1;j++) {
		  double old = lcm*packed[j];
		  int newOne;
		  if (old>0.0)
		    newOne=xInt[j]/thisGcd;
		  else
		    newOne=-xInt[j]/thisGcd;
		  if (fabs(((double) newOne)-old)>
		      1.0e-10*(fabs(newOne)+1.0)) {
		    // say no good - first see if happens
		    printf("Fix this test 456 - just skip\n");
		    abort();
		  }
		} 
#endif		  
#if CGL_DEBUG>1
		printf("The gcd of xInt is %i\n",thisGcd);    
#endif
		
		// construct new cut by dividing through by gcd and 
		double minMultiplier=1.0e100;
		double maxMultiplier=0.0;
		for (j=0; j<number+1;j++) {
		  double old=packed[j];
		  if (old>0.0) {
		    packed[j]=xInt[j]/thisGcd;
		  } else {
		    packed[j]=-xInt[j]/thisGcd;
		  }
#if CGL_DEBUG>1
		  printf("%g => %g   \n",old,packed[j]);
#endif
		  if (packed[j]) {
		    if (fabs(packed[j])>maxMultiplier*fabs(old))
		      maxMultiplier = packed[j]/old;
		    if (fabs(packed[j])<minMultiplier*fabs(old))
		      minMultiplier = packed[j]/old;
		  }
		}
		rhs = packed[number];
#ifdef CGL_DEBUG
		printf("min, max multipliers - %g, %g\n",
		       minMultiplier,maxMultiplier);
#endif
		assert(maxMultiplier/minMultiplier>0.9999&&maxMultiplier/minMultiplier<1.0001);
	      }
	    }
	    // erase cutElement
	    CoinFillN(cutElement,number+1,0.0);
	  } else {
	    // relax rhs a tiny bit
	    rhs += 1.0e-8;
	    // relax if lots of elements for mixed gomory
	    if (number>=20) {
	      rhs  += 1.0e-7*(static_cast<double> (number/20));
	    }
	  }
	  // Take off tiny elements
	  // for first pass reject
#define TINY_ELEMENT 1.0e-12
	  {
	    int i,number2=number;
	    number=0;
	    double largest=0.0;
	    double smallest=1.0e30;
	    for (i=0;i<number2;i++) {
	      double value=fabs(packed[i]);
	      if (value<TINY_ELEMENT) {
		int iColumn = cutIndex[i];
		if (colUpper[iColumn]-colLower[iColumn]<LARGE_BOUND) {
		  // weaken cut
		  if (packed[i]>0.0) 
		    rhs -= value*colLower[iColumn];
		  else
		    rhs += value*colUpper[iColumn];
		} else {
		  // throw away
		  number=limit+1;
		  numberNonInteger=1;
		  break;
		}
	      } else {
		int iColumn = cutIndex[i];
		if (colUpper[iColumn]!=colLower[iColumn]||globalCuts) {
		  largest=CoinMax(largest,value);
		  smallest=CoinMin(smallest,value);
		  cutIndex[number]=cutIndex[i];
		  packed[number++]=packed[i];
		} else {
		  // fixed so subtract out
		  rhs -= packed[i]*colLower[iColumn];
		}
	      }
	    }
	    if (largest>1.0e10*smallest) {
	      number=limit+1; //reject
	      numberNonInteger=1;
	    } else if (largest>1.0e9*smallest) {
#ifdef CLP_INVESTIGATE2
	      printf("WOuld reject %g %g ratio %g\n",smallest,largest,
		     smallest/largest);
#endif
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
	      accurate=false;
#endif
	    }
	  }
	  if (number<limit||!numberNonInteger) {
	    bounds[1]=rhs;
	    if (number>50&&numberNonInteger)
	      bounds[1] = rhs+tolerance6+1.0e-8*fabs(rhs); // weaken
	    double test = CoinMin(largestFactor*largestFactorMultiplier_,
				  relaxation);
	    if (number>5&&numberNonInteger&&test>1.0e-20) {
#ifdef CLP_INVESTIGATE2
	      printf("relaxing rhs by %g - largestFactor was %g, rel %g\n",
	         CoinMin(test*fabs(rhs),tolerance9),largestFactor,relaxation);
#endif
	      //bounds[1] = CoinMax(bounds[1],
	      //		  rhs+CoinMin(test*fabs(rhs),tolerance9)); // weaken
	      bounds[1] = bounds[1]+CoinMin(test*fabs(rhs),tolerance9); // weaken
	    }
#ifdef MORE_GOMORY_CUTS
	    if (accurate) {
#else
	    {
#endif
	      OsiRowCut rc;
	      rc.setRow(number,cutIndex,packed,false);
	      rc.setLb(bounds[0]);
	      rc.setUb(bounds[1]);   
#if MORE_GOMORY_CUTS<2
	      nTotalEls -= number;
	      cs.insert(rc);
#else
	      if(number<saveLimit) {
		nTotalEls -= number;
		cs.insert(rc);
	      } else {
		longCuts.insert(rc);
	      }
#endif
	      //printf("nTot %d kCol %d iCol %d ibasic %d\n",
	      //     nTotalEls,kColumn,iColumn,iBasic);
	      numberAdded++;
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
	    } else if (accurate2) {
	      OsiRowCut rc;
	      rc.setRow(number,cutIndex,packed,false);
	      rc.setLb(bounds[0]);
	      rc.setUb(bounds[1]);   
	      secondaryCuts.insert(rc);
#endif
	    }
	  } else {
#ifdef CGL_DEBUG
	    std::cout<<"cut has "<<number<<" entries - skipped"
		     <<std::endl;
	    if (!number)
	      std::cout<<"******* Empty cut - infeasible"<<std::endl;
#endif
	  }
	} else {
	  // why dropped?
#ifdef CGL_DEBUG
	  std::cout<<"first violation "<<violation<<" now "
		   <<sum-rhs<<" why?, rhs= "<<rhs<<std::endl;
	  
	  for (j=0;j<number;j++) {
	    int jColumn =cutIndex[j];
	    double value=packed[j];
	    std::cout<<"("<<jColumn<<","<<value<<","<<colsol[jColumn]
		     <<") ";
	  }
	  std::cout<<std::endl;
	  //abort();
#endif
	}
      }
    } else {
      // not basic
#if CGL_DEBUG>1
      {
	// put column into array
	array.setVector(columnLength[iColumn],row+columnStart[iColumn],
			columnElements+columnStart[iColumn]);
	// get column in tableau
#ifdef CLP_OSL
	if (!alternateFactorization_)
#endif
	  factorization.updateColumn ( &work, &array );
#ifdef CLP_OSL
	else
	  factorization2->updateColumn ( &work, &array );
#endif
	int numberInArray=array.getNumElements();
	printf("non-basic %d\n",iColumn);
	for (int j=0;j<numberInArray;j++) {
	  int indexValue=arrayRows[j];
	  double value=arrayElements[indexValue];
	  if (fabs(value)>1.0e-6) {
	    printf("%d %g\n",indexValue,value);
	  }
	}
      }
#endif
    }
  }
#ifdef CLP_OSL
  delete factorization2;
#endif

  delete [] rowActivity;
  delete [] swap;
  delete [] rowType;
  delete [] packed;
  delete [] rowIsBasic;
  delete [] columnIsBasic;
#ifdef MORE_GOMORY_CUTS
#if MORE_GOMORY_CUTS==1
  int numberInaccurate = secondaryCuts.sizeRowCuts();
#ifdef CLP_INVESTIGATE2
  int numberOrdinary = numberAdded-numberInaccurate;
  if (!info.inTree&&(infoOptions&512)==0) 
    printf("Gomory has %d ordinary and %d less accurate cuts(%d els)\n",
	   numberOrdinary,numberInaccurate,saveTotalEls-nTotalEls);
#endif
#elif MORE_GOMORY_CUTS==2
  int numberLong = longCuts.sizeRowCuts();
#ifdef CLP_INVESTIGATE2
  int numberOrdinary = numberAdded-numberLong;
  if (!info.inTree&&(infoOptions&512)==0) 
    printf("Gomory has %d ordinary and %d long cuts(%d els)\n",
	   numberOrdinary,numberLong,saveTotalEls-nTotalEls);
#endif
#elif MORE_GOMORY_CUTS==3
  int numberLong = longCuts.sizeRowCuts();
  int numberInaccurate = secondaryCuts.sizeRowCuts();
#ifdef CLP_INVESTIGATE2
  int numberOrdinary = numberAdded-numberLong-numberInaccurate;
  if (!info.inTree&&(infoOptions&512)==0) 
    printf("Gomory has %d ordinary, %d long and %d less accurate cuts(%d els)\n",
	   numberOrdinary,numberLong,numberInaccurate,saveTotalEls-nTotalEls);
#endif
#endif
  if (doSorted&&limit<numberColumns) {
    // Just half
    nTotalEls -= saveTotalEls/2;
#if MORE_GOMORY_CUTS==2||MORE_GOMORY_CUTS==3
    while (nTotalEls>0) {
      for (int i=0;i<numberLong;i++) {
	nTotalEls -= longCuts.rowCutPtr(i)->row().getNumElements();
	cs.insert(longCuts.rowCut(i));
	numberAdded ++;
	if (nTotalEls<=0)
	  break;
      }
      break;
    }
#endif
#if MORE_GOMORY_CUTS==1||MORE_GOMORY_CUTS==3
    while (nTotalEls>0) {
      for (int i=0;i<numberInaccurate;i++) {
	nTotalEls -= secondaryCuts.rowCutPtr(i)->row().getNumElements();
	cs.insert(secondaryCuts.rowCut(i));
	numberAdded ++;
	if (nTotalEls<=0)
	  break;
      }
      break;
    }
#endif
  }
#else
#ifdef CLP_INVESTIGATE2
  if (!info.inTree&&(infoOptions&512)==0) 
    printf("Gomory added %d cuts(%d els)\n",numberAdded,saveTotalEls-nTotalEls);
#endif
#endif
  return numberAdded;
}
// Limit stuff
void CglGomory::setLimit(int limit)
{
  if (limit>=0)
    limit_=limit; 
}
int CglGomory::getLimit() const
{
  return limit_;
}
// Limit stuff at root
void CglGomory::setLimitAtRoot(int limit)
{
  if (limit>=0)
    limitAtRoot_=limit;
}
int CglGomory::getLimitAtRoot() const
{
  return limitAtRoot_;
}
// Return maximum length of cut in tree
int 
CglGomory::maximumLengthOfCutInTree() const
{
  if (limit_)
    return limit_;
  else
    return dynamicLimitInTree_;
}

// Away stuff
void CglGomory::setAway(double value)
{
  if (value>0.0&&value<=0.5)
    away_=value;
}
double CglGomory::getAway() const
{
  return away_;
}

// Away stuff at root
void CglGomory::setAwayAtRoot(double value)
{
  if (value>0.0&&value<=0.5)
    awayAtRoot_=value;
}
double CglGomory::getAwayAtRoot() const
{
  return awayAtRoot_;
}

// ConditionNumberMultiplier stuff
void CglGomory::setConditionNumberMultiplier(double value)
{
  if (value>=0.0)
    conditionNumberMultiplier_=value;
}
double CglGomory::getConditionNumberMultiplier() const
{
  return conditionNumberMultiplier_;
}

// LargestFactorMultiplier stuff
void CglGomory::setLargestFactorMultiplier(double value)
{
  if (value>=0.0)
    largestFactorMultiplier_=value;
}
double CglGomory::getLargestFactorMultiplier() const
{
  return largestFactorMultiplier_;
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglGomory::CglGomory ()
:
CglCutGenerator(),
away_(0.05),
awayAtRoot_(0.05),
conditionNumberMultiplier_(1.0e-18),
largestFactorMultiplier_(1.0e-13),
originalSolver_(NULL),
limit_(50),
limitAtRoot_(0),
dynamicLimitInTree_(-1),
alternateFactorization_(0),
gomoryType_(0)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglGomory::CglGomory (const CglGomory & source) :
  CglCutGenerator(source),
  away_(source.away_),
  awayAtRoot_(source.awayAtRoot_),
  conditionNumberMultiplier_(source.conditionNumberMultiplier_),
  largestFactorMultiplier_(source.largestFactorMultiplier_),
  originalSolver_(NULL),
  limit_(source.limit_),
  limitAtRoot_(source.limitAtRoot_),
  dynamicLimitInTree_(source.dynamicLimitInTree_),
  alternateFactorization_(source.alternateFactorization_),
  gomoryType_(source.gomoryType_)
{ 
  if (source.originalSolver_)
    originalSolver_ = source.originalSolver_->clone();
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglGomory::clone() const
{
  return new CglGomory(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglGomory::~CglGomory ()
{
  delete originalSolver_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglGomory &
CglGomory::operator=(const CglGomory& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    away_=rhs.away_;
    awayAtRoot_=rhs.awayAtRoot_;
    conditionNumberMultiplier_ = rhs.conditionNumberMultiplier_;
    largestFactorMultiplier_ = rhs.largestFactorMultiplier_;
    limit_=rhs.limit_;
    limitAtRoot_=rhs.limitAtRoot_;
    dynamicLimitInTree_ = rhs.dynamicLimitInTree_;
    alternateFactorization_=rhs.alternateFactorization_; 
    gomoryType_ = rhs.gomoryType_;
    delete originalSolver_;
    if (rhs.originalSolver_)
      originalSolver_ = rhs.originalSolver_->clone();
    else
      originalSolver_=NULL;
  }
  return *this;
}
// Pass in a copy of original solver (clone it)
void 
CglGomory::passInOriginalSolver(OsiSolverInterface * solver)
{
  delete originalSolver_;
  if (solver) {
    if (!gomoryType_)
      gomoryType_=1;
    originalSolver_ = solver->clone();
  } else {
    gomoryType_=0;
    originalSolver_=NULL;
  }
}
// Returns true if needs optimal basis to do cuts
bool 
CglGomory::needsOptimalBasis() const
{
  return true;
}
// Does actual work - returns number of cuts
int
CglGomory::generateCuts( const OsiRowCutDebugger * debugger, 
                         OsiCuts & cs,
                         const CoinPackedMatrix & columnCopy,
                         const double * colsol,
                         const double * colLower, const double * colUpper,
                         const double * rowLower, const double * rowUpper,
			 const char * intVar,
                         const CoinWarmStartBasis* warm,
                         const CglTreeInfo info) const
{
  CoinPackedMatrix rowCopy;
  rowCopy.reverseOrderedCopyOf(columnCopy);
  return generateCuts( debugger, cs, columnCopy, rowCopy,
		       colsol, colLower, colUpper,
		       rowLower, rowUpper, intVar, warm, info);
}
// Create C++ lines to get to current state
std::string
CglGomory::generateCpp( FILE * fp) 
{
  CglGomory other;
  fprintf(fp,"0#include \"CglGomory.hpp\"\n");
  fprintf(fp,"3  CglGomory gomory;\n");
  if (limit_!=other.limit_)
    fprintf(fp,"3  gomory.setLimit(%d);\n",limit_);
  else
    fprintf(fp,"4  gomory.setLimit(%d);\n",limit_);
  if (limitAtRoot_!=other.limitAtRoot_)
    fprintf(fp,"3  gomory.setLimitAtRoot(%d);\n",limitAtRoot_);
  else
    fprintf(fp,"4  gomory.setLimitAtRoot(%d);\n",limitAtRoot_);
  if (away_!=other.away_)
    fprintf(fp,"3  gomory.setAway(%g);\n",away_);
  else
    fprintf(fp,"4  gomory.setAway(%g);\n",away_);
  if (awayAtRoot_!=other.awayAtRoot_)
    fprintf(fp,"3  gomory.setAwayAtRoot(%g);\n",awayAtRoot_);
  else
    fprintf(fp,"4  gomory.setAwayAtRoot(%g);\n",awayAtRoot_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  gomory.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  gomory.setAggressiveness(%d);\n",getAggressiveness());
  return "gomory";
}
