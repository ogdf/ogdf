// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <string>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cfloat>

#include "CglPreProcess.hpp"
#include "CglMessage.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglStored.hpp"
#include "CglCutGenerator.hpp"
#include "CoinTime.hpp"
#include "CoinSort.hpp"
#include "CoinBuild.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinWarmStartBasis.hpp"

#include "CglProbing.hpp"
#include "CglDuplicateRow.hpp"
#include "CglClique.hpp"

OsiSolverInterface *
CglPreProcess::preProcess(OsiSolverInterface & model, 
                       bool makeEquality, int numberPasses)
{
  // Tell solver we are in Branch and Cut
  model.setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;
  // Default set of cut generators
  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbeRoot(model.getNumCols());
  generator1.setMaxElements(100);
  generator1.setMaxLookRoot(50);
  generator1.setRowCuts(3);
  // Add in generators
  addCutGenerator(&generator1);
  OsiSolverInterface * newSolver = preProcessNonDefault(model,makeEquality ? 1 : 0,numberPasses);
  // Tell solver we are not in Branch and Cut
  model.setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;
  if (newSolver)
    newSolver->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;
  return newSolver;
}
static void outSingletons(int & nCol, int & nRow,
			  int * startCol, int * row, double * element,
			  int * startRow, int *column)
{
  int iRow,iCol;
  bool singletons=false;
  int * countRow = new int [nRow];
  int * countCol = new int [nCol];
  int * temp = new int[nRow];
  // make row copy
  memset(countRow,0,nRow*sizeof(int));
  memset(countCol,0,nCol*sizeof(int));
  for (iCol=0;iCol<nCol;iCol++) {
    for (int j=startCol[iCol];j<startCol[iCol+1];j++) {
      int iRow = row[j];
      countRow[iRow]++;
      countCol[iCol]++;
    }
  }
  startRow[0]=0;
  for (iRow=0;iRow<nRow;iRow++) {
    int k = countRow[iRow]+startRow[iRow];
    temp[iRow]=startRow[iRow];
    startRow[iRow+1]=k;
  }
  for (iCol=0;iCol<nCol;iCol++) {
    for (int j=startCol[iCol];j<startCol[iCol+1];j++) {
      int iRow = row[j];
      int k=temp[iRow];
      temp[iRow]++;
      column[k]=iCol;
    }
  }
  for (iRow=0;iRow<nRow;iRow++) {
    if (countRow[iRow]<=1)
      singletons=true;
  }
  for (iCol=0;iCol<nCol;iCol++) {
    if (countCol[iCol]<=1)
      singletons=true;
  }
  if (singletons) {
    while (singletons) {
      singletons=false;
      for (iCol=0;iCol<nCol;iCol++) {
	if (countCol[iCol]==1) {
	  singletons=true;
	  countCol[iCol]=0;
	  int iRow = row[startCol[iCol]];
	  int start = startRow[iRow];
	  int end = start+countRow[iRow];
	  countRow[iRow]--;
	  int j;
	  for ( j=start;j<end;j++) {
	    if (column[j]==iCol) {
	      column[j]=column[end-1];
	      break;
	    }
	  }
	  assert (j<end);
	}
      }
      for (iRow=0;iRow<nRow;iRow++) {
	if (countRow[iRow]==1) {
	  singletons=true;
	  countRow[iRow]=0;
	  int iCol = column[startRow[iRow]];
	  int start = startCol[iCol];
	  int end = start+countCol[iCol];
	  countCol[iCol]--;
	  int j;
	  for ( j=start;j<end;j++) {
	    if (row[j]==iRow) {
	      row[j]=row[end-1];
	      if (element)
		element[j]=element[end-1];
	      break;
	    }
	  }
	  assert (j<end);
	}
      }
    }  
    // Pack down
    int newNrow=0;
    for (iRow=0;iRow<nRow;iRow++) {
      if (countRow[iRow]==0) {
	temp[iRow]=-1;
      } else {
	assert (countRow[iRow]>1);
	temp[iRow]=newNrow;
	newNrow++;
      }
    }
    int newNcol=0;
    int nEl=0;
    int iNext=0;
    for (iCol=0;iCol<nCol;iCol++) {
      int start = iNext;
      iNext=startCol[iCol+1];
      if (countCol[iCol]==0) {
	countCol[iCol]=-1;
      } else {
	assert (countCol[iCol]>1);
	int end = start+countCol[iCol];
	countCol[iCol]=newNcol;
	int j;
	for ( j=start;j<end;j++) {
	  int iRow=row[j];
	  iRow = temp[iRow];
	  assert (iRow>=0);
	  row[nEl]=iRow;
	  if (element)
	    element[nEl]=element[j];
	  nEl++;
	}
	newNcol++;
	startCol[newNcol]=nEl;
      }
    }
    newNrow=0;
    nEl=0;
    iNext=0;
    for (iRow=0;iRow<nRow;iRow++) {
      int start = iNext;
      iNext = startRow[iRow+1];
      if (countRow[iRow]>1) {
	int end = start+countRow[iRow];
	int j;
	for ( j=start;j<end;j++) {
	  int iCol = column[j];
	  iCol = countCol[iCol];
	  assert (iCol>=0);
	  column[nEl++]=iCol;
	}
	newNrow++;
	startRow[newNrow]=nEl;
      }
    }
    nRow=newNrow;
    nCol=newNcol;
  }
  delete [] countCol;
  delete [] countRow;
  delete [] temp;
}
static int makeIntegers2(OsiSolverInterface * model,int mode)
{
  // See whether we should make variables integer
  const double *objective = model->getObjCoefficients() ;
  const double *lower = model->getColLower() ;
  const double *upper = model->getColUpper() ;
  const double *rowLower = model->getRowLower() ;
  const double *rowUpper = model->getRowUpper() ;
  int numberRows = model->getNumRows() ;
  double * rhs = new double [numberRows];
  int * count = new int [numberRows];
  int iColumn;
  bool makeAll = (mode>1);
  int numberColumns = model->getNumCols() ;
  // Column copy of matrix
  const double * element = model->getMatrixByCol()->getElements();
  const int * row = model->getMatrixByCol()->getIndices();
  const CoinBigIndex * columnStart = model->getMatrixByCol()->getVectorStarts();
  const int * columnLength = model->getMatrixByCol()->getVectorLengths();
  // Row copy
  CoinPackedMatrix matrixByRow(*model->getMatrixByRow());
  //const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();
  int numberIntegers=1;
  int totalNumberIntegers=0;
  while (numberIntegers) {
    memset(rhs,0,numberRows*sizeof(double));
    memset(count,0,numberRows*sizeof(int));
    int currentNumber=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex start = columnStart[iColumn];
      CoinBigIndex end = start + columnLength[iColumn];
      if (upper[iColumn]==lower[iColumn]) {
	for (CoinBigIndex j=start;j<end;j++) {
	  int iRow = row[j];
	  rhs[iRow] += lower[iColumn]*element[j];
	}
      } else if (model->isInteger(iColumn)) {
	currentNumber++;
	for (CoinBigIndex j=start;j<end;j++) {
	  int iRow = row[j];
	  if (fabs(element[j]-floor(element[j]+0.5))>1.0e-10) 
	    rhs[iRow]  = COIN_DBL_MAX;
	}
      } else {
	for (CoinBigIndex j=start;j<end;j++) {
	  int iRow = row[j];
	  count[iRow]++;
	  if (fabs(element[j])!=1.0)
	    rhs[iRow]  = COIN_DBL_MAX;
	}
      }
    }
#ifdef COIN_DEVELOP
    printf("Current number of integers is %d\n",currentNumber);
#endif
    // now look at continuous
    bool allGood=true;
    double direction = model->getObjSense() ;
    int numberObj=0;
    int numberEq=0;
    int numberEqI=0;
    int numberZero=0;
    int numberNonZero=0;
    if (false) {
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (upper[iColumn]>lower[iColumn]) {
	  if (!model->isInteger(iColumn)) {
	    CoinBigIndex start = columnStart[iColumn];
	    CoinBigIndex end = start + columnLength[iColumn];
	    int nC=0;
	    for (CoinBigIndex j=start;j<end;j++) {
	      int iRow = row[j];
	      if (count[iRow]>1) {
		nC++;
	      }
	    }
	    if (nC>2) {
	      for (CoinBigIndex j=start;j<end;j++) {
		int iRow = row[j];
		if (count[iRow]>1) 
		  count[iRow]=999999;
	      }
	    }
	  }
	}
      }
    }
    int * newInts = new int[numberColumns];
    // Columns to zap
    int nColumnZap=0;
    int * columnZap = new int[numberColumns];
    char * noGoodColumn = new char [numberColumns];
    memset(noGoodColumn,0,numberColumns);
    int nNew=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn]>lower[iColumn]) {
	double objValue = objective[iColumn]*direction;
	bool thisGood=true;
	if ((objValue||makeAll)&&!model->isInteger(iColumn)) {
	  if (objValue) {
	    numberObj++;
	  } else if (columnLength[iColumn]==1) {
	    continue; // don't bother with singletons
	  }
	  CoinBigIndex start = columnStart[iColumn];
	  CoinBigIndex end = start + columnLength[iColumn];
	  if (objValue>=0.0) {
	    // wants to be as low as possible
	    if (lower[iColumn]<-1.0e10||fabs(lower[iColumn]-floor(lower[iColumn]+0.5))>1.0e-10) {
	      thisGood=false;
	    } else if (upper[iColumn]<1.0e10&&fabs(upper[iColumn]-floor(upper[iColumn]+0.5))>1.0e-10) {
	      thisGood=false;
	    }
	    bool singletonRow=true;
	    bool equality=false;
	    int nC=0;
	    int xxxx=0;
	    bool badCount=false;
	    for (CoinBigIndex j=start;j<end;j++) {
	      int iRow = row[j];
	      if (count[iRow]>1) {
		singletonRow=false;
		//printf("col %d row%d element %g - row count %d\n",iColumn,iRow,element[j],count[iRow]);
		if (count[iRow]==999999)
		  badCount=true;
		if (element[j]==1.0) {
		  if ((xxxx&1)==0)
		    xxxx |= 1;
		  else
		    xxxx = 15;
		} else {
		  if ((xxxx&2)==0)
		    xxxx |= 2;
		  else
		    xxxx = 15;
		}
		nC++;
	      } else if (rowLower[iRow]==rowUpper[iRow]) {
		equality=true;
	      }
	      double rhsValue = rhs[iRow];
	      double lowerValue = rowLower[iRow];
	      double upperValue = rowUpper[iRow];
	      if (rhsValue<1.0e20) {
		if(lowerValue>-1.0e20)
		  lowerValue -= rhsValue;
		if(upperValue<1.0e20)
		  upperValue -= rhsValue;
	      }
	      if (fabs(rhsValue)>1.0e20||fabs(rhsValue-floor(rhsValue+0.5))>1.0e-10
		  ||fabs(element[j])!=1.0) {
		// no good
		thisGood=false;
		break;
	      }
	      if (element[j]>0.0) {
		if (lowerValue>-1.0e20&&fabs(lowerValue-floor(lowerValue+0.5))>1.0e-10) {
		  // no good
		  thisGood=false;
		  break;
		}
	      } else {
		if (upperValue<1.0e20&&fabs(upperValue-floor(upperValue+0.5))>1.0e-10) {
		  // no good
		  thisGood=false;
		  break;
		}
	      }
	    }
	    if (!model->isInteger(iColumn)&&false)
	      printf("%d has %d rows with >1 - state network %s interaction %s\n",iColumn,nC,xxxx>3 ? "bad" : "good",
		     badCount ? "too much" : "ok");
	    // If not good here then mark rows
	    if (!thisGood) {
	      for (CoinBigIndex j=start;j<end;j++) {
		int iRow = row[j];
		count[iRow]=999999;
	      }
	    }
	    if (!singletonRow&&end>start+1&&!equality)
	      thisGood=false;
	    // Can we make equality
	    if (end==start+1&&!equality&&false) {
	      numberEq++;
	      int iRow = row[start];
	      if (element[start]>0.0)
		model->setRowUpper(iRow,rowLower[iRow]);
	      else
		model->setRowLower(iRow,rowUpper[iRow]);
	    }
	  } else {
	    // wants to be as high as possible
	    if (upper[iColumn]>1.0e10||fabs(upper[iColumn]-floor(upper[iColumn]+0.5))>1.0e-10) {
	      thisGood=false;
	    } else if (lower[iColumn]>-1.0e10&&fabs(lower[iColumn]-floor(lower[iColumn]+0.5))>1.0e-10) {
	      thisGood=false;
	    }
	    bool singletonRow=true;
	    bool equality=false;
	    for (CoinBigIndex j=start;j<end;j++) {
	      int iRow = row[j];
	      if (count[iRow]>1) {
		singletonRow=false;
		thisGood=false;
	      } else if (rowLower[iRow]==rowUpper[iRow]) {
		equality=true;
	      }
	      double rhsValue = rhs[iRow];
	      double lowerValue = rowLower[iRow];
	      double upperValue = rowUpper[iRow];
	      if (rhsValue<1.0e20) {
		if(lowerValue>-1.0e20)
		  lowerValue -= rhsValue;
		if(upperValue<1.0e20)
		  upperValue -= rhsValue;
	      }
	      if (fabs(rhsValue)>1.0e20||fabs(rhsValue-floor(rhsValue+0.5))>1.0e-10
		  ||fabs(element[j])!=1.0) {
		// no good
		thisGood=false;
		break;
	      }
	      if (element[j]<0.0) {
		if (lowerValue>-1.0e20&&fabs(lowerValue-floor(lowerValue+0.5))>1.0e-10) {
		  // no good
		  thisGood=false;
		  break;
		}
	      } else {
		if (upperValue<1.0e20&&fabs(upperValue-floor(upperValue+0.5))>1.0e-10) {
		  // no good
		  thisGood=false;
		  break;
		}
	      }
	    }
	    if (!singletonRow&&end>start+1&&!equality)
	      thisGood=false;
	    // If not good here then mark rows
	    if (!thisGood) {
	      for (CoinBigIndex j=start;j<end;j++) {
		int iRow = row[j];
		count[iRow]=999999;
	      }
	    }
	    // Can we make equality
	    if (end==start+1&&!equality&&false) {
	      numberEq++;
	      int iRow = row[start];
	      if (element[start]<0.0)
		model->setRowUpper(iRow,rowLower[iRow]);
	      else
		model->setRowLower(iRow,rowUpper[iRow]);
	    }
	  }
	} else if (objValue) {
	  CoinBigIndex start = columnStart[iColumn];
	  CoinBigIndex end = start + columnLength[iColumn];
	  if (end==start+1) {
	    int iRow = row[start];
	    if (rowUpper[iRow]>rowLower[iRow]&&!count[iRow]) {
	      if (fabs(rhs[iRow])>1.0e20||fabs(rhs[iRow]-floor(rhs[iRow]+0.5))>1.0e-10
		  ||fabs(element[start])!=1.0) {
		// no good
	      } else if (false) {
		numberEqI++;
		if (element[start]*objValue>0.0)
		  model->setRowUpper(iRow,rowLower[iRow]);
		else
		  model->setRowLower(iRow,rowUpper[iRow]);
	      }
	    }
	  }
	}
	if (!thisGood) {
	  if (objValue)
	    allGood=false;
	  // look at again
	  columnZap[nColumnZap++]=iColumn;
	} else if (makeAll&&!model->isInteger(iColumn)&&
		   upper[iColumn]-lower[iColumn]<10) {
	  newInts[nNew++] = iColumn;
	}
      }
    }
    // Rows to look at
    int * rowLook = new int[numberRows];
    while (nColumnZap) {
      int nRowLook=0;
      for (int i=0;i<nColumnZap;i++) {
	int iColumn=columnZap[i];
	noGoodColumn[iColumn]=1;
	CoinBigIndex start = columnStart[iColumn];
	CoinBigIndex end = start + columnLength[iColumn];
	for (CoinBigIndex j=start;j<end;j++) {
	  int iRow = row[j];
	  if (count[iRow]!=999999) {
	    count[iRow]=999999;
	    rowLook[nRowLook++]=iRow;
	  }
	}
      }
      nColumnZap=0;
      if (nRowLook) {
	for (int i=0;i<nRowLook;i++) {
	  int iRow=rowLook[i];
	  CoinBigIndex start = rowStart[iRow];
	  CoinBigIndex end = start + rowLength[iRow];
	  for (CoinBigIndex j=start;j<end;j++) {
	    int iColumn = column[j];
	    if (upper[iColumn]>lower[iColumn]&&
		!model->isInteger(iColumn)) {
	      if (!noGoodColumn[iColumn]) {
		noGoodColumn[iColumn]=1;
		columnZap[nColumnZap++]=iColumn;
	      }
	    }
	  }
	}
      }
    }
    delete [] rowLook;
    delete [] noGoodColumn;
    delete [] columnZap;
    // Final look
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn]>lower[iColumn]&&
	  !model->isInteger(iColumn)&&objective[iColumn]) {
	CoinBigIndex start = columnStart[iColumn];
	CoinBigIndex end = start + columnLength[iColumn];
	for (CoinBigIndex j=start;j<end;j++) {
	  int iRow = row[j];
	  if (count[iRow]==999999) 
	    allGood=false;
	}
      }
    }
    // do if some but not too many
    if (nNew&&nNew<currentNumber) {
      for (int i=0;i<nNew;i++) {
	int iColumn = newInts[i];
	double objValue = objective[iColumn];
	bool thisGood=true;
	CoinBigIndex start = columnStart[iColumn];
	CoinBigIndex end = start + columnLength[iColumn];
	for (CoinBigIndex j=start;j<end;j++) {
	  int iRow = row[j];
	  if (count[iRow]==999999) {
	    thisGood=false;
	    break;
	  }
	} 
	if (thisGood&&upper[iColumn]<lower[iColumn]+10.0) {
	  model->setInteger(iColumn);
	  if (objValue)
	    numberNonZero++;
	  else
	    numberZero++; 
	} else if (objValue) {
	  // unable to fix all with obj
	  allGood=false;
	}
      }
    }
    delete [] newInts;
    // Can we look at remainder and make any integer
    if (makeAll&&false) {
      int nLook=0;
      int nEl=0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (upper[iColumn]>lower[iColumn]&&!model->isInteger(iColumn)) {
	  CoinBigIndex start = columnStart[iColumn];
	  CoinBigIndex end = start + columnLength[iColumn];
	  bool possible=true;
	  int n=0;
	  for (CoinBigIndex j=start;j<end;j++) {
	    int iRow = row[j];
	    if (count[iRow]>1) {
	      if (count[iRow]==999999) {
		possible=false;
		break;
	      } else {
		n++;
	      }
	    }
	  }
	  if (possible) {
	    nLook++;
	    nEl+=n;
	  }
	}
      }
      if (nLook) {
	int * startC = new int [nLook+1];
	int * back = new int [nLook];
	int * row2 = new int[nEl];
	double * element2 = new double [nEl];
	int * backRow = new int [numberRows];
	int jRow;
	for (jRow=0;jRow<numberRows;jRow++) {
	  backRow[jRow]=-1;
	}
	int nCol=nLook;
	nLook=0;
	nEl=0;
	startC[0]=0;
	int nRow=0;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (upper[iColumn]>lower[iColumn]&&!model->isInteger(iColumn)) {
	    CoinBigIndex start = columnStart[iColumn];
	    CoinBigIndex end = start + columnLength[iColumn];
	    bool possible=true;
	    int n=0;
	    for (CoinBigIndex j=start;j<end;j++) {
	      int iRow = row[j];
	      if (count[iRow]>1) {
		if (count[iRow]==999999) {
		  possible=false;
		  break;
		} else {
		  n++;
		}
	      }
	    }
	    if (!n)
	      possible=false; // may be done later
	    if (possible) {
	      back[nLook]=iColumn;
	      for (CoinBigIndex j=start;j<end;j++) {
		int iRow = row[j];
		if (count[iRow]>1) {
		  int jRow=backRow[iRow];
		  if (jRow<0) {
		    // new row
		    backRow[iRow]=nRow;
		    jRow=nRow;
		    nRow++;
		  }
		  element2[nEl]=element[j];
		  row2[nEl++]=jRow;
		}
	      }
	      nLook++;
	      startC[nLook]=nEl;
	    }
	  }
	}
	// Redo nCol
	nCol=nLook;
	delete [] backRow;
	int * startRow = new int [nRow+1];
	int * column2 = new int [nEl];
	// take out singletons and do row copy
	outSingletons(nCol, nRow,
		      startC, row2, element2,
		      startRow, column2);
	// Decompose
	int * rowBlock = new int[nRow];
	int * stack = new int [nRow];
	for (int iRow=0;iRow<nRow;iRow++)
	  rowBlock[iRow]=-2;
	int numberBlocks = 0;
	// to say if column looked at
	int * columnBlock = new int[nCol];
	int iColumn;
	for (iColumn=0;iColumn<nCol;iColumn++)
	  columnBlock[iColumn]=-2;
	for (iColumn=0;iColumn<nCol;iColumn++) {
	  int kstart = startC[iColumn];
	  int kend = startC[iColumn+1];
	  if (columnBlock[iColumn]==-2) {
	    // column not allocated
	    int j;
	    int nstack=0;
	    for (j=kstart;j<kend;j++) {
	      int iRow= row2[j];
	      if (rowBlock[iRow]!=-1) {
		assert(rowBlock[iRow]==-2);
		rowBlock[iRow]=numberBlocks; // mark
		stack[nstack++] = iRow;
	      }
	    }
	    if (nstack) {
	      // new block - put all connected in
	      numberBlocks++;
	      columnBlock[iColumn]=numberBlocks-1;
	      while (nstack) {
		int iRow = stack[--nstack];
		int k;
		for (k=startRow[iRow];k<startRow[iRow+1];k++) {
		  int iColumn = column2[k];
		  int kkstart = startC[iColumn];
		  int kkend = startC[iColumn+1];
		  if (columnBlock[iColumn]==-2) {
		    columnBlock[iColumn]=numberBlocks-1; // mark
		    // column not allocated
		    int jj;
		    for (jj=kkstart;jj<kkend;jj++) {
		      int jRow= row2[jj];
		      if (rowBlock[jRow]==-2) {
			rowBlock[jRow]=numberBlocks-1;
			stack[nstack++]=jRow;
		      }
		    }
		  } else {
		    assert (columnBlock[iColumn]==numberBlocks-1);
		  }
		}
	      }
	    } else {
	      // Only in master
	      columnBlock[iColumn]=-1;
	      // empty - should already be integer
	      abort();
	    }
	  }
	}
	// See if each block OK 
	for (int iBlock=0;iBlock<numberBlocks;iBlock++) {
	  // Get block
	  int * startCB = new int [nCol+1];
	  int * row2B = new int[nEl];
	  int * startCC = new int [nCol+1];
	  int * row2C = new int[nEl];
	  int * startRowC = new int [nRow+1];
	  int * column2C = new int [nEl];
	  int * whichRow = new int [nRow];
	  int * whichCol = new int [nCol];
	  int i;
	  int nRowB=0;
	  int nColB=0;
	  int nElB=0;
	  for (i=0;i<nRow;i++) {
	    if (rowBlock[i]==iBlock) {
	      whichRow[i]=nRowB;
	      nRowB++;
	    } else {
	      whichRow[i]=-1;
	    }
	  }
	  bool network=true;
	  // even if not network - take out network columns NO
	  startCB[0]=0;
	  for (i=0;i<nCol;i++) {
	    if (columnBlock[i]==iBlock) {
	      int type=0;
	      whichCol[i]=nColB;
	      for (int j=startC[i];j<startC[i+1];j++) {
		int iRow = row2[j];
		iRow = whichRow[iRow];
		if (iRow>=0) { 
		  if (element2[j]==1.0) {
		    if ((type&1)==0)
		      type |=1;
		    else 
		      type=7;
		  } else {
		    assert (element2[j]==-1.0);
		    if ((type&2)==0)
		      type |=2;
		    else 
		      type=7;
		  }
		  row2B[nElB++]=iRow;
		}
	      }
	      if (type!=3) 
		network=false;
	      nColB++;
	      startCB[nColB]=nElB;
	      assert (startCB[nColB]>startCB[nColB-1]+1);
	    } else {
	      whichCol[i]=-1;
	    }
	  }
	  // See if network
	  bool goodInteger=false;
	  if (!network) {
	    // take out singletons
	    outSingletons(nColB, nRowB,
			  startCB, row2B, NULL,
			  startRowC, column2C);
	    // See if totally balanced;
	    int * split = new int [nRowB];
	    int * best = new int[nRowB];
	    int * current = new int [nRowB];
	    int * size = new int [nRowB];
	    {
	      memset(size,0,nRowB*sizeof(int));
	      for (i=0;i<nColB;i++) {
		int j;
		for (j=startCB[i];j<startCB[i+1];j++) {
		  int iRow = row2B[j];
		  size[iRow]++;
		}
	      }
	      for (i=0;i<nRowB;i++)
		if (size[i]<2)
		  printf("%d entries in row %d\n",size[i],i);
	    }
	    for (i=0;i<nColB;i++)
	      whichCol[i]=i;
	    for (i=0;i<nRowB;i++)
	      whichRow[i]=0;
	    int nLeft=nColB;
	    int nSet=1;
	    size[0]=nRowB;
	    while (nLeft) {
	      // find best column
	      int iBest=-1;
	      memset(best,0,nSet*sizeof(int));
	      memset(current,0,nSet*sizeof(int));
	      for (i=0;i<nColB;i++) {
		if (whichCol[i]<nLeft) {
		  int j;
		  for (j=startCB[i];j<startCB[i+1];j++) {
		    int iRow = row2B[j];
		    int iSet=whichRow[iRow];
		    current[iSet]++;
		  }
		  // See if better - could this be done faster
		  bool better=false;
		  for (j=nSet-1;j>=0;j--) {
		    if (current[j]>best[j]) {
		      better=true;
		      break;
		    } else if (current[j]<best[j]) {
		      break;
		    }
		  }
		  if (better) {
		    iBest=i;
		    memcpy(best,current,nSet*sizeof(int));
		  }
		  for (j=startCB[i];j<startCB[i+1];j++) {
		    int iRow = row2B[j];
		    int iSet=whichRow[iRow];
		    current[iSet]=0;
		  }
		}
	      }
	      assert (iBest>=0);
	      // swap
	      for (i=0;i<nColB;i++) {
		if (whichCol[i]==nLeft-1) {
		  whichCol[i]=whichCol[iBest];
		  whichCol[iBest]=nLeft-1;
		  break;
		}
	      }
	      // See which ones will have to split
	      int nMore=0;
	      for (i=0;i<nSet;i++) {
		current[i]=i+nMore;
		if (best[i]>0&&best[i]<size[i]) {
		  split[i]=i+nMore;
		  nMore++;
		} else {
		  split[i]=-1;
		}
	      }
	      if (nMore) {
		int j;
		for (j=startCB[iBest];j<startCB[iBest+1];j++) {
		  int iRow = row2B[j];
		  int iSet=whichRow[iRow];
		  int newSet = split[iSet];
		  if (newSet>=0) {
		    whichRow[iRow]=newSet+1+nRowB;
		  }
		}
		nSet += nMore;
		memset(size,0,nSet*sizeof(int));
		for (i=0;i<nRowB;i++) {
		  int iSet = whichRow[i];
		  if (iSet>=nRowB) {
		    // has 1 - correct it
		    iSet -= nRowB;
		  } else {
		    // 0 part of split set or not split
		    iSet=current[iSet];
		  }
		  whichRow[i]=iSet;
		  size[iSet]++;
		}
	      }
	      nLeft--;
	    }
	    if (nSet<nRowB) {
	      // ties - need to spread out whichRow
	      memset(split,0,nRowB*sizeof(int));
	      for (i=0;i<nRowB;i++) {
		int iSet = whichRow[i];
		split[iSet]++;
	      }
	      current[0]=0;
	      for (i=0;i<nSet;i++) {
		current[i+1]=current[i]+split[i];
		split[i]=current[i];
	      }
	      for (i=0;i<nRowB;i++) {
		int iSet = whichRow[i];
		int k = split[iSet];
		split[iSet]=k;
		whichRow[i]=k;
	      }
	    }
	    // Get inverse of whichCol
	    for (i=0;i<nColB;i++) {
	      int iColumn = whichCol[i];
	      startCC[iColumn]=i;
	    }
	    memcpy(whichCol,startCC,nColB*sizeof(int));
	    // Permute matrix
	    startCC[0]=0;
	    int nelB=0;
	    memset(split,0,nRowB*sizeof(int));
	    for (i=0;i<nColB;i++) {
	      int iColumn = whichCol[i];
	      int j;
	      for (j=startCB[iColumn];j<startCB[iColumn+1];j++) {
		int iRow = row2B[j];
		int iSet=whichRow[iRow];
		row2C[nelB++]=iSet;
		split[iSet]++;
	      }
	      startCC[i+1]=nelB;
	    }
	    startRowC[0]=0;
	    for (i=0;i<nRowB;i++) {
	      startRowC[i+1]=startRowC[i]+split[i];
	      split[i]=0;
	    }
	    for (i=0;i<nColB;i++) {
	      int j;
	      for (j=startCC[i];j<startCC[i+1];j++) {
		int iRow = row2C[j];
		int k=split[iRow]+startRowC[iRow];
		split[iRow]++;
		column2C[k]=i;
	      }
	    }
	    for (i=0;i<nRowB;i++)
	      split[i]=0;
	    goodInteger=true;
	    for (i=nColB-1;i>0;i--) {
	      int j;
	      for (j=startCC[i];j<startCC[i+1];j++) {
		int iRow = row2C[j];
		split[iRow]=1;
	      }
	      for (j=startCC[i];j<startCC[i+1];j++) {
		int iRow = row2C[j];
		for (int k=startRowC[iRow];k<startRowC[iRow+1];k++) {
		  int iColumn = column2C[k];
		  if (iColumn<i) {
		    for (int jj=startCC[iColumn];jj<startCC[iColumn+1];jj++) {
		      int jRow = row2C[jj];
		      if (jRow>iRow&&!split[jRow]) {
			// bad
			goodInteger=false;
			break;
		      }
		    }
		  }
		}
	      }
	      if (!goodInteger)
		break;
	      for (j=startCC[i];j<startCC[i+1];j++) {
		int iRow = row2C[j];
		split[iRow]=0;
	      }
	    }
	    delete [] split;
	    delete [] best;
	    delete [] current;
	    delete [] size;
	  } else {
	    // was network
	    goodInteger=true;
	  }
	  if (goodInteger) {
	    printf("Block %d can be integer\n",iBlock);
	    for (i=0;i<nCol;i++) {
	      if (columnBlock[i]==iBlock) {
		int iBack = back[i];
		model->setInteger(iBack);
	      }
	    }
	  }
	  delete [] startRowC;
	  delete [] column2C;
	  delete [] startCB;
	  delete [] row2B;
	  delete [] startCC;
	  delete [] row2C;
	  delete [] whichRow;
	  delete [] whichCol;
	}
	delete [] startRow;
	delete [] column2;
	delete [] element2;
	delete [] startC;
	delete [] row2;
	delete [] back;
      }
    }
    numberIntegers=numberNonZero;
    if (allGood&&numberObj) {
#ifdef COIN_DEVELOP
      int numberPossible = 0;
#endif
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if (upper[iColumn]>lower[iColumn]&&objective[iColumn]&&!model->isInteger(iColumn)) {
#ifdef COIN_DEVELOP
	  numberPossible++;
#endif
	  if (upper[iColumn]<=lower[iColumn]+10) {
	    model->setInteger(iColumn);
	    numberIntegers++;
	  }
	}
      }
#ifdef COIN_DEVELOP
      printf("ZZZZYY CglPreProcess analysis says all (%d) continuous with costs could be made integer - %d were\n",numberPossible,numberIntegers-numberNonZero);
#endif
    }
#ifdef COIN_DEVELOP
    if (numberZero)
      printf("ZZZZYY %d continuous with zero cost were made integer\n",numberZero);
#endif
    numberIntegers += numberZero;
#ifdef COIN_DEVELOP
    if (numberEq||numberEqI)
      printf("ZZZZYY %d rows made equality from continuous, %d from integer\n",numberEq,numberEqI);
#endif
    totalNumberIntegers += numberIntegers;
    if (!makeAll)
      numberIntegers=0;
  }
  delete [] rhs;
  delete [] count;
  return (totalNumberIntegers);
}
//#define CGL_WRITEMPS 
#ifdef CGL_WRITEMPS
extern double * debugSolution;
extern int debugNumberColumns;
static int mpsNumber=0;
static void writeDebugMps(const OsiSolverInterface * solver,
			  const char * where,
			  OsiPresolve * pinfo)
{ 
  mpsNumber++;
  char name[20];
  sprintf(name,"presolve%2.2d.mps",mpsNumber);
  printf("saving %s from %s\n",name,where);
  solver->writeMpsNative(name,NULL,NULL,0,1,0);
  if (pinfo&&debugSolution) {
    int n = solver->getNumCols();
    if (n<debugNumberColumns) {
      const int * original = pinfo->originalColumns();
      if (!original) {
	printf("No original columns\n");
	abort();
      }
      for (int i=0;i<n;i++) 
	debugSolution[i]=debugSolution[original[i]];
      debugNumberColumns=n;
    }
  }
  if (debugSolution) {
    OsiSolverInterface * newSolver = solver->clone();
    const double * lower = newSolver->getColLower();
    const double * upper = newSolver->getColUpper();
    for (int i = 0; i<debugNumberColumns;i++) {
      if (newSolver->isInteger(i)) {
	double value = floor(debugSolution[i]+0.5);
	if (value<lower[i]||value>upper[i]) {
	  printf("Bad value %d - %g %g %g\n",i,lower[i],debugSolution[i],
		 upper[i]);
	} else {
	  newSolver->setColLower(i,value);
	  newSolver->setColUpper(i,value);
	}
      }
    }
    printf("Starting solve %d\n",mpsNumber);
    newSolver->resolve();
    printf("Ending solve %d - status %s obj %g\n",mpsNumber,
	   newSolver->isProvenOptimal() ? "ok" : "bad",
	   newSolver->getObjValue());
    delete newSolver;
  }
}
#else
#define writeDebugMps(x,y,z)
#endif
OsiSolverInterface *
CglPreProcess::preProcessNonDefault(OsiSolverInterface & model, 
				    int makeEquality, int numberPasses,
				    int tuning)
{
#if 0
   bool rcdActive = true ;
   std::string modelName ;
   model.getStrParam(OsiProbName,modelName) ;
   std::cout
     << "  Attempting to activate row cut debugger for "
     << modelName << " ... " ;
   writeDebugMps(&model,"IPP:preProcessNonDefault",0) ;
   model.activateRowCutDebugger(modelName.c_str()) ;
   if (model.getRowCutDebugger())
     std::cout << "on optimal path." << std::endl ;
   else if (model.getRowCutDebuggerAlways())
     std::cout << "not on optimal path." << std::endl ;
   else {
     std::cout << "failure." << std::endl ;
     rcdActive = false ;
   }
   if (rcdActive) {
     const OsiRowCutDebugger *debugger = model.getRowCutDebuggerAlways() ;
     std::cout << "  Optimal solution is:" << std::endl ;
     debugger->printOptimalSolution(model) ;
   }
# endif
  originalModel_ = & model;
  numberSolvers_ = numberPasses;
  model_ = new OsiSolverInterface * [numberSolvers_];
  modifiedModel_ = new OsiSolverInterface * [numberSolvers_];
  presolve_ = new OsiPresolve * [numberSolvers_];
  for (int i=0;i<numberSolvers_;i++) {
    model_[i]=NULL;
    modifiedModel_[i]=NULL;
    presolve_[i]=NULL;
  }
  // Put presolve option on
  tuning |= 8;
  // clear original
  delete [] originalColumn_;
  delete [] originalRow_;
  originalColumn_=NULL;
  originalRow_=NULL;
  //startModel_=&model;
  // make clone
  delete startModel_;
  startModel_ = originalModel_->clone();
  CoinPackedMatrix matrixByRow(*originalModel_->getMatrixByRow());
  int numberRows = originalModel_->getNumRows();
  if (rowType_)
    assert (numberRowType_==numberRows);
  int numberColumns = originalModel_->getNumCols();
  //int originalNumberColumns=numberColumns;
  int minimumLength = 5;
  int numberModifiedPasses=10;
  if (numberPasses<=1)
    numberModifiedPasses=1; // lightweight preprocessing
  else if (numberPasses<=2)
    numberModifiedPasses=2; // fairly lightweight preprocessing
  if (tuning>=10000) {
    numberModifiedPasses=(tuning-10000)/10000;
    tuning %= 10000;
    //minimumLength = tuning;
  }
  //bool heavyProbing = (tuning&1)!=0;
  int makeIntegers = (tuning&6)>>1;
  // See if we want to do initial presolve
  int doInitialPresolve = (tuning&8)>>3;
  if (numberSolvers_<2)
    doInitialPresolve=0;
  // We want to add columns
  int numberSlacks=0;
  int * rows = new int[numberRows];
  double * element =new double[numberRows];
  
  int iRow;
  
  int numberCliques=0;
  int * which = new int[numberColumns];

  // Statistics
  int totalP1=0,totalM1=0;
  int numberFixed=0;
  // May just find it is infeasible
  bool feasible=true;
  
  // Row copy
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();
  
  const double * lower = originalModel_->getColLower();
  const double * upper = originalModel_->getColUpper();
  const double * rowLower = originalModel_->getRowLower();
  const double * rowUpper = originalModel_->getRowUpper();
  // Clean bounds
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (originalModel_->isInteger(iColumn)) {
      double lo = CoinMax(lower[iColumn],ceil(lower[iColumn]-1.0e-6));
      if (lo>lower[iColumn])
	originalModel_->setColLower(iColumn,lo);
      double up = CoinMin(upper[iColumn],floor(upper[iColumn]+1.0e-6));
      if (up<upper[iColumn])
	originalModel_->setColUpper(iColumn,up);
      if (lo>up)
	feasible=false;
    }
  }
  bool allToGub = makeEquality==5;
  if (allToGub)
    makeEquality=3;
  // Initialize random seed
  CoinThreadRandom randomGenerator(987654321);
  bool justOnesWithObj=false;
  if (makeEquality==2||makeEquality==3||makeEquality==4) {
    int iRow, iColumn;
    int numberIntegers = 0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (originalModel_->isInteger(iColumn))
        numberIntegers++;
    }
    // Look for possible SOS
    int numberSOS=0;
    int * mark = new int[numberColumns];
    CoinFillN(mark,numberColumns,-1);
    int numberOverlap=0;
    int numberInSOS=0;
    // See if worthwhile creating accumulation variables
    int firstOther=numberRows;
    int * whichRow = new int[numberRows];
    for (iRow=0;iRow<numberRows;iRow++) {
      if (rowUpper[iRow]==1.0) {
        if (rowLength[iRow]<5)
          continue;
        bool goodRow=true;
	bool overlap=false;
        for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
          int iColumn = column[j];
          if (elementByRow[j]!=1.0||!originalModel_->isInteger(iColumn)||lower[iColumn]) {
            goodRow=false;
            break;
          }
          if (mark[iColumn]>=0) {
            overlap=true;
            numberOverlap++;
          }
        }
        if (goodRow) {
	  if (!overlap) {
	    // mark all
	    for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	      int iColumn = column[j];
	      mark[iColumn]=numberSOS;
	    }
	    numberSOS++;
	    numberInSOS += rowLength[iRow];
	  }
	  // May still be interesting even if overlap
	  if (rowLength[iRow]>=5) {
	    firstOther--;
	    whichRow[firstOther]=iRow;
	  }
        }
      }
    }
    if (makeEquality==2&&false) {
      if(numberOverlap||numberIntegers>numberInSOS+1) {
	// try just ones with costs
	CoinFillN(mark,numberColumns,-1);
	numberOverlap=0;
	numberInSOS=0;
	bool allCostsInSOS=true;
	const double *objective = originalModel_->getObjCoefficients() ;
	for (iRow=0;iRow<numberRows;iRow++) {
	  if (rowUpper[iRow]==1.0&&rowLength[iRow]>=5) {
	    bool goodRow=true;
	    bool overlap=false;
	    int nObj=0;
	    for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	      int iColumn = column[j];
	      if (elementByRow[j]!=1.0||!originalModel_->isInteger(iColumn)||lower[iColumn]) {
		goodRow=false;
	      }
	      if (objective[iColumn])
		nObj++;
	      if (mark[iColumn]>=0) {
		overlap=true;
		numberOverlap++;
	      }
	    }
	    if (nObj&&nObj>=rowLength[iRow]-1) {
	      if (goodRow) {
		if (!overlap) {
		  // mark all
		  for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
		    int iColumn = column[j];
		    mark[iColumn]=numberSOS;
		  }
		  numberSOS++;
		  numberInSOS += rowLength[iRow];
		}
	      } else {
		// no good
		allCostsInSOS=false;
	      }
	    }
	  }
	}
	if (numberInSOS&&allCostsInSOS) {
	  int nGoodObj=0;
	  int nBadObj=0;
	  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	    if (objective[iColumn]) {
	      if (mark[iColumn]>=0)
		nGoodObj++;
	      else
		nBadObj++;
	    }
	  }
	  if (nBadObj*10<nGoodObj) {
	    justOnesWithObj=true;
	    makeEquality=3;
#ifdef CLP_INVESTIGATE
	    printf("trying SOS as all costs there\n");
#endif
	  }
	}
      }
    }
    if (firstOther<numberRows&&makeEquality==4) {
      CoinPackedMatrix * matrixByColumn = const_cast<CoinPackedMatrix *>(startModel_->getMatrixByCol());
      // Column copy
      const int * row = matrixByColumn->getIndices();
      const CoinBigIndex * columnStart = matrixByColumn->getVectorStarts();
      const int * columnLength = matrixByColumn->getVectorLengths(); 
      double * columnElements = matrixByColumn->getMutableElements();
      int * rowCount = new int[numberRows];
      memset(rowCount,0,numberRows*sizeof(int));
      double * rowValue = new double [numberRows];
      int numberY=0;
      int numberElements=0;
      int numberSOS=0;
      for (int kRow=firstOther;kRow<numberRows;kRow++) {
	int iRow=whichRow[kRow];
	int n=0;
	int j;
	for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	  int iColumn = column[j];
	  for (int k=columnStart[iColumn];k<columnStart[iColumn]+columnLength[iColumn];
	       k++) {
	    int jRow = row[k];
	    double value = columnElements[k];
	    if (jRow!=iRow) {
	      if (rowCount[jRow]>0) {
		if (value!=rowValue[jRow])
		  rowCount[jRow]=-1; // no good
		else
		  rowCount[jRow]++;
	      } else if (!rowCount[jRow]) {
		whichRow[n++]=jRow;
		rowCount[jRow]=1;
		rowValue[jRow]=value;
	      }
	    }
	  }
	}
	int bestRow=-1;
	int bestCount=4;
	for (j=0;j<n;j++) {
	  int jRow = whichRow[j];
	  int count=rowCount[jRow];
	  rowCount[jRow]=0;
	  if (count>=5) {
	    numberY++;
	    numberElements+=count;
	  }
	  if (count>bestCount) {
	    // possible
	    bestRow=jRow;
	    bestCount=count;
	  }
	}
	if (bestRow>=0) {
	  numberSOS++;
	  numberY++;
	  numberElements+=bestCount;
	}
      }
      if (numberY) {
	// Some may be duplicates
	// make sure ordered
	matrixByRow.orderMatrix();
	elementByRow = matrixByRow.getElements();
	column = matrixByRow.getIndices();
	rowStart = matrixByRow.getVectorStarts();
	rowLength = matrixByRow.getVectorLengths();
	CoinBigIndex * newStart = new CoinBigIndex[numberY+1];
	int * newColumn = new int [numberElements];
	double * newValue = new double [numberElements];
	double * hash = new double [numberY];
	double * hashColumn = new double [numberColumns];
	int i;
	for (i=0;i<numberColumns;i++)
	  hashColumn[i]=randomGenerator.randomDouble();
	double * valueY = new double [3*numberY];
	int * rowY = new int [3*numberY];
	int * columnY = new int[3*numberY];
	// For new solution
	double * newSolution = new double [numberColumns+numberY];
	memcpy(newSolution,startModel_->getColSolution(),numberColumns*sizeof(double));
	memset(rowCount,0,numberRows*sizeof(int));
	// List of SOS entries to zero out
	CoinBigIndex * where = new CoinBigIndex[numberColumns];
	numberY=0;
	numberElements=0;
	int numberElementsY=0;
	newStart[0]=0;
	for (int kRow=firstOther;kRow<numberRows;kRow++) {
	  int iRow=whichRow[kRow];
	  int n=0;
	  int j;
	  int saveNumberY=numberY;
	  for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	    int iColumn = column[j];
	    for (int k=columnStart[iColumn];k<columnStart[iColumn]+columnLength[iColumn];
		 k++) {
	      int jRow = row[k];
	      double value = columnElements[k];
	      if (jRow!=iRow) {
		if (rowCount[jRow]>0) {
		  if (value!=rowValue[jRow])
		    rowCount[jRow]=-1; // no good
		  else
		    rowCount[jRow]++;
		} else if (!rowCount[jRow]) {
		  whichRow[n++]=jRow;
		  rowCount[jRow]=1;
		  rowValue[jRow]=value;
		  assert (value);
		}
	      }
	    }
	  }
	  for (i=0;i<n;i++) {
	    // Sort so fewest first
	    std::sort(whichRow,whichRow+n);
	    int jRow = whichRow[i];
	    int count=rowCount[jRow];
	    rowCount[jRow]=0;
	    if (count>=5) {
	      //assert (count<rowLength[jRow]); // not error - just need to think
	      // mark so not looked at again
	      rowCount[jRow]=-count;
	      // form new row
	      double value=0.0;
	      double hashValue=0.0;
	      int nInSOS=0;
	      double valueOfY=0.0;
	      for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
		int iColumn = column[j];
		for (int k=columnStart[iColumn];k<columnStart[iColumn]+columnLength[iColumn];
		     k++) {
		  if (row[k]==jRow) {
		    value = columnElements[k];
		    newColumn[numberElements]=iColumn;
		    newValue[numberElements++]=1.0;
		    hashValue += hashColumn[iColumn];
		    columnElements[k]=0.0;
		    valueOfY += newSolution[iColumn];
		  } else if (row[k]==iRow) {
		    if (columnElements[k])
		      where[nInSOS++]=k;
		  }
		}
	      }
	      // See if already exists
	      int n=numberElements-newStart[numberY];
	      for (j=0;j<numberY;j++) {
		if (hashValue==hash[j]) {
		  // Double check
		  int offset=newStart[numberY]-newStart[j];
		  if (n==newStart[j+1]-newStart[j]) {
		    int k;
		    for (k=newStart[j];k<newStart[j]+n;k++) {
		      if (newColumn[k]!=newColumn[k+offset])
			break;
		    }
		    if (k==newStart[j+1])
		      break;
		  }
		}
	      }
	      if (j==numberY) {
		// not duplicate
		newSolution[numberY+numberColumns]=valueOfY;
		numberY++;
		newStart[numberY]=numberElements;
		hash[j]=hashValue;
		// Now do -1
		rowY[numberElementsY]=j+numberRows;
		columnY[numberElementsY]=j;
		valueY[numberElementsY++]=-1;
		if (n==nInSOS) {
		  // SOS entry
		  rowY[numberElementsY]=iRow;
		  columnY[numberElementsY]=j;
		  valueY[numberElementsY++]=1;
		  for (int i=0;i<n;i++) {
		    int iEl = where[i];
		    columnElements[iEl]=0.0;
		  }
		}
	      } else {
		// duplicate
		numberElements=newStart[numberY];
	      }
	      // Now do 
	      rowY[numberElementsY]=jRow;
	      columnY[numberElementsY]=j;
	      valueY[numberElementsY++]=value;
	    }
	  }
	  if (numberY>saveNumberY)
	    rowCount[iRow]=-1000;
	}
	delete [] hash;
	delete [] hashColumn;
	matrixByColumn->cleanMatrix();
	// Now add rows
	double * rhs = new double[numberY];
	memset(rhs,0,numberY*sizeof(double));
	startModel_->addRows(numberY,newStart,newColumn,newValue,rhs,rhs);
	delete [] rhs;
	delete [] newStart;
	delete [] newColumn;
	delete [] newValue;
	delete [] where;
	// Redo matrix
	CoinPackedMatrix add(true,rowY,columnY,valueY,numberElementsY);
	delete [] valueY;
	delete [] rowY;
	delete [] columnY;
	const int * row = add.getIndices();
	const CoinBigIndex * columnStart = add.getVectorStarts();
	//const int * columnLength = add.getVectorLengths(); 
	double * columnElements = add.getMutableElements();
	double * lo = new double [numberY];
	double * up = new double [numberY];
	for (i=0;i<numberY;i++) {
	  lo[i]=0.0;
	  up[i]=1.0;
	}
	startModel_->addCols(numberY,columnStart,row,columnElements,lo,up,NULL);
	delete [] lo;
	delete [] up;
	for (i=0;i<numberY;i++) 
	  startModel_->setInteger(i+numberColumns);
	CoinWarmStartBasis* basis =
	  dynamic_cast <CoinWarmStartBasis*>(startModel_->getWarmStart()) ;
	if (basis) {
	  for (i=0;i<numberY;i++) {
	    basis->setArtifStatus(i+numberRows,CoinWarmStartBasis::atLowerBound);
	    basis->setStructStatus(i+numberColumns,CoinWarmStartBasis::basic);
	  }
	  startModel_->setWarmStart(basis);
	  delete basis;
	}
	startModel_->setColSolution(newSolution);
	delete [] newSolution;
	writeDebugMps(startModel_,"start",NULL);
	if (numberElements<10*CoinMin(numberColumns,100*numberY)) {
	  handler_->message(CGL_ADDED_INTEGERS,messages_)
	    <<numberY<<numberSOS<<numberElements
	    <<CoinMessageEol;
	  numberColumns += numberY;
	  bool saveTakeHint;
	  OsiHintStrength saveStrength;
	  startModel_->getHintParam(OsiDoDualInResolve,
			      saveTakeHint,saveStrength);
	  startModel_->setHintParam(OsiDoDualInResolve,false,OsiHintTry);
	  startModel_->resolve();
	  numberIterationsPre_ += startModel_->getIterationCount();
	  startModel_->setHintParam(OsiDoDualInResolve,saveTakeHint,saveStrength);
	} else {
	  // not such a good idea?
	  delete startModel_;
	  startModel_=NULL;
	}
      }
      delete [] rowValue;
      delete [] rowCount;
    }
    if (makeEquality==4) {
      makeEquality=0;
#if 1
      // Try and make continuous variables integer
      // make clone
      if (!startModel_)
	startModel_ = originalModel_->clone();
      makeInteger();
#endif
    }
    delete [] whichRow;
    delete [] mark;
    if (numberSOS) {
      if (makeEquality==2) {
	if(numberOverlap||numberIntegers>numberInSOS+1) {
	  handler_->message(CGL_PROCESS_SOS2,messages_)
	    <<numberSOS<<numberInSOS<<numberIntegers<<numberOverlap
	    <<CoinMessageEol;
	  makeEquality=0;
	}
      }
    } else {
      // no sos
      makeEquality=0;
    }
  }
  if (startModel_) {
    lower = originalModel_->getColLower();
    upper = originalModel_->getColUpper();
    rowLower = originalModel_->getRowLower();
    rowUpper = originalModel_->getRowUpper();
  }
  // See if all + 1
  bool allPlusOnes=true;
  int nPossible=0;
  int numberMadeEquality=0;
  for (iRow=0;iRow<numberRows;iRow++) {
    int numberP1=0, numberM1=0;
    int numberTotal=0;
    int j;
    double upperValue=rowUpper[iRow];
    double lowerValue=rowLower[iRow];
    bool good=true;
    bool possibleSlack=true;
    bool allPlus=true;
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn = column[j];
      double value = elementByRow[j];
      if (upper[iColumn]-lower[iColumn]<1.0e-8) {
        // fixed
        upperValue -= lower[iColumn]*value;
        lowerValue -= lower[iColumn]*value;
        continue;
      } else if (!originalModel_->isBinary(iColumn)) {
        good = false;
	possibleSlack=false;
        //break;
      } else {
	numberTotal++;
      }
      
      if (fabs(value-floor(value+0.5))>1.0e-12)
	possibleSlack=false;;
      if (fabs(value)!=1.0) {
        good=false;
        allPlus=false;
      } else if (value>0.0) {
        which[numberP1++]=iColumn;
      } else {
        numberM1++;
        which[numberColumns-numberM1]=iColumn;
        allPlus=false;
      }
    }
    if (possibleSlack) {
      if(upperValue>1.0e20&&lowerValue>-1.0e12) {
	possibleSlack =  (fabs(lowerValue-floor(lowerValue+0.5))<1.0e-12);
      } else if(lowerValue<-1.0e20&&upperValue<1.0e12) {
	possibleSlack =  (fabs(upperValue-floor(upperValue+0.5))<1.0e-12);
      } else {
	possibleSlack=false;
      }
    }
    if (allPlus)
      nPossible++;
    int iUpper = static_cast<int> (floor(upperValue+1.0e-5));
    int iLower = static_cast<int> (ceil(lowerValue-1.0e-5));
    int state=0;
    if (upperValue<1.0e6) {
      if (iUpper==1-numberM1)
        state=1;
      else if (iUpper==-numberM1)
        state=2;
      else if (iUpper<-numberM1)
        state=3;
      if (fabs((static_cast<double> (iUpper))-upperValue)>1.0e-9)
        state =-1;
    }
    if (!state&&lowerValue>-1.0e6) {
      if (-iLower==1-numberP1)
        state=-1;
      else if (-iLower==-numberP1)
        state=-2;
      else if (-iLower<-numberP1)
        state=-3;
      if (fabs((static_cast<double> (iLower))-lowerValue)>1.0e-9)
        state =-1;
    }
    if (good&&state>0) {
      if (abs(state)==3) {
        // infeasible
        feasible=false;
        break;
      } else if (abs(state)==2) {
        // we can fix all
        numberFixed += numberP1+numberM1;
        int i;
        if (state>0) {
          // fix all +1 at 0, -1 at 1
          for (i=0;i<numberP1;i++)
            originalModel_->setColUpper(which[i],0.0);
          for (i=0;i<numberM1;i++)
            originalModel_->setColLower(which[numberColumns-i-1],1.0);
        } else {
          // fix all +1 at 1, -1 at 0
          for (i=0;i<numberP1;i++)
            originalModel_->setColLower(which[i],1.0);
          for (i=0;i<numberM1;i++)
            originalModel_->setColUpper(which[numberColumns-i-1],0.0);
        }
      } else {
        if (!makeEquality||(makeEquality==-1&&numberM1+numberP1<minimumLength))
          continue;
        if (makeEquality==2||makeEquality==3) {
          if (numberM1||numberP1<minimumLength) 
            continue;
        }
        numberCliques++;
        if (iLower!=iUpper) {
          element[numberSlacks]=state;
          rows[numberSlacks++]=iRow;
        }
        if (state>0) {
          totalP1 += numberP1;
          totalM1 += numberM1;
        } else {
          totalP1 += numberM1;
          totalM1 += numberP1;
        }
      }
    }
    if (possibleSlack&&makeEquality==-2&&(!good||state<=0)) {
      if (numberTotal<minimumLength)
	continue;
      numberMadeEquality++;
      element[numberSlacks] = (upperValue<1.0e10) ? 1.0 : -1.0; 
      rows[numberSlacks++]=iRow+numberRows;
    }
  }
  // allow if some +1's
  allPlusOnes = 10*nPossible>numberRows;
  delete [] which;
  if (!feasible) {
    handler_->message(CGL_INFEASIBLE,messages_)
      <<CoinMessageEol;
    delete [] rows;
    delete [] element;
    return NULL;
  } else {
    if (numberCliques) {
      handler_->message(CGL_CLIQUES,messages_)
        <<numberCliques
        << (static_cast<double>(totalP1+totalM1))/
	(static_cast<double> (numberCliques))
        <<CoinMessageEol;
      //printf("%d of these could be converted to equality constraints\n",
      //     numberSlacks);
    }
    if (numberFixed)
      handler_->message(CGL_FIXED,messages_)
        <<numberFixed
        <<CoinMessageEol;
  }
  if (numberSlacks&&makeEquality&&!justOnesWithObj) {
    handler_->message(CGL_SLACKS,messages_)
      <<numberSlacks
      <<CoinMessageEol;
    // add variables to make equality rows
    // Get new model
    if (!startModel_) {
      assert (!startModel_);
      startModel_ = originalModel_->clone();
    }
    for (int i=0;i<numberSlacks;i++) {
      int iRow = rows[i];
      double value = element[i];
      double lowerValue = 0.0;
      double upperValue = 1.0;
      double objValue  = 0.0;
      if (iRow>=numberRows) {
	// just a slack not a clique
	upperValue=COIN_DBL_MAX;
	iRow -= numberRows;
      }
      CoinPackedVector column(1,&iRow,&value);
      startModel_->addCol(column,lowerValue,upperValue,objValue);
      // set integer
      startModel_->setInteger(numberColumns+i);
      if (value >0)
	startModel_->setRowLower(iRow,rowUpper[iRow]);
      else
	startModel_->setRowUpper(iRow,rowLower[iRow]);
    }
  } else if (!startModel_) {
    // make clone anyway so can tighten bounds
    startModel_ = originalModel_->clone();
  }
  // move objective to integers or to aggregated
  lower = startModel_->getColLower();
  upper = startModel_->getColUpper();
  rowLower = startModel_->getRowLower();
  rowUpper = startModel_->getRowUpper();
  matrixByRow = CoinPackedMatrix(*startModel_->getMatrixByRow());
  elementByRow = matrixByRow.getElements();
  column = matrixByRow.getIndices();
  rowStart = matrixByRow.getVectorStarts();
  rowLength = matrixByRow.getVectorLengths();
  char * marked = new char [numberColumns];
  memset(marked,0,numberColumns);
  numberRows=startModel_->getNumRows();
  numberColumns=startModel_->getNumCols();
  //CoinPackedMatrix * matrixByColumn = const_cast<CoinPackedMatrix *>(startModel_->getMatrixByCol());
  // Column copy
  //const int * row = matrixByColumn->getIndices();
  //const CoinBigIndex * columnStart = matrixByColumn->getVectorStarts();
  //const int * columnLength = startModel_->getMatrixByCol()->getVectorLengths(); 
  //const double * columnElements = matrixByColumn->getElements();
  double * obj = CoinCopyOfArray(startModel_->getObjCoefficients(),numberColumns);
  double offset;
  int numberMoved=0;
  startModel_->getDblParam(OsiObjOffset,offset);
  for (iRow=0;iRow<numberRows;iRow++) {
    //int slack = -1;
    int nPlus=0;
    int nMinus=0;
    int iPlus=-1;
    int iMinus=-1;
    double valuePlus=0;
    double valueMinus=0;
    //bool allInteger=true;
    double rhs = rowLower[iRow];
    if (rhs!=rowUpper[iRow])
      continue;
    //int multiple=0;
    //int iSlack=-1;
    int numberContinuous=0;
    for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn = column[j];
      double value = elementByRow[j];
      if (upper[iColumn]>lower[iColumn]) {
	if (startModel_->isInteger(iColumn)) {
#if 0
	  if (columnLength[iColumn]==1) {
	    if (value==1.0) {
	    }
	  }
	  if (value!=floor(value+0.5))
	    allInteger=false;
	  if (allInteger&&fabs(value)<1.0e8) {
	    if (!multiple)
	      multiple = static_cast<int> (fabs(value));
	    else if (multiple>0)
	      multiple = gcd(multiple,static_cast<int> (fabs(value)));
	  } else {
	    allInteger=false;
	  }
#endif
	} else {
	  numberContinuous++;
	}
	if (value>0.0) {
	  if (nPlus>0&&value!=valuePlus) {
	    nPlus = - numberColumns;
	  } else if (!nPlus) {
	    nPlus=1;
	    iPlus=iColumn;
	    valuePlus=value;
	  } else {
	    nPlus++;
	  }
	} else {
	  if (nMinus>0&&value!=valueMinus) {
	    nMinus = - numberColumns;
	  } else if (!nMinus) {
	    nMinus=1;
	    iMinus=iColumn;
	    valueMinus=value;
	  } else {
	    nMinus++;
	  }
	}
      } else {
	rhs -= lower[iColumn]*value;
      }
    }
    if (((nPlus==1&&startModel_->isInteger(iPlus)&&nMinus>0)||
	 (nMinus==1&&startModel_->isInteger(iMinus)&&nPlus>0))&&numberContinuous&&true) {
      int jColumn;
      double multiplier;
      if (nPlus==1) {
	jColumn = iPlus;
	multiplier = fabs(valuePlus/valueMinus);
	rhs /= -valueMinus;
      } else {
	jColumn = iMinus;
	multiplier = fabs(valueMinus/valuePlus);
	rhs /= valuePlus;
      }
      double smallestPos=COIN_DBL_MAX;
      double smallestNeg=-COIN_DBL_MAX;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	int iColumn = column[j];
	if (iColumn!=jColumn) {
	  double objValue = obj[iColumn];
	  if (upper[iColumn]>lower[iColumn]) {
	    if (objValue>=0.0)
	      smallestPos=CoinMin(smallestPos,objValue);
	    else
	      smallestNeg=CoinMax(smallestNeg,objValue);
	  }
	}
      }
      if (smallestPos>0.0) {
	double move=0.0;
	if(smallestNeg==-COIN_DBL_MAX)
	  move=smallestPos;
	else if (smallestPos==COIN_DBL_MAX)
	  move=smallestNeg;
	if (move) {
	  // can move objective
	  numberMoved++;
#ifdef COIN_DEVELOP
	  if (rhs)
	    printf("ZZZ on col %d move %g offset %g\n",
		   jColumn,move,move*rhs);
#endif
	  offset += move*rhs;
	  for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	    int iColumn = column[j];
	    if (iColumn!=jColumn) {
	      if (upper[iColumn]>lower[iColumn]) {
		obj[iColumn] -= move;
	      }
	    } else {
	      obj[jColumn] += move*multiplier;
	    }
	  }
	}
      }
    }
  }
#ifdef COIN_DEVELOP
  if (numberMoved)
    printf("ZZZ %d costs moved\n",numberMoved);
#endif
  startModel_->setDblParam(OsiObjOffset,offset);
  startModel_->setObjective(obj);
  delete [] obj;
  delete [] marked;
  delete [] rows;
  delete [] element;
  if (makeIntegers) {
    makeIntegers2(startModel_,makeIntegers);
  }
  int infeas=0;
  OsiSolverInterface * startModel2 = startModel_;
  // Do we want initial presolve
  if (doInitialPresolve) {
    assert (doInitialPresolve==1);
    OsiSolverInterface * presolvedModel;
    OsiSolverInterface * oldModel = startModel2;
    OsiPresolve * pinfo = new OsiPresolve();
    int presolveActions=0;
    // Allow dual stuff on integers
    // Allow stuff which may not unroll cleanly
    presolveActions=1+16;
    if ((tuning&32)!=0)
      presolveActions |= 32;
    // Do not allow all +1 to be tampered with
    //if (allPlusOnes)
    //presolveActions |= 2;
    // allow transfer of costs
    // presolveActions |= 4;
    // If trying for SOS don't allow some transfers
    if (makeEquality==2||makeEquality==3)
      presolveActions |= 8;
    pinfo->setPresolveActions(presolveActions);
    if (prohibited_)
      assert (numberProhibited_==oldModel->getNumCols());
    int saveLogLevel = oldModel->messageHandler()->logLevel();
    if (saveLogLevel==1)
      oldModel->messageHandler()->setLogLevel(0);
    std::string solverName;
    oldModel->getStrParam(OsiSolverName,solverName);
    // Extend if you want other solvers to keep solution
    bool keepSolution=solverName=="clp";
    presolvedModel = pinfo->presolvedModel(*oldModel,1.0e-7,true,5,prohibited_,keepSolution,rowType_);
    oldModel->messageHandler()->setLogLevel(saveLogLevel);
    if (presolvedModel) {
      presolvedModel->messageHandler()->setLogLevel(saveLogLevel);
      //presolvedModel->writeMps("new");
      writeDebugMps(presolvedModel,"ordinary",pinfo);
      // update prohibited and rowType
      update(pinfo,presolvedModel);
      if (!presolvedModel->getNumRows()) {
	doInitialPresolve=0;
	delete presolvedModel;
	delete pinfo;
      } else {
	model_[0]=presolvedModel;
	presolve_[0]=pinfo;
	modifiedModel_[0]=presolvedModel->clone();
	startModel2 = modifiedModel_[0];
      }
    } else {
      infeas=1;
      doInitialPresolve=0;
      delete presolvedModel;
      delete pinfo;
    }
  }
  // tighten bounds
/*

  Virtuous solvers may require a refresh via initialSolve if this
  call is ever changed to give a nonzero value to the (default) second
  parameter. Previous actions may have made significant changes to the
  constraint system. Safe as long as tightenPrimalBounds doesn't ask for
  the current solution.
*/
  if (!infeas&&true) {
    // may be better to just do at end
    writeDebugMps(startModel2,"before",NULL);
    infeas = tightenPrimalBounds(*startModel2);
    writeDebugMps(startModel2,"after",NULL);
  }
  if (infeas) {
    handler_->message(CGL_INFEASIBLE,messages_)
      <<CoinMessageEol;
    return NULL;
  }
  OsiSolverInterface * returnModel=NULL;
  int numberChanges;
  if ((tuning&128)!=0) {
    // take out cliques
    OsiSolverInterface * newSolver=cliqueIt(*startModel2,0.0001);
    if (newSolver) {
      if (startModel2 == modifiedModel_[0])
	modifiedModel_[0]=newSolver;
      delete startModel2;
      startModel2=newSolver;
      newSolver->initialSolve();
      assert (newSolver->isProvenOptimal());
      //printf("new size %d rows, %d columns\n",
      //     newSolver->getNumRows(),newSolver->getNumCols());
    }
  }
  {
    // Give a hint to do dual
    bool saveTakeHint;
    OsiHintStrength saveStrength;
    startModel2->getHintParam(OsiDoDualInInitial,
			      saveTakeHint,saveStrength);
    startModel2->setHintParam(OsiDoDualInInitial,true,OsiHintTry);
    startModel2->initialSolve();
    numberIterationsPre_ += startModel2->getIterationCount();
    // double check
    if (!startModel2->isProvenOptimal()) {
      if (!startModel2->isProvenDualInfeasible()) {
	// Do presolves
	bool saveHint;
	OsiHintStrength saveStrength;
	startModel2->getHintParam(OsiDoPresolveInInitial,saveHint,saveStrength);
	startModel2->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
	startModel2->setHintParam(OsiDoDualInInitial,false,OsiHintTry);
	startModel2->initialSolve();
	numberIterationsPre_ += startModel2->getIterationCount();
	if (!startModel2->isProvenDualInfeasible()) {
	  CoinWarmStart * empty = startModel2->getEmptyWarmStart();
	  startModel2->setWarmStart(empty);
	  delete empty;
	  startModel2->setHintParam(OsiDoDualInInitial,true,OsiHintTry);
	  startModel2->initialSolve();
	  numberIterationsPre_ += startModel2->getIterationCount();
	}
	startModel2->setHintParam(OsiDoPresolveInInitial,saveHint,saveStrength);
      }
    }
    startModel2->setHintParam(OsiDoDualInInitial,saveTakeHint,saveStrength);
  }
  if (!startModel2->isProvenOptimal()) {
    if (!startModel2->isProvenDualInfeasible()) {
      handler_->message(CGL_INFEASIBLE,messages_)<< CoinMessageEol ;
#ifdef COIN_DEVELOP
      startModel2->writeMps("infeas");
#endif
    } else {
      handler_->message(CGL_UNBOUNDED,messages_)<< CoinMessageEol ;
    }
    return NULL;
  }
  reducedCostFix(*startModel2);
  if (!numberSolvers_) {
    // just fix
    OsiSolverInterface * newModel = modified(startModel2,false,numberChanges,0,numberModifiedPasses);
    if (startModel_!=originalModel_)
      delete startModel_;
    if (startModel2!=startModel_)
      delete startModel2;
    startModel_=newModel;
    returnModel=startModel_;
  } else {
    OsiSolverInterface * presolvedModel;
    OsiSolverInterface * oldModel = startModel2;
    if (doInitialPresolve)
      oldModel = modifiedModel_[0];
    //CglDuplicateRow dupCuts(oldModel);
    //dupCuts.setLogLevel(1);
    // If +1 try duplicate rows
#define USECGLCLIQUE 512
    if ((options_&8)!=0)
      tuning &= ~USECGLCLIQUE;
    if ((options_&4)!=0)
      allPlusOnes=false;
    if (allPlusOnes||(tuning&USECGLCLIQUE)!=0) {
#if 1
      // put at beginning
      int nAdd= ((tuning&(64+USECGLCLIQUE))==64+USECGLCLIQUE&&allPlusOnes) ? 2 : 1;
      CglCutGenerator ** temp = generator_;
      generator_ = new CglCutGenerator * [numberCutGenerators_+nAdd];
      memcpy(generator_+nAdd,temp,numberCutGenerators_*sizeof(CglCutGenerator *));
      delete[] temp ;
      numberCutGenerators_+=nAdd;
      if (nAdd==2||(tuning&USECGLCLIQUE)!=0) {
	CglClique * cliqueGen=new CglClique(false,true);
	cliqueGen->setStarCliqueReport(false);
	cliqueGen->setRowCliqueReport(false);
	if ((tuning&USECGLCLIQUE)==0)
	  cliqueGen->setMinViolation(-2.0);
	else
	  cliqueGen->setMinViolation(-3.0);
	generator_[0]=cliqueGen;
      }
      if (allPlusOnes) {
	CglDuplicateRow * dupCuts =new CglDuplicateRow(oldModel);
	if ((tuning&256)!=0)
	  dupCuts->setMaximumDominated(numberColumns);
	generator_[nAdd-1]=dupCuts;
      }
#else
      CglDuplicateRow dupCuts(oldModel);
      addCutGenerator(&dupCuts);
#endif
    }
    for (int iPass=doInitialPresolve;iPass<numberSolvers_;iPass++) {
      // Look at Vubs
      {
        const double * columnLower = oldModel->getColLower();
        const double * columnUpper = oldModel->getColUpper();
        const CoinPackedMatrix * rowCopy = oldModel->getMatrixByRow();
        const int * column = rowCopy->getIndices();
        const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
        const int * rowLength = rowCopy->getVectorLengths(); 
        const double * rowElements = rowCopy->getElements();
        const CoinPackedMatrix * columnCopy = oldModel->getMatrixByCol();
        //const int * row = columnCopy->getIndices();
        //const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
        const int * columnLength = columnCopy->getVectorLengths(); 
        //const double * columnElements = columnCopy->getElements();
        const double * rowLower = oldModel->getRowLower();
        const double * rowUpper = oldModel->getRowUpper();
        const double * objective = oldModel->getObjCoefficients();
        double direction = oldModel->getObjSense();
        int numberRows = oldModel->getNumRows();
        for (int iRow=0;iRow<numberRows;iRow++) {
          if (rowLength[iRow]==2&&(rowLower[iRow]<-1.0e20||rowUpper[iRow]>1.0e20)) {
            CoinBigIndex start = rowStart[iRow];
            int iColumn1 = column[start];
            int iColumn2 = column[start+1];
            double value1 = rowElements[start];
            double value2 = rowElements[start+1];
            double upper;
            if (rowLower[iRow]<-1.0e20) {
              if (rowUpper[iRow]<1.0e20)
                upper = rowUpper[iRow];
              else
                continue; // free row
            } else {
              upper = - rowLower[iRow];
              value1=-value1;
              value2=-value2;
            }
            //for now just singletons
            bool integer1 = oldModel->isInteger(iColumn1);
            bool integer2 = oldModel->isInteger(iColumn2);
            int debug=0;
            if (columnLength[iColumn1]==1) {
              if (integer1) {
                debug=0;// no good
              } else if (integer2) {
                // possible
                debug=1;
              }
            } else if (columnLength[iColumn2]==1) {
              if (integer2) {
                debug=-1; // print and skip
              } else if (integer1) {
                // possible
                debug=1;
                double valueD = value1;
                value1 = value2;
                value2 = valueD;
                int valueI = iColumn1;
                iColumn1 = iColumn2;
                iColumn2 = valueI;
                bool valueB = integer1;
                integer1 = integer2;
                integer2 = valueB;
              }
            }
            if (debug&&0) {
              printf("%d %d elements%selement %g and %d %d elements%selement %g <= %g\n",
                     iColumn1,columnLength[iColumn1],integer1 ? " (integer) " : " ",value1,
                     iColumn2,columnLength[iColumn2],integer2 ? " (integer) " : " ",value2,
                     upper);
            }
            if (debug>0) {
              if (value1>0.0&&objective[iColumn1]*direction<0.0) {
                // will push as high as possible so make ==
                // highest effective rhs
                if (value2>0) 
                  upper -= value2 * columnLower[iColumn2];
                else
                  upper -= value2 * columnUpper[iColumn2];
                if (columnUpper[iColumn1]>1.0e20||
                    columnUpper[iColumn1]*value1>=upper) {
                  //printf("looks possible\n");
                  // make equality
                  if (rowLower[iRow]<-1.0e20) 
                    oldModel->setRowLower(iRow,rowUpper[iRow]);
                  else
                    oldModel->setRowUpper(iRow,rowLower[iRow]);
                } else {
                  // may be able to make integer
                  // may just be better to use to see objective integral
                  if (upper==floor(upper)&&value2==floor(value2)&&
                      value1==floor(value1)&&objective[iColumn1]==floor(objective[iColumn1]))
                    oldModel->setInteger(iColumn1);
                  //printf("odd3\n");
                }
              } else if (value1<0.0&&objective[iColumn1]*direction>0.0) {
                //printf("odd4\n");
              } else {
                //printf("odd2\n");
              }
            } else if (debug<0) {
              //printf("odd1\n");
            }
          }
        }
      }
      OsiPresolve * pinfo = new OsiPresolve();
      int presolveActions=0;
      // Allow dual stuff on integers
      // Allow stuff which may not unroll cleanly
      presolveActions=1+16;
      // Do not allow all +1 to be tampered with
      //if (allPlusOnes)
      //presolveActions |= 2;
      // allow transfer of costs
      // presolveActions |= 4;
      // If trying for SOS don't allow some transfers
      if (makeEquality==2||makeEquality==3)
        presolveActions |= 8;
      pinfo->setPresolveActions(presolveActions);
      if (prohibited_)
        assert (numberProhibited_==oldModel->getNumCols());
/*
  VIRTUOUS but possible bad for performance 
  
  At this point, the solution is most likely stale: we may have added cuts as
  we left the previous call to modified(), or we may have changed row bounds
  in VUB analysis just above. Continuous presolve doesn't need a solution
  unless we want it to transform the current solution to match the presolved
  model.
*/
      int saveLogLevel = oldModel->messageHandler()->logLevel();
      if (saveLogLevel==1)
	oldModel->messageHandler()->setLogLevel(0);
      std::string solverName;
      oldModel->getStrParam(OsiSolverName,solverName);
      // Extend if you want other solvers to keep solution
      bool keepSolution=solverName=="clp";
      presolvedModel = pinfo->presolvedModel(*oldModel,1.0e-7,true,5,
					     prohibited_,keepSolution,rowType_);
      oldModel->messageHandler()->setLogLevel(saveLogLevel);
      if (!presolvedModel) {
        returnModel=NULL;
	delete pinfo;
        break;
      }
      presolvedModel->messageHandler()->setLogLevel(saveLogLevel);
      // update prohibited and rowType
      update(pinfo,presolvedModel);
      writeDebugMps(presolvedModel,"ordinary2",pinfo);
      model_[iPass]=presolvedModel;
      presolve_[iPass]=pinfo;
      if (!presolvedModel->getNumRows()) {
        // was returnModel=oldModel;
        returnModel=presolvedModel;
        numberSolvers_=iPass+1;
        break; // model totally solved
      }
      bool constraints = iPass<numberPasses-1;
      // Give a hint to do primal
      bool saveTakeHint;
      OsiHintStrength saveStrength;
      presolvedModel->getHintParam(OsiDoDualInInitial,
                                   saveTakeHint,saveStrength);
      //if (iPass)
      presolvedModel->setHintParam(OsiDoDualInInitial,false,OsiHintTry);
      presolvedModel->initialSolve();
      numberIterationsPre_ += presolvedModel->getIterationCount();
      presolvedModel->setHintParam(OsiDoDualInInitial,saveTakeHint,saveStrength);
      if (!presolvedModel->isProvenOptimal()) {
	writeDebugMps(presolvedModel,"bad2",NULL);
	CoinWarmStartBasis *slack =
	  dynamic_cast<CoinWarmStartBasis *>(presolvedModel->getEmptyWarmStart()) ;
	presolvedModel->setWarmStart(slack);
	delete slack ;
	presolvedModel->resolve();
	if (!presolvedModel->isProvenOptimal()) {
	  returnModel=NULL;
	  //printf("infeasible\n");
	  break;
	} else {
	  //printf("feasible on second try\n");
	}
      }
      // maybe we can fix some
      int numberFixed = 
      reducedCostFix(*presolvedModel);
#ifdef COIN_DEVELOP
      if (numberFixed)
	printf("%d variables fixed on reduced cost\n",numberFixed);
#endif
      OsiSolverInterface * newModel = modified(presolvedModel,constraints,numberChanges,iPass-doInitialPresolve,numberModifiedPasses);
      returnModel=newModel;
      if (!newModel) {
        break;
      }
      modifiedModel_[iPass]=newModel;
      oldModel=newModel;
      writeDebugMps(newModel,"ordinary3",NULL);
      if (!numberChanges&&!numberFixed) {
#ifdef COIN_DEVELOP
	printf("exiting after pass %d of %d\n",iPass,numberSolvers_);
#endif
        numberSolvers_=iPass+1;
        break;
      }
    }
  }
  if (returnModel) {
    if (returnModel->getNumRows()) {
      // tighten bounds
      int infeas = tightenPrimalBounds(*returnModel);
      if (infeas) {
        delete returnModel;
	for (int iPass=0;iPass<numberSolvers_;iPass++) {
	  if (returnModel==modifiedModel_[iPass])
	    modifiedModel_[iPass]=NULL;
	}
	//printf("startModel_ %p startModel2 %p originalModel_ %p returnModel %p\n",
	//     startModel_,startModel2,originalModel_,returnModel);
        if (returnModel==startModel_&&startModel_!=originalModel_)
          startModel_=NULL;
        returnModel=NULL;
      }
    }
  } else {
    handler_->message(CGL_INFEASIBLE,messages_)
      <<CoinMessageEol;
  }
  int numberIntegers=0;
  if (returnModel) {
    int iColumn;
    int numberColumns = returnModel->getNumCols();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (returnModel->isInteger(iColumn))
        numberIntegers++;
    }
  }
  if ((makeEquality==2||makeEquality==3)&&numberCliques&&returnModel) {
    int iRow, iColumn;
    int numberColumns = returnModel->getNumCols();
    int numberRows = returnModel->getNumRows();
    const double * objective = returnModel->getObjCoefficients();
    // get row copy
    const CoinPackedMatrix * matrix = returnModel->getMatrixByRow();
    const double * element = matrix->getElements();
    const int * column = matrix->getIndices();
    const CoinBigIndex * rowStart = matrix->getVectorStarts();
    const int * rowLength = matrix->getVectorLengths();
    const double * rowLower = returnModel->getRowLower();
    const double * rowUpper = returnModel->getRowUpper();
    const double * columnLower = returnModel->getColLower();
    
    // Look for possible SOS
    int numberSOS=0;
    int * mark = new int[numberColumns];
    int * sosRow = new int [numberRows];
    CoinZeroN(sosRow,numberRows);
    CoinFillN(mark,numberColumns,-1);
    int numberOverlap=0;
    int numberInSOS=0;
    for (iRow=0;iRow<numberRows;iRow++) {
      if (rowLower[iRow]==1.0&&rowUpper[iRow]==1.0) {
        if ((rowLength[iRow]<5&&!justOnesWithObj)||(rowLength[iRow]<20&&allToGub))
          continue;
        bool goodRow=true;
	int nObj=0;
        for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
          iColumn = column[j];
          if (element[j]!=1.0||!returnModel->isInteger(iColumn)||columnLower[iColumn]) {
            goodRow=false;
            break;
          }
          if (mark[iColumn]>=0&&!allToGub) {
            goodRow=false;
            numberOverlap++;
          }
	  if (objective[iColumn])
	    nObj++;
        }
	if (goodRow&&justOnesWithObj) {
	  if (!nObj||nObj<rowLength[iRow]-1) 
	    goodRow=false;
	}
        if (goodRow) {
          // mark all
          for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
            int iColumn = column[j];
            mark[iColumn]=numberSOS;
          }
          sosRow[numberSOS++]=iRow;
          numberInSOS += rowLength[iRow];
        }
      }
    }
    if (numberSOS) {
      if (makeEquality==2&&(numberOverlap||numberIntegers>numberInSOS+1)) {
        handler_->message(CGL_PROCESS_SOS2,messages_)
          <<numberSOS<<numberInSOS<<numberIntegers<<numberOverlap
          <<CoinMessageEol;
      } else {
        handler_->message(CGL_PROCESS_SOS1,messages_)
          <<numberSOS<<numberInSOS
          <<CoinMessageEol;
        numberSOS_=numberSOS;
        typeSOS_ = new int[numberSOS_];
        startSOS_ = new int[numberSOS_+1];
        whichSOS_ = new int[numberInSOS];
        weightSOS_ = new double[numberInSOS];
        numberInSOS=0;
        startSOS_[0]=0;
        const CoinPackedMatrix * columnCopy = returnModel->getMatrixByCol();
        const int * row = columnCopy->getIndices();
        const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
        const int * columnLength = columnCopy->getVectorLengths(); 
        const double * columnElements = columnCopy->getElements();
        const double * objective = returnModel->getObjCoefficients();
        int * numberInRow = new int [numberRows];
        double * sort = new double[numberColumns];
        int * which = new int[numberColumns];
        for (int iSOS =0;iSOS<numberSOS_;iSOS++) {
          int n=0;
          int numberObj=0;
          CoinZeroN(numberInRow,numberRows);
	  int kRow = sosRow[iSOS];
          for (int j=rowStart[kRow];j<rowStart[kRow]+rowLength[kRow];j++) {
            int iColumn = column[j];
	    whichSOS_[numberInSOS]=iColumn;
	    weightSOS_[numberInSOS]=n;
	    numberInSOS++;
	    n++;
	    if (objective[iColumn])
	      numberObj++;
	    for (CoinBigIndex j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      if (!sosRow[iRow])
		numberInRow[iRow]++;
	    }
	  }
          // See if any rows look good
          int bestRow=-1;
          int numberDifferent=1;
          int start = startSOS_[iSOS];
          for (int iRow=0;iRow<numberRows;iRow++) {
            if (numberInRow[iRow]>=n-1) {
              // See how many different
              int i;
              for ( i=0;i<n;i++) {
                int iColumn = whichSOS_[i+start];
                sort[i]=0.0;
                which[i]=iColumn;
                for (CoinBigIndex j=columnStart[iColumn];
                     j<columnStart[iColumn]+columnLength[iColumn];j++) {
                  int jRow = row[j];
                  if (jRow==iRow) {
                    sort[i]=columnElements[j];
                    break;
                  }
                }
              }
              // sort
              CoinSort_2(sort,sort+n,which);
              double last = sort[0];
              int nDiff=1;
              for ( i=1;i<n;i++) {
                if (sort[i]>last+CoinMax(fabs(last)*1.0e-8,1.0e-5)) {
                  nDiff++;
                }
                last = sort[i];
              }
              if (nDiff>numberDifferent) {
                numberDifferent = nDiff;
                bestRow=iRow;
              }
            }
          }
          if (numberObj>=n-1||bestRow<0) {
            int i;
            for ( i=0;i<n;i++) {
              int iColumn = whichSOS_[i+start];
              sort[i]=objective[iColumn];
              which[i]=iColumn;
            }
            // sort
            CoinSort_2(sort,sort+n,which);
            double last = sort[0];
            int nDiff=1;
            for ( i=1;i<n;i++) {
              if (sort[i]>last+CoinMax(fabs(last)*1.0e-8,1.0e-5)) {
                nDiff++;
              }
              last = sort[i];
            }
            if (nDiff>numberDifferent) {
              numberDifferent = nDiff;
              bestRow=numberRows;
            }
          }
          if (bestRow>=0) {
            // if not objective - recreate
            if (bestRow<numberRows) {
              int i;
              for ( i=0;i<n;i++) {
                int iColumn = whichSOS_[i+start];
                sort[i]=0.0;
                which[i]=iColumn;
                for (CoinBigIndex j=columnStart[iColumn];
                     j<columnStart[iColumn]+columnLength[iColumn];j++) {
                  int jRow = row[j];
                  if (jRow==bestRow) {
                    sort[i]=columnElements[j];
                    break;
                  }
                }
              }
              // sort
              CoinSort_2(sort,sort+n,which);
            }
            // make sure gaps OK
            double last = sort[0];
            for (int i=1;i<n;i++) {
              double next = last+CoinMax(fabs(last)*1.0e-8,1.0e-5);
              sort[i]=CoinMax(sort[i],next);
              last = sort[i];
            }
            //CoinCopyN(sort,n,weightSOS_+start);
            //CoinCopyN(which,n,whichSOS_+start);
          }
          typeSOS_[iSOS]=1;
          startSOS_[iSOS+1]=numberInSOS;
        }
        delete [] numberInRow;
        delete [] sort;
        delete [] which;
      }
    }
    delete [] mark;
    delete [] sosRow;
  }
  if (returnModel) {
    if (makeIntegers) 
      makeIntegers2(returnModel,makeIntegers);
    handler_->message(CGL_PROCESS_STATS2,messages_)
      <<returnModel->getNumRows()<<returnModel->getNumCols()
      <<numberIntegers<<returnModel->getNumElements()
      <<CoinMessageEol;
    // If can make some cuts then do so
    if (rowType_) {
      int numberRows = returnModel->getNumRows();
      int numberCuts=0;
      for (int i=0;i<numberRows;i++) {
	if (rowType_[i]>0)
	  numberCuts++;
      }
      if (numberCuts) {
	CglStored stored;
	
	int * whichRow = new int[numberRows];
	// get row copy
	const CoinPackedMatrix * rowCopy = returnModel->getMatrixByRow();
	const int * column = rowCopy->getIndices();
	const int * rowLength = rowCopy->getVectorLengths();
	const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
	const double * rowLower = returnModel->getRowLower();
	const double * rowUpper = returnModel->getRowUpper();
	const double * element = rowCopy->getElements();
	int iRow,nDelete=0;
	for (iRow=0;iRow<numberRows;iRow++) {
	  if (rowType_[iRow]==1) {
	    // take out
	    whichRow[nDelete++]=iRow;
	  }
	}
	for (int jRow=0;jRow<nDelete;jRow++) {
	  iRow=whichRow[jRow];
	  int start = rowStart[iRow];
	  stored.addCut(rowLower[iRow],rowUpper[iRow],rowLength[iRow],
			column+start,element+start);
	}
	returnModel->deleteRows(nDelete,whichRow);
	delete [] whichRow;
	cuts_ = stored;
      }
    }
  }
#if 0
  if (returnModel) {
    int numberColumns = returnModel->getNumCols();
    int numberRows = returnModel->getNumRows();
    int * del = new int [CoinMax(numberColumns,numberRows)];
    int * original = new int [numberColumns];
    int nDel=0;
    for (int i=0;i<numberColumns;i++) {
      original[i]=i;
      if (returnModel->isInteger(i))
	del[nDel++]=i;
    }
    int nExtra=0;
    if (nDel&&nDel!=numberColumns&&(options_&1)!=0&&false) {
      OsiSolverInterface * yyyy = returnModel->clone();
      int nPass=0;
      while (nDel&&nPass<10) {
	nPass++;
	OsiSolverInterface * xxxx = yyyy->clone();
	int nLeft=0;
	for (int i=0;i<nDel;i++) 
	  original[del[i]]=-1;
	for (int i=0;i<numberColumns;i++) {
	  int kOrig=original[i];
	  if (kOrig>=0)
	    original[nLeft++]=kOrig;
	}
	assert (nLeft==numberColumns-nDel);
	xxxx->deleteCols(nDel,del);
	numberColumns = xxxx->getNumCols();
	const CoinPackedMatrix * rowCopy = xxxx->getMatrixByRow();
	numberRows = rowCopy->getNumRows();
	const int * column = rowCopy->getIndices();
	const int * rowLength = rowCopy->getVectorLengths();
	const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
	const double * rowLower = xxxx->getRowLower();
	const double * rowUpper = xxxx->getRowUpper();
	const double * element = rowCopy->getElements();
        const CoinPackedMatrix * columnCopy = xxxx->getMatrixByCol();
        const int * columnLength = columnCopy->getVectorLengths(); 
	nDel=0;
	// Could do gcd stuff on ones with costs
	for (int i=0;i<numberRows;i++) {
	  if (!rowLength[i]) {
	    del[nDel++]=i;
	  } else if (rowLength[i]==1) {
	    int k=rowStart[i];
	    int iColumn = column[k];
	    if (!xxxx->isInteger(iColumn)) {
	      double mult =1.0/fabs(element[k]);
	      if (rowLower[i]<-1.0e20) {
		double value = rowUpper[i]*mult;
		if (fabs(value-floor(value+0.5))<1.0e-8) {
		  del[nDel++]=i;
		  if (columnLength[iColumn]==1) {
		    xxxx->setInteger(iColumn);
		    int kOrig=original[iColumn];
		    returnModel->setInteger(kOrig);
		  }
		}
	      } else if (rowUpper[i]>1.0e20) {
		double value = rowLower[i]*mult;
		if (fabs(value-floor(value+0.5))<1.0e-8) {
		  del[nDel++]=i;
		  if (columnLength[iColumn]==1) {
		    xxxx->setInteger(iColumn);
		    int kOrig=original[iColumn];
		    returnModel->setInteger(kOrig);
		  }
		}
	      } else {
		double value = rowUpper[i]*mult;
		if (rowLower[i]==rowUpper[i]&&
		    fabs(value-floor(value+0.5))<1.0e-8) {
		  del[nDel++]=i;
		  xxxx->setInteger(iColumn);
		  int kOrig=original[iColumn];
		  returnModel->setInteger(kOrig);
		}
	      }
	    }
	  } else {
	    // only if all singletons
	    bool possible=false;
	    if (rowLower[i]<-1.0e20) {
	      double value = rowUpper[i];
	      if (fabs(value-floor(value+0.5))<1.0e-8) 
		possible=true;
	    } else if (rowUpper[i]>1.0e20) {
	      double value = rowLower[i];
	      if (fabs(value-floor(value+0.5))<1.0e-8) 
		possible=true;
	    } else {
	      double value = rowUpper[i];
	      if (rowLower[i]==rowUpper[i]&&
		  fabs(value-floor(value+0.5))<1.0e-8)
		possible=true;
	    }
	    if (possible) {
	      for (CoinBigIndex j=rowStart[i];
		   j<rowStart[i]+rowLength[i];j++) {
		int iColumn = column[j];
		if (columnLength[iColumn]!=1||fabs(element[j])!=1.0) {
		  possible=false;
		  break;
		}
	      }
	      if (possible) {
		for (CoinBigIndex j=rowStart[i];
		     j<rowStart[i]+rowLength[i];j++) {
		  int iColumn = column[j];
		  if (!xxxx->isInteger(iColumn)) {
		    xxxx->setInteger(iColumn);
		    int kOrig=original[iColumn];
		    returnModel->setInteger(kOrig);
		  }
		}
		del[nDel++]=i;
	      }
	    }
	  }
	}
	if (nDel) {
	  xxxx->deleteRows(nDel,del);
	}
	if (nDel!=numberRows) {
	  nDel=0;
	  for (int i=0;i<numberColumns;i++) {
	    if (xxxx->isInteger(i)) {
	      del[nDel++]=i;
	      nExtra++;
	    }
	  }
	} 
	delete yyyy;
	yyyy=xxxx->clone();
      }
      numberColumns = yyyy->getNumCols();
      numberRows = yyyy->getNumRows();
      if (!numberColumns||!numberRows) {
	printf("All gone\n");
	int numberColumns = returnModel->getNumCols();
	for (int i=0;i<numberColumns;i++)
	  assert(returnModel->isInteger(i));
      }
      // Would need to check if original bounds integer
      //yyyy->writeMps("noints");
      delete yyyy;
      printf("Creating simplified model with %d rows and %d columns - %d extra integers\n",
	     numberRows,numberColumns,nExtra);
    }
    delete [] del;
    delete [] original;
    //exit(2);
  }
#endif
  return returnModel;
}

/* Tightens primal bounds to make dual and branch and cutfaster.  Unless
   fixed, bounds are slightly looser than they could be.
   Returns non-zero if problem infeasible
   Fudge for branch and bound - put bounds on columns of factor *
   largest value (at continuous) - should improve stability
   in branch and bound on infeasible branches (0.0 is off)
*/
int 
CglPreProcess::tightenPrimalBounds(OsiSolverInterface & model,double factor)
{
  
  // Get a row copy in standard format
  CoinPackedMatrix copy = *model.getMatrixByRow();
  // get matrix data pointers
  const int * column = copy.getIndices();
  const CoinBigIndex * rowStart = copy.getVectorStarts();
  const int * rowLength = copy.getVectorLengths(); 
  double * element = copy.getMutableElements();
  int numberChanged=1,iPass=0;
  double large = model.getInfinity()*0.1; // treat bounds > this as infinite
  int numberInfeasible=0;
  int totalTightened = 0;

  double tolerance;
  model.getDblParam(OsiPrimalTolerance,tolerance);


  int numberColumns=model.getNumCols();
  const double * colLower = model.getColLower();
  const double * colUpper = model.getColUpper();
  // New and saved column bounds
  double * newLower = new double [numberColumns];
  memcpy(newLower,colLower,numberColumns*sizeof(double));
  double * newUpper = new double [numberColumns];
  memcpy(newUpper,colUpper,numberColumns*sizeof(double));
  double * columnLower = new double [numberColumns];
  memcpy(columnLower,colLower,numberColumns*sizeof(double));
  double * columnUpper = new double [numberColumns];
  memcpy(columnUpper,colUpper,numberColumns*sizeof(double));

  int iRow, iColumn;

  // If wanted - tighten column bounds using solution
  if (factor) {
    /*
      Callers need to ensure that the solution is fresh 
    */
    const double * solution =  model.getColSolution();
    double largest=0.0;
    if (factor>0.0) {
      assert (factor>1.0);
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
        if (columnUpper[iColumn]-columnLower[iColumn]>tolerance) {
          largest = CoinMax(largest,fabs(solution[iColumn]));
        }
      }
      largest *= factor;
    } else {
      // absolute
       largest = - factor;
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (columnUpper[iColumn]-columnLower[iColumn]>tolerance) {
        newUpper[iColumn] = CoinMin(columnUpper[iColumn],largest);
        newLower[iColumn] = CoinMax(columnLower[iColumn],-largest);
      }
    }
  }
  int numberRows = model.getNumRows();
  const double * rowLower = model.getRowLower();
  const double * rowUpper = model.getRowUpper();
#ifndef NDEBUG
  double large2= 1.0e10*large;
#endif
#define MAXPASS 10

  // Loop round seeing if we can tighten bounds
  // Would be faster to have a stack of possible rows
  // and we put altered rows back on stack
  int numberCheck=-1;
  while(numberChanged>numberCheck) {

    numberChanged = 0; // Bounds tightened this pass
    
    if (iPass==MAXPASS) break;
    iPass++;
    
    for (iRow = 0; iRow < numberRows; iRow++) {

      if (rowLower[iRow]>-large||rowUpper[iRow]<large) {

	// possible row
	int infiniteUpper = 0;
	int infiniteLower = 0;
	double maximumUp = 0.0;
	double maximumDown = 0.0;
	double newBound;
	CoinBigIndex rStart = rowStart[iRow];
	CoinBigIndex rEnd = rowStart[iRow]+rowLength[iRow];
	CoinBigIndex j;
	// Compute possible lower and upper ranges
      
	for (j = rStart; j < rEnd; ++j) {
	  double value=element[j];
	  iColumn = column[j];
	  if (value > 0.0) {
	    if (newUpper[iColumn] >= large) {
	      ++infiniteUpper;
	    } else {
	      maximumUp += newUpper[iColumn] * value;
	    }
	    if (newLower[iColumn] <= -large) {
	      ++infiniteLower;
	    } else {
	      maximumDown += newLower[iColumn] * value;
	    }
	  } else if (value<0.0) {
	    if (newUpper[iColumn] >= large) {
	      ++infiniteLower;
	    } else {
	      maximumDown += newUpper[iColumn] * value;
	    }
	    if (newLower[iColumn] <= -large) {
	      ++infiniteUpper;
	    } else {
	      maximumUp += newLower[iColumn] * value;
	    }
	  }
	}
	// Build in a margin of error
	maximumUp += 1.0e-8*fabs(maximumUp);
	maximumDown -= 1.0e-8*fabs(maximumDown);
	double maxUp = maximumUp+infiniteUpper*1.0e31;
	double maxDown = maximumDown-infiniteLower*1.0e31;
	if (maxUp <= rowUpper[iRow] + tolerance && 
	    maxDown >= rowLower[iRow] - tolerance) {
	  
	  // Row is redundant - make totally free
	} else {
	  if (maxUp < rowLower[iRow] -100.0*tolerance ||
	      maxDown > rowUpper[iRow]+100.0*tolerance) {
	    // problem is infeasible - exit at once
	    numberInfeasible++;
	    break;
	  }
	  double lower = rowLower[iRow];
	  double upper = rowUpper[iRow];
	  for (j = rStart; j < rEnd; ++j) {
	    double value=element[j];
	    iColumn = column[j];
	    double nowLower = newLower[iColumn];
	    double nowUpper = newUpper[iColumn];
	    if (value > 0.0) {
	      // positive value
	      if (lower>-large) {
		if (!infiniteUpper) {
		  assert(nowUpper < large2);
		  newBound = nowUpper + 
		    (lower - maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumUp);
		} else if (infiniteUpper==1&&nowUpper>large) {
		  newBound = (lower -maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumUp);
		} else {
		  newBound = -COIN_DBL_MAX;
		}
		if (newBound > nowLower + 1.0e-12&&newBound>-large) {
		  // Tighten the lower bound 
		  newLower[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (nowUpper - newBound < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust
		  double now;
		  if (nowLower<-large) {
		    now=0.0;
		    infiniteLower--;
		  } else {
		    now = nowLower;
		  }
		  maximumDown += (newBound-now) * value;
		  nowLower = newBound;
		}
	      } 
	      if (upper <large) {
		if (!infiniteLower) {
		  assert(nowLower >- large2);
		  newBound = nowLower + 
		    (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumDown);
		} else if (infiniteLower==1&&nowLower<-large) {
		  newBound =   (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumDown);
		} else {
		  newBound = COIN_DBL_MAX;
		}
		if (newBound < nowUpper - 1.0e-12&&newBound<large) {
		  // Tighten the upper bound 
		  newUpper[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (newBound - nowLower < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust 
		  double now;
		  if (nowUpper>large) {
		    now=0.0;
		    infiniteUpper--;
		  } else {
		    now = nowUpper;
		  }
		  maximumUp += (newBound-now) * value;
		  nowUpper = newBound;
		}
	      }
	    } else {
	      // negative value
	      if (lower>-large) {
		if (!infiniteUpper) {
		  assert(nowLower < large2);
		  newBound = nowLower + 
		    (lower - maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumUp);
		} else if (infiniteUpper==1&&nowLower<-large) {
		  newBound = (lower -maximumUp) / value;
		  // relax if original was large
		  if (fabs(maximumUp)>1.0e8)
		    newBound += 1.0e-12*fabs(maximumUp);
		} else {
		  newBound = COIN_DBL_MAX;
		}
		if (newBound < nowUpper - 1.0e-12&&newBound<large) {
		  // Tighten the upper bound 
		  newUpper[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (newBound - nowLower < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust
		  double now;
		  if (nowUpper>large) {
		    now=0.0;
		    infiniteLower--;
		  } else {
		    now = nowUpper;
		  }
		  maximumDown += (newBound-now) * value;
		  nowUpper = newBound;
		}
	      }
	      if (upper <large) {
		if (!infiniteLower) {
		  assert(nowUpper < large2);
		  newBound = nowUpper + 
		    (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumDown);
		} else if (infiniteLower==1&&nowUpper>large) {
		  newBound =   (upper - maximumDown) / value;
		  // relax if original was large
		  if (fabs(maximumDown)>1.0e8)
		    newBound -= 1.0e-12*fabs(maximumDown);
		} else {
		  newBound = -COIN_DBL_MAX;
		}
		if (newBound > nowLower + 1.0e-12&&newBound>-large) {
		  // Tighten the lower bound 
		  newLower[iColumn] = newBound;
		  numberChanged++;
		  // check infeasible (relaxed)
		  if (nowUpper - newBound < 
		      -100.0*tolerance) {
		    numberInfeasible++;
		  }
		  // adjust
		  double now;
		  if (nowLower<-large) {
		    now=0.0;
		    infiniteUpper--;
		  } else {
		    now = nowLower;
		  }
		  maximumUp += (newBound-now) * value;
		  nowLower = newBound;
		}
	      }
	    }
	  }
	}
      }
    }
    totalTightened += numberChanged;
    if (iPass==1)
      numberCheck=numberChanged>>4;
    if (numberInfeasible) break;
  }
  if (!numberInfeasible) {
    // Set bounds slightly loose unless integral
    double useTolerance = 1.0e-2;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (columnUpper[iColumn]>columnLower[iColumn]) {
        double lower = newLower[iColumn];
        double upper = newUpper[iColumn];
        if (model.isInteger(iColumn)) {
          if (fabs(lower-floor(lower+0.5))<1.0e-5)
            lower=floor(lower+0.5);
          else
            lower = ceil(lower);
          if (fabs(upper-floor(upper+0.5))<1.0e-5)
            upper=floor(upper+0.5);
          else
            upper = floor(upper);
          if (lower>upper)
            numberInfeasible++;
        } else {
          if (fabs(upper)<1.0e-8&&fabs(lower)<1.0e-8) {
            lower=0.0;
            upper=0.0;
          } else {
            // Relax unless integral
            if (fabs(lower-floor(lower+0.5))>1.0e-9)
              lower -= useTolerance;
            else
              lower = floor(lower+0.5);
            lower=CoinMax(columnLower[iColumn],lower);
            if (fabs(upper-floor(upper+0.5))>1.0e-9)
              upper += useTolerance;
            else
              upper = floor(upper+0.5);
            upper=CoinMin(columnUpper[iColumn],upper);
          }
	}
        model.setColLower(iColumn,lower);
        model.setColUpper(iColumn,upper);
        newLower[iColumn]=lower;
        newUpper[iColumn]=upper;
      }
    }
    if (!numberInfeasible) {
      // check common bad formulations
      int numberChanges=0;
      for (iRow = 0; iRow < numberRows; iRow++) {
        if (rowLower[iRow]>-large||rowUpper[iRow]<large) {
          // possible row
          double sumFixed=0.0;
          int infiniteUpper = 0;
          int infiniteLower = 0;
          double maximumUp = 0.0;
          double maximumDown = 0.0;
          double largest = 0.0;
          CoinBigIndex rStart = rowStart[iRow];
          CoinBigIndex rEnd = rowStart[iRow]+rowLength[iRow];
          CoinBigIndex j;
          int numberInteger=0;
          int whichInteger=-1;
          // Compute possible lower and upper ranges
          for (j = rStart;j < rEnd; ++j) {
            double value=element[j];
            iColumn = column[j];
            if (newUpper[iColumn]>newLower[iColumn]) {
              if (model.isInteger(iColumn)) {
                numberInteger++;
                whichInteger=iColumn;
              }
              largest = CoinMax(largest,fabs(value));
              if (value > 0.0) {
                if (newUpper[iColumn] >= large) {
                  ++infiniteUpper;
                } else {
                  maximumUp += newUpper[iColumn] * value;
                }
                if (newLower[iColumn] <= -large) {
                  ++infiniteLower;
                } else {
                  maximumDown += newLower[iColumn] * value;
                }
              } else if (value<0.0) {
                if (newUpper[iColumn] >= large) {
                  ++infiniteLower;
                } else {
                  maximumDown += newUpper[iColumn] * value;
                }
                if (newLower[iColumn] <= -large) {
                  ++infiniteUpper;
                } else {
                  maximumUp += newLower[iColumn] * value;
                }
              }
            } else {
              // fixed
              sumFixed += newLower[iColumn]*value;
            }
          }
          // Adjust
          maximumUp += sumFixed;
          maximumDown += sumFixed;
          // For moment just when all one sign and ints
          //maximumUp += 1.0e-8*fabs(maximumUp);
          //maximumDown -= 1.0e-8*fabs(maximumDown);
          double gap = 0.0;
          if ((rowLower[iRow]>maximumDown&&largest>rowLower[iRow]-maximumDown)&&
              ((maximumUp<=rowUpper[iRow]&&!infiniteUpper)||rowUpper[iRow]>=1.0e30)) {
            gap = rowLower[iRow]-maximumDown;
            if (infiniteLower)
              gap=0.0; // switch off
          } else if ((maximumUp>rowUpper[iRow]&&largest>maximumUp-rowUpper[iRow])&&
                     ((maximumDown>=rowLower[iRow]&&!infiniteLower)||rowLower[iRow]<=-1.0e30)) {
            gap = -(maximumUp-rowUpper[iRow]);
            if (infiniteUpper)
              gap=0.0; // switch off
          }
          if (fabs(gap)>1.0e-8) {
            for (j = rStart;j < rEnd; ++j) {
              double value=element[j];
              iColumn = column[j];
              double difference = newUpper[iColumn]-newLower[iColumn];
              if (difference>0.0&&difference<=1.0) {
                double newValue=value;
                if (value*gap>0.0&&model.isInteger(iColumn)) {
                  if (fabs(value*difference) > fabs(gap)) {
                    // No need for it to be larger than
                    newValue = gap/difference;
                  }
                  if (fabs(value-newValue)>1.0e-12) {
                    numberChanges++;
		    // BUT variable may have bound
		    double rhsAdjust=0.0;
		    if (gap>0.0) {
		      // rowLower
		      if (value>0.0) {
			// new value is based on going up from lower bound
			if (colLower[iColumn]) 
			  rhsAdjust = colLower[iColumn]*(value-newValue);
		      } else {
			// new value is based on going down from upper bound
			if (colUpper[iColumn]) 
			  rhsAdjust = colUpper[iColumn]*(value-newValue);
		      }
		    } else {
		      // rowUpper
		      if (value<0.0) {
			// new value is based on going up from lower bound
			if (colLower[iColumn]) 
			  rhsAdjust = colLower[iColumn]*(value-newValue);
		      } else {
			// new value is based on going down from upper bound
			if (colUpper[iColumn]) 
			  rhsAdjust = colUpper[iColumn]*(value-newValue);
		      }
		    }
		    if (rhsAdjust) {
#ifdef CLP_INVESTIGATE
		      printf("FFor column %d bounds %g, %g on row %d bounds %g, %g coefficient was changed from %g to %g with rhs adjustment of %g\n",
			     iColumn,colLower[iColumn],colUpper[iColumn],
			     iRow,rowLower[iRow],rowUpper[iRow],
			     value,newValue,rhsAdjust);
#endif
		      if (rowLower[iRow]>-1.0e20)
			model.setRowLower(iRow,rowLower[iRow]-rhsAdjust);
		      if (rowUpper[iRow]<1.0e20)
			model.setRowUpper(iRow,rowUpper[iRow]-rhsAdjust);
#ifdef CLP_INVESTIGATE
		      printf("FFor column %d bounds %g, %g on row %d bounds %g, %g coefficient was changed from %g to %g with rhs adjustment of %g\n",
			     iColumn,colLower[iColumn],colUpper[iColumn],
			     iRow,rowLower[iRow],rowUpper[iRow],
			     value,newValue,rhsAdjust);
#endif
		    }
                    element[j]=newValue;
		    handler_->message(CGL_ELEMENTS_CHANGED2,messages_)
		      <<iRow<<iColumn<<value<<newValue
		      <<CoinMessageEol;
#ifdef CGL_DEBUG
                    const OsiRowCutDebugger * debugger = model.getRowCutDebugger();
                    if (debugger&&debugger->numberColumns()==numberColumns) {
                      const double * optimal = debugger->optimalSolution();
                      double sum=0.0;
                      for (int jj = rStart;jj < rEnd; ++jj) {
                        double value=element[j];
                        int jColumn = column[jj];
                        sum += value*optimal[jColumn];
                      }
                      assert (sum>=rowLower[iRow]-1.0e7&&sum<=rowUpper[iRow]+1.0e-7);
                    }
#endif
                  }
                }
              }
            }
          }
        }
      }
      if (numberChanges) {
	handler_->message(CGL_ELEMENTS_CHANGED1,messages_)
			<<numberChanges<<CoinMessageEol;
        model.replaceMatrixOptional(copy);
      }
    }
  }
  if (numberInfeasible) {
    handler_->message(CGL_INFEASIBLE,messages_)
      <<CoinMessageEol;
    // restore column bounds
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      model.setColLower(iColumn,columnLower[iColumn]);
      model.setColUpper(iColumn,columnUpper[iColumn]);
    }
  }
  delete [] newLower;
  delete [] newUpper;
  delete [] columnLower;
  delete [] columnUpper;
  return (numberInfeasible);
}

void
CglPreProcess::postProcess(OsiSolverInterface & modelIn)
{
  // Do presolves
  bool saveHint;
  OsiHintStrength saveStrength;
  originalModel_->getHintParam(OsiDoPresolveInInitial,saveHint,saveStrength);
  bool saveHint2;
  OsiHintStrength saveStrength2;
  originalModel_->getHintParam(OsiDoDualInInitial,
                        saveHint2,saveStrength2);
  OsiSolverInterface * clonedCopy=NULL;
  double saveObjectiveValue = modelIn.getObjValue();
  if (!modelIn.isProvenOptimal()) {
    CoinWarmStartBasis *slack =
      dynamic_cast<CoinWarmStartBasis *>(modelIn.getEmptyWarmStart()) ;
    modelIn.setWarmStart(slack);
    delete slack ;
    modelIn.resolve();
  }
  if (modelIn.isProvenOptimal()) {
    OsiSolverInterface * modelM = &modelIn;
    // If some cuts add back rows
    if (cuts_.sizeRowCuts()) {
      clonedCopy = modelM->clone();
      modelM=clonedCopy;
      int numberRowCuts = cuts_.sizeRowCuts() ;
      // add in
      CoinBuild build;
      for (int k = 0;k<numberRowCuts;k++) {
	const OsiRowCut * thisCut = cuts_.rowCutPointer(k) ;
	int n=thisCut->row().getNumElements();
	const int * columnCut = thisCut->row().getIndices();
	const double * elementCut = thisCut->row().getElements();
	double lower = thisCut->lb();
	double upper = thisCut->ub();
	build.addRow(n,columnCut,elementCut,lower,upper);
      }
      modelM->addRows(build);
      // basis is wrong
      CoinWarmStartBasis empty;
      modelM->setWarmStart(&empty);
    }
    for (int iPass=numberSolvers_-1;iPass>=0;iPass--) {
      OsiSolverInterface * model = model_[iPass];
      if (model->getNumCols()) {
	CoinWarmStartBasis* basis =
	  dynamic_cast <CoinWarmStartBasis*>(modelM->getWarmStart()) ;
	if (basis) {
	  model->setWarmStart(basis);
	  delete basis;
	}
	int numberColumns = modelM->getNumCols();
	const double * solutionM = modelM->getColSolution();
	const double * columnLower2 = model->getColLower(); 
	const double * columnUpper2 = model->getColUpper();
	const double * columnLower = modelM->getColLower(); 
	const double * columnUpper = modelM->getColUpper();
	int iColumn;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (modelM->isInteger(iColumn)) {
	    double value = solutionM[iColumn];
	    double value2 = floor(value+0.5);
	    // if test fails then empty integer
	    if (fabs(value-value2)<1.0e-3) {
	      model->setColLower(iColumn,value2);
	      model->setColUpper(iColumn,value2);
	    } else {
#ifdef COIN_DEVELOP
	      printf("NPASS=%d, ipass %d var %d values %g %g %g\n",
		     numberSolvers_,iPass,iColumn,model->getColLower()[iColumn],
		     value,model->getColUpper()[iColumn]);
#endif
	    }
	  } else if (columnUpper[iColumn]==columnLower[iColumn]) {
	    if (columnUpper2[iColumn]>columnLower2[iColumn]&&!model->isInteger(iColumn)) {
	      model->setColUpper(iColumn,columnLower[iColumn]);
	    }
	  }
	}
      }
      int numberColumns = modelM->getNumCols();
      const double * solutionM = modelM->getColSolution();
      int iColumn;
      // Give a hint to do primal
      //model->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
      model->setHintParam(OsiDoDualInInitial,false,OsiHintTry);
      // clean
/*
  VIRTUOUS - I am not happy here (JJF) - This was put in for Special Ordered Sets of type 2

  Previous loop has likely made nontrivial bound changes, hence invalidated
  solution. Why do we need this? We're about to do an initialSolve, which
  will overwrite solution. Perhaps belongs in same guarded block with
  following feasibility check? If this is necessary for clp, solution should
  be acquired before bounds changes.
*/
      if (0)
      {
	int numberColumns = model->getNumCols();
	const double * lower = model->getColLower();
	const double * upper = model->getColUpper();
	double * solution = CoinCopyOfArray(model->getColSolution(),numberColumns);
	int i;
	for ( i=0;i<numberColumns;i++) {
	  double value = solution[i];
	  value = CoinMin(value,upper[i]);
	  value = CoinMax(value,lower[i]);
	  solution[i]=value;
	}
	model->setColSolution(solution);
	delete [] solution;
      }
      if (0) {
	int numberColumns = model->getNumCols();
	int numberRows = model->getNumRows();
	const double * lower = model->getColLower();
	const double * upper = model->getColUpper();
	const double * rowLower = model->getRowLower();
	const double * rowUpper = model->getRowUpper();
#ifndef NDEBUG
	double primalTolerance=1.0e-8;
#endif
	// Column copy
	const CoinPackedMatrix * matrix = model->getMatrixByCol();
	const double * element = matrix->getElements();
	const int * row = matrix->getIndices();
	const CoinBigIndex * columnStart = matrix->getVectorStarts();
	const int * columnLength = matrix->getVectorLengths();
	double * rowActivity = new double[numberRows];
	memset(rowActivity,0,numberRows*sizeof(double));
	int i;
	for ( i=0;i<numberColumns;i++) {
	  int j;
	  double value = lower[i];
	  if (value<lower[i]) {
	    value=lower[i];
	  } else if (value>upper[i]) {
	    value=upper[i];
	  }
	  assert (upper[i]>=lower[i]);
	  assert ((fabs(value)<1.0e20));
	  if (value) {
	    for (j=columnStart[i];
		 j<columnStart[i]+columnLength[i];j++) {
	      int iRow=row[j];
	      rowActivity[iRow] += value*element[j];
	    }
	  }
	}
	// check was feasible - if not adjust (cleaning may move)
	for (i=0;i<numberRows;i++) {
	  if(rowActivity[i]<rowLower[i]) {
	    assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
	    rowActivity[i]=rowLower[i];
	  } else if(rowActivity[i]>rowUpper[i]) {
	    assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
	    rowActivity[i]=rowUpper[i];
	  }
	}
      }
      {
	int numberFixed=0;
	int numberColumns = model->getNumCols();
	const double * columnLower = model->getColLower(); 
	const double * columnUpper = model->getColUpper();
	int iColumn;
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (columnLower[iColumn]==columnUpper[iColumn])
	    numberFixed++;
	}
	if (numberColumns>2000&&numberFixed<numberColumns&&
	    numberFixed*5>numberColumns) {
	  model->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
	}
      }
      model->initialSolve();
      numberIterationsPost_ += model->getIterationCount();
      if (!model->isProvenOptimal()) {
#ifdef COIN_DEVELOP
	  model->writeMps("bad2");
	  printf("bad unwind in postprocess\n");
#endif
      } else {
	//model->writeMps("good2");
      }
      presolve_[iPass]->postsolve(true);
      OsiSolverInterface * modelM2;
      if (iPass)
	modelM2 = modifiedModel_[iPass-1];
      else
	modelM2 = startModel_;
      // and fix
      const int * originalColumns = presolve_[iPass]->originalColumns();
      const double * columnLower2 = modelM2->getColLower(); 
      const double * columnUpper2 = modelM2->getColUpper();
      const double * columnLower = modelM->getColLower(); 
      const double * columnUpper = modelM->getColUpper();
      const double * solutionM2 = modelM2->getColSolution();
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	int jColumn = originalColumns[iColumn];
	if (!modelM2->isInteger(jColumn)) {
	  if (columnUpper[iColumn]==columnLower[iColumn]) {
	    if (columnUpper2[jColumn]>columnLower2[jColumn]&&!modelM2->isInteger(jColumn)) {
	      double value = solutionM[iColumn];
	      value = CoinMax(value,columnLower[iColumn]);
	      value = CoinMin(value,columnUpper[iColumn]);
#ifdef COIN_DEVELOP
	      printf ("assuming %d fixed to solution of %g (was %g) - bounds %g %g, old bounds and sol %g %g\n",
		      jColumn,value,solutionM2[jColumn],columnLower2[jColumn],columnUpper2[jColumn],
		      columnLower[iColumn],columnUpper[iColumn]);
#endif
	      modelM2->setColUpper(jColumn,value);
	    }
	  } else {
#if 0
	    if (columnUpper2[jColumn]>columnLower2[jColumn]&&!modelM2->isInteger(jColumn)) {
	      double value = solutionM[iColumn];
	      value = CoinMax(value,columnLower[iColumn]);
	      value = CoinMin(value,columnUpper[iColumn]);
	      printf ("assuming %d not fixed to solution of %g (was %g) - bounds %g %g, old bounds and sol %g %g\n",
		      jColumn,value,solutionM2[jColumn],columnLower2[jColumn],columnUpper2[jColumn],
		      columnLower[iColumn],columnUpper[iColumn]);
	    }
#endif
	  }
	} else {
	  // integer - dupcol bounds may be odd so use solution
	  double value = floor(solutionM2[jColumn]+0.5);
	  if (value<columnLower2[jColumn]) {
	    //printf("changing lower bound for %d from %g to %g to allow feasibility\n",
	    //	   jColumn,columnLower2[jColumn],value);
	    modelM2->setColLower(jColumn,value);
	  } else if (value>columnUpper2[jColumn]) {
	    //printf("changing upper bound for %d from %g to %g to allow feasibility\n",
	    //	   jColumn,columnUpper2[jColumn],value);
	    modelM2->setColUpper(jColumn,value);
	  } 
	}
      }
      delete modifiedModel_[iPass];;
      delete model_[iPass];;
      delete presolve_[iPass];
      modifiedModel_[iPass]=NULL;
      model_[iPass]=NULL;
      presolve_[iPass]=NULL;
      modelM = modelM2;
    }
    // should be back to startModel_;
    OsiSolverInterface * model = originalModel_;
    // Use number of columns in original
    int numberColumns = model->getNumCols();
    const double * solutionM = modelM->getColSolution();
    int iColumn;
    const double * columnLower2 = model->getColLower(); 
    const double * columnUpper2 = model->getColUpper();
    const double * columnLower = modelM->getColLower(); 
    const double * columnUpper = modelM->getColUpper();
    int numberBadValues = 0;
    CoinWarmStartBasis* basis =
      dynamic_cast <CoinWarmStartBasis*>(modelM->getWarmStart()) ;
    if (basis) {
      model->setWarmStart(basis);
      delete basis;
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (modelM->isInteger(iColumn)) {
	double value = solutionM[iColumn];
	double value2 = floor(value+0.5);
	// if test fails then empty integer
	if (fabs(value-value2)<1.0e-3) {
	  value2 = CoinMax(CoinMin(value2,columnUpper[iColumn]),columnLower[iColumn]);
	  model->setColLower(iColumn,value2);
	  model->setColUpper(iColumn,value2);
	} else {
#ifdef COIN_DEVELOP
	  printf("NPASS=%d, ipass end var %d values %g %g %g\n",
		 numberSolvers_,iColumn,model->getColLower()[iColumn],
		 value,model->getColUpper()[iColumn]);
#endif
	  numberBadValues++;
	}
      } else if (columnUpper[iColumn]==columnLower[iColumn]) {
	if (columnUpper2[iColumn]>columnLower2[iColumn]&&!model->isInteger(iColumn)) {
	  model->setColUpper(iColumn,columnLower[iColumn]);
	}
      }
    }
    if(numberBadValues) {
      const CoinPackedMatrix * columnCopy = model->getMatrixByCol();
      const int * row = columnCopy->getIndices();
      const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
      const int * columnLength = columnCopy->getVectorLengths(); 
      const double * element = columnCopy->getElements();
      int numberRows = model->getNumRows();
      double * rowActivity = new double[numberRows];
      memset(rowActivity,0,numberRows*sizeof(double));
      double * solution = CoinCopyOfArray(solutionM,numberColumns);
      for ( iColumn=0;iColumn<numberColumns;iColumn++) {
	double value = solutionM[iColumn];
	if (modelM->isInteger(iColumn)) {
	  double value2 = floor(value+0.5);
	  // if test fails then empty integer
	  if (fabs(value-value2)<1.0e-3) 
	    value = value2;
	}
	solution[iColumn] = value;
	for (CoinBigIndex j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow=row[j];
	  rowActivity[iRow] += value*element[j];
	}
      }
      const double * rowLower = model->getRowLower();
      const double * rowUpper = model->getRowUpper();
      //const double * columnLower = model->getColLower(); 
      //const double * columnUpper = model->getColUpper();
      const double * objective = model->getObjCoefficients();
      double direction = model->getObjSense();
      int numberCheck=0;
      double tolerance;
      model->getDblParam(OsiPrimalTolerance,tolerance);
      tolerance *= 10.0;
      for ( iColumn=0;iColumn<numberColumns;iColumn++) {
	double value = solution[iColumn];
	if (model->isInteger(iColumn)) {
	  double value2 = floor(value);
	  // See if empty integer
	  if (value!=value2) {
	    numberCheck++;
	    int allowed=0;
	    // can we go up
	    double movement = value2+1.0-value;
	    CoinBigIndex j;
	    bool good=true;
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow=row[j];
#ifdef COIN_DEVELOP
	      if (rowLower[iRow]>-1.0e20&&rowUpper[iRow]<1.0e20)
		printf("odd row with both bounds %d %g %g - element %g\n",
		       iRow,rowLower[iRow],rowUpper[iRow],element[j]);
#endif
	      double newActivity = rowActivity[iRow] + movement*element[j];
	      if (newActivity>rowUpper[iRow]+tolerance||
		  newActivity<rowLower[iRow]-tolerance)
		good=false;
	    }
	    if (good)
	      allowed=1;
	    // can we go down
	    movement = value2-value;
	    good=true;
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow=row[j];
	      double newActivity = rowActivity[iRow] + movement*element[j];
	      if (newActivity>rowUpper[iRow]+tolerance||
		  newActivity<rowLower[iRow]-tolerance)
		good=false;
	    }
	    if (good)
	      allowed |= 2;
	    if (allowed) {
	      if (allowed==3) {
		if (direction*objective[iColumn]>0.0)
		  allowed=2;
		else
		  allowed=1;
	      }
	      if(allowed==1) value2++;
	      movement =  value2-value;
	      solution[iColumn]=value2;
	      model->setColLower(iColumn,value2);
	      model->setColUpper(iColumn,value2);
	      for (j=columnStart[iColumn];
		   j<columnStart[iColumn]+columnLength[iColumn];j++) {
		int iRow=row[j];
		rowActivity[iRow] += movement*element[j];
	      }
	    } else {
#ifdef COIN_DEVELOP
	      printf("Unable to move %d\n",iColumn);
#endif
	    }
	  }
	}
      }
      assert (numberCheck==numberBadValues);
      model->setColSolution(solution);
      delete [] rowActivity;
      delete [] solution;
    }
  } else {
    // infeasible 
    for (int iPass=numberSolvers_-1;iPass>=0;iPass--) {
      delete modifiedModel_[iPass];;
      delete model_[iPass];;
      delete presolve_[iPass];
      modifiedModel_[iPass]=NULL;
      model_[iPass]=NULL;
      presolve_[iPass]=NULL;
    }
    // Back to startModel_;
    OsiSolverInterface * model = originalModel_;
    // Use number of columns in original
    int numberColumns = model->getNumCols();
    const double * columnLower = model->getColLower();
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (model->isInteger(iColumn)) 
	model->setColUpper(iColumn,columnLower[iColumn]);
    }
  }
  delete clonedCopy;
  originalModel_->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
  originalModel_->setHintParam(OsiDoDualInInitial,false,OsiHintTry);
  //printf("dumping ss.mps.gz in postprocess\n");
  //originalModel_->writeMps("ss");
  {
    int numberFixed=0;
    int numberColumns = originalModel_->getNumCols();
    const double * columnLower = originalModel_->getColLower(); 
    const double * columnUpper = originalModel_->getColUpper();
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (columnLower[iColumn]==columnUpper[iColumn])
	numberFixed++;
    }
    //printf("XX %d columns, %d fixed\n",numberColumns,numberFixed);
    //double time1 = CoinCpuTime();
    //originalModel_->initialSolve();
    //printf("Time with basis %g seconds, %d iterations\n",CoinCpuTime()-time1,originalModel_->getIterationCount());
    if (numberColumns<10000||numberFixed==numberColumns) {
      CoinWarmStart * empty = originalModel_->getEmptyWarmStart();
      originalModel_->setWarmStart(empty);
      delete empty;
    }
  }
  //double time1 = CoinCpuTime();
  originalModel_->initialSolve();
  numberIterationsPost_ += originalModel_->getIterationCount();
  //printf("Time without basis %g seconds, %d iterations\n",CoinCpuTime()-time1,originalModel_->getIterationCount());
  double objectiveValue = originalModel_->getObjValue();
  double testObj = 1.0e-8*CoinMax(fabs(saveObjectiveValue),
				  fabs(objectiveValue))+1.0e-4;
  if (!originalModel_->isProvenOptimal()) {
#ifdef COIN_DEVELOP
    originalModel_->writeMps("bad3");
    printf("bad end unwind in postprocess\n");
#endif
    handler_->message(CGL_POST_INFEASIBLE,messages_)
      <<CoinMessageEol;
  } else if (fabs(saveObjectiveValue-objectiveValue)>testObj) {
    handler_->message(CGL_POST_CHANGED,messages_)
      <<saveObjectiveValue<<objectiveValue
      <<CoinMessageEol;
  }
  originalModel_->setHintParam(OsiDoDualInInitial,saveHint2,saveStrength2);
  originalModel_->setHintParam(OsiDoPresolveInInitial,saveHint,saveStrength);
}
//-------------------------------------------------------------------
// Returns the greatest common denominator of two 
// positive integers, a and b, found using Euclid's algorithm 
//-------------------------------------------------------------------
static int gcd(int a, int b) 
{
  int remainder = -1;
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
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
}
/* Return model with useful modifications.  
   If constraints true then adds any x+y=1 or x-y=0 constraints
   If NULL infeasible
*/
OsiSolverInterface * 
CglPreProcess::modified(OsiSolverInterface * model,
			bool constraints,
			int & numberChanges,
                        int iBigPass,
			int numberPasses)
{
  OsiSolverInterface * newModel = model->clone();
  OsiCuts twoCuts;
  int numberRows = newModel->getNumRows();
  int numberColumns = newModel->getNumCols();
  int number01Integers=0;
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (newModel->isBinary(iColumn))
      number01Integers++;
  }
  OsiRowCut ** whichCut = new OsiRowCut * [numberRows+1];
  numberChanges=0;
  CoinThreadRandom randomGenerator;
  CglTreeProbingInfo info(model);
  info.level = 0;
  info.pass = 0;
  info.formulation_rows = numberRows;
  info.inTree = false;
  info.options = !numberProhibited_ ? 0 : 2;
  info.randomNumberGenerator=&randomGenerator;
  info.strengthenRow= whichCut;
  // See if user asked for heavy probing
  bool heavyProbing=false;
  for (int iGenerator=0;iGenerator<numberCutGenerators_;iGenerator++) {
    CglProbing * probingCut = dynamic_cast<CglProbing *> (generator_[iGenerator]);
    if (probingCut&&probingCut->getMaxPassRoot()>1) {
      heavyProbing=true;
      break;
    }
  }
  bool feasible=true;
  int firstGenerator=0;
  int lastGenerator=numberCutGenerators_;
  for (int iPass=0;iPass<numberPasses;iPass++) {
    // Statistics
    int numberFixed=0;
    int numberTwo=twoCuts.sizeRowCuts();
    int numberStrengthened=0;
    info.pass = iPass;
    info.options=0;
    int numberChangedThisPass=0;
    /*
      needResolve    solution is stale
      rebuilt   constraint system deleted and recreated (implies initialSolve)
    */
    for (int iGenerator=firstGenerator;iGenerator<lastGenerator;iGenerator++) {
      bool needResolve = false ;
      bool rebuilt = false ;
      OsiCuts cs;
      CoinZeroN(whichCut,numberRows);
      CglProbing * probingCut=NULL;
      int numberFromCglDuplicate=0;
      const int * duplicate=NULL;
      CglDuplicateRow * dupRow = NULL;
      CglClique * cliqueGen = NULL;
      if (iGenerator>=0) {
        //char name[20];
        //sprintf(name,"prex%2.2d.mps",iGenerator);
        //newModel->writeMpsNative(name, NULL, NULL,0,1,0);
        // refresh as model may have changed
        generator_[iGenerator]->refreshSolver(newModel);
        // skip duplicate rows except once
        dupRow = dynamic_cast<CglDuplicateRow *> (generator_[iGenerator]);
	cliqueGen = dynamic_cast<CglClique *> (generator_[iGenerator]);
        if ((dupRow||cliqueGen)&&(iPass/*||iBigPass*/))
            continue;
        probingCut = dynamic_cast<CglProbing *> (generator_[iGenerator]);
	if (!probingCut) {
	  generator_[iGenerator]->generateCuts(*newModel,cs,info);
	} else {
	  info.options=64;
	  probingCut->generateCutsAndModify(*newModel,cs,&info);
	}
#if 1 //def CLIQUE_ANALYSIS
	if (probingCut) {
	  //printf("ordinary probing\n");
	  info.analyze(*newModel);
	} 
#endif
        // If CglDuplicate may give us useless rows
        if (dupRow) {
          numberFromCglDuplicate = dupRow->numberOriginalRows();
          duplicate = dupRow->duplicate();
        }
	if (cliqueGen&&cs.sizeRowCuts()) {
	  int n = cs.sizeRowCuts();
	  printf("%d clique cuts\n",n);
	  OsiSolverInterface * copySolver = newModel->clone();
	  numberRows=copySolver->getNumRows();
	  copySolver->applyCuts(cs);
	  //static int kk=0;
	  //char name[20];
	  //kk++;
	  //sprintf(name,"matrix%d",kk);
	  //printf("writing matrix %s\n",name);
	  //copySolver->writeMps(name);
	  CglDuplicateRow dupCuts(copySolver);
	  dupCuts.setMode(8);
	  OsiCuts cs2;
	  dupCuts.generateCuts(*copySolver,cs2,info);
	  // -1 not used, -2 delete, -3 not clique
	  const int * duplicate = dupCuts.duplicate();
	  // -1 not used, >=0 earliest row affected
	  const int * used = duplicate+numberRows+n;
	  int numberDrop=0;
	  int * drop = new int[numberRows];
	  for (int iRow=0;iRow<numberRows;iRow++) {
	    if (duplicate[iRow]==-2) 
	      drop[numberDrop++]=iRow;
	  }
	  int nOther=0;
	  for (int iRow=numberRows+n-1;iRow>=numberRows;iRow--) {
#if 1
	    int earliest = used[iRow];
	    while (earliest>=numberRows) {
	      if (duplicate[earliest]==-2) 
		earliest = used[earliest];
	      else
		break;
	    }
#else
	    int earliest=0;
#endif
	    if (duplicate[iRow]==-2||earliest==-1||earliest>=numberRows) { 
	      cs.eraseRowCut(iRow-numberRows);
	      nOther++;
	    }
	  }
	  n -= nOther;
	  int newNumberRows = numberRows-numberDrop+n;
	  bool special = (cliqueGen->getMinViolation()==-3.0);
	  printf("could drop %d rows - current nrows %d other %d - new nrows %d\n",
		 numberDrop,numberRows,nOther,newNumberRows);
	  if (n<=numberDrop||special) {
	    printf("Dropping rows current nrows %d - new nrows %d\n",
		   numberRows,newNumberRows);
	    if (newNumberRows>numberRows) {
	      // need new array
	      delete [] whichCut;
	      whichCut = new OsiRowCut * [newNumberRows+1];
	      CoinZeroN(whichCut,newNumberRows);
	      info.strengthenRow= whichCut;
	    }
	    newModel->deleteRows(numberDrop,drop);
	    // may be able to delete some added cliques
	    newModel->applyCuts(cs);
	    numberRows=newModel->getNumRows();
	    newModel->resolve();
#if 0
	    int numberRows2=copySolver->getNumRows();
	    const double * rowLower = copySolver->getRowLower();
	    const double * rowUpper = copySolver->getRowUpper();
	    const CoinPackedMatrix * matrixByRow = copySolver->getMatrixByRow();
	    // Row copy
	    const double * elementByRow = matrixByRow->getElements();
	    const int * column = matrixByRow->getIndices();
	    const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
	    const int * rowLength = matrixByRow->getVectorLengths();
	    const double * solution = newModel->getColSolution();
	    for (int iRow=0;iRow<numberRows2;iRow++) {
	      double sum=0.0;
	      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
		int iColumn = column[j];
		double value = elementByRow[j];
		sum += value*solution[iColumn];
	      }
	      assert (sum>rowLower[iRow]-1.0e-4&&sum<rowUpper[iRow]+1.0e-4);
	    }
#endif
	  }
	  delete copySolver;
	  delete [] drop;
	  continue;
	  //for (int i=0;i<n;i++) {
	  //OsiRowCut & thisCut = cs.rowCut(i);
	  //thisCut.print();
	  //}
	}
      } else {
#if 0
        // special probing
        CglProbing generator1;
        probingCut=&generator1;
        generator1.setUsingObjective(false);
        generator1.setMaxPass(1);
        generator1.setMaxPassRoot(1);
        generator1.setMaxProbeRoot(100);
        generator1.setMaxLook(100);
        generator1.setRowCuts(3);
	if (heavyProbing) {
	  generator1.setMaxElements(300);
	  generator1.setMaxProbeRoot(model->getNumCols());
	}
	// out for now - think about cliques
        if(!generator1.snapshot(*newModel,NULL,false)) {
          generator1.createCliques(*newModel,2,1000,true);
          generator1.setMode(0);
          // To get special stuff
          info.pass=4;
          CoinZeroN(whichCut,numberRows);
          generator1.generateCutsAndModify(*newModel,cs,&info);
#ifdef CLIQUE_ANALYSIS
	  printf("special probing\n");
	  info.analyze(*newModel);
#endif
        } else {
          feasible=false;
        }
#endif
      }
      // check changes
      // first are any rows strengthened by cuts
      int iRow;
      for (iRow=0;iRow<numberRows;iRow++) {
        if(whichCut[iRow])
          numberStrengthened++;
      }
      // Also can we get rid of duplicate rows
      int numberDrop=0;
      for (iRow=0;iRow<numberFromCglDuplicate;iRow++) {
        if (duplicate[iRow]==-2||duplicate[iRow]>=0) {
          numberDrop++;
          newModel->setRowBounds(iRow,-COIN_DBL_MAX,COIN_DBL_MAX);
        }
      }
      const double * columnLower = newModel->getColLower();
      const double * columnUpper = newModel->getColUpper();
      if ((numberStrengthened||numberDrop)&&feasible) {
	/*
	  
	Deleting all rows and rebuilding invalidates everything, initialSolve will
	be required.
	*/
	needResolve = true ;
	rebuilt = true ;
        // Easier to recreate entire matrix
        const CoinPackedMatrix * rowCopy = newModel->getMatrixByRow();
        const int * column = rowCopy->getIndices();
        const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
        const int * rowLength = rowCopy->getVectorLengths(); 
        const double * rowElements = rowCopy->getElements();
        const double * rowLower = newModel->getRowLower();
        const double * rowUpper = newModel->getRowUpper();
        CoinBuild build;
	// For basis
	char * keepRow = new char [numberRows];
        for (iRow=0;iRow<numberRows;iRow++) {
	  keepRow[iRow]=0;
          OsiRowCut * thisCut = whichCut[iRow];
          whichCut[iRow]=NULL;
          if (rowLower[iRow]>-1.0e20||rowUpper[iRow]<1.0e20) {
#if 0
	    if (thisCut) {
	      printf("Cut on row %d\n",iRow);
	      thisCut->print();
	    }
#endif
            if (!thisCut) {
              // put in old row
              int start=rowStart[iRow];
	      int kInt=-1;
	      double newValue=0.0;
	      if (!iPass&&!iBigPass) {
		// worthwhile seeing if odd gcd
		int end = start + rowLength[iRow];
		double rhsAdjustment=0.0;
		int nPosInt=0;
		int nNegInt=0;
		// Find largest integer coefficient
		int k;
		for ( k = start; k < end; ++k) {
		  CoinBigIndex j = column[k];
		  if (columnUpper[j]>columnLower[j]) {
		    if (newModel->isInteger(j)) {
		      if (rowElements[k]>0.0)
			nPosInt++;
		      else
			nNegInt++;
		    } else {
		      break; // no good
		    }
		  } else {
		    rhsAdjustment += columnLower[j]*rowElements[k];
		  }
		}
		if (k==end) {
		  // see if singleton coefficient can be strengthened
		  if ((nPosInt==1&&nNegInt>1)||(nNegInt==1&&nPosInt>1)) {
		    double lo;
		    double up;
		    if (rowLower[iRow]>-1.0e20)
		      lo = rowLower[iRow]-rhsAdjustment;
		    else
		      lo=-COIN_DBL_MAX;
		    if(rowUpper[iRow]<1.0e20)
		      up = rowUpper[iRow]-rhsAdjustment;
		    else
		      up=COIN_DBL_MAX;
		    double multiplier=1.0;
		    if (nNegInt==1) {
		      // swap signs
		      multiplier=lo;
		      lo=-up;
		      up=-multiplier;
		      multiplier=-1.0;
		    }
		    bool possible=true;
		    double singletonValue=0;
		    double scale = 4.0*5.0*6.0;
		    int kGcd=-1;
		    double smallestSum = 0.0;
		    double largestSum = 0.0;
		    for ( k = start; k < end; ++k) {
		      CoinBigIndex j = column[k];
		      double value=multiplier*rowElements[k];
		      if (columnUpper[j]>columnLower[j]) {
			if (value>0.0) {
			  // singleton
			  kInt=j;
			  if (columnUpper[j]-columnLower[j]!=1.0) {
			    possible = false;
			    break;
			  } else {
			    singletonValue=value;
			  }
			} else {
			  if (columnLower[j]>-1.0e10)
			    smallestSum += value*columnLower[j];
			  else
			    smallestSum = -COIN_DBL_MAX;
			  if (columnUpper[j]<-1.0e10)
			    largestSum += value*columnUpper[j];
			  else
			    largestSum = COIN_DBL_MAX;
			  value *=-scale;
			  if (fabs(value-floor(value+0.5))>1.0e-12) {
			    possible=false;
			    break;
			  } else {
			    int kVal = static_cast<int> (floor(value+0.5));
			    if (kGcd>0) 
			      kGcd = gcd(kGcd,kVal);
			    else
			      kGcd=kVal;
			  }
			}
		      }
		    }
		    if (possible) {
		      double multiple = (static_cast<double> (kGcd))/scale;
		      int interesting=0;
		      double saveLo=lo;
		      double saveUp=up;
		      double nearestLo0=lo;
		      double nearestUp0=up;
		      double nearestLo1=lo;
		      double nearestUp1=up;
		      // adjust rhs for singleton 
		      if (lo!=-COIN_DBL_MAX) {
			// singleton at lb
			lo -= columnLower[kInt]*singletonValue;
			double exact = lo/multiple;
			if (fabs(exact-floor(exact+0.5))>1.0e-4) {
			  interesting +=1;
			  nearestLo0 = ceil(exact)*multiple;
			} 
			// singleton at ub
			lo -= singletonValue;
			exact = lo/multiple;
			if (fabs(exact-floor(exact+0.5))>1.0e-4) {
			  interesting +=2;
			  nearestLo1 = ceil(exact)*multiple;
			}
		      }
		      if (up!=COIN_DBL_MAX) {
			// singleton at lb
			up -= columnLower[kInt]*singletonValue;
			double exact = up/multiple;
			if (fabs(exact-floor(exact+0.5))>1.0e-4) {
			  interesting +=4;
			  nearestUp0 = floor(exact)*multiple;
			} 
			// singleton at ub
			up -= singletonValue;
			exact = up/multiple;
			if (fabs(exact-floor(exact+0.5))>1.0e-4) {
			  interesting +=8;
			  nearestUp1 = floor(exact)*multiple;
			} 
		      }
		      if (interesting) {
#ifdef CLP_INVESTIGATE
			printf("Row %d interesting %d lo,up %g,%g singleton %d value %g bounds %g,%g - gcd %g\n",iRow,interesting,saveLo,saveUp,kInt,singletonValue,
			       columnLower[kInt],columnUpper[kInt],multiple);
			printf("Orig lo,up %g,%g %d\n",rowLower[iRow],rowUpper[iRow],end-start);
			for ( k = start; k < end; ++k) {
			  CoinBigIndex j = column[k];
			  double value=multiplier*rowElements[k];
			  printf(" (%d, %g - bds %g, %g)",j,value,
				 columnLower[j],columnUpper[j]);
			}
			printf("\n");
#endif
			if(columnLower[kInt]) {
#ifdef CLP_INVESTIGATE
			  printf("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ\n"); //think
#endif
			  interesting=0;
			}
			newValue = singletonValue;
			double newLoRhs = rowLower[iRow];
			double newUpRhs = rowUpper[iRow];
			if ((interesting&3)!=0) {
			  newLoRhs = nearestLo0;
			  newValue = nearestLo0-nearestLo1;
			}
			if (saveLo==saveUp&&((interesting&5)==5||(interesting&10)==10)) {
#ifdef CLP_INVESTIGATE
			  printf("INFEAS? ");
#endif
			  interesting=0; //ninfeas++;
			}
			if ((interesting&12)) {
#ifdef CLP_INVESTIGATE
			  double value2 = newValue;
#endif
			  newUpRhs = nearestUp0;
			  newValue = nearestUp0-nearestUp1;
#ifdef CLP_INVESTIGATE
			  if (newValue!=value2) {
			    printf("??? old newvalue %g ",newValue);
			  }
#endif
			}
#ifdef CLP_INVESTIGATE
			printf("guess is new lo %g, new up %g, new value %g\n",
			       newLoRhs,newUpRhs,newValue);
#endif
		      }
		      // Just do mzzv11 case
		      double exact = singletonValue/multiple;
		      if (fabs(exact-floor(exact+0.5))<1.0e-5)
			interesting &= ~2;
		      if (!smallestSum&&interesting==2&&!saveLo&&saveUp>1.0e20) {
			newValue = multiple*floor(exact);
			newValue *= multiplier;
#ifdef CLP_INVESTIGATE
			printf("New coefficient for %d will be %g\n",kInt,newValue);
#endif
		      } else {
			// don't do
			kInt=-1;
		      }
		    } else {
		      kInt=-1;
		    }
		  }
		}
		// endgcd
	      }
	      int length = rowLength[iRow];
	      if (kInt<0) {
		build.addRow(length,column+start,rowElements+start,
			     rowLower[iRow],rowUpper[iRow]);
	      } else {
		double * els = CoinCopyOfArray(rowElements+start,length);
		if (fabs(newValue)>1.0e-13) {
		  for ( int k = 0; k < length; ++k) {
		    int j = column[k+start];
		    if (j==kInt) {
		      els[k] = newValue;
		      break;
		    }
		  }
		  build.addRow(length,column+start,els,
			       rowLower[iRow],rowUpper[iRow]);
		} else {
		  // strengthened to zero!
#ifdef CLP_INVESTIGATE
		  printf("CglPreProcess - element strenthened to zero!\n");
#endif
		  int * cols = CoinCopyOfArray(column+start,length);
		  int n=0;
		  for ( int k = 0; k < length; ++k) {
		    int j = cols[k];
		    if (j!=kInt) {
		      els[n] = els[k];
		      cols[n++]=j;
		    }
		  }
		  build.addRow(n,cols,els,
			       rowLower[iRow],rowUpper[iRow]);
		  delete [] cols;
		}
		delete [] els;
	      }
	      keepRow[iRow]=1;
            } else {
              // strengthens this row (should we check?)
              // may be worth going round again
              numberChangedThisPass++;
              int n=thisCut->row().getNumElements();
              const int * columnCut = thisCut->row().getIndices();
              const double * elementCut = thisCut->row().getElements();
              double lower = thisCut->lb();
              double upper = thisCut->ub();
              if (probingCut) {
                int i;
                int nFree=0;
                for ( i=0;i<n;i++) {
                  int iColumn = columnCut[i];
                  if (columnUpper[iColumn]>columnLower[iColumn]+1.0e-12)
                    nFree++;
                }
                bool good = (n==nFree);
		nFree=0;
                int n1=rowLength[iRow];
                int start=rowStart[iRow];
                for ( i=0;i<n1;i++) {
                  int iColumn = column[start+i];
                  if (columnUpper[iColumn]>columnLower[iColumn]+1.0e-12)
                    nFree++;
                }
                if (n1!=nFree) 
		  good=false;
                if (good) {
#if 0
                  printf("Original row %.8d %g <= ",iRow,rowLower[iRow]);
                  for ( i=0;i<n1;i++) 
		    printf("%g * x%d ",rowElements[start+i],column[start+i]);
                  printf("<= %g\n",rowUpper[iRow]);
                  printf("New                   %g <= ",lower);
                  for ( i=0;i<n;i++) 
		    printf("%g * x%d ",elementCut[i],columnCut[i]);
                  printf("<= %g\n",upper);
#endif
                } else {
                  // can't use
                  n=-1;
		  numberStrengthened--;
                  // put in old row
                  int start=rowStart[iRow];
                  build.addRow(rowLength[iRow],column+start,rowElements+start,
                               rowLower[iRow],rowUpper[iRow]);
		  keepRow[iRow]=1;
                }
              }
              if (n>0) {
                build.addRow(n,columnCut,elementCut,lower,upper);
		keepRow[iRow]=1;
              } else if (!n) {
                // Either infeasible or redundant
                if (lower<=0.0&&upper>=0.0) {
                  // redundant - row will go
                } else {
                  // infeasible!
                  feasible=false;
                  break;
                }
              }
            }
          }
        }
	if (rowType_) {
	  assert (numberRowType_==numberRows);
	  int numberRowType_=0;
	  for (iRow=0;iRow<numberRows;iRow++) {
	    if (keepRow[iRow]) {
	      rowType_[numberRowType_++]=rowType_[iRow];
	    }
	  }
	}
	// recreate
        int * del = new int[numberRows];
	// save and modify basis
	CoinWarmStartBasis* basis =
	  dynamic_cast <CoinWarmStartBasis*>(newModel->getWarmStart()) ;
	if (basis) {
	  int n=0;
	  for (iRow=0;iRow<numberRows;iRow++) {
	    if (!keepRow[iRow]) {
	      del[n++]=iRow;
	    }
	  }
	  basis->deleteRows(n,del);
	}
        for (iRow=0;iRow<numberRows;iRow++) 
          del[iRow]=iRow;
        newModel->deleteRows(numberRows,del);
        newModel->addRows(build);
        numberRows = newModel->getNumRows();
	if (basis) {
	  assert (numberRows==basis->getNumArtificial());
	  newModel->setWarmStart(basis);
	  delete basis;
	}
        if (dupRow&&cs.sizeRowCuts())
	  newModel->applyCuts(cs);
	delete [] keepRow;
        delete [] del;
/*
  VIRTUOUS

  A solver is not required to hold these pointers constant across complete
  row deletion and rebuild.
*/
	columnLower = newModel->getColLower();
	columnUpper = newModel->getColUpper();
      }
      if (!feasible)
        break;
      // now see if we have any x=y x+y=1
      if (constraints) {
        int numberRowCuts = cs.sizeRowCuts() ;
        for (int k = 0;k<numberRowCuts;k++) {
          OsiRowCut * thisCut = cs.rowCutPtr(k) ;
          int n=thisCut->row().getNumElements();
          double lower = thisCut->lb();
          double upper = thisCut->ub();
          if (n==2&&lower==upper) {
            twoCuts.insert(*thisCut);
          }
        }
      }
      if (probingCut) {
	// we could look harder for infeasibilities 
	assert (info.numberVariables()==numberColumns);
	int number01 = info.numberIntegers();
	const cliqueEntry * entry = info.fixEntries();
	const int * toZero = info.toZero();
	const int * toOne = info.toOne();
	const int * which = info.integerVariable();
	int numberBounds=0;
	char * markLB = new char [number01];
	memset(markLB,0,number01);
	char * markUB = new char [number01];
	memset(markUB,0,number01);
	CoinRelFltEq equality;
	for (int k=0;k<number01;k++) {
	  int start = toZero[k];
	  int end = toOne[k];
	  // to zero
	  int j;
	  for (j=start;j<end;j++) {
	    int goingToOne = oneFixesInCliqueEntry(entry[j]);
	    int v = sequenceInCliqueEntry(entry[j]);
	    if (v>=number01)
	      continue;
	    if (goingToOne) {
	      // If v -> 1 means k -> 0 we have k+v==1 
	      int startV = toOne[v];
	      int endV = toZero[v+1];
	      for (int jv=startV;jv<endV;jv++) {
		if(k==static_cast<int> (sequenceInCliqueEntry(entry[jv]))) {
		  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
		  double lo,up;
		  if (!goingToOneV) {
		    lo=1.0;
		    up=1.0;
		    OsiRowCut thisCut;
		    thisCut.setLb(lo);
		    thisCut.setUb(up);
		    double values[]={1.0,1.0};
		    int columns[2];
		    columns[0]=which[k];
		    columns[1]=which[v];
		    thisCut.setRow(2,columns,values,false);
		    twoCuts.insertIfNotDuplicate(thisCut,equality);
		  } else {
		    // means infeasible for k to go to 0
		    markLB[k]=1;
		    numberBounds++;
		  }
		  break;
		}
	      }
	    } else {
	      // If v -> 0 means k -> 0 we have k==v 
	      int startV = toZero[v];
	      int endV = toOne[v];
	      for (int jv=startV;jv<endV;jv++) {
		if(k==static_cast<int> (sequenceInCliqueEntry(entry[jv]))) {
		  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
		  double lo,up;
		  if (!goingToOneV) {
		    lo=0.0;
		    up=0.0;
		    OsiRowCut thisCut;
		    thisCut.setLb(lo);
		    thisCut.setUb(up);
		    double values[]={1.0,-1.0};
		    int columns[2];
		    columns[0]=which[k];
		    columns[1]=which[v];
		    thisCut.setRow(2,columns,values,false);
		    twoCuts.insertIfNotDuplicate(thisCut,equality);
		  } else {
		    // means infeasible for k to go to 0
		    markLB[k]=1;
		    numberBounds++;
		  }
		  break;
		}
	      }
	    }
	  }
	  start = toOne[k];
	  end = toZero[k+1];
	  // to one
	  for (j=start;j<end;j++) {
	    int goingToOne = oneFixesInCliqueEntry(entry[j]);
	    int v = sequenceInCliqueEntry(entry[j]);
	    if (v>=number01)
	      continue;
	    if (goingToOne) {
	      // If v -> 1 means k -> 1 we have k==v 
	      int startV = toOne[v];
	      int endV = toZero[v+1];
	      for (int jv=startV;jv<endV;jv++) {
		if(k==static_cast<int> (sequenceInCliqueEntry(entry[jv]))) {
		  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
		  double lo,up;
		  if (goingToOneV) {
		    lo=0.0;
		    up=0.0;
		    OsiRowCut thisCut;
		    thisCut.setLb(lo);
		    thisCut.setUb(up);
		    double values[]={1.0,-1.0};
		    int columns[2];
		    columns[0]=which[k];
		    columns[1]=which[v];
		    thisCut.setRow(2,columns,values,false);
		    twoCuts.insertIfNotDuplicate(thisCut,equality);
		  } else {
		    // means infeasible for k to go to 1
		    markUB[k]=1;
		    numberBounds++;
		  }
		  break;
		}
	      }
	    } else {
	      // If v -> 0 means k -> 1 we have k+v==1 
	      int startV = toZero[v];
	      int endV = toOne[v];
	      for (int jv=startV;jv<endV;jv++) {
		if(k==static_cast<int> (sequenceInCliqueEntry(entry[jv]))) {
		  int goingToOneV = oneFixesInCliqueEntry(entry[jv]);
		  double lo,up;
		  if (goingToOneV) {
		    lo=1.0;
		    up=1.0;
		    OsiRowCut thisCut;
		    thisCut.setLb(lo);
		    thisCut.setUb(up);
		    double values[]={1.0,1.0};
		    int columns[2];
		    columns[0]=which[k];
		    columns[1]=which[v];
		    thisCut.setRow(2,columns,values,false);
		    twoCuts.insertIfNotDuplicate(thisCut,equality);
		  } else {
		    // means infeasible for k to go to 1
		    markUB[k]=1;
		    numberBounds++;
		  }
		  break;
		}
	      }
	    }
	  }
	}
	if (numberBounds) {
	  CoinPackedVector lbs;
	  CoinPackedVector ubs;
	  for (int k=0;k<number01;k++) {
	    if (markLB[k]&&markUB[k]) {
	      // infeasible
	      feasible=false;
	      break;
	    } else if (markLB[k]) {
	      lbs.insert(which[k],1.0);
	    } else if (markUB[k]) {
	      ubs.insert(which[k],0.0);
	    }
	  }
	  OsiColCut cc;
	  cc.setUbs(ubs);
	  cc.setLbs(lbs);
	  cc.setEffectiveness(1.0e-5);
	  cs.insert(cc);
	}
	delete [] markLB;
	delete [] markUB;
      }
      // see if we have any column cuts
      int numberColumnCuts = cs.sizeColCuts() ;
      int numberBounds=0;
      for (int k = 0;k<numberColumnCuts;k++) {
        OsiColCut * thisCut = cs.colCutPtr(k) ;
	/*
	  Nontrivial bound changes will invalidate current solution.
	*/
	if (thisCut->effectiveness() > 1.0) {
	  needResolve = true ;
	}
	const CoinPackedVector & lbs = thisCut->lbs() ;
	const CoinPackedVector & ubs = thisCut->ubs() ;
	int j ;
	int n ;
	const int * which ;
	const double * values ;
	n = lbs.getNumElements() ;
	which = lbs.getIndices() ;
	values = lbs.getElements() ;
	for (j = 0;j<n;j++) {
	  int iColumn = which[j] ;
          if (values[j]>columnLower[iColumn]&&values[j]>-1.0e20) {
            //printf("%d lower from %g to %g\n",iColumn,columnLower[iColumn],values[j]);
            newModel->setColLower(iColumn,values[j]) ;
            if (false) {
              OsiSolverInterface * xx = newModel->clone();
              xx->initialSolve();
              assert (xx->isProvenOptimal());
              delete xx;
            }
            numberChangedThisPass++;
            if (columnLower[iColumn]==columnUpper[iColumn]) {
              numberFixed++;
            } else {
              numberBounds++;
	      if (columnLower[iColumn]>columnUpper[iColumn]) 
		feasible=false;
	    }
          }
	}
	n = ubs.getNumElements() ;
	which = ubs.getIndices() ;
	values = ubs.getElements() ;
	for (j = 0;j<n;j++) {
	  int iColumn = which[j] ;
          if (values[j]<columnUpper[iColumn]&&values[j]<1.0e20) {
            //printf("%d upper from %g to %g\n",iColumn,columnUpper[iColumn],values[j]);
            newModel->setColUpper(iColumn,values[j]) ;
            if (false) {
              OsiSolverInterface * xx = newModel->clone();
              xx->initialSolve();
              assert (xx->isProvenOptimal());
              delete xx;
            }
            numberChangedThisPass++;
            if (columnLower[iColumn]==columnUpper[iColumn]) {
              numberFixed++;
            } else {
              numberBounds++;
	      if (columnLower[iColumn]>columnUpper[iColumn]) 
		feasible=false;
	    }
          }
        }
      }
      numberTwo = twoCuts.sizeRowCuts()-numberTwo;
      numberChanges += numberTwo + numberStrengthened/10;
      if (numberFixed||numberTwo||numberStrengthened||numberBounds)
        handler_->message(CGL_PROCESS_STATS,messages_)
          <<numberFixed<<numberBounds<<numberStrengthened<<numberTwo
          <<CoinMessageEol;
      if (!feasible)
        break;
      /*
	If solution needs to be refreshed, do resolve or initialSolve as appropriate.
      */
      if (needResolve) {
	if (rebuilt) {
	  // basis shot to bits?
	  //CoinWarmStartBasis *slack =
	  //dynamic_cast<CoinWarmStartBasis *>(newModel->getEmptyWarmStart()) ;
	  //newModel->setWarmStart(slack);
	  //delete slack ;
	  newModel->initialSolve() ;
	} else {
	  newModel->resolve() ;
	}
	numberIterationsPre_ += newModel->getIterationCount();
	feasible = newModel->isProvenOptimal();
	if (!feasible&&getCutoff()>1.0e20) {
	  // Better double check
	  CoinWarmStartBasis *slack =
	    dynamic_cast<CoinWarmStartBasis *>(newModel->getEmptyWarmStart()) ;
	  newModel->setWarmStart(slack);
	  delete slack ;
	  newModel->resolve() ;
	  numberIterationsPre_ += newModel->getIterationCount();
	  feasible = newModel->isProvenOptimal();
	  //if (!feasible)
	  //newModel->writeMpsNative("infeas.mps",NULL,NULL,2,1);
	}
      }
      if (!feasible)
        break;
    }
    if (!feasible)
      break;
    numberChanges +=  numberChangedThisPass;
    if (iPass<numberPasses-1) {
      if ((!numberFixed&&numberChangedThisPass<1000*(numberRows+numberColumns))||iPass==numberPasses-2) {
        // do special probing at end - but not if very last pass
	if (iBigPass<numberSolvers_-1) {
	  firstGenerator=-1;
	  lastGenerator=0;
	}
        iPass=numberPasses-2;
      }
    }
  }
  delete [] whichCut;
  int numberRowCuts = twoCuts.sizeRowCuts() ;
  if (numberRowCuts) {
    // add in x=y etc
    CoinBuild build;
    for (int k = 0;k<numberRowCuts;k++) {
      OsiRowCut * thisCut = twoCuts.rowCutPtr(k) ;
      int n=thisCut->row().getNumElements();
      const int * columnCut = thisCut->row().getIndices();
      const double * elementCut = thisCut->row().getElements();
      double lower = thisCut->lb();
      double upper = thisCut->ub();
      build.addRow(n,columnCut,elementCut,lower,upper);
    }
    newModel->addRows(build);
    if (rowType_) {
      // adjust
      int numberRows=newModel->getNumRows();
      char * temp = CoinCopyOfArrayPartial(rowType_,numberRows,numberRowType_);
      delete [] rowType_;
      rowType_ = temp;
      for (int iRow=numberRowType_;iRow<numberRows;iRow++)
	rowType_[iRow]=-1;
      numberRowType_=numberRows;
    }
  }
  if (!feasible) {
    delete newModel;
    newModel=NULL;
  }
  return newModel;
}

/* Default Constructor

*/
CglPreProcess::CglPreProcess() 

:
  originalModel_(NULL),
  startModel_(NULL),
  numberSolvers_(0),
  model_(NULL),
  modifiedModel_(NULL),
  presolve_(NULL),
  handler_(NULL),
  defaultHandler_(true),
  appData_(NULL),
  originalColumn_(NULL),
  originalRow_(NULL),
  numberCutGenerators_(0),
  generator_(NULL),
  numberSOS_(0),
  typeSOS_(NULL),
  startSOS_(NULL),
  whichSOS_(NULL),
  weightSOS_(NULL),
  numberProhibited_(0),
  numberIterationsPre_(0),
  numberIterationsPost_(0),
  prohibited_(NULL),
  numberRowType_(0),
  options_(0),
  rowType_(NULL)
{
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  messages_ = CglMessage();
}

// Copy constructor.

CglPreProcess::CglPreProcess(const CglPreProcess & rhs)
:
  numberSolvers_(rhs.numberSolvers_),
  defaultHandler_(rhs.defaultHandler_),
  appData_(rhs.appData_),
  originalColumn_(NULL),
  originalRow_(NULL),
  numberCutGenerators_(rhs.numberCutGenerators_),
  numberProhibited_(rhs.numberProhibited_),
  numberIterationsPre_(rhs.numberIterationsPre_),
  numberIterationsPost_(rhs.numberIterationsPost_),
  numberRowType_(rhs.numberRowType_),
  options_(rhs.options_)
{
  if (defaultHandler_) {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(rhs.handler_->logLevel());
  } else {
    handler_ = rhs.handler_;
  }
  messages_ = rhs.messages_;
  if (numberCutGenerators_) {
    generator_ = new CglCutGenerator * [numberCutGenerators_];
    for (int i=0;i<numberCutGenerators_;i++) {
      generator_[i]=rhs.generator_[i]->clone();
    }
  } else {
    generator_=NULL;
  }
  if (rhs.originalModel_) {
    originalModel_ = rhs.originalModel_;
    // If no make equality then solvers are same
    if (rhs.originalModel_!=rhs.startModel_) {
      startModel_=rhs.startModel_->clone();
    } else {
      startModel_=originalModel_;
    }
  } else {
    originalModel_=NULL;
    startModel_=NULL;
  }
  if (numberSolvers_) {
    model_ = new OsiSolverInterface * [numberSolvers_];
    modifiedModel_ = new OsiSolverInterface * [numberSolvers_];
    presolve_ = new OsiPresolve * [numberSolvers_];
    for (int i=0;i<numberSolvers_;i++) {
      model_[i]=rhs.model_[i]->clone();
      modifiedModel_[i]=rhs.modifiedModel_[i]->clone();
      presolve_[i]=new OsiPresolve(*rhs.presolve_[i]);
    }
  } else {
    model_=NULL;
    presolve_=NULL;
  }
  numberSOS_=rhs.numberSOS_;
  if (numberSOS_) {
    int numberTotal = rhs.startSOS_[numberSOS_];
    typeSOS_= CoinCopyOfArray(rhs.typeSOS_,numberSOS_);
    startSOS_= CoinCopyOfArray(rhs.startSOS_,numberSOS_+1);
    whichSOS_= CoinCopyOfArray(rhs.whichSOS_,numberTotal);
    weightSOS_= CoinCopyOfArray(rhs.weightSOS_,numberTotal);
  } else {
    typeSOS_ = NULL;
    startSOS_ = NULL;
    whichSOS_ = NULL;
    weightSOS_ = NULL;
  }
  prohibited_ = CoinCopyOfArray(rhs.prohibited_,numberProhibited_);
  rowType_ = CoinCopyOfArray(rhs.rowType_,numberRowType_);
  cuts_ = rhs.cuts_;
}
  
// Assignment operator 
CglPreProcess & 
CglPreProcess::operator=(const CglPreProcess& rhs)
{
  if (this!=&rhs) {
    gutsOfDestructor();
    numberSolvers_=rhs.numberSolvers_;
    defaultHandler_=rhs.defaultHandler_;
    appData_=rhs.appData_;
    numberCutGenerators_=rhs.numberCutGenerators_;
    numberProhibited_ = rhs.numberProhibited_;
    numberIterationsPre_ = rhs.numberIterationsPre_;
    numberIterationsPost_ = rhs.numberIterationsPost_;
    numberRowType_ = rhs.numberRowType_;
    options_ = rhs.options_;
    if (defaultHandler_) {
      handler_ = new CoinMessageHandler();
      handler_->setLogLevel(rhs.handler_->logLevel());
    } else {
      handler_ = rhs.handler_;
    }
    messages_ = rhs.messages_;
    if (numberCutGenerators_) {
      generator_ = new CglCutGenerator * [numberCutGenerators_];
      for (int i=0;i<numberCutGenerators_;i++) {
        generator_[i]=rhs.generator_[i]->clone();
      }
    }
    if (rhs.originalModel_) {
      originalModel_ = rhs.originalModel_;
      // If no make equality then solvers are same
      if (rhs.originalModel_!=rhs.startModel_) {
        startModel_=rhs.startModel_->clone();
      } else {
        startModel_=originalModel_;
      }
    } else {
      originalModel_=NULL;
      startModel_=NULL;
    }
    if (numberSolvers_) {
      model_ = new OsiSolverInterface * [numberSolvers_];
      modifiedModel_ = new OsiSolverInterface * [numberSolvers_];
      presolve_ = new OsiPresolve * [numberSolvers_];
      for (int i=0;i<numberSolvers_;i++) {
        model_[i]=rhs.model_[i]->clone();
        modifiedModel_[i]=rhs.modifiedModel_[i]->clone();
        presolve_[i]=new OsiPresolve(*rhs.presolve_[i]);
      }
    } else {
      model_=NULL;
      presolve_=NULL;
    }
    numberSOS_=rhs.numberSOS_;
    if (numberSOS_) {
      int numberTotal = rhs.startSOS_[numberSOS_];
      typeSOS_= CoinCopyOfArray(rhs.typeSOS_,numberSOS_);
      startSOS_= CoinCopyOfArray(rhs.startSOS_,numberSOS_+1);
      whichSOS_= CoinCopyOfArray(rhs.whichSOS_,numberTotal);
      weightSOS_= CoinCopyOfArray(rhs.weightSOS_,numberTotal);
    } else {
      typeSOS_ = NULL;
      startSOS_ = NULL;
      whichSOS_ = NULL;
      weightSOS_ = NULL;
    }
    prohibited_ = CoinCopyOfArray(rhs.prohibited_,numberProhibited_);
    rowType_ = CoinCopyOfArray(rhs.rowType_,numberRowType_);
    cuts_ = rhs.cuts_;
  }
  return *this;
}
  
// Destructor 
CglPreProcess::~CglPreProcess ()
{
  gutsOfDestructor();
}
// Clears out as much as possible (except solver)
void 
CglPreProcess::gutsOfDestructor()
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  if (startModel_!=originalModel_) 
    delete startModel_;
  startModel_=NULL;
  //delete originalModel_;
  originalModel_=NULL;
  int i;
  for (i=0;i<numberCutGenerators_;i++) {
    delete generator_[i];
  }
  delete [] generator_;
  generator_=NULL;
  for (i=0;i<numberSolvers_;i++) {
    delete model_[i];
    delete modifiedModel_[i];
    delete presolve_[i];
  }
  delete [] model_;
  delete [] modifiedModel_;
  delete [] presolve_;
  model_=NULL;
  presolve_=NULL;
  delete [] originalColumn_;
  delete [] originalRow_;
  originalColumn_=NULL;
  originalRow_=NULL;
  delete [] typeSOS_;
  delete [] startSOS_;
  delete [] whichSOS_;
  delete [] weightSOS_;
  typeSOS_ = NULL;
  startSOS_ = NULL;
  whichSOS_ = NULL;
  weightSOS_ = NULL;
  delete [] prohibited_;
  prohibited_=NULL;
  numberProhibited_=0;
  numberIterationsPre_=0;
  numberIterationsPost_=0;
  delete [] rowType_;
  rowType_=NULL;
  numberRowType_=0;
}
// Add one generator
void 
CglPreProcess::addCutGenerator(CglCutGenerator * generator)
{
  CglCutGenerator ** temp = generator_;
  generator_ = new CglCutGenerator * [numberCutGenerators_+1];
  memcpy(generator_,temp,numberCutGenerators_*sizeof(CglCutGenerator *));
  delete[] temp ;
  generator_[numberCutGenerators_++]=generator->clone(); 
}
//#############################################################################
// Set/Get Application Data
// This is a pointer that the application can store into and retrieve
// This field is the application to optionally define and use.
//#############################################################################

void CglPreProcess::setApplicationData(void * appData)
{
  appData_ = appData;
}
//-----------------------------------------------------------------------------
void * CglPreProcess::getApplicationData() const
{
  return appData_;
}
/* Set cutoff bound on the objective function.
   
When using strict comparison, the bound is adjusted by a tolerance to
avoid accidentally cutting off the optimal solution.
*/
void 
CglPreProcess::setCutoff(double value) 
{
  // Solvers know about direction
  double direction = originalModel_->getObjSense();
  originalModel_->setDblParam(OsiDualObjectiveLimit,value*direction); 
}

// Get the cutoff bound on the objective function - always as minimize
double 
CglPreProcess::getCutoff() const
{ 
  double value ;
  originalModel_->getDblParam(OsiDualObjectiveLimit,value) ;
  return value * originalModel_->getObjSense() ;
}
// Pass in Message handler (not deleted at end)
void 
CglPreProcess::passInMessageHandler(CoinMessageHandler * handler)
{
  if (defaultHandler_)
    delete handler_;
  defaultHandler_=false;
  handler_=handler;
}
// Set language
void 
CglPreProcess::newLanguage(CoinMessages::Language language)
{
  messages_ = CglMessage(language);
}
// Return a pointer to the original columns (without clique slacks)
const int * 
CglPreProcess::originalColumns() const
{
  if (!originalColumn_) 
    createOriginalIndices();
  return originalColumn_;
}
// Return a pointer to the original rows
const int * 
CglPreProcess::originalRows() const
{
  if (!originalRow_)
    createOriginalIndices();
  return originalRow_;
}
// create original columns and rows
void 
CglPreProcess::createOriginalIndices() const
{
  // Find last model and presolve
  int iPass;
  for (iPass=numberSolvers_-1;iPass>=0;iPass--) {
    if (presolve_[iPass])
      break;
  }
  int nRows,nColumns;
  if (iPass>=0) {
    nRows=model_[iPass]->getNumRows();
    nColumns=model_[iPass]->getNumCols();
  } else {
    nRows=originalModel_->getNumRows();
    nColumns=originalModel_->getNumCols();
  }
  delete [] originalColumn_;
  originalColumn_=new int [nColumns];
  delete [] originalRow_;
  originalRow_ = new int[nRows];
  if (iPass>=0) {
    memcpy(originalColumn_,presolve_[iPass]->originalColumns(),
           nColumns*sizeof(int));
    memcpy(originalRow_,presolve_[iPass]->originalRows(),
           nRows*sizeof(int));
    iPass--;
    for (;iPass>=0;iPass--) {
      const int * originalColumns = presolve_[iPass]->originalColumns();
      int i;
      for (i=0;i<nColumns;i++)
        originalColumn_[i]=originalColumns[originalColumn_[i]];
      const int * originalRows = presolve_[iPass]->originalRows();
      int nRowsNow=model_[iPass]->getNumRows();
      for (i=0;i<nRows;i++) {
	int iRow=originalRow_[i];
	if (iRow>=0&&iRow<nRowsNow)
	  originalRow_[i]=originalRows[iRow];
	else
	  originalRow_[i]=-1;
      }
    }
    std::sort(originalColumn_,originalColumn_+nColumns);
  } else {
    int i;
    for (i=0;i<nColumns;i++)
      originalColumn_[i]=i;
    for (i=0;i<nRows;i++)
      originalRow_[i]=i;
  }
}
// Update prohibited and rowType
void 
CglPreProcess::update(const OsiPresolve * pinfo,
		      const OsiSolverInterface * solver)
{
  if (prohibited_) {
    const int * original = pinfo->originalColumns();
    int numberColumns = solver->getNumCols();
    // number prohibited must stay constant
    int n=0;
    int i;
    for (i=0;i<numberProhibited_;i++) {
      if(prohibited_[i])
	n++;
    }
    int last=-1;
    int n2=0;
    for (i=0;i<numberColumns;i++) {
      int iColumn = original[i];
      assert (iColumn>last);
      last=iColumn;
      char p = prohibited_[iColumn];
      if (p)
	n2++;
      prohibited_[i]=p;
    }
    assert (n==n2);
    numberProhibited_=numberColumns;
  }
  if (rowType_) {
    const int * original = pinfo->originalRows();
    int numberRows = solver->getNumRows();
#ifdef COIN_DEVELOP
    int nMarked1=0;
    for (int i=0;i<pinfo->getNumRows();i++) {
      if (rowType_[i])
	nMarked1++;
    }
    int nMarked2=0;
    int k=-1;
    for (int i=0;i<numberRows;i++) {
      int iRow = original[i];
      if (iRow<i)
	abort();
      if (iRow<=k)
	abort();
      k=iRow;
      if (rowType_[iRow])
	nMarked2++;
    }
    if (nMarked1>nMarked2)
      printf("Marked rows reduced from %d to %d\n",
	     nMarked1,nMarked2);
#endif
    for (int i=0;i<numberRows;i++) {
      int iRow = original[i];
      rowType_[i]=rowType_[iRow];
    }
    numberRowType_=numberRows;
  }
}
/* Fix some of problem - returning new problem.
   Uses reduced costs.
   Optional signed character array
   1 always keep, -1 always discard, 0 use djs
   
*/
OsiSolverInterface * 
CglPreProcess::someFixed(OsiSolverInterface & model, 
                                 double fractionToKeep,
                                 bool fixContinuousAsWell,
                                 char * keep) const
{
  model.resolve();
  int numberColumns = model.getNumCols();
  OsiSolverInterface * newModel = model.clone();
  int i;
  const double * lower = model.getColLower();
  const double * upper = model.getColUpper();
  const double * solution = model.getColSolution();
  double * dj = CoinCopyOfArray(model.getReducedCost(),numberColumns);
  int * sort = new int[numberColumns];
  int number=0;
  int numberThrow=0;
  int numberContinuous=0;
  for (i=0;i<numberColumns;i++) {
    if (!model.isInteger(i)&&upper[i]>lower[i])
      numberContinuous++;
    if (model.isInteger(i)||fixContinuousAsWell) {
      if (keep) {
        if (keep[i]==1) {
          continue; // always keep
        } else if (keep[i]==-1) {
          // always fix
          dj[number]=-1.0e30;
          numberThrow++;
          sort[number++]=i;
          continue;
        }
      }
      double value = solution[i];
      if (value<lower[i]+1.0e-8) {
        dj[number]=-dj[i];
        sort[number++]=i;
      } else if (value>upper[number]-1.0e-8) {
        dj[number]=-dj[i];
        sort[number++]=i;
      }
    }
  }
  CoinSort_2(dj,dj+number,sort);
  int numberToFix = static_cast<int> (numberColumns *(1.0-fractionToKeep));
  if (!fixContinuousAsWell)
    numberToFix = static_cast<int> ((numberColumns-numberContinuous) *(1.0-fractionToKeep));
  numberToFix = CoinMax(numberToFix,numberThrow);
  numberToFix = CoinMin(number,numberToFix);
  printf("%d columns fixed\n",numberToFix);
  for (i=0;i<numberToFix;i++) {
    int iColumn = sort[i];
    double value = solution[iColumn];
    if (value<lower[iColumn]+1.0e-8) {
      newModel->setColUpper(iColumn,lower[iColumn]);
    } else if (value>upper[number]-1.0e-8) {
      newModel->setColLower(iColumn,lower[iColumn]);
    } else {
      // must be a throw away on - go to lower
      newModel->setColUpper(iColumn,lower[iColumn]);
    }
  }
  return newModel;
}
// If we have a cutoff - fix variables
int 
CglPreProcess::reducedCostFix(OsiSolverInterface & model)
{
  double cutoff ;
  model.getDblParam(OsiDualObjectiveLimit,cutoff) ;
  double direction = model.getObjSense() ;
  cutoff *= direction;
  double gap = cutoff - model.getObjValue()*direction ;
  double tolerance;
  model.getDblParam(OsiDualTolerance,tolerance) ;
  if (gap<=0.0||fabs(cutoff)>1.0e20)
    return 0;
  gap += 100.0*tolerance;
  // not really but thats all we can get
  double integerTolerance;
  model.getDblParam(OsiPrimalTolerance,integerTolerance) ;
  int numberColumns = model.getNumCols();

  const double *lower = model.getColLower() ;
  const double *upper = model.getColUpper() ;
  const double *solution = model.getColSolution() ;
  const double *reducedCost = model.getReducedCost() ;

  int numberFixed = 0 ;
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (model.isInteger(iColumn)&&upper[iColumn]>lower[iColumn]) {
      double djValue = direction*reducedCost[iColumn] ;
      if (solution[iColumn] < lower[iColumn]+integerTolerance && djValue > gap) {
	model.setColUpper(iColumn,lower[iColumn]) ;
	numberFixed++ ; 
      } else if (solution[iColumn] > upper[iColumn]-integerTolerance && -djValue > gap) {
	model.setColLower(iColumn,upper[iColumn]) ;
	numberFixed++ ;
      }
    }
  }
  return numberFixed;
}
// Pass in prohibited columns 
void 
CglPreProcess::passInProhibited(const char * prohibited,int numberColumns)
{
  delete [] prohibited_;
  prohibited_ = CoinCopyOfArray(prohibited,numberColumns);
  numberProhibited_ = numberColumns;
}
/* Pass in row types
   0 normal
   1 cut rows - will be dropped if remain in
   At end of preprocess cut rows will be dropped
   and put into cuts
*/
void 
CglPreProcess::passInRowTypes(const char * rowTypes,int numberRows)
{
  delete [] rowType_;
  rowType_ = CoinCopyOfArray(rowTypes,numberRows);
  numberRowType_ = numberRows;
  cuts_ = CglStored();
}
// Make continuous variables integer
void 
CglPreProcess::makeInteger()
{
  // First check if we need to
  int numberInteger=0;
  {
    const double *lower = startModel_->getColLower() ;
    const double *upper = startModel_->getColUpper() ;
    int numberColumns = startModel_->getNumCols() ;
    int iColumn;
    int numberContinuous=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8) {
	if(startModel_->isInteger(iColumn)) 
	  numberInteger++;
	else
	  numberContinuous++;
      }
    }
    if (!numberContinuous)
      return;
  }
  OsiSolverInterface * solver = startModel_->clone();
  const double *objective = solver->getObjCoefficients() ;
  const double *lower = solver->getColLower() ;
  const double *upper = solver->getColUpper() ;
  int numberColumns = solver->getNumCols() ;
  int numberRows = solver->getNumRows();
  double direction = solver->getObjSense();
  int iRow,iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix  matrixByCol(*solver->getMatrixByCol());
  const double * element = matrixByCol.getElements();
  const int * row = matrixByCol.getIndices();
  const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
  const int * columnLength = matrixByCol.getVectorLengths();

  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();

  char * ignore = new char [numberRows];
  int * changed = new int[numberColumns];
  int * which = new int[numberRows];
  double * changeRhs = new double[numberRows];
  memset(changeRhs,0,numberRows*sizeof(double));
  memset(ignore,0,numberRows);
  int numberChanged=0;
  bool finished=false;
  while (!finished) {
    int saveNumberChanged = numberChanged;
    for (iRow=0;iRow<numberRows;iRow++) {
      int numberContinuous=0;
      double value1=0.0,value2=0.0;
      bool allIntegerCoeff=true;
      double sumFixed=0.0;
      int jColumn1=-1,jColumn2=-1;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (upper[jColumn] > lower[jColumn]+1.0e-8) {
          if (!solver->isInteger(jColumn)) {
            if (numberContinuous==0) {
              jColumn1=jColumn;
              value1=value;
            } else {
              jColumn2=jColumn;
              value2=value;
            }
            numberContinuous++;
          } else {
            if (fabs(value-floor(value+0.5))>1.0e-12)
              allIntegerCoeff=false;
          }
        } else {
          sumFixed += lower[jColumn]*value;
        }
      }
      double low = rowLower[iRow];
      if (low>-1.0e20) {
        low -= sumFixed;
        if (fabs(low-floor(low+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      double up = rowUpper[iRow];
      if (up<1.0e20) {
        up -= sumFixed;
        if (fabs(up-floor(up+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      if (!allIntegerCoeff)
        continue; // can't do
      if (numberContinuous==1) {
        // see if really integer
        // This does not allow for complicated cases
        if (low==up) {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // integer
              changed[numberChanged++]=jColumn1;
              solver->setInteger(jColumn1);
              if (upper[jColumn1]>1.0e20)
                solver->setColUpper(jColumn1,1.0e20);
              if (lower[jColumn1]<-1.0e20)
                solver->setColLower(jColumn1,-1.0e20);
            }
          }
        } else {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // This constraint will not stop it being integer
              ignore[iRow]=1;
            }
          }
        }
      } else if (numberContinuous==2) {
        if (low==up) {
          /* need general theory - for now just look at 2 cases -
             1 - +- 1 one in column and just costs i.e. matching objective
             2 - +- 1 two in column but feeds into G/L row which will try and minimize
          */
          if (fabs(value1)==1.0&&value1*value2==-1.0&&!lower[jColumn1]
              &&!lower[jColumn2]) {
            int n=0;
            int i;
            double objChange=direction*(objective[jColumn1]+objective[jColumn2]);
            double bound = CoinMin(upper[jColumn1],upper[jColumn2]);
            bound = CoinMin(bound,1.0e20);
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                which[n++]=jRow;
                changeRhs[jRow]=value;
              }
            }
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                if (!changeRhs[jRow]) {
                  which[n++]=jRow;
                  changeRhs[jRow]=value;
                } else {
                  changeRhs[jRow]+=value;
                }
              }
            }
            if (objChange>=0.0) {
              // see if all rows OK
              bool good=true;
              for (i=0;i<n;i++) {
                int jRow = which[i];
                double value = changeRhs[jRow];
                if (value) {
                  value *= bound;
                  if (rowLength[jRow]==1) {
                    if (value>0.0) {
                      double rhs = rowLower[jRow];
                      if (rhs>0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    } else {
                      double rhs = rowUpper[jRow];
                      if (rhs<0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    }
                  } else if (rowLength[jRow]==2) {
                    if (value>0.0) {
                      if (rowLower[jRow]>-1.0e20)
                        good=false;
                    } else {
                      if (rowUpper[jRow]<1.0e20)
                        good=false;
                    }
                  } else {
                    good=false;
                  }
                }
              }
              if (good) {
                // both can be integer
                changed[numberChanged++]=jColumn1;
                solver->setInteger(jColumn1);
                if (upper[jColumn1]>1.0e20)
                  solver->setColUpper(jColumn1,1.0e20);
                if (lower[jColumn1]<-1.0e20)
                  solver->setColLower(jColumn1,-1.0e20);
                changed[numberChanged++]=jColumn2;
                solver->setInteger(jColumn2);
                if (upper[jColumn2]>1.0e20)
                  solver->setColUpper(jColumn2,1.0e20);
                if (lower[jColumn2]<-1.0e20)
                  solver->setColLower(jColumn2,-1.0e20);
              }
            }
            // clear
            for (i=0;i<n;i++) {
              changeRhs[which[i]]=0.0;
            }
          }
        }
      }
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8&&!solver->isInteger(iColumn)) {
        double value;
        value = upper[iColumn];
        if (value<1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        value = lower[iColumn];
        if (value>-1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        bool integer=true;
        for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          if (!ignore[iRow]) {
            integer=false;
            break;
          }
        }
        if (integer) {
          // integer
          changed[numberChanged++]=iColumn;
          solver->setInteger(iColumn);
          if (upper[iColumn]>1.0e20)
            solver->setColUpper(iColumn,1.0e20);
          if (lower[iColumn]<-1.0e20)
            solver->setColLower(iColumn,-1.0e20);
        }
      }
    }
    finished = numberChanged==saveNumberChanged;
  }
  delete [] which;
  delete [] changeRhs;
  delete [] ignore;
  //increment=0.0;
  if (numberChanged) {
    handler_->message(CGL_MADE_INTEGER,messages_)
      <<numberChanged
      <<CoinMessageEol;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (solver->isInteger(iColumn)&&objective[iColumn])
	startModel_->setInteger(iColumn);
    }
  }
  delete solver;
  delete [] changed;
}
//#define BRON_TIMES
#ifdef BRON_TIMES
static int numberTimesX=0;
#endif
/* Replace cliques by more maximal cliques
   Returns NULL if rows not reduced by greater than cliquesNeeded*rows
   
*/
OsiSolverInterface * 
CglPreProcess::cliqueIt(OsiSolverInterface & model,
			double cliquesNeeded) const
{
  /*
    Initial arrays
    * Candidate nodes (columns)
    First nIn already in
    Next nCandidate are candidates
    numberColumns-1 back to nNot are Nots
    * Starts
    * Other node
    * Original row (paired with other node)
    * CliqueIn expanded array with 1 in, 2 not, 3 out, 0 possible, -1 never
    * Type (for original row)
    */
  const double *lower = model.getColLower() ;
  const double *upper = model.getColUpper() ;
  const double *rowLower = model.getRowLower() ;
  const double *rowUpper = model.getRowUpper() ;
  int numberRows = model.getNumRows() ;
  //int numberColumns = model.getNumCols() ;
  // Column copy of matrix
  //const double * element = model.getMatrixByCol()->getElements();
  //const int * row = model.getMatrixByCol()->getIndices();
  //const CoinBigIndex * columnStart = model.getMatrixByCol()->getVectorStarts();
  //const int * columnLength = model.getMatrixByCol()->getVectorLengths();
  // Row copy
  CoinPackedMatrix matrixByRow(*model.getMatrixByRow());
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();
  char * type = new char [numberRows];
  int numberElements=0;
  int numberCliques=0;
  for (int i=0;i<numberRows;i++) {
    type[i]=-1;
    if (rowUpper[i]!=1.0||
	(rowLower[i]>0.0&&rowLower[i]!=1.0))
      continue;
    bool possible = true;
    CoinBigIndex start = rowStart[i];
    CoinBigIndex end = start + rowLength[i];
    int n=0;
    for (CoinBigIndex j=start;j<end;j++) {
      int iColumn = column[j];
      if (upper[iColumn]==1.0&&lower[iColumn]==0.0&&
	  model.isInteger(iColumn)&&elementByRow[j]==1.0) {
	n++;
      } else {
	possible=false;
	break;
      }
    }
    // temp fix to get working faster for client
    if (rowLower[i]>0.0||n!=2)
      possible=false;
    if (possible) {
      numberElements+=n;
      numberCliques++;
      if (rowLower[i]>0.0)
	type[i]=1;
      else
	type[i]=0;
    }
  }
  OsiSolverInterface * newSolver = NULL;
  if (numberCliques>CoinMax(1,static_cast<int>(cliquesNeeded*numberRows))) {
#ifdef BRON_TIMES
    double time1 = CoinCpuTime();
#endif
    CglBK bk(model,type,numberElements);
    bk.bronKerbosch();
    newSolver = bk.newSolver(model);
#ifdef BRON_TIMES
    printf("Time %g - bron called %d times\n",CoinCpuTime()-time1,numberTimesX);
#endif
  }
  delete [] type;
  return newSolver;
}
// Default constructor
CglBK::CglBK()
{
  candidates_=NULL;
  mark_=NULL;
  start_=NULL;
  otherColumn_=NULL;
  originalRow_=NULL;
  dominated_=NULL;
  cliqueMatrix_=NULL;
  rowType_=NULL;
  numberColumns_=0;
  numberRows_=0;
  numberPossible_=0;
  numberCandidates_=0;
  firstNot_=0;
  numberIn_=0;
  left_=0;
  lastColumn_=0;
} 
  
// Useful constructor
CglBK::CglBK(const OsiSolverInterface & model, const char * rowType,
	     int numberElements)
{
  const double *lower = model.getColLower() ;
  const double *upper = model.getColUpper() ;
  const double *rowLower = model.getRowLower() ;
  const double *rowUpper = model.getRowUpper() ;
  numberRows_ = model.getNumRows() ;
  numberColumns_ = model.getNumCols() ;
  // Column copy of matrix
#ifndef NDEBUG
  const double * element = model.getMatrixByCol()->getElements();
#endif
  const int * row = model.getMatrixByCol()->getIndices();
  const CoinBigIndex * columnStart = model.getMatrixByCol()->getVectorStarts();
  const int * columnLength = model.getMatrixByCol()->getVectorLengths();
  start_ = new CoinBigIndex[numberColumns_+1];
  otherColumn_ = new int [numberElements];
  candidates_ = new int [2*numberColumns_];
  CoinZeroN(candidates_,2*numberColumns_); // for valgrind
  originalRow_ = new int [numberElements];
  dominated_ = new int [numberRows_];
  CoinZeroN(dominated_,numberRows_);
  numberElements=0;
  numberPossible_=0;
  rowType_=rowType;
  // Row copy
  CoinPackedMatrix matrixByRow(*model.getMatrixByRow());
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();
#if 1
  // take out duplicate doubleton rows
  double * sort = new double[numberRows_];
  int * which = new int [numberRows_];
  double * randomValues = new double [numberColumns_];
  // Initialize random seed
  CoinThreadRandom randomGenerator(987654321);
  for (int i=0;i<numberColumns_;i++)
    randomValues[i]=randomGenerator.randomDouble();
  int nSort=0;
  for (int i=0;i<numberRows_;i++) {
    if (rowLength[i]==2&&rowUpper[i]==1.0) {
      int first = rowStart[i];
      int last = first+1;
      if (column[first]>column[last]) {
	first=last;
	last=rowStart[i];
      }
      int iColumn1 = column[first];
      int iColumn2 = column[last];
      double value = elementByRow[first]*randomValues[iColumn1]+
	elementByRow[last]*randomValues[iColumn2];
      sort[nSort]=value;
      which[nSort++]=i;
    }
  }
  CoinSort_2(sort,sort+nSort,which);
  double value=sort[0];
  int nDominated=0;
  for (int i=1;i<nSort;i++) {
    if (sort[i]==value) {
      int i1=which[i-1];
      int i2=which[i];
      if (rowLower[i1]==rowLower[i2]) {
	int first1 = rowStart[i1];
	int last1 = first1+1;
	if (column[first1]>column[last1]) {
	  first1=last1;
	  last1=rowStart[i1];
	}
	int iColumn11 = column[first1];
	int iColumn12 = column[last1];
	int first2 = rowStart[i2];
	int last2 = first2+1;
	if (column[first2]>column[last2]) {
	  first2=last2;
	  last2=rowStart[i2];
	}
	int iColumn21 = column[first2];
	int iColumn22 = column[last2];
	if (iColumn11==iColumn21&&
	    iColumn12==iColumn22&&
	    elementByRow[first1]==elementByRow[first2]&&
	    elementByRow[last1]==elementByRow[last2]) {
	  dominated_[i2]=1;
	  nDominated++;
	}
      }
    }
    value=sort[i];
  }
  //if (nDominated)
  //printf("%d duplicate doubleton rows!\n",nDominated);
  delete [] randomValues;
  delete [] sort;
  delete [] which;
#endif
  for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
    start_[iColumn]=numberElements;
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex end = start + columnLength[iColumn];
    if (upper[iColumn]==1.0&&lower[iColumn]==0.0&&
	model.isInteger(iColumn)) {
      for (CoinBigIndex j=start;j<end;j++) {
	int iRow = row[j];
	if (rowType[iRow]>=0&&!dominated_[iRow]) {
	  assert(element[j]==1.0);
	  CoinBigIndex r=rowStart[iRow];
	  assert (rowLength[iRow]==2);
	  int kColumn = column[r];
	  if (kColumn==iColumn)
	    kColumn=column[r+1];
	  originalRow_[numberElements]=iRow;
	  otherColumn_[numberElements++]=kColumn;
	}
      }
      if (numberElements>start_[iColumn]) {
	candidates_[numberPossible_++]=iColumn;
      }
    }
  }
  start_[numberColumns_]=numberElements;
  numberCandidates_=numberPossible_;
  numberIn_=0;
  firstNot_ = numberPossible_;
  left_=numberPossible_;
  lastColumn_=-1;
  mark_ = new char [numberColumns_];
  memset(mark_,0,numberColumns_);
  cliqueMatrix_=new CoinPackedMatrix(false,0.5,0.0);
  int n=0;
  for (int i=0;i<numberRows_;i++) {
    if (rowType[i]>=0)
      n++;
  }
  cliqueMatrix_->reserve(CoinMin(100,n),5*numberPossible_);
} 
  
// Copy constructor .
CglBK::CglBK(const CglBK & rhs)
{
  // This only copies data in candidates_
  // rest just points
  candidates_ = CoinCopyOfArray(rhs.candidates_,2*rhs.numberPossible_);
  mark_=rhs.mark_;
  start_=rhs.start_;
  otherColumn_=rhs.otherColumn_;
  originalRow_=rhs.originalRow_;
  dominated_=rhs.dominated_;
  cliqueMatrix_=rhs.cliqueMatrix_;
  rowType_=rhs.rowType_;
  numberColumns_=rhs.numberColumns_;
  numberRows_=rhs.numberRows_;
  numberPossible_=rhs.numberPossible_;
  numberCandidates_=rhs.numberCandidates_;
  firstNot_=rhs.firstNot_;
  numberIn_=rhs.numberIn_;
  left_=rhs.left_;
  lastColumn_=rhs.lastColumn_;
} 

// Assignment operator 
CglBK & CglBK::operator=(const CglBK& rhs)
{
  if (this!=&rhs) {
    delete [] candidates_;
    // This only copies data in candidates_
    // rest just points
    candidates_ = CoinCopyOfArray(rhs.candidates_,2*numberPossible_);
    mark_=rhs.mark_;
    start_=rhs.start_;
    otherColumn_=rhs.otherColumn_;
    originalRow_=rhs.originalRow_;
    dominated_=rhs.dominated_;
    cliqueMatrix_=rhs.cliqueMatrix_;
    rowType_=rhs.rowType_;
    numberColumns_=rhs.numberColumns_;
    numberRows_=rhs.numberRows_;
    numberPossible_=rhs.numberPossible_;
    numberCandidates_=rhs.numberCandidates_;
    firstNot_=rhs.firstNot_;
    numberIn_=rhs.numberIn_;
    left_=rhs.left_;
    lastColumn_=rhs.lastColumn_;
  }
  return *this;
} 

// Destructor 
CglBK::~CglBK ()
{
  delete [] candidates_;
  // only deletes if left_==-1
  if (left_==-1) {
    delete [] mark_;
    delete [] start_;
    delete [] otherColumn_;
    delete [] originalRow_;
    delete [] dominated_;
    delete cliqueMatrix_;
  }
}
// For Bron-Kerbosch
void 
CglBK::bronKerbosch()
{
#ifdef BRON_TIMES
  numberTimesX++;
  if ((numberTimesX%1000)==0)
    printf("times %d - %d candidates left\n",numberTimesX,numberCandidates_);
#endif
  if (!numberCandidates_&&firstNot_==numberPossible_) {
    // mark original rows which are dominated
    // save if clique size >2
    if (numberIn_>2) {
      double * elements = new double [numberIn_];
      int * column = candidates_+numberPossible_;
      // mark those in clique
      for (int i=0;i<numberIn_;i++) {
	int iColumn=column[i];
	mark_[iColumn]=1;
      }
      for (int i=0;i<numberIn_;i++) {
	elements[i]=1.0;
	int iColumn=column[i];
	for (int j=start_[iColumn];j<start_[iColumn+1];j++) {
	  int jColumn = otherColumn_[j];
	  if (mark_[jColumn]) {
	    int iRow=originalRow_[j];
	    dominated_[iRow]++;
	  }
	}
      }
      for (int i=0;i<numberIn_;i++) {
	int iColumn=column[i];
	mark_[iColumn]=0;
      }
      cliqueMatrix_->appendRow(numberIn_,column,elements);
      delete [] elements;
    }
  } else {
#if 0
    int nCplusN=numberCandidates_+(numberPossible_-firstNot_);
    int iChoose = CoinDrand48()*nCplusN;
    iChoose=CoinMin(0,nCplusN-1);
    if (iChoose>=numberCandidates_) {
      iChoose -= numberCandidates_;
      iChoose += firstNot_;
    }
#else
    for (int i=0;i<numberCandidates_;i++) {
      int jColumn = candidates_[i];
      mark_[jColumn]=1;
    }
    int nMax=0;
    int iChoose=0;
    for (int i=numberPossible_-1;i>=firstNot_;i--) {
      int iColumn = candidates_[i];
      int n=0;
      for (int j=start_[iColumn];j<start_[iColumn+1];j++) {
	int jColumn = otherColumn_[j];
	n += mark_[jColumn];
      }
      if (n>nMax) {
	nMax=n;
	iChoose=i;
      } 
    }
    if (nMax<numberCandidates_-1||!nMax) {
      for (int i=0;i<numberCandidates_;i++) {
	int iColumn = candidates_[i];
	int n=0;
	for (int j=start_[iColumn];j<start_[iColumn+1];j++) {
	  int jColumn = otherColumn_[j];
	  n += mark_[jColumn];
	}
	if (n>nMax) {
	  nMax=n;
	  iChoose=i;
	}
      }
    }
    for (int i=0;i<numberCandidates_;i++) {
      int jColumn = candidates_[i];
      mark_[jColumn]=0;
    }
#endif
    iChoose = candidates_[iChoose];
    int * temp = candidates_+numberPossible_+numberIn_;
    int nTemp=0;
    if (nMax<numberCandidates_) {
      // Neighborhood of iColumn
      for (int j=start_[iChoose];j<start_[iChoose+1];j++) {
	int jColumn = otherColumn_[j];
	mark_[jColumn]=1;
      }
      for (int i=0;i<numberCandidates_;i++) {
	int jColumn = candidates_[i];
	if (!mark_[jColumn])
	  temp[nTemp++]=jColumn;
      }
      for (int j=start_[iChoose];j<start_[iChoose+1];j++) {
	int jColumn = otherColumn_[j];
	mark_[jColumn]=0;
      }
    }
    //if (nMax==numberCandidates_)
    //assert (!nTemp);
    for (int kk=0;kk<nTemp;kk++) {
      int iColumn=temp[kk];
      // move up
      int put=0;
      for (int i=0;i<numberCandidates_;i++) {
	if (candidates_[i]!=iColumn)  
	  candidates_[put++]=candidates_[i];
      }
      numberCandidates_--;
      CglBK bk2(*this);
      int * newCandidates=bk2.candidates_;
#if 0
      printf("%p (next %p) iColumn %d, %d candidates %d not %d in\n",
	     this,&bk2,iColumn,numberCandidates_,
	     numberPossible_-firstNot_,numberIn_);
      for (int i=0;i<numberCandidates_;i++) {
	printf(" %d",candidates_[i]);
      }
      printf("\n");
#endif
      newCandidates[numberPossible_+numberIn_]=iColumn;
      bk2.numberIn_=numberIn_+1;
      // Neighborhood of iColumn
      for (int j=start_[iColumn];j<start_[iColumn+1];j++) {
	int jColumn = otherColumn_[j];
	mark_[jColumn]=1;
      }
      // Intersection of candidates with neighborhood
      int numberCandidates=0;
      for (int i=0;i<bk2.numberCandidates_;i++) {
	int jColumn = newCandidates[i];
	if (mark_[jColumn])
	  newCandidates[numberCandidates++]=jColumn;
      }
      bk2.numberCandidates_=numberCandidates;
      // Intersection of not with neighborhood
      int firstNot=numberPossible_;
      for (int i=numberPossible_-1;i>=bk2.firstNot_;i--) {
	int jColumn = newCandidates[i];
	if (mark_[jColumn])
	  newCandidates[--firstNot]=jColumn;
      }
      bk2.firstNot_=firstNot;
      for (int j=start_[iColumn];j<start_[iColumn+1];j++) {
	int jColumn = otherColumn_[j];
	mark_[jColumn]=0;
      }
      //for (int i=0;i<numberColumns_;i++)
      //assert (!mark_[i]);
      bk2.bronKerbosch();
      candidates_[--firstNot_]=iColumn;
    }
  }
}
// Creates strengthened smaller model
OsiSolverInterface * 
CglBK::newSolver(const OsiSolverInterface & model)
{
  // See how many rows can be deleted
  int nDelete=0;
  int * deleted = new int [numberRows_];
  for (int i=0;i<numberRows_;i++) {
    if (dominated_[i]) {
      deleted[nDelete++]=i;
    }
  }
  int nAdd=cliqueMatrix_->getNumRows();
  printf ("%d rows can be deleted with %d new cliques\n",
	  nDelete,nAdd);

  OsiSolverInterface * newSolver = NULL;
  if (nDelete>nAdd) {
    newSolver = model.clone();
    newSolver->deleteRows(nDelete,deleted);
    double * lower = new double [nAdd];
    double * upper = new double [nAdd];
    for (int i=0;i<nAdd;i++) {
      lower[i]=-COIN_DBL_MAX;
      upper[i]=1.0;
    }
    const double * elementByRow = cliqueMatrix_->getElements();
    const int * column = cliqueMatrix_->getIndices();
    const CoinBigIndex * rowStart = cliqueMatrix_->getVectorStarts();
    //const int * rowLength = cliqueMatrix_->getVectorLengths();
    assert (cliqueMatrix_->getNumElements()==rowStart[nAdd]);
    newSolver->addRows(nAdd,rowStart,column,elementByRow,lower,upper);
    delete [] lower;
    delete [] upper;
  }
  delete [] deleted;
  // mark so everything will be deleted
  left_=-1;
  return newSolver;
}
