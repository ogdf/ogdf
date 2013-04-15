/*
 * $Revision: 3386 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-10 14:00:02 +0200 (Mi, 10. Apr 2013) $
 ***************************************************************/

/*!\file
 * \author Frank Baumann
 *
 * \par License:
 * This file is part of ABACUS - A Branch And CUt System
 * Copyright (C) 1995 - 2003
 * University of Cologne, Germany
 *
 * \par
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * \par
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * \see http://www.gnu.org/copyleft/gpl.html
 */

#include <ogdf/abacus/osiif.h>
#include <ogdf/abacus/master.h>
#include <ogdf/abacus/lpmasterosi.h>
#include <ogdf/abacus/lpsubosi.h>
#include <ogdf/abacus/column.h>
#include <ogdf/abacus/row.h>

namespace abacus {


OsiIF::OsiIF(Master *master)
	:
	LP(master),
	osiLP_(NULL),
	value_(0.),
	xVal_(NULL),
	barXVal_(NULL),
	reco_(NULL),
	yVal_(NULL),
	cStat_(NULL),
	rStat_(NULL),
	rhs_(NULL),
	rowactivity_(NULL),
	rowsense_(NULL),
	colupper_(NULL),
	collower_(NULL),
	objcoeff_(NULL),
	ws_(NULL)
{
	lpMasterOsi_ = master->lpMasterOsi();
}


OsiIF::OsiIF(
	Master *master,
	OptSense sense,
	int nRow,
	int maxRow,
	int nCol,
	int maxCol,
	Array<double> &obj,
	Array<double> &lb,
	Array<double> &ub,
	Array<Row*> &rows)
	:
LP(master),
	osiLP_(NULL),
	value_(0.),
	xVal_(NULL),
	barXVal_(NULL),
	reco_(NULL),
	yVal_(NULL),
	cStat_(NULL),
	rStat_(NULL),
	rhs_(NULL),
	rowactivity_(NULL),
	rowsense_(NULL),
	colupper_(NULL),
	collower_(NULL),
	objcoeff_(NULL),
	ws_(NULL)
{
	lpMasterOsi_ = master->lpMasterOsi();

	_initialize(sense, nRow, maxRow, nCol, maxCol, obj, lb, ub, rows);
}


OsiIF::~OsiIF()
{
	delete ws_;
	delete osiLP_;

	freeDouble(xVal_);
	freeDouble(yVal_);
	freeDouble(reco_);
	freeDouble(rowactivity_);
	freeChar(cStat_);
	freeChar(rStat_);
}


void OsiIF::_initialize(
	OptSense sense,
	int nRow,
	int maxRow,
	int nCol,
	int maxCol,
	Array<double> &obj,
	Array<double> &lBound,
	Array<double> &uBound,
	Array<Row*> &rows)
{
	osiLP_ = getDefaultInterface();
	currentSolverType_ = Exact;

	// switch off output from the solver
	// can be reset in setSolverParameters
	osiLP_->setHintParam(OsiDoReducePrint, true, OsiHintDo);
	osiLP_->messageHandler()->setLogLevel(0);
	master_->setSolverParameters(osiLP_, currentSolverType() == Approx);

	numRows_ = nRow;
	numCols_ = nCol;
	double *lbounds = new double[numCols_];
	double *ubounds = new double[numCols_];
	double *objectives = new double[numCols_];

	CoinPackedVector *coinrow = new CoinPackedVector();
	CoinPackedMatrix *matrix =  new CoinPackedMatrix(false,0,0);
	matrix->setDimensions(0, numCols_);

	for (int i = 0; i < numCols_; i++){
		lbounds[i] = lBound[i];
		ubounds[i] = uBound[i];
		objectives[i] = obj[i];
	}

	if (currentSolverType() == Exact && numRows_ == 0 && master_->defaultLpSolver() == Master::CPLEX) {
		loadDummyRow(osiLP_, lbounds, ubounds, objectives);
	}
	else {
		char *senses = new char[numRows_];
		double *rhs = new double[numRows_];
		double *ranges = new double[numRows_];

		for (int i = 0; i < numRows_; i++){
			coinrow->clear();
			for (int j = 0; j < rows[i]->nnz(); j++){
				coinrow->insert(rows[i]->support(j), rows[i]->coeff(j));
			}
			matrix->appendRow(*coinrow);
			senses[i] = csense2osi(rows[i]->sense());
			rhs[i] = rows[i]->rhs();
			ranges[i] = 0.0;
		}
		lpSolverTime_.start();
		osiLP_->loadProblem(*matrix, lbounds, ubounds, objectives, senses, rhs, ranges);
		lpSolverTime_.stop();

		freeChar(senses);
		freeDouble(rhs);
		freeDouble(ranges);
	}

	// set the sense of the optimization
	_sense(sense);

	// get the pointers to the solution, reduced costs etc.
	lpSolverTime_.start();
	numRows_ = osiLP_->getNumRows();
	numCols_ = osiLP_->getNumCols();
	rhs_ = osiLP_->getRightHandSide();
	rowsense_ = osiLP_->getRowSense();
	colupper_ = osiLP_->getColUpper();
	collower_ = osiLP_->getColLower();
	objcoeff_ = osiLP_->getObjCoefficients();
	if( ws_ != NULL )
		delete ws_;
	//ws_ = dynamic_cast<CoinWarmStartBasis *>(osiLP_->getWarmStart());
	ws_=0;

	xValStatus_ = recoStatus_ = yValStatus_ = slackStatus_ = basisStatus_ = Missing;
	lpSolverTime_.stop();

	delete coinrow;
	delete matrix;
	freeDouble(lbounds);
	freeDouble(ubounds);
	freeDouble(objectives);
}


void OsiIF::_loadBasis(
	Array<LPVARSTAT::STATUS> &lpVarStat,
	Array<SlackStat::STATUS> &slackStat)
{
	int lps = lpVarStat.size();
	int sls = slackStat.size();

	CoinWarmStartBasis *ws = NULL;
	ws = new CoinWarmStartBasis();
	ws->setSize(numCols_, numRows_);

	if (osiLP_->getNumCols() > lps) {
		Logger::ifout() << "OsiIF::_loadBasis: mismatch in number of columns: OSI " << osiLP_->getNumCols() << ", Abacus: " << lps << "\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
	for (int i = 0; i < numCols_; i++)
		ws->setStructStatus(i, lpVarStat2osi(lpVarStat[i]));

	if (osiLP_->getNumRows() > sls) {
		Logger::ifout() << "OsiIF::_loadBasis: mismatch in number of rows: OSI " << osiLP_->getNumCols() << ", Abacus: " << sls << "\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
	for (int i = 0; i < numRows_; i++)
		ws->setArtifStatus(i, slackStat2osi(slackStat[i]));

	lpSolverTime_.start();
	slackStatus_ = basisStatus_ = Missing;
	int status = 0;
	// FIXME loadBasis
	// better test whether the number of basic structurals is correct?
	if (ws->numberBasicStructurals() > 0) {
		status = osiLP_->setWarmStart(dynamic_cast<CoinWarmStart *> (ws));
		if (ws_ != NULL) delete ws_;
		ws_ = dynamic_cast<CoinWarmStartBasis*> (osiLP_->getWarmStart());
		if (ws_ != NULL) {
			delete[] cStat_;
			int nStructBytes = (int) ceil( ws_->getNumStructural() / 4.0);
			cStat_ = new char[nStructBytes];
			for(int i = 0; i < nStructBytes; i++) {
				cStat_[i] = ws_->getStructuralStatus()[i];
			}

			delete[] rStat_;
			int nArtBytes = (int) ceil( ws_->getNumArtificial() / 4.0 );
			rStat_ = new char[nArtBytes];
			for(int i = 0; i < nArtBytes; i++) {
				rStat_[i] = ws_->getArtificialStatus()[i];
			}

			basisStatus_ = Available;
		} else
			basisStatus_ = Missing;
	} else
		status = 2;
	lpSolverTime_.stop();

	delete ws;

	if (status == 0) {
		Logger::ifout()
			<< "OsiIF::_loadBasis(): loading the new basis has failed. Status "
			<< status << endl;
		// FIXME loadBasis
		return;
	} else
		return;
}


int OsiIF::_getInfeas(int &infeasRow, int &infeasCol, double *bInvRow) const
{
	// This is only implemented in CPLEX and rarely used.
	Logger::ifout() << "OsiIF::_getInfeas(): currently not available\n";
	OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
}


void OsiIF::_remRows(ArrayBuffer<int> &ind)
{
	const int n = ind.size();
	int *indices = new int[n];

	for (int i = 0; i < n; i++) {
		indices[i] = ind[i];
	}

	lpSolverTime_.start();
	osiLP_->deleteRows(n, indices);
	numRows_ = osiLP_->getNumRows();
	rowsense_ = osiLP_->getRowSense();
	rhs_ = osiLP_->getRightHandSide();
	lpSolverTime_.stop();

	freeInt(indices);
}


void OsiIF::_addRows(ArrayBuffer<Row*> &rows)
{
	CoinPackedVector *coinrow = new CoinPackedVector();

	for (int i = 0; i < rows.size(); i++) {
		coinrow->clear();
		for (int j = 0; j < rows[i]->nnz(); j++){
			coinrow->insert(rows[i]->support(j), rows[i]->coeff(j));
		}
		lpSolverTime_.start();
		osiLP_->addRow(*coinrow, csense2osi(rows[i]->sense()), rows[i]->rhs(), 0.0);
		lpSolverTime_.stop();
	}

	delete coinrow;
	lpSolverTime_.start();
	numRows_ = osiLP_->getNumRows();
	rhs_ = osiLP_->getRightHandSide();
	numCols_ = osiLP_->getNumCols();
	collower_ = osiLP_->getColLower();
	colupper_ = osiLP_->getColUpper();
	lpSolverTime_.stop();

}


void OsiIF::_remCols(ArrayBuffer<int> &vars)
{
	int num = vars.size();
	int *indices = new int[num];

	for (int i = 0; i < num; i++)
		indices[i] = vars[i];

	lpSolverTime_.start();
	osiLP_->deleteCols(num, indices);
	numCols_ = osiLP_->getNumCols();
	collower_ = osiLP_->getColLower();
	colupper_ = osiLP_->getColUpper();
	objcoeff_ = osiLP_->getObjCoefficients();
	lpSolverTime_.stop();

	freeInt(indices);
}


void OsiIF::_addCols(ArrayBuffer<Column*> &newCols)
{
	int num;
	double ub, lb, obj;
	int  *supports; //!< supports of added rows
	double  *coeffs; //!< coefficients of added rows
	CoinPackedVector *newcol = new CoinPackedVector;

	for (int i = 0; i < newCols.size(); i++) {
		num = newCols[i]->nnz();
		ub =  newCols[i]->uBound();
		lb =  newCols[i]->lBound();
		obj =  newCols[i]->obj();
		supports = new int[num];
		coeffs = new double[num];

		for (int j = 0; j < num; j++) {
			supports[j] = newCols[i]->support(j);
			coeffs[j] = newCols[i]->coeff(j);
		}

		newcol->setVector(num, supports, coeffs);
		lpSolverTime_.start();
		osiLP_->addCol(*newcol, lb, ub, obj);
		lpSolverTime_.stop();

		freeInt(supports);
		freeDouble(coeffs);
	}

	lpSolverTime_.start();
	numCols_ = osiLP_->getNumCols();
	collower_ = osiLP_->getColLower();
	colupper_ = osiLP_->getColUpper();
	objcoeff_ = osiLP_->getObjCoefficients();
	lpSolverTime_.stop();
	delete newcol;
}


void OsiIF::_changeRhs(Array<double> &newRhs)
{
	lpSolverTime_.start();

	for (int i = 0; i < newRhs.size(); i++)
		osiLP_->setRowType(i, rowsense_[i], newRhs[i], 0);

	rowsense_ = osiLP_->getRowSense();
	rhs_ = osiLP_->getRightHandSide();
	lpSolverTime_.stop();
}


void OsiIF::_changeLBound(int i, double newLb)
{
	lpSolverTime_.start();

	osiLP_->setColLower(i, newLb);
	collower_ = osiLP_->getColLower();

	lpSolverTime_.stop();
}


void OsiIF::_changeUBound(int i, double newUb)
{
	lpSolverTime_.start();

	osiLP_->setColUpper(i, newUb);
	colupper_ = osiLP_->getColUpper();

	lpSolverTime_.stop();
}


int OsiIF::_pivotSlackVariableIn(ArrayBuffer<int> &rows)
{
	Logger::ifout() << "OsiIF:::_pivotSlackVariableIn(): currently not implemented" << endl;
	return 1;
}


LP::OPTSTAT OsiIF::_primalSimplex()
{
	lpSolverTime_.start();

	// switch the interface, if necessary
	if (currentSolverType_ != Exact){
		currentSolverType_ = Exact;
		osiLP_ = switchInterfaces(Exact);
	}
	osiLP_->initialSolve();

	lpSolverTime_.stop();

	// check for solver statuses
	if (osiLP_->isAbandoned()){
		Logger::ifout() << "OsiIF::_primalSimplex():\nWarning: solver Interface reports status isAbandoned\nThere have been numerical difficulties, aborting...\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}

	// get information about the solution
	getSol();

	// The order is important here
	if (osiLP_->isProvenOptimal()) return Optimal;
	if (osiLP_->isProvenPrimalInfeasible()) return Infeasible;
	if (osiLP_->isProvenDualInfeasible()) return Unbounded;
	if (osiLP_->isIterationLimitReached()) return LimitReached;
	else {
		Logger::ifout() << "OsiIF::_primalSimplex():\nunable to determine status of LP, aborting...\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


LP::OPTSTAT OsiIF::_dualSimplex()
{
	lpSolverTime_.start();

	// switch the interface, if necessary
	if (currentSolverType_ != Exact){
		currentSolverType_ = Exact;
		osiLP_ = switchInterfaces(Exact);
	}
	osiLP_->resolve();

	lpSolverTime_.stop();

	// check for solver statuses
	if (osiLP_->isAbandoned()){
		Logger::ifout() << "OsiIF::_dualSimplex():\nWarning: solver Interface reports staus isAbandoned\nThere have been numerical difficulties, aborting...\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}

	// get information about the solution
	getSol();

	// The order is important here
	if (osiLP_->isProvenOptimal()) return Optimal;
	if (osiLP_->isProvenPrimalInfeasible()) return Infeasible;
	if (osiLP_->isProvenDualInfeasible()) return Unbounded;
	if (osiLP_->isIterationLimitReached()) return LimitReached;
	else {
		Logger::ifout() << "OsiIF::_dualSimplex():\nunable to determine status of LP, aborting...\n";
		//FIXME what about strong branching?
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


LP::OPTSTAT OsiIF::_barrier(bool doCrossover)
{
	// The barrier method is not implemented in Osi.
	// We use the primal simplex method instead
	Logger::ifout() << "OsiIF::_barrier: Sorry, Osi implements no barrier method." << endl;
	Logger::ifout() << "Using primal simplex method instead." << endl;
	return _primalSimplex();
}


LP::OPTSTAT OsiIF::_approx()
{
	lpSolverTime_.start();

	// switch the interface, if necessary
	if (currentSolverType() != Approx){
		currentSolverType_ = Approx;
		osiLP_ = switchInterfaces(Approx);
	}
	osiLP_->resolve();

	lpSolverTime_.stop();

	// check for solver statuses
	if (osiLP_->isAbandoned()){
		Logger::ifout() << "OsiIF::_approx():\nWarning: solver Interface reports staus isAbandoned\nThere have been numerical difficulties, aborting...\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}

	// get information about the solution
	getSol();

	// The order is important here
	if (osiLP_->isProvenOptimal()) return Optimal;
	if (osiLP_->isProvenPrimalInfeasible()) return Infeasible;
	if (osiLP_->isProvenDualInfeasible()) return Unbounded;
	if (osiLP_->isIterationLimitReached()) return LimitReached;
	else {
		Logger::ifout() << "OsiIF::_approx(): ";
		Logger::ifout() << "unable to determine status of LP, assume the solution is optimal..." << endl;
		return Optimal;
	}
}


OptSense OsiIF::_sense() const
{
	OptSense s;

	double osiSense = osiLP_->getObjSense();

	if (osiSense == 1.)
		s.sense(OptSense::Min);
	else
		s.sense(OptSense::Max);

	return s;
}

void OsiIF::_sense(const OptSense &newSense)
{
	if (newSense.unknown()) {
		Logger::ifout() << "OsiIF::_sense: The objective sense can not be set to 'unknown' with OSI.\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
	else if (newSense.max())
		osiLP_->setObjSense(-1.);
	else
		osiLP_->setObjSense(1.);
}


void OsiIF::_row(int i, Row &r) const
{
	const CoinPackedMatrix* coinMatrix;
	CoinPackedVector coinVector;
	int coinNumEl;
	const int* coinIndices;
	const double* coinElements;

	coinMatrix = osiLP_->getMatrixByRow();
	coinVector = coinMatrix->getVector(i);
	coinNumEl = coinVector.getNumElements();
	coinIndices = coinVector.getIndices();
	coinElements = coinVector.getElements();

	r.clear();

	for (int j = 0; j < coinNumEl; j++)
		r.insert(coinIndices[j], coinElements[j]);
	r.sense(osi2csense(rowsense_[i]));
	r.rhs(_rhs(i));
}


double OsiIF::_barXVal(int i) const
{
	// The barrier algorithm is only supported by cplex
	// support may be added later
	Logger::ifout() << "OsiIF::_barXVal: The barrier algorithm is currently not supported\n";
	OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	//return barXValStatus_ ? 0. : barXVal_[i];
}


double OsiIF::_slack(int i) const
{
	double rhs = rhs_[i];
	double ract = rowactivity_[i];
	switch (rowsense_[i]){
	case 'L':
		return rhs - ract;
	case 'G':
		return ract - rhs;
	case 'E':
		return 0.0;
	default:
		Logger::ifout() << "OsiIF::_slack : slack not defined for sense " << rowsense_[i] << " for row " << i << " of " << osiLP_->getNumRows() << " osiLP_->getNumRows()\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


void OsiIF::getSol()
{
#ifdef OSIDEBUG
	Logger::ifout() << "<> isProvenOptimal: " << osiLP_->isProvenOptimal() << endl
	 << "<> isProvenDualInfeasible: " << osiLP_->isProvenDualInfeasible() << endl
	 << "<> isProvenPrimalInfeasible: " << osiLP_->isProvenPrimalInfeasible() << endl
	 << "<> isProvenDualInfeasible: " << osiLP_->isProvenDualInfeasible() << endl
	 << "<> isAbandoned: " << osiLP_->isAbandoned() << endl
	 << "<> isPrimalObjetiveLimitReached: " << osiLP_->isPrimalObjectiveLimitReached() << endl
	 << "<> isDualObjetiveLimitReached: " << osiLP_->isDualObjectiveLimitReached() << endl
	 << "<> isIterationLimitReached: " << osiLP_->isIterationLimitReached() << endl;
	double POL;
	osiLP_->getDblParam(OsiPrimalObjectiveLimit, POL);
	Logger::ifout() << "<> PrimalObjectiveLimit: " << POL << endl
	 << "<> ObjectiveSense: " << osiLP_->getObjSense() << endl
	 << "<> ObjectiveValue: " << osiLP_->getObjValue() << endl;
#endif

	lpSolverTime_.start();
	// get the solution
	xValStatus_ = recoStatus_ = yValStatus_ = slackStatus_ = basisStatus_ = Missing;

	numCols_ = osiLP_->getNumCols();
	numRows_ = osiLP_->getNumRows();
	collower_ = osiLP_->getColLower();
	colupper_ = osiLP_->getColUpper();
	objcoeff_ = osiLP_->getObjCoefficients();
	rhs_ = osiLP_->getRightHandSide();
	rowsense_ = osiLP_->getRowSense();

	if( !osiLP_->isProvenPrimalInfeasible() && !osiLP_->isAbandoned() ) {
		delete[] xVal_;
		xVal_ = new double[numCols_];
		for(int i=0; i < numCols_; i++){
			xVal_[i] = osiLP_->getColSolution()[i];
		}
		xValStatus_ = Available;

		delete[] rowactivity_;
		rowactivity_ = new double[numRows_];
		for(int j=0; j < numRows_; j++){
			rowactivity_[j] = osiLP_->getRowActivity()[j];
		}
	}

	if( !osiLP_->isProvenDualInfeasible() && !osiLP_->isAbandoned() ) {
		delete[] yVal_;
		yVal_ = new double[numRows_];
		for(int j=0; j < numRows_; j++){
			yVal_[j] = osiLP_->getRowPrice()[j];
		}
		yValStatus_ = Available;
	}

	if (osiLP_->isProvenOptimal() || osiLP_->isIterationLimitReached()) {
		value_ = osiLP_->getObjValue();
		delete[] reco_;
		reco_ = new double[numCols_];
		for(int i=0; i < numCols_; i++){
			reco_[i] = osiLP_->getReducedCost()[i];
		}
		recoStatus_ = Available;
		// get information about the basis
		if (currentSolverType() != Approx) {
			if (ws_ != NULL) delete ws_;
			ws_ = dynamic_cast<CoinWarmStartBasis*> (osiLP_->getWarmStart());
			delete[] cStat_;
			int nStructBytes = (int) ceil( ws_->getNumStructural() / 4.0);
			cStat_ = new char[nStructBytes];
			for(int i=0; i < nStructBytes; i++){
				cStat_[i] = ws_->getStructuralStatus()[i];
			}

			delete[] rStat_;
			int nArtBytes = (int) ceil( ws_->getNumArtificial() / 4.0 );
			rStat_ = new char[nArtBytes];
			for(int i=0; i < nArtBytes; i++){
				rStat_[i] = ws_->getArtificialStatus()[i];
			}

			basisStatus_ = Available;
			slackStatus_ = Available;
		} else {
			// when the solver is not exact all variables are assumed to be non-basic
			// this makes all variables candidates for fixing, that are at one of their bounds
			if (ws_ != NULL) delete ws_;
			ws_ = new CoinWarmStartBasis();
			ws_->setSize(numCols_, numRows_);
			for (int i = 0; i < numCols_; i++) {
				if (_uBound(i) - _xVal(i) < master_->eps()){
					ws_->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
				} else if (_xVal(i) - _lBound(i) < master_->eps()){
					ws_->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
				} else {
					ws_->setStructStatus(i, CoinWarmStartBasis::isFree);
				}
			}
			delete[] cStat_;
			int nStatBytes = (int) ceil(numCols_ / 4.0);
			cStat_ = new char[nStatBytes];
			for(int i=0; i < nStatBytes; i++){
				cStat_[i] = ws_->getStructuralStatus()[i];
			}
			basisStatus_ = Available;
		}
	}

	lpSolverTime_.stop();
}


void OsiIF::_rowRealloc(int newSize)
{
	// Memory management is completely handled by Osi
	//Logger::ifout() << "OsiIF::_rowRealloc : This function does not do anything. Memory management is completely handled by Osi." << endl;
}


void OsiIF::_colRealloc(int newSize)
{
	// Memory management is completely handled by Osi
	//Logger::ifout() << "OsiIF::_colRealloc : This function does not do anything. Memory management is completely handled by Osi." << endl;
}


LPVARSTAT::STATUS OsiIF::osi2lpVarStat(CoinWarmStartBasis::Status stat) const
{
	switch (stat) {
	case CoinWarmStartBasis::isFree:
		return LPVARSTAT::NonBasicFree;
	case CoinWarmStartBasis::basic:
		return LPVARSTAT::Basic;
	case CoinWarmStartBasis::atUpperBound:
		return LPVARSTAT::AtUpperBound;
	case CoinWarmStartBasis::atLowerBound:
		return LPVARSTAT::AtLowerBound;
	default:
		Logger::ifout() << "OsiIF::osi2lpVarStat( " << stat << " ) unknown status\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


CoinWarmStartBasis::Status OsiIF::slackStat2osi(SlackStat::STATUS stat) const
{
	switch(stat) {
	case SlackStat::NonBasicZero:
		return CoinWarmStartBasis::atLowerBound;
	case SlackStat::Basic:
		return CoinWarmStartBasis::basic;
	case SlackStat::NonBasicNonZero:
		return CoinWarmStartBasis::atUpperBound;
	case SlackStat::Unknown:
		return CoinWarmStartBasis::atLowerBound;
	default:
		Logger::ifout() << "OsiIF::slackStat2osi( " << stat << " ) corresponding OSI status unknown\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


SlackStat::STATUS OsiIF::osi2slackStat(CoinWarmStartBasis::Status stat) const
{
	switch (stat) {
	case CoinWarmStartBasis::atLowerBound:
		return SlackStat::NonBasicZero;
	case CoinWarmStartBasis::atUpperBound:
		return SlackStat::NonBasicNonZero;
	case CoinWarmStartBasis::basic:
		return SlackStat::Basic;
	case CoinWarmStartBasis::isFree:
	default:
		Logger::ifout() << "OsiIF::osi2slackStat( " << stat << " ) unknown status\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


LPVARSTAT::STATUS OsiIF::_lpVarStat(int i) const
{
	return basisStatus_ ? LPVARSTAT::Unknown : osi2lpVarStat(getStatus(cStat_, i));
}


SlackStat::STATUS OsiIF::_slackStat(int i) const
{
	return slackStatus_ ? SlackStat::Unknown : osi2slackStat(getStatus(rStat_, i));
}


char OsiIF::csense2osi(CSense *sense) const
{
	switch(sense->sense()) {
	case CSense::Less:  return 'L';
	case CSense::Equal:  return 'E';
	case CSense::Greater: return 'G';
	default:
		Logger::ifout() << "OsiIF::csense2osi unknown sense\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


CSense::SENSE OsiIF::osi2csense(char sense) const
{
	switch(sense) {
	case 'L': return CSense::Less;
	case 'E': return CSense::Equal;
	case 'G': return CSense::Greater;
	default:
		Logger::ifout() << "OsiIF::osi2csense( " << sense << " ) unknown sense";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


CoinWarmStartBasis::Status OsiIF::lpVarStat2osi(LPVARSTAT::STATUS stat) const
{
	// Status Unknown <--> atlowerBound
	switch (stat) {
	case LPVARSTAT::AtLowerBound:
		return CoinWarmStartBasis::atLowerBound;
	case LPVARSTAT::AtUpperBound:
		return CoinWarmStartBasis::atUpperBound;
	case LPVARSTAT::Basic:
		return CoinWarmStartBasis::basic;
	case LPVARSTAT::NonBasicFree:
		return CoinWarmStartBasis::isFree;
	case LPVARSTAT::Unknown:
		return CoinWarmStartBasis::atLowerBound;
	default:
		Logger::ifout() << "OsiIF::lpVarStat2osi( " << stat << " ) unknown status\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}
}


OsiSolverInterface* OsiIF::getDefaultInterface()
{
	OsiSolverInterface *interface=NULL;
	switch(master_->defaultLpSolver()) {
#ifdef OSI_CBC
	case Master::Cbc:
		interface = new OsiCbcSolverInterface;
		break;
#endif
#ifdef OSI_CLP
	case Master::Clp:
		interface = new OsiClpSolverInterface;
		break;
#endif
#ifdef OSI_CPX
	case Master::CPLEX:
		interface = new OsiCpxSolverInterface;
		break;
#endif
#ifdef OSI_DYLP
	case Master::DyLP:
		interface = new OsiDylpSolverInterface;
		break;
#endif
#ifdef OSI_FORTMP
	case Master::FortMP:
		interface = new OsiFmpSolverInterface;
		break;
#endif
#ifdef OSI_GLPK
	case Master::GLPK:
		interface =  new OsiGlpkSolverInterface;
		break;
#endif
#ifdef OSI_MOSEK
	case Master::MOSEK:
		interface =  new OsiMskSolverInterface;
		break;
#endif
#ifdef OSI_OSL
	case Master::OSL:
		interface = new OsiOslSolverInterface;
		break;
#endif
#ifdef OSI_SOPLEX
	case Master::SoPlex:
		interface = new OsiSpxSolverInterface;
		break;
#endif
#ifdef OSI_SYM
	case Master::SYMPHONY:
		interface = new OsiSymSolverInterface;
		break;
#endif
#ifdef OSI_XPRESS
	case Master::XPRESS_MP:
		interface = new OsiXprSolverInterface;
		break;
#endif
#ifdef OSI_GRB
	case Master::Gurobi:
		interface = new OsiGrbSolverInterface;
		break;
#endif
#ifdef OSI_CSDP
	case Master::Csdp:
		interface = new OsiCsdpSolverInterface;
		break;
#endif
	default:
		Logger::ifout() << "No support for solver " << Master::OSISOLVER_[master_->defaultLpSolver()] << " in Coin-Osi! (see defaultLP-Solver)\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
	}

	if(interface != NULL){
		interface->setHintParam(OsiDoDualInInitial, false, OsiHintDo);
		interface->setHintParam(OsiDoDualInResolve, true, OsiHintDo);
	}

	return interface;
}


OsiSolverInterface* OsiIF::switchInterfaces(SOLVERTYPE newMethod)
{
	OsiSolverInterface *s2 = NULL;
	if( newMethod == Exact )
	{
		s2 = getDefaultInterface();
	}
	else
	{
#ifdef OSI_VOL
		// TODO switchInterfaces
		// s2 = getApproxInterface()
		s2 = new OsiVolSolverInterface;
#else
		Logger::ifout() << "ABACUS has not been compiled with support for the Volume Algorithm, cannot switch to approximate solver.\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcOsiIf);
#endif
	}

	s2->setHintParam(OsiDoReducePrint, true, OsiHintDo);
	s2->messageHandler()->setLogLevel(0);
	master_->setSolverParameters(s2, currentSolverType() == Approx);

	if (currentSolverType() == Exact && numRows_ == 0 && master_->defaultLpSolver() == Master::CPLEX) {
		loadDummyRow(s2, osiLP_->getColLower(), osiLP_->getColUpper(), osiLP_->getObjCoefficients());
	}
	else {
		s2->loadProblem(*osiLP_->getMatrixByCol(), osiLP_->getColLower(),
			osiLP_->getColUpper(), osiLP_->getObjCoefficients(),
			osiLP_->getRowLower(), osiLP_->getRowUpper());
	}

	s2->setObjSense(osiLP_->getObjSense());

	delete osiLP_;

	// get the pointers to the solution, reduced costs etc.
	rhs_ = s2->getRightHandSide();
	rowsense_ = s2->getRowSense();
	colupper_ = s2->getColUpper();
	collower_ = s2->getColLower();
	objcoeff_ = s2->getObjCoefficients();
	if( ws_ != NULL )
		delete ws_;
	ws_ = dynamic_cast<CoinWarmStartBasis* >(s2->getWarmStart());

	xValStatus_ = recoStatus_= yValStatus_ = slackStatus_ = basisStatus_ = Missing;
	return s2;
}


void OsiIF::loadDummyRow(OsiSolverInterface* s2, const double* lbounds, const double* ubounds, const double* objectives)
{
	CoinPackedVector *coinrow = new CoinPackedVector();
	CoinPackedMatrix *matrix =  new CoinPackedMatrix(false,0,0);
	matrix->setDimensions(0, numCols_);
	ArrayBuffer<int> dummy(1,false);
	dummy.push(0);
	char *senses = new char[1];
	double *rhs = new double[1];
	double *ranges = new double[1];
	coinrow->insert(0, 1.);
	matrix->appendRow(*coinrow);
	senses[0] = 'E';
	rhs[0] = 1.;
	ranges[0] = 0.;

	lpSolverTime_.start();
	s2->loadProblem(*matrix, lbounds, ubounds, objectives, senses, rhs, ranges);
	lpSolverTime_.stop();
	_remRows(dummy);

	delete coinrow;
	delete matrix;
	freeChar(senses);
	freeDouble(rhs);
	freeDouble(ranges);

	return;
}
} //namespace abacus
