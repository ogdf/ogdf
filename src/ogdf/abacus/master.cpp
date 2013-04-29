/*
 * $Revision: 3452 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-24 13:50:34 +0200 (Mi, 24. Apr 2013) $
 ***************************************************************/

/*!\file
* \author Matthias Elf
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
*
* $Id: master.cc,v 2.19 2009-05-13 14:17:37 baumann Exp $
*/

#include <ogdf/abacus/master.h>
#include <ogdf/abacus/sub.h>
#include <ogdf/abacus/lpsub.h>
#include <ogdf/abacus/lpvarstat.h>
#include <ogdf/abacus/history.h>
#include <ogdf/abacus/opensub.h>
#include <ogdf/abacus/fixcand.h>
#include <ogdf/abacus/setbranchrule.h>
#include <ogdf/abacus/standardpool.h>

namespace abacus {

const char* Master::STATUS_[] = {
	"Optimal", "Error", "OutOfMemory", "Unprocessed", "Processing",
	"Guaranteed", "MaxLevel", "MaxCpuTime", "MaxNSub",
	"MaxCowTime", "ExceptionFathom"
};

const char * Master::ENUMSTRAT_[] = {
	"BestFirst", "BreadthFirst",
	"DepthFirst", "DiveAndBest"
};

const char * Master::BRANCHINGSTRAT_[] = {
	"CloseHalf", "CloseHalfExpensive"
};

const char* Master::PRIMALBOUNDMODE_[] = {
	"None", "Optimum", "OptimumOne"
};

const char* Master::SKIPPINGMODE_[] = {
	"SkipByNode", "SkipByLevel"
};

const char* Master::CONELIMMODE_[] = {
	"None", "NonBinding", "Basic"
};

const char* Master::VARELIMMODE_[] = {
	"None", "ReducedCost"
};

const char* Master::VBCMODE_[] = {
	"None", "File", "Pipe"
};

// Possible values for the defaultLpSolver_ parameter
const char* Master::OSISOLVER_[] = {
	"Cbc", "Clp", "CPLEX", "DyLP", "FortMP",
	"GLPK", "MOSEK", "OSL", "SoPlex",
	"SYMPHONY", "XPRESS_MP", "Gurobi", "Csdp"
};


Master::Master(
	const char *problemName,
	bool cutting,
	bool pricing,
	OptSense::SENSE optSense,
	double eps,
	double machineEps,
	double infinity,
	bool readParamFromFile)
	:
	AbacusGlobal(eps, machineEps, infinity),
	problemName_(problemName),
	readParamFromFile_(readParamFromFile),
	optSense_(optSense),
	root_(0),
	rRoot_(0),
	openSub_(0),
	history_(0),
	enumerationStrategy_(BestFirst),
	branchingStrategy_(CloseHalfExpensive),
	nBranchingVariableCandidates_(1),
	nStrongBranchingIterations_(50),
	defaultLpSolver_(
#if defined(COIN_OSI_CPX)
	CPLEX
#elif defined(COIN_OSI_SYM)
	SYMPHONY
#elif defined(COIN_OSI_GRB)
	Gurobi
#else
	Clp
#endif
	),
	lpMasterOsi_(0),
	conPool_(0),
	cutPool_(0),
	varPool_(0),
	fixCand_(0),
	cutting_(cutting),
	pricing_(pricing),
	solveApprox_(false),
	nSubSelected_(0),
	VbcLog_(NoVbc),
	treeStream_(0),
	requiredGuarantee_(0.0),
	maxLevel_(numeric_limits<int>::max()),
	maxNSub_(numeric_limits<int>::max()),
	maxCpuTime_(__int64(999999)*3600+59*60+59),
	maxCowTime_(__int64(999999)*3600+59*60+59),
	objInteger_(false),
	tailOffNLp_(0),
	tailOffPercent_(0.000001),
	dbThreshold_(0),
	minDormantRounds_(1),
	pbMode_(NoPrimalBound),
	pricingFreq_(0),
	skipFactor_(1),
	skippingMode_(SkipByNode),
	fixSetByRedCost_(true),
	printLP_(false),
	maxConAdd_(100),
	maxConBuffered_(100),
	maxVarAdd_(100),
	maxVarBuffered_(100),
	maxIterations_(numeric_limits<int>::max()),
	eliminateFixedSet_(false),
	newRootReOptimize_(false),
	showAverageCutDistance_(false),
	conElimMode_(NoConElim),
	varElimMode_(NoVarElim),
	conElimEps_(0.001),
	varElimEps_(0.001),
	conElimAge_(1),
	varElimAge_(1),
	status_(Unprocessed),
	nSub_(0),
	nLp_(0),
	highestLevel_(0),
	nFixed_(0),
	nAddCons_(0),
	nRemCons_(0),
	nAddVars_(0),
	nRemVars_(0),
	nNewRoot_(0)
{
	_createLpMasters();
	// Master::Master(): allocate some members
	history_ = new History(this);
	openSub_ = new OpenSub(this);
	fixCand_ = new FixCand(this);


	_initializeParameters();

	// start the timers
	totalCowTime_.start();
	totalTime_.start();
}


Master::~Master()
{
	if (treeStream_ != &cout)
		delete treeStream_;
	delete history_;
	delete conPool_;
	delete cutPool_;
	delete varPool_;
	delete openSub_;
	delete fixCand_;
	_deleteLpMasters();
}


Master::STATUS Master::optimize()
{
	// startup the optimization
	/* We call \a _initializeParameters() again, because the call in the
	*   constructor of Master will not call a possibly redefined virtual
	*   function in a derived class.
	*
	*   The default implementation of the virtual function \a initializeOptimization()
	*   does nothing. A reimplementation can be used for problem specific
	*   initializations.
	*/
	_initializeParameters();

	// initialize the tree-log file
	/* If \a VbcLog_ is \a File, then we generate a log-file which can be read
	*   by the VBC-tool to visualize the implicit enumeration. The name
	*   of the log-file is composed of the problem name, the process ID, and
	*   the extension {\tt .tree}. The length of the array storing the name
	*   of the log-file is sufficient for process ids represented by integers with
	*   up to 64 bit (a 64 bit integer has most 19 decimal digits).
	*   If \a VbcLog_ is \a Pipe then we write the instructions to standard out.
	*/
	if (VbcLog_ == File) {
		string treeStreamName = problemName_ + "." + to_string(ogdf::System::getProcessID()) + ".tree";
		treeStream_ = new ofstream(OGDF_STRING_OPEN(treeStreamName));
	}
	else if (VbcLog_ == Pipe)
	{
		string pipeName;
		if (getParameter("VbcPipeName", pipeName)==0)
			treeStream_ = new ofstream(pipeName.c_str());
		else
			treeStream_ = &cout;
	}

	// initialize the primal and the dual bound
	/* The initialization of the primal and dual bound cannot be performed
	*   in any case in the constructor, as the sense of optimization
	*   might be unknown, if it is read from the constructor of a derived class.
	*   Therefore, we initialize these bounds now. If the initialization of the
	*   sense of the optimization has been forgotten, then the program terminates
	*   with an error.
	*/

	// initialize the bounds according to sense of optimization
	switch (optSense_.sense()) {
	case OptSense::Min:
		primalBound_ =  infinity();
		dualBound_   = -infinity();
		break;
	case OptSense::Max:
		primalBound_ = -infinity();
		dualBound_   =  infinity();
		break;
	case OptSense::Unknown:
		Logger::ifout() << "Master::optimize(): optimization sense unknown.\n"
			<< "Specify optimization sense in the constructor or use function initializeOptSense().\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}

	// depending on \a PRIMALBOUNDMODE reinitialize \a primalBound
	double opt;
	if (knownOptimum(opt)) {
		switch (pbMode_) {
		case Optimum: primalBound_= opt;
			break;
		case OptimumOne:
			if (optSense_.max()) primalBound_ = opt - 1.0;
			else                 primalBound_ = opt + 1.0;
		default: break;
		}
	}

	initializeOptimization();
	// print the parameters
	/* The parameters are only not printed if the output level or the log
	*   level is \a Silent.
	*/

	printParameters();

	Logger::ilout(Logger::LL_MINOR) << endl << "   #sub   #open   current   #iter         LP       dual     primal" << endl;

	// perform the branch and bound algorithm
	// initialize the root node
	/* The function \a firstSub() is a pure virtual function returning the problem
	*   specific subproblem derived from the class Sub, which is associated with
	*   the root node of the \bab\ tree.
	*/
	status_ = Processing;
	root_   = firstSub();
	rRoot_  = root_;

	openSub_->insert(root_);

	treeInterfaceNewNode(root_);

	// process all subproblems
	/* The following lines of code are the core of the \bab\ algorithm.
	*   The function \a select() returns subproblems which still have to be
	*   processed until either the set of open subproblem becomes empty or
	*   one of the resource limits (e.g., maximal cpu time) is exceeded.
	*   If the optimization of a subproblem fails we quit the optimization
	*   immediately..
	*/
	Sub *current;

	while ((current = select())) {
		++nSubSelected_;

		if (current->optimize()) {
			status_ = Error;
			break;
		}
	}

	if (status_ == Processing) status_ = Optimal;

	// output history and statistics
	/* The virtual dummy function \a output() can be redefined in derived classes
	*   for problem specific output.
	*/

	// output the solution history
	Logger::ilout(Logger::LL_DEFAULT) << endl << *history_ << endl;

	// output information on the tree, variables, constraints, etc.
	const int w = 6;

	Logger::ilout(Logger::LL_DEFAULT) << "Miscellaneous Statistics" << endl << endl
	 << "  Dual bound of the root node       : "
	 << setw(w) << rootDualBound_     << endl
	 << "  Number of subproblems             : "
	 << setw(w) << nSubSelected_     << endl
	 << "  Number of solved LPs              : "
	 << setw(w) << nLp_      << endl
	 << "  Highest level in tree             : "
	 << setw(w) << highestLevel_   << endl
	 << "  Number of fixed variables         : "
	 << setw(w) << nFixed_   << endl
	 << endl
	 << "  Number of added constraints       : "
	 << setw(w) << nAddCons_ << endl
	 << "  Number of removed constraints     : "
	 << setw(w) << nRemCons_ << endl
	 << "  Number of added variables         : "
	 << setw(w) << nAddVars_ << endl
	 << "  Number of removed variables       : "
	 << setw(w) << nRemVars_ << endl
	 << endl
	 << "  Number of root changes            : "
	 << setw(w) << nNewRoot_ << endl;

	_outputLpStatistics();

	output();

	Logger::ilout(Logger::LL_DEFAULT) << endl << endl;

	// output the timing statistics
	/* The cpu time for branching may include time for linear programming.
	*   Therefore, we should not subract it from the total time. Otherwise,
	*   the LP time for LP ranking of branching rules would be subtracted twice.
	*/
	double lpPercent;   //!< percentage of lp time of total cpu time
	double sepPercent;  //!< percentage of separation time of total cpu time
	double heuPercent;  //!< percentage of heuristic time of total cpu time
	double priPercent;  //!< percentage of pricing time of total cpu time
	double misPercent;  //!< percentage of miscellaneous time of total cpu time

	double totcsecs = (double) totalTime_.centiSeconds();

	// miscellaneous time in centi-seconds
	__int64 misTime = totalTime_.centiSeconds()
		- lpTime_.centiSeconds()
		- separationTime_.centiSeconds()
		- improveTime_.centiSeconds()
		- pricingTime_.centiSeconds();

	ogdf::StopwatchCPU misTimer(10*misTime);  //!< required for nice output of the value

	if (totcsecs < eps()) totcsecs = 1.0;  //!< avoid zero division in following lines

	lpPercent  = (double)lpTime_.centiSeconds() / totcsecs*100.0;
	sepPercent = (double)separationTime_.centiSeconds() / totcsecs*100.0;
	heuPercent = (double)improveTime_.centiSeconds() / totcsecs*100.0;
	priPercent = (double)pricingTime_.centiSeconds() / totcsecs*100.0;
	misPercent = (double)misTime / totcsecs*100.0;

	const int wpc = 7;

	Logger::ilout(Logger::LL_DEFAULT) << "Timing Statistics" << endl << endl
	 << "  Elapsed time           : " << totalCowTime_ << endl
	 << "  Total cpu time         : " << totalTime_   << endl
	 << "  LP cpu time            : " << lpTime_
	 << "  (" << setw(wpc) <<  lpPercent << "%)" << endl
	 << "  LP solver cpu time     : " << lpSolverTime_ << endl
	 << "  Separation cpu time    : " << separationTime_
	 << "  (" << setw(wpc) <<  sepPercent << "%)" << endl
	 << "  Heuristics cpu time    : " << improveTime_
	 << "  (" << setw(wpc) << heuPercent << "%)" << endl
	 << "  Pricing cpu time       : " << pricingTime_
	 << "  (" << setw(wpc) << priPercent << "%)" << endl
	 << "  Branching cpu time     : " << branchingTime_ << endl
	 << "  Miscellaneous cpu time : " << misTimer
	 << "  (" << setw(wpc) << misPercent << "%)" << endl << endl;

	// output the value of the best solution
	if (feasibleFound())
		Logger::ilout(Logger::LL_DEFAULT) << "Best solution: " << primalBound_ << endl;
	else
		Logger::ilout(Logger::LL_DEFAULT) << "No feasible solution found." << endl;

	// output the status of the optimization
	Logger::ilout(Logger::LL_DEFAULT) << endl << "ABACUS optimization terminated with status "
	 << STATUS_[status_] << "." << endl;

	// Master::optimize(): clean up and return
	/* Before cleaning up we call the virtual function \a terminateOptimization().
	*   The default implementation of this function does nothing but it
	*   can be used as entrance point for problem specific actions
	*   (e.g., visualization of the best solution).
	*
	*   By deleting the \a root_ of the branch and cut tree, recursively the
	*   destructors of all subproblems in the \bac\ tree are called.
	*
	*   After the optimization we delete the candidates for fixing such
	*   that the variable pools can be safely deleted, as it is an error
	*   to delete a still referenced slot of a pool.
	*/

	terminateOptimization();

	delete root_;

	fixCand_->deleteAll();

	if (VbcLog_ == Pipe) {
		writeTreeInterface("#END_OF_OUTPUT");
	}

	return status_;
}


void Master::initializePools(
	ArrayBuffer<Constraint*>  &constraints,
	ArrayBuffer<Variable*>    &variables,
	int varPoolSize,
	int cutPoolSize,
	bool dynamicCutPool)
{
	// initialize the constraint pool of the Master
	const int nCons = constraints.size();

	conPool_ = new StandardPool<Constraint, Variable>(this, nCons);

	for (int i = 0; i < nCons; i++)
		conPool_->insert(constraints[i]);

	// initialize the variable pool of the Master
	/* The size of the variable pool must be at least the number of variables
	*   which are inserted in this function. A variable pool is always an
	*   automatically growing pool, as the addition of a variable must be
	*   always possible for the correctness of the algorithm. Therefore, the
	*   third argument of the constructor of StandardPool is \a true.
	*/
	const int nVars = variables.size();

	if (varPoolSize < nVars) varPoolSize = nVars;

	varPool_ = new StandardPool<Variable, Constraint>(this, varPoolSize, true);

	for (int i = 0; i < nVars; i++)
		varPool_->insert(variables[i]);

	// set up the pool for dynamically generated constraints
	/* If the third argument of the constructor of StandardPool is missing
	*   the size of the pool is fixed, of the third argument is \a true then
	*   the pool is automatically increased if it is full and an insertion
	*   is performed.
	*/
	if (cutPoolSize > 0)
		cutPool_ = new StandardPool<Constraint, Variable>(this, cutPoolSize, dynamicCutPool);
}


void Master::initializePools(
	ArrayBuffer<Constraint*> &constraints,
	ArrayBuffer<Constraint*> &cuts,
	ArrayBuffer<Variable*> &variables,
	int varPoolSize,
	int cutPoolSize,
	bool dynamicCutPool)
{
	initializePools(constraints, variables, varPoolSize, cutPoolSize, dynamicCutPool);

	// is there enough space for all the cuts in the cut pool
	/* Even if the cut pool is dynamic it is an error if the cut pool is
	*   too small. Because before the cut pool is enlarged, inactive, i.e.\
	*   just inserted, constraints would be removed.
	*/
	const int nCuts = cuts.size();

	if (nCuts > cutPoolSize) {
		//char *_error=new char[255];
		Logger::ifout() << "Master::initializePools(): size of cutting plane too small for all initialize cuts.\nsize of cut pool: " << cutPoolSize << "\n number of cuts: " << nCuts << "\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}

	for (int i = 0; i < nCuts; i++)
		cutPool_->insert(cuts[i]);
}


Sub* Master::select()
{
	// check if we should terminate the optimization
	/* If one of the criteria for early termination is satisfied then
	*   we fathom all subproblems of the tree, in order to perform a
	*   correct cleaning up.
	*
	*   The maximal level of the enumeration tree is no termination criterion
	*   in this sense, because it only prevents the generation of further
	*   sons of subproblems having this maximal level, but does not
	*   stop the optimization.
	*/
	if (totalTime_.exceeds(maxCpuTime())) {
		Logger::ilout(Logger::LL_DEFAULT) << "Maximal CPU time " << maxCpuTimeAsString() << " exceeded." << endl
		 << "Stop optimization." << endl;
		root_->fathomTheSubTree();
		status_ = MaxCpuTime;
		return 0;
	}

	if (totalCowTime_.exceeds(maxCowTime())) {
		Logger::ilout(Logger::LL_DEFAULT) << "Maximal elapsed time " << maxCowTimeAsString() << " exceeded." << endl
		 << "Stop optimization." << endl;
		root_->fathomTheSubTree();
		status_ = MaxCowTime;
		return 0;
	}

	if (guaranteed()) {
		Logger::ilout(Logger::LL_DEFAULT) << endl
		 << "Guarantee " << requiredGuarantee() << " % reached." << endl
		 << "Terminate optimization." << endl;
		status_ = Guaranteed;
		root_->fathomTheSubTree();
		return 0;
	}

	if (nSubSelected_ >= maxNSub()) {
		Logger::ilout(Logger::LL_DEFAULT) << endl
		 << "Maximal number of subproblems reached: " << maxNSub() << endl
		 << "Terminate optimization." << endl;
		status_ = MaxNSub;
		root_->fathomTheSubTree();
		return 0;
	}

	return openSub_->select();
}


int Master::enumerationStrategy(const Sub *s1, const Sub *s2)
{
	switch (enumerationStrategy_) {
	case BestFirst:    return bestFirstSearch(s1, s2);
	case BreadthFirst: return breadthFirstSearch(s1, s2);
	case DepthFirst:   return depthFirstSearch(s1, s2);
	case DiveAndBest:  return diveAndBestFirstSearch(s1, s2);
	default:
		Logger::ifout() << "Master::enumerationStrategy(): Unknown enumeration strategy\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
}


int Master::bestFirstSearch(const Sub *s1, const Sub *s2) const
{
	double dual1 = s1->dualBound();
	double dual2 = s2->dualBound();

	if (optSense()->max()) {
		if (dual1 > dual2) return  1;
		if (dual1 < dual2) return -1;
		return equalSubCompare(s1, s2);
	}
	else {
		if (dual1 > dual2) return -1;
		if (dual1 < dual2) return  1;
		return equalSubCompare(s1, s2);
	}
}


int Master::equalSubCompare(const Sub *s1, const Sub *s2) const
{
	if (!s1->branchRule() || !s2->branchRule())
		return 0;

	if (s1->branchRule()->branchOnSetVar() && s2->branchRule()->branchOnSetVar())
	{
		if (((SetBranchRule*) s1->branchRule())->setToUpperBound())
		{
			if (((SetBranchRule*) s2->branchRule())->setToUpperBound())
				return 0;
			else
				return 1;
		}
		else if (((SetBranchRule*) s2->branchRule())->setToUpperBound())
			return -1;
		return 0;
	}
	else return 0;
}


int Master::depthFirstSearch(const Sub* s1, const Sub* s2) const
{
	if(s1->level() > s2->level()) return  1;
	if(s1->level() < s2->level()) return -1;
	return equalSubCompare(s1, s2);
}


int Master::breadthFirstSearch(const Sub* s1, const Sub* s2) const
{
	if (s1->level() > s2->level()) return -1;
	if (s1->level() < s2->level()) return  1;

	if (s1->id() < s2->id()) return 1;
	return -1;
}


int Master::diveAndBestFirstSearch(const Sub *s1, const Sub* s2) const
{
	if (feasibleFound()) return bestFirstSearch(s1, s2);
	else                 return depthFirstSearch(s1, s2);
}


void Master::primalBound(double x)
{
	if (optSense()->max()) {
		if (x < primalBound_) {
			Logger::ifout() << "Error: Master::primalBound(): got worse\nold bound: " << primalBound_ << "\nnew bound: " << x << "\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcPrimalBound);
		}
	}
	else
		if (x > primalBound_) {
			Logger::ifout() << "Error: Master::primalBound(): got worse\nold bound: " << primalBound_ << "\nnew bound: " << x << "\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcPrimalBound);
		}

		// make sure that this is an integer value for \a objInteger_
		/* After the test for integrality, we round the value to the next
		*   integer. The reason is that a new primal bound coming obtained
		*   by the solution of the relaxation in a subproblem might have some
		*   garbage in it (e.g., \a 10.000000000034 or \a 9.999999999987).
		*   With this garbage the
		*   fathoming of node might fail in the subproblem optimization, although
		*   it is already correct.
		*/
		if (objInteger_) {
			if (!isInteger(x, eps())) {
				Logger::ifout() << "Master::primalBound(): value " << x << " is not integer, but feasible solutions with integer objective function values are expected.\n";
				OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcNotInteger);
			}

			x = floor(x + eps());
		}

		primalBound_ = x;

		// update the primal bound in the Tree Interface
		if (optSense()->max()) treeInterfaceLowerBound(x);
		else                   treeInterfaceUpperBound(x);

		history_->update();
}


void Master::dualBound(double x)
{
	if (optSense()->max()) {
		if (x > dualBound_) {
			Logger::ifout() << "Error: Master::dualBound(): got worse\nold bound: " << dualBound_ << "\nnew bound: " << x << "\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcDualBound);
		}
	}
	else
		if (x < dualBound_) {
			Logger::ifout() << "Error: Master::dualBound(): got worse\nold bound: " << dualBound_ << "\nnew bound: " << x << "\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcDualBound);
		}

		dualBound_ = x;

		//! update the dual bound in the Tree Interface
		if (optSense()->max()) treeInterfaceUpperBound(x);
		else                   treeInterfaceLowerBound(x);

		history_->update();
}


bool Master::betterDual(double x) const
{
	if (optSense_.max()) {
		return x < dualBound_;
	}
	else {
		return x > dualBound_;
	}
}


bool Master::primalViolated(double x) const
{
	if (optSense_.max()) {
		if (objInteger_) {
			return x <= primalBound_;
		}
		else {
			return x + eps() <= primalBound_;
		}
	}
	else {
		if (objInteger_) {
			return x >= primalBound_;
		}
		else {
			return x - eps() >= primalBound_;
		}
	}
}


bool Master::betterPrimal(double x) const
{
	if (optSense_.max()) {
		return x > primalBound_;
	}
	else {
		return x < primalBound_;
	}
}


bool Master::feasibleFound() const
{
	if (optSense_.max()) {
		return primalBound_ > -infinity();
	}
	else {
		return primalBound_ <  infinity();
	}
}


void Master::rRoot(Sub *newRoot, bool /* reoptimize */)
{
	if (rRoot_ == newRoot) return;

	rRoot_ = newRoot;

	Logger::ilout(Logger::LL_DEFAULT) << "\t" << "subproblem " << newRoot->id() << " is now root of remaining tree" << endl;

	if ((newRoot->status() == Sub::Processed ||
		newRoot->status() == Sub::Dormant     ) && newRootReOptimize_)
		newRoot->reoptimize();

	++nNewRoot_;
}


bool Master::guaranteed() const
{
	if (fabs(lowerBound()) < machineEps() && fabs(upperBound()) > machineEps())
		return false;

	if (guarantee() + machineEps() < requiredGuarantee_)
		return true;
	else
		return false;
}


double Master::guarantee() const
{
	double lb = lowerBound();

	if (fabs(lb) < machineEps()) {
		if (fabs(upperBound()) < machineEps())
			return 0.0;
		else {
			Logger::ifout() << "Master::guarantee(): cannot compute guarantee with lower bound 0\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
		}
	}
	return fabs((upperBound() - lb)/lb * 100.0);
}


void Master::printGuarantee() const // warning: this should only be called if it has already been ensured that logging allows this output!
{
	double lb = lowerBound();
	double ub = upperBound();

	if (lb == -infinity() || ub == infinity() ||
		((fabs(lb) < machineEps()) && (fabs(ub) > machineEps())))
		Logger::ifout() << "---";
	else
		Logger::ifout() << guarantee() << '%';
}


bool Master::check() const
{
	double optVal;

	if (!knownOptimum(optVal)) return false;

	if (optVal - eps() < primalBound() && primalBound() < optVal + eps())
		return true;
	else
		return false;
}


bool Master::knownOptimum(double &optVal) const
{
	// open the file storing the optimum solutions
	ifstream optimumFile(optimumFileName_.c_str(), ios::in );

	if (!optimumFile) return false;

	// go through the lines of the file until the problem is found
	/* As the problem name might be preceded by a path we check if the
	*   \a name read from the file, matches the last letters of the problem name.
	*/
	char name[256];
	double value;

	while (!optimumFile.eof()) {
		optimumFile >> name >> value;

		if (endsWith(problemName_, name)) {
			optVal = value;
			return true;
		}
		optimumFile >> ws;
	}

	return false;
}


void Master::writeTreeInterface(const string &info, bool time) const
{
	if (VbcLog_ == NoVbc) return;

	if (VbcLog_ == Pipe) *treeStream_ << '$';
	if (VbcLog_ == File && time) *treeStream_ << totalTime_ << ' ';
	*treeStream_ << info << endl;
}


void Master::treeInterfaceNewNode(Sub *sub) const
{
	if (VbcLog_ != NoVbc) {
		int fatherId =  (sub == root_) ? 0 : sub->father()->id();

		string info = string("N ") + to_string(fatherId) + " " + to_string(sub->id()) + " 5";
		writeTreeInterface(info);
	}
}


void Master::treeInterfacePaintNode(int id, int color) const
{
	if (VbcLog_ == NoVbc) return;

	string info = string("P ") + to_string(id) + " " + to_string(color);
	writeTreeInterface(info);
}


void Master::treeInterfaceLowerBound(double lb) const
{
	if (VbcLog_ == NoVbc) return;

	string info = string("L ") + to_string(lb);
	writeTreeInterface(info);
}


void Master::treeInterfaceUpperBound(double ub) const
{
	if (VbcLog_ == NoVbc) return;

	string info = string("U ") + to_string(ub);
	writeTreeInterface(info);
}


void Master::treeInterfaceNodeBounds(int id, double lb, double ub)
{
	if (VbcLog_ == NoVbc) return;

	//char string[256];
	ostringstream info;

	if (isInfinity(fabs(lb))) {
		if (isInfinity(fabs(ub)))
			info << "I " << id << " \\iLower Bound: ---\\nUpper Bound:  ---\\i";
		else
			info << "I " << id << " \\iLower Bound: ---\\nUpper Bound:  "
				<< ios::fixed << std::setprecision(2) << std::setw(6) << ub << "\\i";
	}
	else {
		if (isInfinity(fabs(ub)))
			info << "I " << id << " \\iLower Bound: "
				<< ios::fixed << std::setprecision(2) << std::setw(6) << lb
				<< "\\nUpper Bound:  ---\\i";
		else
			info << "I " << id << " \\iLower Bound: "
				<< ios::fixed << std::setprecision(2) << std::setw(6) << lb
				<< "\\nUpper Bound:  "
				<< ios::fixed << std::setprecision(2) << std::setw(6) << ub
				<< "\\i";
	}

	writeTreeInterface(info.str());
}


void Master::newSub(int level)
{
	++nSub_;
	if(level > highestLevel_)
		highestLevel_ = level;
}


void Master::_initializeParameters()
{
	// Master::_initializeParameters(): local variables
	/* The function \a _initializeLpParameters() initializes the LP solver
	*   specific Parameters. It is called after the parameter-file
	*   {\tt .abacus} was read and all general parameters were processed.
	*
	*   This function is implemented in the file \a lpif.cc.
	*/

	if(readParamFromFile_) {
		// set up the name of the configuration file
		// The location of the configuration file .abacus is given by the environment variable ABACUS_DIR.
		char *abacusDir = getenv("ABACUS_DIR");
		if (!abacusDir) {
			Logger::ifout() << "environment variable ABACUS_DIR not found\n";
			OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
		}

#ifdef OGDF_SYSTEM_UNIX
		string configFileName = string(abacusDir) + "/.abacus";
#else
		string configFileName = string(abacusDir) + "\\.abacus";
#endif

		readParameters(configFileName);

	} else {

		// set default values for abacus parameters
		insertParameter("EnumerationStrategy","BestFirst");
		insertParameter("BranchingStrategy","CloseHalfExpensive");
		insertParameter("NBranchingVariableCandidates","1");
		insertParameter("NStrongBranchingIterations","50");
		insertParameter("Guarantee","0.0");
		insertParameter("MaxLevel","999999");
		insertParameter("MaxNSub","9999999");
		insertParameter("MaxCpuTime","99999:59:59");
		insertParameter("MaxCowTime","99999:59:59");
		insertParameter("ObjInteger","false");
		insertParameter("TailOffNLps","0");
		insertParameter("TailOffPercent","0.0001");
		insertParameter("DelayedBranchingThreshold","0");
		insertParameter("MinDormantRounds","1");
		insertParameter("PrimalBoundInitMode","None");
		insertParameter("PricingFrequency","0");
		insertParameter("SkipFactor","1");
		insertParameter("SkippingMode","SkipByNode");
		insertParameter("FixSetByRedCost","true");
		insertParameter("PrintLP","false");
		insertParameter("MaxConAdd","100");
		insertParameter("MaxConBuffered","100");
		insertParameter("MaxVarAdd","500");
		insertParameter("MaxVarBuffered","500");
		insertParameter("MaxIterations","-1");
		insertParameter("EliminateFixedSet","false");
		insertParameter("NewRootReOptimize","false");
		insertParameter("ShowAverageCutDistance","false");
		insertParameter("ConstraintEliminationMode","Basic");
		insertParameter("ConElimEps","0.001");
		insertParameter("ConElimAge","1");
		insertParameter("VariableEliminationMode","ReducedCost");
		insertParameter("VarElimEps","0.001");
		insertParameter("VarElimAge","1");
		insertParameter("VbcLog","None");
#if defined(COIN_OSI_CPX)
	insertParameter("DefaultLpSolver","CPLEX");
#elif defined(COIN_OSI_SYM)
	insertParameter("DefaultLpSolver","SYMPHONY");
#elif defined(COIN_OSI_GRB)
	insertParameter("DefaultLpSolver","Gurobi");
#else
		insertParameter("DefaultLpSolver","Clp");
#endif
		insertParameter("SolveApprox","false");

		_setDefaultLpParameters();
	}

	// let the user set/overwrite parameters in 4 steps.
	// 1) Assign parameters that were read from the global abacus config file
	assignParameters();

	// 2) Allow user to read parameters from custom config file
	initializeParameters();

	// 3) Assign parameters from custom config file.
	// this assigns all parameters that were read from the abacus config file
	// in the abacus directory and those that the user read from a file
	// in initializeParameters(). However, this will overwrite any parameter
	// that was set directly.
	assignParameters();

	// 4) Allow user to set parameters directly
	// call initializeParameters again so that parameters that were set directly
	// (and not via a config file) can be set (undo the overwrite from step 2).
	// Since we do *not* call assignParameters again, those parameters are not
	// overwritten this time.
	initializeParameters();

	_initializeLpParameters();
}


void Master::assignParameters()
{
	string stringVal;   //!< auxiliary variable for reading strings

	// get the enumeration strategy
	enumerationStrategy_= (ENUMSTRAT)
		findParameter("EnumerationStrategy",4,ENUMSTRAT_);

	// get the branching strategy

	branchingStrategy_=(BRANCHINGSTRAT)
		findParameter("BranchingStrategy",2,BRANCHINGSTRAT_);

	// get the number of tested candidates for branching variables
	assignParameter(nBranchingVariableCandidates_, "NBranchingVariableCandidates",0, numeric_limits<int>::max());

	// get the number of simplex iterations to perform when testing
	// branching candidates within strong branching
	assignParameter(nStrongBranchingIterations_, "NStrongBranchingIterations",-1, numeric_limits<int>::max());

	// get the solution guarantee
	assignParameter(requiredGuarantee_,"Guarantee",0.0,infinity());

	// get the maximal level in the enumeration tree
	assignParameter(maxLevel_, "MaxLevel", 1, numeric_limits<int>::max());

	// get the maximal level in the enumeration tree
	assignParameter(maxNSub_, "MaxNSub", 1, numeric_limits<int>::max());

	// get the maximal cpu time
	assignParameter(stringVal,"MaxCpuTime",0);
	maxCpuTime(stringVal);

	// get the maximal wall clock time
	assignParameter(stringVal,"MaxCowTime",0);
	maxCowTime(stringVal);

	// get the integrality condition on the objective function
	assignParameter(objInteger_,"ObjInteger" );

	// get the number of linear programs for the tailing off analysis
	assignParameter(tailOffNLp_, "TailOffNLps", numeric_limits<int>::min(), numeric_limits<int>::max());

	// get the minimal improvement for the tailing off analysis
	assignParameter(tailOffPercent_, "TailOffPercent", 0.0, infinity());

	// get the threshold for delayed branching
	assignParameter(dbThreshold_, "DelayedBranchingThreshold", 0, numeric_limits<int>::max());

	// get the minimal number of rounds a subproblem should stay dormant
	assignParameter(minDormantRounds_,  "MinDormantRounds", 1, numeric_limits<int>::max());

	//initializeParameters();

	// get the primal bound initialization mode
	pbMode_ = (PRIMALBOUNDMODE) findParameter("PrimalBoundInitMode", 3, PRIMALBOUNDMODE_);

	// get the pricing frequency
	assignParameter(pricingFreq_, "PricingFrequency", 0, numeric_limits<int>::max());

	// get the skip factor
	assignParameter(skipFactor_, "SkipFactor", 0, numeric_limits<int>::max());

	// get the skipping mode

	skippingMode_ = (SKIPPINGMODE) findParameter("SkippingMode", 2, SKIPPINGMODE_);

	// is fixing/setting by reduced costs turned on?
	assignParameter(fixSetByRedCost_, "FixSetByRedCost");

	// should the LP be output every iteration?
	assignParameter(printLP_, "PrintLP");

	// get the maximal number of added constraints
	assignParameter(maxConAdd_, "MaxConAdd", 0, numeric_limits<int>::max());

	// get the size of the buffer for constraints
	assignParameter(maxConBuffered_, "MaxConBuffered", 0, numeric_limits<int>::max());

	// get the maximal number of added variables
	assignParameter(maxVarAdd_, "MaxVarAdd", 0, numeric_limits<int>::max());

	// get the size of the buffer for variables
	assignParameter(maxVarBuffered_, "MaxVarBuffered", 0, numeric_limits<int>::max());

	// get the maximal number of iterations in the cutting plane phase
	assignParameter(maxIterations_, "MaxIterations", -1, numeric_limits<int>::max());

	// should fixed or set variables be eliminated from the LP?
	assignParameter(eliminateFixedSet_, "EliminateFixedSet");

	// should a new root be reoptimized?
	assignParameter(newRootReOptimize_, "NewRootReOptimize");

	// get the name of the file storing the optimal solutions
	/* The name of the file storing the optimal solutions is only an optional
	*   parameter which can be missing.
	*/
	getParameter("OptimumFileName", optimumFileName_);

	// should the average distance of the cuts per iteration be output?
	assignParameter(showAverageCutDistance_, "ShowAverageCutDistance");

	// get the constraint elimination mode
	conElimMode_ = (CONELIMMODE) findParameter("ConstraintEliminationMode", 3, CONELIMMODE_);

	// get the variable elimination mode
	varElimMode_ = (VARELIMMODE) findParameter("VariableEliminationMode", 2, VARELIMMODE_);

	// get the tolerance for variable elimination
	assignParameter(varElimEps_, "VarElimEps", 0.0, infinity());

	// get the tolerance for constraint elimination
	assignParameter(conElimEps_, "ConElimEps", 0.0, infinity());

	// get the age for constraint elimination
	assignParameter(conElimAge_, "ConElimAge", 1, numeric_limits<int>::max());

	// get the age for variable elimination
	assignParameter(varElimAge_, "VarElimAge", 1, numeric_limits<int>::max());

	// should a log-file of the enumeration tree be generated?
	VbcLog_ = (VBCMODE) findParameter("VbcLog", 3, VBCMODE_);

	//! get the default LP-solver
	defaultLpSolver_ = (OSISOLVER) findParameter("DefaultLpSolver", 12, OSISOLVER_);
	assignParameter(solveApprox_, "SolveApprox", false);
}


void Master::printParameters() const
{
	Logger::ilout(Logger::LL_DEFAULT) << "Branch and Cut Parameters:" << endl << endl

	 << "  Enumeration strategy                   : "
	 << ENUMSTRAT_[enumerationStrategy_]
	 << endl

	 << "  Branching Strategy                     : "
	 << BRANCHINGSTRAT_[branchingStrategy_]
	 << endl
	 << "  Tested candidates for branching var.   : "
	 << nBranchingVariableCandidates_ << endl
	 << "  Simplex iterations when testing" << endl
	 << "         candidates for branching var.   : "
	 << nStrongBranchingIterations_ << endl
	 << "  Guarantee                              : "
	 << requiredGuarantee_ << " %" << endl
	 << "  Maximal enumeration level              : "
	 << maxLevel_ << endl
	 << "  Maximal number of subproblems          : "
	 << maxNSub_ << endl
	 << "  CPU time limit                         : "
	 << maxCpuTimeAsString() << endl
	 << "  Wall-clock time limit                  : "
	 << maxCowTimeAsString() << endl
	 << "  Objective function values integer      : "
	 << onOff(objInteger_) << endl
	 << "  Tailing Off Parameters" << endl
	 << "                    Number of LPs        : "
	 << tailOffNLp_ << endl
	 << "                    Minimal improvement  :    "
	 << tailOffPercent_ << '%' << endl
	 << "  Delayed branching threshold            : "
	 << dbThreshold_ << endl
	 << "  Maximal number of dormant rounds       : "
	 << minDormantRounds_ << endl

	 << "  Primal Bound Initialization            : "
	 << PRIMALBOUNDMODE_[pbMode_]
	 << endl
	 << "  Frequency of additional pricing        : "
	 << pricingFreq_ << " LPs" << endl
	 << "  Cutting skip factor                    : "
	 << skipFactor_ << endl
	 << "  Skipping mode                          : "
	 << ((skippingMode_ == SkipByNode)? "by node": "by tree") << endl

	<< "  Fix/set by reduced costs               : "
	 << onOff(fixSetByRedCost_) << endl
	 << "  Output of the linear program           : "
	 << onOff(printLP_) << endl
	 << "  Maximal number of added constraints    : "
	 << maxConAdd_ << endl
	 << "  Maximal number of buffered constraints : "
	 << maxConBuffered_ << endl
	 << "  Maximal number of added variables      : "
	 << maxVarAdd_ << endl
	 << "  Maximal number of buffered variables   : "
	 << maxVarBuffered_ << endl
	 << "  Maximal number of iterations per" << endl
	 << "                     cutting plane phase : "
	 << maxIterations_ << endl
	 << "  Elimination of fixed and set variables : "
	 << onOff(eliminateFixedSet_) << endl
	 << "  Reoptimization after a root change     : "
	 << onOff(newRootReOptimize_) << endl
	 << "  File storing optimum solutions         : "
	 << optimumFileName_ << endl
	 << "  Show average distance of added cuts    : "
	 << onOff(showAverageCutDistance_) << endl
	 << "  Elimination of constraints             : "
	 << CONELIMMODE_[conElimMode_] << endl
	 << "  Elimination of variables               : "
	 << VARELIMMODE_[varElimMode_] << endl
	 << "  Tolerance for constraint elimination   : "
	 << conElimEps_ << endl
	 << "  Tolerance for variable elimination     : "
	 << varElimEps_ << endl
	 << "  Age for constraint elimination         : "
	 << conElimAge_ << endl
	 << "  Age for variable elimination           : "
	 << varElimAge_ << endl
	 << "  Default LP-solver                      : "
	 << OSISOLVER_[defaultLpSolver_] << endl
	 << "  Usage of approximate solver            : "
	 << onOff(solveApprox_) << endl;
	_printLpParameters();
}


void Master::maxCpuTime(int hour, int min, int sec)
{

	if(sec >59||min >59){
		Logger::ifout() << "Master::setCpuTime() invalid argument \n - correct value: sec,min <=60\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	maxCpuTime_ = __int64(hour)*3600 + min*60 + sec;
}


void Master::nBranchingVariableCandidates(int n)
{
	if (n < 1) {
		Logger::ifout() << "Master::nBranchingVariableCandidates() invalid argument\ncorrect value: positive integer number\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	nBranchingVariableCandidates_ = n;
}


void Master::nStrongBranchingIterations(int n)
{
	if (n < 1) {
		Logger::ifout() << "Master::nStrongBranchingIterations() invalid argument\ncorrect value: positive integer number\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	nStrongBranchingIterations_ = n;
}


void Master::requiredGuarantee(double g)
{
	if (g < 0.0) {
		Logger::ifout() << "Master::guarantee: " << g << "\nchoose nonnegative value.";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}

	requiredGuarantee_ = g;
}


void Master::maxLevel(int max)
{
	if (max < 1) {
		Logger::ifout() << "Master::maxLevel " << max << ", only positive integers are valid\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	maxLevel_ = max;
}


void Master::maxNSub(int max)
{
	if (max < 1) {
		Logger::ifout() << "Master::maxNSubl " << max << ", only positive integers are valid\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	maxNSub_ = max;
}


void Master::tailOffPercent(double p)
{
	if (p < 0.0) {
		Logger::ifout() << "Master::tailing_off(p): choose nonnegative value\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	tailOffPercent_ = p;
}


bool Master::delayedBranching(int nOpt) const
{
	if (nOpt >= dbThreshold_ + 1) return false;
	else                          return true;
}


void Master::pricingFreq(int f)
{
	if (f < 0) {
		Logger::ifout() << "Master::pricingFreq(): nonnegative frequency expected\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	pricingFreq_ = f;
}


void Master::skipFactor(int f)
{
	if (f < 0) {
		Logger::ifout() << "Master::skipFactor(): nonnegative value expected\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcIllegalParameter);
	}
	skipFactor_ = f;
}


void Master::rootDualBound(double x)
{
	rootDualBound_ = x;
}


bool Master::setSolverParameters(OsiSolverInterface* /* interface */, bool /* solverIsApprox */)
{
	return false;
}


static __int64 getSecondsFromString(const string &str)
{
	// convert time string in seconds
	__int64 s  = 0;   //!< seconds in \a maxTime
	__int64 m  = 0;   //!< minutes in \a maxTime
	__int64 h  = 0;   //!< hours in \a maxTime

	// begin at the end of \a maxTime and read the seconds
	/* We start at the end of \a buf and go left until a \a ':' is found
	*   or we have reached the first character of the string. The string
	*   on the right side of the current position are the seconds.
	*/
	int pos = (int)str.size() - 1;
	while (pos >= 0 && str[pos] != ':')
		pos--;
	s = stoi(str.substr(pos+1));

	int posStop = pos;

	// proceed with the minutes
	/* If we have not reached the beginning of \a buf we replace the \a ':'
	*   by \a '\0' to mark the new end of the string. We go left
	*   until we reach the beginning of the string or the next \a ':'.
	*   The string on the right side of the current position are the
	*   minutes.
	*/
	if(pos >= 0) {
		pos--;
		while (pos >= 0 && str[pos] != ':')
			pos--;
		m = stoi(str.substr(pos+1, posStop-(pos+1)));

		posStop = pos;
	}

	// proceed with the hours
	/* If we still have not reached the beginning of the string, we again replace
	*   the \a ':' by a \a '\0' to mark the new end of the string. If there
	*   are still characters on the left side from the current position,
	*   then this string is the hours.
	*/
	if (pos >= 0) {
		if (--pos >= 0)
			h = stoi(str.substr(0, posStop));
	}

	return s + 60*m + 3600*h;
}


void Master::maxCpuTime(const string &t)
{
	maxCpuTime_ = getSecondsFromString(t);
}


void Master::maxCowTime(const string &t)
{
	maxCowTime_ = getSecondsFromString(t);
}


string Master::maxCowTimeAsString() const
{
	string str;

	__int64 rest = maxCowTime_;
	__int64 sec  = rest % 60;
	rest /= 60;
	__int64 min  = rest % 60;
	rest /= 60;

	str = to_string(rest) + ":";
	if(min < 10) str += '0';
	str += to_string(min);
	str += ':';
	if(sec < 10) str += '0';
	str += to_string(sec);

	return str;
}


string Master::maxCpuTimeAsString() const
{
	string str;

	__int64 rest = maxCpuTime_;
	__int64 sec  = rest % 60;
	rest /= 60;
	__int64 min  = rest % 60;
	rest /= 60;

	str = to_string(rest) + ":";
	if(min < 10) str += '0';
	str += to_string(min);
	str += ':';
	if(sec < 10) str += '0';
	str += to_string(sec);

	return str;
}
} //namespace abacus
