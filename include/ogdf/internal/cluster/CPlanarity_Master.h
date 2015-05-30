/** \file
 * \brief Declaration of the CPlanarityMaster class for the Branch&Cut algorithm
 * for c-planarity testing via an extension to complete connectivity.
 *
 * This class is managing the optimization.
 * Variables and initial constraints are generated and pools are initialized.
 *
 * \author Karsten Klein, Markus Chimani
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifndef OGDF_CPLANARITY_MASTER_H
#define OGDF_CPLANARITY_MASTER_H

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/ArrayBuffer.h>

#include <ogdf/internal/cluster/basics.h>
#include <ogdf/internal/cluster/CPlanar_Edge.h>
#include <ogdf/internal/cluster/CP_MasterBase.h>
#include <ogdf/cluster/ClusterAnalysis.h>


namespace ogdf {



class CPlanarityMaster : public CP_MasterBase {

	friend class CPlanaritySub;

	// Pointers to the given Clustergraph and underlying Graph are stored.
	//const ClusterGraph *m_C;
	//const Graph *m_G;


	// Each time the primal bound is improved, the integer solution induced Graph is built.
	// \a m_solutionGraph is a pointer to the currently best solution induced Graph.
	//GraphCopy *m_solutionGraph;

	//List<nodePair> m_connectionOneEdges;  //<! Contains connection nodePairs whose variable is set to 1.0

public:

	// Construction and default values
	CPlanarityMaster(const ClusterGraph &C,
			//Check what we really still need here
			int heuristicLevel=0,
			int heuristicRuns=2,
			double heuristicOEdgeBound=0.3,
			int heuristicNPermLists=5,
			int kuratowskiIterations=3,
			int subdivisions=10,
			int kSupportGraphs=3,
			double kuratowskiHigh=0.75,
			double kuratowskiLow=0.3,
			bool perturbation=false,
			double branchingGap=0.4,
			const char *time="00:20:00", //maximum computation time
			bool dopricing = false,
			int numAddVariables = 15,
			double strongConstraintViolation = 0.3,
			double strongVariableViolation = 0.3);

	// Destruction
	virtual ~CPlanarityMaster();

	// Initialization of the first Subproblem
	virtual abacus::Sub *firstSub();

	// Returns the number of variables
	int nMaxVars() const {return m_nMaxVars;}

	// Returns a pointer to the underlying Graph
	const Graph *getGraph() const {return m_G;}

	// Returns a pointer to the given Clustergraph.
	const ClusterGraph *getClusterGraph() const {return m_C;}

	// Returns a pointer on the search space graph, which is
	// a copy of the input graph with all edges added that
	// correspond to variables add at initialization.
	// Be aware that this is not dynamically updated, e.g. for pricing.

	const GraphCopy *searchSpaceGraph() const {return m_ssg; }
	// Updates the "best" Subgraph \a m_solutionGraph found so far and fills edge lists with
	// corresponding edges (nodePairs).
	void updateBestSubGraph(List<nodePair> &connection);

	// Returns the optimal solution induced Clustergraph
	Graph *solutionInducedGraph() {return static_cast<Graph*>(m_solutionGraph);}

	// Returns nodePairs of connecting optimal solution edges in list \a edges.
	void getConnectionOptimalSolutionEdges(List<nodePair> &egdes) const;

	// Get parameters
	int getKIterations() const {return m_nKuratowskiIterations;}
	int getNSubdivisions() const {return m_nSubdivisions;}
	int getNKuratowskiSupportGraphs() const {return m_nKuratowskiSupportGraphs;}
	int getHeuristicLevel() const {return m_heuristicLevel;}
	int getHeuristicRuns() const {return m_nHeuristicRuns;}
	double getKBoundHigh() const {return m_kuratowskiBoundHigh;}
	double getKBoundLow() const {return m_kuratowskiBoundLow;}
	bool perturbation() const {return m_usePerturbation;}
	double branchingOEdgeSelectGap() const {return m_branchingGap;}
	double getHeuristicFractionalBound() const {return m_heuristicFractionalBound;}
	int numberOfHeuristicPermutationLists() const {return m_nHeuristicPermutationLists;}
	bool getMPHeuristic() const {return m_mpHeuristic;}
	int getNumAddVariables() const {return m_numAddVariables;}
	double getStrongConstraintViolation() const {return m_strongConstraintViolation;}
	double getStrongVariableViolation() const {return m_strongVariableViolation;}

	// Read global constraint counter, i.e. the number of added constraints of specific type.
	//int addedKConstraints() const {return m_nKConsAdded;}
	//int addedCConstraints() const {return m_nCConsAdded;}


	// Set parameters
	void setKIterations(int n) {m_nKuratowskiIterations = n;}
	void setNSubdivisions(int n) {m_nSubdivisions = n;}
	void setNKuratowskiSupportGraphs(int n) {m_nKuratowskiSupportGraphs = n;}
	void setNHeuristicRuns(int n) {m_nHeuristicRuns = n;}
	void setKBoundHigh(double n) {m_kuratowskiBoundHigh = ((n>0.0 && n<1.0) ? n : 0.8);}
	void setKBoundLow(double n) {m_kuratowskiBoundLow = ((n>0.0 && n<1.0) ? n : 0.2);}
	void heuristicLevel(int level) {m_heuristicLevel = level;}
	void setHeuristicRuns(int n) {m_nHeuristicRuns = n;}
	void setPertubation(bool b) {m_usePerturbation = b;}
	void setHeuristicFractionalBound(double b) {m_heuristicFractionalBound = b;}
	void setHeuristicPermutationLists(int n) {m_nHeuristicPermutationLists = n;}
	void setMPHeuristic(bool b) {m_mpHeuristic = b;}//!< Switches use of lower bound heuristic
	void setNumAddVariables(int i) {m_numAddVariables=i;}
	void setStrongConstraintViolation(double d) { m_strongConstraintViolation=d;}
	void setStrongVariableViolation(double d) { m_strongVariableViolation=d;}
	void setSearchSpaceShrinking(bool b) {m_shrink = b;}//!< Toggles reduction of search space (using only some bag/satchel connections) on/off

	// Updating global constraint counter
	//void updateAddedCCons(int n) {m_nCConsAdded += n;}
	//void updateAddedKCons(int n) {m_nKConsAdded += n;}

	// Returns global primal and dual bounds.
	//double getPrimalBound() {return globalPrimalBound;}
	//double getDualBound() {return globalDualBound;}

	// Cut pools for connectivity and planarity
	//! Returns cut pool for connectivity
	//StandardPool<Constraint, Variable> *getCutConnPool() {return m_cutConnPool;}
	//! Returns cut pool for planarity
	//StandardPool<Constraint, Variable> *getCutKuraPool() {return m_cutKuraPool;}

	//! Returns true if default cut pool is used. Otherwise, separate
	//! connectivity and Kuratowski pools are generated and used.
	//bool &useDefaultCutPool() { return m_useDefaultCutPool;}

#ifdef OGDF_DEBUG
	bool m_solByHeuristic; //solution computed by heuristic or ILP
		// Simple output function to print the given graph to the console.
	// Used for debugging only.
	void printGraph(const Graph &G);
#endif

	//! The name of the file that contains the standard, i.e., non-cut,
	//! constraints (may be deleted by ABACUS and shouldn't be stored twice)
	const char* getStdConstraintsFileName()
	{
		return "StdConstraints.txt";
	}

	int getNumInactiveVars() { return m_inactiveVariables.size();}

	//! Returns reference to cluster nodes member list for \a c.
	const List<node> &getClusterNodes(cluster c) const {
		return m_cNodes[c];
	}
	//! Copies cluster nodes from member list to parameter list.
	//! Should be used if node list needs to be altered.
	void getClusterNodes(cluster c, List<node> &nodeList) const {
		ListConstIterator<node> it = m_cNodes[c].begin();
		while (it.valid())
		{
			nodeList.pushBack(*it);
			++it;
		}
	}

protected:

	// Initializes constraints and variables and an initial dual bound.
	virtual void initializeOptimization();

	// Function that is invoked at the end of the optimization
	virtual void terminateOptimization();

	double heuristicInitialLowerBound();

	//! All variables that have to be present at start of optimization
	//! are created here. Their number is returned.
	void createInitialVariables(List<CPlanarEdgeVar*>& initVars);

	//! Adds inner cluster connection variables in bag-reduced search space.
	void addInnerConnections(cluster c, List<CPlanarEdgeVar*>& connectVars);

	//! Create variables for external cluster connections in case we search
	//! only in the bag-reduced search space.
	void addExternalConnections(cluster c, List<CPlanarEdgeVar*>& connectVars);


	bool isCP() {return feasibleFound();}//(dualBound()!=-infinity());}

	//! Node pair is potential candidate for new edge variable
	bool goodVar(node a, node b);

private:

	// Computes a dual bound for the optimal solution.
	// Tries to find as many edge-disjoint Kuratowski subdivisions as possible.
	// If k edge-disjoint groups of subdivisions are found, the upper bound can be
	// initialized with number of edges in underlying graph minus k.
	double heuristicInitialUpperBound();

	//! Is invoked by heuristicInitialLowerBound()
	double clusterConnection(cluster c, GraphCopy &GC);

	//! Creates variables for complete connectivity
	void createCompConnVars(List<CPlanarEdgeVar*>& initVars);

	//! Computes the graphtheoretical distances of edges incident to node \a u.
	void nodeDistances(node u, NodeArray<NodeArray<int> > &dist);


	// Parameters
	/*int m_nKuratowskiSupportGraphs; 	// Maximal number of times the Kuratowski support graph is computed
	int m_nKuratowskiIterations; 		// Maximal number of times BoyerMyrvold is invoked
	int m_nSubdivisions; 				// Maximal number of extracted Kuratowski subdivisions
	int m_nMaxVars; 					// Max Number of variables
	int m_heuristicLevel; 				// Indicates if primal heuristic shall be used or not
	int m_nHeuristicRuns; 				// Counts how often the primal heuristic has been called

	bool m_usePerturbation; 			// Indicates whether C-variables should be perturbated or not
	double m_branchingGap;				// Modifies the branching behaviour
	double m_heuristicFractionalBound;
	int m_nHeuristicPermutationLists;	// The number of permutation lists used in the primal heuristic
	bool m_mpHeuristic;                 //!< Indicates if simple max planar subgraph heuristic
										// should be used to derive lower bound if only root cluster exists

	double m_kuratowskiBoundHigh;		// Upper bound for deterministic edge addition in computation of the Supportgraph
	double m_kuratowskiBoundLow;		// Lower bound for deterministic edge deletion in computation of the Supportgraph

	int m_numAddVariables;				// how many variables should i add maximally per pricing round?
	double m_strongConstraintViolation; // when do i consider a constraint strongly violated -> separate in first stage
	double m_strongVariableViolation;   // when do i consider a variable strongly violated (red.cost) -> separate in first stage
	*/
	//AbaString *m_maxCpuTime;			// Time threshold for optimization

	// If perturbation is used, this variable stores the largest occuring coeff,
	// i.e. the one closest to 0. Otherwise it corresponds to \a m_epsilon
	//double m_largestConnectionCoeff;

	// Counters for the number of added constraints
	/*int m_nCConsAdded;
	int m_nKConsAdded;
	int m_solvesLP;
	int m_varsInit;
	int m_varsAdded;
	int m_varsPotential;
	int m_varsMax;
	int m_varsCut;
	int m_varsKura;
	int m_varsPrice;
	int m_varsBranch;
	int m_activeRepairs;
	ArrayBuffer<int> m_repairStat;
	inline void clearActiveRepairs() {
		if(m_activeRepairs) {
			m_repairStat.push(m_activeRepairs);
			m_activeRepairs = 0;
		}
	}

	double globalPrimalBound;
	double globalDualBound;

	inline double getDoubleTime(const Stopwatch* act) {
		long tempo = act->centiSeconds()+100*act->seconds()+6000*act->minutes()+360000*act->hours();
		return  ((double) tempo)/ 100.0;
	}
	*/
	//number of calls of the fast max planar subgraph heuristic
	//const int m_fastHeuristicRuns;

	//! Cut pools for connectivity and Kuratowski constraints
	//StandardPool< Constraint, Variable > *m_cutConnPool; //!< Connectivity Cuts
	//StandardPool< Constraint, Variable > *m_cutKuraPool; //!< Kuratowski Cuts

	//! Defines if the ABACUS default cut pool or the separate Connectivity
	//! and Kuratowski constraint pools are used
	//bool m_useDefaultCutPool;

	//!
	//double m_delta;
	//double m_deltaCount;
	//Switch to minimization of additional edges, no delta necessary
	virtual double nextConnectCoeff() {return 1.0;}
	//double nextConnectCoeff() { return  -1  + m_deltaCount--*m_delta; };
	//! Variable creation for nodePair
	virtual CPlanarEdgeVar* createVariable(ListIterator<nodePair>& it) {
		++m_varsAdded;
		CPlanarEdgeVar* v = new CPlanarEdgeVar(this, nextConnectCoeff(), (*it).v1, (*it).v2);
		v->printMe(Logger::slout());
		m_inactiveVariables.del(it);
		//we don't need to check symmetry
		m_varCreated[(*it).v1][(*it).v2] = true;
		return v;
	}
	//! Variable creation for pair of nodes with lower bound
	virtual CPlanarEdgeVar* createVariable(node a, node b, double lbound) {
				OGDF_ASSERT(!(m_varCreated[a][b] || m_varCreated[b][a]));
				++m_varsAdded;
				CPlanarEdgeVar* v = new CPlanarEdgeVar(this, nextConnectCoeff(), lbound, a, b);
				v->printMe(Logger::slout());
				//we don't need to check symmetry
				m_varCreated[a][b] = true;
				return v;
	}
	//! Variable creation for pair of nodes which is not stored in m_inactiveVariables.
	virtual CPlanarEdgeVar* createVariable(node a, node b) {
			OGDF_ASSERT(!(m_varCreated[a][b] || m_varCreated[b][a]));
			++m_varsAdded;
			CPlanarEdgeVar* v = new CPlanarEdgeVar(this, nextConnectCoeff(), a, b);
			v->printMe(Logger::slout());
			//we don't need to check symmetry
			m_varCreated[a][b] = true;
			return v;
	}
	//List<nodePair> m_inactiveVariables;
	//used in initialization
	virtual void generateVariablesForFeasibility(const List<ChunkConnection*>& ccons, List<CPlanarEdgeVar*>& connectVars);
	// Keeps track of created variables
	// NodeArray< NodeArray<bool> > m_varCreated;

	//! writes coefficients of all orig and connect variables in constraint con into
	//! emptied list coeffs
	virtual void getCoefficients(abacus::Constraint* con, const List<CPlanarEdgeVar* > & connect,
			List<double> & coeffs);

	//! Used to check if variables are truly needed wrt to search space
	//! reduction (Chimani/Klein)
	ClusterAnalysis* m_ca;

	bool m_shrink; //!< If set to true, search space reduction is performed. Reduction is
				   //!< only feasible when only a single independent bag exists, which
				   //!< has to be assured by external partitioning.
	GraphCopy* m_ssg; //Search space graph, input graph plus edges modelled by initial variables.
	int m_nSep; //!< Stores number of separation calls
	ClusterArray<List<node> > m_cNodes; //!< Static storage of cluster node lists to avoid repeated computation.
};

}//end namespace

#endif
