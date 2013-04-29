/*
 * $Revision: 3472 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-29 15:52:12 +0200 (Mo, 29. Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implements class SubgraphPlanarizerUML.
 *
 * \author Markus Chimani, Carsten Gutwenger
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

#include <ogdf/uml/SubgraphPlanarizerUML.h>
#include <ogdf/uml/VariableEmbeddingInserterUML.h>
#include <ogdf/planarity/MaximalPlanarSubgraphSimple.h>
#include <ogdf/basic/CriticalSection.h>
#include <ogdf/basic/Thread.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/internal/planarity/CrossingStructure.h>


#ifdef OGDF_HAVE_CPP11
using std::minstd_rand;
#endif


namespace ogdf
{

class SubgraphPlanarizerUML::ThreadMaster {

	CrossingStructure *m_pCS;
	int                m_bestCR;

	const PlanRep     &m_pr;
	int                m_cc;

	const EdgeArray<int> *m_pCost;
	const List<edge>     &m_delEdges;

	int m_seed;
	__int32               m_perms;
	__int64               m_stopTime;
	CriticalSection       m_criticalSection;

public:
	ThreadMaster(
		const PlanRep &pr,
		int cc,
		const EdgeArray<int>  *pCost,
		const List<edge> &delEdges,
		int seed,
		int perms,
		__int64 stopTime);

	~ThreadMaster() { delete m_pCS; }

	const PlanRep &planRep() const { return m_pr; }
	int currentCC() const { return m_cc; }

	const EdgeArray<int> *cost() const { return m_pCost; }
	const List<edge> &delEdges() const { return m_delEdges; }

	int rseed(long id) const { return (int)id * m_seed; }

	int queryBestKnown() const { return m_bestCR; }
	CrossingStructure *postNewResult(CrossingStructure *pCS);
	bool getNextPerm();

	void restore(PlanRep &PG, int &cr);
};


class SubgraphPlanarizerUML::Worker : public Thread {

	ThreadMaster *m_pMaster;
	UMLEdgeInsertionModule *m_pInserter;

public:
	Worker(ThreadMaster *pMaster, UMLEdgeInsertionModule *pInserter) : m_pMaster(pMaster), m_pInserter(pInserter) { }
	~Worker() { delete m_pInserter; }

protected:
	virtual void doWork();
};


SubgraphPlanarizerUML::ThreadMaster::ThreadMaster(
	const PlanRep &pr,
	int cc,
	const EdgeArray<int>  *pCost,
	const List<edge> &delEdges,
	int seed,
	int perms,
	__int64 stopTime)
	:
	m_pCS(0), m_bestCR(numeric_limits<int>::max()), m_pr(pr), m_cc(cc),
	m_pCost(pCost),
	m_delEdges(delEdges), m_seed(seed), m_perms(perms), m_stopTime(stopTime)
{ }


CrossingStructure *SubgraphPlanarizerUML::ThreadMaster::postNewResult(CrossingStructure *pCS)
{
	int newCR = pCS->numberOfCrossings();

	m_criticalSection.enter();

	if(newCR < m_bestCR) {
		std::swap(pCS, m_pCS);
		m_bestCR = newCR;
	}

	m_criticalSection.leave();

	return pCS;
}


bool SubgraphPlanarizerUML::ThreadMaster::getNextPerm()
{
	if(m_stopTime >= 0 && System::realTime() >= m_stopTime)
		return false;
	return atomicDec(&m_perms) >= 0;
}


void SubgraphPlanarizerUML::ThreadMaster::restore(PlanRep &PG, int &cr)
{
	m_pCS->restore(PG, m_cc);
	cr = m_bestCR;
}


bool SubgraphPlanarizerUML::doSinglePermutation(
	PlanRepLight &PG,
	int cc,
	const EdgeArray<int>  *pCost,
	Array<edge> &deletedEdges,
	UMLEdgeInsertionModule &inserter,
#ifdef OGDF_HAVE_CPP11
	minstd_rand &rng,
#endif
	int &crossingNumber
	)
{
	PG.initCC(cc);

	const int nG = PG.numberOfNodes();
	const int high = deletedEdges.high();

	for(int j = 0; j <= high; ++j)
		PG.delEdge(PG.copy(deletedEdges[j]));

	// permute
#ifdef OGDF_HAVE_CPP11
	std::uniform_int_distribution<int> dist(0,high);

	for(int j = 0; j <= high; ++j)
		deletedEdges.swap(j, dist(rng));
#else
	for(int j = 0; j <= high; ++j)
		deletedEdges.swap(j, randomNumber(0,high));
#endif

	ReturnType ret = inserter.callEx(PG, deletedEdges, pCost, 0);

	if(isSolution(ret) == false)
		return false; // no solution found, that's bad...

	if(pCost == 0)
		crossingNumber = PG.numberOfNodes() - nG;
	else {
		crossingNumber = 0;
		node n;
		forall_nodes(n, PG) {
			if(PG.original(n) == 0) { // dummy found -> calc cost
				edge e1 = PG.original(n->firstAdj()->theEdge());
				edge e2 = PG.original(n->lastAdj()->theEdge());
				crossingNumber += (*pCost)[e1] * (*pCost)[e2];
			}
		}
	}

	return true;
}

void SubgraphPlanarizerUML::doWorkHelper(ThreadMaster &master, UMLEdgeInsertionModule &inserter
#ifdef OGDF_HAVE_CPP11
	, minstd_rand &rng
#endif
)
{
	const List<edge> &delEdges = master.delEdges();

	const int m = delEdges.size();
	Array<edge> deletedEdges(m);
	int j = 0;
	for(ListConstIterator<edge> it = delEdges.begin(); it.valid(); ++it)
		deletedEdges[j++] = *it;

	PlanRepLight PG(master.planRep());
	int cc = master.currentCC();

	const EdgeArray<int>  *pCost = master.cost();

	do {
		int crossingNumber;
		if(doSinglePermutation(PG, cc, pCost, deletedEdges, inserter,
#ifdef OGDF_HAVE_CPP11
			rng,
#endif
			crossingNumber)
			&& crossingNumber < master.queryBestKnown())
		{
			CrossingStructure *pCS = new CrossingStructure;
			pCS->init(PG, crossingNumber);
			pCS = master.postNewResult(pCS);
			delete pCS;
		}

	} while(master.getNextPerm());
}


void SubgraphPlanarizerUML::Worker::doWork()
{
#ifdef OGDF_HAVE_CPP11
	minstd_rand rng(m_pMaster->rseed(threadID())); // different seeds per thread
	doWorkHelper(*m_pMaster, *m_pInserter, rng);
#else
	doWorkHelper(*m_pMaster, *m_pInserter);
#endif
}


// default constructor
SubgraphPlanarizerUML::SubgraphPlanarizerUML()
{
	m_subgraph.set(new MaximalPlanarSubgraphSimple);
	m_inserter.set(new VariableEmbeddingInserterUML);

	m_permutations = 1;
	m_setTimeout = true;

#ifdef OGDF_MEMORY_POOL_NTS
	m_maxThreads = 1;
#else
	m_maxThreads = System::numberOfProcessors();
#endif
}


// copy constructor
SubgraphPlanarizerUML::SubgraphPlanarizerUML(const SubgraphPlanarizerUML &planarizer)
	: UMLCrossingMinimizationModule(planarizer), Logger()
{
	m_subgraph.set(planarizer.m_subgraph.get().clone());
	m_inserter.set(planarizer.m_inserter.get().clone());

	m_permutations = planarizer.m_permutations;
	m_setTimeout   = planarizer.m_setTimeout;
	m_maxThreads   = planarizer.m_maxThreads;
}


// clone method
UMLCrossingMinimizationModule *SubgraphPlanarizerUML::clone() const {
	return new SubgraphPlanarizerUML(*this);
}


// assignment operator
SubgraphPlanarizerUML &SubgraphPlanarizerUML::operator=(const SubgraphPlanarizerUML &planarizer)
{
	m_timeLimit = planarizer.m_timeLimit;
	m_subgraph.set(planarizer.m_subgraph.get().clone());
	m_inserter.set(planarizer.m_inserter.get().clone());

	m_permutations = planarizer.m_permutations;
	m_setTimeout   = planarizer.m_setTimeout;
	m_maxThreads   = planarizer.m_maxThreads;

	return *this;
}


Module::ReturnType SubgraphPlanarizerUML::doCall(
	PlanRepUML           &pr,
	int                   cc,
	const EdgeArray<int> *pCostOrig,
	int                  &crossingNumber)
{
	OGDF_ASSERT(m_permutations >= 1);

	PlanarSubgraphModule   &subgraph = m_subgraph.get();
	UMLEdgeInsertionModule &inserter = m_inserter.get();

	int nThreads = min(m_maxThreads, m_permutations);

	__int64 startTime;
	System::usedRealTime(startTime);
	__int64 stopTime = (m_timeLimit >= 0) ? (startTime + __int64(1000.0*m_timeLimit)) : -1;

	//
	// Compute subgraph
	//
	if(m_setTimeout)
		subgraph.timeLimit(m_timeLimit);

	pr.initCC(cc);

	// gather generalization edges, which should all be in the planar subgraph
	List<edge> preferedEdges;
	edge e;
	forall_edges(e,pr) {
		if (pr.typeOf(e) == Graph::generalization)
			preferedEdges.pushBack(e);
	}

	List<edge> delEdges;
	ReturnType retValue;

	if(pCostOrig) {
		EdgeArray<int> costPG(pr);
		forall_edges(e,pr)
			costPG[e] = (*pCostOrig)[pr.original(e)];

		retValue = subgraph.call(pr, costPG, preferedEdges, delEdges);

	} else
		retValue = subgraph.call(pr, preferedEdges, delEdges);

	if(isSolution(retValue) == false)
		return retValue;

	const int m = delEdges.size();
	for(ListIterator<edge> it = delEdges.begin(); it.valid(); ++it)
		*it = pr.original(*it);

	//
	// Permutation phase
	//

	int seed = rand();
#ifdef OGDF_HAVE_CPP11
	minstd_rand rng(seed);
#endif

	if(nThreads > 1) {
		//
		// Parallel implementation
		//
		ThreadMaster master(
			pr, cc,
			pCostOrig,
			delEdges,
			seed,
			m_permutations - nThreads,
			stopTime);

		Array<Worker *> thread(nThreads-1);
		for(int i = 0; i < nThreads-1; ++i) {
			thread[i] = new Worker(&master, inserter.clone());
			thread[i]->start();
		}

#ifdef OGDF_HAVE_CPP11
		doWorkHelper(master, inserter, rng);
#else
		doWorkHelper(master, inserter);
#endif

		for(int i = 0; i < nThreads-1; ++i) {
			thread[i]->join();
			delete thread[i];
		}

		master.restore(pr, crossingNumber);

	} else {
		//
		// Sequential implementation
		//
		PlanRepLight prl(pr);

		Array<edge> deletedEdges(m);
		int j = 0;
		for(ListIterator<edge> it = delEdges.begin(); it.valid(); ++it)
			deletedEdges[j++] = *it;

		bool foundSolution = false;
		CrossingStructure cs;
		for(int i = 1; i <= m_permutations; ++i)
		{
			int cr;
			bool ok = doSinglePermutation(prl, cc, pCostOrig, deletedEdges, inserter,
#ifdef OGDF_HAVE_CPP11
				rng,
#endif
				cr);

			if(ok && (foundSolution == false || cr < cs.weightedCrossingNumber())) {
				foundSolution = true;
				cs.init(prl, cr);
			}

			if(stopTime >= 0 && System::realTime() >= stopTime) {
				if(foundSolution == false)
					return retTimeoutInfeasible; // not able to find a solution...
				break;
			}
		}

		cs.restore(pr,cc); // restore best solution in PG
		crossingNumber = cs.weightedCrossingNumber();

		OGDF_ASSERT(isPlanar(pr) == true);
	}

	return retFeasible;
}


} // namspace ogdf
