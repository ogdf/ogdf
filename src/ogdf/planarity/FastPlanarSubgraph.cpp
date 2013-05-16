/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:39 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the FastPlanarSubgraph.
 *
 * Class is derived from base class PlanarSubgraphModule
 * and implements the interface for the Planarization algorithm
 * based on PQ-trees.
 *
 * \author Sebastian Leipert
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


#include <ogdf/basic/basic.h>
#include <ogdf/basic/Array.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/CriticalSection.h>
#include <ogdf/basic/Thread.h>
#include <ogdf/internal/planarity/PlanarSubgraphPQTree.h>
#include <ogdf/internal/planarity/PlanarLeafKey.h>
#include <ogdf/planarity/FastPlanarSubgraph.h>


namespace ogdf {

// default constructor
FastPlanarSubgraph::FastPlanarSubgraph()
{
	m_nRuns = 10;
}


// copy constructor
FastPlanarSubgraph::FastPlanarSubgraph(const FastPlanarSubgraph &fps) : PlanarSubgraphModule(fps)
{
	m_nRuns = fps.m_nRuns;
}


// clone method
PlanarSubgraphModule *FastPlanarSubgraph::clone() const
{
	return new FastPlanarSubgraph(*this);
}


// assignment operator
FastPlanarSubgraph &FastPlanarSubgraph::operator=(const FastPlanarSubgraph &fps)
{
	maxThreads(fps.maxThreads());

	m_timeLimit = fps.m_timeLimit;
	m_nRuns = fps.m_nRuns;
	return *this;
}


class FastPlanarSubgraph::ThreadMaster {

	Array<int>          m_bestSolution;  // value of best solution for block
	Array<List<edge> *> m_bestDelEdges;  // best solution for block

	int m_nBlocks;  // number of blocks
	const Array<BlockType> &m_block;  // the blocks (graph and edge mapping)
	const EdgeArray<int>   *m_pCost;  // edge cost (may be 0)
	__int32 m_runs;

	CriticalSection     m_criticalSection; // thread synchronization

public:
	ThreadMaster(const Array<BlockType> &block, const EdgeArray<int> *pCost, int runs);

	int numBlocks() const { return m_nBlocks; }
	const Graph &block(int i) const { return *m_block[i].first; }

	bool considerBlock(int i) const {
		return (m_bestSolution[i] > 1);
	}

	List<edge> *postNewResult(int i, List<edge> *pNewDelEdges);
	void buildSolution(List<edge> &delEdges);

	bool getNextRun() {
		return atomicDec(&m_runs) >= 0;
	}
};


class FastPlanarSubgraph::Worker : public Thread {

	ThreadMaster *m_pMaster;  // associated master

public:
	Worker(ThreadMaster *pMaster) : m_pMaster(pMaster) { }
	~Worker() { }

protected:
	virtual void doWork() {
		doWorkHelper(*m_pMaster);
	}
};


FastPlanarSubgraph::ThreadMaster::ThreadMaster(const Array<BlockType> &block, const EdgeArray<int> *pCost, int runs)
	: m_bestSolution(block.size()), m_bestDelEdges(block.size()), m_nBlocks(block.size()), m_block(block), m_pCost(pCost), m_runs(runs)
{
	for(int i = 0; i < m_nBlocks; ++i) {
		m_bestDelEdges[i] = 0;
		m_bestSolution[i] = (m_block[i].first != 0) ? numeric_limits<int>::max() : 0;
	}
}


// post a new solution and update best if better
List<edge> *FastPlanarSubgraph::ThreadMaster::postNewResult(int i, List<edge> *pNewDelEdges)
{
	int newSolution;
	if(m_pCost == 0) {
		newSolution = pNewDelEdges->size();

	} else {
		const EdgeArray<edge> &origEdge = *m_block[i].second;

		newSolution = 0;
		for(ListConstIterator<edge> it = pNewDelEdges->begin(); it.valid(); ++it)
			newSolution += (*m_pCost)[origEdge[*it]];
	}

	m_criticalSection.enter();

	if(newSolution < m_bestSolution[i]) {
		std::swap(pNewDelEdges, m_bestDelEdges[i]);
		m_bestSolution[i] = newSolution;
	}

	m_criticalSection.leave();

	return pNewDelEdges;
}


// creates a solution for the original graph from best solutions per block
// also deletes (intermediate) solution lists
void FastPlanarSubgraph::ThreadMaster::buildSolution(List<edge> &delEdges)
{
	for(int i = 0; i < m_nBlocks; ++i) {
		if(m_bestDelEdges[i] != 0) {
			const EdgeArray<edge> &origEdge = *m_block[i].second;
			for(ListConstIterator<edge> it = m_bestDelEdges[i]->begin(); it.valid(); ++it)
				delEdges.pushBack(origEdge[*it]);
			delete m_bestDelEdges[i];
		}
	}
}


// performs actual work
void FastPlanarSubgraph::doWorkHelper(ThreadMaster &master)
{
	const int nBlocks = master.numBlocks();

	do {
		for(int i = 0; i < nBlocks; ++i) {
			if(master.considerBlock(i)) {
				const Graph &B = master.block(i);

				// compute (randomized) st-numbering
				NodeArray<int> numbering(B,0);
				stNumber(B,numbering,0,0,true);

				List<edge> *pCurrentDelEdges = new List<edge>;
				planarize(B,numbering,*pCurrentDelEdges);

				//int currentSolution = master.costOf(i, pCurrentDelEdges);
				pCurrentDelEdges = master.postNewResult(i, pCurrentDelEdges);
				delete pCurrentDelEdges;
			}
		}

	} while(master.getNextRun());
}


// Prepares the planarity test and the planar embedding
Module::ReturnType FastPlanarSubgraph::doCall(
	const Graph &G,
	const List<edge> & /*preferedEdges*/,
	List<edge> &delEdges,
	const EdgeArray<int>  *pCost,
	bool /*preferedImplyPlanar*/)
{
	delEdges.clear();

	if (G.numberOfEdges() < 9)
		return retOptimal;

	// Determine Biconnected Components
	EdgeArray<int> componentID(G);
	int nBlocks = biconnectedComponents(G,componentID);

	// Determine edges per biconnected component
	Array<SList<edge> > blockEdges(0,nBlocks-1);
	edge e;
	forall_edges(e,G) {
		if (!e->isSelfLoop())
			blockEdges[componentID[e]].pushFront(e);
	}

	// Build non-trivial blocks
	Array<BlockType> block(nBlocks);
	NodeArray<node> copyV(G,0);

	for(int i = 0; i < nBlocks; i++)
	{
		if(blockEdges[i].size() < 9) {
			block[i] = BlockType((Graph*)0, (EdgeArray<edge>*)0); // explicit casts required for VS2010
			continue;
		}

		Graph *bc = new Graph;
		EdgeArray<edge> *origE = new EdgeArray<edge>(*bc,0);
		block[i] = BlockType(bc,origE);

		SList<node> marked;
		SListIterator<edge> it;
		for (it = blockEdges[i].begin(); it.valid(); ++it)
		{
			e = *it;
			if (copyV[e->source()] == 0) {
				copyV[e->source()] = bc->newNode();
				marked.pushBack(e->source());
			}
			if (copyV[e->target()] == 0) {
				copyV[e->target()] = bc->newNode();
				marked.pushBack(e->target());
			}

			(*origE)[bc->newEdge(copyV[e->source()],copyV[e->target()])] = e;
		}

		SListIterator<node> itn;
		for (itn = marked.begin(); itn.valid(); ++itn)
			copyV[*itn] = 0;
	}
	copyV.init();

	int nRuns = max(1, m_nRuns);
	int nThreads = min(maxThreads(), nRuns);

	if(nThreads == 1)
		seqCall(block, pCost, nRuns, (m_nRuns == 0), delEdges);
	else
		parCall(block, pCost, nRuns, nThreads, delEdges);

	// clean-up
	for(int i = 0; i < nBlocks; i++) {
		delete block[i].first;
		delete block[i].second;
	}

	return retFeasible;
}


//
// sequential implementation
//
void FastPlanarSubgraph::seqCall(const Array<BlockType> &block, const EdgeArray<int> *pCost, int nRuns, bool randomize, List<edge> &delEdges)
{
	const int nBlocks = block.size();

	Array<int>          bestSolution(nBlocks);
	Array<List<edge> *> bestDelEdges(nBlocks);

	for(int i = 0; i < nBlocks; ++i) {
		bestDelEdges[i] = 0;
		bestSolution[i] = (block[i].first != 0) ? numeric_limits<int>::max() : 0;
	}

	for(int run = 0; run < nRuns; ++run)
	{
		for(int i = 0; i < nBlocks; ++i) {
			if(bestSolution[i] > 1) {
				const Graph           &B        = *block[i].first;
				const EdgeArray<edge> &origEdge = *block[i].second;

				// compute (randomized) st-numbering
				NodeArray<int> numbering(B,0);
				stNumber(B,numbering,0,0,randomize);

				List<edge> *pCurrentDelEdges = new List<edge>;
				planarize(B,numbering,*pCurrentDelEdges);

				int currentSolution;
				if(pCost == 0) {
					currentSolution = pCurrentDelEdges->size();
				} else {
					currentSolution = 0;
					for(ListConstIterator<edge> it = pCurrentDelEdges->begin(); it.valid(); ++it)
						currentSolution += (*pCost)[origEdge[*it]];
				}

				if(currentSolution < bestSolution[i]) {
					delete bestDelEdges[i];
					bestDelEdges[i] = pCurrentDelEdges;
					bestSolution[i] = currentSolution;
				} else
					delete pCurrentDelEdges;
			}
		}
	}

	// build final solution from block solutions
	for(int i = 0; i < nBlocks; ++i) {
		if(bestDelEdges[i] != 0) {
			const EdgeArray<edge> &origEdge = *block[i].second;
			for(ListConstIterator<edge> it = bestDelEdges[i]->begin(); it.valid(); ++it)
				delEdges.pushBack(origEdge[*it]);
			delete bestDelEdges[i];
		}
	}
}


//
// parallel implementation
//
void FastPlanarSubgraph::parCall(const Array<BlockType> &block, const EdgeArray<int> *pCost, int nRuns, int nThreads, List<edge> &delEdges)
{
	ThreadMaster master(block, pCost, nRuns-nThreads);

	Array<Worker *> thread(nThreads-1);
	for(int i = 0; i < nThreads-1; ++i) {
		thread[i] = new Worker(&master);
		thread[i]->start();
	}

	doWorkHelper(master);

	for(int i = 0; i < nThreads-1; ++i) {
		thread[i]->join();
		delete thread[i];
	}

	master.buildSolution(delEdges);
}


// Performs a planarity test on a biconnected component
// of G. numbering contains an st-numbering of the component.
void FastPlanarSubgraph::planarize(
	const Graph &G,
	NodeArray<int> &numbering,
	List<edge> &delEdges)
{
	node v;

	NodeArray<SListPure<PlanarLeafKey<whaInfo*>* > > inLeaves(G);
	NodeArray<SListPure<PlanarLeafKey<whaInfo*>* > > outLeaves(G);
	Array<node> table(G.numberOfNodes()+1);

	forall_nodes(v,G)
	{
		edge e;
		forall_adj_edges(e,v)
		{
			if (numbering[e->opposite(v)] > numbering[v])
				// sideeffect: ignores selfloops
			{
				PlanarLeafKey<whaInfo*>* L = OGDF_NEW PlanarLeafKey<whaInfo*>(e);
				inLeaves[v].pushFront(L);
			}
		}
		table[numbering[v]] = v;
	}

	forall_nodes(v,G)
	{
		SListIterator<PlanarLeafKey<whaInfo*>* > it;
		for (it = inLeaves[v].begin(); it.valid(); ++it)
		{
			PlanarLeafKey<whaInfo*>* L = *it;
			outLeaves[L->userStructKey()->opposite(v)].pushFront(L);
		}
	}

	SList<PQLeafKey<edge,whaInfo*,bool>*> totalEliminatedKeys;

	PlanarSubgraphPQTree T;
	T.Initialize(inLeaves[table[1]]);
	for (int i = 2; i < G.numberOfNodes(); i++)
	{
		SList<PQLeafKey<edge,whaInfo*,bool>*> eliminatedKeys;
		T.Reduction(outLeaves[table[i]],eliminatedKeys);

		totalEliminatedKeys.conc(eliminatedKeys);
		T.ReplaceRoot(inLeaves[table[i]]);
		T.emptyAllPertinentNodes();
	}


	SListIterator<PQLeafKey<edge,whaInfo*,bool>* > it;
	for (it = totalEliminatedKeys.begin(); it.valid(); ++it)
	{
		edge e = (*it)->userStructKey();
		delEdges.pushBack(e);
	}

	//cleanup
	forall_nodes(v,G)
	{
		while (!inLeaves[v].empty())
		{
			PlanarLeafKey<whaInfo*>* L = inLeaves[v].popFrontRet();
			delete L;
		}
	}

	T.Cleanup();	// Explicit call for destructor necessary. This allows to call virtual
					// funtion CleanNode for freeing node information class.
}


}
