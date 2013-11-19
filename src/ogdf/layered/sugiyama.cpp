/*
 * $Revision: 3832 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 11:16:27 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of Sugiyama algorithm (classes Hierarchy,
 * Level, SugiyamaLayout)
 *
 * \author Carsten Gutwenger
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


#include <ogdf/layered/Hierarchy.h>
#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/LongestPathRanking.h>
#include <ogdf/layered/BarycenterHeuristic.h>
#include <ogdf/layered/SplitHeuristic.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/OptimalHierarchyClusterLayout.h>
#include <ogdf/packing/TileToRowsCCPacker.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/Thread.h>
#include <ogdf/basic/CriticalSection.h>


#ifdef OGDF_HAVE_CPP11
#include <random>
using std::minstd_rand;
#endif


namespace ogdf {


//---------------------------------------------------------
// GraphCopyAttributes
// manages access on copy of an attributed graph
//---------------------------------------------------------

void GraphCopyAttributes::transform()
{
	node v;
	forall_nodes(v,*m_pGC)
	{
		node vG = m_pGC->original(v);
		if(vG) {
			m_pAG->x(vG) = m_x[v];
			m_pAG->y(vG) = m_y[v];
		}
	}

	edge e;
	forall_edges(e,*m_pGC)
	{
		edge eG = m_pGC->original(e);
		if(eG == 0 || e != m_pGC->chain(eG).front())
			continue;
		// current implementation does not layout self-loops;
		// they are simply ignored
		//if (e->isSelfLoop()) continue;

		DPolyline &dpl = m_pAG->bends(eG);
		dpl.clear();

		ListConstIterator<edge> itE = m_pGC->chain(eG).begin();
		node v      = (*itE)->source();
		node vAfter = (*itE)->target();

		for (++itE; itE.valid(); ++itE)
		{
			node vBefore = v;
			v      = vAfter;
			vAfter = (*itE)->target();

			// filter real bend points
			if((m_x[v] != m_x[vBefore] || m_x[v] != m_x[vAfter]) &&
				(m_y[v] != m_y[vBefore] || m_y[v] != m_y[vAfter]))
				dpl.pushBack(DPoint(m_x[v],m_y[v]));
		}

		if (m_pGC->isReversed(eG))
			dpl.reverse();
	}
}


//---------------------------------------------------------
// ClusterGraphCopyAttributes
// manages access on copy of an attributed cluster graph
//---------------------------------------------------------

void ClusterGraphCopyAttributes::transform()
{
	node v;
	forall_nodes(v,*m_pH)
	{
		node vG = m_pH->origNode(v);
		if(vG) {
			m_pACG->x(vG) = m_x[v];
			m_pACG->y(vG) = m_y[v];
		}
	}

	edge e;
	forall_edges(e,*m_pH)
	{
		edge eG = m_pH->origEdge(e);
		if(eG == 0 || e != m_pH->chain(eG).front())
			continue;
		// current implementation does not layout self-loops;
		// they are simply ignored
		//if (e->isSelfLoop()) continue;

		DPolyline &dpl = m_pACG->bends(eG);
		dpl.clear();

		ListConstIterator<edge> itE = m_pH->chain(eG).begin();
		node v      = (*itE)->source();
		node vAfter = (*itE)->target();

		for (++itE; itE.valid(); ++itE)
		{
			node vBefore = v;
			v      = vAfter;
			vAfter = (*itE)->target();

			// filter real bend points
			if((m_x[v] != m_x[vBefore] || m_x[v] != m_x[vAfter]) &&
				(m_y[v] != m_y[vBefore] || m_y[v] != m_y[vAfter]))
				dpl.pushBack(DPoint(m_x[v],m_y[v]));
		}

		if (m_pH->isReversed(eG))
			dpl.reverse();
	}
}


//---------------------------------------------------------
// Level
// representation of levels in hierarchies
//---------------------------------------------------------

const Array<node> &Level::adjNodes(node v) const
{
	return m_pLevels->adjNodes(v);
}


void Level::swap(int i, int j)
{
	m_nodes.swap(i,j);
	m_pLevels->m_pos[m_nodes[i]] = i;
	m_pLevels->m_pos[m_nodes[j]] = j;
}


void Level::recalcPos()
{
	NodeArray<int> &pos = m_pLevels->m_pos;

	for(int i = 0; i <= high(); ++i)
		pos[m_nodes[i]] = i;

	m_pLevels->buildAdjNodes(m_index);
}


void Level::getIsolatedNodes(SListPure<Tuple2<node,int> > &isolated) const
{
	for (int i = 0; i <= high(); ++i)
		if (adjNodes(m_nodes[i]).high() < 0)
			isolated.pushBack(Tuple2<node,int>(m_nodes[i],i));
}


void Level::setIsolatedNodes(SListPure<Tuple2<node,int> > &isolated)
{
	const int sizeL = size();
	Array<node> sortedNodes(sizeL);
	isolated.pushBack(Tuple2<node,int>(0,sizeL));
	SListConstIterator<Tuple2<node,int> > itIsolated = isolated.begin();

	int nextPos = (*itIsolated).x2();
	for( int iNodes = 0, iSortedNodes = 0; nextPos <= sizeL; ) {
		if( iSortedNodes == nextPos ) {
			if( iSortedNodes == sizeL )
				break;
			sortedNodes[iSortedNodes++] = (*itIsolated).x1();
			nextPos = (*(++itIsolated)).x2();
		} else {
			node v = m_nodes[iNodes++];
			if( adjNodes(v).size() > 0 )
				sortedNodes[iSortedNodes++] = v;
		}
	}

	for( int i = 0; i < sizeL; ++i)
		m_nodes[i] = sortedNodes[i];
}


class WeightBucket : public BucketFunc<node> {
	const NodeArray<int> *m_pWeight;

public:
	WeightBucket(const NodeArray<int> *pWeight) : m_pWeight(pWeight) { }

	int getBucket(const node &v) { return (*m_pWeight)[v]; }
};


void Level::sort(NodeArray<double> &weight)
{
	SListPure<Tuple2<node,int> > isolated;
	getIsolatedNodes(isolated);

	WeightComparer<> cmp(&weight);
	std::stable_sort(&m_nodes[0], &m_nodes[0]+m_nodes.size(), cmp);

	if (!isolated.empty()) setIsolatedNodes(isolated);
	recalcPos();
}


void Level::sortByWeightOnly(NodeArray<double> &weight)
{
	WeightComparer<> cmp(&weight);
	std::stable_sort(&m_nodes[0], &m_nodes[0]+m_nodes.size(), cmp);
	recalcPos();
}


void Level::sort(NodeArray<int> &weight, int minBucket, int maxBucket)
{
	SListPure<Tuple2<node,int> > isolated;
	getIsolatedNodes(isolated);

	WeightBucket bucketFunc(&weight);
	bucketSort(m_nodes,minBucket,maxBucket,bucketFunc);

	if (!isolated.empty()) setIsolatedNodes(isolated);
	recalcPos();
}



//---------------------------------------------------------
// Hierarchy
// representation of proper hierarchies used by Sugiyama
//---------------------------------------------------------

Hierarchy::Hierarchy(const Graph &G, const NodeArray<int> &rank) :
	m_GC(G), m_rank(m_GC)
{
	doInit(rank);
}


void Hierarchy::createEmpty(const Graph &G)
{
	m_GC.createEmpty(G);
	m_rank.init(m_GC);
}


void Hierarchy::initByNodes(const List<node> &nodes,
	EdgeArray<edge> &eCopy,
	const NodeArray<int> &rank)
{
	m_GC.initByNodes(nodes,eCopy);

	doInit(rank);
}


void Hierarchy::doInit(const NodeArray<int> &rank)
{
	makeLoopFree(m_GC);

	int maxRank = 0;

	node v;
	forall_nodes(v,m_GC) {
		int r = m_rank[v] = rank[m_GC.original(v)];
		OGDF_ASSERT(r >= 0)
		if (r > maxRank) maxRank = r;
	}

	SListPure<edge> edges;
	m_GC.allEdges(edges);
	SListConstIterator<edge> it;
	for(it = edges.begin(); it.valid(); ++it)
	{
		edge e = *it;

		int rankSrc = m_rank[e->source()], rankTgt = m_rank[e->target()];

		if (rankSrc > rankTgt) {
			m_GC.reverseEdge(e); std::swap(rankSrc,rankTgt);
		}

		if (rankSrc == rankTgt) {
			e = m_GC.split(e);
			m_GC.reverseEdge(e);
			if ((m_rank[e->target()] = rankSrc+1) > maxRank)
				maxRank = rankSrc+1;

		} else {
			for(++rankSrc; rankSrc < rankTgt; ++rankSrc)
				m_rank[(e = m_GC.split(e))->source()] = rankSrc;
		}
	}

	m_size.init(0,maxRank,0);
	forall_nodes(v,m_GC)
		m_size[m_rank[v]]++;
}


//---------------------------------------------------------
// HierarchyLevels
// stores and maintains the levels in a hierarchy
//---------------------------------------------------------

HierarchyLevels::HierarchyLevels(const Hierarchy &H) : m_H(H), m_pLevel(0,H.maxRank()), m_pos(H), m_lowerAdjNodes(H), m_upperAdjNodes(H), m_nSet(H,0)
{
	const GraphCopy &GC = m_H;
	int maxRank = H.maxRank();

	for(int i = 0; i <= maxRank; ++i)
		m_pLevel[i] = new Level(this,i,H.size(i));

	Array<int> next(0,maxRank,0);

	node v;
	forall_nodes(v,GC) {
		int r = H.rank(v), pos = next[r]++;
		m_pos[(*m_pLevel[r])[pos] = v] = pos;

		m_lowerAdjNodes[v].init(v->indeg());
		m_upperAdjNodes[v].init(v->outdeg());
	}

	buildAdjNodes();
}


HierarchyLevels::~HierarchyLevels()
{
	for(int i = 0; i <= high(); ++i)
		delete m_pLevel[i];
}


void HierarchyLevels::buildAdjNodes()
{
	for(int i = 0; i <= high(); ++i)
		buildAdjNodes(i);
}


void HierarchyLevels::buildAdjNodes(int i)
{
	if (i > 0) {
		const Level &lowerLevel = *m_pLevel[i-1];

		for(int j = 0; j <= lowerLevel.high(); ++j)
			m_nSet[lowerLevel[j]] = 0;
	}

	if (i < high()) {
		const Level &upperLevel = *m_pLevel[i+1];

		for(int j = 0; j <= upperLevel.high(); ++j)
			m_nSet[upperLevel[j]] = 0;
	}

	const Level &level = *m_pLevel[i];

	for(int j = 0; j <= level.high(); ++j) {
		node v = level[j];
		edge e;
		forall_adj_edges(e,v) {
			if (e->source() == v) {
				(m_lowerAdjNodes[e->target()])[m_nSet[e->target()]++] = v;
			} else {
				(m_upperAdjNodes[e->source()])[m_nSet[e->source()]++] = v;
			}
		}
	}
}


void HierarchyLevels::storePos (NodeArray<int> &oldPos) const
{
	oldPos = m_pos;
}


void HierarchyLevels::restorePos (const NodeArray<int> &newPos)
{
	const GraphCopy &GC = m_H;

	m_pos = newPos;

	node v;
	forall_nodes(v,GC) {
		(*m_pLevel[m_H.rank(v)])[m_pos[v]] = v;
	}

	//check();

	buildAdjNodes();
}


void HierarchyLevels::permute()
{
	for(int i = 0; i < m_pLevel.high(); ++i) {
		Level &level = *m_pLevel[i];
		level.m_nodes.permute();
		for(int j = 0; j <= level.high(); ++j)
			m_pos[level[j]] = j;
	}

	//check();

	buildAdjNodes();
}


void HierarchyLevels::separateCCs(int numCC, const NodeArray<int> &component)
{
	Array<SListPure<node> > table(numCC);

	for(int i = 0; i < m_pLevel.high(); ++i) {
		Level &level = *m_pLevel[i];
		for(int j = 0; j <= level.high(); ++j) {
			node v = level[j];
			table[component[v]].pushBack(v);
		}
	}

	Array<int> count(0, m_pLevel.high(), 0);
	for(int c = 0; c < numCC; ++c) {
		SListConstIterator<node> it;
		for(it = table[c].begin(); it.valid(); ++it)
			m_pos[*it] = count[m_H.rank(*it)]++;
	}

	const GraphCopy &GC = m_H;

	node v;
	forall_nodes(v,GC) {
		(*m_pLevel[m_H.rank(v)])[m_pos[v]] = v;
	}

	//check();

	buildAdjNodes();
}


int HierarchyLevelsBase::calculateCrossings() const
{
	int nCrossings = 0;

	for(int i = 0; i < this -> high(); ++i) {
		nCrossings += calculateCrossings(i);
	}

	return nCrossings;
}


// calculation of edge crossings between level i and i+1
// implementation by Michael Juenger, Decembre 2000, adapted by Carsten Gutwenger
// implements the algorithm by Barth, Juenger, Mutzel

int HierarchyLevelsBase::calculateCrossings(int i) const
{
	// const Level &L = *m_pLevel[i];             // level i
	// const int nUpper = m_pLevel[i+1]->size();  // number of nodes on level i+1

	const LevelBase &L = (*this)[i];             // level i
	const int nUpper = (*this)[i+1].size();  // number of nodes on level i+1

	int nc = 0; // number of crossings

	int fa = 1;
	while (fa < nUpper)
		fa *= 2;

	int nTreeNodes = 2*fa - 1; // number of tree nodes
	fa -= 1;         // "first address:" indexincrement in tree

	Array<int> nin(0,nTreeNodes-1,0);

	for(int j = 0; j < L.size(); ++j)
	{
		// const Array<node> &adjNodes = m_upperAdjNodes[L[j]];
		const Array<node> &adjNodes = this->adjNodes(L[j], upward); // m_upperAdjNodes[L[j]];

		for(int k = 0; k < adjNodes.size(); ++k)
		{
			// index of tree node for vertex adjNode[k]
			//int index = m_pos[adjNodes[k]] + fa;
			int index = pos(adjNodes[k]) + fa;
			nin[index]++;

			while (index>0) {
				if (index % 2)
					nc += nin[index+1]; // new crossing
				index = (index - 1) / 2;
				nin[index]++;
			}

		}
	}

	return nc;
}


int HierarchyLevels::calculateCrossingsSimDraw(const EdgeArray<__uint32> *edgeSubGraphs) const
{
	int nCrossings = 0;

	for(int i = 0; i < m_pLevel.high(); ++i) {
		nCrossings += calculateCrossingsSimDraw(i, edgeSubGraphs);
	}

	return nCrossings;
}


// naive calculation of edge crossings between level i and i+1
// for SimDraw-calculation by Michael Schulz

int HierarchyLevels::calculateCrossingsSimDraw(int i, const EdgeArray<__uint32> *edgeSubGraphs) const
{
	const int maxGraphs = 32;

	const Level &L = *m_pLevel[i];             // level i
	const GraphCopy &GC = m_H;

	int nc = 0; // number of crossings

	for(int j = 0; j < L.size(); ++j)
	{
		node v = L[j];
		edge e;
		forall_adj_edges(e,v) {
			if (e->source() == v){
				int pos_adj_e = pos(e->target());
				for (int k = j+1; k < L.size(); k++) {
					node w = L[k];
					edge f;
					forall_adj_edges(f,w) {
						if (f->source() == w) {
							int pos_adj_f = pos(f->target());
							if(pos_adj_f < pos_adj_e)
							{
								int graphCounter = 0;
								for(int numGraphs = 0; numGraphs < maxGraphs; numGraphs++)
									if((1 << numGraphs) & (*edgeSubGraphs)[GC.original(e)] & (*edgeSubGraphs)[GC.original(f)])
										graphCounter++;
								nc += graphCounter;
							}
						}
					}
				}
			}
		}
	}

	return nc;
}


int HierarchyLevels::transposePart(
	const Array<node> &adjV,
	const Array<node> &adjW)
{
	const int vSize = adjV.size();
	int iW = 0, iV = 0, sum = 0;

	for(; iW <= adjW.high(); ++iW) {
		int p = m_pos[adjW[iW]];
		while(iV < vSize && m_pos[adjV[iV]] <= p) ++iV;
		sum += vSize - iV;
	}

	return sum;
}


bool HierarchyLevels::transpose(node v)
{
	int rankV = m_H.rank(v), posV = m_pos[v];
	node w = (*m_pLevel[rankV])[posV+1];

	int d = 0;
	d += transposePart(m_upperAdjNodes[v],m_upperAdjNodes[w]);
	d -= transposePart(m_upperAdjNodes[w],m_upperAdjNodes[v]);
	d += transposePart(m_lowerAdjNodes[v],m_lowerAdjNodes[w]);
	d -= transposePart(m_lowerAdjNodes[w],m_lowerAdjNodes[v]);

	if (d > 0) {
		m_pLevel[rankV]->swap(posV,posV+1);
		return true;
	}

	return false;
}


void HierarchyLevels::print(ostream &os) const
{
	for(int i = 0; i <= m_pLevel.high(); ++i) {
		os << i << ": ";
		const Level &level = *m_pLevel[i];
		for(int j = 0; j < level.size(); ++j)
			os << level[j] << " ";
		os << endl;
	}

	os << endl;

	const GraphCopy &GC = m_H;

	node v;
	forall_nodes(v,GC) {
		os << v << ": lower: " << (m_lowerAdjNodes[v]) <<
			", upper: " << (m_upperAdjNodes[v]) << endl;
	}

}


void HierarchyLevels::check() const
{
	int i, j;
	for(i = 0; i <= high(); ++i) {
		Level &L = *m_pLevel[i];
		for(j = 0; j <= L.high(); ++j) {
			if (m_pos[L[j]] != j) {
				cerr << "m_pos[" << L[j] << "] wrong!" << endl;
			}
			if (m_H.rank(L[j]) != i) {
				cerr << "m_rank[" << L[j] << "] wrong!" << endl;
			}
		}
	}
}


//---------------------------------------------------------
// LayerByLayerSweep
// Crossing reduction algorithm using 2-level heuristic.
//---------------------------------------------------------

// LayerByLayerSweep::CrossMinMaster

class LayerByLayerSweep::CrossMinMaster {

	NodeArray<int>  *m_pBestPos;
	int              m_bestCR;

	const SugiyamaLayout &m_sugi;
	const Hierarchy      &m_H;

	int              m_seed;
	__int32          m_runs;
	CriticalSection  m_criticalSection;

public:
	CrossMinMaster(
		const SugiyamaLayout &sugi,
		const Hierarchy &H,
		int seed,
		int runs);

	const Hierarchy &hierarchy() const { return m_H; }
	int rseed(long id) const { return (int)id * m_seed; }

	void restore(HierarchyLevels &levels, int &cr);

	void doWorkHelper(
		LayerByLayerSweep        *pCrossMin,
		TwoLayerCrossMinSimDraw *pCrossMinSimDraw,
		HierarchyLevels         &levels,
		NodeArray<int>          &bestPos,
		bool                     permuteFirst
#ifdef OGDF_HAVE_CPP11
		, std::minstd_rand      &rng
#endif
		);

private:
	const EdgeArray<__uint32> *subgraphs() const { return m_sugi.subgraphs(); }
	int fails() const { return m_sugi.fails(); }
	bool transpose() { return m_sugi.transpose(); }

	bool arrangeCCs() const { return m_sugi.arrangeCCs(); }
	int arrange_numCC() const { return m_sugi.numCC(); }
	const NodeArray<int> &arrange_compGC() const { return m_sugi.compGC(); }

	bool transposeLevel(int i, HierarchyLevels &levels, Array<bool> &levelChanged);
	void doTranspose(HierarchyLevels &levels, Array<bool> &levelChanged);
	void doTransposeRev(HierarchyLevels &levels, Array<bool> &levelChanged);

	int traverseTopDown(
		HierarchyLevels &levels,
		LayerByLayerSweep *pCrossMin,
		TwoLayerCrossMinSimDraw *pCrossMinSimDraw,
		Array<bool>             *pLevelChanged);

	int traverseBottomUp(
		HierarchyLevels &levels,
		LayerByLayerSweep *pCrossMin,
		TwoLayerCrossMinSimDraw *pCrossMinSimDraw,
		Array<bool>             *pLevelChanged);

	int queryBestKnown() const { return m_bestCR; }
	bool postNewResult(int cr, NodeArray<int> *pPos);
	bool getNextRun();
};



LayerByLayerSweep::CrossMinMaster::CrossMinMaster(
	const SugiyamaLayout &sugi,
	const Hierarchy &H,
	int seed,
	int runs)
	: m_pBestPos(0), m_bestCR(numeric_limits<int>::max()), m_sugi(sugi), m_H(H), m_seed(seed), m_runs(runs) { }


bool LayerByLayerSweep::CrossMinMaster::postNewResult(int cr, NodeArray<int> *pPos)
{
	bool storeResult = false;
	m_criticalSection.enter();

	if(cr < m_bestCR) {
		m_bestCR = cr;
		m_pBestPos = pPos;
		storeResult = true;

		if(cr == 0)
			m_runs = 0;
	}

	m_criticalSection.leave();
	return storeResult;
}


bool LayerByLayerSweep::CrossMinMaster::getNextRun()
{
	return atomicDec(&m_runs) >= 0;
}


void LayerByLayerSweep::CrossMinMaster::restore(HierarchyLevels &levels, int &cr)
{
	levels.restorePos(*m_pBestPos);
	cr = m_bestCR;
}


bool LayerByLayerSweep::CrossMinMaster::transposeLevel(int i, HierarchyLevels &levels, Array<bool> &levelChanged)
{
	bool improved = false;

	if (levelChanged[i] || levelChanged[i-1] || levelChanged[i+1]) {
		Level &L = levels[i];

		for (int j = 0; j < L.high(); j++) {
			if (levels.transpose(L[j])) improved = true;
		}
	}

	if (improved) levels.buildAdjNodes(i);
	return (levelChanged[i] = improved);
}


void LayerByLayerSweep::CrossMinMaster::doTranspose(HierarchyLevels &levels, Array<bool> &levelChanged)
{
	levelChanged.fill(true);

	bool improved;
	do {
		improved = false;

		for (int i = 0; i <= levels.high(); ++i)
			improved |= transposeLevel(i,levels,levelChanged);
	} while (improved);
}


void LayerByLayerSweep::CrossMinMaster::doTransposeRev(HierarchyLevels &levels, Array<bool> &levelChanged)
{
	levelChanged.fill(true);

	bool improved;
	do {
		improved = false;

		for (int i = levels.high(); i >= 0 ; --i)
			improved |= transposeLevel(i,levels,levelChanged);
	} while (improved);
}


int LayerByLayerSweep::CrossMinMaster::traverseTopDown(
	HierarchyLevels           &levels,
	LayerByLayerSweep          *pCrossMin,
	TwoLayerCrossMinSimDraw   *pCrossMinSimDraw,
	Array<bool>               *pLevelChanged)
{
	levels.direction(HierarchyLevels::downward);

	for (int i = 1; i <= levels.high(); ++i) {
		if(pCrossMin != 0)
			pCrossMin->call(levels[i]);
		else
			pCrossMinSimDraw->call(levels[i], subgraphs());
	}

	if(pLevelChanged != 0)
		doTranspose(levels, *pLevelChanged);
	if(arrangeCCs() == false)
		levels.separateCCs(arrange_numCC(), arrange_compGC());

	return (pCrossMin != 0) ? levels.calculateCrossings() : levels.calculateCrossingsSimDraw(subgraphs());
}


int LayerByLayerSweep::CrossMinMaster::traverseBottomUp(
	HierarchyLevels           &levels,
	LayerByLayerSweep          *pCrossMin,
	TwoLayerCrossMinSimDraw   *pCrossMinSimDraw,
	Array<bool>               *pLevelChanged)
{
	levels.direction(HierarchyLevels::upward);

	for (int i = levels.high()-1; i >= 0; i--) {
		if(pCrossMin != 0)
			pCrossMin->call(levels[i]);
		else
			pCrossMinSimDraw->call(levels[i], subgraphs());
	}

	if (pLevelChanged != 0)
		doTransposeRev(levels, *pLevelChanged);
	if(arrangeCCs() == false)
		levels.separateCCs(arrange_numCC(), arrange_compGC());

	return (pCrossMin != 0) ? levels.calculateCrossings() : levels.calculateCrossingsSimDraw(subgraphs());
}


void LayerByLayerSweep::CrossMinMaster::doWorkHelper(
	LayerByLayerSweep        *pCrossMin,
	TwoLayerCrossMinSimDraw *pCrossMinSimDraw,
	HierarchyLevels         &levels,
	NodeArray<int>          &bestPos,
	bool                     permuteFirst
#ifdef OGDF_HAVE_CPP11
	, minstd_rand &rng
#endif
	)
{
	if(permuteFirst)
#ifdef OGDF_HAVE_CPP11
		levels.permute(rng);
#else
		levels.permute();
#endif

	int nCrossingsOld = (pCrossMin != 0) ? levels.calculateCrossings() : levels.calculateCrossingsSimDraw(subgraphs());
	if(postNewResult(nCrossingsOld, &bestPos) == true)
		levels.storePos(bestPos);

	if(queryBestKnown() == 0)
		return;

	if(pCrossMin != 0)
		pCrossMin->init(levels);
	else
		pCrossMinSimDraw->init(levels);

	Array<bool> *pLevelChanged = 0;
	if(transpose()) {
		pLevelChanged = new Array<bool>(-1,levels.size());
		(*pLevelChanged)[-1] = (*pLevelChanged)[levels.size()] = false;
	}

	int maxFails = fails();
	for( ; ; ) {

		int nFails = maxFails+1;
		do {

			// top-down traversal
			int nCrossingsNew = traverseTopDown(levels, pCrossMin, pCrossMinSimDraw, pLevelChanged);
			if(nCrossingsNew < nCrossingsOld) {
				if(nCrossingsNew < queryBestKnown() && postNewResult(nCrossingsNew, &bestPos) == true)
					levels.storePos(bestPos);

				nCrossingsOld = nCrossingsNew;
				nFails = maxFails+1;
			} else
				--nFails;

			// bottom-up traversal
			nCrossingsNew = traverseBottomUp(levels, pCrossMin, pCrossMinSimDraw, pLevelChanged);
			if(nCrossingsNew < nCrossingsOld) {
				if(nCrossingsNew < queryBestKnown() && postNewResult(nCrossingsNew, &bestPos) == true)
					levels.storePos(bestPos);

				nCrossingsOld = nCrossingsNew;
				nFails = maxFails+1;
			} else
				--nFails;

		} while(nFails > 0);

		if(getNextRun() == false)
			break;

#ifdef OGDF_HAVE_CPP11
		levels.permute(rng);
#else
		levels.permute();
#endif

		nCrossingsOld = (pCrossMin != 0) ? levels.calculateCrossings() : levels.calculateCrossingsSimDraw(subgraphs());
		if(nCrossingsOld < queryBestKnown() && postNewResult(nCrossingsOld, &bestPos) == true)
			levels.storePos(bestPos);
	}

	delete pLevelChanged;

	if(pCrossMin != 0)
		pCrossMin->cleanup();
	else
		pCrossMinSimDraw->cleanup();
}


// LayerByLayerSweep::CrossMinWorker

class LayerByLayerSweep::CrossMinWorker : public Thread {

	LayerByLayerSweep::CrossMinMaster &m_master;
	LayerByLayerSweep        *m_pCrossMin;
	TwoLayerCrossMinSimDraw *m_pCrossMinSimDraw;

	NodeArray<int>   m_bestPos;

public:
	CrossMinWorker(LayerByLayerSweep::CrossMinMaster &master, LayerByLayerSweep *pCrossMin, TwoLayerCrossMinSimDraw *pCrossMinSimDraw)
		: m_master(master), m_pCrossMin(pCrossMin), m_pCrossMinSimDraw(pCrossMinSimDraw)
	{
		OGDF_ASSERT( (pCrossMin != 0 && pCrossMinSimDraw == 0) || (pCrossMin == 0 && pCrossMinSimDraw != 0));
	}

	~CrossMinWorker() { delete m_pCrossMin; }

protected:
	virtual void doWork();
};

void LayerByLayerSweep::CrossMinWorker::doWork()
{
	HierarchyLevels levels(m_master.hierarchy());

#ifdef OGDF_HAVE_CPP11
	minstd_rand rng(m_master.rseed(threadID())); // different seeds per thread
	m_master.doWorkHelper(m_pCrossMin, m_pCrossMinSimDraw, levels, m_bestPos, true, rng);
#else
	m_master.doWorkHelper(m_pCrossMin, m_pCrossMinSimDraw, levels, m_bestPos, true);
#endif
}


//---------------------------------------------------------
// SugiyamaLayout
// Sugiyama drawing algorithm for hierarchical graphs
//---------------------------------------------------------

// SugiyamaLayout

SugiyamaLayout::SugiyamaLayout()
{
	m_ranking        .set(new LongestPathRanking);
	m_crossMin       .set(new BarycenterHeuristic);
	m_crossMinSimDraw.set(new SplitHeuristic);
	m_layout         .set(new FastHierarchyLayout);
	m_clusterLayout  .set(new OptimalHierarchyClusterLayout);
	m_packer         .set(new TileToRowsCCPacker);

	m_fails = 4;
	m_runs = 15;
	m_transpose = true;
	m_permuteFirst = false;

	m_arrangeCCs = true;
	m_minDistCC = LayoutStandards::defaultCCSeparation();
	m_pageRatio = 1.0;

#ifdef OGDF_MEMORY_POOL_NTS
	m_maxThreads = 1;
#else
	m_maxThreads = System::numberOfProcessors();
#endif

	m_alignBaseClasses = false;
	m_alignSiblings = false;

	m_subgraphs = 0;

	m_maxLevelSize = -1;
	m_numLevels = -1;
	m_timeReduceCrossings = 0.0;
}


void SugiyamaLayout::call(GraphAttributes &AG)
{
	doCall(AG,false);
}


void SugiyamaLayout::call(GraphAttributes &AG, NodeArray<int> &rank)
{
	doCall(AG,false,rank);
}


void SugiyamaLayout::doCall(GraphAttributes &AG, bool umlCall)
{
	NodeArray<int> rank;
	doCall(AG, umlCall, rank);
}


void SugiyamaLayout::doCall(GraphAttributes &AG, bool umlCall, NodeArray<int> &rank)
{
	const Graph &G = AG.constGraph();
	if (G.numberOfNodes() == 0)
		return;

	// compute connected component of G
	NodeArray<int> component(G);
	m_numCC = connectedComponents(G,component);

	const bool optimizeHorizEdges = (umlCall || rank.valid());
	if(!rank.valid())
	{
		if(umlCall)
		{
			LongestPathRanking ranking;
			ranking.alignBaseClasses(m_alignBaseClasses);
			ranking.alignSiblings(m_alignSiblings);

			ranking.callUML(AG,rank);

		} else {
			m_ranking.get().call(AG.constGraph(),rank);
		}
	}

	if(m_arrangeCCs) {
		// intialize the array of lists of nodes contained in a CC
		Array<List<node> > nodesInCC(m_numCC);

		node v;
		forall_nodes(v,G)
			nodesInCC[component[v]].pushBack(v);

		Hierarchy H;
		H.createEmpty(G);

		EdgeArray<edge> auxCopy(G);
		Array<DPoint> boundingBox(m_numCC);
		Array<DPoint> offset1(m_numCC);

		m_numLevels = m_maxLevelSize = 0;

		int totalCrossings = 0;
		for(int i = 0; i < m_numCC; ++i)
		{
			// adjust ranks in cc to start with 0
			int minRank = numeric_limits<int>::max();
			ListConstIterator<node> it;
			for(it = nodesInCC[i].begin(); it.valid(); ++it)
				if(rank[*it] < minRank)
					minRank = rank[*it];

			if(minRank != 0) {
				for(it = nodesInCC[i].begin(); it.valid(); ++it)
					rank[*it] -= minRank;
			}
			H.createEmpty(G);
			H.initByNodes(nodesInCC[i],auxCopy,rank);
			//HierarchyLevels levels(H);
			//reduceCrossings(levels);
			const HierarchyLevelsBase *pLevels = reduceCrossings(H);
			const HierarchyLevelsBase &levels = *pLevels;
			totalCrossings += m_nCrossings;

			const GraphCopy &GC = H;
			NodeArray<bool> mark(GC);

			m_layout.get().call(levels,AG);

			double
				minX =  numeric_limits<double>::max(),
				maxX = -numeric_limits<double>::max(),
				minY =  numeric_limits<double>::max(),
				maxY = -numeric_limits<double>::max();

			node vCopy;
			forall_nodes(vCopy,GC)
			{
				mark[vCopy] = false;
				node v = GC.original(vCopy);
				if(v == 0) continue;

				if(AG.x(v)-AG.width (v)/2 < minX) minX = AG.x(v)-AG.width(v) /2;
				if(AG.x(v)+AG.width (v)/2 > maxX) maxX = AG.x(v)+AG.width(v) /2;
				if(AG.y(v)-AG.height(v)/2 < minY) minY = AG.y(v)-AG.height(v)/2;
				if(AG.y(v)+AG.height(v)/2 > maxY) maxY = AG.y(v)+AG.height(v)/2;
			}

			if(optimizeHorizEdges)
			{
				for(int i = 0; i < levels.size(); ++i) {
					const LevelBase &L = levels[i];
					for(int j = 0; j < L.size(); ++j) {
						node v = L[j];
						if(!GC.isDummy(v)) continue;
						edge e = GC.original(v->firstAdj()->theEdge());
						if(e == 0) continue;
						node src = GC.copy(e->source());
						node tgt = GC.copy(e->target());

						if(H.rank(src) == H.rank(tgt)) {
							int minPos = levels.pos(src), maxPos = levels.pos(tgt);
							if(minPos > maxPos) std::swap(minPos,maxPos);

							bool straight = true;
							const LevelBase &L_e = levels[H.rank(src)];
							for(int p = minPos+1; p < maxPos; ++p) {
								if(!H.isLongEdgeDummy(L_e[p]) && mark[L_e[p]] == false) {
									straight = false; break;
								}
							}
							if(straight) {
								AG.bends(e).clear();
								mark[v] = true;
							}
						}
					}
				}
			}

			edge eCopy;
			forall_edges(eCopy,GC)
			{
				edge e = GC.original(eCopy);
				if(e == 0 || eCopy != GC.chain(e).front()) continue;

				const DPolyline &dpl = AG.bends(e);
				ListConstIterator<DPoint> it;
				for(it = dpl.begin(); it.valid(); ++it)
				{
					if((*it).m_x < minX) minX = (*it).m_x;
					if((*it).m_x > maxX) maxX = (*it).m_x;
					if((*it).m_y < minY) minY = (*it).m_y;
					if((*it).m_y > maxY) maxY = (*it).m_y;
				}
			}

			minX -= m_minDistCC;
			minY -= m_minDistCC;

			boundingBox[i] = DPoint(maxX - minX, maxY - minY);
			offset1    [i] = DPoint(minX,minY);

			m_numLevels = max(m_numLevels, levels.size());
			for(int i = 0; i <= levels.high(); i++) {
				const LevelBase &l = levels[i];
				m_maxLevelSize = max(m_maxLevelSize, l.size());
			}
			delete pLevels;
		}

		m_nCrossings = totalCrossings;

		// call packer
		Array<DPoint> offset(m_numCC);
		m_packer.get().call(boundingBox,offset,m_pageRatio);

		// The arrangement is given by offset to the origin of the coordinate
		// system. We still have to shift each node and edge by the offset
		// of its connected component.

		for(int i = 0; i < m_numCC; ++i)
		{
			const List<node> &nodes = nodesInCC[i];

			const double dx = offset[i].m_x - offset1[i].m_x;
			const double dy = offset[i].m_y - offset1[i].m_y;

			// iterate over all nodes in ith CC
			ListConstIterator<node> it;
			for(it = nodes.begin(); it.valid(); ++it)
			{
				node v = *it;

				AG.x(v) += dx;
				AG.y(v) += dy;

				edge e;
				forall_adj_edges(e,v)
				{
					if(e->isSelfLoop() || e->source() != v) continue;

					DPolyline &dpl = AG.bends(e);
					ListIterator<DPoint> it;
					for(it = dpl.begin(); it.valid(); ++it)
					{
						(*it).m_x += dx;
						(*it).m_y += dy;
					}
				}
			}
		}

	} else {
		int minRank = numeric_limits<int>::max();
		node v;
		forall_nodes(v,G)
			if(rank[v] < minRank)
				minRank = rank[v];

		if(minRank != 0) {
			forall_nodes(v,G)
				rank[v] -= minRank;
		}

		Hierarchy H(G,rank);

		{  // GC's scope is limited to allow reassignment after crossing reduction phase
			const GraphCopy &GC = H;

			m_compGC.init(GC);
			forall_nodes(v,GC) {
				node vOrig = GC.original(v);
				if(vOrig == 0)
					vOrig = GC.original(v->firstAdj()->theEdge())->source();

				m_compGC[v] = component[vOrig];
			}
		}

		const HierarchyLevelsBase *pLevels = reduceCrossings(H);
		const HierarchyLevelsBase &levels = *pLevels;
		//HierarchyLevels levels(H);
		//reduceCrossings(levels);
		m_compGC.init();

		const GraphCopy &GC = H;

		m_layout.get().call(levels,AG);

		if(optimizeHorizEdges)
		{
			NodeArray<bool> mark(GC,false);
			for(int i = 0; i < levels.size(); ++i) {
				const LevelBase &L = levels[i];
				for(int j = 0; j < L.size(); ++j) {
					node v = L[j];
					if(!GC.isDummy(v)) continue;
					edge e = GC.original(v->firstAdj()->theEdge());
					if(e == 0) continue;
					node src = GC.copy(e->source());
					node tgt = GC.copy(e->target());

					if(H.rank(src) == H.rank(tgt)) {
						int minPos = levels.pos(src), maxPos = levels.pos(tgt);
						if(minPos > maxPos) std::swap(minPos,maxPos);

						bool straight = true;
						const LevelBase &L_e = levels[H.rank(src)];
						for(int p = minPos+1; p < maxPos; ++p) {
							if(!H.isLongEdgeDummy(L_e[p]) && mark[L_e[p]] == false) {
								straight = false; break;
							}
						}
						if(straight) {
							AG.bends(e).clear();
							mark[v] = true;
						}
					}
				}
			}
		}

		m_numLevels = levels.size();
		m_maxLevelSize = 0;
		for(int i = 0; i <= levels.high(); i++) {
			const LevelBase &l = levels[i];
			if (l.size() > m_maxLevelSize)
				m_maxLevelSize = l.size();
		}
		delete pLevels;
	}
}


void SugiyamaLayout::callUML(GraphAttributes &AG)
{
	doCall(AG,true);
}

/*
void SugiyamaLayout::reduceCrossings(HierarchyLevels &levels)
{
	OGDF_ASSERT(m_runs >= 1);


	__int64 t;
	System::usedRealTime(t);

	LayerByLayerSweep          *pCrossMin = 0;
	TwoLayerCrossMinSimDraw   *pCrossMinSimDraw = 0;

	if(useSubgraphs() == false)
		pCrossMin = &m_crossMin.get();
	else
		pCrossMinSimDraw = &m_crossMinSimDraw.get();

	int nThreads = min(m_maxThreads, m_runs);

	int seed = rand();
#ifdef OGDF_HAVE_CPP11
	minstd_rand rng(seed);
#endif

	LayerByLayerSweep::CrossMinMaster master(*this, levels.hierarchy(), seed, m_runs - nThreads);

	Array<LayerByLayerSweep::CrossMinWorker *> thread(nThreads-1);
	for(int i = 0; i < nThreads-1; ++i) {
		thread[i] = new LayerByLayerSweep::CrossMinWorker(master,
			(pCrossMin        != 0) ? pCrossMin       ->clone() : 0,
			(pCrossMinSimDraw != 0) ? pCrossMinSimDraw->clone() : 0);
		thread[i]->start();
	}

	NodeArray<int> bestPos;
	master.doWorkHelper(pCrossMin, pCrossMinSimDraw, levels, bestPos, m_permuteFirst
#ifdef OGDF_HAVE_CPP11
		, rng
#endif
		);

	for(int i = 0; i < nThreads-1; ++i)
		thread[i]->join();

	master.restore(levels, m_nCrossings);

	for(int i = 0; i < nThreads-1; ++i)
		delete thread[i];

	t = System::usedRealTime(t);
	m_timeReduceCrossings = double(t) / 1000;
}
//*/
const HierarchyLevels *LayerByLayerSweep::reduceCrossings(const SugiyamaLayout &sugi, const Hierarchy &H)
{
	HierarchyLevels *levels = new HierarchyLevels(H);
	OGDF_ASSERT(sugi.runs() >= 1);

	int nThreads = min( sugi.maxThreads(), sugi.runs() );

	int seed = rand();
#ifdef OGDF_HAVE_CPP11
	minstd_rand rng(seed);
#endif

	LayerByLayerSweep::CrossMinMaster master(sugi, levels->hierarchy(), seed, sugi.runs() - nThreads);

	Array<LayerByLayerSweep::CrossMinWorker *> thread(nThreads-1);
	for ( int i = 0; i < nThreads - 1; ++i ) {
		thread[i] = new LayerByLayerSweep::CrossMinWorker(master,
			clone(), 0);
		thread[i] -> start();
	}

	NodeArray<int> bestPos;
	master.doWorkHelper(this, 0, *levels, bestPos, sugi.permuteFirst()
#ifdef OGDF_HAVE_CPP11
		, rng
#endif
		);

	for (int i = 0; i < nThreads - 1; ++i )
		thread[i] -> join();

	// ??
	int x = 0;
	master.restore(*levels, x );

	for ( int i = 0; i < nThreads - 1; ++i )
		delete thread[i];

	return levels;
}

const HierarchyLevelsBase *SugiyamaLayout::reduceCrossings(Hierarchy &H)
{
	OGDF_ASSERT(m_runs >= 1);


	if (useSubgraphs() == false) {
		__int64 t;
		System::usedRealTime(t);
		const HierarchyLevelsBase *levels = m_crossMin.get().reduceCrossings(*this,H);
		t = System::usedRealTime(t);
		m_timeReduceCrossings = double(t) / 1000;
		m_nCrossings = levels -> calculateCrossings();
		return levels;
	}



	//unchanged crossing reduction of subgraphs
	HierarchyLevels *pLevels = new HierarchyLevels(H);
	HierarchyLevels levels = *pLevels;

	__int64 t;
	System::usedRealTime(t);

	LayerByLayerSweep          *pCrossMin = 0;
	TwoLayerCrossMinSimDraw   *pCrossMinSimDraw = 0;

	//if(useSubgraphs() == false)
	//	pCrossMin = &m_crossMin.get();
	//else
		pCrossMinSimDraw = &m_crossMinSimDraw.get();

	int nThreads = min(m_maxThreads, m_runs);

	int seed = rand();
#ifdef OGDF_HAVE_CPP11
	minstd_rand rng(seed);
#endif

	LayerByLayerSweep::CrossMinMaster master(*this, levels.hierarchy(), seed, m_runs - nThreads);

	Array<LayerByLayerSweep::CrossMinWorker *> thread(nThreads-1);
	for(int i = 0; i < nThreads-1; ++i) {
		thread[i] = new LayerByLayerSweep::CrossMinWorker(master,
			(pCrossMin        != 0) ? pCrossMin       ->clone() : 0,
			(pCrossMinSimDraw != 0) ? pCrossMinSimDraw->clone() : 0);
		thread[i]->start();
	}

	NodeArray<int> bestPos;
	master.doWorkHelper(pCrossMin, pCrossMinSimDraw, levels, bestPos, m_permuteFirst
#ifdef OGDF_HAVE_CPP11
		, rng
#endif
		);

	for(int i = 0; i < nThreads-1; ++i)
		thread[i]->join();

	master.restore(levels, m_nCrossings);

	for(int i = 0; i < nThreads-1; ++i)
		delete thread[i];

	t = System::usedRealTime(t);
	m_timeReduceCrossings = double(t) / 1000;

	return pLevels;
}

} // end namespace ogdf
