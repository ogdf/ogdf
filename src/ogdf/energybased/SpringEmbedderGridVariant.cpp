/** \file
 * \brief Implementation of Spring-Embedder algorithm
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

#include <ogdf/internal/energybased/SEGV_ForceModel.h>

#include <ogdf/packing/TileToRowsCCPacker.h>
#include <ogdf/basic/GraphCopyAttributes.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/Barrier.h>
#include <ogdf/basic/Thread.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/Math.h>
#include <random>


using std::minstd_rand;
using std::uniform_real_distribution;


namespace ogdf {

template<typename T>
inline void updateMin(T &var, T val)
{
	if(val < var) var = val;
}

template<typename T>
inline void updateMax(T &var, T val)
{
	if(val > var) var = val;
}


SpringEmbedderGridVariant::SpringEmbedderGridVariant()
{
	// default parameters
	m_iterations        = 400;
	m_iterationsImprove = 200;
	m_coolDownFactor    = 0.999;
	m_forceLimitStep    = 0.5;

	m_xleft = m_ysmall = 0.0;
	m_xright = m_ybig = 250.0;
	m_noise = true;

	m_forceModel        = SpringForceModel::FruchtermanReingold;
	m_forceModelImprove = SpringForceModel::FruchtermanReingoldModRep;

	m_avgConvergenceFactor = 0.1;
	m_maxConvergenceFactor = 0.2;

	m_scaling = Scaling::scaleFunction;
	m_scaleFactor = 4.0;
	m_bbXmin = 0.0;
	m_bbXmax = 100.0;
	m_bbYmin = 0.0;
	m_bbYmax = 100.0;

	double def_nw = LayoutStandards::defaultNodeWidth();
	double def_nh = LayoutStandards::defaultNodeHeight();
	m_idealEdgeLength = LayoutStandards::defaultNodeSeparation() + sqrt( def_nw*def_nw + def_nh*def_nh);
	m_minDistCC = LayoutStandards::defaultCCSeparation();
	m_pageRatio = 1.0;

	m_maxThreads = System::numberOfProcessors();
}


class SpringEmbedderGridVariant::Master {

	const SpringEmbedderGridVariant &m_spring;
	const GraphCopy                 &m_gc;
	GraphAttributes                 &m_ga;
	DPoint                          &m_boundingBox;

	NodeArray<int>         m_index;
	Array<NodeInfo>        m_vInfo;
	Array<DPoint>          m_disp;
	Array<int>             m_adjLists;
	Array2D<ListPure<int>> m_gridCell;

	ForceModelBase *m_forceModel;
	ForceModelBase *m_forceModelImprove;

	Array<Worker*>  m_worker;
	Barrier        *m_barrier;

	double m_idealEdgeLength;

	double m_tNull;
	double m_cF;
	double m_t;
	double m_coolingFactor;

	double m_k2;

	double m_xleft;
	double m_xright;
	double m_ysmall;
	double m_ybig;

	double m_avgDisplacement;
	double m_maxDisplacement;
	double m_scaleFactor;

public:
	Master(const SpringEmbedderGridVariant &spring, const GraphCopy &gc, GraphAttributes &ga, DPoint &boundingBox);
	~Master() {
		delete m_barrier;
		delete m_forceModel;
		delete m_forceModelImprove;
	}

	int numberOfNodes() const { return m_vInfo.size(); }
	int numberOfIterations() const { return m_spring.m_iterations; }
	int numberOfIterationsImprove() const { return m_spring.m_iterationsImprove; }

	double xleft () const { return m_xleft;  }
	double ysmall() const { return m_ysmall; }

	void initUnfoldPhase();
	void initImprovementPhase();
	void coolDown();
	double maxForceLength() const { return m_t; }
	double coolingFactor() const { return m_coolingFactor; }

	double idealEdgeLength() const { return m_idealEdgeLength; }
	double boxLength() const { return m_k2; }
	bool noise() const { return m_spring.m_noise; }

	const GraphCopy &getGraph() const { return m_gc; }
	GraphAttributes &getAttributes() { return m_ga; }

	const NodeArray<int> &index() const { return m_index; }
	Array<NodeInfo> &vInfo() { return m_vInfo; }
	Array<DPoint> &disp() { return m_disp; }
	Array<int> &adjLists() { return m_adjLists; }
	Array2D<ListPure<int>> &gridCell() { return m_gridCell; }

	const ForceModelBase &forceModel() const { return *m_forceModel; }
	const ForceModelBase &forceModelImprove() const { return *m_forceModelImprove; }

	void syncThreads() {
		if(m_barrier)
			m_barrier->threadSync();
	}

	void computeGrid(double wsum, double hsum, double xmin, double xmax, double ymin, double ymax);
	void updateGridAndMoveNodes();
	void scaleLayout(double sumLengths);
	double scaleFactor() const { return m_scaleFactor; }
	void computeFinalBB();

	bool hasConverged() const {
		return m_avgDisplacement <= m_spring.m_avgConvergenceFactor * m_idealEdgeLength
			&& m_maxDisplacement <= m_spring.m_maxConvergenceFactor * m_idealEdgeLength;
	}

	double avgDisplacement() const { return m_avgDisplacement; }
	double maxDisplacement() const { return m_maxDisplacement; }
};


class SpringEmbedderGridVariant::Worker {

	friend class Master;

	unsigned int m_id;
	Master  &m_master;

	int  m_vStartIndex;
	int  m_vStopIndex;
	node m_vStart;
	node m_vStop;
	int  m_eStartIndex;

	double m_wsum;
	double m_hsum;
	double m_xmin;
	double m_xmax;
	double m_ymin;
	double m_ymax;

	double m_sumForces;
	double m_maxForce;
	double m_sumLengths;

public:
	Worker(unsigned int id, Master &master, int vStartIndex, int vStopIndex, node vStart, node vStop, int eStartIndex)
		: m_id(id), m_master(master), m_vStartIndex(vStartIndex), m_vStopIndex(vStopIndex), m_vStart(vStart), m_vStop(vStop), m_eStartIndex(eStartIndex) { }

	void operator()();

private:
	double sumUpLengths(Array<NodeInfo> &vInfo, const Array<int> &adjLists);
	void scaling(Array<NodeInfo> &vInfo, const Array<int> &adjLists);
	void finalScaling(Array<NodeInfo> &vInfo, const Array<int> &adjLists);
};



SpringEmbedderGridVariant::Master::Master(const SpringEmbedderGridVariant &spring, const GraphCopy &gc, GraphAttributes &ga, DPoint &boundingBox) :
	m_spring(spring), m_gc(gc), m_ga(ga), m_boundingBox(boundingBox),
	m_index(gc), m_vInfo(gc.numberOfNodes()), m_disp(gc.numberOfNodes()), m_adjLists(2*gc.numberOfEdges()), m_forceModel(nullptr), m_forceModelImprove(nullptr),
	m_avgDisplacement(numeric_limits<double>::max()), m_maxDisplacement(numeric_limits<double>::max())
{
	const unsigned int minNodesPerThread = 64;
	const unsigned int n = gc.numberOfNodes();

	unsigned int nThreads = max( 1u, min(spring.m_maxThreads, (n/4) / (minNodesPerThread/4)) );
	m_worker.init(nThreads);

	if(nThreads == 1) {

		m_barrier = nullptr;

		int nextIndex = 0;
		for(node v : gc.nodes)
			m_index[v] = nextIndex++;

		m_worker[0] = new Worker(0, *this, 0, n, gc.firstNode(), nullptr, 0);
		(*m_worker[0])();

	} else {

		unsigned int nodesPerThread = 4*((n/4) / nThreads);

		Array<node> startNode(nThreads+1);
		Array<int>  startIndex(nThreads+1);
		Array<int>  eStartIndex(nThreads+1);

		int nextIndex = 0, j = 0, t = 0;
		for(node v : gc.nodes) {
			if(nextIndex % nodesPerThread == 0) {
				startNode [t] = v;
				startIndex[t] = nextIndex;
				eStartIndex [t] = j;
				++t;
			}
			m_index[v] = nextIndex++;
			j += v->degree();
		}

		startNode [nThreads] = nullptr;
		startIndex[nThreads] = gc.numberOfNodes();

		m_barrier = new Barrier(nThreads);

		Array<Thread> thread(nThreads-1);

		for(unsigned int i = 1; i < nThreads; ++i) {
			m_worker[i] =
				new Worker(i, *this, startIndex[i], startIndex[i+1], startNode[i], startNode[i+1], eStartIndex[i]);
			thread [i-1] = Thread(*m_worker[i]);
		}

		m_worker[0] = new Worker(0, *this, startIndex[0], startIndex[1], startNode[0], startNode[1], eStartIndex[0]);
		(*m_worker[0])();

		for(unsigned int i = 1; i < nThreads; ++i) {
			thread[i-1].join();
			delete m_worker[i];
		}
	}

	delete m_worker[0];
}


void SpringEmbedderGridVariant::Master::computeGrid(double wsum, double hsum, double xmin, double xmax, double ymin, double ymax)
{
	const int n = m_gc.numberOfNodes();

	for(int t = 1; t <= m_worker.high(); ++t) {
		updateMin(xmin, m_worker[t]->m_xmin);
		updateMax(xmax, m_worker[t]->m_xmax);
		updateMin(ymin, m_worker[t]->m_ymin);
		updateMax(ymax, m_worker[t]->m_ymax);
		wsum += m_worker[t]->m_wsum;
		hsum += m_worker[t]->m_hsum;
	}

	Scaling scaling = m_spring.m_scaling;

	// handle special case of zero area bounding box
	if(xmin == xmax || ymin == ymax) {
		if(scaling == Scaling::userBoundingBox) {
			m_xleft  = m_spring.m_bbXmin;
			m_xright = m_spring.m_bbXmax;
			m_ysmall = m_spring.m_bbYmin;
			m_ybig   = m_spring.m_bbYmax;

		} else {
			m_idealEdgeLength = max(1e-3, m_spring.m_idealEdgeLength);
			m_xleft  = m_ysmall = 0;
			m_xright = m_ybig   = m_idealEdgeLength * sqrt(double(n));
		}

		minstd_rand rng(randomSeed());
		uniform_real_distribution<> rand_x(m_xleft ,m_xright);
		uniform_real_distribution<> rand_y(m_ysmall,m_ybig  );

		for(int j = 0; j < n; ++j) {
			m_vInfo[j].m_pos.m_x = rand_x(rng);
			m_vInfo[j].m_pos.m_y = rand_y(rng);
		}

	} else {
		double scaleFactor = m_spring.m_scaleFactor;

		switch(scaling) {
		case Scaling::input:
			m_xleft  = xmin;
			m_xright = xmax;
			m_ysmall = ymin;
			m_ybig   = ymax;
			break;

		case Scaling::userBoundingBox:
		case Scaling::scaleFunction:
		case Scaling::useIdealEdgeLength:

			if(scaling == Scaling::userBoundingBox) {
				m_xleft  = m_spring.m_bbXmin;
				m_xright = m_spring.m_bbXmax;
				m_ysmall = m_spring.m_bbYmin;
				m_ybig   = m_spring.m_bbYmax;

			} else if(scaling == Scaling::scaleFunction) {
				double sqrt_n = sqrt((double)n);
				m_xleft  = m_ysmall = 0;
				m_xright = (wsum > 0) ? scaleFactor * wsum / sqrt_n : 1;
				m_ybig   = (hsum > 0) ? scaleFactor * hsum / sqrt_n : 1;
			} else {
				m_idealEdgeLength = max(1e-3, m_spring.m_idealEdgeLength);
				double w = xmax - xmin;
				double h = ymax - ymin;
				double r = (w > 0) ? h / w : 1.0;
				m_xleft = m_ysmall = 0;
				m_xright = m_idealEdgeLength * sqrt(double(n) / r);
				m_ybig   = r * m_xright;
			}
			// Compute scaling such that layout coordinates fit into used bounding box
			double fx = (xmax == xmin) ? 1.0 : m_xright / (xmax - xmin);
			double fy = (ymax == ymin) ? 1.0 : m_ybig   / (ymax - ymin);

			// Adjust coordinates accordingly
			for(int j = 0; j < n; ++j) {
				m_vInfo[j].m_pos.m_x = m_xleft  + (m_vInfo[j].m_pos.m_x - xmin) * fx;
				m_vInfo[j].m_pos.m_y = m_ysmall + (m_vInfo[j].m_pos.m_y - ymin) * fy;
			}
		}
	}

	double width  = m_xright - m_xleft;
	double height = m_ybig - m_ysmall;

	OGDF_ASSERT(width >= 0);
	OGDF_ASSERT(height >= 0);

	initUnfoldPhase();

	double k = sqrt(width*height / n);
	m_k2 = max(1e-3, 2*k);

	if(scaling != Scaling::useIdealEdgeLength)
		m_idealEdgeLength = k;

	switch(m_spring.m_forceModel) {
	case SpringForceModel::FruchtermanReingold:
		m_forceModel = new ForceModelFR(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::FruchtermanReingoldModAttr:
		m_forceModel = new ForceModelFRModAttr(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::FruchtermanReingoldModRep:
		m_forceModel = new ForceModelFRModRep(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::Eades:
		m_forceModel = new ForceModelEades(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::Hachul:
		m_forceModel = new ForceModelHachul(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::Gronemann:
		m_forceModel = new ForceModelGronemann(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	}

	switch(m_spring.m_forceModelImprove) {
	case SpringForceModel::FruchtermanReingold:
		m_forceModelImprove = new ForceModelFR(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::FruchtermanReingoldModAttr:
		m_forceModelImprove = new ForceModelFRModAttr(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::FruchtermanReingoldModRep:
		m_forceModelImprove = new ForceModelFRModRep(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::Eades:
		m_forceModelImprove = new ForceModelEades(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::Hachul:
		m_forceModelImprove = new ForceModelHachul(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	case SpringForceModel::Gronemann:
		m_forceModelImprove = new ForceModelGronemann(m_vInfo, m_adjLists, m_gridCell, m_idealEdgeLength);
		break;
	}

	// build  grid cells
	int xA = int(width / m_k2 + 2);
	int yA = int(height / m_k2 + 2);
	m_gridCell.init(-1,xA,-1,yA);

	for(int j = 0; j < n; ++j)
	{
		NodeInfo &vj = m_vInfo[j];

		vj.m_gridX = int((vj.m_pos.m_x - m_xleft ) / m_k2);
		vj.m_gridY = int((vj.m_pos.m_y - m_ysmall) / m_k2);

		OGDF_ASSERT(vj.m_gridX < xA);
		OGDF_ASSERT(vj.m_gridX > -1);
		OGDF_ASSERT(vj.m_gridY < yA);
		OGDF_ASSERT(vj.m_gridY > -1);

		m_vInfo[j].m_lit = m_gridCell(vj.m_gridX,vj.m_gridY).pushFront(j);
	}
}


void SpringEmbedderGridVariant::Master::updateGridAndMoveNodes()
{
	const Worker &worker = *m_worker[0];

	double xmin = worker.m_xmin;
	double xmax = worker.m_xmax;
	double ymin = worker.m_ymin;
	double ymax = worker.m_ymax;
	double maxForce  = worker.m_maxForce;
	double sumForces = worker.m_sumForces;

	for(int t = 1; t <= m_worker.high(); ++t) {
		updateMin(xmin, m_worker[t]->m_xmin);
		updateMax(xmax, m_worker[t]->m_xmax);
		updateMin(ymin, m_worker[t]->m_ymin);
		updateMax(ymax, m_worker[t]->m_ymax);
		updateMax(maxForce, m_worker[t]->m_maxForce);
		sumForces += m_worker[t]->m_sumForces;
	}

	m_avgDisplacement = sumForces / numberOfNodes();

	const int xA = m_gridCell.high1();
	const int yA = m_gridCell.high2();

	// prevent drawing area from getting too small
	double hMargin = 0.5 * max(0.0, m_idealEdgeLength * xA - (xmax-xmin));
	double vMargin = 0.5 * max(0.0, m_idealEdgeLength * yA - (ymax-ymin));

	// set new values
	m_xleft  = xmin - hMargin;
	m_xright = xmax + hMargin;
	m_ysmall = ymin - vMargin;
	m_ybig   = ymax + vMargin;

	m_k2 = max( (m_xright-m_xleft) / (xA-1), (m_ybig-m_ysmall) / (yA-1) );

	// move nodes
	for(int j = 0; j <= m_vInfo.high(); ++j) {
		NodeInfo &vj = m_vInfo[j];

		// new position
		vj.m_pos += m_disp[j];

		// new cell
		int grid_x = int((vj.m_pos.m_x - m_xleft ) / m_k2);
		int grid_y = int((vj.m_pos.m_y - m_ysmall) / m_k2);

		OGDF_ASSERT(grid_x >= 0);
		OGDF_ASSERT(grid_x < m_gridCell.high1());
		OGDF_ASSERT(grid_y >= 0);
		OGDF_ASSERT(grid_y < m_gridCell.high2());

		// move to different cell?
		if( (grid_x != vj.m_gridX) || (grid_y != vj.m_gridY) ) {
			m_gridCell(vj.m_gridX,vj.m_gridY).moveToFront(vj.m_lit, m_gridCell(grid_x,grid_y));
			vj.m_gridX = grid_x;
			vj.m_gridY = grid_y;
		}
	}
}


void SpringEmbedderGridVariant::Master::computeFinalBB()
{
	double xmin = m_worker[0]->m_xmin;
	double xmax = m_worker[0]->m_xmax;
	double ymin = m_worker[0]->m_ymin;
	double ymax = m_worker[0]->m_ymax;

	for(int t = 1; t <= m_worker.high(); ++t) {
		updateMin(xmin, m_worker[t]->m_xmin);
		updateMax(xmax, m_worker[t]->m_xmax);
		updateMin(ymin, m_worker[t]->m_ymin);
		updateMax(ymax, m_worker[t]->m_ymax);
	}

	xmin -= m_spring.m_minDistCC;
	ymin -= m_spring.m_minDistCC;
	m_boundingBox = DPoint(xmax - xmin, ymax - ymin);

	m_xleft  = xmin;
	m_ysmall = ymin;
}


void SpringEmbedderGridVariant::Master::coolDown()
{
	m_cF += m_spring.forceLimitStep();
	m_t = m_tNull / Math::log2(m_cF);

	m_coolingFactor *= m_spring.coolDownFactor();
}


void SpringEmbedderGridVariant::Master::initUnfoldPhase()
{
	// cool down
	m_t = m_tNull = 0.25 * m_idealEdgeLength * sqrt(numberOfNodes());
	m_cF = 2.0;
	m_coolingFactor = m_spring.coolDownFactor();

	// convergence
	m_avgDisplacement = numeric_limits<double>::max();
	m_maxDisplacement = numeric_limits<double>::max();
}


void SpringEmbedderGridVariant::Master::initImprovementPhase()
{
	// cool down
	m_t  = m_tNull;
	m_cF = 2.0;
	m_coolingFactor = m_spring.coolDownFactor();

	// convergence
	m_avgDisplacement = numeric_limits<double>::max();
	m_maxDisplacement = numeric_limits<double>::max();
}


void SpringEmbedderGridVariant::Master::scaleLayout(double sumLengths)
{
	for(int t = 1; t <= m_worker.high(); ++t)
		sumLengths += m_worker[t]->m_sumLengths;

	const int m = m_gc.numberOfEdges();
	m_scaleFactor = m_idealEdgeLength / sumLengths * m;

	// set new values
	m_xleft  *= m_scaleFactor;
	m_xright *= m_scaleFactor;
	m_ysmall *= m_scaleFactor;
	m_ybig   *= m_scaleFactor;

	m_k2 = max( (m_xright-m_xleft) / (m_gridCell.high1()-1), (m_ybig-m_ysmall) / (m_gridCell.high2()-1) );
}


void SpringEmbedderGridVariant::call(GraphAttributes &AG)
{
	//static int i = 0;

	//string inputfilename = "input_" + to_string(i) + ".gml";
	//string outputfilename = "output_" + to_string(i++) + ".gml";

	//GraphIO::writeGML(AG, inputfilename);


	const Graph &G = AG.constGraph();
	if(G.empty())
		return;

	// all edges straight-line
	AG.clearAllBends();

	GraphCopy GC;
	GC.createEmpty(G);

	// compute connected component of G
	NodeArray<int> component(G);
	int numCC = connectedComponents(G,component);

	// intialize the array of lists of nodes contained in a CC
	Array<List<node> > nodesInCC(numCC);

	for(node v : G.nodes)
		nodesInCC[component[v]].pushBack(v);

	EdgeArray<edge> auxCopy(G);
	Array<DPoint> boundingBox(numCC);

	for(int i = 0; i < numCC; ++i)
	{
		GC.initByNodes(nodesInCC[i],auxCopy);
		makeSimpleUndirected(GC);

		const int n = GC.numberOfNodes();

		// special case: just one node
		if(n == 1) {
			node vOrig = GC.original(GC.firstNode());
			AG.x(vOrig) = AG.y(vOrig) = 0;
			boundingBox[i] = DPoint(0,0);
			continue;
		}

		Master master(*this, GC, AG, boundingBox[i]);
	}

	Array<DPoint> offset(numCC);
	TileToRowsCCPacker packer;
	packer.call(boundingBox,offset,m_pageRatio);

	// The arrangement is given by offset to the origin of the coordinate
	// system. We still have to shift each node and edge by the offset
	// of its connected component.

	for(int i = 0; i < numCC; ++i)
	{
		const List<node> &nodes = nodesInCC[i];

		const double dx = offset[i].m_x;
		const double dy = offset[i].m_y;

		// iterate over all nodes in ith CC
		ListConstIterator<node> it;
		for(node v : nodes) {
			AG.x(v) += dx;
			AG.y(v) += dy;
		}
	}

	//GraphIO::writeGML(AG, outputfilename);
}


double SpringEmbedderGridVariant::Worker::sumUpLengths(Array<NodeInfo> &vInfo, const Array<int> &adjLists)
{
	double sumLengths = 0.0;
	for(int j = m_vStartIndex; j < m_vStopIndex; ++j)
	{
		const NodeInfo &vj = vInfo[j];
		for(int l = vj.m_adjBegin; l != vj.m_adjStop; ++l) {
			int u = adjLists[l];
			if(u < j) {
				DPoint dist = vj.m_pos - vInfo[u].m_pos;
				sumLengths += dist.norm();
			}
		}
	}

	return sumLengths;
}


void SpringEmbedderGridVariant::Worker::scaling(Array<NodeInfo> &vInfo, const Array<int> &adjLists)
{
	m_sumLengths = sumUpLengths(vInfo, adjLists);

	m_master.syncThreads();

	if(m_id == 0)
		m_master.scaleLayout(m_sumLengths);

	m_master.syncThreads();

	double s = m_master.scaleFactor();
	for(int j = m_vStartIndex; j < m_vStopIndex; ++j)
		vInfo[j].m_pos *= s;

	if(m_id == 0)
		m_master.initImprovementPhase();

	m_master.syncThreads();
}


void SpringEmbedderGridVariant::Worker::finalScaling(Array<NodeInfo> &vInfo, const Array<int> &adjLists)
{
	m_sumLengths = sumUpLengths(vInfo, adjLists);

	m_master.syncThreads();

	if(m_id == 0)
		m_master.scaleLayout(m_sumLengths);

	m_master.syncThreads();

	double s = m_master.scaleFactor();

	const GraphCopy &gc = m_master.getGraph();
	GraphAttributes &ga = m_master.getAttributes();

	double xmin = numeric_limits<double>::max(), xmax = -numeric_limits<double>::max();
	double ymin = numeric_limits<double>::max(), ymax = -numeric_limits<double>::max();

	node v = m_vStart;
	for(int j = m_vStartIndex; j < m_vStopIndex; ++j) {
		node vOrig = gc.original(v);
		NodeInfo &vj = vInfo[j];

		double xv = s * vj.m_pos.m_x;
		double yv = s * vj.m_pos.m_y;
		vj.m_pos.m_x = xv;
		vj.m_pos.m_y = yv;

		double wv = ga.width(vOrig);
		double hv = ga.height(vOrig);

		updateMin(xmin, xv-0.5*wv);
		updateMax(xmax, xv+0.5*wv);
		updateMin(ymin, yv-0.5*hv);
		updateMax(ymax, yv+0.5*hv);
	}

	m_xmin = xmin; m_xmax = xmax;
	m_ymin = ymin; m_ymax = ymax;

	m_master.syncThreads();
}


void SpringEmbedderGridVariant::Worker::operator()()
{
	const double forceScaleFactor = 0.1; //0.05;

	const GraphCopy &gc = m_master.getGraph();
	GraphAttributes &ga = m_master.getAttributes();

	const NodeArray<int>  &index    = m_master.index();
	Array<NodeInfo>       &vInfo    = m_master.vInfo();
	Array<int>            &adjLists = m_master.adjLists();

	//---------------------------------
	// Initialize
	//---------------------------------

	double wsum = 0, hsum = 0;
	double xmin = numeric_limits<double>::max(), xmax = -numeric_limits<double>::max();
	double ymin = numeric_limits<double>::max(), ymax = -numeric_limits<double>::max();

	int adjCounter = m_eStartIndex;
	int j = m_vStartIndex;
	for(node v = m_vStart; v != m_vStop; v = v->succ()) {
		node vOrig = gc.original(v);

		double x = ga.x(vOrig), y = ga.y(vOrig);
		wsum += ga.width(vOrig);
		hsum += ga.height(vOrig);

		vInfo[j].m_pos.m_x = x;
		vInfo[j].m_pos.m_y = y;
		if(x < xmin) xmin = x;
		if(x > xmax) xmax = x;
		if(y < ymin) ymin = y;
		if(y > ymax) ymax = y;

		vInfo[j].m_adjBegin = adjCounter;
		for(adjEntry adj : v->adjEntries)
			adjLists[adjCounter++] = index[adj->twinNode()];
		vInfo[j].m_adjStop = adjCounter;

		++j;
	}

	m_xmin = xmin; m_xmax = xmax;
	m_ymin = ymin; m_ymax = ymax;
	m_wsum = wsum; m_hsum = hsum;

	m_master.syncThreads();

	if(m_id == 0)
		m_master.computeGrid(wsum, hsum, xmin, xmax, ymin, ymax);

	m_master.syncThreads();


	//---------------------------------
	// Main step
	//---------------------------------

	const bool noise = m_master.noise();
	Array<DPoint> &disp = m_master.disp();

	// random number generator for adding noise
	minstd_rand rng(randomSeed());
	uniform_real_distribution<> rand(0.75,1.25);

	//Array2D<ListPure<int>> &gridCell = m_master.gridCell();

	// --- Unfold Phase ---

	const ForceModelBase &forceModel = m_master.forceModel();

	const int numIter = m_master.numberOfIterations();
	for(int iter = 1; m_master.hasConverged() == false && iter <= numIter; ++iter)
	{
		double boxLength = m_master.boxLength();

		xmin = numeric_limits<double>::max(); xmax = -numeric_limits<double>::max();
		ymin = numeric_limits<double>::max(); ymax = -numeric_limits<double>::max();

		double sumForces = 0.0;
		double maxForce  = 0.0;

		double t = m_master.maxForceLength();
		double f = m_master.coolingFactor() * forceScaleFactor;

		for(int j = m_vStartIndex; j < m_vStopIndex; ++j) {
			NodeInfo &vj = vInfo[j];

			DPoint dp = forceModel.computeDisplacement(j, boxLength);

			// noise
			if(noise) {
				dp.m_x *= rand(rng);
				dp.m_y *= rand(rng);
			}

			double length = dp.norm();
			sumForces += length;
			updateMax(maxForce, length);

			//double s = (length <= t) ? forceScaleFactor : forceScaleFactor * t / length;
			double s = (length <= t) ? f : f * t / length;
			dp *= s;

			DPoint newPos = vj.m_pos + dp;

			// update new bounding box
			updateMin(xmin, newPos.m_x);
			updateMax(xmax, newPos.m_x);
			updateMin(ymin, newPos.m_y);
			updateMax(ymax, newPos.m_y);

			// store displacement
			disp[j] = dp;
		}

		m_xmin = xmin; m_xmax = xmax;
		m_ymin = ymin; m_ymax = ymax;
		m_sumForces = sumForces;
		m_maxForce  = maxForce;

		m_master.syncThreads();

		if(m_id == 0) {
			m_master.updateGridAndMoveNodes();
			m_master.coolDown();
		}

		m_master.syncThreads();
	}


	// --- Improvement Phase ---

	const int numIterImp = m_master.numberOfIterationsImprove();
	if(numIterImp > 0)
	{
		// scaling to ideal edge length
		scaling(vInfo,adjLists);

		const ForceModelBase &forceModelImprove = m_master.forceModelImprove();

		for(int iter = 1; !m_master.hasConverged() && iter <= numIterImp; ++iter)
		{
			double boxLength = m_master.boxLength();

			xmin = numeric_limits<double>::max(); xmax = -numeric_limits<double>::max();
			ymin = numeric_limits<double>::max(); ymax = -numeric_limits<double>::max();

			double sumForces = 0.0;
			double maxForce  = 0.0;

			double t = m_master.maxForceLength();
			double f = m_master.coolingFactor() * forceScaleFactor;

			for(int j = m_vStartIndex; j < m_vStopIndex; ++j) {
				NodeInfo &vj = vInfo[j];

				DPoint dp = forceModelImprove.computeDisplacement(j, boxLength);

				// noise
				if(noise) {
					dp.m_x *= rand(rng);
					dp.m_y *= rand(rng);
				}

				double length = dp.norm();
				sumForces += length;
				updateMax(maxForce, length);

				//double s = (length <= t) ? forceScaleFactor : forceScaleFactor * t / length;
				double s = (length <= t) ? f : f * t / length;
				dp *= s;

				DPoint newPos = vj.m_pos + dp;

				// update new bounding box
				updateMin(xmin, newPos.m_x);
				updateMax(xmax, newPos.m_x);
				updateMin(ymin, newPos.m_y);
				updateMax(ymax, newPos.m_y);

				// store displacement
				disp[j] = dp;
			}

			m_xmin = xmin; m_xmax = xmax;
			m_ymin = ymin; m_ymax = ymax;
			m_sumForces = sumForces;
			m_maxForce  = maxForce;

			m_master.syncThreads();

			// last iteration?
			if(iter == numIterImp) {

				for(int j = m_vStartIndex; j < m_vStopIndex; ++j) {
					NodeInfo &vj = vInfo[j];
					vj.m_pos += disp[j];
				}

			} else if(m_id == 0) {
				m_master.updateGridAndMoveNodes();
				m_master.coolDown();
			}

			m_master.syncThreads();
		}
	}


	//---------------------------------
	// Compute final layout
	//---------------------------------

	// scale layout to ideal edge length and compute bounding box
	finalScaling(vInfo, adjLists);

	if(m_id == 0)
		m_master.computeFinalBB();

	m_master.syncThreads();

	xmin = m_master.xleft();
	ymin = m_master.ysmall();

	node v = m_vStart;
	for(int j = m_vStartIndex; j < m_vStopIndex; v = v->succ(), ++j) {
		node vOrig = gc.original(v);

		ga.x(vOrig) = vInfo[j].m_pos.m_x - xmin;
		ga.y(vOrig) = vInfo[j].m_pos.m_y - ymin;
	}
}


} // end namespace ogdf
