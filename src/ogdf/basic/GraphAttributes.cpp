/*
 * $Revision: 3521 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-05-31 14:52:33 +0200 (Fr, 31. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of class GraphAttributes.
 *
 * Class GraphAttributes extends a graph by graphical attributes like
 * node position, color, etc.
 *
 * \author Carsten Gutwenger, Karsten Klein
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


#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/fileformats/GmlParser.h>


namespace ogdf {

//---------------------------------------------------------
// GraphAttributes
// graph topology + graphical attributes
//---------------------------------------------------------

GraphAttributes::GraphAttributes() : m_pGraph(0), m_directed(true) { }



GraphAttributes::GraphAttributes(const Graph &G, long initAttr) :
	m_pGraph(&G), m_directed(true), m_attributes(0)
{
	initAttributes(m_attributes = initAttr);
}


void GraphAttributes::initAttributes(long attr)
{
	m_attributes |= attr;

	// no node style without graphics
	OGDF_ASSERT( (m_attributes & nodeGraphics) != 0 || (m_attributes & nodeStyle) == 0);
	// no edge style without graphics
	OGDF_ASSERT( (m_attributes & edgeGraphics) != 0 || (m_attributes & edgeStyle) == 0);

	if (attr & nodeGraphics) {
		m_x        .init( *m_pGraph, 0.0 );
		m_y        .init( *m_pGraph, 0.0 );
		m_width    .init( *m_pGraph, LayoutStandards::defaultNodeWidth () );
		m_height   .init( *m_pGraph, LayoutStandards::defaultNodeHeight() );
		m_nodeShape.init( *m_pGraph, LayoutStandards::defaultNodeShape () );
	}

	if (attr & threeD) {
		m_z.init(*m_pGraph, 0.0);
	}

	if (attr & nodeStyle)
	{
		m_nodeStroke.init( *m_pGraph, LayoutStandards::defaultNodeStroke() );
		m_nodeFill  .init( *m_pGraph, LayoutStandards::defaultNodeFill  () );
	}

	if (attr & edgeGraphics) {
		m_bends.init( *m_pGraph, DPolyline() );
	}

	if (attr & edgeStyle)
	{
		m_edgeStroke.init( *m_pGraph, LayoutStandards::defaultEdgeStroke() );
	}

	if (attr & nodeWeight) {
		m_nodeIntWeight.init( *m_pGraph, 0 );
	}
	if (attr & edgeIntWeight) {
		m_intWeight.init( *m_pGraph, 1 );
	}
	if (attr & edgeDoubleWeight) {
		m_doubleWeight.init( *m_pGraph, 1.0 );
	}
	if (attr & nodeLabel) {
		m_nodeLabel.init(*m_pGraph);
	}
	if (attr & edgeLabel) {
		m_edgeLabel.init(*m_pGraph);
	}
	if (attr & edgeType) {
		m_eType.init( *m_pGraph, Graph::association ); //should be Graph::standard and explicitly set
	}
	if (attr & nodeType) {
		m_vType.init( *m_pGraph, Graph::vertex );
	}
	if (attr & nodeId) {
		m_nodeId.init( *m_pGraph, -1 );
	}
	if (attr & edgeArrow) {
		m_edgeArrow.init( *m_pGraph, LayoutStandards::defaultEdgeArrow() );
	}
	if (attr & nodeTemplate) {
		m_nodeTemplate.init(*m_pGraph);
	}
	if (attr & edgeSubGraphs) {
		m_subGraph.init( *m_pGraph, 0 );
	}
}


void GraphAttributes::destroyAttributes(long attr)
{
	m_attributes &= ~attr;

	if (attr & nodeGraphics) {
		m_x     .init();
		m_y     .init();
		m_width .init();
		m_height.init();
		m_nodeShape.init();
		if (attr & nodeStyle)
		{
			m_nodeStroke.init();
			m_nodeFill  .init();
		}
	}

	if (attr & threeD) {
		m_z.init();
	}

	if (attr & edgeGraphics) {
		m_bends.init();
	}
	if (attr & edgeStyle)
	{
		m_edgeStroke.init();
	}

	if (attr & nodeWeight) {
		m_nodeIntWeight.init();
	}
	if (attr & edgeIntWeight) {
		m_intWeight.init();
	}
	if (attr & edgeDoubleWeight) {
		m_doubleWeight.init();
	}
	if (attr & nodeLabel) {
		m_nodeLabel.init();
	}
	if (attr & edgeLabel) {
		m_edgeLabel.init();
	}
	if (attr & nodeId) {
		m_nodeId.init();
	}
	if (attr & edgeArrow) {
		m_edgeArrow.init();
	}
	if (attr & nodeTemplate) {
		m_nodeTemplate.init();
	}
	if (attr & edgeSubGraphs) {
		m_subGraph.init();
	}
}


void GraphAttributes::init(const Graph &G, long initAttr)
{
	m_pGraph = &G;
	destroyAttributes(m_attributes);
	m_attributes = 0;
	initAttributes(m_attributes = initAttr);
}

void GraphAttributes::setAllWidth(double w)
{
	node v;
	forall_nodes(v,*m_pGraph)
		m_width[v] = w;
}


void GraphAttributes::setAllHeight(double h)
{
	node v;
	forall_nodes(v,*m_pGraph)
		m_height[v] = h;
}


void GraphAttributes::clearAllBends()
{
	edge e;
	forall_edges(e,*m_pGraph)
		m_bends[e].clear();
}


//
// calculates the bounding box of the graph
const DRect GraphAttributes::boundingBox() const
{
	double minx, maxx, miny, maxy;
	const Graph           &G  = constGraph();
	const GraphAttributes &AG = *this;
	node v = G.firstNode();

	if (v == 0) {
		minx = maxx = miny = maxy = 0.0;
	}
	else {
		minx = AG.x(v) - AG.width(v)/2;
		maxx = AG.x(v) + AG.width(v)/2;
		miny = AG.y(v) - AG.height(v)/2;
		maxy = AG.y(v) + AG.height(v)/2;

		forall_nodes(v, G) {
			double x1 = AG.x(v) - AG.width(v)/2;
			double x2 = AG.x(v) + AG.width(v)/2;
			double y1 = AG.y(v) - AG.height(v)/2;
			double y2 = AG.y(v) + AG.height(v)/2;

			if (x1 < minx) minx = x1;
			if (x2 > maxx) maxx = x2;
			if (y1 < miny) miny = y1;
			if (y2 > maxy) maxy = y2;
		}
	}

	edge e;
	forall_edges(e, G) {
		const DPolyline &dpl = AG.bends(e);
		ListConstIterator<DPoint> iter;
		for (iter = dpl.begin(); iter.valid(); ++iter) {
			if ((*iter).m_x < minx) minx = (*iter).m_x;
			if ((*iter).m_x > maxx) maxx = (*iter).m_x;
			if ((*iter).m_y < miny) miny = (*iter).m_y;
			if ((*iter).m_y > maxy) maxy = (*iter).m_y;
		}
	}

	return DRect(minx, miny, maxx, maxy);
}


//
// returns a list of all hierachies in the graph (a hierachy consists of a set of nodes)
// at least one list is returned, which is the list of all nodes not belonging to any hierachy
// this is always the first list
// the return-value of this function is the number of hierachies
int GraphAttributes::hierarchyList(List<List<node>* > &list) const
{
	// list must be empty during startup
	OGDF_ASSERT(list.empty());

	const Graph &G = constGraph();
	Array<bool> processed(0, G.maxNodeIndex(), false);
	node v;
	edge e;

	// initialize the first list of all single nodes
	List<node> *firstList = OGDF_NEW List<node>;
	list.pushBack(firstList);

	forall_nodes(v, G) { // scan all nodes

		// skip, if already processed
		if (processed[v->index()])
			continue;

		List<node> nodeSet;                    // set of nodes in this hierachy,
		// whose neighbours have to be processed
		List<node> *hierachy = OGDF_NEW List<node>; // holds all nodes in this hierachy

		nodeSet.pushBack(v);           // push the unprocessed node to the list
		processed[v->index()] = true;  // and mark it as processed

		do { // scan all neighbours of nodes in 'nodeSet'
			node v = nodeSet.popFrontRet();
			hierachy->pushBack(v); // push v to the list of nodes in this hierachy

			// process all the neighbours of v, e.g. push them into 'nodeSet'
			forall_adj_edges(e, v) {
				if (type(e) == Graph::generalization) {
					node w = e->source() == v ? e->target() : e->source();
					if (!processed[w->index()]) {
						nodeSet.pushBack(w);
						processed[w->index()] = true;
					}
				}
			}
		} while (!nodeSet.empty());

		// skip adding 'hierachy', if it contains only one node
		if (hierachy->size() == 1) {
			firstList->conc(*hierachy);
			delete hierachy;
		}
		else
			list.pushBack(hierachy);
	}

	return list.size() - 1 + (*list.begin())->size();
}


//
// returns a list of all hierarchies in the graph (in this case, a hierarchy consists of a set of edges)
// list may be empty, if no generalizations are used
// the return-value of this function is the number of hierarchies with generalizations
int GraphAttributes::hierarchyList(List<List<edge>* > &list) const
{
	// list must be empty during startup
	OGDF_ASSERT(list.empty());

	const Graph &G = constGraph();
	Array<bool> processed(0, G.maxNodeIndex(), false);
	node v;
	edge e;

	forall_nodes(v, G) { // scan all nodes

		// skip, if already processed
		if (processed[v->index()])
			continue;

		List<node> nodeSet;                    // set of nodes in this hierarchy,
		// whose neighbours have to be processed
		List<edge> *hierarchy = OGDF_NEW List<edge>; // holds all edges in this hierarchy

		nodeSet.pushBack(v);           // push the unprocessed node to the list
		processed[v->index()] = true;  // and mark it as processed

		do { // scan all neighbours of nodes in 'nodeSet'
			node v = nodeSet.popFrontRet();

			// process all the neighbours of v, e.g. push them into 'nodeSet'
			forall_adj_edges(e, v) {
				if (type(e) == Graph::generalization) {
					node w = e->source() == v ? e->target() : e->source();
					if (!processed[w->index()]) {
						nodeSet.pushBack(w);
						processed[w->index()] = true;
						hierarchy->pushBack(e); // push e to the list of edges in this hierarchy
					}
				}
			}
		} while (!nodeSet.empty());

		// skip adding 'hierarchy', if it contains only one node
		if (hierarchy->empty())
			delete hierarchy;
		else
			list.pushBack(hierarchy);
	}

	return list.size();
}



void GraphAttributes::removeUnnecessaryBendsHV()
{
	edge e;
	forall_edges(e,*m_pGraph)
	{
		DPolyline &dpl = m_bends[e];

		if(dpl.size() < 3)
			continue;

		ListIterator<DPoint> it1, it2, it3;

		it1 = dpl.begin();
		it2 = it1.succ();
		it3 = it2.succ();

		do {
			if(((*it1).m_x == (*it2).m_x && (*it2).m_x == (*it3).m_x) ||
				((*it1).m_y == (*it2).m_y && (*it2).m_y == (*it3).m_y))
			{
				dpl.del(it2);
				it2 = it3;
			} else {
				it1 = it2;
				it2 = it3;
			}

			it3 = it2.succ();
		} while(it3.valid());
	}
}


void GraphAttributes::addNodeCenter2Bends(int mode)
{
	edge e;
	forall_edges(e, *m_pGraph) {
		node v = e->source();
		node w = e->target();
		DPolyline &bendpoints = bends(e);
		switch (mode) {
		case 0 : // push center to the bends and return
			bendpoints.pushFront(DPoint(x(v), y(v)));
			bendpoints.pushBack (DPoint(x(w), y(w)));
			break;
		case 1 : // determine intersection with node and [center, last-bend-point]
			bendpoints.pushFront(DPoint(x(v), y(v)));
			bendpoints.pushBack (DPoint(x(w), y(w)));
		case 2 : // determine intersection between node and last bend-segment
			{
				DPoint sp1(x(v) - width(v)/2, y(v) - height(v)/2);
				DPoint sp2(x(v) - width(v)/2, y(v) + height(v)/2);
				DPoint sp3(x(v) + width(v)/2, y(v) + height(v)/2);
				DPoint sp4(x(v) + width(v)/2, y(v) - height(v)/2);
				DLine sourceRect[4] = {
					DLine(sp1, sp2),
					DLine(sp2, sp3),
					DLine(sp3, sp4),
					DLine(sp4, sp1)
				};

				DPoint tp1(x(w) - width(w)/2, y(w) - height(w)/2);
				DPoint tp2(x(w) - width(w)/2, y(w) + height(w)/2);
				DPoint tp3(x(w) + width(w)/2, y(w) + height(w)/2);
				DPoint tp4(x(w) + width(w)/2, y(w) - height(w)/2);
				DLine targetRect[4] = {
					DLine(tp1, tp2),
					DLine(tp2, tp3),
					DLine(tp3, tp4),
					DLine(tp4, tp1)
				};

				DRect source(sp1, sp3);
				DRect target(tp1, tp3);

				DPoint c1 = bendpoints.popFrontRet();
				DPoint c2 = bendpoints.popBackRet();

				while (!bendpoints.empty() && source.contains(bendpoints.front()))
					c1 = bendpoints.popFrontRet();
				while (!bendpoints.empty() && target.contains(bendpoints.back()))
					c2 = bendpoints.popBackRet();

				DPoint a1, a2;
				int i;
				if (bendpoints.size() == 0) {
					DLine cross(c1, c2);
					for (i = 0; i < 4; i++)
						if (cross.intersection(sourceRect[i], a1)) break;
					for (i = 0; i < 4; i++)
						if (cross.intersection(targetRect[i], a2)) break;
				}
				else {
					DLine cross1(c1, bendpoints.front());
					for (i = 0; i < 4; i++)
						if (cross1.intersection(sourceRect[i], a1)) break;
					DLine cross2(bendpoints.back(), c2);
					for (i = 0; i < 4; i++)
						if (cross2.intersection(targetRect[i], a2)) break;
				}
				bendpoints.pushFront(a1);
				bendpoints.pushBack(a2);
				break;
			}
			OGDF_NODEFAULT
		}
		bendpoints.normalize();
	}
}



//void GraphAttributes::writeSVG(const char *fileName, int fontSize, const string &fontColor) const
//	{
//		ofstream os(fileName);
//		writeSVG(os, fontSize, fontColor);
//	}
//
//void GraphAttributes::writeSVG(ostream &os, int fontSize, const string &fontColor) const
//{
//	os.setf(ios::showpoint);
//	os.precision(10);
//
//	os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
//	os << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" version=\"1.1\" baseProfile=\"full\" ";
//
//	// determine bounding box of svg
//	OGDF_ASSERT((*m_pGraph).numberOfNodes() > 0);
//	double maxX = x((*m_pGraph).firstNode());
//	double maxY = y((*m_pGraph).firstNode());
//	double minX = x((*m_pGraph).firstNode());
//	double minY = y((*m_pGraph).firstNode());
//	double nodeStrokeWidth;
//
//	node v;
//	forall_nodes(v, *m_pGraph) {
//		if (m_attributes & nodeStyle) {
//			nodeStrokeWidth = strokeWidth(v);
//		} else {
//			nodeStrokeWidth = 1.0;
//		}
//		maxX = max(maxX, x(v) + m_width[v]/2 + nodeStrokeWidth);
//		maxY = max(maxY, y(v) + m_height[v]/2 + nodeStrokeWidth);
//		minX = min(minX, x(v) - m_width[v]/2 - nodeStrokeWidth);
//		minY = min(minY, y(v) - m_height[v]/2 - nodeStrokeWidth);
//	}
//
//	edge e;
//	ListConstIterator<DPoint> it;
//	double edgeStrokeWidth;
//	forall_edges(e, *m_pGraph) {
//		if (m_attributes & edgeGraphics) {
//			if (attributes() & GraphAttributes::edgeStyle) {
//				edgeStrokeWidth = strokeWidth(e);
//			} else {
//				edgeStrokeWidth = 1.0;
//			}
//			const DPolyline &dpl = m_bends[e];
//			if (!dpl.empty()) {
//				for(it = dpl.begin(); it.valid(); ++it) {
//					maxX = max(maxX, (*it).m_x + edgeStrokeWidth);
//					maxY = max(maxY, (*it).m_y + edgeStrokeWidth);
//					minX = min(minX, (*it).m_x - edgeStrokeWidth);
//					minY = min(minY, (*it).m_y - edgeStrokeWidth);
//				}
//			}
//		}
//	}
//
//	os << "width=\"" << (maxX - minX) << "px\" ";
//	os << "height=\"" << (maxY - minY) << "px\" ";
//	os << "viewBox=\"" << 0 << " " << 0 << " " << (maxX - minX) << " " << (maxY - minY) << "\">\n";
//
//	forall_edges(e, *m_pGraph) {
//
//		const DPolyline &dpl = m_bends[e];
//		if (m_attributes & edgeGraphics) {
//			if (!dpl.empty()) { //polyline
//				os << "<polyline fill=\"none\" ";
//
//				if ((m_attributes & edgeColor) && (m_edgeColor[e].length() != 0)) {
//					os << "stroke=\"" << m_edgeColor[e] << "\" ";
//				}
//
//				if (attributes() & GraphAttributes::edgeStyle) {
//					os << "stroke-width=\"" << strokeWidth(e) << "px\" ";
//				} else {
//					os << "stroke=\"#000000\" ";
//				}
//
//				os << "points=\"";
//				node v = e->source();
//				if(dpl.front().m_x < m_x[v] - m_width[v]/2 ||
//						dpl.front().m_x > m_x[v] + m_width[v]/2 ||
//						dpl.front().m_y < m_y[v] - m_height[v]/2 ||
//						dpl.front().m_y > m_y[v] + m_height[v]/2)
//				{
//					os << (m_x[e->source()] - minX) << "," << (m_y[e->source()] - minY) << " ";
//				}
//
//				for(it = dpl.begin(); it.valid(); ++it)
//				os << ((*it).m_x - minX) << "," << ((*it).m_y - minY) << " ";
//
//				v = e->target();
//				if(dpl.back().m_x < m_x[v] - m_width[v]/2 ||
//						dpl.back().m_x > m_x[v] + m_width[v]/2 ||
//						dpl.back().m_y < m_y[v] - m_height[v]/2 ||
//						dpl.back().m_y > m_y[v] + m_height[v]/2)
//				{
//					os << (m_x[e->target()] - minX) << "," << (m_y[e->target()] - minY) << " ";
//				}
//
//				os << "\"/>\n";
//			} else { // single line
//				os << "<line ";
//				os << "x1=\"" << x(e->source()) - minX << "\" ";
//				os << "y1=\"" << y(e->source()) - minY << "\" ";
//				os << "x2=\"" << x(e->target()) - minX << "\" ";
//				os << "y2=\"" << y(e->target()) - minY<< "\" ";
//
//				if ((m_attributes & edgeColor) && (m_edgeColor[e].length() != 0)) {
//					os << "stroke=\"" << m_edgeColor[e] << "\" ";
//				} else {
//					os << "stroke=\"#000000\" ";
//				}
//
//				if (attributes() & GraphAttributes::edgeStyle) {
//					os << "stroke-width=\"" << strokeWidth(e) << "px\" ";
//				}
//
//				os << "/>\n";
//			}
//		}
//	}
//
//	forall_nodes(v,*m_pGraph) {
//		if (m_attributes & nodeGraphics) {
//			switch (m_nodeShape[v])
//			{
//				case shRect:
//				os << "<rect ";
//				os << "x=\"" << m_x[v] - minX - m_width[v]/2 << "\" ";
//				os << "y=\"" << m_y[v] - minY - m_height[v]/2 << "\" ";
//				os << "width=\"" << m_width[v] << "\" ";
//				os << "height=\"" << m_height[v] << "\" ";
//				break;
//				case shEllipse:
//				os << "<ellipse ";
//				os << "cx=\"" << m_x[v] - minX << "\" ";
//				os << "cy=\"" << m_y[v] - minY << "\" ";
//				os << "rx=\"" << m_width[v]/2 << "\" ";
//				os << "ry=\"" << m_height[v]/2 << "\" ";
//				break;
//			}
//
//			if (m_attributes & nodeColor) {
//				os << "fill=\"" << m_nodeColor[v] << "\" ";
//				os << "stroke=\"" << m_nodeLine[v] << "\" ";
//			}
//
//			if (m_attributes & nodeStyle)
//			{
//				os << "stroke-width=\"" << strokeWidth(v) << "px\" ";
//			}
//
//			os << "/>\n";
//
//			if(m_attributes & nodeLabel){
//				os << "<text x=\"" << m_x[v] - minX - m_width[v]/2 << "\" y=\"" << m_y[v] - minY << "\" textLength=\"" << m_width[v] << "\" font-size=\"" << fontSize << "\" fill=\"" << fontColor << "\" lengthAdjust=\"spacingAndGlyphs\">" << m_nodeLabel[v] << "</text>\n";
//			}
//		}
//	}
//
//	os << "</svg>\n";
//}

} // end namespace ogdf
