/*
 * $Revision: 3960 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2014-03-13 11:36:28 +0100 (Thu, 13 Mar 2014) $
 ***************************************************************/

/** \file
 * \brief Declaration of class TopologyModule.
 *
 * The TopologyModule transports the layout information from
 * GraphAttributes to PlanRep on that Graph, i.e., it computes a
 * combinatorial embedding for the input.
 *
 * \author Karsten Klein
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


#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_TOPOLOGYMODULE_H
#define OGDF_TOPOLOGYMODULE_H



#include <ogdf/planarity/PlanRep.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/EdgeComparer.h>

namespace ogdf {

class EdgeLeg;

//===============================================
//main function(s):
//		setEmbeddingFromGraph(PlanRep&, GraphAttributes&)
//	    assumes that PG(AG) without bend nodes in PG
//
//		sortEdgesFromLayout(GraphAttributes &AG)
//      sort the edges in AG, no crossing insertion
//===============================================

//! Constructs embeddings from given layout.
/**
 * This class comprises functions for constructing the combinatorial embedding of a graph
 * or a planarized representation from a given layout.
 *
 * The main functions of the class are the following:
 *   - setEmbeddingFromGraph(PlanRep &PG, GraphAttributes &GA)
 *   - sortEdgesFromLayout(GraphAttributes &GA)
 */
class OGDF_EXPORT TopologyModule
{
public:
	TopologyModule() : m_options(opDegOneCrossings | opGenToAss |
		opCrossFlip | opLoop | opFlipUML) {}
	virtual ~TopologyModule() {}

	//! The (pre/post)processing options
	/**
	 * opCrossFlip increases running time by constant * n,
	 * opLoop increases running time by constant * m
	 */
	enum Options {
		opDegOneCrossings = 0x0001, //!< should degree one node's edge be crossed
		opGenToAss        = 0x0002, //!< should generalizations be turned into associations
		opCrossFlip       = 0x0004, //!< if there is a crossing between two edges with the same start or end point,
		                            //!< should their position at the node be flipped and the crossing be skipped?
		opFlipUML         = 0x0010, //!< only flip if same edge type
		opLoop            = 0x0008  //!< should loops between crossings (consecutive on both crossing edges) be deleted
		                            //!< (we dont check for enclosed CC's; therefore it is safe to remove the crossing).
	};

	void setOptions(int i) { m_options = i; }
	void addOption(TopologyModule::Options o)  { m_options = (m_options | o); }

	//! Uses the layout \a GA to determine an embedding for \a PG.
	/**
	 * Non-constness of GA in the following methods is only used when resolving problems,
	 * e.g., setting edge types if two generalizations cross in the input layout
	 *
	 * @param PG  is the input graph.
	 * @param GA  is the input layout.
	 * @param adjExternal  is assigned the external face (if \a setExternal is true).
	 * @param setExternal if true, we run over faces to compute the external face.
	 * @param reuseGAEmbedding If true, the call only checks for a correct embedding of \a PG
	 *                         and tries to insert crossings detected in the given layout otherwise.
	 *                         This allows to assign already sorted UMLGraphs.
	 *                         NOTE: if the sorting of the edges does not correspond to the layout given in \a GA,
	 *                         this cannot work correctly
	 * @return false if planarization fails; true otherwise.
	 */
	bool setEmbeddingFromGraph(
		PlanRep &PG,
		GraphAttributes &GA,
		adjEntry &adjExternal,
		bool setExternal = true,
		bool reuseGAEmbedding = false);

	//! Sorts the edges around all nodes of \a GA corresponding to the layout given in \a GA.
	/**
	 * There is no check of the embedding afterwards because this method could be used as a first step of a planarization
	 *
	 * @param G  is the input graph whose adjacency lists get sorted.
	 * @param GA is the input layout.
	 */
	void sortEdgesFromLayout(Graph &G, GraphAttributes &GA);

	face getExternalFace(PlanRep &PG, const GraphAttributes &AG);

	double faceSum(PlanRep &PG, const GraphAttributes &AG, face f);

protected:
	//compute a planarization, i.e. insert crossing vertices,
	//corresponding to the AG layout
	void planarizeFromLayout(PlanRep &PG,
							 GraphAttributes &AG);
	//compute crossing point and return if existing
	bool hasCrossing(EdgeLeg* legA, EdgeLeg* legB, DPoint& xp);
	//check if node v is a crossing of two edges with a common
	//endpoint adjacent to v, crossing is removed if flip is set
	bool checkFlipCrossing(PlanRep &PG,node v, bool flip = true);

	void postProcess(PlanRep &PG);
	void handleImprecision(PlanRep &PG);
	bool skipable(EdgeLeg* legA, EdgeLeg* legB);

private:
	//compare vectors for sorting
	int compare_vectors(const double& x1, const double& y1,
						const double& x2, const double& y2);
	//compute and return the angle defined by p-q,p-r
	double angle(DPoint p, DPoint q, DPoint r);
	//we have to save the position of the inserted crossing vertices
	//in order to compute the external face
	NodeArray<DPoint> m_crossPosition;

	//we save a list of EdgeLegs for all original edges in
	//AG
	EdgeArray< List<EdgeLeg*> > m_eLegs;


	//option settings as bits
	int m_options;

};//TopologyModule


//! Sorts EdgeLegs according to their xp distance to a reference point.
class PointComparer {
public:
	PointComparer(DPoint refPoint) : m_refPoint(refPoint) {}
	//compares the crossing points stored in the EdgeLeg
	int compare(const ListIterator<EdgeLeg*> &ep1,
		const ListIterator<EdgeLeg*> &ep2) const;
	OGDF_AUGMENT_COMPARER(ListIterator<EdgeLeg*>)

	// undefined methods to avoid automatic creation
	PointComparer &operator=(const PointComparer &);

private:
	const DPoint m_refPoint;
};//PointComparer

//! Helper class for the computation of crossings
/**
 * Represents a part of the edge between two consecutive bends (in the layout, there are no bends
 *  allowed in the representation) or crossings. There can be multiple EdgeLegs associated to one
 * copy edge in the PlanRep because of bends.
 */
class EdgeLeg
{
public:
	EdgeLeg() : m_topDown(false)
		{m_copyEdge = 0; m_xp = DPoint(0.0, 0.0);}
	EdgeLeg(edge e, int number, DPoint p1, DPoint p2) :
		m_xp(DPoint(0.0,0.0)), m_topDown(false),
		m_copyEdge(e), m_p1(p1), m_p2(p2), m_number(number)
		{}

	DPoint& start() {return m_p1;}
	DPoint& end()   {return m_p2;}
	int& number()    {return m_number;}
	edge& copyEdge() { return m_copyEdge;}

	//to avoid sorting both edgelegs and crossing points,
	//do not store a pair of them, but allow the xp to be
	//stored in the edgeleg
	DPoint m_xp;
	//we store the direction of the crossed EdgeLeg, too
	//if crossingEdgeLeg is horizontally left to right
	bool m_topDown;

	//each edgeLeg holds an entry with a ListIterator pointing to
	//its entry in a <edgeLeg*> List for an original edge
	ListIterator< EdgeLeg* > m_eIterator;

private:
	edge m_copyEdge; //the edge in the PlanRep copy corresponding to the EdgeLeg
	DPoint m_p1, m_p2;  //"Starting" and "End" point of the EdgeLeg
	int    m_number;    //the order nuumber on the edge, starting at 0

};//edgeleg


} // end namespace ogdf

#endif
