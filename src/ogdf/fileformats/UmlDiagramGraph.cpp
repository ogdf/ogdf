/*
 * $Revision: 2966 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-11-05 21:26:11 +0100 (Mo, 05. Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Implementation of the class UmlDiagramGraph
 *
 * \author Dino Ahr
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


#include <ogdf/fileformats/UmlDiagramGraph.h>


namespace ogdf {

	//
	// C o n s t r u c t o r
	//
	UmlDiagramGraph::UmlDiagramGraph(const UmlModelGraph &umlModelGraph,
											 UmlDiagramType diagramType,
											 const string &diagramName):
		m_modelGraph(umlModelGraph),
		m_diagramName(diagramName),
		m_diagramType(diagramType)
	{
	}

	//
	// D e s t r u c t o r
	//
	UmlDiagramGraph::~UmlDiagramGraph()
	{
		// Remove elements from lists
		m_containedNodes.clear();
		m_containedEdges.clear();
		m_x.clear();
		m_y.clear();
		m_w.clear();
		m_h.clear();
	}

	//
	// a d d N o d e W i t h G e o m e t r y
	//
	void UmlDiagramGraph::addNodeWithGeometry(
		NodeElement* node,
		double x, double y,
		double w, double h)
	{
		// Append node to the end of the list
		m_containedNodes.pushBack(node);

		// Dito with coordinates
		m_x.pushBack(x);
		m_y.pushBack(y);
		m_w.pushBack(w);
		m_h.pushBack(h);

	}

	//
	// a d d E d g e
	//
	void UmlDiagramGraph::addEdge(EdgeElement* edge)
	{
		// Append edge to the end of the list
		m_containedEdges.pushBack(edge);
	}

	//
	// g e t D i a g r a m T y p e S t r i n g
	//
	const char *UmlDiagramGraph::getDiagramTypeString() const
	{
		switch(m_diagramType){

		case (classDiagram):
			return "Class diagram";
			break;
		case (moduleDiagram):
			return "Module diagram";
			break;
		case (sequenceDiagram):
			return "Sequence diagram";
			break;
		case (collaborationDiagram):
			return "Collaboration diagram";
			break;
		case (componentDiagram):
			return "Component diagram";
			break;
		case (unknownDiagram):
			return "Unknown type diagram";
			break;
		default:
			return "";
		}

	} // getDiagramTypeString



	//
	// o u t p u t O p e r a t o r  for UmlDiagramGraph
	//
	ostream &operator<<(ostream &os, const UmlDiagramGraph &diagramGraph)
	{
		// Header with diagram name and type
		os << "\n--- " << diagramGraph.getDiagramTypeString()
			<< " \"" << diagramGraph.m_diagramName << "\" ---\n" << endl;

		// Nodes

		// Initialize iterators
		SListConstIterator<NodeElement*> nodeIt = diagramGraph.m_containedNodes.begin();
		SListConstIterator<double> xIt = diagramGraph.m_x.begin();
		SListConstIterator<double> yIt = diagramGraph.m_y.begin();
		SListConstIterator<double> wIt = diagramGraph.m_w.begin();
		SListConstIterator<double> hIt = diagramGraph.m_h.begin();

		// Traverse lists
		while (nodeIt.valid()){

			os << "Node " << diagramGraph.m_modelGraph.getNodeLabel(*nodeIt)
				<< " with geometry ("
				<< *xIt << ", "
				<< *yIt << ", "
				<< *wIt << ", "
				<< *hIt << ")." << endl;

			++nodeIt;
			++xIt;
			++yIt;
			++wIt;
			++hIt;

		} // while

		// Edges

		// Traverse lists
		SListConstIterator<EdgeElement*> edgeIt = diagramGraph.m_containedEdges.begin();
		for (edgeIt = diagramGraph.m_containedEdges.begin();
			 edgeIt.valid();
			 ++edgeIt)
		{
			os << "Edge between "
				<< diagramGraph.m_modelGraph.getNodeLabel((*edgeIt)->source())
				<< " and "
				<< diagramGraph.m_modelGraph.getNodeLabel((*edgeIt)->target())
				<< endl;
		}

		return os;

	} // <<


} // namespace ogdf
