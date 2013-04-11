/*
 * $Revision: 2977 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-11-06 14:33:34 +0100 (Di, 06. Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Contains the class UmlDiagramGraph which represents one
 * particular diagram of the complete UML Model.
 *
 * Each diagram refers to the node and edge information of
 * UmlModelGraph. Essentially a diagram contains selected nodes
 * and edges of the model provides with additional geometric
 * information.
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_DINO_UML_DIAGRAM_GRAPH_H
#define OGDF_DINO_UML_DIAGRAM_GRAPH_H

#include <ogdf/fileformats/UmlModelGraph.h>
#include <ogdf/basic/SList.h>


namespace ogdf {

	//---------------------------------------------------------
	// U m l D i a g r a m G r a p h
	//---------------------------------------------------------
	/** Contains the class UmlDiagramGraph which represents one
	 *  particular diagram of the complete UML Model. Each diagram refers
	 *  to the node and edge information of UmlModelGraph. Essentially
	 *  a diagram contains selected nodes and edges of the model provides
	 *  with additional geometric information.
	 */
	class OGDF_EXPORT UmlDiagramGraph {

		friend ostream &operator<<(ostream&, const UmlDiagramGraph &);

	public:

		//---------------------------------------------------------
		// U m l D i a g r a m T y p e
		//---------------------------------------------------------
		/** This enum type represents the different diagram types of UML.
		 */
		enum UmlDiagramType{
			classDiagram,
			moduleDiagram,
			sequenceDiagram,
			collaborationDiagram,
			componentDiagram,
			unknownDiagram

		}; // enum UmlDiagramType

	private:

		/** Reference to the model graph. */
		const UmlModelGraph &m_modelGraph;

		/** The name of the diagram. */
		string m_diagramName;

		/** The type of diagram. */
		UmlDiagramType m_diagramType;

		/** This list holds pointer to the nodes contained in
		 * the represented diagram.
		 */
		SList<NodeElement*> m_containedNodes;

		/** This list holds pointer to the edges contained in
		 * the represented diagram.
		 */
		SList<EdgeElement*> m_containedEdges;

		/** This list contains the x-coordinates of the nodes
		 * contained in the represented diagram.
		 */
		SList<double> m_x;

		/** This list contains the y-coordinates of the nodes
		 * contained in the represented diagram.
		 */
		SList<double> m_y;

		/** This list contains the width of the nodes
		 * contained in the represented diagram.
		 */
		SList<double> m_w;

		/** This list contains the height of the nodes
		 * contained in the represented diagram.
		 */
		SList<double> m_h;

	public:

		/** Constructor. */
		UmlDiagramGraph(const UmlModelGraph &umlModelGraph,
							UmlDiagramType diagramType,
							const string &diagramName);

		/** Destructor. */
		~UmlDiagramGraph();

		/** Adds a node with the given coordinates. */
		void addNodeWithGeometry(NodeElement* node,
								 double x, double y, double w, double h);

		/** Adds an edge. */
		void addEdge(EdgeElement* edge);

		/** Returns the name of the diagram. */
		const string &getDiagramName() const{
			return m_diagramName;
		}

		/** Returns the type of the diagram as string. */
		const char *getDiagramTypeString() const;

		/** Access to contained nodes. */
		const SList<NodeElement*> &getNodes() const{
			return m_containedNodes;
		}

		/** Access to contained edges. */
		const SList<EdgeElement*> &getEdges() const{
			return m_containedEdges;
		}

		/** Access to x-coordinates. */
		const SList<double> &getX() const{
			return m_x;
		}

		/** Access to y-coordinates. */
		const SList<double> &getY() const{
			return m_y;
		}

		/** Access to width. */
		const SList<double> &getWidth() const{
			return m_w;
		}

		/** Access to height. */
		const SList<double> &getHeight() const{
			return m_h;
		}

	}; // class UmlDiagramGraph

	/** Output operator for UmlDiagramGraph. */
	ostream &operator<<(ostream &os, const UmlDiagramGraph &diagramGraph);


} // end namespace ogdf

#endif
