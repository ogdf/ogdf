/** \file
 * \brief In the Geometric Edge Insertion Approach we iteratively insert a given set of edges into a drawing of a graph.
 * In each insertion the end vertices are moved to their optimal position.
 * The optimal position can be for example a position that minimizes the number of crossings for this vertex.
 *
 * Based on the implementation and techniques of the following papers:
 *
 * Marcel Radermacher, Klara Reichard, Ignaz Rutter, and Dorothea Wagner.
 * Geometric Heuristics for Rectilinear Crossing Minimization.
 * Journal of Experimental Algorithmics 24:1, pages 1.12:1–1.12:21, 2019. doi:10.1145/3325861 .
 *
 * Marcel Radermacher and Ignaz Rutter.
 * Geometric Crossing Minimization - A Scalable Randomized Approach.
 * In Proceedings of the 27th Annual European Symposium on Algorithms (ESA’19).
 * Ed. by Michael A. Bender, Ola Svensson, and Grzegorz Herman. Leibniz International Proceedings in Informatics (LIPIcs), pages 76:1–76:16.
 * Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2019. doi: 10.4230/LIPIcs.ESA
 *
 * \author Marcel Radermacher
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
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
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#pragma once

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/LayoutModule.h>
#include <ogdf/geometric/VertexPositionModule.h>

namespace ogdf {

class OGDF_EXPORT GeometricEdgeInsertion : public LayoutModule {
public:
	//! Constructor, sets options to default values.
	GeometricEdgeInsertion(Graph& _g) : g(_g) { (void)g; }

	~GeometricEdgeInsertion() { }

	//! The main call to the algorithm. GA should have nodeGraphics attributes enabled.
	virtual void call(GraphAttributes& GA) override;

	//!sets the routine the compute the new position of a vertex
	void setVertexPosition(VertexPositionModule* opt_pos) { m_pos = opt_pos; }

	//! sets the set of edges that have to reinserted
	void setHiddenEdgeSet(List<edge>* edge_set) {
		m_edge_set = edge_set;
		//m_hidden_edges = hidden_edge_set;
	}

	//! sets the method to compute the initial layout of the computed (planar) subgraph
	void setInitialLayouter(ogdf::LayoutModule* initial_layout_module) {
		m_initial_layout_module = initial_layout_module;
	}

protected:
private:
	Graph& g;
	VertexPositionModule* m_pos = nullptr;
	LayoutModule* m_initial_layout_module = nullptr;

	List<edge>* m_edge_set = nullptr;


	OGDF_NEW_DELETE
};

}
