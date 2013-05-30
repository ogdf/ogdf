/*
 * $Revision: 3505 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:47 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief A declaration of EdgeStandardRep class representing a graph
 *        representation of a hypergraph in the edge standard form.
 *
 * This class provides a kind of an intermediate repsenetation of a
 * hypergraph between Hypergraph and edge standard based layout classes.
 * It is derived from HypergraphObserver and provides some additional
 * functionality to handle edge standard representation of hypergraph.
 * It follows Observer design pattern.
 *
 * \author Ondrej Moris
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

#ifndef OGDF_EDGE_STANDARD_REP_H
#define OGDF_EDGE_STANDARD_REP_H

#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/hypergraph/Hypergraph.h>
#include <ogdf/hypergraph/HypergraphArray.h>
#include <ogdf/hypergraph/HypergraphObserver.h>

namespace ogdf {

//!< Enumeration class of possible edge standard representations.
/**
 * There are the following possibilities:
 *
 * a) clique - no new dummy nodes are introduced, for every hyperedge
 *             e = (v_1, ..., v_l), we add a cliqie K_l connecting
 *             hypernodes incident with e,
 * b) star   - for every hyperedge e = {v_1, ..., v_l} a single new dummy
 *             node v_e is introduced, moreover, v_e becomes the center of
 *             a new star connecting all hypernodes incident with e (ie.
 *             {v_1, v_e}, ..., {v_l, v_e} are added)
 * c) tree   - for every hyperedge e a minimal subcubic tree connecting all
 *             hypernodes incident with e together is added with all its
 *             nodes and edges, leaves of tree are hypernodes, all non-leaf
 *             nodes are newly introduced dummy nodes.
 */
class OGDF_EXPORT EdgeStandardType
{
public:

	enum Type {
		clique = 0x0001,
		star   = 0x0002,
		tree   = 0x0003
	};
};

//! Edge standard representation of hypergraphs.
class OGDF_EXPORT EdgeStandardRep : public HypergraphObserver
{
private:

	//!< The type of edge standard representation.
	EdgeStandardType::Type m_type;

	//!< The reference to the original hypergraph.
	const Hypergraph *m_hypergraph;

	//!< Edge standard representation of the hypergraph.
	Graph m_graphRep;

	//!< The map from representation nodes to hypernodes.
	NodeArray<hypernode> m_hypernodeMap;

	//!< The map from representation hypernodes to nodes.
	HypernodeArray<node> m_nodeMap;

	//!< The map from representation edge to hyperedges.
	EdgeArray<hyperedge> m_hyperedgeMap;

	//!< The map from representation hyperedge to edges.
	HyperedgeArray<List<edge> > m_edgeMap;

	//!< The list of all newly created nodes.
	List<node> m_dummyNodes;

public:

	//!< Creates an edge standard representation.
	EdgeStandardRep();

	//!< Creates an edge standard rep. of a given type associated with \a H.
	EdgeStandardRep(const Hypergraph &pH, EdgeStandardType::Type pType);

	//!< Desctructor.
	virtual ~EdgeStandardRep();

	//!< Clears all cluster data.
	void clear();

	//! Conversion to original hypergraph reference.
	const Hypergraph & hypergraph() const
	{
		return *m_hypergraph;
	}

	//! Returns a reference to the representation graph.
	const Graph & constGraph() const
	{
		return m_graphRep;
	}

	//!< Returns the type of edge standard representation.
	EdgeStandardType::Type type() const {
		return m_type;
	}

	//!< Returns the node associated with the hypernode.
	node nodeMap(hypernode v)
	{
		return m_nodeMap[v];
	}

	//!< Returns the hypernode associated with the node (if any).
	hypernode hypernodeMap(node v)
	{
		return m_hypernodeMap[v];
	}

	//!< Returns the list of edges associated with the hyperedge.
	const List<edge> & edgeMap(hyperedge e)
	{
		return m_edgeMap[e];
	}

	//!< Returns the hyperedge associated with the edge.
	hyperedge hyperedgeMap(edge e)
	{
		return m_hyperedgeMap[e];
	}

	//!< Returns the list of dummy nodes.
	const List<node> & dummyNodes() const
	{
		return m_dummyNodes;
	}

protected:

	//!< Hypernode removal reaction.
	virtual void hypernodeDeleted(hypernode v);

	//!< Hypernode addition reaction.
	virtual void hypernodeAdded(hypernode v);

	//!< Hyperedge removal reaction.
	virtual void hyperedgeDeleted(hyperedge e);

	//!< Hyperedge addition reaction.
	virtual void hyperedgeAdded(hyperedge e);

	//!< Hypergraph clean-up reaction.
	virtual void cleared();

private:

	void constructCliqueRep();

	void constructStarRep();

	void constructTreeRep();

	void hyperedgeToTree(hyperedge e, int degree);

	void hyperedgeToClique(hyperedge e);

	void cloneHypernodes();

};

} // end namespace ogdf

#endif
