/** \file
 * \brief Declaration and implementation of NodeSet, EdgeSet, and AdjEntrySet classes.
 *
 * \author Simon D. Fink, Matthias Pfretzschner
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h> // IWYU pragma: keep
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/basic.h>

#include <iterator>
#include <utility>

namespace ogdf {

//! Node sets.
/**
 * @ingroup graph-containers
 *
 * Maintains a subset of nodes contained in an associated graph.
 *
 * Provides efficient operations for testing membership,
 * iteration, insertion, and deletion of elements, as well as clearing the set.
 *
 * \tparam SupportFastSizeQuery Whether this set supports querying its #size in
 * constant instead of linear time (in the size).
 */
template<bool SupportFastSizeQuery = true>
class NodeSet : public RegisteredSet<internal::GraphNodeRegistry, SupportFastSizeQuery> {
	using RS = RegisteredSet<internal::GraphNodeRegistry, SupportFastSizeQuery>;

public:
	//! Creates a new node set associated with \p graph.
	explicit NodeSet(const Graph& graph) : RS((const internal::GraphNodeRegistry&)graph) { }

	OGDF_REGSET_CONSTR(NodeSet, RS)

	//! Returns a reference to the list of nodes contained in this set.
	const typename RS::list_type& nodes() { return RS::elements(); }

	//! Returns the associated graph.
	const Graph& graphOf() const {
		OGDF_ASSERT(RS::registeredAt());
		return *RS::registeredAt();
	}
};

//! Edge sets.
/**
 * @ingroup graph-containers
 *
 * Maintains a subset of edges contained in an associated graph.
 *
 * Provides efficient operations for testing membership,
 * iteration, insertion, and deletion of elements, as well as clearing the set.
 *
 * Note that an EdgeSet is not notified if an edge is hidden using HiddenEdgeSet. Thus, hidden
 * edges will stay part of EdgeSets even if they are (temporarily) no longer accessible as entry in
 * the edge list of their corresponding Graph.
 *
 * \tparam SupportFastSizeQuery Whether this set supports querying its #size in
 * constant instead of linear time (in the size).
 */
template<bool SupportFastSizeQuery = true>
class EdgeSet : public RegisteredSet<internal::GraphEdgeRegistry, SupportFastSizeQuery> {
	using RS = RegisteredSet<internal::GraphEdgeRegistry, SupportFastSizeQuery>;

public:
	//! Creates a new edge set associated with \p graph.
	explicit EdgeSet(const Graph& graph) : RS((const internal::GraphEdgeRegistry&)graph) { }

	OGDF_REGSET_CONSTR(EdgeSet, RS)

	//! Returns a reference to the list of edges contained in this set.
	const typename RS::list_type& edges() { return RS::elements(); }

	//! Returns the associated graph.
	const Graph& graphOf() const {
		OGDF_ASSERT(RS::registeredAt());
		return *RS::registeredAt();
	}
};

//! AdjEntry sets.
/**
 * @ingroup graph-containers
 *
 * Maintains a subset of adjEntries contained in an associated graph.
 *
 * Provides efficient operations for testing membership,
 * iteration, insertion, and deletion of elements, as well as clearing the set.
 *
 * \tparam SupportFastSizeQuery Whether this set supports querying its #size in
 * constant instead of linear time (in the size).
 */
template<bool SupportFastSizeQuery = true>
class AdjEntrySet : public RegisteredSet<internal::GraphAdjRegistry, SupportFastSizeQuery> {
	using RS = RegisteredSet<internal::GraphAdjRegistry, SupportFastSizeQuery>;

public:
	//! Creates a new adjEntry set associated with \p graph.
	explicit AdjEntrySet(const Graph& graph) : RS((const internal::GraphAdjRegistry&)graph) { }

	OGDF_REGSET_CONSTR(AdjEntrySet, RS)

	//! Returns a reference to the list of adjEntries contained in this set.
	const typename RS::list_type& adjEntries() { return RS::elements(); }

	//! Returns the associated graph.
	const Graph& graphOf() const {
		OGDF_ASSERT(RS::registeredAt());
		return *RS::registeredAt();
	}
};

template<OGDF_NODE_LIST NL>
std::pair<int, int> Graph::insert(const NL& nodeList, const EdgeSet<true>& edgeSet,
		NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap) {
	using std::size;
	m_regNodeArrays.reserveSpace(size(nodeList));
	m_regEdgeArrays.reserveSpace(size(edgeSet));
	m_regAdjArrays.reserveSpace(size(edgeSet));
	using std::begin;
	using std::end;
	return insert(begin(nodeList), end(nodeList), edgeSet, nodeMap, edgeMap);
}

template<OGDF_NODE_LIST NL>
std::pair<int, int> Graph::insert(const NL& nodeList, const EdgeSet<false>& edgeSet,
		NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap) {
	using std::size;
	m_regNodeArrays.reserveSpace(size(nodeList));
	using std::begin;
	using std::end;
	return insert(begin(nodeList), end(nodeList), edgeSet, nodeMap, edgeMap);
}

}
