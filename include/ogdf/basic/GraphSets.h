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
#include <ogdf/basic/internal/copy_move.h>

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
 */
class OGDF_EXPORT NodeSet : public RegisteredSet<internal::GraphNodeRegistry> {
public:
	//! Creates a new node set associated with \p graph.
	explicit NodeSet(const Graph& graph)
		: RegisteredSet((const internal::GraphNodeRegistry&)graph) { }

	using RegisteredSet::RegisteredSet;
	OGDF_DEFAULT_COPY(NodeSet);

	//! Returns a reference to the list of nodes contained in this set.
	const list_type& nodes() { return elements(); }

	//! Returns the associated graph.
	const Graph& graphOf() const {
		OGDF_ASSERT(registeredAt());
		return *registeredAt()->graphOf();
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
 */
class OGDF_EXPORT EdgeSet : public RegisteredSet<internal::GraphEdgeRegistry> {
public:
	//! Creates a new edge set associated with \p graph.
	explicit EdgeSet(const Graph& graph)
		: RegisteredSet((const internal::GraphEdgeRegistry&)graph) { }

	using RegisteredSet::RegisteredSet;
	OGDF_DEFAULT_COPY(EdgeSet);

	//! Returns a reference to the list of edges contained in this set.
	const list_type& edges() { return elements(); }

	//! Returns the associated graph.
	const Graph& graphOf() const {
		OGDF_ASSERT(RegisteredSet::registeredAt());
		return *registeredAt()->graphOf();
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
 */
class OGDF_EXPORT AdjEntrySet : public RegisteredSet<internal::GraphAdjRegistry> {
public:
	//! Creates a new adjEntry set associated with \p graph.
	explicit AdjEntrySet(const Graph& graph)
		: RegisteredSet((const internal::GraphAdjRegistry&)graph) { }

	using RegisteredSet::RegisteredSet;
	OGDF_DEFAULT_COPY(AdjEntrySet);

	//! Returns a reference to the list of adjEntries contained in this set.
	const list_type& adjEntries() { return elements(); }

	//! Returns the associated graph.
	const Graph& graphOf() const {
		OGDF_ASSERT(RegisteredSet::registeredAt());
		return *registeredAt()->graphOf();
	}
};

template<OGDF_NODE_LIST NL>
std::pair<int, int> Graph::insert(const NL& nodeList, const EdgeSet& edgeSet,
		NodeArray<node>& nodeMap, EdgeArray<edge>& edgeMap) {
	using std::size;
	m_regNodeArrays.reserveSpace(size(nodeList));
	m_regEdgeArrays.reserveSpace(size(edgeSet));
	m_regAdjArrays.reserveSpace(size(edgeSet));
	using std::begin;
	using std::end;
	return insert(begin(nodeList), end(nodeList), edgeSet, nodeMap, edgeMap);
}

}
