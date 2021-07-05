/** \file
 * \brief Declaration and implementation of NodeSet, EdgeSet, and AdjEntrySet classes.
 *
 * \author Niko Fink, Matthias Pfretzschner
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

#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/RegisteredSet.h>

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
 * \tparam SupportFastSizeQuery Whether this set supports querying it's #size in
 * constant instead of linear time (in the size).
 */
template<bool SupportFastSizeQuery = true>
class NodeSet : public RegisteredSet<node, GraphRegistry<NodeElement>, SupportFastSizeQuery> {
	using RS = RegisteredSet<node, GraphRegistry<NodeElement>, SupportFastSizeQuery>;

public:
	explicit NodeSet(const Graph& graph) : RS((const GraphRegistry<NodeElement>&)graph) {};

	//! Returns a reference to the list of nodes contained in this set.
	const typename RS::ListType& nodes() { return RS::elements(); }

	//! Returns the associated graph
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
 * \tparam SupportFastSizeQuery Whether this set supports querying it's #size in
 * constant instead of linear time (in the size).
 */
template<bool SupportFastSizeQuery = true>
class EdgeSet : public RegisteredSet<edge, GraphRegistry<EdgeElement>, SupportFastSizeQuery> {
	using RS = RegisteredSet<edge, GraphRegistry<EdgeElement>, SupportFastSizeQuery>;

public:
	explicit EdgeSet(const Graph& graph) : RS((const GraphRegistry<EdgeElement>&)graph) {};

	//! Returns a reference to the list of edges contained in this set.
	const typename RS::ListType& edges() { return RS::elements(); }

	//! Returns the associated graph
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
 * \tparam SupportFastSizeQuery Whether this set supports querying it's #size in
 * constant instead of linear time (in the size).
 */
template<bool SupportFastSizeQuery = true>
class AdjEntrySet
	: public RegisteredSet<adjEntry, GraphRegistry<AdjElement, GraphAdjIterator>, SupportFastSizeQuery> {
	using RS =
			RegisteredSet<adjEntry, GraphRegistry<AdjElement, GraphAdjIterator>, SupportFastSizeQuery>;

public:
	explicit AdjEntrySet(const Graph& graph)
		: RS((const GraphRegistry<AdjElement, GraphAdjIterator>&)graph) {};

	//! Returns a reference to the list of adjEntries contained in this set.
	const typename RS::ListType& adjEntries() { return RS::elements(); }

	//! Returns the associated graph
	const Graph& graphOf() const {
		OGDF_ASSERT(RS::registeredAt());
		return *RS::registeredAt();
	}
};

}
