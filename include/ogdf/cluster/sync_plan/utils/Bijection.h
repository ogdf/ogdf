/** \file
 * \brief Utilities for working with the bijections between the edges incident to the two endpoints of a Pipe.
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/comparer.h>
#include <ogdf/cluster/sync_plan/basic/Iterators.h>

#include <functional>
#include <utility>

namespace ogdf::sync_plan {

using PipeBijIterator = ZipIterator<ogdf::internal::GraphObjectContainer<AdjElement>::iterator,
		ogdf::internal::GraphObjectContainer<AdjElement>::reverse_iterator>;
using PipeBijRange = Range<PipeBijIterator>;
using PipeBijPair = std::pair<adjEntry, adjEntry>;
using FrozenPipeBijPair = std::pair<int, int>;
using PipeBij = List<PipeBijPair>;
using FrozenPipeBij = List<FrozenPipeBijPair>;

OGDF_DECLARE_COMPARER(PipeBijCmp, PipeBijPair, int, x.first->theEdge()->index());

OGDF_DECLARE_COMPARER(FrozenPipeBijCmp, FrozenPipeBijPair, int, x.first);

OGDF_EXPORT PipeBijRange getPipeBijection(node u, node v);

OGDF_EXPORT void getPipeBijection(node u, node v, PipeBij& out);

OGDF_EXPORT void getPipeBijection(node u, node v, AdjEntryArray<adjEntry>& out);

OGDF_EXPORT void getPipeBijection(node u, node v, EdgeArray<edge>& out);

OGDF_EXPORT void getFrozenPipeBijection(node u, node v, FrozenPipeBij& out);

OGDF_EXPORT void freezePipeBijection(const PipeBij& in, FrozenPipeBij& out);

OGDF_EXPORT std::pair<node, node> split(Graph& G, sync_plan::PipeBij& bij,
		const EdgeArray<edge>* new_edges = nullptr, const EdgeArray<bool>* reverse_edges = nullptr,
		node src = nullptr, node tgt = nullptr);

OGDF_EXPORT void join(Graph& G, node u, node v, sync_plan::PipeBij& bij,
		List<bool>* reverse_v = nullptr);

OGDF_EXPORT void join(Graph& G, node u, node v, sync_plan::PipeBij& bij,
		const std::function<void(node)>& deleteNode, const std::function<void(edge)>& deleteEdge,
		List<bool>* reverse_v = nullptr);

}
