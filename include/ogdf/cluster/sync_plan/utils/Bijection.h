/** \file
 * \brief TODO Document
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
#include <ogdf/cluster/sync_plan/basic/Iterators.h>

using namespace ogdf;

using PipeBijIterator = ZipIterator<internal::GraphObjectContainer<AdjElement>::iterator,
		internal::GraphObjectContainer<AdjElement>::reverse_iterator>;
using PipeBijRange = Range<PipeBijIterator>;
using PipeBijPair = std::pair<adjEntry, adjEntry>;
using FrozenPipeBijPair = std::pair<int, int>;
using PipeBij = List<PipeBijPair>; // TODO convert to vector?
using FrozenPipeBij = List<FrozenPipeBijPair>;

OGDF_DECLARE_COMPARER(PipeBijCmp, PipeBijPair, int, x.first->theEdge()->index());

OGDF_DECLARE_COMPARER(FrozenPipeBijCmp, FrozenPipeBijPair, int, x.first);

PipeBijRange getPipeBijection(node u, node v);

void getPipeBijection(node u, node v, PipeBij& out);

void getPipeBijection(node u, node v, AdjEntryArray<adjEntry>& out);

void getPipeBijection(node u, node v, EdgeArray<edge>& out);

void getFrozenPipeBijection(node u, node v, FrozenPipeBij& out);

void freezePipeBijection(const PipeBij& in, FrozenPipeBij& out);
