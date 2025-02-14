/** \file
 * \brief Tests for ogdf::HypernodeArray and ogdf::HyperedgeArray
 *
 * \author Matthias Pfretzschner
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
#include <ogdf/basic/List.h>
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/hypergraph/Hypergraph.h>

#include <functional>
#include <string>

#include "array_helper.h"
#include <testing.h>

using HypernodeSet = RegisteredSet<HypergraphRegistry<HypernodeElement>>;
using HyperedgeSet = RegisteredSet<HypergraphRegistry<HyperedgeElement>>;

go_bandit([]() {
	auto chooseHypernode = [](const Hypergraph& H) { return H.randomHypernode(); };

	auto allHypernodes = [](const Hypergraph& H, List<hypernode>& list) { H.allHypernodes(list); };

	auto createHypernode = [](Hypergraph& H) { return H.newHypernode(); };

	auto deleteHypernode = [](Hypergraph& H, hypernode n) { H.delHypernode(n); };
	auto deleteHyperedge = [](Hypergraph& H, hyperedge e) { H.delHyperedge(e); };

	auto chooseHyperedge = [](const Hypergraph& H) { return H.randomHyperedge(); };

	auto allHyperedges = [](const Hypergraph& H, List<hyperedge>& list) { H.allHyperedges(list); };

	auto createHyperedge = [](Hypergraph& H) {
		List<hypernode> newNodes;
		newNodes.pushBack(H.randomHypernode());
		newNodes.pushBack(H.newHypernode());
		newNodes.pushBack(H.newHypernode());

		return H.newHyperedge(newNodes);
	};

	auto init = [&](Hypergraph& H) {
		H.clear();
		H.newHypernode();
		for (int i = 0; i < 42; ++i) {
			createHyperedge(H);
		}
	};
	auto clear = [&](Hypergraph& H) { H.clear(); };

	runBasicArrayTests<Hypergraph, HypernodeArray, hypernode>( //
			"HypernodeArray", init, chooseHypernode, allHypernodes, createHypernode);

	runBasicArrayTests<Hypergraph, HyperedgeArray, hyperedge>( //
			"HyperedgeArray", init, chooseHyperedge, allHyperedges, createHyperedge);

	runBasicSetTests<Hypergraph, HypernodeSet, hypernode>("HypergraphRegistry<HypernodeElement>",
			init, chooseHypernode, allHypernodes, createHypernode, deleteHypernode, clear);
	runBasicSetTests<Hypergraph, HyperedgeSet, hyperedge>("HypergraphRegistry<HyperedgeElement>",
			init, chooseHyperedge, allHyperedges, createHyperedge, deleteHyperedge, clear);
});
