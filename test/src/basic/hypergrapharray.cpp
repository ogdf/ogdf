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
#include <ogdf/hypergraph/Hypergraph.h>

#include "array_helper.h"

go_bandit([]() {
	auto chooseHypernode = [](const Hypergraph& H) { return H.randomHypernode(); };

	auto allHypernodes = [](const Hypergraph& H, List<hypernode>& list) { H.allHypernodes(list); };

	auto createHypernode = [](Hypergraph& H) { return H.newHypernode(); };

	auto chooseHyperedge = [](const Hypergraph& H) { return H.randomHyperedge(); };

	auto allHyperedges = [](const Hypergraph& H, List<hyperedge>& list) { H.allHyperedges(list); };

	auto createHyperedge = [](Hypergraph& H) {
		List<hypernode> newNodes;
		newNodes.pushBack(H.newHypernode());
		newNodes.pushBack(H.newHypernode());
		newNodes.pushBack(H.newHypernode());

		return H.newHyperedge(newNodes);
	};

	auto init = [&](Hypergraph& H) {
		H.clear();
		for (int i = 0; i < 42; ++i) {
			createHyperedge(H);
		}
	};

	describeArray<Hypergraph, HypernodeArray, hypernode, int>("HypernodeArray filled with ints", 42,
			43, init, chooseHypernode, allHypernodes, createHypernode);
	describeArray<Hypergraph, HypernodeArray, hypernode, List<int>>(
			"HypernodeArray filled with lists of ints", {1, 2, 3}, {42}, init, chooseHypernode,
			allHypernodes, createHypernode);

	describeArrayWithoutDefault<Hypergraph, HypernodeArrayWithoutDefault, hypernode,
			std::unique_ptr<int>>("HypernodeArray filled with unique pointers", init,
			chooseHypernode, allHypernodes, createHypernode);
	describeArrayWithoutDefault<Hypergraph, HypernodeArrayWithoutDefault, hypernode,
			std::vector<std::unique_ptr<int>>>("HypernodeArray filled with vectors of unique pointers",
			init, chooseHypernode, allHypernodes, createHypernode);

	describeArray<Hypergraph, HyperedgeArray, hyperedge, int>("HyperedgeArray filled with ints", 42,
			43, init, chooseHyperedge, allHyperedges, createHyperedge);
	describeArray<Hypergraph, HyperedgeArray, hyperedge, List<int>>(
			"HyperedgeArray filled with lists of ints", {1, 2, 3}, {42}, init, chooseHyperedge,
			allHyperedges, createHyperedge);

	describeArrayWithoutDefault<Hypergraph, HyperedgeArrayWithoutDefault, hyperedge,
			std::unique_ptr<int>>("HyperedgeArray filled with unique pointers", init,
			chooseHyperedge, allHyperedges, createHyperedge);
	describeArrayWithoutDefault<Hypergraph, HyperedgeArrayWithoutDefault, hyperedge,
			std::vector<std::unique_ptr<int>>>("HyperedgeArray filled with vectors of unique pointers",
			init, chooseHyperedge, allHyperedges, createHyperedge);
});
