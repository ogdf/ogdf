/** \file
 * \brief Tests for ogdf::AdjEntryArray
 *
 * \author Mirko Wagner, Tilo Wiedera
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
#include <ogdf/basic/AdjEntryArray.h>

#include "array_helper.h"

using namespace ogdf;
using namespace bandit;

go_bandit([]() {
	auto chooseAdjEntry = [](const Graph& graph) {
		edge e = graph.chooseEdge();
		return randomNumber(0, 1) ? e->adjSource() : e->adjTarget();
	};

	auto allAdjEntries = [](const Graph& graph, List<adjEntry>& list) {
		list.clear();

		for (edge e : graph.edges) {
			list.pushBack(e->adjSource());
			list.pushBack(e->adjTarget());
		}
	};

	auto createAdjEntry = [](Graph& graph) {
		edge e = graph.newEdge(graph.chooseNode(), graph.chooseNode());
		return randomNumber(0, 1) ? e->adjSource() : e->adjTarget();
	};

	auto init = [](Graph& graph) { randomGraph(graph, 42, 168); };

	describeArray<Graph, AdjEntryArray, adjEntry, int>( //
			"AdjEntryArray filled with ints", //
			42, 43, //
			init, chooseAdjEntry, allAdjEntries, createAdjEntry);
	describeArray<Graph, AdjEntryArray, adjEntry, List<int>>( //
			"AdjEntryArray filled with lists of ints", //
			{1, 2, 3}, {42}, //
			init, chooseAdjEntry, allAdjEntries, createAdjEntry);
	describeArray<Graph, AdjEntryArray, adjEntry, bool>( //
			"AdjEntryArray filled with bools", //
			false, true, //
			init, chooseAdjEntry, allAdjEntries, createAdjEntry);

	describeArrayWithoutDefault<Graph, AdjEntryArray, adjEntry, std::unique_ptr<int>>( //
			"AdjEntryArray filled with unique pointers", //
			init, chooseAdjEntry, allAdjEntries, createAdjEntry);
	describeArrayWithoutDefault<Graph, AdjEntryArray, adjEntry, std::vector<std::unique_ptr<int>>>( //
			"AdjEntryArray filled with vectors of unique pointers", //
			init, chooseAdjEntry, allAdjEntries, createAdjEntry);

	describe("AdjEntryArray", [&]() {
		it("keeps the correct values when splitting/unsplitting edges", [&]() {
			Graph G;
			AdjEntryArray<int> R(G, 0);
			node n1 = G.newNode();
			node n2 = G.newNode();
			edge e = G.newEdge(n1, n2);
			R[e->adjSource()] = 1;
			R[e->adjTarget()] = 2;

			edge e2 = G.split(e);

			AssertThat(R[e->adjSource()], Equals(1));
			AssertThat(R[e2->adjSource()], Equals(1));
			AssertThat(R[e->adjTarget()], Equals(2));
			AssertThat(R[e2->adjTarget()], Equals(2));

			R[e->adjTarget()] = 3;
			R[e2->adjSource()] = 4;

			G.unsplit(e, e2);

			AssertThat(R[e->adjSource()], Equals(1));
			AssertThat(R[e->adjTarget()], Equals(2));
		});
	});
});
