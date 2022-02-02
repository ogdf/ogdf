/** \file
 * \brief Tests for algorithms computing Maximum Adjacency Orderings
 *
 * \author Max Ilsen, Ivo Hedtke
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

#include <ogdf/graphalg/MaxAdjOrdering.h>
#include <ogdf/basic/graph_generators.h>

#include <graphs.h>
#include <testing.h>

void testAllMAOs(const Graph &G)
{
	Graph P;
	emptyGraph(P, G.numberOfNodes());
	MaxAdjOrdering perms;
	ListPure<ListPure<node>> allPerms;
	perms.calcAll(&P,&allPerms);

	MaxAdjOrdering m;
	ListPure<ListPure<node>> MAOs;
	m.calcAll(&G,&MAOs);

	AssertThat(m.testIfAllMAOs(&G,&MAOs,&allPerms), IsTrue());
}

void testMAOBfs(const Graph &G)
{
	MaxAdjOrdering m;
	ListPure<node> MAO;
	m.calcBfs(&G,&MAO);

	AssertThat(m.testIfMAO(&G,&MAO), IsTrue());
	AssertThat(m.testIfMAOBfs(&G,&MAO), IsTrue());
}

go_bandit([](){
	describe("Maximum Adjacency Orderings", [](){
		describe("calculate exactly all MAOs", [](){
			constexpr int MIN_N = 4;
			constexpr int MAX_N = 8;

			forEachGraphItWorks({GraphProperty::simple}, testAllMAOs,
				GraphSizes(MIN_N, MAX_N, 1), 0, MAX_N
			);
		});

		describe("calculate MAOs with correct lex-bfs tie breaking", [](){
			constexpr int MIN_N = 10;
			constexpr int MAX_N = 20;

			forEachGraphItWorks({GraphProperty::simple}, testMAOBfs,
				GraphSizes(MIN_N, MAX_N, 1), 0, MAX_N
			);
		});
	});
});
