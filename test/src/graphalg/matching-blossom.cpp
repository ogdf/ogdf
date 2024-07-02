/** \file
 * \brief Tests for the blossom matching algorithm.
 *
 * \author Joshua Sangmeister
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/graphalg/Matching.h>
#include <ogdf/graphalg/MatchingBlossom.h>
#include <ogdf/graphalg/MatchingBlossomV.h>

#include <algorithm>
#include <functional>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include <testing.h>

namespace ogdf {
template<typename TWeight>
class MatchingModule;
} // namespace ogdf

using namespace ogdf::Matching;

const double DELTA = 1e-8;

template<class TWeight>
void createGraph(Graph& g, EdgeArray<TWeight>& weights,
		std::vector<std::tuple<int, int, TWeight>>& edges) {
	List<std::pair<int, int>> rawEdges;
	int n = -1;
	for (auto e : edges) {
		int u = std::get<0>(e);
		int v = std::get<1>(e);
		n = std::max({n, u, v});
		rawEdges.pushBack({u, v});
	}
	customGraph(g, n + 1, rawEdges);
	auto it = g.edges.begin();
	for (auto e : edges) {
		weights[*it] = std::get<2>(e);
		it++;
	}
}

template<class TWeight>
void call(MatchingModule<TWeight>& m, ogdf::Graph& graph, EdgeArray<TWeight>& weights,
		TWeight expected) {
	std::unordered_set<edge> matching;
	bool result = m.minimumWeightPerfectMatching(graph, weights, matching);
	if (expected != -1) {
		AssertThat(result, IsTrue());
		AssertThat(isPerfectMatching(graph, matching), IsTrue());
		AssertThat(m.matchingWeight(matching, weights), EqualsWithDelta(expected, DELTA));
	} else {
		AssertThat(result, IsFalse());
	}
}

template<template<typename> class MatchingImpl, class TWeight>
void testWithGraph(ogdf::Graph& graph, EdgeArray<TWeight>& weights, TWeight expected) {
	MatchingImpl<TWeight> m;
	call(m, graph, weights, expected);
}

template<template<typename> class MatchingImpl, class TWeight>
void test(std::vector<std::tuple<int, int, TWeight>> edges, TWeight expected) {
	Graph graph;
	EdgeArray<TWeight> weights(graph);
	createGraph(graph, weights, edges);
	testWithGraph<MatchingImpl, TWeight>(graph, weights, expected);
}

template<template<typename> class MatchingImpl>
void runAllTests() {
	it("calculates the matching on an empty graph correctly", [&] { test<MatchingImpl>({}, 0); });
	it("finds no matching for a single node", [&] {
		Graph g;
		EdgeArray<int> weights(g);
		emptyGraph(g, 1);
		testWithGraph<MatchingImpl>(g, weights, -1);
	});
	it("finds no matching for an odd number of vertices", [&] {
		test<MatchingImpl>({{0, 1, 1}, {0, 2, 1}}, -1);
	});
	it("finds a simple matching with a single zero-weight edge", [&] {
		test<MatchingImpl>({{0, 1, 0}}, 0);
	});
	it("finds no matching with two zero-weight edges", [&] {
		test<MatchingImpl>({{0, 1, 0}, {1, 2, 0}}, -1);
	});
	it("finds no matching without edges", [&] {
		Graph g;
		EdgeArray<int> weights(g);
		emptyGraph(g, 4);
		testWithGraph<MatchingImpl>(g, weights, -1);
	});
	it("finds a simple matching with three zero-weight edges", [&] {
		test<MatchingImpl>({{0, 1, 0}, {1, 2, 0}, {2, 3, 0}}, 0);
	});
	it("finds a zero-cost matching with a single odd-length cycle", [&] {
		test<MatchingImpl>({{0, 1, 0}, {1, 2, 0}, {2, 0, 0}, {2, 3, 0}}, 0);
	});
	it("finds a simple matching with a single edge", [&] { test<MatchingImpl>({{0, 1, 1}}, 1); });
	it("finds a simple matching with multiple edges", [&] {
		test<MatchingImpl>({{0, 1, 1}, {0, 2, 5}, {2, 3, 2}}, 3);
	});
	it("calculates the matching on an example graph correctly", [&] {
		test<MatchingImpl>(
				{{0, 1, 3}, {0, 4, 5}, {0, 5, 4}, {1, 2, 7}, {1, 5, 3}, {1, 6, 2}, {2, 3, 1},
						{2, 6, 3}, {2, 7, 2}, {3, 7, 1}, {4, 5, 2}, {5, 6, 1}, {6, 7, 4}},
				9);
	});
	it("calculates the matching on a bipartite graph correctly", [&] {
		test<MatchingImpl>({{0, 4, 1}, {1, 4, 1}, {1, 5, 1}, {1, 2, 4}, {5, 6, 3}, {2, 6, 9},
								   {3, 6, 2}, {3, 7, 7}},
				15);
	});
	it("handles multiple shrinkings correctly", [&] {
		test<MatchingImpl>(
				{{0, 1, 1}, {0, 2, 2}, {1, 2, 1}, {2, 3, 3}, {2, 4, 4}, {3, 4, 3}, {4, 5, 5}}, 9);
	});
	it("calculates a matching on a small odd cycle correctly", [&] {
		test<MatchingImpl>({{0, 1, 8}, {1, 2, 8}, {2, 3, 8}, {3, 0, 8}, {1, 3, 2}}, 16);
	});
	it("calculates a matching with zero weight correctly", [&] {
		test<MatchingImpl>({{0, 1, 3}, {0, 2, 2}, {0, 3, 8}, {0, 4, 10}, {0, 5, 8}, {1, 2, 1},
								   {1, 3, 4}, {1, 4, 8}, {1, 5, 0}, {2, 3, 8}, {2, 4, 5}, {2, 5, 0},
								   {3, 4, 7}, {3, 5, 0}, {4, 5, 0}},
				6);
	});
	it("calculates the matching on a triangulation of 10 points correctly (1)", [&] {
		test<MatchingImpl>(
				{{5, 6, 11}, {5, 4, 15}, {4, 6, 10}, {7, 2, 11}, {2, 5, 11}, {6, 0, 14}, {5, 7, 16},
						{7, 4, 17}, {8, 7, 10}, {7, 1, 10}, {1, 4, 15}, {9, 8, 1}, {8, 1, 4},
						{3, 9, 20}, {9, 1, 4}, {4, 0, 7}, {3, 4, 8}, {1, 3, 16}, {3, 0, 3}},
				35);
	});
	it("calculates the matching on a triangulation of 10 points correctly (2)", [&] {
		test<MatchingImpl>({{1, 9, 4}, {7, 8, 9}, {9, 7, 7}, {8, 9, 13}, {0, 6, 8}, {4, 0, 10},
								   {6, 4, 3}, {7, 0, 12}, {4, 7, 7}, {1, 7, 10}, {4, 8, 6},
								   {0, 1, 20}, {6, 2, 8}, {2, 4, 7}, {2, 3, 12}, {3, 0, 8},
								   {5, 0, 2}, {3, 5, 5}, {2, 8, 5}, {3, 6, 7}, {6, 5, 7}},
				25);
	});
	it("calculates the matching on a triangulation of 10 points correctly (3)", [&] {
		test<MatchingImpl>(
				{{2, 5, 21}, {6, 2, 14}, {5, 6, 7}, {4, 6, 11}, {3, 4, 15}, {6, 3, 18}, {0, 6, 17},
						{5, 0, 18}, {0, 3, 8}, {8, 1, 6}, {7, 8, 7}, {1, 7, 9}, {1, 3, 14},
						{0, 1, 22}, {8, 3, 8}, {3, 7, 12}, {7, 4, 9}, {4, 2, 6}, {9, 7, 3},
						{1, 9, 12}, {9, 2, 11}, {4, 9, 7}, {2, 1, 23}},
				30);
	});
	it("calculates the matching on a triangulation of 10 points correctly (4)", [&] {
		test<MatchingImpl>({{9, 1, 9}, {6, 1, 19}, {7, 6, 19}, {1, 7, 14}, {5, 0, 3}, {2, 5, 14},
								   {0, 2, 11}, {4, 9, 8}, {1, 3, 5}, {9, 7, 8}, {4, 2, 6},
								   {2, 9, 9}, {5, 4, 19}, {7, 2, 3}, {0, 7, 10}, {8, 5, 9},
								   {0, 6, 17}, {6, 3, 17}, {8, 0, 10}, {6, 8, 11}},
				30);
	});
	it("calculates the matching on a triangulation of 10 points correctly (5)", [&] {
		test<MatchingImpl>({{6, 1, 7}, {4, 6, 7}, {1, 4, 7}, {5, 0, 8}, {1, 5, 7}, {0, 1, 3},
								   {7, 8, 4}, {9, 7, 3}, {8, 9, 3}, {5, 7, 10}, {9, 5, 11},
								   {8, 4, 19}, {6, 5, 1}, {5, 2, 12}, {2, 0, 8}, {3, 9, 20},
								   {7, 4, 16}, {7, 6, 9}, {9, 2, 19}, {0, 3, 9}, {3, 2, 1}},
				24);
	});
	it("calculates the matching on a triangulation of 10 points correctly (6)", [&] {
		test<MatchingImpl>({{7, 3, 7}, {1, 9, 12}, {4, 1, 18}, {9, 4, 11}, {7, 8, 7}, {5, 7, 4},
								   {8, 5, 3}, {3, 4, 5}, {3, 0, 4}, {0, 4, 8}, {7, 0, 4},
								   {2, 4, 16}, {0, 2, 10}, {5, 0, 5}, {5, 2, 5}, {8, 6, 11},
								   {2, 8, 4}, {6, 2, 8}, {2, 1, 18}, {1, 6, 16}},
				32);
	});
	it("calculates the matching on a triangulation of 10 points correctly (7)", [&] {
		test<MatchingImpl>({{1, 7, 15}, {5, 0, 13}, {4, 5, 14}, {0, 4, 12}, {4, 6, 11}, {0, 6, 4},
								   {6, 2, 15}, {7, 5, 14}, {4, 7, 2}, {5, 2, 10}, {2, 0, 11},
								   {1, 9, 6}, {5, 1, 1}, {9, 5, 7}, {8, 3, 2}, {3, 1, 9}, {3, 9, 7},
								   {2, 8, 9}, {8, 9, 6}, {2, 9, 4}},
				13);
	});
	it("calculates the matching on a triangulation of 10 points correctly (8)", [&] {
		test<MatchingImpl>({{7, 6, 2}, {6, 1, 14}, {3, 6, 22}, {1, 3, 8}, {9, 8, 6}, {5, 9, 4},
								   {8, 5, 8}, {3, 4, 19}, {6, 5, 4}, {5, 1, 15}, {7, 5, 3},
								   {2, 1, 15}, {8, 2, 14}, {1, 8, 16}, {0, 9, 9}, {2, 3, 17},
								   {9, 7, 6}, {0, 8, 7}, {4, 2, 4}, {0, 4, 20}, {2, 0, 18}},
				25);
	});
	it("calculates the matching on a triangulation of 10 points correctly (9)", [&] {
		test<MatchingImpl>({{6, 4, 9}, {6, 0, 13}, {8, 6, 11}, {0, 8, 4}, {5, 8, 12}, {0, 5, 16},
								   {9, 5, 8}, {2, 9, 2}, {5, 2, 8}, {9, 8, 9}, {2, 1, 8}, {7, 2, 4},
								   {1, 7, 6}, {3, 5, 7}, {1, 4, 3}, {4, 7, 6}, {3, 2, 4}, {1, 3, 7},
								   {7, 9, 5}, {9, 6, 6}, {7, 6, 4}},
				20);
	});
	it("calculates the matching on a triangulation of 10 points correctly (10)", [&] {
		test<MatchingImpl>({{3, 4, 16}, {9, 5, 8}, {1, 9, 7}, {5, 1, 5}, {0, 1, 1}, {8, 0, 8},
								   {1, 8, 9}, {0, 9, 7}, {9, 8, 11}, {2, 9, 7}, {8, 2, 13}, {9, 7, 4},
								   {7, 5, 11}, {4, 2, 15}, {8, 4, 25}, {2, 7, 5}, {6, 7, 2},
								   {2, 6, 5}, {5, 3, 18}, {3, 6, 6}, {2, 3, 9}, {6, 5, 13}},
				38);
	});
	it("calculates the matching on a triangulation of 20 points correctly", [&] {
		test<MatchingImpl>(
				{{6, 5, 41}, {4, 6, 37}, {5, 4, 59}, {16, 6, 16}, {4, 16, 29}, {13, 5, 47},
						{10, 13, 35}, {5, 10, 19}, {17, 1, 20}, {9, 17, 33}, {1, 9, 14}, {5, 15, 55},
						{13, 2, 41}, {7, 13, 21}, {2, 7, 43}, {17, 4, 44}, {13, 17, 46}, {4, 13, 51},
						{12, 11, 88}, {14, 0, 16}, {4, 14, 17}, {0, 4, 26}, {15, 12, 83},
						{0, 15, 69}, {16, 0, 40}, {15, 16, 30}, {6, 15, 24}, {19, 14, 26},
						{4, 19, 38}, {18, 13, 32}, {7, 18, 28}, {2, 10, 24}, {5, 2, 42},
						{18, 9, 19}, {1, 18, 12}, {8, 11, 2}, {9, 8, 3}, {11, 9, 2}, {8, 18, 20},
						{7, 8, 41}, {18, 17, 27}, {12, 0, 14}, {3, 12, 12}, {0, 3, 13}, {19, 17, 15},
						{9, 19, 47}, {11, 19, 49}, {12, 19, 39}, {14, 3, 14}, {3, 19, 28}},
				200);
	});
	it("calculates the matching on a triangulation of 50 points correctly", [&] {
		test<MatchingImpl>(
				{{37, 17, 72}, {49, 37, 67}, {17, 49, 46}, {39, 48, 141}, {21, 39, 155},
						{48, 21, 146}, {19, 11, 217}, {48, 13, 83}, {11, 48, 92}, {13, 11, 105},
						{13, 39, 139}, {9, 48, 46}, {11, 9, 89}, {46, 20, 71}, {21, 46, 109},
						{20, 21, 64}, {14, 28, 387}, {42, 14, 258}, {28, 42, 129}, {33, 24, 149},
						{34, 33, 121}, {24, 34, 52}, {33, 7, 200}, {7, 24, 129}, {38, 32, 75},
						{35, 38, 59}, {32, 35, 107}, {45, 3, 65}, {37, 35, 117}, {35, 17, 113},
						{15, 5, 65}, {29, 15, 44}, {5, 29, 42}, {35, 30, 145}, {30, 17, 45},
						{8, 4, 142}, {2, 8, 123}, {4, 2, 159}, {32, 26, 192}, {26, 35, 182},
						{9, 7, 158}, {21, 9, 143}, {7, 21, 133}, {44, 43, 120}, {6, 44, 98},
						{43, 6, 48}, {18, 23, 301}, {2, 18, 186}, {23, 2, 115}, {25, 46, 191},
						{21, 25, 169}, {20, 31, 169}, {39, 20, 150}, {31, 39, 285}, {24, 9, 68},
						{1, 24, 98}, {9, 1, 93}, {11, 1, 28}, {6, 25, 159}, {32, 6, 157},
						{25, 32, 189}, {46, 31, 128}, {25, 26, 34}, {25, 7, 65}, {7, 26, 53},
						{30, 49, 66}, {7, 30, 191}, {30, 26, 172}, {1, 34, 79}, {19, 1, 189},
						{19, 33, 81}, {34, 19, 140}, {19, 14, 211}, {47, 19, 89}, {14, 47, 148},
						{47, 49, 81}, {30, 47, 93}, {47, 33, 75}, {12, 37, 40}, {49, 12, 42},
						{30, 33, 81}, {49, 41, 7}, {41, 12, 36}, {47, 41, 82}, {14, 41, 147},
						{5, 6, 166}, {32, 5, 67}, {14, 37, 170}, {12, 14, 136}, {38, 14, 328},
						{42, 38, 97}, {38, 37, 159}, {0, 27, 6}, {38, 0, 44}, {27, 38, 39},
						{27, 32, 53}, {28, 0, 84}, {38, 28, 105}, {0, 32, 55}, {40, 10, 33},
						{28, 40, 67}, {10, 28, 77}, {10, 32, 63}, {0, 10, 31}, {15, 28, 94},
						{29, 10, 32}, {40, 29, 9}, {15, 40, 45}, {3, 15, 112}, {10, 5, 49},
						{3, 5, 116}, {3, 6, 176}, {36, 45, 4}, {45, 6, 150}, {18, 36, 52},
						{36, 6, 149}, {22, 16, 64}, {2, 22, 113}, {16, 2, 111}, {18, 44, 57},
						{6, 18, 129}, {22, 31, 62}, {31, 16, 78}, {43, 2, 157}, {4, 43, 16},
						{43, 25, 134}, {8, 46, 69}, {46, 4, 172}, {4, 25, 124}, {46, 16, 94},
						{16, 8, 38}, {44, 2, 135}, {31, 23, 75}, {23, 22, 21}},
				1601);
	});
	it("calculates the matching on a triangulation of 100 points correctly", [&] {
		test<MatchingImpl>(
				{{20, 79, 1125}, {56, 20, 147}, {79, 56, 1254}, {56, 42, 538}, {29, 56, 641},
						{42, 29, 399}, {40, 27, 834}, {26, 40, 866}, {27, 26, 250}, {98, 72, 1504},
						{40, 98, 326}, {72, 40, 1326}, {56, 48, 1271}, {94, 56, 308},
						{48, 94, 1579}, {22, 33, 69}, {12, 22, 136}, {33, 12, 156}, {93, 57, 909},
						{74, 93, 275}, {57, 74, 1184}, {9, 0, 919}, {17, 9, 1276}, {0, 17, 1090},
						{35, 40, 728}, {26, 35, 210}, {25, 89, 51}, {63, 25, 366}, {89, 63, 417},
						{9, 91, 1955}, {74, 9, 1524}, {91, 74, 432}, {67, 27, 661}, {77, 67, 1476},
						{27, 77, 928}, {40, 67, 398}, {77, 40, 1427}, {68, 83, 992}, {27, 40, 834},
						{35, 27, 394}, {12, 9, 1022}, {74, 12, 701}, {40, 60, 88}, {60, 67, 481},
						{60, 36, 1108}, {93, 60, 648}, {36, 93, 519}, {76, 72, 485}, {63, 76, 149},
						{72, 63, 340}, {3, 52, 600}, {73, 3, 743}, {52, 73, 1298}, {32, 92, 1485},
						{56, 32, 1411}, {92, 56, 434}, {41, 28, 498}, {6, 33, 1007}, {49, 6, 868},
						{33, 49, 218}, {17, 95, 1565}, {69, 17, 1552}, {95, 69, 476}, {20, 56, 147},
						{56, 79, 1254}, {41, 86, 638}, {86, 28, 909}, {56, 92, 434}, {92, 79, 1375},
						{92, 56, 434}, {1, 32, 1334}, {56, 1, 425}, {86, 62, 456}, {7, 86, 274},
						{62, 7, 186}, {0, 39, 231}, {39, 17, 1067}, {0, 24, 240}, {53, 0, 167},
						{24, 53, 398}, {24, 95, 2330}, {95, 53, 2017}, {39, 11, 529}, {11, 17, 611},
						{95, 18, 940}, {0, 95, 2166}, {18, 0, 2266}, {14, 0, 758}, {53, 14, 789},
						{51, 46, 1515}, {8, 51, 1000}, {46, 8, 613}, {27, 98, 963}, {48, 72, 1773},
						{67, 48, 1376}, {72, 67, 1213}, {67, 98, 325}, {67, 74, 757},
						{91, 67, 1147}, {40, 16, 1123}, {16, 60, 1210}, {66, 16, 837},
						{40, 66, 1919}, {66, 60, 2008}, {58, 76, 50}, {63, 58, 126}, {89, 76, 444},
						{58, 89, 474}, {65, 56, 1148}, {92, 65, 1350}, {25, 76, 397}, {17, 48, 1842},
						{67, 17, 466}, {63, 76, 149}, {28, 95, 664}, {24, 28, 2382}, {88, 10, 97},
						{18, 88, 438}, {10, 18, 358}, {18, 95, 940}, {17, 18, 1286}, {95, 88, 1288},
						{95, 10, 1191}, {0, 53, 167}, {76, 48, 2075}, {60, 76, 1837}, {48, 60, 921},
						{2, 23, 167}, {54, 2, 159}, {23, 54, 97}, {48, 63, 1954}, {13, 22, 353},
						{12, 13, 357}, {80, 36, 1651}, {60, 80, 581}, {67, 93, 554}, {12, 93, 520},
						{93, 9, 1250}, {60, 72, 1391}, {29, 34, 798}, {34, 56, 602}, {20, 94, 455},
						{94, 56, 308}, {94, 64, 192}, {64, 56, 220}, {90, 32, 1174}, {1, 90, 1895},
						{56, 60, 1331}, {32, 56, 1411}, {21, 95, 396}, {28, 21, 942}, {56, 94, 308},
						{28, 95, 664}, {96, 8, 927}, {46, 96, 1436}, {8, 46, 613}, {44, 60, 931},
						{48, 44, 897}, {60, 48, 921}, {30, 44, 1061}, {48, 30, 279}, {50, 28, 619},
						{81, 50, 825}, {28, 81, 259}, {48, 56, 1271}, {69, 95, 476}, {28, 69, 353},
						{99, 95, 228}, {69, 99, 649}, {75, 28, 183}, {50, 75, 664}, {75, 81, 166},
						{28, 75, 183}, {28, 52, 1733}, {75, 69, 268}, {86, 30, 450}, {48, 86, 654},
						{62, 9, 729}, {86, 62, 456}, {9, 86, 421}, {48, 28, 359}, {28, 69, 353},
						{69, 48, 466}, {41, 62, 301}, {9, 5, 367}, {5, 86, 58}, {44, 5, 1125},
						{9, 44, 1452}, {86, 44, 1068}, {49, 15, 436}, {15, 6, 477}, {34, 60, 1711},
						{15, 33, 625}, {29, 42, 399}, {60, 29, 934}, {42, 60, 793}, {60, 9, 603},
						{93, 60, 648}, {15, 93, 555}, {93, 33, 446}, {42, 78, 113}, {78, 29, 482},
						{48, 42, 950}, {13, 49, 418}, {49, 22, 222}, {57, 12, 802}, {93, 12, 520},
						{93, 37, 632}, {37, 12, 236}, {12, 49, 357}, {93, 12, 520}, {49, 93, 296},
						{80, 93, 1224}, {9, 67, 938}, {40, 60, 88}, {90, 56, 2189}, {0, 17, 1090},
						{11, 0, 489}, {51, 0, 304}, {0, 8, 768}, {8, 24, 921}, {96, 24, 33},
						{96, 28, 2396}, {28, 38, 1694}, {46, 28, 1392}, {38, 46, 567},
						{19, 52, 1872}, {85, 19, 201}, {52, 85, 1759}, {46, 28, 1392}, {19, 71, 259},
						{70, 19, 713}, {71, 70, 625}, {46, 97, 1078}, {19, 46, 89}, {97, 19, 1020},
						{71, 55, 862}, {55, 70, 366}, {46, 23, 767}, {23, 97, 862}, {47, 46, 1172},
						{51, 47, 1334}, {23, 71, 444}, {19, 23, 679}, {71, 70, 625}, {47, 23, 589},
						{83, 47, 1008}, {46, 43, 938}, {43, 8, 926}, {71, 45, 513}, {23, 71, 444},
						{45, 23, 290}, {19, 28, 1478}, {46, 19, 89}, {59, 46, 713}, {19, 59, 668},
						{19, 52, 1872}, {19, 73, 983}, {19, 71, 259}, {71, 73, 800}, {31, 71, 431},
						{23, 31, 203}, {2, 31, 345}, {54, 73, 947}, {71, 54, 524}, {54, 31, 191},
						{54, 84, 867}, {84, 73, 280}, {54, 73, 947}, {3, 87, 1838}, {82, 3, 1799},
						{87, 82, 105}, {73, 4, 1060}, {87, 73, 1096}, {4, 87, 324}, {61, 83, 229},
						{68, 61, 954}, {3, 68, 766}, {73, 47, 1435}, {47, 4, 600}, {54, 47, 494},
						{3, 4, 1789}, {4, 68, 1084}, {82, 4, 221}, {83, 4, 426}, {61, 4, 216}},
				29923);
	});
	it("can handle floating point weights", [&] {
		test<MatchingImpl>({{0, 1, 8.7}, {1, 2, 8.3}, {2, 3, 8.15}, {3, 0, 8.33}, {1, 3, 2.123}},
				16.63);
	});
	it("can be reused", [&] {
		Graph graph;
		completeGraph(graph, 2);
		EdgeArray<int> weights(graph, 1);
		MatchingImpl<int> m;
		call(m, graph, weights, 1);
		graph.clear();
		completeGraph(graph, 4);
		weights.init(graph, 2);
		call(m, graph, weights, 4);
	});
	it("can handle expands during the algorithm", [&] {
		test<MatchingImpl>({{0, 1, 1}, {1, 2, 0}, {2, 3, 0}, {1, 3, 0}, {1, 4, 1}, {4, 5, 10},
								   {5, 6, 1}, {6, 7, 0}, {7, 8, 0}, {6, 8, 0}, {6, 9, 1}},
				12);
	});
	it("calculates a maximum weight perfect matching correctly", [&] {
		Graph graph;
		EdgeArray<int> weights(graph);
		std::vector<std::tuple<int, int, int>> edges = {{0, 1, 3}, {0, 4, 5}, {0, 5, 4}, {1, 2, 7},
				{1, 5, 3}, {1, 6, 2}, {2, 3, 1}, {2, 6, 3}, {2, 7, 2}, {3, 7, 1}, {4, 5, 2},
				{5, 6, 1}, {6, 7, 4}};
		createGraph(graph, weights, edges);
		std::unordered_set<edge> matching;
		MatchingImpl<int> m;
		bool result = m.maximumWeightPerfectMatching(graph, weights, matching);
		AssertThat(result, IsTrue());
		AssertThat(isPerfectMatching(graph, matching), IsTrue());
		AssertThat(m.matchingWeight(matching, weights), Equals(14));
	});
	it("calculates a maximum weight matching correctly", [&] {
		Graph graph;
		EdgeArray<int> weights(graph);
		std::vector<std::tuple<int, int, int>> edges = {{0, 1, 3}, {0, 4, 5}, {0, 5, 4}, {1, 2, 7},
				{1, 5, 3}, {1, 6, 2}, {2, 3, 1}, {2, 6, 3}, {2, 7, 2}, {3, 7, 1}, {4, 5, 2},
				{5, 6, 1}, {6, 7, 4}};
		createGraph(graph, weights, edges);
		std::unordered_set<edge> matching;
		MatchingImpl<int> m;
		m.maximumWeightMatching(graph, weights, matching);
		AssertThat(isMatching(graph, matching), IsTrue());
		AssertThat(m.matchingWeight(matching, weights), Equals(16));
	});
}

go_bandit([] {
	describe("Blossom I", [&] { runAllTests<MatchingBlossom>(); });
	describe("Blossom V", [&] { runAllTests<MatchingBlossomV>(); });
});
