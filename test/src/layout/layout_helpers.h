/** \file
 * \brief Test helpers for layout algorithms
 *
 * \author Carsten Gutwenger, Tilo Wiedera
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

#include <set>
#include <random>
#include <bandit/bandit.h>
#include <resources.h>

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/module/LayoutModule.h>
#include <ogdf/module/GridLayoutModule.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>

namespace ogdf {

constexpr int MAX_NODES = 200;
constexpr int MIN_NODES = 25;
constexpr int STEP_SIZE = 50;

enum class GraphRequirement {
	planar,
	triconnected,
	connected
};

inline void insertGraph(Graph &g, const Graph &g2)
{
	NodeArray<node> map(g2);

	for(node v : g2.nodes)
		map[v] = g.newNode();

	for(edge e: g2.edges)
		g.newEdge(map[e->source()], map[e->target()]);
}

inline void createDisconnectedGraph(
	Graph &G,
	int n_max,
	double density_min,
	double density_max,
	int cc)
{
	int n = randomNumber(3,n_max);
	int m = randomNumber((int)(density_min*n),(int)(density_max*n));
	int b = n/25+1;
	planarCNBGraph(G, n/b+1, m/b+1, b);

	for(int i = cc-1; i > 0; --i) {
		Graph g2;
		n = randomNumber(3,n);
		m = randomNumber((int)(density_min*n),(int)(density_max*n));
		b = n/25+1;
		planarCNBGraph(g2, n/b+1, m/b+1, b);
		insertGraph(G,g2);
	}
}

inline void createAlmostPlanarGraph(Graph &G, int n, int m, int add_m)
{
	planarBiconnectedGraph(G, n, m);

	Array<node> table(n);

	int i = 0;
	for(node v : G.nodes) {
		OGDF_ASSERT(i < n);
		table[i++] = v;
	}

	for(i = 0; i < add_m; ++i) {
		G.newEdge(table[randomNumber(0, n-1)], table[randomNumber(0, n-1)]);
	}

	makeSimpleUndirected(G);
}

inline void getRandomLayout(GraphAttributes &GA)
{
	const Graph &G = GA.constGraph();
	double max_x = 2.0 * sqrt(G.numberOfNodes());
	double max_y = max_x;

	std::minstd_rand rng(randomSeed());
	std::uniform_real_distribution<> rand_x(0.0,max_x);
	std::uniform_real_distribution<> rand_y(0.0,max_y);

	for(node v : G.nodes) {
		GA.x(v) = rand_x(rng);
		GA.y(v) = rand_y(rng);
	}
}

inline int64_t callLayout(const Graph &G, LayoutModule &L, bool isGridLayout, long extraAttributes)
{
	int64_t T;

	if(isGridLayout) {
		GridLayout gl;
		System::usedRealTime(T);
		((GridLayoutModule*)&L)->callGrid(G, gl);
	}
	else {
		extraAttributes |= GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics;
		GraphAttributes GA(G, extraAttributes);
		getRandomLayout(GA);
		System::usedRealTime(T);
		L.call(GA);
	}
	return System::usedRealTime(T);
}

inline bool doesNotInclude(const std::set<GraphRequirement> &req, std::initializer_list<GraphRequirement> values) {
	for(auto v : values) {
		if(req.find(v) != req.end()) {
			return false;
		}
	}

	return true;
}

/**
 * Runs several tests for a given layout module.
 * The layout algorithm is executed for different graphs.
 * There are no assertions yet.
 *
 * \param name
 * 	the name to be used for describing this module
 * \param L
 * 	the module to be tested
 * \param extraAttributes
 *  init attributes for GraphAttributes
 * \param req
 * 	the requirements for graphs to be drawn, see GraphRequirementFlags for details
 * \param maxNodes
 * 	the maximum number of nodes the algorithm should be tested on
 * \param isGridLayout
 * 	set this to true if the layout module is a grid layout module
 * \param skipTreeWithProbablyNegativeCoordinates
 *  set this to true if the file negative_coordinates_tree.gml should be skipped
 * \param skipMe
 *  set this to true to skip the entire test
 */
inline void describeLayoutModule(
  const std::string name,
  LayoutModule &L,
  long extraAttributes = 0,
  std::set<GraphRequirement> req = {},
  int maxNodes = MAX_NODES,
  bool isGridLayout = false,
  bool skipMe = false)
{
	maxNodes = max(maxNodes, MIN_NODES+1);
	int steps = static_cast<int>(ceil((maxNodes-MIN_NODES) / static_cast<double>(STEP_SIZE)));
	auto performTest = [&](){
		auto call = [&](const Graph &G) {
			return callLayout(G, L, isGridLayout, extraAttributes);
		};

		if(doesNotInclude(req, {GraphRequirement::triconnected})) {
			bandit::it("works on trees", [&](){
				int64_t time = 0;
				for(int n = MIN_NODES; n < maxNodes; n += STEP_SIZE) {
					Graph G;
					randomTree(G, n);
					time += call(G);
				}
				cout << endl << "      average time was " << time/steps << "ms" << endl;
			});

			for_each_graph_it("works on a tree with probably negative coordinates",
			                  { "misc/negative_coordinates_tree.gml" },
			                  [&](Graph &G, const string &filename) { call(G); });

			bandit::it("works on planar connected graphs", [&](){
				int64_t time = 0;
				for(int n = MIN_NODES; n < maxNodes; n += STEP_SIZE) {
					Graph G;
					planarConnectedGraph(G, n, randomNumber(n, 2*n));
					makeSimpleUndirected(G);
					time += call(G);
				}
				cout << endl << "      average time was " << time/steps << "ms" << endl;
			});

			bandit::it("works on planar biconnected graphs", [&](){
				int64_t time = 0;
				for(int n = MIN_NODES; n < maxNodes; n += STEP_SIZE) {
					Graph G;
					planarBiconnectedGraph(G, n, randomNumber(3*n/2, 2*n));
					time += call(G);
				}
				cout << endl << "      average time was " << time/steps << "ms" << endl;
			});
		}

		bandit::it("works on planar triconnected graphs", [&](){
			int64_t time = 0;
			for(int n = MIN_NODES; n < maxNodes; n += STEP_SIZE) {
				Graph G;
				planarTriconnectedGraph(G, n, .5, .5);
				time += call(G);
			}
			cout << endl << "      average time was " << time/steps << "ms" << endl;
		});

		if(doesNotInclude(req, {GraphRequirement::planar, GraphRequirement::triconnected})) {
			bandit::it("works on almost planar graphs", [&](){
				int64_t time = 0;
				for(int n = MIN_NODES; n < MAX_NODES; n += STEP_SIZE) {
					Graph G;
					createAlmostPlanarGraph(G, n, 2*n, 10);
					time += call(G);
				}
				cout << endl << "      average time was " << time/steps << "ms" << endl;
			});

			bandit::it("works on biconnected graphs", [&](){
				int64_t time = 0;
				for(int n = MIN_NODES; n < maxNodes; n += STEP_SIZE) {
					Graph G;
					randomBiconnectedGraph(G, n, n*(n-1)/6);
					makeSimpleUndirected(G);
					time += call(G);
				}
				cout << endl << "      average time was " << time/steps << "ms" << endl;
			});

			if(doesNotInclude(req, {GraphRequirement::connected})) {
				bandit::it("works on disconnected graphs", [&](){
					int64_t time = 0;
					for(int n = MIN_NODES; n < maxNodes; n += STEP_SIZE) {
						Graph G;
						createDisconnectedGraph(G, n/7, 1.4, 2.6, 7);
						time += call(G);
					}
					cout << endl << "      average time was " << time/steps << "ms" << endl;
				});
			}
		}
	};

	if(skipMe) {
		bandit::describe_skip(name, performTest);
	} else {
		bandit::describe(name, performTest);
	}
}

/**
 * Runs several tests for a given layout module.
 * The layout algorithm is executed for different graphs.
 * There are no assertions yet.
 *
 * \param name
 * 	the name to be used for describing this module
 * \param L
 * 	the module to be tested
 * \param req
 * 	the requirements for graphs to be drawn, see GraphRequirementFlags for details
 * \param maxNodes
 * 	the maximum number of nodes the algorithm should be tested on
 * \param isGridLayout
 * 	set this to true if the layout module is a grid layout module
 **/
inline void describeGridLayoutModule(
  const std::string name,
  GridLayoutModule &L,
  std::initializer_list<GraphRequirement> req = {},
  int maxNodes = MAX_NODES,
  bool skipMe = false)
{
	describeLayoutModule(name, L, 0, req, maxNodes, true, skipMe);
}

}
