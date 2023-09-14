/** \file
 * \brief Implements simple maximum density subgraph algorithm.
 *
 * \author Finn Stutzenstein
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

#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/Stopwatch.h>
#include <ogdf/graphalg/MaximumDensitySubgraph.h>
#include <ogdf/graphalg/MinSTCutMaxFlow.h>

namespace ogdf {

bool maximumDensitySubgraph(Graph& G, NodeSet<true>& subgraphNodes,
		std::function<node(node)> resultNodeMap, int64_t timelimit) {
	StopwatchCPU watch;
	watch.start();

	subgraphNodes.clear();
	int n = G.numberOfNodes();
	int m = G.numberOfEdges();

	// special cases
	if (m == 0) {
		return true;
	}
	if (m == 1) {
		edge e = G.firstEdge();
		subgraphNodes.insert(resultNodeMap(e->source()));
		subgraphNodes.insert(resultNodeMap(e->target()));
		return true;
	}

	double l = 0;
	double u = G.numberOfEdges();
	EpsilonTest eps;

	// save degrees, so the source and target (added below) do not interfere.
	NodeArray<int> degree(G);
	for (node v : G.nodes) {
		degree[v] = v->degree();
	}

	// construct network
	node s = G.newNode();
	node t = G.newNode();
	EdgeArray<double> weights(G);
	List<edge> targetEdges;
	for (edge e : G.edges) {
		weights[e] = 1.0;
	}
	for (node v : G.nodes) {
		if (v != s && v != t) {
			edge sourceEdge = G.newEdge(s, v);
			weights[sourceEdge] = static_cast<double>(m);
			targetEdges.pushBack(G.newEdge(v, t));
		}
	}

	const double limit = 1.0 / (n * (n - 1));
	while (eps.geq((u - l), limit)) {
		double g = (u + l) / 2;

		// check timelimit
		if (timelimit >= 0 && watch.milliSeconds() >= timelimit) {
			return false;
		}

		// set network weights for target edges (all other edges have constant (preset) weights)
		for (edge e : targetEdges) {
			weights[e] = static_cast<double>(m) + 2 * g - degree[e->source()];
			OGDF_ASSERT(weights[e] >= 0);
		}

		MinSTCutMaxFlow<double> mstc(true, new MaxFlowGoldbergTarjan<double>(), true, false,
				new EpsilonTest());

		List<edge> cutEdges;
		mstc.call(G, weights, s, t, cutEdges);

		List<node> sPartition;
		for (node v : G.nodes) {
			if (mstc.isInFrontCut(v)) {
				sPartition.pushBack(v);
			}
		}

#ifdef OGDF_DEBUG
		OGDF_ASSERT(sPartition.size() >= 1);
		bool sInSPartition = false;
		for (node v : sPartition) {
			OGDF_ASSERT(v != t);
			if (v == s) {
				sInSPartition = true;
			}
		}
		OGDF_ASSERT(sInSPartition);
#endif

		if (sPartition.size() == 1) {
			// only source is in A
			u = g;
		} else {
			l = g;

			subgraphNodes.clear();
			for (node v : sPartition) {
				if (v != s) {
					subgraphNodes.insert(resultNodeMap(v));
				}
			}
		}
	}

	return true;
}

}
