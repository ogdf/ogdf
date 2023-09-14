/** \file
 * \brief Implementation of the random k-spanner algorithm
 * of Elkin and Neiman 2018.
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

#pragma once

#include <ogdf/basic/Queue.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/SpannerIteratedWrapper.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Randomized multiplicative spanner calculation by propagating random messages
 * through the graph.
 *
 * M. Elkin und O. Neiman. Efficient Algorithms for Constructing Very Sparse Spanners
 * and Emulators. ACM Trans. Algorithms 15.1 (Nov. 2018). doi: 10.1145/3274651.
 *
 * Conditions for the graph:
 * - undirected
 * - unweighted
 *
 * The stretch \f$s\f$ must satisfy: \f$s\in\{1,3,5,7,\ldots\}\f$.
 *
 * The preconditions can be checked with SpannerElkinNeiman::preconditionsOk.
 *
 * Calculates a \f$(2k-1)\f$-spanner with size \f$\mathcal{O}(n^{1+1/k}/\epsilon)\f$,
 * \f$0<\epsilon<1\f$. The runtime is \f$\mathcal{O}(m)\f$ with a success probability
 * of at least \f$1-\epsilon\f$.
 *
 * Note that in practice, \f$\epsilon\f$ can be larger than 1 which can result in even better
 * results, but has a larger probability of failing. Remember, that the success probability
 * is a lower bound, so larger \f$\epsilon\f$ does not imply that it is impossible to generate
 * a solution.
 *
 * @note
 * Implementation details:
 *
 * Since the paper does not state much detail about the algorithm, some things must be explained
 * here that are not trivial:
 * - "Each vertex x that received a message originated at u, stores [...] a neighbor p_u(x)
 * 	 that lies on a shortest path from x to u": This also includes, that the message must be send
 *   via a shortest path. In the implementation a breadth first search is used since the graph is
 *   unweighted (no need for dijkstra). With marking visited notes it is ensured that when visiting
 *   a node the first time, it *is* a shortest path.
 * - Continuing: "(this neighbor sent x the message from u, breaking ties arbitrarily if there is
 *   more than one)". To not make the code more complicated to search for all shortest paths, saving
 *   them into a list and in a second step choosing one random element from the list, just *the first*
 *   shortest path is taken (see above). It works and the results do not seem to be much worse.
 * - There are two feasibility checks below with the second one modified. See the comment at the very
 *   end of the execute function.
 * - The paper describes the algorithm, that u broadcasts a value to all nodes within distance k. Then,
 *   each receiving node x is handled (e.g. building m(x) and C(x) only with respect to x). If it is
 *   implemented this way, there has to be a n*n matrix for m_u(x) and p_u(x) which is very memory consuming.
 *   Tests with some large instances (n about 40000) resulted in out of memory issues during the
 *   initialization of the arrays.
 *   To lower the memory footprint to just O(n), the logic is switched around: The algorithm always fixes
 *   a receiving node x, which is fine for the later parts. So now all nodes u are searched, that sends
 *   a message to x. This can also be done by the same BFS with little modification. In summary, this allows
 *   for just storing the value and edge for each node, not a matrix of nodes.
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerElkinNeiman : public SpannerModule<TWeight> {
	//! The original graph
	const Graph* m_G;

	EpsilonTest m_eps;

	double m_c; //!< the parameter for the exponential distribution. It depends on epsilon
	double m_beta; //!< The parameter for the exponential distribution
	int m_k; //!< the parameter k derived from the stretch
	int m_intStretch; //!< m_stretch, but as an integer

	/**
	 * The default value for epsilon.
	 * This is a result of some experimental tests on arbitrary graphs: This epsilon provides a good compromise of
	 * - Around 75% chance of success
	 * - low mean value of the resulting spanner sizes
	 * - medium statistical dispersion -> The results mainly disperse evenly, so a high dispersion
	 *   leads to a higher probability of finding better minimum spanners
	 */
	const double DEFAULT_EPSILON = 0.8;

public:
	SpannerElkinNeiman() : SpannerModule<TWeight>() { setEpsilon(DEFAULT_EPSILON); }

	SpannerElkinNeiman(double epsilon) : SpannerModule<TWeight>() { setEpsilon(epsilon); }

	/**
	 * Sets the epsilon for this algorithm. Note that it must be >0 and
	 * bounded with <1 in the paper. The upper bound can be exceeded, so values like 2 or 3
	 * can also be used.
	 */
	void setEpsilon(double epsilon) {
		OGDF_ASSERT(epsilon > 0);
		m_c = 3.0 / epsilon; // see corollary 8.
	}

	//! @copydoc ogdf::SpannerModule::preconditionsOk
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		if (GA.directed()) {
			error = "The graph must be undirected";
			return false;
		}
		double integralPart;
		if (std::modf(stretch, &integralPart) != 0.0) {
			error = "The stretch is required to be an integer, not " + to_string(m_stretch);
			return false;
		}
		int intStretch = static_cast<int>(stretch);
		if (intStretch < 1) {
			error = "The stretch must be >= 1.0";
			return false;
		}
		if (intStretch % 2 == 0) {
			error = "The stretch must be odd";
			return false;
		}
		if (GA.has(GraphAttributes::edgeDoubleWeight) || GA.has(GraphAttributes::edgeIntWeight)) {
			error = "The graph must be unweighted";
			return false;
		}
		return true;
	}

private:
	//! @copydoc ogdf::SpannerModule::init
	virtual void init(const GraphAttributes& GA, double stretch, GraphCopySimple& spanner,
			EdgeArray<bool>& inSpanner) override {
		SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
		m_intStretch = static_cast<int>(stretch);
		m_k = (m_intStretch + 1) / 2;
		m_G = &GA.constGraph();
		m_beta = log(m_c * m_G->numberOfNodes()) / m_k; // log is logarithm naturale
	}

	struct BfsNode {
		node n; //!< The current node to visit all neighbors from
		int depth; //!< The depth of the node n. It is equal to the distance from the source node
		edge p; //!< The first edge of the path

		BfsNode(node _n, double _depth, edge _p) : n(_n), depth(_depth), p(_p) { }
	};

	//! @copydoc ogdf::SpannerModule::execute
	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		NodeArray<double> r(*m_G);
		// First, assign each node the random variable. This is done here, so the first
		// feasibility check can be done fast (This is the one that mostly fails).
		for (node u : m_G->nodes) {
			r[u] = randomDoubleExponential(m_beta);
			// Feasibility check 1: See proof of Theorem 1 and "Standard Centralized Model" in 2.1.1
			if (m_eps.geq(r[u], static_cast<double>(m_k))) {
				return SpannerModule<TWeight>::ReturnType::NoFeasibleSolution;
			}
		}

		// Paper: u broadcasts a message x within distance k
		// Here: Fix a receiving node x and find all nodes u within distance k, that send messages to x
		for (node x : m_G->nodes) {
			// values[u]: x recieves the value from u (m_u(x))
			NodeArray<double> values(*m_G, std::numeric_limits<double>::lowest());

			// firstEdge[u]: The first edge from the path x->u (p_u(x))
			// This is the last edge from the path u->x as described in the paper
			NodeArray<edge> firstEdge(*m_G, nullptr);

			QueuePure<BfsNode> queue;
			NodeArray<bool> marked(*m_G, false); // Just visit a node once.

			queue.emplace(x, 0, nullptr); // p==nullptr: we do not have a first edge
			marked[x] = true; // Do not revisit x
			while (!queue.empty()) {
				BfsNode bfsNode = queue.pop();
				int distance = bfsNode.depth + 1;

				for (adjEntry adj : bfsNode.n->adjEntries) {
					node u = adj->twinNode();
					if (marked[u]) {
						continue;
					}

					edge p = bfsNode.p;
					if (p == nullptr) {
						p = adj->theEdge(); // Set the first edge: Only happens in the first expansion step of the bfs
					}

					values[u] = r[u] - distance;
					firstEdge[u] = p;

					if (distance < m_k) { // distance is <= m_k when a dfs node is popped from the stack.
						queue.emplace(u, distance, p);
						marked[u] = true;
					}
				}
			}

			assertTimeLeft();

			// find edges that must be added to the spanner:
			// m(x) = max{m_u(x) | u in V}
			double maxValue = std::numeric_limits<double>::lowest();
			for (node u : m_G->nodes) {
				Math::updateMax(maxValue, values[u]);
			}
			maxValue -= 1.0; // m(x)-1

			assertTimeLeft();

			// Add the edges (sets C(x) in the paper) to the spanner.
			for (node u : m_G->nodes) {
				if (values[u] != std::numeric_limits<double>::lowest()
						&& m_eps.geq(values[u], maxValue)) {
					edge e = firstEdge[u];
					if (!(*m_inSpanner)[e]) {
						m_spanner->newEdge(e);
						(*m_inSpanner)[e] = true;
					}
				}
			}

			assertTimeLeft();
		}

		// Feasibility check 2: See proof of Theorem 1 and "Standard Centralized Model" in 2.1.1
		// Note: The paper states, that the spanner must have >= n-1 edges. This does not work, for
		// all graphs, e.g. a forest as an input. The paper never states, that the graph must be connected
		// nor that it can be disconnected. Using n-<amount components> does work and solves the issue above.
		int components = connectedComponents(*m_G);
		if (m_spanner->numberOfEdges() < m_G->numberOfNodes() - components) {
			return SpannerModule<TWeight>::ReturnType::NoFeasibleSolution;
		}

		return SpannerModule<TWeight>::ReturnType::Feasible;
	}

	using SpannerModule<TWeight>::assertTimeLeft;
	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_stretch;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

/**
 * Use the ogdf::SpannerIteratedWrapper to execute the ogdf::SpannerElkinNeiman algorithm
 * up to 200 times.
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerElkinNeimanIterated : public SpannerIteratedWrapper<TWeight> {
public:
	SpannerElkinNeimanIterated()
		: SpannerIteratedWrapper<TWeight>(new SpannerElkinNeiman<TWeight>(), 200) {};
};

}
