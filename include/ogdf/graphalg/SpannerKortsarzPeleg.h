/** \file
 * \brief Implementation of the 2-spanner approximation algorithm of
 * Kortsarz and Peleg 1994.
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

#include <ogdf/basic/NodeSet.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/MaximumDensitySubgraph.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Approximation multiplicative 2-spanner calculation.
 *
 * G. Kortsarz und D. Peleg. Generating Sparse 2-Spanners. Journal of Algorithms 17.2
 * (1994), S. 222–236. doi: https://doi.org/10.1006/jagm.1994.1032.
 *
 * Conditions for the graph:
 * - simple
 * - undirected
 * - unweighted
 *
 * The stretch \f$s\f$ must satisfy: \f$s=2\f$.
 *
 * The preconditions can be checked with SpannerKortsarzPeleg::preconditionsOk.
 *
 * Calculates a 2-spanner with an approximation ratio \f$\mathcal{O}(|E|/|V|)\f$. The
 * runtime is \f$\mathcal{O}(m^2n^2\log(n^2/m))\f$ (given in the paper).
 *
 * The implemented runtime can be bound by \f$\mathcal{O}(nm\cdot MDS(n, m))\f$ with
 * \f$MDS(n, m)\f$ being the runtime for the maximumDensitySubgraph implementation.
 * The used implementation has a runtime of \f$\mathcal{O}(mn^2\log n)\f$ so this
 * implementation has a runtime of \f$\mathcal{O}(m^2n^3\log n)\f$.
 *
 * Note that the practical runtime is not as high since the MDS-problem is only called
 * for neighbor-graphs of a vertex. For sparse graphs, these subgraphs typically are not
 * as huge as the input graph. But remember: For complete graphs, the running time will
 * approach the upper bound given in the O-notation.
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerKortsarzPeleg : public SpannerModule<TWeight> {
	//! The copy of GA.constGraph() to work on. Edges will be removed in each iteration.
	GraphCopySimple m_G;
	EpsilonTest m_eps;

public:
	//! @copydoc ogdf::SpannerModule::preconditionsOk
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		if (GA.directed()) {
			error = "The graph must be undirected";
			return false;
		}
		if (stretch != 2.0) {
			error = "The stretch must be 2.0";
			return false;
		}
		if (!isSimple(GA.constGraph())) {
			error = "The graph is not simple";
			return false;
		}
		if (GA.has(GraphAttributes::edgeDoubleWeight) || GA.has(GraphAttributes::edgeIntWeight)) {
			error = "The graph must be unweighted";
			return false;
		}
		return true;
	}

private:
	virtual void init(const GraphAttributes& GA, double stretch, GraphCopySimple& spanner,
			EdgeArray<bool>& inSpanner) override {
		SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
		m_G.clear();
		m_G.init(GA.constGraph());
	}

	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		while (true) {
			node maxNode = nullptr;
			double maxDensity = 0.0;
			NodeSet<true> maxDenseSubset(m_GA->constGraph());
			List<edge> maxE_U; // holds the edges of the maximum dense subgraph

			for (node v : m_G.nodes) {
				// Create neighbor graph
				GraphCopySimple neighborGraph;
				// Note: v is not part of the neighbor graph.
				neighborGraph.createEmpty(m_G);
				for (adjEntry adj : v->adjEntries) {
					neighborGraph.newNode(adj->twinNode());
				}
				for (edge e : m_G.edges) {
					if (neighborGraph.copy(e->target()) != nullptr
							&& neighborGraph.copy(e->source()) != nullptr) {
						neighborGraph.newEdge(e);
					}
				}

				// Calculate dense subgraph and find all edges E_U of this subgraph
				NodeSet<true> denseSubset(m_G);
				int64_t timelimit = -1;
				if (isTimelimitEnabled()) {
					timelimit = max(static_cast<int64_t>(0), getTimeLeft());
				}
				bool success = maximumDensitySubgraph(
						neighborGraph, denseSubset,
						[&](node n) {
							OGDF_ASSERT(neighborGraph.original(n) != nullptr);
							return neighborGraph.original(n);
						},
						timelimit);
				if (!success) {
					throw typename SpannerModule<TWeight>::TimeoutException();
				}

				List<edge> E_U;
				for (edge e : m_G.edges) {
					if (denseSubset.isMember(e->source()) && denseSubset.isMember(e->target())) {
						E_U.pushBack(e);
					}
				}

				// Did we found a new maximum dense subgraph?
				double density = denseSubset.size() == 0
						? 0.0
						: static_cast<double>(E_U.size()) / denseSubset.size();
				if (m_eps.greater(density, maxDensity)) {
					maxDensity = density;
					maxNode = v;
					maxDenseSubset = denseSubset;
					maxE_U = E_U;
				}
			}

			if (m_eps.less(maxDensity, 1.0)) {
				break;
			}

			OGDF_ASSERT(maxNode != nullptr);

			// add edges to spanner
			for (node u : maxDenseSubset.nodes()) {
				edge e = m_G.searchEdge(maxNode, u); // e is part of m_G
				addToSpanner(e);

				// remove from m_G
				m_G.delEdge(e);
			}

			// Remove R(H^u, maxDenseSubset) from m_G
			// Improvement: H^c (see paper) is actually not needed. Since
			// R(H^u, maxDenseSubset) is the intersection of E_U and the current spanner
			// edges, but E_U only consists of current spanner edges (Note that the loop
			// above only contains edges not in E_U since maxNode is not part of the
			// neighbor graph), the intersection is not needed. One can directly delete
			// all edges in E_U.
			for (edge e : maxE_U) {
				m_G.delEdge(e);
			}
		}

		// Add uncovered edges to spanner
		for (edge e : m_G.edges) {
			addToSpanner(e);
		}

		return SpannerModule<TWeight>::ReturnType::Feasible;
	}

	/**
	 * Adds \p e from \e m_G to the spanner and sets inSpanner
	 */
	inline void addToSpanner(edge e) {
		m_spanner->newEdge(m_G.original(e));
		(*m_inSpanner)[m_G.original(e)] = true;
	}

	using SpannerModule<TWeight>::assertTimeLeft;
	using SpannerModule<TWeight>::isTimelimitEnabled;
	using SpannerModule<TWeight>::getTimeLeft;
	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

}
