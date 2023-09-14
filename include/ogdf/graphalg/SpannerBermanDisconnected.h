/** \file
 * \brief Implementation of an k-spanner approximation algorithm from
 * Berman et al.
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

#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/SpannerBerman.h>
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Wrapper around SpannerBerman: For each component of the
 * graph, the algorithm will be called. This allows to use the
 * algorithm on graphs with more than one component.
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerBermanDisconnected : public SpannerModule<TWeight> {
public:
	//! @copydoc ogdf::SpannerModule::preconditionsOk
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		if (!isSimple(GA.constGraph())) {
			error = "The graph is not simple";
			return false;
		}
		if (stretch < 1.0) {
			error = "The stretch must be >= 1.0";
			return false;
		}
		return true;
	}

private:
	//! @copydoc ogdf::SpannerModule::execute
	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		const Graph& G = m_GA->constGraph();

		if (G.numberOfNodes() == 0) {
			return SpannerModule<TWeight>::ReturnType::Feasible;
		}

		NodeArray<int> components(G);
		int amountComponenets = connectedComponents(G, components);

		SpannerBerman<TWeight> spannerBerman;

		for (int c = 0; c < amountComponenets; c++) {
			assertTimeLeft();

			List<node> nodes;
			for (node n : G.nodes) {
				if (components[n] == c) {
					nodes.pushBack(n);
				}
			}

			GraphCopySimple GC;
			GC.createEmpty(G);
			inducedSubGraph(G, nodes.begin(), GC);

			GraphCopySimple spanner(GC);
			EdgeArray<bool> inSpanner;
			GraphAttributes GCA(GC, 0);
			GCA.directed() = m_GA->directed();
			if (m_GA->has(GraphAttributes::edgeIntWeight)) {
				GCA.addAttributes(GraphAttributes::edgeIntWeight);
				for (edge e : GC.edges) {
					GCA.intWeight(e) = m_GA->intWeight(GC.original(e));
				}
			}
			if (m_GA->has(GraphAttributes::edgeDoubleWeight)) {
				GCA.addAttributes(GraphAttributes::edgeDoubleWeight);
				for (edge e : GC.edges) {
					GCA.doubleWeight(e) = m_GA->doubleWeight(GC.original(e));
				}
			}

			if (isTimelimitEnabled()) {
				int64_t timeLeft = max(getTimeLeft(), static_cast<int64_t>(0));
				spannerBerman.setTimelimit(timeLeft);
			}
			typename SpannerModule<TWeight>::ReturnType r =
					spannerBerman.call(GCA, m_stretch, spanner, inSpanner);
			if (r != SpannerModule<TWeight>::ReturnType::Feasible) {
				return r;
			}

			// copy used edges to m_spanner and m_inSpanner
			for (edge e : spanner.edges) {
				edge eOrig = GC.original(spanner.original(e));
				(*m_inSpanner)[eOrig] = true;
				m_spanner->newEdge(eOrig);
			}
		}

		return SpannerModule<TWeight>::ReturnType::Feasible;
	}

	using SpannerModule<TWeight>::getWeight;
	using SpannerModule<TWeight>::assertTimeLeft;
	using SpannerModule<TWeight>::getTimeLeft;
	using SpannerModule<TWeight>::isTimelimitEnabled;
	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_stretch;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

}
