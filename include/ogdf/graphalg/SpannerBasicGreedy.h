/** \file
 * \brief Implementation of the basic greedy (2k-1)-spanner
 * algorithm of Althöfer et al. 2007.
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
#include <ogdf/graphalg/SpannerModule.h>

namespace ogdf {

/**
 * Multiplicative spanner by greedily adding edges.
 *
 * I. Althöfer, G. Das, D. Dobkin und D. Joseph. On Sparse Spanners of Weighted Graphs.
 * Discrete Comput Geom 9 (1993), S. 81–100. doi: https://doi.org/10.1007/BF02189308.
 *
 * Conditions for the graph:
 * - undirected
 * - weighted
 *
 * The stretch \f$s\f$ must satisfy: \f$s\geq1\in\mathbb{R}\f$.
 *
 * The preconditions can be checked with SpannerBasicGreedy::preconditionsOk.
 *
 * Calculates a \f$(2k-1)\f$-spanner with size \f$<n^{(1+1/k)}\f$ (\f$\mathcal{O}(n^{1+1/k})\f$)
 * and lightness \f$<1+n/k\f$ (\f$\mathcal{O}(n/k)\f$). The runtime is
 * \f$\mathcal{O}(mn^{1+1/k}\log n)\f$
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerBasicGreedy : public SpannerModule<TWeight> {
	EdgeArray<TWeight> m_spannerWeights;
	EpsilonTest m_eps;

public:
	//! @copydoc ogdf::SpannerModule::preconditionsOk
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		if (GA.directed()) {
			error = "The graph must be undirected";
			return false;
		}
		if (stretch < 1.0) {
			error = "The stretch must be >= 1.0";
			return false;
		}
		return true;
	}

private:
	virtual void init(const GraphAttributes& GA, double stretch, GraphCopySimple& spanner,
			EdgeArray<bool>& inSpanner) override {
		SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);
		m_spannerWeights.init(spanner);
	}

	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		Array<edge> edges(m_GA->constGraph().numberOfEdges());
		int i = 0;
		for (edge e : m_GA->constGraph().edges) {
			edges[i++] = e;
		}
		edges.quicksort(EdgeWeightComparator(*m_GA));

		assertTimeLeft();

		for (edge e : edges) {
			node u = m_spanner->copy(e->source());
			node v = m_spanner->copy(e->target());
			double maxDistance = m_stretch * weight(e);
			double currentSpannerDistance = distanceInSpanner(u, v, maxDistance);
			if (m_eps.greater(currentSpannerDistance, maxDistance)) {
				edge newEdge = m_spanner->newEdge(e);
				m_spannerWeights[newEdge] = weight(e);
				(*m_inSpanner)[e] = true;
			}
			assertTimeLeft();
		}

		return SpannerModule<TWeight>::ReturnType::Feasible;
	}

	OGDF_EXPORT double distanceInSpanner(node s, node t, double maxLookupDist);

	/**
	 * \returns the weights of an edge \p e from m_G
	 */
	inline TWeight weight(edge e) { return getWeight(*m_GA, e); }

	struct EdgeWeightComparator {
		EdgeWeightComparator(const GraphAttributes& GA) : m_GA(GA) {};
		EdgeWeightComparator() = delete;
		EdgeWeightComparator(const EdgeWeightComparator&) = delete;
		EdgeWeightComparator& operator=(const EdgeWeightComparator&) = delete;

		bool less(edge a, edge b) const {
			return m_eps.less(getWeight(m_GA, a), getWeight(m_GA, b));
		}

	private:
		const GraphAttributes& m_GA;
		const EpsilonTest m_eps;
	};

	using SpannerModule<TWeight>::getWeight;
	using SpannerModule<TWeight>::assertTimeLeft;
	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_stretch;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

}
