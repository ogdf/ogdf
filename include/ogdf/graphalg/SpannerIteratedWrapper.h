/** \file
 * \brief A wrapper class for iterating calls to spanner algorithms.
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

#include <ogdf/graphalg/SpannerModule.h>

#include <memory>

namespace ogdf {

/**
 * A implementation-independed wrapper class to execute a spanner algorithm multiple times.
 * The amount of maximum iterations can be set. If this is reached or the optional
 * timelimit is exceeded, the best solution so far is returned. The quality of the solution is
 * determined by the sparsity of the spanner, namely the spanner with the least amount of
 * edges regardless of whether the spanner is weighted or not.
 *
 * If no valid solution was found and the timelimit is exceeded, SpannerModule::ReturnType::TimeoutInfeasible
 * is returned. If all iterations are completed and no valid solution was found,
 * SpannerModule::ReturnType::NoFeasibleSolution is returned. If an error happens, it
 * is directly returned and the iterative execution aborted. If there is at least one valid
 * result, the return type will be SpannerModule::ReturnType::Feasible and the spanner
 * is valid.
 *
 * @tparam TWeight The type of weights to get from \p GA
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerIteratedWrapper : public SpannerModule<TWeight> {
public:
	/**
	 * Initializes the wrapper.
	 *
	 * @param module The algorithm implementation to use. It is internally put into a unique
	 * pointer, so this class will take care about the deletion of the algorithm.
	 * @param maxIterations The maximum amount of alls to the algorithm.
	 */
	SpannerIteratedWrapper(SpannerModule<TWeight>* module, int maxIterations)
		: m_module(module), m_maxIterations(maxIterations) { }

	//! @copydoc ogdf::SpannerModule::preconditionsOk
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		return m_module->preconditionsOk(GA, stretch, error);
	}

	/**
	 * @returns the amount of iterations during the last call of the wrapper.
	 */
	int getExecutedIterations() { return m_iterations; }

private:
	std::unique_ptr<SpannerModule<TWeight>> m_module;
	const int m_maxIterations;
	int m_iterations;

	//! @copydoc ogdf::SpannerModule::execute
	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		const Graph& G = m_GA->constGraph();
		if (isTimelimitEnabled()) {
			int64_t timeLeft = max(getTimeLeft(), static_cast<int64_t>(0));
			m_module->setTimelimit(timeLeft);
		}

		int bestSize = std::numeric_limits<int>::max();

		for (m_iterations = 1; m_iterations <= m_maxIterations; m_iterations++) {
			// No assertTimeLeft here: If there is an timeout, check, if we had a previous solution
			if (isTimelimitEnabled() && getTimeLeft() <= 0) {
				if (bestSize == std::numeric_limits<int>::max()) {
					return SpannerModule<TWeight>::ReturnType::TimeoutInfeasible;
				} else {
					return SpannerModule<TWeight>::ReturnType::Feasible;
				}
			}

			GraphCopySimple spanner(G);
			EdgeArray<bool> inSpanner(G);

			typename SpannerModule<TWeight>::ReturnType r =
					m_module->call(*m_GA, m_stretch, spanner, inSpanner);
			if (r == SpannerModule<TWeight>::ReturnType::NoFeasibleSolution) {
				continue;
			}
			if (r == SpannerModule<TWeight>::ReturnType::TimeoutInfeasible) {
				// do we found at least one solution?
				if (bestSize == std::numeric_limits<int>::max()) {
					return SpannerModule<TWeight>::ReturnType::TimeoutInfeasible;
				} else {
					return SpannerModule<TWeight>::ReturnType::Feasible;
				}
			}
			if (r != SpannerModule<TWeight>::ReturnType::Feasible) {
				return r; // return errors directly
			}

			if (spanner.numberOfEdges() < bestSize) {
				// found a new best solution
				bestSize = spanner.numberOfEdges();

				// Copy solution to the members
				m_spanner->clear();
				m_spanner->createEmpty(G);
				for (node n : G.nodes) {
					m_spanner->newNode(n);
				}
				m_inSpanner->init(G, false);

				for (edge e : spanner.edges) {
					edge eOrig = spanner.original(e);
					(*m_inSpanner)[eOrig] = true;
					m_spanner->newEdge(eOrig);
				}
			}
		}
		m_iterations--; // remove the overflow iteration, which is not executed

		if (bestSize == std::numeric_limits<int>::max()) {
			return SpannerModule<TWeight>::ReturnType::NoFeasibleSolution;
		} else {
			return SpannerModule<TWeight>::ReturnType::Feasible;
		}
	}

	using SpannerModule<TWeight>::getTimeLeft;
	using SpannerModule<TWeight>::isTimelimitEnabled;
	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_stretch;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

}
