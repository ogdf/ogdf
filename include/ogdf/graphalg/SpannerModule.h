/** \file
 * \brief Basic module for spanner algorithms
 *
 * Includes some useful functions dealing with min-cost flow
 * (problem generator, problem checker).
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

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/Stopwatch.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/graphalg/ShortestPathAlgorithms.h>

namespace ogdf {

/**
 * \brief Interface for spanner algorithms.
 *
 * All spanner implementations must:
 * - implement SpannerModule::execute and SpannerModule::preconditionsOk
 * - be timelimited: Call SpannerModule::assertTimeLeft e.g. in the main loop
 * - create a spanner (m_spanner) and a set of used edges regarding the original graph (m_inSpanner)
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerModule : public Module {
public:
	//! Initializes a spanner module.
	SpannerModule() : m_GA(nullptr) { }

	// destruction
	virtual ~SpannerModule() { }

	/**
	 * \brief Executes the algorithm.
	 *
	 * Example:
	 * \code
	 * GraphAttributes GA = ....
	 * int stretch = 2;
	 * GraphCopySimple spanner;
	 * EdgeArray<bool> inSpanner;
	 * std::unique_ptr<SpannerModule<int>> sm = std::make_unique<BaswanaSenClusterSpanner<int>>();
	 * typename SpannerModule<int>::Result result = sm->call(GA, stretch, spanner, inSpanner);
	 * \endcode
	 *
	 * @tparam TWeight the type of weights to get from \p GA
	 *
	 * @param GA The graph used. GraphAttributes can be used to pass in weights
	 * @param stretch The multiplicative shortest-path-length factor
	 * @param spanner The resulting spanner. Must be a copy of GA.constGraph()
	 * @param inSpanner Maps true to an edge iff the edge is in the spanner
	 * @return The state in which the algorithm returns. Only for OK the spanner is valid.
	 */
	virtual ReturnType call(const GraphAttributes& GA, double stretch, GraphCopySimple& spanner,
			EdgeArray<bool>& inSpanner) {
		init(GA, stretch, spanner, inSpanner);
#ifdef OGDF_DEBUG
		std::string error;
		OGDF_ASSERT(preconditionsOk(GA, stretch, error));
#endif

		if (isTimelimitEnabled()) {
			m_watch.start(true);
		}
		ReturnType result;
		try {
			result = execute();
		} catch (const TimeoutException&) {
			result = ReturnType::TimeoutInfeasible;
		}
		if (isTimelimitEnabled()) {
			m_watch.stop();
		}
		return result;
	}

	/**
	 * @returns true, if the given \p GA and \p stretch are valid for a specific algorithm. If not,
	 * an error message is provided via \p error
	 */
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch, std::string& error) = 0;

	/**
	 * Sets the timelimit for the algorithm in milliseconds. Use -1 to disable the timelimit.
	 * On timelimit, the algorithm will be aborted. Note that 0 is a valid timelimit. Do not
	 * change this during the execution of the algorithm.
	 */
	void setTimelimit(int64_t milliseconds) { m_timelimit = milliseconds; }

	/**
	 * @returns the time in milliseconds needed for executing the algorithm (see the execute method).
	 */
	int64_t getTimeNeeded() { return m_watch.milliSeconds(); }

protected:
	const GraphAttributes* m_GA;
	double m_stretch;
	GraphCopySimple* m_spanner;
	EdgeArray<bool>* m_inSpanner;

	/**
	 * Initializes members and create an empty spanner.
	 */
	virtual void init(const GraphAttributes& GA, double stretch, GraphCopySimple& spanner,
			EdgeArray<bool>& inSpanner) {
		m_GA = &GA;
		m_stretch = stretch;
		m_spanner = &spanner;
		m_inSpanner = &inSpanner;

		const Graph& G = GA.constGraph();
		m_spanner->clear();
		m_spanner->createEmpty(G);
		for (node n : G.nodes) {
			m_spanner->newNode(n);
		}
		m_inSpanner->init(G, false);
	}

	/**
	 * Executes the core algorithm. Called after initialization. This method is used for the
	 * timelimit, so do not forget to call assertTimeLeft from time to time.
	 */
	virtual ReturnType execute() = 0;

	/**
	 * @returns true, if there should be a timelimit.
	 */
	bool isTimelimitEnabled() { return m_timelimit != -1; }

	/**
	 * @returns the amount of time in milliseconds left for execution.
	 */
	int64_t getTimeLeft() { return m_timelimit - getTimeNeeded(); }

	/**
	 * Assert, that time is left. If there is no time left (and the timelimit is active),
	 * an exception will be thwron to abort the execution.
	 */
	void assertTimeLeft() {
		if (isTimelimitEnabled() && getTimeLeft() <= 0) {
			throw TimeoutException();
		}
	}

private:
	int64_t m_timelimit = -1; //!< Timelimit in milliseconds. -1 disables the timeout.
	StopwatchCPU m_watch; //!< Used for keeping track of time.

protected:
	struct TimeoutException {
	}; //!< A simple exception used to exit from the execution, if the timelimit is reached.

public:
	/**
	 * Validates a spanner
	 *
	 * \returns true iff the \p spanner is a k-multiplicative spanner of \p GA
	 */
	static bool isMultiplicativeSpanner(const GraphAttributes& GA, const GraphCopySimple& spanner,
			double stretch);

	/**
	 * Calculates an all-pair shortest-path on \p spanner with the weights given by \p GA
	 */
	static void apspSpanner(const GraphAttributes& GA, const GraphCopySimple& spanner,
			NodeArray<NodeArray<TWeight>>& shortestPathMatrix);

protected:
	/**
	 * \returns the weight of an edge \p e from \p GA
	 */
	static TWeight getWeight(const GraphAttributes& GA, edge e);
};

template<typename TWeight>
void SpannerModule<TWeight>::apspSpanner(const GraphAttributes& GA, const GraphCopySimple& spanner,
		NodeArray<NodeArray<TWeight>>& shortestPathMatrix) {
	EdgeArray<TWeight> weights(spanner);
	for (edge e : spanner.edges) {
		weights[e] = getWeight(GA, spanner.original(e));
	}
	dijkstra_SPAP(spanner, shortestPathMatrix, weights);
}

template<typename TWeight>
bool SpannerModule<TWeight>::isMultiplicativeSpanner(const GraphAttributes& GA,
		const GraphCopySimple& spanner, double stretch) {
	const Graph& G = GA.constGraph();
	OGDF_ASSERT(&(spanner.original()) == &G);

	if (G.numberOfNodes() != spanner.numberOfNodes()) {
		return false;
	}

	NodeArray<NodeArray<TWeight>> originalDistances(G);
	EdgeArray<TWeight> weights(G);
	for (edge e : G.edges) {
		weights[e] = getWeight(GA, e);
	}
	dijkstra_SPAP(G, originalDistances, weights);

	NodeArray<NodeArray<TWeight>> newDistances(spanner);
	apspSpanner(GA, spanner, newDistances);

	EpsilonTest m_eps;
	for (edge e : G.edges) {
		node u = e->source();
		node v = e->target();
		TWeight originalDistance = originalDistances[u][v];
		TWeight newDistance = newDistances[spanner.copy(u)][spanner.copy(v)];
		if (m_eps.greater(static_cast<double>(newDistance), stretch * originalDistance)) {
			return false;
		}
	}

	return true;
}

template<>
inline int SpannerModule<int>::getWeight(const GraphAttributes& GA, edge e) {
	if (GA.has(GraphAttributes::edgeIntWeight)) {
		return GA.intWeight(e);
	} else {
		return 1;
	}
}

template<>
inline double SpannerModule<double>::getWeight(const GraphAttributes& GA, edge e) {
	if (GA.has(GraphAttributes::edgeDoubleWeight)) {
		return GA.doubleWeight(e);
	} else {
		return 1.0;
	}
}

}
