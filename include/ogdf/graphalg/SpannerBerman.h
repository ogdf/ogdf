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

#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/external/coin.h>
#include <ogdf/graphalg/SpannerModule.h>

#include <iomanip>

namespace ogdf {

/**
 * Approximation algorithm for calculating spanners.
 *
 * P. Berman, A. Bhattacharyya, K. Makarychev, S. Raskhodnikova und G. Yaroslavtsev.
 * Approximation algorithms for spanner problems and Directed Steiner Forest. Information
 * and Computation 222 (2013). 38th International Colloquium on Automata, Languages
 * and Programming (ICALP 2011), S. 93–107. doi: https://doi.org/10.1016/j.ic.
 * 2012.10.007.
 *
 * Conditions for the graph:
 * - simple
 * - weakly connected: The graph must be connected, if acrs are seen as edges.
 *   To calculate a spanner on disconnected graphs use SpannerBermanDisconnected.
 *
 * The stretch \f$s\f$ must satisfy: \f$s\geq1\in\mathbb{R}\f$.
 *
 * The preconditions can be checked with SpannerBerman::preconditionsOk.
 *
 * Calculates a k-spanner with an approximation ratio of \f$\mathcal{O}(n^{1/2}\log n)\f$ for
 * the size. The graph can be directed and weighted. Undirected and/or unweighted graphs
 * still preserve the approximation ratio.
 *
 * To enable logging you have to set the log level of the algorithm's logger to
 * Logger::Level::Default or below:
 * \code
 * SpannerBerman<int>::logger.localLogLevel(Logger::Level::Default);
 * // or
 * SpannerBerman<double>::logger.localLogLevel(Logger::Level::Default);
 * \endcode
 *
 * @ingroup ga-spanner
 */
template<typename TWeight>
class SpannerBerman : public SpannerModule<TWeight> {
public:
	static Logger logger;

	SpannerBerman() {
		m_OPT = 0;
		m_osi = nullptr;
	}

	virtual ~SpannerBerman() { resetLP(); }

	//! @copydoc ogdf::SpannerModule::preconditionsOk
	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		if (!isSimple(GA.constGraph())) {
			error = "The graph is not simple";
			return false;
		}
		if (!isConnected(GA.constGraph())) {
			error = "The graph is not connected";
			return false;
		}
		if (stretch < 1.0) {
			error = "The stretch must be >= 1.0";
			return false;
		}
		return true;
	}

private:
	/**
	 * Used to indicate the result of the separation method
	 */
	enum class SeparationResult { NewConstraint, Fail, Solution };

	const Graph* m_G; //!< const reference to the original graph
	EdgeArray<TWeight> m_weight; //!< weights of each edge from m_G
	EdgeArray<TWeight> m_spannerWeight; //!< weights for m_spanner
	//! Caches, if an edge is thick (or thin). Use isThinEdge and
	//! isThickEdge for access. It is fully cached after calling calculateThickEdges
	EdgeArray<bool> m_isThickEdge;

	// Some precalculated values that are needed all over the place.
	int m_nSquared; //!< n^2
	double m_beta; //!< sqrt(n)
	double m_sqrtlog; //!< sqrt(n) * ln(n)
	int m_thickEdgeNodeAmountLimit; //!< n/beta

	EdgeArray<bool> m_E1; //!< Holds the set E1 from the first part of the algorithm
	EdgeArray<bool> m_E2; //!< Holds the set E2 from the second part of the algorithm

	NodeArray<NodeArray<TWeight>> m_inDistance; //!< m_inDistance[v][w]: distance from v in a v-rooted in-arborescense to w
	NodeArray<NodeArray<TWeight>> m_outDistance; //!< m_outDistance[v][w]: distance from v in a v-rooted out-arborescense to w

	// Solving the LP
	OsiSolverInterface* m_osi; //!< the solver. Initial nullptr since it is initialized in the second phase.
	std::vector<CoinPackedVector*> m_constraints; //!< Holds all constraints so they can be freed at destruction.
	int m_OPT; //!< The current guess for our optimal LP value.

	EpsilonTest m_eps;

	inline bool isThinEdge(edge e) { return !m_isThickEdge(e); }

	inline bool isThickEdge(edge e) { return m_isThickEdge(e); }

	/**
	 * Resets the LP defining variables.
	 *
	 * Deletes m_osi as well as every entry in m_constraints. Sets both variables to nullptr/empty list.
	 */
	void resetLP() {
		for (auto c : m_constraints) {
			delete c;
		}
		m_constraints.clear();

		delete m_osi;
		m_osi = nullptr;
	}

	//! @copydoc ogdf::SpannerModule::init
	virtual void init(const GraphAttributes& GA, double stretch, GraphCopySimple& spanner,
			EdgeArray<bool>& inSpanner) override {
		resetLP();
		SpannerModule<TWeight>::init(GA, stretch, spanner, inSpanner);

		m_G = &GA.constGraph();
		m_weight.init(*m_G);
		for (edge e : m_G->edges) {
			m_weight[e] = getWeight(*m_GA, e);
		}
		m_spannerWeight.init(*m_spanner);
		m_isThickEdge.init(*m_G, false);

		m_nSquared = m_G->numberOfNodes() * m_G->numberOfNodes();
		m_beta = sqrt(m_G->numberOfNodes()),
		m_thickEdgeNodeAmountLimit = ceil(m_G->numberOfNodes() / m_beta);
		m_sqrtlog = sqrt(m_G->numberOfNodes()) * log(m_G->numberOfNodes());

		m_E1.init(*m_G, false);
		m_E2.init(*m_G, false);

		m_inDistance.init(*m_G);
		m_outDistance.init(*m_G);
	}

	//! @copydoc ogdf::SpannerModule::execute
	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		if (m_G->numberOfNodes() == 0 || m_G->numberOfEdges() == 0) {
			return SpannerModule<TWeight>::ReturnType::Feasible;
		}
		firstPart();
		printStats();

		logger.lout() << "\n2. Part\n" << std::endl;
		bool ok = secondPart();

		if (!ok) {
			return SpannerModule<TWeight>::ReturnType::Error;
		}

		printStats(true);

		// Fill m_inSpanner since only E1, E2 and m_spanner are correctly filled.
		for (edge e : m_G->edges) {
			(*m_inSpanner)[e] = m_E1[e] || m_E2[e];
		}
		return SpannerModule<TWeight>::ReturnType::Feasible;
	}

	/**
	 * First part of the algorithm: Settle all thick edges.
	 */
	void firstPart() {
		NodeArray<NodeArray<edge>> inPredecessor(*m_G);
		NodeArray<NodeArray<edge>> outPredecessor(*m_G);
		for (node n : m_G->nodes) {
			inArborescence(*m_GA, n, inPredecessor[n], m_inDistance[n]);
			outArborescence(*m_GA, n, outPredecessor[n], m_outDistance[n]);
			assertTimeLeft();
		}

		int max = ceil(m_beta * log(m_G->numberOfNodes())); // log is ln
		logger.lout() << "max=" << max << std::endl;
		for (int i = 0; i < max; i++) {
			node v = m_G->chooseNode();
			for (node n : m_G->nodes) {
				edge e1 = inPredecessor[v][n];
				if (e1) {
					addEdgeToSpanner(e1);
				}

				edge e2 = outPredecessor[v][n];
				if (e2) {
					addEdgeToSpanner(e2);
				}
			}
		}
		assertTimeLeft();
		logger.lout() << "thickEdgeNodeAmountLimit=" << m_thickEdgeNodeAmountLimit << std::endl;
		printStats();

		logger.lout() << "add unsettled thick edges" << std::endl;
		calculateThickEdges();
		addUnsettledThickEdges();
	}

	/**
	 * Calculates an in-arborescense rooted at \p root
	 */
	void inArborescence(const GraphAttributes& GA, node root, NodeArray<edge>& predecessor,
			NodeArray<TWeight>& distance) {
		Dijkstra<TWeight> dijkstra;
		dijkstra.callUnbound(*m_G, m_weight, root, predecessor, distance, true, true);
	}

	/**
	 * Calculates an out-arborescense rooted at \p root
	 */
	void outArborescence(const GraphAttributes& GA, node root, NodeArray<edge>& predecessor,
			NodeArray<TWeight>& distance) {
		Dijkstra<TWeight> dijkstra;
		dijkstra.callUnbound(*m_G, m_weight, root, predecessor, distance, true, false);
	}

	/**
	 * Add an edge to the spanner. Only used in the first phase!
	 */
	void addEdgeToSpanner(edge e) {
		if (!m_E1[e]) {
			m_E1[e] = true;
			m_spanner->newEdge(e);
			m_spannerWeight[m_spanner->copy(e)] = m_weight[e];
		}
	}

	void printStats(bool assert = false) {
		logger.lout() << "\nStats:" << std::endl;
		logger.lout() << m_spanner->numberOfEdges() << " covered of " << m_G->numberOfEdges()
					  << std::endl;

		int E1amount = 0;
		int E1amountThick = 0;
		int E1amountThinNotInE2 = 0;
		int E2amount = 0;
		int E2amountThin = 0;
		int E2amountThickNotInE1 = 0;
		int edgesCovered = 0;
		int duplicateCovered = 0;
		for (edge e : m_G->edges) {
			if (m_E1[e]) {
				E1amount++;
				if (isThickEdge(e)) {
					E1amountThick++;
				} else if (!m_E2[e]) {
					E1amountThinNotInE2++;
				}
			}
			if (m_E2[e]) {
				E2amount++;
				if (isThinEdge(e)) {
					E2amountThin++;
				} else if (!m_E1[e]) {
					E2amountThickNotInE1++;
				}
			}
			if (m_E1[e] || m_E2[e]) {
				edgesCovered++;
			}
			if (m_E1[e] && m_E2[e]) {
				duplicateCovered++;
			}
		}

		logger.lout() << "covered: " << edgesCovered << ", duplicate edges: " << duplicateCovered
					  << std::endl;
		logger.lout() << "E1: " << E1amount << " edges, " << E1amountThick << " thick, "
					  << (E1amount - E1amountThick) << " thin, " << E1amountThinNotInE2
					  << " thin edges not in E2" << std::endl;
		logger.lout() << "E2: " << E2amount << " edges, " << E2amountThin << " thin, "
					  << (E2amount - E2amountThin) << " thick, " << E2amountThickNotInE1
					  << " thick edges not in E1" << std::endl;

#if defined(OGDF_DEBUG)
		if (assert) {
			OGDF_ASSERT((E1amountThick + E2amountThin + E1amountThinNotInE2 + E2amountThickNotInE1)
					== m_spanner->numberOfEdges());
		}
#endif

		OGDF_ASSERT(edgesCovered == m_spanner->numberOfEdges());
	}

	void calculateThickEdges() {
		for (edge e : m_G->edges) {
			m_isThickEdge[e] = _isThickEdge(e);
		}
	}

	/**
	 * Actually calculates whether \p e is thick or not.
	 */
	bool _isThickEdge(edge e) {
		const TWeight maxDistance = m_stretch * m_weight[e];

		int amountNodesInShortestPaths = 0;
		for (node n : m_G->nodes) {
			TWeight sd = m_outDistance[e->source()][n];
			TWeight td = m_inDistance[e->target()][n];
			if (sd == std::numeric_limits<TWeight>::max()
					|| td == std::numeric_limits<TWeight>::max()) {
				continue;
			}

			OGDF_ASSERT(m_eps.geq(sd + td, static_cast<TWeight>(0)));
			OGDF_ASSERT(m_eps.less(sd + td, std::numeric_limits<TWeight>::max()));

			if (m_eps.leq((sd + td), maxDistance)) {
				amountNodesInShortestPaths++;
			}

			if (amountNodesInShortestPaths >= m_thickEdgeNodeAmountLimit) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Shortcut for calling isSettledEdge with m_spanner and the current spanner weights.
	 */
	bool isSettledEdge(const edge e) { return isSettledEdge(e, *m_spanner, m_spannerWeight); }

	/**
	 * @returns true, if the edge is settles in the given spanner.
	 */
	bool isSettledEdge(const edge e, const GraphCopySimple& _spanner,
			const EdgeArray<TWeight>& _spannerWeight) {
		const TWeight maxDistance = m_stretch * m_weight[e];
		node u = _spanner.copy(e->source());
		node v = _spanner.copy(e->target());
		TWeight currentSpannerDistance = distance(_spanner, _spannerWeight, u, v, ceil(maxDistance));
		return m_eps.leq(currentSpannerDistance, maxDistance);
	}

	TWeight distance(const GraphCopySimple& G, const EdgeArray<TWeight>& weights, const node s,
			const node t, int maxLookupDist) {
		NodeArray<TWeight> distances;
		NodeArray<edge> predecessor;
		Dijkstra<TWeight> dijkstra;
		dijkstra.callBound(G, weights, s, predecessor, distances, m_GA->directed(),
				false, // arcs are not reversed
				t, maxLookupDist);
		return distances[t];
	}

	void addUnsettledThickEdges() {
		for (edge e : m_G->edges) {
			if (m_E1[e]) {
				continue;
			}

			if (isThickEdge(e)) {
				if (!isSettledEdge(e)) {
					addEdgeToSpanner(e);
				}
				assertTimeLeft();
			}
		}
	}

	/**
	 * The second part: Settling all thin edges.
	 */
	bool secondPart() {
		logger.lout() << "Using LP-Solver: " << Configuration::whichLPSolver() << std::endl;
		m_osi = CoinManager::createCorrectOsiSolverInterface();
		m_osi->messageHandler()->setLogLevel(0); // 0=nothing .. 4=verbose
		m_osi->setObjSense(1); // minimize

		// One column per edge (x_e)
		CoinPackedVector zero;
		EdgeArray<int> indices(*m_G);
		int i = 0;
		for (edge e : m_G->edges) {
			m_osi->addCol(zero, 0, 1, 1); // vec, lower, upper, objective
			indices[e] = i++;
		}


		CoinPackedVector* optConstraint = new CoinPackedVector;
		for (edge e : m_G->edges) {
			optConstraint->insert(indices[e], 1);
		}
		// sense: 'E' ==   'G' >=   'L' <=
		m_osi->addRow(*optConstraint, 'L', 0, 0); // constraint, sense, rhs (will be set below), range
		m_constraints.push_back(optConstraint);

		// Set the initial value of the OPT constraint rhs to n-1
		setOpt(m_G->numberOfNodes() - 1);

		int amountSolverCalls = 0;
		while (true) {
			if (amountSolverCalls == 0) {
				m_osi->initialSolve();
			} else {
				m_osi->resolve();
			}
			amountSolverCalls++;

			assertTimeLeft();

			logger.lout() << amountSolverCalls << ". solve... ";
			if (m_osi->isProvenOptimal()) {
				logger.lout() << "-> optimal (" << m_osi->getObjValue() << ")" << std::endl;
				logger.lout() << "  separation... ";
				SeparationResult result = separation(m_osi->getColSolution(), indices);

				if (result == SeparationResult::Solution) {
					for (edge e : m_G->edges) {
						if (m_E2[e] && !m_E1[e]) {
							m_spanner->newEdge(e);
							m_spannerWeight[m_spanner->copy(e)] = m_weight[e];
						}
					}
					logger.lout() << "-> Found solution\n" << std::endl;
					logger.lout() << "solver calls: " << amountSolverCalls << std::endl;
					logger.lout() << "cols: " << m_osi->getNumCols() << std::endl;
					logger.lout() << "rows: " << m_osi->getNumRows() << std::endl;
					return true;
				} else if (result == SeparationResult::NewConstraint) {
					logger.lout() << "-> New constraint" << std::endl;
				} else { // Fail
					logger.lout() << "-> Failed" << std::endl;
					if (!setOpt(m_OPT + 1)) {
						return false;
					}
				}
			} else if (m_osi->isProvenPrimalInfeasible()) {
				logger.lout() << "-> Infeasible" << std::endl;
				if (!setOpt(m_OPT + 1)) {
					return false;
				}
			} else if (m_osi->isProvenDualInfeasible()) {
				logger.lout() << "-> Unbounded" << std::endl;
				return false;
			} else {
				logger.lout() << "-> No solution found" << std::endl;
				return false;
			}
		};
	}

	/**
	 * Set a new value for m_OPT. Logs about the new value and modifies the LP.
	 *
	 * @returns false, if the new value is too large. If so, the calculation can be aborted.
	 */
	bool setOpt(int opt) {
		m_OPT = opt;
		if (m_OPT > m_nSquared) {
			logger.lout() << "  OPT is too large. Abort." << std::endl;
			return false;
		}
		float percentage = static_cast<float>(m_OPT) / m_nSquared * 100.0;
		logger.lout() << "  Set OPT to " << m_OPT << " (" << std::fixed << std::setprecision(2)
					  << percentage << "% of " << m_nSquared << ")" << std::endl;

		m_osi->setRowBounds(0, 0.0,
				m_OPT); // this is the opt constraint; it is 0 <= sum(x_e) <= m_OPT
		return true;
	}

	SeparationResult separation(const double* solution, const EdgeArray<int>& indices) {
		EdgeArray<bool> out(*m_G, false);
		randomizedSelection(solution, out);

		GraphCopySimple copy;
		copy.createEmpty(*m_G);
		for (node n : m_G->nodes) {
			copy.newNode(n);
		}
		EdgeArray<TWeight> copyWeight(copy);
		int amountCoveredEdges = 0;
		for (edge e : m_G->edges) {
			if (out[e]) {
				amountCoveredEdges++;
				copy.newEdge(e);
				copyWeight[copy.copy(e)] = m_weight[e];
			}
		}

		bool settlesAllThinEdges = true;
		edge unsettledThinEdge = nullptr;
		for (edge e : m_G->edges) {
			if (isThinEdge(e) && !isSettledEdge(e, copy, copyWeight)) {
				settlesAllThinEdges = false;
				unsettledThinEdge = e;
				break;
			}
		}

		if (settlesAllThinEdges) {
			if (amountCoveredEdges <= ceil(2.0 * m_OPT * m_sqrtlog)) {
				m_E2 = out;
				return SeparationResult::Solution;
			} else {
				return SeparationResult::Fail;
			}
		} else {
			EdgeArray<bool> antispanner(*m_G, false);
			createAntispanner(unsettledThinEdge, out, antispanner);

			double rowValue = 0.0;
			int i = 0;
			for (edge e : m_G->edges) {
				if (antispanner[e]) {
					rowValue += solution[i];
				}
				i++;
			}
			if (m_eps.less(rowValue, 1.0)) {
				// create new constraint
				CoinPackedVector* c = new CoinPackedVector;
				for (edge e : m_G->edges) {
					if (antispanner[e]) {
						c->insert(indices[e], 1);
					}
				}
				m_osi->addRow(*c, 'G', 1, 0);
				m_constraints.push_back(c);
				return SeparationResult::NewConstraint;
			} else {
				return SeparationResult::Fail;
			}
		}
	}

	void randomizedSelection(const double* fractions, EdgeArray<bool>& out) {
		int i = 0;
		for (edge e : m_G->edges) {
			double p = m_sqrtlog * fractions[i++]; // P may not be limited to 1: if it is >1 the
					// result will always be true, which is ok.
			out[e] = randomDouble(0.0, 1.0) <= p;
		}
	}

	void createAntispanner(const edge unsettledThinEdge, const EdgeArray<bool>& out,
			EdgeArray<bool>& antispanner) {
		GraphCopySimple copy;
		copy.createEmpty(*m_G);
		for (node n : m_G->nodes) {
			copy.newNode(n);
		}
		EdgeArray<TWeight> copyWeight(copy);

		const TWeight maxDistance = m_stretch * m_weight[unsettledThinEdge];
		for (edge e : m_G->edges) {
			// s->e->t
			TWeight s_e = m_outDistance[unsettledThinEdge->source()][e->source()];
			TWeight e_t = m_inDistance[unsettledThinEdge->target()][e->target()];
			bool isInLocalSubgraph = (s_e != std::numeric_limits<TWeight>::max()
					&& e_t != std::numeric_limits<TWeight>::max()
					&& m_eps.leq(s_e + m_weight[e] + e_t, maxDistance));

			// Graph to maximize is E_st\(E_st\E2) = E_st intersection E2
			// intersection[e] = isInLocalSubgraph && out[e];
			if (isInLocalSubgraph && out[e]) {
				copy.newEdge(e);
				copyWeight[copy.copy(e)] = m_weight[e];
			}

			// e in antispanner <=> e in E_st and not in E2 (=out)
			antispanner[e] = isInLocalSubgraph && !out[e];
		}
		assertTimeLeft();

		// minimize antispanner
		bool changed;
		do {
			changed = false;
			for (edge e : m_G->edges) {
				if (!antispanner[e]) {
					continue;
				}

				// try to remove e and check, if still no path <=stretch*max exists
				// -> removing e from antispanner is equivalent to adding e to the graph
				copy.newEdge(e);
				copyWeight[copy.copy(e)] = m_weight[e];
				antispanner[e] = false;

				if (!isSettledEdge(unsettledThinEdge, copy, copyWeight)) {
					// we've found an unnecessary edge.
					changed = true;
				} else {
					antispanner[e] = true;
					copy.delEdge(copy.copy(e));
				}
			}

			assertTimeLeft();
		} while (changed);
	}

	using SpannerModule<TWeight>::getWeight;
	using SpannerModule<TWeight>::assertTimeLeft;
	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_stretch;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

template<typename TWeight>
Logger SpannerBerman<TWeight>::logger(Logger::LogMode::Log,
		Logger::Level::High); // do not log by default

}
