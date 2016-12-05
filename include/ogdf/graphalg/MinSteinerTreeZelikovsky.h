/** \file
 * \brief Implementation of Zelikovsky's 11/6-approximation algorithm
 * 	      for the minimum Steiner tree problem.
 *
 * \author Matthias Woste, Stephan Beyer
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

#include <ogdf/basic/List.h>
#include <ogdf/graphalg/Voronoi.h>
#include <ogdf/internal/steinertree/SaveStatic.h>
#include <ogdf/internal/steinertree/SaveEnum.h>
#include <ogdf/internal/steinertree/SaveDynamic.h>
#include <ogdf/internal/steinertree/Triple.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/module/MinSteinerTreeModule.h>
#include <ogdf/internal/steinertree/common_algorithms.h>

namespace ogdf {

/*!
 * \brief This class implements the 11/6-approximation algorithm by Zelikovsky
 * for the minimum Steiner tree problem along with variants and practical improvements.
 *
 * @ingroup ga-steiner
 *
 * This implementation is based on:
 *
 * (A. Zelikovsky, An 11/6-Approximation Algorithm for the Network Steiner Problem,
 * Algorithmica, volume 9, number 5, pages 463-470, Springer, 1993)
 *
 * (A. Zelikovsky, A faster approximation algorithm for the Steiner problem in graphs,
 * Information Processing Letters, volume 46, number 2, pages 79-83, 1993)
 *
 * (A. Zelikovsky, Better approximation bound for the network and euclidean Steiner
 * tree problems, Technical Report, 2006)
 */
template<typename T>
class MinSteinerTreeZelikovsky: public MinSteinerTreeModule<T> {
public:
	template<typename TYPE> using Save = steinertree::Save<TYPE>;
	template<typename TYPE> using Triple = steinertree::Triple<TYPE>;

	//! Choice of objective function
	enum WinCalculation {
		absolute, //!< win=gain-cost
		relative  //!< win=gain/cost
	};

	//! Choice of triple generation
	enum TripleGeneration {
		exhaustive, //!< try all possibilities
		voronoi, //!< use voronoi regions
		ondemand  //!< generate triples "on the fly", only usable with WinCalculation::absolute
	};

	//! Switches immediate triple dropping
	enum TripleReduction {
		on, //!< removes triples as soon as their gain is known to be non positive
		off //!< keeps triples all the time
	};

	//! Different methods for obtaining save edges
	enum SaveCalculation {
		//! Stores explicitly the save edge for every pair of terminals.
		//! Needs O(n^2) space but has fast query times
		staticEnum,
		//! Builds a "weight tree" (save edges are inner nodes, terminals are leaves
		//! and searches save edges via LCA calculation of two nodes
		staticLCATree,
		//! Same as staticLCATree but each time a triple has been contracted
		//! the "weight tree" is updated dynamically rather than completely
		//! new from scratch.  Has the fastest update time
		dynamicLCATree,
		//! Uses staticEnum for the triple generation phase (many queries)
		//! and dynamicLCATree during the contraction phase (few updates)
		hybrid
	};

	//! Enables a heuristic version (for TG exhaustive and voronoi only)
	enum Pass {
		//! heuristic: evaluate all triples, sort them descending by gain,
		//! traverse sorted triples once, contract when possible
		one,
		//! normal, greedy version
		multi
	};

	MinSteinerTreeZelikovsky(WinCalculation wc = absolute, TripleGeneration tg = voronoi, SaveCalculation sc = hybrid, TripleReduction tr = on, Pass pass = multi)
	 : m_winCalculation(wc)
	 , m_tripleGeneration(tg)
	 , m_saveCalculation(sc)
	 , m_tripleReduction(tr)
	 , m_pass(pass)
	 , m_ssspDistances(true)
	{
	}

	virtual ~MinSteinerTreeZelikovsky() { }

	virtual T call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree) override
	{
		m_triplesGenerated = 0;
		m_tripleLookUps = 0;
		m_triplesContracted = 0;
		return MinSteinerTreeModule<T>::call(G, terminals, isTerminal, finalSteinerTree);
	}

	/*!
	 * \brief For the 3-restricted case, it is sufficient to compute an SSSP from every terminal
	 *  instead of doing a full APSP. In case a full APSP is faster, use this method.
	 * @param force True to force APSP instead of SSSP.
	 */
	void forceAPSP(bool force = true)
	{
		m_ssspDistances = !force;
	}

	//! Sets type of gain calculation \see MinSteinerTreeZelikovsky::WinCalculation
	void winCalculation(WinCalculation wc)
	{
		m_winCalculation = wc;
	}

	//! Returns type of gain calculation currently in use \see MinSteinerTreeZelikovsky::WinCalculation
	WinCalculation winCalculation() const
	{
		return m_winCalculation;
	}

	//! Sets type of triple generation \see MinSteinerTreeZelikovsky::TripleGeneration
	void tripleGeneration(TripleGeneration tg)
	{
		m_tripleGeneration = tg;
	}

	//! Returns type of triple generation currently in use \see MinSteinerTreeZelikovsky::TripleGeneration
	TripleGeneration tripleGeneration() const
	{
		return m_tripleGeneration;
	}

	//! Sets type of triple reduction \see MinSteinerTreeZelikovsky::TripleReduction
	void tripleReduction(TripleReduction tr)
	{
		m_tripleReduction = tr;
	}

	//! Returns type of triple reduction currently in use \see MinSteinerTreeZelikovsky::TripleReduction
	TripleReduction tripleReduction() const
	{
		return m_tripleReduction;
	}

	//! Sets type of save calculation \see MinSteinerTreeZelikovsky::SaveCalculation
	void saveCalculation(SaveCalculation sv)
	{
		m_saveCalculation = sv;
	}

	//! Returns type of save calculation currently in use \see MinSteinerTreeZelikovsky::SaveCalculation
	SaveCalculation saveCalculation() const
	{
		return m_saveCalculation;
	}

	//! Sets type of pass \see MinSteinerTreeZelikovsky::Pass
	void pass(Pass p) {
		m_pass = p;
	}

	//! Returns type of pass currently in use \see MinSteinerTreeZelikovsky::Pass
	Pass pass() const
	{
		return m_pass;
	}

	//! Returns the number of generated triples
	long numberOfGeneratedTriples() const
	{
		return m_triplesGenerated;
	}

	//! Returns the number of contracted triples
	long numberOfContractedTriples() const
	{
		return m_triplesContracted;
	}

	//! Returns the number of triple lookups during execution time
	long numberOfTripleLookUps() const
	{
		return m_tripleLookUps;
	}

protected:
	/*!
	 * \brief Builds a minimum Steiner tree given a weighted graph and a list of terminals \see MinSteinerTreeModule::call
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param finalSteinerTree The final Steiner tree
	 * @return The objective value (sum of edge costs) of the final Steiner tree
	 */
	virtual T computeSteinerTree(
		const EdgeWeightedGraph<T> &G,
		const List<node> &terminals,
		const NodeArray<bool> &isTerminal,
		EdgeWeightedGraphCopy<T> *&finalSteinerTree) override;

	/*!
	 * \brief Computes the distance matrix for the original graph
	 */
	void computeDistanceMatrix();

	/*!
	 * \brief Mark unreachable nodes in distance matrix by setting their distance to infinity (numeric_limits<T>::max())
	 */
	inline void markUnreachablesInDistanceMatrix(NodeArray<T> &dist, const NodeArray<edge> &pred) {
		for (node u : m_originalGraph->nodes) {
			if (pred[u] == nullptr) {
				dist[u] = numeric_limits<T>::max();
			}
		}
	}

	/*!
	 * \brief Update center node if it is the best so far. (Just a helper to avoid code duplication.)
	 * @param x The node to test
	 * @param center The returned best center node
	 * @param minCost The returned cost of the component with that node
	 * @param dist1 SSSP distance vector of the first terminal
	 * @param dist2 SSSP distance vector of the second terminal
	 * @param dist3 SSSP distance vector of the third terminal
	 */
	inline void updateBestCenter(node x, node &center, T &minCost, const NodeArray<T> &dist1, const NodeArray<T> &dist2, const NodeArray<T> &dist3) const
	{
		const T tmp = dist1[x] + dist2[x] + dist3[x];
		if (minCost > tmp
		 && dist1[x] != numeric_limits<T>::max()
		 && dist2[x] != numeric_limits<T>::max()
		 && dist3[x] != numeric_limits<T>::max()) {
			center = x;
			minCost = tmp;
		}
	}

	/*!
	 * \brief Add a found triple to the triples list. (Just a helper to avoid code duplication.)
	 */
	inline void generateTriple(node u, node v, node w, node center, const T &minCost, const Save<T> &save)
	{
		if (center // center found
		 && !(*m_isTerminal)[center]) { // cheapest center is not a terminal
			const double gain = save.gain(u, v, w);
			const double win = calcWin(gain, minCost);
			if (tripleReduction() == off
			 || win > 0) {
				++m_triplesGenerated;
				OGDF_ASSERT(center);
				Triple<T> triple(u, v, w, center, minCost, win);
				m_triples.pushBack(triple);
			}
		}
	}

	/*!
	 * \brief Generates triples according to voronoi regions
	 * @param save data structure for calculation save edges
	 */
	void generateVoronoiTriples(const Save<T> &save);

	/*!
	 * \brief Generates all possible triples
	 * @param save data structure for determining save edges
	 */
	void generateExhaustiveTriples(const Save<T> &save);

	/*!
	 * \brief Generates triples according to the chosen option \see TripleGeneration
	 * @param save data structure for calculation save edges
	 */
	void generateTriples(const Save<T> &save)
	{
		if (tripleGeneration() == voronoi) {
			generateVoronoiTriples(save);
		} else
		if (tripleGeneration() == exhaustive) {
			generateExhaustiveTriples(save);
		}
	}

	/*!
	 * \brief Contracts a triple and updates the save data structure
	 * @param triple triple to be contracted
	 * @param save save data structure
	 * @param isNewTerminal true for nodes to be interpreted as terminals
	 */
	void contractTriple(const Triple<T> &triple, Save<T> &save, NodeArray<bool> &isNewTerminal)
	{
		++m_triplesContracted;
		save.update(triple);
		isNewTerminal[triple.z()] = true;
	}

	/*!
	 * \brief Contraction phase for algorithm generating triples on demand \see MinSteinerTreeZelikovsky::ondemand
	 * @param save save data structure
	 * @param isNewTerminal true for nodes to be interpreted as terminals
	 */
	void tripleOnDemand(Save<T> &save, NodeArray<bool> &isNewTerminal);

	/*!
	 * \brief Contraction phase for the original version of the algorithm \see MinSteinerTreeZelikovsky::multi
	 * @param save save data structure
	 * @param isNewTerminal true for nodes to be interpreted as terminals
	 */
	void multiPass(Save<T> &save, NodeArray<bool> &isNewTerminal);

	/*!
	 * \brief Contraction phase for the one pass heuristic \see MinSteinerTreeZelikovsky::one
	 * @param save save data structure
	 * @param isNewTerminal true for nodes to be interpreted as terminals
	 */
	void onePass(Save<T> &save, NodeArray<bool> &isNewTerminal);

	/*!
	 * \brief Calculate the win
	 */
	double calcWin(double gain, T cost) const
	{
		switch (winCalculation()) {
		case relative:
			return gain / cost - 1.0;
		case absolute:
		default:
			return gain - cost;
		}
	}

	void generateInitialTerminalSpanningTree(EdgeWeightedGraphCopy<T> &steinerTree)
	{
		// generate complete graph
		for (node v : *m_terminals) {
			steinerTree.newNode(v);
		}
		for (node u : steinerTree.nodes) {
			const NodeArray<T> &dist = m_distance[steinerTree.original(u)];
			for (node v = u->succ(); v; v = v->succ()) {
				steinerTree.newEdge(u, v, dist[steinerTree.original(v)]);
			}
		}
		// compute MST
		makeMinimumSpanningTree(steinerTree, steinerTree.edgeWeights());
	}

private:
	WinCalculation m_winCalculation; //!< Chosen option for gain calculation \see WinCalculation
	TripleGeneration m_tripleGeneration; //!< Chosen option for triple generation \see TripleGeneration
	SaveCalculation m_saveCalculation; //!< Chosen option for save calculation \see SaveCalculation
	TripleReduction m_tripleReduction; //!< Chosen option for triple reduction \see TripleReduction
	Pass m_pass; //!< Chosen option for pass \see Pass
	bool m_ssspDistances; //!< True iff we only compute SSSP from terminals instead of APSP for full component construction

	const EdgeWeightedGraph<T> *m_originalGraph; //!< The original edge-weighted graph
	const NodeArray<bool> *m_isTerminal; //!< Incidence vector for terminal nodes
	const List<node> *m_terminals; //!< List of terminal nodes
	NodeArray<NodeArray<T>> m_distance; //!< The distance matrix
	List<Triple<T>> m_triples; //!< The list of triples during the algorithm

	long m_triplesGenerated; //!< Number of generated triples
	long m_triplesContracted; //!< Number of contracted triples
	long m_tripleLookUps; //!< Number of triple lookups

};

} // end namespace ogdf

// ============= Implementation =================

namespace ogdf {

template<typename T>
T MinSteinerTreeZelikovsky<T>::computeSteinerTree(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree)
{
	OGDF_ASSERT(tripleGeneration() != ondemand // combinations that only work with ondemand:
	 || (winCalculation() == absolute
	  && saveCalculation() != hybrid
	  && pass() != one));

	m_originalGraph = &G;
	m_terminals = &terminals;
	m_isTerminal = &isTerminal;

	// build the complete terminal graph and a distance matrix
	NodeArray<bool> isNewTerminal(G, false);
	for (node v : *m_terminals) {
		isNewTerminal[v] = true;
	}

	if (terminals.size() >= 3) {
		m_distance.init(G); // the distance matrix
		computeDistanceMatrix();

		// init terminal-spanning tree and its save-edge data structure
		EdgeWeightedGraphCopy<T> steinerTree; // the terminal-spanning tree to be modified
		steinerTree.createEmpty(G);
		generateInitialTerminalSpanningTree(steinerTree);

		Save<T> *save = nullptr;
		switch (saveCalculation()) {
		case staticEnum:
			save = new steinertree::SaveEnum<T>(steinerTree);
			break;
		case staticLCATree:
			save = new steinertree::SaveStatic<T>(steinerTree);
			break;
		case dynamicLCATree:
		case hybrid:
			save = new steinertree::SaveDynamic<T>(steinerTree);
			break;
		}
		OGDF_ASSERT(save);

		if (tripleGeneration() == ondemand) { // ondemand triple generation
			tripleOnDemand(*save, isNewTerminal);
		} else { // exhaustive or voronoi
			// triple generation phase
			if (saveCalculation() == hybrid) {
				steinertree::SaveEnum<T> saveTriple(steinerTree);
				generateTriples(saveTriple);
			} else {
				generateTriples(*save);
			}
			// contraction phase
			if (pass() == multi) {
				multiPass(*save, isNewTerminal);
			} else { // pass() == one
				onePass(*save, isNewTerminal);
			}
		}
		delete save;

		// cleanup
		m_triples.clear();
	}

	// obtain final Steiner Tree using (MST-based) Steiner tree approximation algorithm
	return steinertree::obtainFinalSteinerTree(G, isNewTerminal, isTerminal, finalSteinerTree);
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::computeDistanceMatrix()
{
	if (m_ssspDistances) {
		NodeArray<edge> pred;
		for (node u : *m_terminals) {
			MinSteinerTreeModule<T>::singleSourceShortestPathsStrict(*m_originalGraph, u, *m_isTerminal, m_distance[u], pred);
			markUnreachablesInDistanceMatrix(m_distance[u], pred);
		}
	} else {
		NodeArray<NodeArray<edge>> pred(*m_originalGraph);
		MinSteinerTreeModule<T>::allPairShortestPathsStrict(*m_originalGraph, *m_isTerminal, m_distance, pred);
		for (node u : *m_terminals) {
			markUnreachablesInDistanceMatrix(m_distance[u], pred[u]);
		}
	}
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::generateVoronoiTriples(const Save<T> &save)
{
	Voronoi<T> voronoi(*m_originalGraph, m_originalGraph->edgeWeights(), *m_terminals);
	for (ListConstIterator<node> it_u = m_terminals->begin(); it_u.valid(); ++it_u) {
		const node u = *it_u;
		const NodeArray<T> &uDistance = m_distance[u];
		for (ListConstIterator<node> it_v = it_u.succ(); it_v.valid(); ++it_v) {
			const node v = *it_v;
			const NodeArray<T> &vDistance = m_distance[v];
			for (ListConstIterator<node> it_w = it_v.succ(); it_w.valid(); ++it_w) {
				const node w = *it_w;
				const NodeArray<T> &wDistance = m_distance[w];

				node center = nullptr;
				T minCost = numeric_limits<T>::max();
				// look in all Voronoi regions for the best center node
				for (node x : voronoi.nodesInRegion(u)) {
					updateBestCenter(x, center, minCost, uDistance, vDistance, wDistance);
				}
				for (node x : voronoi.nodesInRegion(v)) {
					updateBestCenter(x, center, minCost, uDistance, vDistance, wDistance);
				}
				for (node x : voronoi.nodesInRegion(w)) {
					updateBestCenter(x, center, minCost, uDistance, vDistance, wDistance);
				}
				generateTriple(u, v, w, center, minCost, save);
			}
		}
	}
}


template<typename T>
void MinSteinerTreeZelikovsky<T>::generateExhaustiveTriples(const Save<T> &save)
{
	for (ListConstIterator<node> it_u = m_terminals->begin(); it_u.valid(); ++it_u) {
		const node u = *it_u;
		const NodeArray<T> &uDistance = m_distance[u];
		for (ListConstIterator<node> it_v = it_u.succ(); it_v.valid(); ++it_v) {
			const node v = *it_v;
			const NodeArray<T> &vDistance = m_distance[v];
			for (ListConstIterator<node> it_w = it_v.succ(); it_w.valid(); ++it_w) {
				const node w = *it_w;
				const NodeArray<T> &wDistance = m_distance[w];

				node center = nullptr;
				T minCost = numeric_limits<T>::max();
				for (node x : m_originalGraph->nodes) {
					updateBestCenter(x, center, minCost, uDistance, vDistance, wDistance);
				}
				generateTriple(u, v, w, center, minCost, save);
			}
		}
	}
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::tripleOnDemand(Save<T> &save, NodeArray<bool> &isNewTerminal)
{
	Triple<T> maxTriple;
	do {
		maxTriple.win(0);
		for (node v : m_originalGraph->nodes) {
			if (!(*m_isTerminal)[v]) { // for each non-terminal v
				// find s0, nearest terminal to v
				T best = numeric_limits<T>::max();
				node s0 = nullptr;
				for (node s : *m_terminals) {
					T tmp = m_distance[s][v];
					if (best > tmp) {
						best = tmp;
						s0 = s;
					}
				}
				OGDF_ASSERT(s0);

				// find s1 maximizing save(s0, s1) - d(v, s1)
				node s1 = nullptr;
				T save1Val(0);
				for (node s : *m_terminals) {
					if (s == s0
					 || m_distance[s][v] == numeric_limits<T>::max()) {
						continue;
					}
					T tmpVal = save.saveWeight(s, s0);
					T tmp = tmpVal - m_distance[s][v];
					if (!s1 || best < tmp) {
						best = tmp;
						s1 = s;
						save1Val = tmpVal;
					}
				}
				if (!s1) { // it may happen that s1 does not exist
					continue;
				}

				node s2 = nullptr;
				T save2Val(0);
				const edge save1 = save.saveEdge(s0, s1);
				for (node s : *m_terminals) {
					if (s == s0
					 || s == s1
					 || m_distance[s][v] == numeric_limits<T>::max()) {
						continue;
					}
					const edge tmp = save.saveEdge(s0, s);
					if (tmp == save1) {
						save2Val = save.saveWeight(s1, s);
					} else {
						save2Val = save.saveWeight(s0, s);
					}
					T tmpWin = save1Val + save2Val - m_distance[s0][v] - m_distance[s1][v] - m_distance[s][v];
					if (!s2 || best < tmpWin) {
						best = tmpWin;
						s2 = s;
					}
				}
				if (!s2) { // it may happen that s2 does not exist
					continue;
				}

				if (s2 // it may happen that s2 does not exist
				 && best > maxTriple.win()) { // best win is better than previous best; also positive
					++m_triplesGenerated;
					maxTriple.s0(s0);
					maxTriple.s1(s1);
					maxTriple.s2(s2);
					maxTriple.z(v);
					maxTriple.win(best);
					//maxTriple.cost(save1Val + save2Val - win); not needed
				}
			}
		}

		if (maxTriple.win() > 0) {
			contractTriple(maxTriple, save, isNewTerminal);
		}
	} while (maxTriple.win() > 0);
}


template<typename T>
void MinSteinerTreeZelikovsky<T>::onePass(Save<T> &save, NodeArray<bool> &isNewTerminal)
{
	class TripleComparer {
	public:
		static double compare(const Triple<T> &x1, const Triple<T> &x2) {
			return x2.win() - x1.win();
		}
		OGDF_AUGMENT_STATICCOMPARER(Triple<T>)
	} tc;
	m_triples.quicksort(tc);

	for (const Triple<T> &t : m_triples) {
		++m_tripleLookUps;
		if (calcWin(double(save.gain(t.s0(), t.s1(), t.s2())), t.cost()) > 0) {
			contractTriple(t, save, isNewTerminal);
		}
	}
}


template<typename T>
void MinSteinerTreeZelikovsky<T>::multiPass(Save<T> &save, NodeArray<bool> &isNewTerminal)
{
	double win = 0;
	ListIterator<Triple<T>> maxTriple;

	do {
		win = 0;
		ListIterator<Triple<T>> nextIt;
		for (ListIterator<Triple<T>> it = m_triples.begin(); it.valid(); it = nextIt) {
			nextIt = it.succ();
			++m_tripleLookUps;
			Triple<T> &t = *it;
			t.win(calcWin(double(save.gain(t.s0(), t.s1(), t.s2())), t.cost()));
			if (t.win() > win) {
				win = t.win();
				maxTriple = it;
			} else {
				if (tripleReduction() == on
				 && t.win() <= 0) {
					m_triples.del(it);
				}
			}
		}

		if (win > 0) {
			contractTriple(*maxTriple, save, isNewTerminal);
			m_triples.del(maxTriple);
		}
	} while (win > 0);
}


}
// end namespace ogdf
