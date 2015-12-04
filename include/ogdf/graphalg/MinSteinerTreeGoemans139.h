/** \file
 * \brief Implementation of an LP-based 1.39+epsilon Steiner tree
 * approximation algorithm by Goemans et al.
 *
 * \author Stephan Beyer
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#pragma once

#include <coin/CoinPackedVector.hpp>
#include <coin/OsiClpSolverInterface.hpp>
#include <ogdf/external/coin.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/Stopwatch.h>
#include <ogdf/internal/steinertree/FullComponentStore.h>
#include <ogdf/internal/steinertree/DreyfusWagnerFullComponentGenerator.h>
#include <ogdf/internal/steinertree/common_algorithms.h>
#include <ogdf/graphalg/MinCostFlowReinelt.h>
#include <ogdf/graphalg/MaxFlowGoldbergTarjan.h>
#include <ogdf/graphalg/MinSTCut.h>

#define OGDF_STEINERTREE_GOEMANS139_SEPARATE_CONNECTED_COMPONENTS // this is faster
#define OGDF_STEINERTREE_GOEMANS139_SEPARATE_YVAR_CONSTRAINTS // if not defined: generate all yvar constraints in the beginning
#define OGDF_STEINERTREE_GOEMANS139_LOSS_ON_METRIC_CLOSURE // faster and better

namespace ogdf {

/*!
 * \brief This class implements the (1.39+epsilon)-approximation algorithm
 * for the Steiner tree problem by Goemans et. al.
 *
 * @ingroup ga-steiner
 *
 * This implementation is based on:
 *
 * M.X. Goemans, N. Olver, T. Rothvoß, R. Zenklusen:
 * Matroids and Integrality Gaps for Hypergraphic Steiner Tree Relaxations.
 * STOC 2012, pages 1161-1176, 2012
 *
 * and
 *
 * S. Beyer, M. Chimani: Steiner Tree 1.39-Approximation in Practice.
 * MEMICS 2014, LNCS 8934, 60-72, Springer, 2014
 */
template<typename T>
class MinSteinerTreeGoemans139 : public MinSteinerTreeModule<T>
{
private:
	class UFCR;

protected:
	int m_restricted;
	bool m_use2approx;
	bool m_forceAPSP;
	bool m_separateCycles;
	int m_seed;

public:
	MinSteinerTreeGoemans139()
	  : m_restricted(3)
	  , m_use2approx(false)
	  , m_forceAPSP(false)
	  , m_separateCycles(false)
	  , m_seed(1337)
	{
	}

	virtual ~MinSteinerTreeGoemans139() { }

	/*!
	 * \brief Sets the maximal number of terminals in a full component
	 * @param k the maximal number of terminals in a full component
	 */
	void setMaxComponentSize(int k)
	{
		m_restricted = k;
	}

	/*!
	 * \brief Set seed for the random number generation.
	 * @param seed The seed
	 */
	void setSeed(int seed)
	{
		m_seed = seed;
	}

	/*!
	 * \brief Use Takahashi-Matsuyama 2-approximation as upper bounds
	 * \note not recommended to use in general
	 * @param use2approx True to apply the bound
	 */
	void use2Approximation(bool use2approx = true)
	{
		m_use2approx = use2approx;
	}

	/*! \brief Force full APSP algorithm even if consecutive SSSP algorithms may work
	 *
	 * For the 3-restricted case, it is sufficient to compute an SSSP from every terminal
	 *  instead of doing a full APSP. In case a full APSP is faster, use this method.
	 * @param force True to force APSP instead of SSSP.
	 */
	void forceAPSP(bool force = true)
	{
		m_forceAPSP = force;
	}

	/*!
	 * \brief Use stronger LP relaxation (not recommended in general)
	 * @param separateCycles True to turn the stronger LP relaxation on
	 */
	void separateCycles(bool separateCycles = true)
	{
		m_separateCycles = separateCycles;
	}

protected:
	/*!
	 * \brief Builds a minimum Steiner tree for a given weighted graph with terminals \see MinSteinerTreeModule::computeSteinerTree
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
};

//! \brief Class managing LP and LP-based approximation
//! \todo should be refactored, done this way for historical reasons
template<typename T>
class MinSteinerTreeGoemans139<T>::UFCR
{
	const EdgeWeightedGraph<T> &m_G;
	const NodeArray<bool> &m_isTerminal;
	List<node> m_terminals; //!< List of terminals
	List<node> m_nonterminals; //!< List of non-terminals
	steinertree::FullComponentWithExtraStore<T, double> m_fullCompStore; //!< all enumerated full components, with solution

	NodeArray< NodeArray<T> > m_distance;
	NodeArray< NodeArray<edge> > m_predAPSP;

	OsiSolverInterface *m_osiSolver;
	CoinPackedMatrix *m_matrix;
	double *m_objective;
	double *m_lowerBounds;
	double *m_upperBounds;

	int m_restricted;
	enum class Approx2State {
		Off,
		On,
		JustUseIt,
	};
	Approx2State m_use2approx;
	bool m_ssspDistances;
	bool m_separateCycles;
#ifdef OGDF_STEINERTREE_GOEMANS139_SEPARATE_CONNECTED_COMPONENTS
	int m_separationStage;
#endif

	double m_eps; // epsilon for double operations

	std::minstd_rand m_rng;

	EdgeWeightedGraphCopy<T> *m_approx2SteinerTree;
	T m_approx2Weight;

	//! Manage witness set and core edges of a blowup graph
	class CoreWitness
	{
		std::minstd_rand &m_rng;

		List<node> m_coreEdges; //!< the core edges as nodes

		/* witness set (for some component and specified K (core edge set))
		 *
		 *  W(e) = { e } if e is a core edge
		 *  W(e) =
		 *     take the path P of loss edges from u to e in the component
		 *     add every incident core edge of u
		 *  that is: if we root the loss connected component at the terminal,
		 *   W(e) = all core edges incident to the subtree below e
		 *
		 * Need fast lookups for:
		 *  (1) |W(f)|
		 *  (2) all loss edges f such that given core edge e is in W(f)
		 *  (3) all loss edges f such that W(f) is a subset of a given basis -> for each e \in B: do (2)
		 *
		 * Data Structures:
		 *  - EdgeArray<int> witnessCard with witnessCard[e] = |W(e)|
		 *  - NodeArray< List<edge> > witness with witness[v_e] = { f | e \in W(f) }
		 * (Note that core edges are given as nodes.)
		 * Also note that we can save it all much sparser.
		 */
		EdgeArray<int> witnessCard;
		NodeArray< List<edge> > witness;

		//! \brief Add a core edge
		//! \remark Note that core edges are implemented by nodes in the blowup graph
		void addCore(node e)
		{
			m_coreEdges.pushBack(e);
		}

		//! Add e to W(f)
		void addWitness(node e, edge f)
		{
			++witnessCard[f];
			witness[e].pushBack(f);
		}

		//! Compute a random set of core edges
		void computeRandomCoreEdges(const EdgeWeightedGraphCopy<T> &G, const List<node> &terminals, EdgeArray<bool> &isInTree)
		{
			// Let's do Kruskal's algorithm without weights but on a randomly permuted edge list.
			// We virtually contract all terminals in the union-find data structure.
			NodeArray<int> setID(G, -1);
			isInTree.init(G, false);
			DisjointSets<> uf(G.numberOfNodes() - terminals.size() + 1);

			int contractedID = uf.makeSet();
			OGDF_ASSERT(contractedID >= 0);
			for (node t : terminals) {
				setID[t] = contractedID;
			}
			for (node v = G.firstNode(); v; v = v->succ()) {
				if (setID[v] < 0) {
					setID[v] = uf.makeSet();
				}
			}

			// obtain a random edge permutation
			List<edge> edgePermutation;
			for (edge e : G.edges) {
				edgePermutation.pushBack(e);
			}
			edgePermutation.permute(m_rng);

			// add edges if we do not close a cycle
			for (edge e : edgePermutation) {
				const int v = setID[e->source()];
				const int w = setID[e->target()];
				if (uf.find(v) != uf.find(w)) {
					isInTree[e] = true;
					uf.link(uf.find(v), uf.find(w));
				}
			}
		}

	public:
		CoreWitness(std::minstd_rand &rng)
		  : m_rng(rng)
		{
		}

		/** \brief Find a "random" set of core edges and "replace" found edges by nodes,
		  *   also find the witness sets for the core edges
		  * \todo Derandomization: find best set of core edges using dynamic programming
		  * \param blowupGraph the blowup graph
		  * \param capacity the capacities of the edges in the blowup graph
		  * \param terminals the terminals in the blowup graph
		  */
		void init(EdgeWeightedGraphCopy<T> &blowupGraph, EdgeArray<int> &capacity, const List<node> &terminals)
		{
			witnessCard.init(blowupGraph, 0);
			witness.init(blowupGraph);

			// compute set of core edges
			EdgeArray<bool> isLossEdge;
			computeRandomCoreEdges(blowupGraph, terminals, isLossEdge);

			// add nodes for core edges and be able to map them
			EdgeArray<node> splitMap(blowupGraph, nullptr);
			ArrayBuffer<edge> coreEdges;
			for (edge e = blowupGraph.firstEdge(); e; e = e->succ()) {
				if (!isLossEdge[e]) {
					splitMap[e] = blowupGraph.newNode();
					coreEdges.push(e);
				}
			}

			// traverse losses from all terminals to find witness sets
			NodeArray<adjEntry> pred(blowupGraph, nullptr);
			for (node t : terminals) {
				ArrayBuffer<node> stack;
				stack.push(t);
				while (!stack.empty()) {
					// for each node v "below" an edge e in the traversal:
					//   add all incident core edges vw to W(e)
					const node v = stack.popRet();
					for (adjEntry adj : v->adjEntries) {
						const edge e = adj->theEdge();
						const node w = adj->twinNode();
						if (!pred[v] || w != pred[v]->theNode()) { // do not look at backward arcs in the tree
							if (isLossEdge[e]) {
								stack.push(w);
								pred[w] = adj;
							} else {
								for (node x = v; pred[x]; x = pred[x]->theNode()) {
									addWitness(splitMap[e], pred[x]->theEdge());
								}
							}
						}
					}
				}
			}

			// finally replace core edges by nodes
			for (edge e : coreEdges) {
				const T w = blowupGraph.weight(e);
				const node x = splitMap[e];
				const int cap = capacity[e];
				OGDF_ASSERT(x);
				capacity[blowupGraph.newEdge(e->source(), x, w)] = cap;
				capacity[blowupGraph.newEdge(x, e->target(), w)] = cap;
				// the cost of a core edge node is hence the weight of one incident edge; also keep capacity
				blowupGraph.delEdge(e);
				addCore(x);
			}
		}

		//! Copy witness sets and core edges for a given copy map
		void makeCopy(const HashArray<edge,edge> &edgeMap)
		{
			for (HashConstIterator<edge,edge> pair = edgeMap.begin(); pair.valid(); ++pair) {
				const edge eO = pair.key();
				const edge eC = pair.info();
				const node vO = eO->target();
				const node vC = eC->target();
				witnessCard[eC] = witnessCard[eO]; // copy witness cardinality
				if (vC == vO) { // target is a terminal, so it cannot be a core edge
					continue;
				}
				for (ListIterator<node> it = m_coreEdges.begin(); it.valid(); ++it) {
					if (*it == vO) { // vO is a core edge
						m_coreEdges.insertAfter(vC, it); // make vC a core edge

						// copy witness sets
						// XXX: do we need this or are witness sets computed after the loop again?!
						for (edge e : witness[vO]) {
							witness[vC].pushBack(edgeMap[e]);
						}
						break;
					}
				}
			}
		}

		//! Return list of core edges (implemented by nodes)
		const List<node> &core() const
		{
			return m_coreEdges;
		}

		//! \brief Remove a core edge
		//! \note the blowup graph is not affected
		void delCore(node e)
		{
			/* What happens when we remove a core edge?
			 *  - loss edges are not affected
			 *  - we have to remove a core edge e from W(f) for all f, which means:
			 *    for all elements f of witness[v_e], decrease witnessCard[f], then remove witness[v_e]
			 */
			for (edge f : witness[e]) {
				--witnessCard[f];
			}
			// witness[e] is removed by removing the node from the graph
			m_coreEdges.removeFirst(e);
		}

		//! Return the number of witnesses of an edge
		int numberOfWitnesses(edge e) const
		{
			return witnessCard[e];
		}

		//! Return list of loss edges f in W(e)
		const List<edge> &witnessList(node e) const
		{
			return witness[e];
		}
	};

	//! Obtain and provides information about components in a given blowup graph
	class BlowupComponents
	{
	protected:
		// To represent Gamma(X) [the set of all components in the blowup graph], we give
		//  - terminals, source and target the component id 0
		//  - all other nodes a component id > 0. Nodes with same id belong to the same component.
		// Note that this is also fine for 2-components with only one edge, since this
		// is a core edge and hence a dummy node is inserted.
		NodeArray<int> componentId;
		// We also save lists of terminals for each component in the specified array buffer.
		ArrayBuffer< List<node> > componentTerminals;
		// To efficiently traverse the found component, we save the edge from the root of the component
		ArrayBuffer<edge> componentRootEdge;
		// Finally a component -> cost array.
		ArrayBuffer<T> componentCost;

		int maxId; // == the size of the arraybuffers

		//! Initialize all information about the component starting with \a rootEdge in the blowup graph.
		void initializeComponent(edge rootEdge, const EdgeWeightedGraphCopy<T> &blowupGraph, const NodeArray<node> &bgOriginal)
		{
			List<node> queue;
			const node start = rootEdge->target();
			queue.pushBack(start);
			componentRootEdge.push(rootEdge);
			componentTerminals.push(List<node>());
			List<node> &terms = componentTerminals[maxId];
			if (!bgOriginal[start]) { // start node is a core edge
				componentCost.push(blowupGraph.weight(rootEdge));
			} else {
				componentCost.push(0);
			}
			T &cost = componentCost[maxId];
			++maxId;
			while (!queue.empty()) {
				const node v = queue.popBackRet();
				componentId[v] = maxId;
				for (adjEntry adj : v->adjEntries) {
					const node w = adj->twinNode();
					if (componentId[w] < 0) {
						// count coreEdge cost only once
						if (bgOriginal[v]) { // v is no core edge
							cost += blowupGraph.weight(adj->theEdge());
						}
						if (blowupGraph.original(w)) { // is terminal?
							terms.pushBack(w);
						} else {
							queue.pushBack(w);
						}
					}
				}
			}
		}

	public:
		//! Find all components in the blowup graph and initialize information them
		BlowupComponents(const EdgeWeightedGraphCopy<T> &blowupGraph, const List<node> &bgTerminals, const std::vector<node> &ignores, const NodeArray<node> &bgOriginal)
		  : componentId(blowupGraph, -1)
		  , maxId(0)
		{
			for (node v : ignores) {
				componentId[v] = 0;
			}

			for (node t : bgTerminals) {
				for (adjEntry rootAdj : t->adjEntries) {
					const edge rootEdge = rootAdj->theEdge();
					if (rootEdge->source() != t) {
						continue;
					}
					const node r = rootAdj->twinNode();
					OGDF_ASSERT(r == rootEdge->target());
					if (componentId[r] < 0) {
						initializeComponent(rootEdge, blowupGraph, bgOriginal);
					}
				}
			}

			// now set terminals to id 0   [XXX: why dont we keep -1?]
			for (node t : bgTerminals) {
				componentId[t] = 0;
			}
		}

		//! Return list of terminals for a given component
		const List<node> &terminals(int id) const
		{
			OGDF_ASSERT(id > 0);
			return componentTerminals[id-1];
		}

		//! Return the component id a given node is contained in
		int id(node v) const
		{
			return componentId[v];
		}

		//! Return total cost of a given component
		const T &cost(int id) const
		{
			OGDF_ASSERT(id > 0);
			return componentCost[id-1];
		}

		//! Return number of components
		int size() const
		{
			return maxId;
		}

		//! Return the edge coming from the root of a given component
		edge rootEdge(int id) const
		{
			OGDF_ASSERT(id > 0);
			return componentRootEdge[id-1];
		}

		//! Set the edge coming from the root for a given component
		void setRootEdge(int id, edge e) // beware of using!
		{
			OGDF_ASSERT(id > 0);
			componentRootEdge[id-1] = e;
			OGDF_ASSERT(componentTerminals[id-1].search(e->source()).valid());
		}
	};

	//! \name Finding full components
	//! @{

	//! Find full components of size 2
	void findFull2Components(SubsetEnumerator<node> &terminalSubset);
	//! Find full components of size 3
	void findFull3Components(SubsetEnumerator<node> &terminalSubset);
	//! Find full components
	void findFullComponents();

	//! Check if a given (sub)graph is a valid full component
	bool isValidComponent(const EdgeWeightedGraphCopy<T> &graph)
	{
		for (edge e : graph.edges) {
			const node u = graph.original(e->source());
			const node v = graph.original(e->target());
			if (m_predAPSP[u][v] == nullptr) {
				return false;
			}
		}
		for (node v : graph.nodes) {
			if (m_isTerminal[graph.original(v)] // is a terminal
			 && v->degree() > 1) { // but not a leaf
				return false;
			}
		}
		return true;
	}

	//! @}
	//! \name LP model generation and separation
	//! @{

	//! Generate the LP model
	void generateProblem(bool perturb = false);
	//! Generate the LP objective
	void generateObjective(bool perturb);
	//! Add terminal cover constraints to the LP
	void addTerminalCoverConstraint();
	//! Add subset cover constraints to the LP for a given subset of terminals
	bool addSubsetCoverConstraint(const List<node> &subset);
	//! Add constraint that the sum of x_C over all components C spanning terminal \a t is at least 1 to ensure y_t >= 0
	void addYConstraint(const node t);

	//! \brief Separate to ensure that the solution is connected
	//! \return True iff new constraints have been introduced
	bool separateConnected(const ArrayBuffer<int> &activeComponents);

	//! \brief Perform the general cut-based separation algorithm
	//! \return True iff new constraints have been introduced
	bool separateMinCut(const ArrayBuffer<int> &activeComponents);

	//! \brief Perform the separation algorithm for cycle constraints (to obtain stronger LP solution)
	//! \return True iff new constraints have been introduced
	bool separateCycles(const ArrayBuffer<int> &activeComponents);

	//! \brief Perform all available separation algorithms
	//! \return True iff new constraints have been introduced
	bool separate();

	//! \brief Generates an auxiliary multi-graph for separation (during LP solving):
	//!  directed, with special source and target, without Steiner vertices of degree 2
	double generateMinCutSeparationGraph(const ArrayBuffer<int> &activeComponents, node &source, node &target, GraphCopy &G, EdgeArray<double> &capacity, int &cutsFound)
	{
		G.createEmpty(m_G);
		capacity.init(G);
		source = G.newNode();
		for (node t : m_terminals) { // generate all terminals
			G.newNode(t);
		}
		target = G.newNode();
		for (int j = 0; j < activeComponents.size(); ++j) {
			const int i = activeComponents[j];
			const double cap = m_fullCompStore.extra(i);
			const Array<node> &terminals = m_fullCompStore.terminals(i);
			// take the first terminal as root
			// XXX: we may generate parallel edges but it's ok
			const auto it0 = terminals.begin();
			const node rC = G.copy(*it0);
			capacity[G.newEdge(source, rC)] = cap;
			if (terminals.size() > 2) {
				const node v = G.newNode();
				capacity[G.newEdge(rC, v)] = cap;
				for (auto it = it0 + 1; it != terminals.end(); ++it) {
					const node w = G.copy(*it);
					capacity[G.newEdge(v, w)] = cap;
				}
			} else { // exactly two terminals: we do not need the inner Steiner node
				capacity[G.newEdge(rC, G.copy(*terminals.rbegin()))] = cap;
			}
		}
		double y_R = 0;
		// TODO: perhaps better to compute y_v before
		// add edges to target and compute y_R
		for (node t : m_terminals) {
			const node v = G.copy(t);
			OGDF_ASSERT(v);
			// compute y_v, simply the sum of all x_C where C contains v and then - 1
			double y_v(-1);
			for (adjEntry adj : v->adjEntries) {
				if (adj->twinNode() != source) {
					y_v += capacity[adj->theEdge()];
				}
			}

#ifdef OGDF_STEINERTREE_GOEMANS139_SEPARATE_YVAR_CONSTRAINTS
			if (y_v < -m_eps) {
				addYConstraint(t);
				++cutsFound;
			}
			else
#endif
			if (y_v > 0) {
				capacity[G.newEdge(v, target)] = y_v;
				y_R += y_v;
			}
		}
#if 0
		// just for output of blow-up graph
		for (edge e : G.edges) {
			if (G.original(e->source())) {
				cout << " T:" << G.original(e->source());
			} else {
				cout << " " << e->source();
			}
			cout << " -> ";
			if (G.original(e->target())) {
				cout << "T:" << G.original(e->target());
			} else {
				cout << e->target();
			}
			cout << " " << capacity[e] << "\n";
		}
#endif

		return y_R;
	}

	//! @}
	//! \name Preliminaries and preprocessing for the approximation algorithm
	//! @{

	//! Remove inactive components from m_fullCompStore (since we do not need them any longer)
	void removeInactiveComponents()
	{
		// XXX: is it faster to do this backwards? (less copying)
		int k = 0;
		while (k < m_fullCompStore.size()) {
			if (m_fullCompStore.extra(k) > m_eps) {
				++k;
			} else {
				m_fullCompStore.remove(k);
			}
		}
	}

	//! Remove the full components with the given ids
	void removeComponents(ArrayBuffer<int> &ids)
	{
		ids.quicksort();
		for (int i = ids.size() - 1; i >= 0; --i) {
			m_fullCompStore.remove(ids[i]);
		}
	}

	//! Add a full component to the final solution (by changing nonterminals to terminals)
	void addComponent(NodeArray<bool> &isNewTerminal, int id)
	{
		m_fullCompStore.foreachNode(id, m_predAPSP, [&](node v) {
			isNewTerminal[v] = true;
		});
	}

	//! \brief Preprocess LP solution
	//! \pre every terminal is covered with >= 1
	void preprocess(NodeArray<bool> &isNewTerminal)
	{
		Graph H; // a graph where each component is a star
		NodeArray<int> id(H); // ids each center of the star to the component id
		NodeArray<node> copy(m_G, nullptr); // ids orig in m_G -> copy in H

		List<node> centers; // all centers
		for (int i = 0; i < m_fullCompStore.size(); ++i) {
			const node center = H.newNode();
			centers.pushBack(center);
			id[center] = i;

			for (node vG : m_fullCompStore.terminals(i)) {
				node vH = copy[vG];
				if (!vH) {
					vH = H.newNode();
					copy[vG] = vH;
				}
				H.newEdge(vH, center); // target is always center
			}
		}

		// find components to be inserted into the steinerTree and insert them
		ArrayBuffer<int> inactive; // ids of components we insert; we have to remove them from the set of active components afterwards
		bool changed;
		do {
			changed = false;
			ListIterator<node> it2;
			for (ListIterator<node> it = centers.begin(); it.valid(); it = it2) {
				it2 = it.succ();
				node c = *it;
				int innerNodes = 0; // count inner nodes
				for (adjEntry adj : c->adjEntries) {
					innerNodes += (adj->twinNode()->degree() != 1);
				}
				if (innerNodes <= 1) { // this center represents a component to add to steinerTree
					// insert component into steinerTree
					addComponent(isNewTerminal, id[c]);

					// remove center from H (adjacent leaves can remain being isolated nodes)
					inactive.push(id[c]);
					H.delNode(c);
					centers.del(it);

					changed = true;
				}
			}
		} while (changed);

		removeComponents(inactive);
	}

	//! @}
	//! \name The approximation algorithm
	//! @{

	//! \brief Inserts (based on m_distance and m_predAPSP) a shortest path from a component, directed, into
	//! the blowupGraph and sets weights and capacity accordingly
	void insertShortestPathIntoBlowupGraph(EdgeWeightedGraphCopy<T> &blowupGraph, NodeArray<node> &bgOriginal, EdgeArray<int> &capacity, node s, node t, const node vC, const node wC, int cap)
	{
		if (m_ssspDistances
		 && !m_G.isTerminal(s)) {
			swap(s, t);
			OGDF_ASSERT(m_G.isTerminal(s));
		}
		node tC = wC;
		for (edge e = m_predAPSP[s][t]; e; t = e->opposite(t), e = m_predAPSP[s][t]) {
			const node x = e->opposite(t);
			node xC = vC;
			if (x != s) {
				xC = blowupGraph.newNode();
				OGDF_ASSERT(!bgOriginal[xC]);
				bgOriginal[xC] = x;
			}
			capacity[blowupGraph.newEdge(xC, tC, m_G.weight(e))] = cap;
			tC = xC;
		}
	}

	//! \brief Generates a special-purpose blowup graph for gammoid computation:
	//!  directed, with special source and target, with core edges (implemented as nodes)
	int generateGammoidGraph(node &source, node &target, CoreWitness &cw, EdgeWeightedGraphCopy<T> &blowupGraph, NodeArray<node> &bgOriginal, List<node> &bgTerminals, EdgeArray<int> &capacity, int &y_R) const
	{
		List<int> denominators;

		for (int i = 0; i < m_fullCompStore.size(); ++i) {
			OGDF_ASSERT(m_fullCompStore.extra(i) <= 1.0 + m_eps && m_fullCompStore.extra(i) >= m_eps);
			int num, denom;
			Math::getFraction(m_fullCompStore.extra(i), num, denom);
			OGDF_ASSERT(Math::gcd(num, denom) == 1);
			denominators.pushBack(denom);
		}
		int lcm = 1;
		for (int denom : denominators) {
			lcm = Math::lcm(lcm, denom);
		}

		blowupGraph.createEmpty(m_G);
		bgOriginal.init(blowupGraph, nullptr);
		for (node t : m_terminals) { // generate all terminals
			const node v = blowupGraph.newNode(t);
			bgTerminals.pushBack(v);
			bgOriginal[v] = t;
		}
		y_R = 0;
		List<node> roots;
		List<int> rootCap;
		capacity.init(blowupGraph);
		for (int i = 0; i < m_fullCompStore.size(); ++i) {
			int cap = int(lcm * m_fullCompStore.extra(i) + m_eps);
			// do a bfs of the component tree to add *directed* components
			// with the first terminal as root
			adjEntry start = m_fullCompStore.start(i);
			node v = m_fullCompStore.original(start->theNode());
			roots.pushBack(blowupGraph.copy(v));
			rootCap.pushBack(cap);
			List<adjEntry> queueT;
			List<node> queueC;
			queueT.pushBack(start);
			queueC.pushBack(blowupGraph.copy(v));
			while (!queueT.empty()) {
				const adjEntry inAdj = queueT.popFrontRet();
				const node wT = inAdj->twinNode();
				const node vC = queueC.popFrontRet();

				const node wO = m_fullCompStore.original(wT);
				if (m_isTerminal[wO]) {
#ifdef OGDF_STEINERTREE_GOEMANS139_LOSS_ON_METRIC_CLOSURE
					capacity[blowupGraph.newEdge(vC, blowupGraph.copy(wO), m_fullCompStore.graph().weight(inAdj->theEdge()))] = cap;
#else
					insertShortestPathIntoBlowupGraph(blowupGraph, bgOriginal, capacity, m_fullCompStore.original(vT), wO, vC, blowupGraph.copy(wO), cap);
#endif
				} else { // not a terminal
					const node wC = blowupGraph.newNode();
					bgOriginal[wC] = wO;
#ifdef OGDF_STEINERTREE_GOEMANS139_LOSS_ON_METRIC_CLOSURE
					capacity[blowupGraph.newEdge(vC, wC, m_fullCompStore.graph().weight(inAdj->theEdge()))] = cap;
#else
					insertShortestPathIntoBlowupGraph(blowupGraph, bgOriginal, capacity, m_fullCompStore.original(vT), wO, vC, wC, cap);
#endif
					const adjEntry back = inAdj->twin();
					for (adjEntry adj = back->cyclicSucc(); adj != back; adj = adj->cyclicSucc()) {
						queueT.pushBack(adj);
						queueC.pushBack(wC);
					}
				}
			}
		}

		// remove isolated terminals (can exist by preprocessing)
		ListIterator<node> it2;
		for (ListIterator<node> it = bgTerminals.begin(); it.valid(); it = it2) {
			it2 = it.succ();
			if ((*it)->degree() == 0) {
				blowupGraph.delNode(*it);
				bgTerminals.del(it);
			}
		}

		// compute core edges (and replace these edges by nodes)
		// and witness sets
		cw.init(blowupGraph, capacity, bgTerminals);

		// connect source to component roots
		source = blowupGraph.newNode();
		while (!roots.empty()) {
			capacity[blowupGraph.newEdge(source, roots.popFrontRet(), 0)] = rootCap.popFrontRet();
		}

		// connect target
		target = blowupGraph.newNode();
		for (node v : bgTerminals) {
			OGDF_ASSERT(v);
			int y_v = -lcm;
			// compute y_v, the number of components containing v in the blow up graph - N
			// NOTE: for the non-blowup variant this is simply the sum of all x_C where C contains v ... - 1
			for (adjEntry adj : v->adjEntries) {
				if (adj->twinNode() != source) {
					y_v += capacity[adj->theEdge()];
				}
			}
			OGDF_ASSERT(y_v >= 0);

			if (y_v > 0) {
				capacity[blowupGraph.newEdge(v, target, 0)] = y_v;
				y_R += y_v;
			}
		}
		return lcm;
	}

	//! Computes the rank of the gammoid (given by the blowup graph)
	int gammoidGetRank(EdgeWeightedGraphCopy<T> &blowupGraph, EdgeArray<int> &capacity, const node source, const node target) const
	{
		MaxFlowGoldbergTarjan<int> maxFlow(blowupGraph);
		return maxFlow.computeValue(capacity, source, target);
	}

	//! \brief Finds the best component and its maximum-weight basis in the given blowup graph with given core and witness set
	//! \return The component id of the best component
	int findComponentAndMaxBasis(List<std::pair<node,int>> *&maxBasis, EdgeWeightedGraphCopy<T> &blowupGraph, const List<node> &bgTerminals, const CoreWitness &cw, EdgeArray<int> &capacity, const BlowupComponents &gamma, const int N, const int y_R, const node source, const node target)
	{
		// there should always be saturated flow to the component roots
		// (contracted matroid)
		EdgeArray<int> lB(blowupGraph, 0);
		for (adjEntry adj : source->adjEntries) {
			const edge e = adj->theEdge();
			lB[e] = capacity[e];
		}

		// compute weights of core edges and add source->core edges
		List<edge> sourceCoreEdges;
		EdgeArray<double> cost(blowupGraph, 0);
		for (ListConstIterator<node> it = cw.core().rbegin(); it.valid(); --it) { // XXX: why backwards?
			// compute weight of core edge
			const node v = *it;
			edge tmp = v->firstAdj()->theEdge();
			double weight = (double)blowupGraph.weight(tmp);
			for (edge e : cw.witnessList(v)) {
				OGDF_ASSERT(cw.numberOfWitnesses(e) > 0);
				weight += (double)blowupGraph.weight(e) / cw.numberOfWitnesses(e);
			}

			// add edges from source to core edges v
			edge e = blowupGraph.newEdge(source, v, 0);
			capacity[e] = capacity[tmp];
			cost[e] = -weight;
			sourceCoreEdges.pushBack(e);
		}

		maxBasis = nullptr;
		NodeArray<int> supply(blowupGraph, 0);
		EdgeArray<int> flow(blowupGraph);
		MinCostFlowReinelt<double> mcf;

		for (int id = 1; id <= gamma.size(); ++id) {
			/* We want to find maximum-weight basis B \in B^K_Q
			 * B_Q = minimal edge set to remove after contracting Q to obtain feasible solution (blowup graph)
			 *     = bases of gammoid M_Q
			 * B^K_Q = { B \in B_Q | B \subseteq K }  [K is the set of core edges]
			 *       = bases of gammoid M^K_Q
			 * M^K_Q is gammoid "obtained by restricting M_Q to K"
			 * M_Q = M'_Q / X'  (M'_Q contracted by X') is a gammoid
			 * M'_Q = gammoid from   X \cup X'  to  Y  in  D'  with arc-disjointness instead of vertex-disjointness
			 *    D' is D like for separation but without source node and with interim node v_e for each edge e
			 *    X  = inserted nodes v_e in each edge in blowup graph
			 *    X' = the (arbitrary) roots of each component
			 *    Y  = Q \cup {t}
			 * that means:
			 *   I (subset of X) is an independent set of M_Q
			 *  <=> (if X' is always an independent set in M'_Q ... which we assume here)
			 *   I \cup X'  is an independent set of M'_Q
			 * Restricting to K means:
			 *   only put v_e in X with e \in K (that is, we only need to *generate* v_e if e \in K)
			 * Hence:
			 * - we generate D'^K (this is blowupGraph)
			 * - compute the max flow from X^K \cup X' to Q \cup {t}
			 * - ASSERT that X' is "saturated"!
			 * - check which subset of X is saturated -> these are the nodes representing the edge set we need
			 */
			// add edges from component's terminals to target
			List<edge> Q_to_target;
			for (node t : gamma.terminals(id)) {
				const edge e = blowupGraph.newEdge(t, target, 0);
				capacity[e] = N * y_R; // upper bound for capacity
				Q_to_target.pushBack(e);
			}

			List<std::pair<node,int>> *basis = new List<std::pair<node,int>>;

			int rank = gammoidGetRank(blowupGraph, capacity, source, target);
			supply[source] = rank;
			supply[target] = -rank;

			// find maximum weight basis
#ifdef OGDF_DEBUG
			bool feasible = mcf.call(blowupGraph, lB, capacity, cost, supply, flow);
			OGDF_ASSERT(feasible);
			double minCostVal;
			OGDF_ASSERT(mcf.checkComputedFlow(blowupGraph, lB, capacity, cost, supply, flow, minCostVal));
#else
			mcf.call(blowupGraph, lB, capacity, cost, supply, flow);
#endif
			double weight(0);
			for (edge e : sourceCoreEdges) {
				if (flow[e] > 0) {
					basis->pushBack(std::pair<node,int>(e->target(), flow[e]));
					weight -= flow[e] * cost[e];
				}
			}

			// remove temporary edges for multi-target max-flow again
			for (edge e : Q_to_target) {
				blowupGraph.delEdge(e);
			}
			// XXX/TODO: we could also keep target edges from all terminals to the target
			// and just change the capacity

			// we choose the component with max cost*N <= weight
			if (gamma.cost(id) * N <= weight + m_eps) {
				maxBasis = basis;
				for (edge e : sourceCoreEdges) { // clean up (XXX: check if necessary)
					blowupGraph.delEdge(e);
				}
				return id;
			}
			delete basis;
		}
		return 0; // no component chosen, fail
	}

	//! \brief For the end of the algorithm: find cheapest component and choose all remaining core edges as basis
	//! \return The component id of the cheapest component
	int findCheapestComponentAndRemainingBasis(List<std::pair<node,int>> *&maxBasis, const CoreWitness &cw, EdgeArray<int> &capacity, const BlowupComponents &gamma)
	{
		// find cheapest component
		int compId = 0;
		double cost = 0;
		for (int id = 1; id <= gamma.size(); ++id) {
			if (gamma.cost(id) > cost) {
				cost = gamma.cost(id);
				compId = id;
			}
		}
		// use all core edges as basis
		maxBasis = new List<std::pair<node,int>>();
		for (node v : cw.core()) {
			edge tmp = v->lastAdj()->theEdge();
			maxBasis->pushBack(std::pair<node,int>(v, capacity[tmp]));
		}
		return compId;
	}

	//! Add a component of the blowup graph to the final solution (by changing nonterminals to terminals)
	void addComponent(NodeArray<bool> &isNewTerminal, const EdgeWeightedGraphCopy<T> &blowupGraph, const NodeArray<node> &bgOriginal, const edge rootEdge)
	{
		OGDF_ASSERT(blowupGraph.original(rootEdge->source()));
		List<node> stack;
		stack.pushBack(rootEdge->target());
		while (!stack.empty()) {
			const node v = stack.popBackRet();
			if (blowupGraph.original(v)) { // v is terminal
				continue;
			}
			const node vO = bgOriginal[v];
			if (vO) {
				isNewTerminal[vO] = true;
			}
			for (adjEntry adj : v->adjEntries) {
				const node w = adj->theEdge()->target();
				if (v != w) { // outgoing edge
					stack.pushBack(w);
				}
			}
		}
	}

	//! Copy a component in the blowup graph and set original capacity to \a origCap and capacity of copy to \a copyCap
	void copyComponent(EdgeWeightedGraphCopy<T> &blowupGraph, EdgeArray<int> &capacity, CoreWitness &cw, NodeArray<node> &bgOriginal, const edge origEdge, const int origCap, const int copyCap)
	{
		if (copyCap == 0) {
			return;
		}
		List<edge> todo;
		List<node> origin;
		HashArray<edge,edge> edgeMap;
		todo.pushBack(origEdge);
		origin.pushBack(origEdge->source());
		while (!todo.empty()) {
			edge eO = todo.popFrontRet();
			node vC = origin.popFrontRet();
			node wO = eO->target();
			node wC = wO;
			if (!blowupGraph.original(wO)) { // is not a terminal
				wC = blowupGraph.newNode();
				bgOriginal[wC] = bgOriginal[wO]; // copy
			}
			edge eC = blowupGraph.newEdge(vC, wC, blowupGraph.weight(eO));
			capacity[eC] = copyCap;
			capacity[eO] = origCap;
			edgeMap[eO] = eC;
			if (!blowupGraph.original(wO)) { // is not a terminal
				for (adjEntry adj = eO->adjTarget()->cyclicSucc(); adj != eO->adjTarget(); adj = adj->cyclicSucc()) {
					OGDF_ASSERT(adj->theEdge()->target() != eO->target()); // outgoing edges
					origin.pushBack(wC);
					todo.pushBack(adj->theEdge());
				}
			}
		}
		cw.makeCopy(edgeMap);
	}

	//! Remove basis (given by \a v, a core edge node of the maximum basis) and cleanup
	void removeBasisAndCleanup(EdgeWeightedGraphCopy<T> &blowupGraph, BlowupComponents &gamma, node v, const node source, const node pseudotarget, const int N)
	{
		List<node> cleanup;
		cleanup.pushBack(v->firstAdj()->twinNode());
		cleanup.pushBack(v->lastAdj()->twinNode());
		OGDF_ASSERT(v->degree() == 2 && v->firstAdj()->twinNode() != v->lastAdj()->twinNode());
		blowupGraph.delNode(v);

		while (!cleanup.empty()) {
			v = cleanup.popBackRet();
			if (!blowupGraph.original(v)) { // is not a terminal
				OGDF_ASSERT(v->degree() >= 1);
				if (v->degree() == 1) { // v is a pendant node, delete
					cleanup.pushBack(v->firstAdj()->twinNode());
					blowupGraph.delNode(v);
				} else
				if (v->indeg() == 0) { // v has no incoming edge, fix
					const node w = v->firstAdj()->twinNode();
					const edge e = v->firstAdj()->theEdge();
					blowupGraph.reverseEdge(e);
					OGDF_ASSERT(e->source() == w);
					if (!blowupGraph.original(w)) { // is not a terminal
						cleanup.pushBack(w);
						// when w is cleaned up, it must not go back to v
						if (w->firstAdj()->theEdge() == e) {
							// move first adjacency entries of w away (w->v is not first anymore)
							blowupGraph.moveAdjAfter(w->firstAdj(), w->lastAdj());
						}
					} else { // w is a terminal
						// in this case we have just fixed the direction of the component,
						// so w is the new root of the component -> update gamma
						gamma.setRootEdge(gamma.id(v), e);
					}
				}
			}
		}
	}

	//! Remove a given basis and cleanup, the basis may be given fractionally
	void removeFractionalBasisAndCleanup(List<std::pair<node,int>> &basis, EdgeWeightedGraphCopy<T> &blowupGraph, NodeArray<node> &bgOriginal, EdgeArray<int> &capacity, BlowupComponents &gamma, CoreWitness &cw, const node source, const node pseudotarget, const int N)
	{
		// remove B from K (K := K \ B) and from blowup graph (X := X - B)
		// and, while at it, remove cleanup edges from blowup graph (X := X - F)
		// and fix components that have no incoming edges
		List<Prioritized<node, int>> fractionalCoreEdges; // we defer fractional basis elements
		for (std::pair<node,int> p : basis) {
			const node v = p.first;
			const int count = p.second;
			OGDF_ASSERT(v->degree() == 2);
			int origCap = capacity[v->firstAdj()->theEdge()];
			OGDF_ASSERT(count <= origCap);
			if (count < origCap) { // only remove a fraction?
				fractionalCoreEdges.pushBack(Prioritized<node,int>(v, -count));
			} else {
				// we are deleting the core edge from the whole component
				cw.delCore(v);
				removeBasisAndCleanup(blowupGraph, gamma, v, source, pseudotarget, N);
			}
		}
		fractionalCoreEdges.quicksort(); // sort decreasing by flow value
		for (Prioritized<node,int> p : fractionalCoreEdges) {
			const node v = p.item();
			const int count = -p.priority();
			OGDF_ASSERT(v->degree() == 2);
			int origCap = capacity[v->firstAdj()->theEdge()];
			OGDF_ASSERT(count <= origCap);
			// copy (split) the component
			copyComponent(blowupGraph, capacity, cw, bgOriginal, gamma.rootEdge(gamma.id(v)), count, origCap - count);
			// we are deleting the core edge from the whole component
			cw.delCore(v);
			removeBasisAndCleanup(blowupGraph, gamma, v, source, pseudotarget, N);
		}
	}

	//! Update arc capacities s->v and v->t
	int updateSourceAndTargetArcCapacities(EdgeWeightedGraphCopy<T> &blowupGraph, EdgeArray<int> &capacity, const node v, const node source, const node pseudotarget, const int N)
	{
		int delta = 0;
		int capSource = 0;
		int capTarget = -N;
		adjEntry adj2;
		for (adjEntry adj = v->firstAdj(); adj; adj = adj2) {
			adj2 = adj->succ();
			if (adj->twinNode() == source) {
				// remove arcs from the source
				blowupGraph.delEdge(adj->theEdge());
			}
			else
			if (adj->twinNode() == pseudotarget) {
				// remove arcs to the pseudotarget
				delta -= capacity[adj->theEdge()];
				blowupGraph.delEdge(adj->theEdge());
			}
			else {
				// compute y_v for the contraction node
				capTarget += capacity[adj->theEdge()];
				// compute s->v capacity
				if (v != adj->theEdge()->target()) { // outgoing edge
					capSource += capacity[adj->theEdge()];
				}
			}
		}
		OGDF_ASSERT(capTarget >= 0);
		if (capTarget > 0) {
			capacity[blowupGraph.newEdge(v, pseudotarget, 0)] = capTarget;
		}
		if (capSource > 0) {
			capacity[blowupGraph.newEdge(source, v, 0)] = capSource;
		}

		return delta + capTarget;
	}

	//! Perform the actual approximation algorithm on the LP solution
	void doGoemansApproximation(NodeArray<bool> &isNewTerminal);

	//! @}

public:
	//! Initialize all attributes, sort the terminal list
	UFCR(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal)
	  : m_G(G)
	  , m_isTerminal(isTerminal)
	  , m_terminals(terminals)
	  , m_nonterminals()
	  , m_fullCompStore(G, m_terminals, isTerminal)
	  , m_distance()
	  , m_predAPSP()
	  , m_osiSolver(CoinManager::createCorrectOsiSolverInterface())
	  , m_matrix(nullptr)
	  , m_objective(nullptr)
	  , m_lowerBounds(nullptr)
	  , m_upperBounds(nullptr)
	  , m_restricted(3)
	  , m_use2approx(Approx2State::Off)
	  , m_ssspDistances(true)
	  , m_separateCycles(false)
#ifdef OGDF_STEINERTREE_GOEMANS139_SEPARATE_CONNECTED_COMPONENTS
	  , m_separationStage(0)
#endif
	  , m_eps(1e-8)
	  , m_approx2SteinerTree(nullptr)
	{
		steinertree::sortTerminals(m_terminals);
	}

	~UFCR()
	{
		if (m_objective) {
			delete[] m_objective;
		}
		if (m_matrix) {
			delete m_matrix;
		}
		if (m_lowerBounds) {
			delete[] m_lowerBounds;
		}
		if (m_upperBounds) {
			delete[] m_upperBounds;
		}
		delete m_osiSolver;
	}

	//! Perform basic initialization algorithms: find shortest paths, find the full components, generate the LP model
	void init(bool perturb = false)
	{
		if (m_use2approx != Approx2State::Off) { // add upper bound by 2-approximation
			MinSteinerTreeTakahashi<T> mstT;
			m_approx2Weight = mstT.call(m_G, m_terminals, m_isTerminal, m_approx2SteinerTree, m_G.firstNode());
		}

		if (m_restricted > m_terminals.size()) {
			m_restricted = m_terminals.size();
		}

		for (node u : m_G.nodes) {
			if (!m_isTerminal[u]) {
				m_nonterminals.pushBack(u);
			}
		}

		m_distance.init(m_G);
		m_predAPSP.init(m_G);
		if (m_ssspDistances
		 && m_restricted <= 3) {
			// for 2- and 3-restricted computations, it is ok to use SSSP from all terminals
			for (node u : m_terminals) {
#ifndef DETOUR
				MinSteinerTreeModule<T>::singleSourceShortestPathsStrict(m_G, u, m_isTerminal, m_distance[u], m_predAPSP[u]);
#else
				MinSteinerTreeModule<T>::singleSourceShortestPathsDetour(m_G, u, m_isTerminal, m_distance[u], m_predAPSP[u]);
#endif
			}
		} else {
			m_ssspDistances = false;
#ifndef DETOUR
			MinSteinerTreeModule<T>::allPairShortestPathsStrict(m_G, m_isTerminal, m_distance, m_predAPSP);
#else
			MinSteinerTreeModule<T>::allPairShortestPathsDetour(m_G, m_nonterminals, m_distance, m_predAPSP);
#endif
		}

		// fill m_fullCompStore
		findFullComponents();

		generateProblem(perturb);
		addTerminalCoverConstraint();

#ifndef OGDF_STEINERTREE_GOEMANS139_SEPARATE_YVAR_CONSTRAINTS
		for (node t : m_terminals) {
			addYConstraint(t);
		}
#endif
	}

	//! Set the maximum full component size to generate
	void setMaxComponentSize(int restricted)
	{
		m_restricted = restricted;
	}

	//! Define if we want to use a 2-approximation as bound
	void use2Approximation(bool use2approx)
	{
		m_use2approx = use2approx ? Approx2State::On : Approx2State::Off;
	}

	//! Define if we want to use cycle constraints for a stronger LP
	void addCycleConstraints(bool add = true)
	{
		m_separateCycles = add;
	}

	//! Define if APSP shortest path algorithm should be used even if SSSP can be used
	void forceAPSP(bool force = true)
	{
		m_ssspDistances = !force;
	}

	//! Solve the LP
	void solve();

	//! Obtain an (1.39+epsilon)-approximation based on the LP solution
	T getApproximation(EdgeWeightedGraphCopy<T> *&finalSteinerTree, const std::minstd_rand &rng, const bool doPreprocessing = true);
};

template<typename T>
void
MinSteinerTreeGoemans139<T>::UFCR::doGoemansApproximation(NodeArray<bool> &isNewTerminal)
{
	node source;
	node pseudotarget;
	EdgeWeightedGraphCopy<T> blowupGraph;
	int y_R;

	CoreWitness cw(m_rng);
	EdgeArray<int> capacity;
	List<node> bgTerminals;
	NodeArray<node> bgOriginal; // original node of a node in blowupGraph (if contracted: just one of it); if nullptr: core edge or special node like source, pseudotarget or target
	const int N = generateGammoidGraph(source, pseudotarget, cw, blowupGraph, bgOriginal, bgTerminals, capacity, y_R);
	// Since we contract terminal nodes in the blow-up graph, the blowup graph
	// needs its own terminal data structures.
	// However, we do not need bgIsTerminal[]. We can just exploit the fact that
	// blowupGraph is a GraphCopy and the only "copied" nodes are terminals.

	node target = blowupGraph.newNode();
	capacity[blowupGraph.newEdge(pseudotarget, target, 0)] = y_R;

	while (bgTerminals.size() > 1) { // T is not a Steiner tree
		// TODO: maybe we should initially compute the blowup components when we *build*
		//       the gammoid graph and update it on each delEdge/delNode
		BlowupComponents gamma(blowupGraph, bgTerminals, {source, pseudotarget, target}, bgOriginal); // Gamma(X)

		OGDF_ASSERT(isLoopFree(blowupGraph));

		// take a component Q in Gamma(X)
		List<std::pair<node,int>> *maxBasis;
		int compId = y_R > 0
		     ? findComponentAndMaxBasis(maxBasis, blowupGraph, bgTerminals, cw,
		         capacity, gamma, N, y_R, source, target)
		     : findCheapestComponentAndRemainingBasis(maxBasis, cw, capacity, gamma);
		OGDF_ASSERT(compId);

		// add component Q to T
		addComponent(isNewTerminal, blowupGraph, bgOriginal, gamma.rootEdge(compId));

		// remove (maybe fractional) basis and do all the small things necessary for update
		removeFractionalBasisAndCleanup(*maxBasis, blowupGraph, bgOriginal, capacity, gamma, cw, source, pseudotarget, N);
		delete maxBasis;

		// contract (X := X / Q)
		ListConstIterator<node> it0 = gamma.terminals(compId).begin();
		node v = *it0;
		for (ListConstIterator<node> it = it0.succ(); it.valid(); ++it) {
			node t = *it;
			if (v->degree() == 0) {
				std::swap(v, t);
			}
			bgTerminals.removeFirst(t);
			if (t->degree() > 0) {
				v = blowupGraph.contract(blowupGraph.newEdge(v, t, 0));
				// the contract method ensures that capacities, weights, and everything is kept
			} else { // v->degree() == 0
				blowupGraph.delNode(t);
			}
		}

		if (bgTerminals.size() > 1) {
			//y_R += updateSourceAndTargetArcCapacities(blowupGraph, capacity, v, source, pseudotarget, N);
			// update capacities from source to terminals and terminals to pseudotarget
			for (node t : bgTerminals) {
				y_R += updateSourceAndTargetArcCapacities(blowupGraph, capacity, t, source, pseudotarget, N);
			}
			// XXX: doing it for v and all terminals we have met during cleanup would be sufficient
			OGDF_ASSERT(target->degree() == 1);
			capacity[target->firstAdj()->theEdge()] = y_R;
		}
	}
}

template<typename T>
T
MinSteinerTreeGoemans139<T>::UFCR::getApproximation(EdgeWeightedGraphCopy<T> *&finalSteinerTree, const std::minstd_rand &rng, const bool doPreprocessing)
{
	if (m_use2approx == Approx2State::JustUseIt) {
		// no remaining components
		finalSteinerTree = m_approx2SteinerTree;
		return m_approx2Weight;
	}
	m_rng = rng;

	const double *constSol = m_osiSolver->getColSolution();
	const int numberOfColumns = m_osiSolver->getNumCols();

	for (int i = 0; i < numberOfColumns; i++) {
		m_fullCompStore.extra(i) = constSol[i];
	}

	removeInactiveComponents();

	NodeArray<bool> isNewTerminal(m_G, false);
	for (node v : m_terminals) {
		isNewTerminal[v] = true;
	}

	if (doPreprocessing) {
		preprocess(isNewTerminal);
	}

	if (!m_fullCompStore.isEmpty()) {
		doGoemansApproximation(isNewTerminal);
	}

	T cost = steinertree::obtainFinalSteinerTree(m_G, isNewTerminal, m_isTerminal, finalSteinerTree);
	if (m_use2approx != Approx2State::Off) {
		if (m_approx2Weight < cost) {
			delete finalSteinerTree;
			finalSteinerTree = m_approx2SteinerTree;
		} else {
			delete m_approx2SteinerTree;
		}
	}

	return cost;
}

template<typename T>
void
MinSteinerTreeGoemans139<T>::UFCR::generateProblem(bool perturb)
{
	int n = m_fullCompStore.size();

	m_matrix = new CoinPackedMatrix;
	m_matrix->setDimensions(0, n);

	m_lowerBounds = new double[n];
	m_upperBounds = new double[n];
	for (int i = 0; i < n; ++i) {
		m_lowerBounds[i] = 0;
		m_upperBounds[i] = 1;
	}

	generateObjective(perturb);
	m_osiSolver->loadProblem(*m_matrix, m_lowerBounds, m_upperBounds, m_objective, nullptr, nullptr);

	if (m_use2approx != Approx2State::Off) { // add upper bound by 2-approximation
		CoinPackedVector row(m_objective);
		m_osiSolver->addRow(row, 0, m_approx2Weight);
	}
}

template<typename T>
void
MinSteinerTreeGoemans139<T>::UFCR::generateObjective(bool perturb)
{
	int n = m_fullCompStore.size();
	m_objective = new double[n];

	for (int i = 0; i < n; ++i) {
		double w = m_fullCompStore.cost(i);

		if (perturb) {
			w += (rand() % 1000)*1e-9;
		}

		m_objective[i] = w;
	}
}

template<typename T>
void
MinSteinerTreeGoemans139<T>::UFCR::addYConstraint(const node t)
{
	CoinPackedVector row;

	for (int i = 0; i < m_fullCompStore.size(); ++i) {
		if (m_fullCompStore.isTerminal(i, t)) { // component spans terminal
			row.insert(i, 1);
		}
	}

	m_osiSolver->addRow(row, 1, m_osiSolver->getInfinity());
}

// add constraint if necessary and return if necessary
template<typename T>
bool
MinSteinerTreeGoemans139<T>::UFCR::addSubsetCoverConstraint(const List<node> &subset)
{
	CoinPackedVector row;
	double test = 0;

	for (int i = 0; i < m_fullCompStore.size(); ++i) {
		// compute the intersection cardinality (linear time because terminal sets are sorted by index)
		int intersectionCard = 0;
		auto terminals = m_fullCompStore.terminals(i);
		node *it1 = terminals.begin();
		auto it2 = subset.begin();
		while (it1 != terminals.end()
		    && it2.valid()) {
			if ((*it1)->index() < (*it2)->index()) {
				++it1;
			} else
			if ((*it1)->index() > (*it2)->index()) {
				++it2;
			} else { // ==
				++intersectionCard;
				++it1;
				++it2;
			}
		}
		// and use it as a coefficient
		if (intersectionCard > 1) {
			row.insert(i, intersectionCard - 1);
			test += (intersectionCard - 1) * m_fullCompStore.extra(i);
		}
	}
	if (test > (subset.size() - 1)) {
		m_osiSolver->addRow(row, 0, subset.size() - 1);
		return true;
	}
	return false;
}

template<typename T>
void
MinSteinerTreeGoemans139<T>::UFCR::addTerminalCoverConstraint()
{
	// we use the sum over all components C of (|C| - 1) * x_C = |R| - 1 constraint
	// from the paper
	CoinPackedVector row;

	for (int i = 0; i < m_fullCompStore.size(); ++i) {
		row.insert(i, m_fullCompStore.terminals(i).size() - 1);
	}

	int value = m_terminals.size() - 1;
	m_osiSolver->addRow(row, value, value);
}

template<typename T>
void MinSteinerTreeGoemans139<T>::UFCR::findFull2Components(SubsetEnumerator<node> &terminalSubset)
{
	// compute special case of 2-components (shortest paths between terminals)
	for (terminalSubset.begin(2); terminalSubset.valid(); terminalSubset.next()) {
		const node s = terminalSubset[0];
		const node t = terminalSubset[1];
		if (m_predAPSP[s][t]) {
			EdgeWeightedGraphCopy<T> minComp;
			minComp.createEmpty(m_G);
			minComp.newEdge(minComp.newNode(s), minComp.newNode(t), m_distance[s][t]);
			m_fullCompStore.insert(minComp);
		}
	}
}

template<typename T>
void MinSteinerTreeGoemans139<T>::UFCR::findFull3Components(SubsetEnumerator<node> &terminalSubset)
{
	// compute special case of 3-components (3-stars)
	for (terminalSubset.begin(3); terminalSubset.valid(); terminalSubset.next()) {
		const node t0 = terminalSubset[0];
		const node t1 = terminalSubset[1];
		const node t2 = terminalSubset[2];
		OGDF_ASSERT(m_distance[t0][t1] < numeric_limits<T>::max());
		OGDF_ASSERT(m_distance[t0][t2] < numeric_limits<T>::max());
		OGDF_ASSERT(m_distance[t1][t2] < numeric_limits<T>::max());
		T minCost = // bounded by the sum of the 2-components
		   min(m_distance[t0][t1] + m_distance[t0][t2],
		   min(m_distance[t1][t0] + m_distance[t1][t2],
		       m_distance[t2][t0] + m_distance[t2][t1]));
		node minCenter = nullptr;
		for (node center : m_nonterminals) {
			OGDF_ASSERT(m_distance[t0][center] < numeric_limits<T>::max());
			OGDF_ASSERT(m_distance[t1][center] < numeric_limits<T>::max());
			OGDF_ASSERT(m_distance[t2][center] < numeric_limits<T>::max());
			const T cost = m_distance[t0][center] + m_distance[t1][center] + m_distance[t2][center];
			if (cost < minCost) {
				minCenter = center;
				minCost = cost;
			}
		}
		if (minCenter
		 && m_predAPSP[t0][minCenter]
		 && m_predAPSP[t1][minCenter]
		 && m_predAPSP[t2][minCenter]) {
			// create a full 3-component
			EdgeWeightedGraphCopy<T> minComp;
			minComp.createEmpty(m_G);
			node minCenterC = minComp.newNode(minCenter);
			minComp.newEdge(minComp.newNode(t0), minCenterC, m_distance[t0][minCenter]);
			minComp.newEdge(minComp.newNode(t1), minCenterC, m_distance[t1][minCenter]);
			minComp.newEdge(minComp.newNode(t2), minCenterC, m_distance[t2][minCenter]);
			m_fullCompStore.insert(minComp);
		}
	}
}

template<typename T>
void MinSteinerTreeGoemans139<T>::UFCR::findFullComponents()
{
	SubsetEnumerator<node> terminalSubset(m_terminals);
	if (m_restricted >= 4) { // use Dreyfus-Wagner based full component generation
		steinertree::DreyfusWagnerFullComponentGenerator<T> fcg(m_G, m_terminals, m_distance);
		fcg.call(m_restricted);
		for (terminalSubset.begin(2, m_restricted); terminalSubset.valid(); terminalSubset.next()) {
			EdgeWeightedGraphCopy<T> component;
			List<node> terminals;
			terminalSubset.list(terminals);
			fcg.getSteinerTreeFor(terminals, component);
			if (isValidComponent(component)) {
				m_fullCompStore.insert(component);
			}
		}
	} else {
		findFull2Components(terminalSubset);
		if (m_restricted == 3) {
			findFull3Components(terminalSubset);
		}
	}
}

template<typename T>
void
MinSteinerTreeGoemans139<T>::UFCR::solve()
{
	bool initialIteration = true;

	do {
		if (initialIteration) {
			m_osiSolver->initialSolve();
			initialIteration = false;
		} else {
			m_osiSolver->resolve();
		}

		if (!m_osiSolver->isProvenOptimal()) {
			// infeasible:
			if (m_use2approx != Approx2State::Off) {
				m_use2approx = Approx2State::JustUseIt;
				break;
			} else {
				cerr << "Failed to optimize LP!" << endl;
				throw(-1);
			}
		}
	} while (separate());
}

template<typename T>
bool
MinSteinerTreeGoemans139<T>::UFCR::separate()
{
	const double *constSol = m_osiSolver->getColSolution();

	ArrayBuffer<int> activeComponents;
	for (int i = 0; i < m_fullCompStore.size(); ++i) {
		m_fullCompStore.extra(i) = constSol[i];
		if (m_fullCompStore.extra(i) > m_eps) {
			activeComponents.push(i);
		}
	}

#ifdef OGDF_STEINERTREE_GOEMANS139_SEPARATE_CONNECTED_COMPONENTS
	if (!m_separationStage) {
		if (separateConnected(activeComponents)) {
			return true;
		}
		m_separationStage = 1;
	}
#endif
	if (separateMinCut(activeComponents)) {
		return true;
	}
	return (m_separateCycles ? separateCycles(activeComponents) : false);
}

template<typename T>
bool
MinSteinerTreeGoemans139<T>::UFCR::separateConnected(const ArrayBuffer<int> &activeComponents)
{
	NodeArray<int> setID(m_G, -1); // XXX: NodeArray over terminals only would be better
	DisjointSets<> uf(m_terminals.size());
	for (node t : m_terminals) {
		setID[t] = uf.makeSet();
	}

	// union all nodes within one component
	for (int j = 0; j < activeComponents.size(); ++j) {
		auto terminals = m_fullCompStore.terminals(activeComponents[j]);
		auto it = terminals.begin();
		const int s1 = setID[*it];
		for (++it; it != terminals.end(); ++it) {
			uf.link(uf.find(s1), uf.find(setID[*it]));
		}
	}

	if (uf.getNumberOfSets() == 1) { // solution is connected
		return false;
	}

	Array< List<node> > components(m_terminals.size());
	List<int> usedComp;
	for (node t : m_terminals) {
		const int k = uf.find(setID[t]);
		if (components[k].empty()) {
			usedComp.pushBack(k);
		}
		components[k].pushBack(t);
	}
	int cutsFound = 0;
	for (ListConstIterator<int> it = usedComp.begin(); it.valid(); ++it) {
		const int k = *it;
		cutsFound += addSubsetCoverConstraint(components[k]);
	}
	return true;
}

template<typename T>
bool
MinSteinerTreeGoemans139<T>::UFCR::separateMinCut(const ArrayBuffer<int> &activeComponents)
{
	int cutsFound = 0;
	node source;
	node pseudotarget;
	GraphCopy auxG;
	EdgeArray<double> capacity;
	double y_R = generateMinCutSeparationGraph(activeComponents, source, pseudotarget, auxG, capacity, cutsFound);

#ifdef OGDF_STEINERTREE_GOEMANS139_SEPARATE_YVAR_CONSTRAINTS
	if (cutsFound > 0) {
		return true;
	}
#endif

	node target = auxG.newNode();
	capacity[auxG.newEdge(pseudotarget, target)] = y_R;

	EdgeArray<double> flow;
	MaxFlowGoldbergTarjan<double> maxFlow;
	MinSTCut<double> minSTCut;
	for (node t : m_terminals) {
		const node v = auxG.copy(t);

		edge v_to_target = auxG.newEdge(v, target);
		capacity[v_to_target] = numeric_limits<double>::max(); // XXX: smaller is better

		maxFlow.init(auxG, &flow);

		const double cutVal = maxFlow.computeValue(capacity, source, target);
		if (cutVal - y_R < 1 - m_eps) {
			minSTCut.call(capacity, flow, source, target);
			List<node> subset;
			for (node tOrig : m_terminals) {
				const node tCopy = auxG.copy(tOrig);
				if (tCopy && minSTCut.isInBackCut(tCopy)) {
					subset.pushBack(tOrig);
				}
			}

			cutsFound += addSubsetCoverConstraint(subset);
		}

		auxG.delEdge(v_to_target);
	}
	return cutsFound != 0;
}

template<typename T>
bool
MinSteinerTreeGoemans139<T>::UFCR::separateCycles(const ArrayBuffer<int> &activeComponents)
{
	int count = 0;

	// generate auxiliary graph
	Graph G;
	NodeArray<int> id(G);
	for (int i : activeComponents) {
		id[G.newNode()] = i;
	}
	for (node u1 : G.nodes) {
		const int i1 = id[u1];
		const Array<node> &terminals1 = m_fullCompStore.terminals(i1);
		for (node u2 = u1->succ(); u2; u2 = u2->succ()) {
			const int i2 = id[u2];
			const Array<node> &terminals2 = m_fullCompStore.terminals(i2);
			// compute intersection cardinality (linear time because terminal sets are sorted by index)
			int intersectionCard = 0;
			const node *it1 = terminals1.begin();
			const node *it2 = terminals2.begin();
			while (it1 != terminals1.end()
			    && it2 != terminals2.end()) {
				if ((*it1)->index() < (*it2)->index()) {
					++it1;
				} else
				if ((*it1)->index() > (*it2)->index()) {
					++it2;
				} else { // ==
					++intersectionCard;
					if (intersectionCard == 2) {
						G.newEdge(u1, u2);
						break;
					}
					++it1;
					++it2;
				}
			}
		}
	}

	if (G.numberOfEdges() == 0) {
		return false;
	}

	const int maxClique = m_restricted + 1;
	// now find cliques
	Array<List<node>> degrees(maxClique);
	for (node v : G.nodes) {
		int idx = v->degree();
		if (idx == 0) { // ignore isolated nodes
			continue;
		}
		if (idx >= maxClique) {
			idx = maxClique - 1;
		}
		--idx;
		degrees[idx].pushBack(v);
	}
	NodeArray<bool> test(G, false);
	for (int k = degrees.size(); k >= 2; --k) {
		degrees[k-2].conc(degrees[k-1]);
		if (degrees[k-2].size() >= k) {
			SubsetEnumerator<node> nodeSubset(degrees[k-2]);
			for (nodeSubset.begin(k); nodeSubset.valid(); nodeSubset.next()) {
				int countEdges = (k * (k-1)) / 2;
				for (int j = 0; j < nodeSubset.size(); ++j) {
					test[nodeSubset[j]] = true;
				}
				for (edge e : G.edges) {
					if (test[e->source()]
					 && test[e->target()]) {
						countEdges -= 1;
					}
				}
				OGDF_ASSERT(countEdges >= 0);
				if (countEdges == 0) {
					// found clique, add constraint
					double val(0);
					CoinPackedVector row;

					for (int j = 0; j < nodeSubset.size(); ++j) {
						int i = id[nodeSubset[j]];
						val += m_fullCompStore.extra(i);
						row.insert(i, 1);
					}
					if (val >= 1 + m_eps) {
						m_osiSolver->addRow(row, 0, 1);
						++count;
					}
				}
				for (int j = 0; j < nodeSubset.size(); ++j) {
					test[nodeSubset[j]] = false;
				}
			}
		}
	}
	return (count > 0);
}

template<typename T>
T MinSteinerTreeGoemans139<T>::computeSteinerTree(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree)
{
	std::minstd_rand rng(m_seed);
	UFCR ufcr(G, terminals, isTerminal);
	ufcr.setMaxComponentSize(m_restricted);
	ufcr.use2Approximation(m_use2approx);
	ufcr.addCycleConstraints(m_separateCycles);
	if (m_forceAPSP) {
		ufcr.forceAPSP();
	}
	ufcr.init(false);
	ufcr.solve();
	return ufcr.getApproximation(finalSteinerTree, rng, true);
}

} // end namespace
