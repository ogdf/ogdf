/** \file
 * \brief Template of base class of node coloring algorithms.
 *
 * \author Jan-Niklas Buckow
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeSet.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This is the base class.
 */
class NodeColoringModule {
public:
	/**
	 * Data type of the node colors
	 */
	using NodeColor = unsigned int;

	/**
	 * Declares the search procedures
	 */
	enum class SearchProcedure {
		linearSearch, ///< Use a linear search
		binarySearch, ///< Use a binary search
		wigdersonSearch ///< Use the search procedure specified by Wigderson
	};

	/**
	 * Declares the procedure of finding nodes
	 * in Ramsey's algorithm.
	 */
	enum class RamseyProcedure {
		smallestIndex, ///< Uses the node with the smallest index
		smallestDegree, ///< Uses the node with the smallest degree
		largestDegree, ///< Uses the node with the largest degree
		extremalDegree, ///< Uses the node with the extremal degree
	};

	/**
	 * Default constructor.
	 * Initializes the Ramsey-procedure with the smallest index procedure.
	 */
	NodeColoringModule() : m_ramseyProcedure(RamseyProcedure::smallestIndex) { }

	/**
	 * The actual algorithm call.
	 *
	 * @param graph The graph for which the node coloring will be calculated.
	 * @param colors The resulting node coloring.
	 * @param start The first color index which will be used.
	 * @return The number of colors used by the node coloring.
	 */
	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors, NodeColor start = 0) = 0;

	/**
	 * Destructor.
	 */
	virtual ~NodeColoringModule() { }

	/**
	 * Checks if the given node coloring is valid.
	 *
	 * @param graph The graph corresponding to the given coloring.
	 * @param colors The given node coloring.
	 * @return True, iff the coloring is valid, false otherwise.
	 */
	virtual bool checkColoring(const Graph& graph, const NodeArray<NodeColor>& colors) const {
		OGDF_ASSERT(colors.graphOf() == &graph);
		for (edge e : graph.edges) {
			if ((colors[e->source()] == colors[e->target()]) && !e->isSelfLoop()) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Preprocesses a graph so that a coloring algorithm can be applied.
	 * The preprocessing removes self-loop edges and redundant edges.
	 *
	 * @param graph The graph to be preprocessed.
	 */
	virtual void preprocessGraph(Graph& graph) const { makeSimpleUndirected(graph); }

protected:
	/// The RamseyProcedure to be used to select nodes.
	RamseyProcedure m_ramseyProcedure;

	/**
	 * Struct to iterate over all node subsets of a given size.
	 */
	struct SubsetIterator {
		/**
		 * Data type for the indices.
		 */
		using Index = int;

		Array<node>& m_set;
		Array<Index> m_indices;
		Index m_numElements;

		/**
		 * Creates the SubsetIterator with a given set to iterate
		 * and the size of the subsets.
		 *
		 * @param set The given node set as array
		 * @param subsetSize The size of each subset
		 */
		SubsetIterator(Array<node>& set, Index subsetSize)
			: m_set(set), m_indices(subsetSize), m_numElements(subsetSize) {
			OGDF_ASSERT(subsetSize <= set.size());
			for (Index i = 0; i < m_numElements; i++) {
				m_indices[i] = i;
			}
		}

		/**
		 * Returns if the iterator, i.e. the next subset to extract is ok.
		 *
		 * @return true, iff the iterator is ok
		 */
		bool isOk() {
			for (Index i = 0; i < m_numElements; i++) {
				if (m_indices[i] > (m_set.size() - m_numElements + i)) {
					return false;
				}
			}
			return true;
		}

		/**
		 * Returns the current subset.
		 *
		 * @return The current subset
		 */
		List<node> currentSubset() {
			List<node> result;
			for (Index i = 0; i < m_numElements; i++) {
				result.emplaceBack(m_set[m_indices[i]]);
			}
			return result;
		}

		/**
		 * Advances the iterator so that the next subset can
		 * be queried.
		 *
		 * @return true, iff the advancing was successful
		 */
		bool advance() {
			for (Index i = m_numElements - 1; i >= 0; i--) {
				m_indices[i]++;
				if (m_indices[i] <= (m_set.size() + i - m_numElements)) {
					for (int j = i + 1; j < m_numElements; j++) {
						m_indices[j] = m_indices[i] + j - i;
					}
					return true;
				}
			}
			return false;
		}
	};

	/**
	 * Wraps the search for the minimum parameter k so that
	 * the same code can be reused for all algorithms.
	 */
	struct SearchWrapper {
		/**
		 * Performs a step in the search procedure.
		 * It gives feedback if the search was successful or not
		 * with the given parameter.
		 * @param k The given parameter
		 * @return True, iff the search was successful
		 */
		virtual bool step(int k) = 0;
	};

	/**
	 * Checks if subgraph induced by the given nodes forms an independent set.
	 *
	 * @tparam CONTAINER The type of container for the node set
	 * @param graph The graph
	 * @param nodes The subset of nodes
	 * @return true, iff the subset is independent
	 */
	template<class CONTAINER>
	bool checkIndependentSet(const Graph& graph, const CONTAINER& nodes) const {
		NodeSet<false> nodeSet(graph);
		for (node v : nodes) {
			nodeSet.insert(v);
		}
		for (node v : nodes) {
			for (adjEntry adj : v->adjEntries) {
				if (nodeSet.isMember(adj->twinNode()) && adj->twinNode() != v) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Creates a partitioning of the node set into buckets
	 * of a given size.
	 *
	 * @param graph The input graph
	 * @param size The bucket size
	 * @param buckets The resulting bucket partitioning of the nodes
	 */
	virtual void createBuckets(const Graph& graph, int size, Array<Array<node>>& buckets) const {
		OGDF_ASSERT(size > 0);
		OGDF_ASSERT(size <= graph.numberOfNodes());
		unsigned int numBuckets = std::floor(graph.numberOfNodes() / size);
		buckets = Array<Array<node>>(numBuckets);
		Array<unsigned int> num_elems(numBuckets);
		for (auto& bucketElems : num_elems) {
			bucketElems = 0;
		}
		// Prepare the bucket sizes
		for (unsigned int bucketIdx = 0; bucketIdx < numBuckets - 1; bucketIdx++) {
			buckets[bucketIdx] = Array<node>(size);
		}
		buckets[numBuckets - 1] = Array<node>(size + graph.numberOfNodes() % size);
		// Fill the buckets
		unsigned int nodeIdx = 0;
		for (auto v : graph.nodes) {
			unsigned int bucketIdx = min(numBuckets - 1, nodeIdx++ / size);
			buckets[bucketIdx][num_elems[bucketIdx]++] = v;
		}
	}

	/**
	 * Calculates the set of neighbors Y for a given set of nodes X.
	 * Thereby, Y is the union of all neighbors for each node x in X.
	 *
	 * @tparam LISTITERATOR The type of iterator for the node set
	 * @param graph The input graph
	 * @param nodes The nodes of whose neighbors will be determined
	 * @param neighbors The resulting list of neighbors
	 */
	template<class LISTITERATOR>
	void getNeighbors(const Graph& graph, LISTITERATOR nodes, List<node>& neighbors) const {
		neighbors.clear();
		NodeArray<bool> isInserted(graph, false);

		for (LISTITERATOR its = nodes; its.valid(); its++) {
			node v = (*its);
			OGDF_ASSERT(v != nullptr);
			OGDF_ASSERT(v->graphOf() == &graph);

			for (adjEntry adj : v->adjEntries) {
				node twin = adj->twinNode();
				if (!isInserted[twin]) {
					neighbors.emplaceBack(twin);
					isInserted[twin] = true;
				}
			}
		}
	}

	/**
	 * Calculates the sum of neighbor nodes degrees
	 * for a given node v.
	 *
	 * @param v The given node
	 * @return The sum of the neighbor's degrees
	 */
	int getNeighborDegrees(const node& v) const {
		int result = 0;
		for (adjEntry adj : v->adjEntries) {
			node twin = adj->twinNode();
			result += twin->degree();
		}
		return result;
	}

	/**
	 * Calculates the set of complement neighbors Y for a given set of nodes X.
	 * Thereby, Y contains every node which is not in X and is not a neighbor
	 * of a node in X.
	 *
	 * @tparam LISTITERATOR The type of iterator for the node set
	 * @param graph The input graph
	 * @param nodes The nodes of whose complement neighbors will be determined
	 * @param complementNeighbors The resulting list of complement neighbors
	 */
	template<class LISTITERATOR>
	void getNeighborsComplement(const Graph& graph, LISTITERATOR nodes,
			List<node>& complementNeighbors) const {
		// Preparations step
		complementNeighbors.clear();
		NodeArray<bool> isComplement(graph, true);

		// Mark all given nodes as false
		for (LISTITERATOR its = nodes; its.valid(); its++) {
			node v = (*its);
			OGDF_ASSERT(v != nullptr);
			OGDF_ASSERT(v->graphOf() == &graph);
			isComplement[v] = false;
		}

		// Mark all neighbors of the given nodes as false
		List<node> neighbors;
		getNeighbors<LISTITERATOR>(graph, nodes, neighbors);
		for (node& v : neighbors) {
			isComplement[v] = false;
		}

		// Store all nodes which are part of the complement
		for (node v : graph.nodes) {
			if (isComplement[v]) {
				complementNeighbors.emplaceBack(v);
			}
		}
	}

	/**
	 * Merges two lists of nodes and deletes the duplicates.
	 *
	 * @tparam LISTITERATOR Type of list iterator
	 * @param graph The corresponding graph
	 * @param firstList Iterator of the first list
	 * @param secondList Iterator of the second list
	 * @param mergedList Reference to the resulting merged list
	 */
	template<class LISTITERATOR>
	void mergeNodeLists(const Graph& graph, LISTITERATOR firstList, LISTITERATOR secondList,
			List<node>& mergedList) const {
		mergedList.clear();
		NodeArray<bool> isInserted(graph, false);
		List<LISTITERATOR> iterators = {firstList, secondList};
		for (auto& iterator : iterators) {
			for (LISTITERATOR its = iterator; its.valid(); its++) {
				node v = (*its);
				OGDF_ASSERT(v != nullptr);
				OGDF_ASSERT(v->graphOf() == &graph);
				if (!isInserted[v]) {
					mergedList.emplaceBack(v);
					isInserted[v] = true;
				}
			}
		}
	}

	/**
	 * Searches for a maximum degree node in the graph.
	 *
	 * @param graph The input graph.
	 * @param maxDegreeNode The resulting maximum degree node.
	 * @return The maximum degree of the graph.
	 */
	virtual int getMaximumDegreeNode(const Graph& graph, node& maxDegreeNode) const {
		OGDF_ASSERT(!graph.empty());
		maxDegreeNode = graph.firstNode();
		int maxDegree = maxDegreeNode->degree();
		for (node v : graph.nodes) {
			int degree = v->degree();
			if (degree > maxDegree) {
				maxDegree = degree;
				maxDegreeNode = v;
			}
		}
		return maxDegree;
	}

	/**
	 * Searches for all nodes with maximum degree in the graph.
	 *
	 * @param graph The input graph.
	 * @param maxDegreeNodes List of all nodes with maximum degree.
	 * @return The maximum degree of the graph.
	 */
	virtual int getMaximumDegreeNodes(const Graph& graph, List<node>& maxDegreeNodes) const {
		OGDF_ASSERT(!graph.empty());
		maxDegreeNodes.clear();
		maxDegreeNodes.pushBack(graph.firstNode());
		int maxDegree = graph.firstNode()->degree();
		for (node v : graph.nodes) {
			int degree = v->degree();
			if (degree > maxDegree) {
				maxDegree = degree;
				maxDegreeNodes.clear();
				maxDegreeNodes.pushBack(v);
			} else if (degree == maxDegree) {
				maxDegreeNodes.pushBack(v);
			}
		}
		return maxDegree;
	}

	/**
	 * Searches for a minimum degree node in the graph.
	 *
	 * @param graph The input graph.
	 * @param minDegreeNode The resulting minimum degree node.
	 * @return The minimum degree of the graph.
	 */
	virtual int getMinimumDegreeNode(const Graph& graph, node& minDegreeNode) const {
		OGDF_ASSERT(!graph.empty());
		minDegreeNode = graph.firstNode();
		int minDegree = minDegreeNode->degree();
		for (node v : graph.nodes) {
			int degree = v->degree();
			if (degree < minDegree) {
				minDegree = degree;
				minDegreeNode = v;
			}
		}
		return minDegree;
	}

	/**
	 * Searches for all nodes with minimum degree in the graph.
	 *
	 * @param graph The input graph.
	 * @param minDegreeNodes List of all nodes with minimum degree.
	 * @return The minimum degree of the graph.
	 */
	virtual int getMinimumDegreeNodes(const Graph& graph, List<node>& minDegreeNodes) const {
		OGDF_ASSERT(!graph.empty());
		minDegreeNodes.clear();
		minDegreeNodes.pushBack(graph.firstNode());
		int minDegree = graph.firstNode()->degree();
		for (node v : graph.nodes) {
			int degree = v->degree();
			if (degree < minDegree) {
				minDegree = degree;
				minDegreeNodes.clear();
				minDegreeNodes.pushBack(v);
			} else if (degree == minDegree) {
				minDegreeNodes.pushBack(v);
			}
		}
		return minDegree;
	}

	/**
	 * Calculates the maximum node color index used in a certain node coloring.
	 * @param colors The node coloring.
	 * @return Maximum node color index.
	 */
	virtual NodeColor getMaximumNodeColor(NodeArray<NodeColor>& colors) {
		auto maxColor = NodeColor(0);
		for (node v : colors.graphOf()->nodes) {
			if (colors[v] > maxColor) {
				maxColor = colors[v];
			}
		}
		return maxColor;
	}

	/**
	 * Reverses the mapping between the nodes sets of a graph
	 * and a subgraph.
	 *
	 * @param graphOrig The original graph
	 * @param graphNew The new (sub-)graph
	 * @param orig2New Existing node table from the original graph to the new graph
	 * @param new2Orig Resulting node table from the new graph to the original graph
	 */
	virtual void reverseNodeTable(const Graph& graphOrig, const Graph& graphNew,
			const NodeArray<node>& orig2New, NodeArray<node>& new2Orig) const {
		OGDF_ASSERT(&graphOrig == orig2New.graphOf());
		new2Orig = NodeArray<node>(graphNew);
		for (node v : graphOrig.nodes) {
			if (orig2New[v]) {
				new2Orig[orig2New[v]] = v;
			}
		}
	}

	/**
	 * Performs the Ramsey algorithm for finding heuristically large
	 * cliques and independents sets in a graph.
	 *
	 * @param graph The input graph
	 * @param clique List with the resulting clique
	 * @param independentSet List with the resulting independent set
	 */
	virtual void ramseyAlgorithm(const Graph& graph, List<node>& clique,
			List<node>& independentSet) const {
		// Preparation
		clique.clear();
		independentSet.clear();

		// End of recursion
		if (graph.empty()) {
			return;
		}

		// Select arbitrary node
		node v = graph.firstNode();
		if (m_ramseyProcedure == RamseyProcedure::smallestDegree) {
			getMinimumDegreeNode(graph, v);
		} else if (m_ramseyProcedure == RamseyProcedure::largestDegree) {
			getMaximumDegreeNode(graph, v);
		} else if (m_ramseyProcedure == RamseyProcedure::extremalDegree) {
			int n = graph.numberOfNodes();
			int extremalDegree = n / 2;
			v = graph.firstNode();
			for (node w : graph.nodes) {
				int degree = w->degree();
				int extremal = min(degree, n - degree);
				if (extremal < extremalDegree) {
					extremalDegree = extremal;
					v = w;
				}
			}
		}
		List<node> nodes = {v};

		// Recursive step with neighbors
		List<node> neighbors;
		getNeighbors<ListIterator<node>>(graph, nodes.begin(), neighbors);
		Graph subGraphNeighbors;
		NodeArray<node> nodeTableOrig2NewNeighbors(graph, nullptr);
		EdgeArray<edge> edgeTableOrig2NewNeighbors(graph, nullptr);
		subGraphNeighbors.insert(neighbors, graph.edges, nodeTableOrig2NewNeighbors,
				edgeTableOrig2NewNeighbors);
		List<node> cliqueNeighbors;
		List<node> independentSetNeighbors;
		ramseyAlgorithm(subGraphNeighbors, cliqueNeighbors, independentSetNeighbors);

		// Recursive step with complement-neighbors
		List<node> complementNeighbors;
		getNeighborsComplement<ListIterator<node>>(graph, nodes.begin(), complementNeighbors);
		Graph subGraphComplement;
		NodeArray<node> nodeTableOrig2NewComplement(graph, nullptr);
		EdgeArray<edge> edgeTableOrig2NewComplement(graph, nullptr);
		subGraphComplement.insert(complementNeighbors, graph.edges, nodeTableOrig2NewComplement,
				edgeTableOrig2NewComplement);
		List<node> cliqueComplement;
		List<node> independentSetComplement;
		ramseyAlgorithm(subGraphComplement, cliqueComplement, independentSetComplement);

		// Reverse the node tables
		NodeArray<node> nodeTableNew2OrigNeighbors;
		reverseNodeTable(graph, subGraphNeighbors, nodeTableOrig2NewNeighbors,
				nodeTableNew2OrigNeighbors);
		NodeArray<node> nodeTableNew2OrigComplement;
		reverseNodeTable(graph, subGraphComplement, nodeTableOrig2NewComplement,
				nodeTableNew2OrigComplement);

		// Select the biggest cliques / independent sets
		if (cliqueNeighbors.size() + 1 >= cliqueComplement.size()) {
			clique.emplaceBack(v);
			for (node& w : cliqueNeighbors) {
				clique.emplaceBack(nodeTableNew2OrigNeighbors[w]);
			}
		} else {
			for (node& w : cliqueComplement) {
				clique.emplaceBack(nodeTableNew2OrigComplement[w]);
			}
		}
		if (independentSetNeighbors.size() >= independentSetComplement.size() + 1) {
			for (node& w : independentSetNeighbors) {
				independentSet.emplaceBack(nodeTableNew2OrigNeighbors[w]);
			}
		} else {
			independentSet.emplaceBack(v);
			for (node& w : independentSetComplement) {
				independentSet.emplaceBack(nodeTableNew2OrigComplement[w]);
			}
		}
	}

	/**
	 * Applies the clique removal algorithm for finding a large
	 * independent set in a given graph.
	 *
	 * @param graph The given graph
	 * @param independentSet The resulting independent set
	 */
	virtual void cliqueRemoval(const Graph& graph, List<node>& independentSet) const {
		// Copy the graph
		GraphCopy graphMain = GraphCopy(graph);
		List<List<node>> multipleIndependentSets;

		// Main loop: Delete large cliques until the graph is empty
		do {
			List<node> localIndependentSet;
			List<node> localClique;
			ramseyAlgorithm(graphMain, localClique, localIndependentSet);
			for (node& v : localIndependentSet) {
				v = graphMain.original(v);
			}
			multipleIndependentSets.emplaceBack(localIndependentSet);
			for (node v : localClique) {
				graphMain.delNode(v);
			}
		} while (!graphMain.empty());

		// Search for the biggest independent set
		independentSet.clear();
		for (List<node>& set : multipleIndependentSets) {
			if (set.size() > independentSet.size()) {
				independentSet = set;
			}
		}
		OGDF_ASSERT(checkIndependentSet(graph, independentSet));
	}

	/**
	 * Performs a linear search in a specified range with a given oracle.
	 * The oracle tells if the current value is valid.
	 * @param searchWrapper The oracle
	 * @param start Start interval of the range
	 * @param end End interval of the range
	 * @return The smallest valid parameter determined
	 */
	int searchLinear(SearchWrapper* searchWrapper, int start, int end) const {
		OGDF_ASSERT(start <= end);
		for (int k = start; k <= end; k++) {
			if (searchWrapper->step(k)) {
				return k;
			}
		}
		return end;
	}

	/**
	 * Performs a binary search in a specified range with a given oracle.
	 * The oracle tells if the current value is valid.
	 * @param searchWrapper The oracle
	 * @param start Start interval of the range
	 * @param end End interval of the range
	 * @return The smallest valid parameter determined
	 */
	int searchBinary(SearchWrapper* searchWrapper, int start, int end) const {
		OGDF_ASSERT(start <= end);
		int left = start;
		int right = end;
		int middle;
		while (left < right) {
			middle = (left + right) / 2;
			if (searchWrapper->step(middle)) {
				right = middle;
			} else {
				left = middle + 1;
			}
		}
		return right;
	}

	/**
	 * Performs a the search procedure specified by Wigderson.
	 * The oracle tells if the current value is valid.
	 * @param searchWrapper The oracle
	 * @return The smallest valid parameter determined
	 */
	int searchWigderson(SearchWrapper* searchWrapper) const {
		int k = 2;
		int l;
		for (l = 1; l <= 30; l++) {
			if (searchWrapper->step(k)) {
				break;
			}
			k *= 2;
		}
		OGDF_ASSERT(l <= 30);
		return searchBinary(searchWrapper, k / 2, k);
	}
};
}
