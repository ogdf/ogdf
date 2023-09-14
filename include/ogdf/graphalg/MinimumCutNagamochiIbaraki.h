/** \file
 * \brief Calculate minimum cut value for a given Graph.
 *
 * \author Sascha Alder
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

#define OGDF_MINCUTNI_MAXLISTSIZE 100
#define OGDF_MINCUTNI_PRTHR 100
#define OGDF_MINCUTNI_CLUSTERSIZE 10

#include <ogdf/basic/Array.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/MinimumCutModule.h>

#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace ogdf {

/**
 * \brief  Calculate minimum cut value for a given Graph.
 *
 * This class implements a version of Nagamochi-Ibaraki with Padberg-Rinaldi heuristics.
 * Implemented to be faster than the standard mincut algorithm.
 */
class OGDF_EXPORT MinimumCutNagamochiIbaraki : public Logger, public MinimumCutModule<int> {
public:
	/**
	 * \brief Standard constructor of the class
	 * @param p If true use Padberg Rinaldi heuristics
	 * @param preprocessing If true use preprocessing
	 * @param logLevel Level of debug messages
	 */
	MinimumCutNagamochiIbaraki(bool p = true, bool preprocessing = false,
			Logger::Level logLevel = Logger::Level::Default);

	/**
	 * \brief Standard destructor of the class
	 */
	virtual ~MinimumCutNagamochiIbaraki() override;

	/**
	 * \brief Compute a minimum cut value for the given weighted graph
	 * @param G input graph
	 * @param capacity Vector representing the capacity values of edges
	 */
	const int& minCutWeighted(const Graph& G, const std::vector<int>& capacity);

	/**
	 * \brief Compute a minimum cut value for the given unweighted graph
	 * @param G input graph
	 */
	const int& minCutUnweighted(const Graph& G);

	//! Computes a minimum cut on graph \p G.
	virtual inline int call(const Graph& G) override { return minCutUnweighted(G); }

	//! Computes a minimum cut on graph \p G with non-negative \p weights on edges.
	virtual int call(const Graph& G, const EdgeArray<int>& weights) override {
		init(G);
		for (edge e : G.edges) {
			edgeCapacity[m_GC.copy(e)->index()] = weights[e];
		}
		const int& value {minCutWeighted()};
		delete hiddenEdges;
		return value;
	}

	//! Output value of last minimum cut computation
	int value() const override { return barLambda; };

	// @warning Currently this algorithm does not support returning the cut edges.
	virtual const ArrayBuffer<edge>& edges() override { return m_cutEdges; };

	// @warning Currently this algorithm does not support returning a min cut partition.
	virtual const ArrayBuffer<node>& nodes() override { return m_partition; };

	/**
	 * \brief Output the number of Padberg-Rinaldi rounds
	 */
	const unsigned int& getPrRounds() const { return prRounds; };

	/**
	 * \brief Output the number of Nagamochi-Ibaraki rounds (equals to MAO-Computations)
	 */
	const unsigned int& getNIRounds() const { return NIRounds; };

private:
	struct BoundedList {
		explicit BoundedList(ListPure<node> l_init = ListPure<node>(), int nodesInList_init = 0,
				int size_init = 0) {
			list = l_init;
			nodesInList = nodesInList_init;
			size = size_init;
		}

		ListPure<node> list;
		int nodesInList;
		int size;
	};

	struct adjInfo {
		int adj = 0;
		ListIterator<node> placeInL = nullptr;
	};

	struct clusterstruct {
		ListPure<node> clusterNodes;
		node clusterhead;
	};

	//! Initializes member variables.
	void init(const Graph& G);

	/**
	 * \brief Compute a minimum cut value for the given weighted graph, assuming
	 * that #m_GC and #edgeCapacity have already been initialized.
	 */
	const int& minCutWeighted();

	/**
	 * \brief Underlying function for minCut computation
	 */
	void minCut();

	/**
	 * \brief Contracts all Clusters trough contraction of each cluster using contractClusters
	 * @param clusters Vector containing all clusters to contract
	 */
	void contractClusters(const std::vector<clusterstruct>& clusters);

	/**
	 * \brief Contracts one clusters
	 * @param s Clusterhead of the current cluster
	 * @param toContract Nodes of the current cluster
	 * @param clusterlevel Integer representing the cluster in setid
	 * @param clusters Vector containing all clusters to contract
	 */
	void contract(const node& s, const ListPure<node>& toContract, const int& clusterlevel,
			const std::vector<clusterstruct>& clusters);

	/**
	 * \brief Tests rule 1 and 2 of Padberg-Rinaldi of all incident edges of lastContracted and contracts on success
	 * @param lastContracted Node to check
	 */
	void PRPass1_2(const node& lastContracted);

	/**
	 * \brief Test for rule 1 of Padberg-Rinaldi
	 * @param eIndex Index of the edge e to test
	 */
	inline bool PRTest1(const unsigned int& eIndex) { return (edgeCapacity[eIndex] >= barLambda); };

	/**
	 * \brief Test for rule 2 of Padberg-Rinaldi
	 * @param eIndex Index of the edge e to test
	 * @param uIndex Index of node incident to e
	 * @param vIndex Index of node incident to e (uIndex != vIndex)
	 */
	inline bool PRTest2(const unsigned int& eIndex, const unsigned int& uIndex,
			const unsigned int& vIndex) {
		return 2 * edgeCapacity[eIndex] >= min(degree[uIndex], degree[vIndex]);
	};

	/**
	 * \brief Add new node to an existing cluster of head or create new cluster with head as head
	 * @param head Adjacent node with adjacency of currentNode >= barlambda
	 * @param currentNode Node to add to Cluster
	 * @param clusters Needed for adding node to or adding new cluster
	 * @param level Setid for the next new cluster
	 */
	void updateClusters(const node& head, const node& currentNode,
			std::vector<clusterstruct>& clusters, int& level);

	/**
	 * \brief Compute a MAO and contract clusters in it
	 * @param s Starting node
	 */
	node MAOComputation(const node& s);

	/**
	 * \brief Deletes the node given through placeInL from L
	 * @param L List which could contain node (only possible if node was added to list )
	 * @param placeInL Information of node (If null, the node wasn't added to the list)
	 */
	void deleteFromL(BoundedList& L, ListIterator<node>& placeInL);

	/**
	 * \brief Refills L, if it's empty but nodes with the same adjacency exists
	 * @param maxAdj Nodes with adjacency maxAdj are added to L
	 * @param unviewed Remaining nodes which are candidates to add to L
	 * @param L List to add nodes to
	 * @param adjToViewed Adjacency values of nodes
	 */
	void fillL(const int& maxAdj, ListPure<node>& unviewed, BoundedList& L,
			std::vector<adjInfo>& adjToViewed);

	/**
	 * \brief Return first node of L and adjust values of L
	 * @param L Not empty list containing nodes
	 */
	node getFirstNode(BoundedList& L);

	/**
	 * \brief Returns edge corresponding to adj
	 * @param adj Given adjEntry
	 * @param s Given node incident to edge
	 * @param opposite For saving opposite node to s for the edge
	 */
	static edge getAdjEdge(const adjEntry& adj, const node& s, node& opposite) {
		auto e = adj->theEdge();
		opposite = e->opposite(s);
		return e;
	}

	/**
	 * \brief Checks for new upper bound
	 * @param nodeDegree Degree of node to check if lesser than barlambda
	 */
	inline void updateLambda(const int nodeDegree) {
		if (nodeDegree < barLambda && size >= 2) {
			barLambda = nodeDegree;
		}
	}

	bool m_preprocess; //if preprocessing should be used

	bool pr; //if pr should be used

	int size; //G.numberOfNodes()

	unsigned int NIRounds = 0;
	unsigned int prRounds = 0; //successful rounds of PR1_2
	//unsigned int pr3Rounds = 0;
	int barLambda = 0; //mincut value

	GraphCopy m_GC; //Copy of the Graph for changing nodes and edges
	int N; //original size
	int M; //original |E|
	std::vector<int> edgeCapacity; //capacity in graphcopy
	std::vector<int> degree;
	std::vector<int> setid;

	Graph::HiddenEdgeSet* hiddenEdges;

	std::unordered_set<node> allNodes;

	ArrayBuffer<node> m_partition;

	ArrayBuffer<edge> m_cutEdges;
};
}
