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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/graphalg/NodeColoringModule.h>

#include <algorithm>
#include <cmath>

namespace ogdf {

using NColor = NodeColoringModule::NodeColor;

bool NodeColoringModule::checkColoring(const Graph& graph, const NodeArray<NColor>& colors) const {
	OGDF_ASSERT(colors.graphOf() == &graph);
	for (edge e : graph.edges) {
		if ((colors[e->source()] == colors[e->target()]) && !e->isSelfLoop()) {
			return false;
		}
	}
	return true;
}

void NodeColoringModule::createBuckets(const Graph& graph, int size,
		Array<Array<node>>& buckets) const {
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

int NodeColoringModule::getNeighborDegrees(const node& v) const {
	int result = 0;
	for (adjEntry adj : v->adjEntries) {
		node twin = adj->twinNode();
		result += twin->degree();
	}
	return result;
}

int NodeColoringModule::getMaximumDegreeNode(const Graph& graph, node& maxDegreeNode) const {
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

int NodeColoringModule::getMaximumDegreeNodes(const Graph& graph, List<node>& maxDegreeNodes) const {
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

int NodeColoringModule::getMinimumDegreeNode(const Graph& graph, node& minDegreeNode) const {
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

int NodeColoringModule::getMinimumDegreeNodes(const Graph& graph, List<node>& minDegreeNodes) const {
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

NColor NodeColoringModule::getMaximumNodeColor(NodeArray<NColor>& colors) {
	auto maxColor = NColor(0);
	for (node v : colors.graphOf()->nodes) {
		if (colors[v] > maxColor) {
			maxColor = colors[v];
		}
	}
	return maxColor;
}

void NodeColoringModule::reverseNodeTable(const Graph& graphOrig, const Graph& graphNew,
		const NodeArray<node>& orig2New, NodeArray<node>& new2Orig) const {
	OGDF_ASSERT(&graphOrig == orig2New.graphOf());
	new2Orig = NodeArray<node>(graphNew);
	for (node v : graphOrig.nodes) {
		if (orig2New[v]) {
			new2Orig[orig2New[v]] = v;
		}
	}
}

void NodeColoringModule::ramseyAlgorithm(const Graph& graph, List<node>& clique,
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

void NodeColoringModule::cliqueRemoval(const Graph& graph, List<node>& independentSet) const {
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

int NodeColoringModule::searchLinear(SearchWrapper* searchWrapper, int start, int end) const {
	OGDF_ASSERT(start <= end);
	for (int k = start; k <= end; k++) {
		if (searchWrapper->step(k)) {
			return k;
		}
	}
	return end;
}

int NodeColoringModule::searchBinary(SearchWrapper* searchWrapper, int start, int end) const {
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

int NodeColoringModule::searchWigderson(SearchWrapper* searchWrapper) const {
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

}
