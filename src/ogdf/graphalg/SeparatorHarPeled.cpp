/** \file
 * \brief Implementation of class SeparatorHarPeled.
 *
 * \author Thomas Klein
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

#include <ogdf/graphalg/SeparatorHarPeled.h>

namespace ogdf {

void BFSTreeHP::construct() {
	SListPure<node> bfs; // current queue of BFS
	bfs.pushBack(root);

	int currentLevel = 0;

	while (!bfs.empty()) {
		// store all nodes of the current level
		List<node> next;
		for (const node no : bfs) {
			next.pushBack(no);
			levelOfNode[no] = currentLevel;
		}

		bfs.clear();

		// collect all nodes of the next level, record their parents
		for (const node w : next) {
			for (const adjEntry adj : w->adjEntries) {
				node v = adj->twinNode();
				if (!mark[v]) {
					mark[v] = true;
					bfs.pushBack(v);
					parentOfNode[v] = w;
					childrenOfNode[w].pushBack(v);
					edgeToParent[v] = adj->twin();
					inTree[adj->theEdge()] = true;
				}
			}
		}
		currentLevel++;
	}
#ifdef OGDF_DEBUG
	calculateDescendants();
#endif
}

void BFSTreeHP::reconstruct(const Cycle& cycle) {
	node newRoot = cycle.getRoot();
	mark.fill(false); // marks which nodes were visited
	inTree.fill(false); // contains for every edge whether it's in the tree or not
	mark[newRoot] = true;
	int rootLevel = levelOfNode[newRoot]; // storing old level of root

	root = newRoot;

	// 1. keep the subtree of the new root as it is
	List<node> subtree;
	for (node child : childrenOfNode[newRoot]) {
		subtree.pushBack(child);
	}

	while (!subtree.empty()) {
		node child = subtree.popFrontRet();

		levelOfNode[child] = levelOfNode[child] - rootLevel;
		inTree[edgeToParent[child]] = true;
		mark[child] = true;

		for (node grandchild : childrenOfNode[child]) {
			subtree.pushBack(grandchild);
		}
	}

	// 2. start regular BFS at newRoot, ignoring marked nodes (therefore ignoring subtree)
	SListPure<node> bfs; // current queue of BFS
	bfs.pushBack(root);

	int currentLevel = 0;

	while (!bfs.empty()) {
		// store all nodes of the current level
		List<node> next;
		for (const node no : bfs) {
			next.pushBack(no);
			levelOfNode[no] = currentLevel;
		}

		bfs.clear();

		// collect all nodes of the next level, record their parents
		for (const node w : next) {
			for (const adjEntry adj : w->adjEntries) {
				node v = adj->twinNode();
				if (!mark[v]) {
					mark[v] = true;
					bfs.pushBack(v);
					parentOfNode[v] = w;
					childrenOfNode[w].pushBack(v);
					edgeToParent[v] = adj->twin();
					inTree[adj->theEdge()] = true;
				}
			}
		}
		currentLevel++;
	}

	calculateDescendants();
}

void BFSTreeHP::calculateDescendants() {
	// recording the number of descendants in a possibly stupid and inefficient manner
	descendantsOfNode.init(*pGraph, 1);
	for (node n : pGraph->nodes) {
		while (n != root) {
			descendantsOfNode[parentOfNode[n]]++;
			n = parentOfNode[n];
		}
	}
}

void SeparatorHarPeled::reset() {
	psi = nullptr;

	faceLevels.init();
	border.init();
	ringIn.init();
	ringOut.init();
	mainSeparator.init();
	rings.init();

	faceFrontiers.clear();
}

#ifdef OGDF_DEBUG

void SeparatorHarPeled::verifyRing(const Ring& ring) const {
	FaceArray<bool> debugMarker(E, false);
	EdgeArray<bool> debugBorder(*graph, false);
	for (adjEntry adj : ring.entries) {
		debugBorder[adj->theEdge()] = true;
	}
	adjEntry firstAdj = ring.entries.front();
	face startFace = E.rightFace(firstAdj);
	debugMarker[startFace] = true;
	SListPure<face> stack;
	stack.pushBack(startFace);

	int counter = 0;
	while (!stack.empty()) {
		face face1 = stack.popFrontRet();
		counter++;
		for (adjEntry adj : face1->entries) {
			if (!debugBorder[adj->theEdge()]) {
				face neighbour = E.leftFace(adj);
				if (!debugMarker[neighbour]) {
					stack.pushBack(neighbour);
					debugMarker[neighbour] = true;
				}
			}
		}
	}
	OGDF_ASSERT(counter == ring.faces);
}

#endif

bool SeparatorHarPeled::doSeparate(const Graph& G, List<node>& separator, List<node>& first,
		List<node>& second) {
	int delta = ceil(sqrt(graph->numberOfNodes() / 2.0));
	triangulate(*graph);
	E.init(*graph);

	makeTree();

	edge edgeSep = findSeparatorEdge();
	node psiSrc = edgeSep->source();
	node psiTrg = edgeSep->target();
	psi = tree->getLevelOfNode(psiSrc) > tree->getLevelOfNode(psiTrg) ? psiSrc : psiTrg;

	Cycle cycle(tree.get(), edgeSep);

#ifdef OGDF_DEBUG
	OGDF_ASSERT(cycle.getInsideCost() <= 2.0 / 3.0 * graph->numberOfNodes());
#endif

	if (cycle.getSize() <= 2 * delta) {
		return finalize("short_main_separator", cycle.getNodes(), separator, first, second);
	}

	tree->reconstruct(cycle);

	// === storage prep === //
	faceLevels.init(E, -1);
	border.init(*graph, -1);
	ringIn.init(*graph);
	ringOut.init(*graph);
	// === storage prep === //

	findFaceLevels(tree->getRoot());

	buildRings(cycle);

	// finding i0
	int i0 = find_i0(delta);

	// find heavy ring
	int f = E.numberOfFaces();
	int innerRingLevel = i0;
	for (int i = i0; i < rings.size(); i += delta) {
		Ring current = rings[i];
		int size = current.faces;

		if (size < f / 3.0) {
			innerRingLevel = i;
		}

		if (size >= floor(f / 3.0) && size <= ceil(f * 2.0 / 3.0)) {
			return finalize("acceptable_ring", current.nodes, separator, first, second);
		}
	}
	int outerRingLevel = innerRingLevel + delta;
	Ring& innerRing = rings[innerRingLevel];

	// if we got here, innerRing is the last ring that contained less than a third of the faces
	// find border of region R1 / R2
	List<node> region;
	if (findRegions(region, cycle, innerRing, outerRingLevel)) {
		return finalize("region_R", region, separator, first, second);
	}

	// worst case, we actually need to build the region K now
	Ring outer(psi);
	if (outerRingLevel < rings.size()) {
		outer = rings[outerRingLevel];
	}
	if (constructK(region, cycle, innerRing, outer)) {
		return finalize("region_K", region, separator, first, second);
	}

	return false; // if we got here, I messed up or the algorithm doesn't actually work
}

/**
 * Auxiliary method to calculate the weight of a subtree.
 *
 * @param x the root of the subtree
 * @param weights stores the weights of subtrees starting at each node
 * @param marked stores which edges have been visited
 * @return the weight of the node x
 */
int calculateWeight(const node x, NodeArray<int>& weights, EdgeArray<bool>& marked) {
	OGDF_ASSERT(weights[x] == 0);

	int sum = 1;
	for (adjEntry adj : x->adjEntries) {
		if (!marked[adj->theEdge()]) {
			marked[adj->theEdge()] = true;
			sum += calculateWeight(adj->twinNode(), weights, marked);
		}
	}
	weights[x] = sum;
	return sum;
}

/**
 * Auxiliary method to descend into the dual tree and find the first separator edge that
 * divides the tree into two halves that are greater than a third of the nodes.
 *
 * @param x the startnode
 * @param weights stores the weights of each subtree
 * @param marked stores which edges have been visited
 * @param size the total size of the graph
 * @return the edge that separates the dual tree
 */
edge findHeavyEdge(const node x, const NodeArray<int>& weights, EdgeArray<bool>& marked,
		const int size) {
	adjEntry max = nullptr;
	int maxW = 0;
	for (adjEntry adj : x->adjEntries) {
		if (!marked[adj->theEdge()] && weights[adj->twinNode()] >= maxW) {
			marked[adj->theEdge()] = true;
			maxW = weights[adj->twinNode()];
			max = adj;
		}
	}

	OGDF_ASSERT(max != nullptr);

	if (maxW <= size) {
		return max->theEdge();
	}

	return findHeavyEdge(max->twinNode(), weights, marked, size);
}

void SeparatorHarPeled::createDual(Graph& Dual, EdgeArray<edge>& oldEdge) const {
	// 0. compute Dual of the graph - not using ogdf methods because they are too slow
	FaceArray<node> dualNode(E);
	oldEdge.init(Dual, nullptr);

	// create a node for each face
	for (face f : E.faces) {
		node x = Dual.newNode();
		dualNode[f] = x;
	}

	// create an edge for each edge, ignoring tree edges
	for (edge e : graph->edges) {
		if (!tree->isInTree(e)) {
			face left = E.leftFace(e->adjSource());
			face right = E.rightFace(e->adjSource());
			edge newEdge = Dual.newEdge(dualNode[left], dualNode[right]);
			oldEdge[newEdge] = e;
		}
	}

	OGDF_ASSERT(isTree(Dual));
}

edge SeparatorHarPeled::findSeparatorEdge() const {
	Graph Dual;
	EdgeArray<edge> oldEdge;
	createDual(Dual, oldEdge);

	// 1. compute number of nodes in the subtree of each node
	node root = Dual.chooseNode();
	NodeArray<int> subTreeWeight(Dual, 0);
	EdgeArray<bool> marked(Dual, false);
	int total = calculateWeight(root, subTreeWeight, marked);

	OGDF_ASSERT(total == E.numberOfFaces());

	marked.fill(false);

	// 2. climb down in tree until we found the separator edge
	edge heavyEdge = findHeavyEdge(root, subTreeWeight, marked, ceil(2.0 / 3.0 * total));
	return oldEdge[heavyEdge];
}

void SeparatorHarPeled::findFaceLevels(const node root) {
	EdgeArray<bool> marked(*graph, false); // marks which edges have been used for search
	isMultiNode.init(*graph, false); // marks which node has multiple incoming ring-edges

	List<adjEntry> frontier; // the current frontier of the expansion search = adjEntries whose right faces are next
	for (const adjEntry& adj : root->adjEntries) {
		frontier.pushBack(adj);
		marked[adj->theEdge()] = true;
	}

	int level = 1; // current level

	while (!frontier.empty()) {
		List<adjEntry> nextFrontier; // new frontier once we are done with this one
		List<face> nextLining; // the inner lining of the next ring (faces directly on the inside of the ring)

		while (!frontier.empty()) {
			// grab an adj, set the level of the right face
			adjEntry adj = frontier.popFrontRet();
			face rightFace = E.rightFace(adj);

			if (faceLevels[rightFace] == -1) { // is unvisited face

				faceLevels[rightFace] = level;
				nextLining.pushBack(rightFace);

				// update border array
				for (adjEntry neighbourAdj : rightFace->entries) {
					face neighbor = E.leftFace(neighbourAdj);

					// does the adjEntry neighbourAdj lead to a face of level "level-1"?
					if (faceLevels[neighbor] == level - 1) {
						border[neighbourAdj->theEdge()] = level - 1;

						// This is a border edge, and its two adjEntries need to be connected to their nodes
						// there is probably a more elegant way of doing this
						node outNode, inNode;
						adjEntry outAdj, inAdj;
						if (faceLevels[E.rightFace(neighbourAdj)] == level - 1) {
							outNode = neighbourAdj->theNode();
							outAdj = neighbourAdj;
							inNode = neighbourAdj->twinNode();
							inAdj = neighbourAdj->twin();
						} else {
							outNode = neighbourAdj->twinNode();
							outAdj = neighbourAdj->twin();
							inNode = neighbourAdj->theNode();
							inAdj = neighbourAdj;
						}

						if (!ringOut[outNode].empty()) {
							isMultiNode[outNode] = true;
						}
						ringOut[outNode].pushBack(outAdj);

						if (!ringIn[inNode].empty()) {
							isMultiNode[inNode] = true;
						}
						ringIn[inNode].pushBack(inAdj);
					}
				}
			}

			// collect adjEntries for the next level of the search
			for (adjEntry out : adj->twinNode()->adjEntries) {
				if (!marked[out->theEdge()]) {
					nextFrontier.pushBack(out);
					marked[out->theEdge()] = true;
				}
			}
		}

		frontier = nextFrontier;
		level++;

		faceFrontiers.pushBack(nextLining);
	}

#ifdef OGDF_DEBUG
	for (face f : E.faces) {
		OGDF_ASSERT(faceLevels[f] > -1);
	}
#endif
}

/**
 * Auxiliary struct to hold all the information we need during ring creation.
 * This is necessary to maintain data about the inside of each ring, even if rings degenerate.
 */
struct ExpansionData {
	FaceArray<bool> counted; // holds for each face whether it has been counted for cycle size calculation
	EdgeArray<int> ringBorder; // like border, but only with those edges that actually form the rings
	SListPure<face> searchOrigin; // starting points for the next round of BFS search over faces
	int numberOfFaces; // total number of faces counted so far
	int borderValue; // which value in ringBorder is the current border

	ExpansionData(const ConstCombinatorialEmbedding& E, const Graph& G) {
		counted.init(E, false);
		ringBorder.init(G, 0);
		numberOfFaces = 0;
		borderValue = 1;
	}
};

/**
 * Performs one step of the expansion: For all faces of the current level, we expand until we
 * run into the border defined by the current ring everywhere.
 *
 * @param E the embedding of the graph
 * @param faces the inner lining of the current region border
 * @param data auxiliary struct that holds data
 * @return the (cumulative) number of faces inside the current ring.
 */
int expandSearch(const ConstCombinatorialEmbedding& E, const List<face>& faces, ExpansionData& data) {
	// push all faces into searchOrigin
	for (face f : faces) {
		if (!data.counted[f]) {
			data.searchOrigin.pushBack(f);
			data.counted[f] = true;
		}
	}

	while (!data.searchOrigin.empty()) {
		face f = data.searchOrigin.popFrontRet();
		++data.numberOfFaces;

		for (adjEntry adj : f->entries) {
			if (data.ringBorder[adj->theEdge()] < data.borderValue) {
				face leftFace = E.leftFace(adj);
				if (!data.counted[leftFace]) {
					data.searchOrigin.pushBack(leftFace);
					data.counted[leftFace] = true;
				}
			}
		}
	}
	return data.numberOfFaces;
}

void SeparatorHarPeled::buildRings(const Cycle& cycle) {
	// the two legs of the main separator: nodesA leads to psi, nodesB leads back
	List<node> nodesA;
	List<node> nodesB;
	List<adjEntry> entryPath; // all adjEntries along separator to psi
	List<adjEntry> entryPathBackwards; // all adjEntries back to root

	// 1. find path to psi
	bool foundRoot = false;
	auto edgeIt = cycle.getEdges().cbegin();
	for (const node& x : cycle.getNodes()) {
		if (x == cycle.getRoot()) {
			foundRoot = true;
			entryPath.pushFront((*edgeIt)->twin());
			++edgeIt;
			continue;
		}

		if (foundRoot) {
			nodesB.pushBack(x);
			entryPathBackwards.pushFront((*edgeIt)->twin());
			++edgeIt;
		} else {
			nodesA.pushFront(x);
			entryPath.pushFront((*edgeIt)->twin());
			++edgeIt;
		}
	}

	entryPath.conc(entryPathBackwards);

	// mS always has the form of: root->along longer leg to psi -> backwards along shorter leg to root
	// (regardless whether it is clockwise or not)
	mainSeparator.init(*graph, nullptr);
	for (adjEntry entry : entryPath) {
		mainSeparator[entry->theNode()] = entry;
	}
	int numberOfRings = min(nodesA.size(), nodesB.size());
	if (nodesA.size() == nodesB.size()) {
		numberOfRings--;
	}
	rings.init(numberOfRings);

	// last node in nodesA is psi, the deeper one of u,v (by construction of cycle)

	// 2. creating a data storage unit for the expansion
	ExpansionData data(E, *graph);

	// 3. walk over path both forwards and backwards, creating only as many rings as can be completed
	auto adjIt = entryPath.cbegin();
	adjIt++; // adjIt starts at root, but we don't construct a ring from root
	auto liningIt = faceFrontiers.cbegin();
	auto it = nodesA.cbegin(), revIt = nodesB.cbegin();

#ifdef OGDF_DEBUG
	NodeArray<bool> assertArray(*graph, false);
#endif
	for (int idx = 0; idx < numberOfRings && it != nodesA.cend() && revIt != nodesB.cend();
			++it, ++revIt, ++adjIt, ++liningIt, ++idx) {
		// build ring
		Ring ring(*it, *revIt, *adjIt, *this);

		// build fence around the area that this ring covers
		for (adjEntry entry : ring.entries) {
			data.ringBorder[entry->theEdge()] = data.borderValue;
		}

		ring.faces = expandSearch(E, *liningIt, data); // counting faces on the inside

		data.borderValue++;
		rings[idx] = ring;

#ifdef OGDF_DEBUG
		verifyRing(ring);
		for (node n : ring.nodes) {
			OGDF_ASSERT(!assertArray[n]);
			assertArray[n] = true;
		}
#endif
	}
}

int SeparatorHarPeled::find_i0(int delta) const {
	for (int i0 = 0; i0 < rings.size(); i0++) {
		int sum = 0;
		for (int i = i0; i < rings.size(); i += delta) {
			sum += rings[i].getSize();
		}
		if (sum <= graph->numberOfNodes() / delta) {
			return i0;
		}
	}
	OGDF_ASSERT(false); // we should never get here
	return 0;
}

bool SeparatorHarPeled::findRegions(List<node>& region, const Cycle& cycle, const Ring& inner,
		int outerIdx) const {
	// main complication: is the outer ring an actual simple cycle, or just one node?
	Ring outer(psi);
	if (outerIdx < rings.size()) {
		outer = rings[outerIdx];
	}

	// try to find R1 first, if that doesn't work, try to find R2
	return (findRegion(region, cycle, inner, outer, true)
			|| findRegion(region, cycle, inner, outer, false));
}

void SeparatorHarPeled::walkAlongSeparator(node startNode, node endNode,
		EdgeArray<bool>& regionBorder, List<node>& region) const {
	node n = startNode;
	while (n != endNode) {
		adjEntry nextAdj = mainSeparator[n];
		region.pushBack(n);
		regionBorder[nextAdj->theEdge()] = true;
		n = nextAdj->twinNode();
	}
}

void SeparatorHarPeled::walkAlongRing(const Ring& ring, bool firstSection, bool invert,
		EdgeArray<bool>& regionBorder, List<node>& region) const {
	// Assumption: Ring walks clockwise around root of tree

	// 1. get all adjEntries from the inside of the ring
	List<adjEntry> entries = ring.getSectionInSeparator(firstSection);

	// 2. if we want inverse direction, invert all adjEntries
	if (invert) {
		List<adjEntry> tmp;
		for (adjEntry adj : entries) {
			tmp.pushFront(adj->twin());
		}
		entries = tmp;
	}

	// walk along the ring, store the node of each entry and mark its edge
	for (adjEntry adj : entries) {
		regionBorder[adj->theEdge()] = true;
		region.pushBack(adj->theNode());
	}
}

bool SeparatorHarPeled::findRegion(List<node>& region, const Cycle& cycle, const Ring& inner,
		const Ring& outer, bool inside) const {
	region.clear(); // important to do this initially, because of the way we try both regions sequentially

	EdgeArray<bool> regionBorder(*graph, false);

	// walk from inner in to outer in (along main separator)
	walkAlongSeparator(inner.in, outer.in, regionBorder, region);

	// walk along outer ring from outer.in to outer.out - if outer ring is degenerate, nothing happens
	walkAlongRing(outer, (cycle.getClockwise() != inside), !(cycle.getClockwise() != inside),
			regionBorder, region);

	// walk from outer out to inner out (along main separator)
	walkAlongSeparator(outer.out, inner.out, regionBorder, region);

	// walk along inner ring from inner.out to inner.in
	walkAlongRing(inner, (cycle.getClockwise() != inside), (cycle.getClockwise() != inside),
			regionBorder, region);

	return testRegionSize(inner.in, regionBorder, (cycle.getClockwise() != inside), region.size());
}

bool SeparatorHarPeled::testRegionSize(node startNode, const EdgeArray<bool>& regionBorder,
		bool inside, int regionSize) const {
	// get one of the faces inside of the region R1
	adjEntry startAdj = mainSeparator[startNode];
	face startFace = inside ? E.rightFace(startAdj) : E.leftFace(startAdj);

	// start BFS at startFace
	FaceArray<bool> marked(E, false);
	marked[startFace] = true;
	int counter = 0;
	SListPure<face> stack;
	stack.pushBack(startFace);

	while (!stack.empty()) {
		face f = stack.popFrontRet();
		counter++;

		for (adjEntry adj : f->entries) {
			if (!regionBorder[adj->theEdge()]) {
				face neigh = E.leftFace(adj);
				if (!marked[neigh]) {
					stack.pushBack(neigh);
					marked[neigh] = true;
				}
			}
		}
	}
	// not doing exactly what HP proposed, but this works better
	int nodesInsideRegion = (counter - regionSize) / 2 + 1;
	return nodesInsideRegion + regionSize >= ceil(graph->numberOfNodes() / 3.0);
}

bool SeparatorHarPeled::constructK(List<node>& region, const Cycle& cycle, const Ring& inner,
		const Ring& outer) const {
	region.clear(); // important to do this initially, because of the way we try both regions sequentially

	EdgeArray<bool> regionBorder(*graph, false);

	// walk from inner in to outer in (along main separator)
	walkAlongSeparator(inner.in, outer.in, regionBorder, region);

	// walk along the outer ring
	walkAlongRing(outer, !cycle.getClockwise(), cycle.getClockwise(), regionBorder, region);

	// walk from outer out to inner out (along main separator)
	walkAlongSeparator(outer.out, inner.out, regionBorder, region);

	// walk along the inner ring
	walkAlongRing(inner, !cycle.getClockwise(), !cycle.getClockwise(), regionBorder, region);

	return true;
}

void SeparatorHarPeled::makeTree() {
	tree = std::make_shared<BFSTreeHP>(*graph, getStartNode(*graph));
}

bool SeparatorHarPeled::finalize(std::string exit, const List<node>& region, List<node>& separator,
		List<node>& first, List<node>& second) {
	exitPoint = exit;
	for (node n : region) {
		separator.pushBack(graph->original(n));
		graph->delNode(n);
	}
	return separateComponents(*graph, separator, first, second, true);
}

}
