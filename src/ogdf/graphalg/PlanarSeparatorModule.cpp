/** \file
 * \brief Implementation of class PlanarSeparatorModule and
 * 		various auxiliary data structures for planar separator algorithms.
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

#include <ogdf/graphalg/Matching.h>
#include <ogdf/graphalg/PlanarSeparatorModule.h>

namespace ogdf {

BFSTreeClassical::BFSTreeClassical(GraphCopy& G, node rootNode, unsigned int iterations,
		bool findLevelsSimple)
	: ArrayBFSTree(G, rootNode), heightMaxIterations(iterations), simple {findLevelsSimple} {
	m_ratio = 2.0 * sqrt(2.0);

	construct(rootNode, heightMaxIterations);
}

int BFSTreeClassical::getSizeOfLevel(int level) const {
	if (level < 0 || level >= levels.size()) {
		return 0;
	}
	return (*levels.get(level)).size();
}

List<node> BFSTreeClassical::getLevel(int level) const {
	if (level < 0 || level >= levels.size()) {
		List<node> empty;
		return empty;
	}
	return *levels.get(level);
}

List<node> BFSTreeClassical::getNodesFromTo(int start, int end) const {
	List<node> list;
	for (int i = start; i < end; i++) {
		for (const node no : *levels.get(i)) {
			list.pushBack(no);
		}
	}
	return list;
}

List<node> BFSTreeClassical::getNodesFrom(int start) const {
	return getNodesFromTo(start, levels.size());
}

void BFSTreeClassical::construct(node rootNode, unsigned int numIterations) {
	for (unsigned int iteration = 0; iteration < numIterations; iteration++) {
		// stores the size and the index of the smallest level
		int smallestSize = pGraph->numberOfNodes();
		int idxOfSmallestLevel = -1;

		root = rootNode;
		init(); // fresh setup to support multiple iterations - not all necessary, but easier this way

		// clear everything
		levels.clear();

		SListPure<node> bfs; // current queue of BFS
		bfs.pushBack(root);

		// main BFS
		currentLevel = 0;
		belowMiddle = true;

		int sum = 0; // stores number of nodes (maybe later node costs -> template) of all nodes so far

		while (!bfs.empty()) {
			// store all nodes of the current level
			List<node> next;
			for (const node no : bfs) {
				sum += 1;
				next.pushBack(no);
				levelOfNode[no] = currentLevel;
			}
			levels.pushBack(next);

			// check if the current level happens to be a suitable separator
			if (next.size() < m_ratio * sqrt(pGraph->numberOfNodes())
					&& sum > 1 / 3.0 * pGraph->numberOfNodes()
					&& sum < 2 / 3.0 * pGraph->numberOfNodes()) {
				// check if we are beating the current best known level
				if (next.size() < smallestSize) {
					smallestSize = next.size();
					idxOfSmallestLevel = currentLevel;
				}
			}

			// step 4: check if all the nodes so far were more than half of total
			if (belowMiddle && sum > pGraph->numberOfNodes() / 2.0) {
				t1 = currentLevel;
				k = sum;
				belowMiddle = false;
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
						edgeToParent[v] = adj->twin();
						inTree[adj->theEdge()] = true;
					}
				}
			}
			currentLevel++;
		}

		// only return for last iteration
		if (iteration == numIterations - 1) {
			if (idxOfSmallestLevel > 0) {
				currentLevel = idxOfSmallestLevel;
				return; // we found a suitable level
			}
		}

		// step 5: find the two separator levels
		if (simple) {
			findLevelsSimple();
		} else {
			findLevels();
		}

		// record for each node between levels t0 and t2 the number of descendants
		for (int i = t2 - 1; i > t0; i--) {
			for (const node x : getLevel(i)) {
				descendantsOfNode[parentOfNode[x]] += descendantsOfNode[x];
				childrenOfNode[parentOfNode[x]].pushBack(x);
			}
		}

		// tree is constructed, now select new root from highest level and start fresh
		rootNode = getLevel(currentLevel - 1).chooseElement();
	}

	currentLevel = -1; // finalize
}

void BFSTreeClassical::visit(node v, node parent, adjEntry adj, SListPure<node>& bfs) {
	if (!mark[v]) {
		mark[v] = true;
		bfs.pushBack(v);
		parentOfNode[v] = parent;
		edgeToParent[v] = adj->twin();
		inTree[adj->theEdge()] = true;
	}
}

void BFSTreeClassical::reconstruct() {
	levels.clear();

	SListPure<node> bfs; // current queue of BFS
	bfs.pushBack(root);

	currentLevel = 0;

	while (!bfs.empty()) {
		// store all nodes of the current level
		List<node> next;
		for (const node no : bfs) {
			next.pushBack(no);
			levelOfNode[no] = currentLevel;
		}
		levels.pushBack(next);

		bfs.clear();

		// collect all nodes of the next level, record their parents
		for (const node w : next) {
			std::set<node> neighbors; // storing neighbors of w for quick access
			List<adjEntry> oldAdjEntries;
			for (const adjEntry adj : w->adjEntries) {
				neighbors.insert(adj->twinNode());
				oldAdjEntries.pushBack(adj);
			}

			for (const adjEntry adj : oldAdjEntries) {
				node v = adj->twinNode();
				visit(v, w, adj, bfs);

				adjEntry adjBefore = adj;
				adjEntry nextInFace = adj->faceCycleSucc();

				while (nextInFace->twinNode() != w) {
					if (neighbors.find(nextInFace->twinNode()) == neighbors.end()) {
						// w is not connected to nextInFace->twinNode, but they are in the same face, so make connection
						nextInFace = nextInFace->faceCycleSucc();
						edge newEdge = pGraph->newEdge(adjBefore, nextInFace, Direction::after);
						neighbors.insert(nextInFace->theNode());
						adjBefore = newEdge->adjSource();
						visit(nextInFace->theNode(), w, newEdge->adjSource(), bfs);

					} else {
						break; // next one would end loop anyway
					}
				}
			}
		}
		currentLevel++;
	}

	// fill up descendant information from the bottom
	for (const auto& lvl : reverse(levels)) {
		for (const node x : lvl) {
			if (x != root) {
				descendantsOfNode[parentOfNode[x]] += descendantsOfNode[x];
				childrenOfNode[parentOfNode[x]].pushBack(x);
			}
		}
	}
}

void BFSTreeClassical::findLevels() {
	t0 = t1;
	while (getSizeOfLevel(t0) + 2 * (t1 - t0) > 2 * sqrt(k)) {
		t0--;
	}

	t2 = t1 + 1;
	while (getSizeOfLevel(t2) + 2 * (t2 - t1 - 1) > 2 * sqrt(pGraph->nodes.size() - k)) {
		t2++;
	}
}

void BFSTreeClassical::findLevelsSimple() {
	t0 = t1;
	while (getSizeOfLevel(t0) > sqrt(pGraph->numberOfNodes())) {
		t0--;
	}

	t2 = t1 + 1;
	while (getSizeOfLevel(t2) > sqrt(pGraph->numberOfNodes())) {
		t2++;
	}
}

void BFSTreeClassical::restructure(List<node>& separator, List<node>& second, bool useTriBFS) {
	removeSeparatorLevels(separator, second);
	createNewRoot(useTriBFS);
}

// removes the nodes below t0 and above t2 and stores them in second
// removes the nodes in t0 and t2 and stores them in the separator
void BFSTreeClassical::removeSeparatorLevels(List<node>& separator, List<node>& second) {
	// 1. remove all nodes from 0 to t0
	for (const node removeNode : getNodesFromTo(0, t0)) {
		second.pushBack(pGraph->original(removeNode));
		pGraph->delNode(removeNode);
	}

	// 2. remove levels t2 and above
	for (const node removeNode : getNodesFrom(t2 + 1)) {
		second.pushBack(pGraph->original(removeNode));
		pGraph->delNode(removeNode);
	}

	// 3. put all nodes from level t0 in the separator
	for (const node sepNode : getLevel(t0)) {
		separator.pushBack(pGraph->original(sepNode));
		pGraph->delNode(sepNode);
	}

	// 4. put all nodes from level t2 in the separator
	for (const node sepNode : getLevel(t2)) {
		separator.pushBack(pGraph->original(sepNode));
		pGraph->delNode(sepNode);
	}
}

void BFSTreeClassical::createNewRoot(bool useTriBFS) {
	root = pGraph->newNode();
	parentOfNode[root] = nullptr;
	descendantsOfNode[root] = 1;
	edgeToParent[root] = nullptr;
	levelOfNode[root] = t0;

	for (const node end : getLevel(t0 + 1)) {
		edge e = pGraph->newEdge(end, root);
		parentOfNode[end] = root;
		edgeToParent[end] = e->adjSource();
		inTree[e] = true;
		descendantsOfNode[root] += descendantsOfNode[end];
		childrenOfNode[root].pushBack(end);
	}

	OGDF_ASSERT(isPlanar(*pGraph));
	planarEmbedPlanarGraph(*pGraph); // TODO there should be a way of avoiding this

	if (useTriBFS) {
		init();
		reconstruct();
		// if any of these asserts fail, I messed up tri-BFS
		OGDF_ASSERT(isPlanar(*pGraph));
		OGDF_ASSERT(isSimple(*pGraph));
		OGDF_ASSERT(pGraph->representsCombEmbedding());
		OGDF_ASSERT(pGraph->numberOfEdges() == 3 * pGraph->numberOfNodes() - 6);
	}
}

bool NodeExpulsor::apply(const Graph& G, List<node>& separator, List<node>& first,
		List<node>& second) {
	// we need constant time access to whether a node is in one of the partitions
	NodeArray<short> assignments;
	assignments.init(G);

	for (node n : separator) {
		assignments[n] = 0;
	}
	for (node n : first) {
		assignments[n] = 1;
	}
	for (node n : second) {
		assignments[n] = 2;
	}

	double ratio = 2.0 / 3.0 * G.nodes.size();

	List<node> removeList;
	for (node n : separator) {
		bool hasNeighbourInFirst = false;
		bool hasNeighbourInSecond = false;
		for (adjEntry entry : n->adjEntries) {
			if (assignments[entry->twinNode()] == 1) {
				hasNeighbourInFirst = true;
			}
			if (assignments[entry->twinNode()] == 2) {
				hasNeighbourInSecond = true;
			}
		}
		// if a node has neighbours in neither of the two, assign it to the smaller component
		if (!(hasNeighbourInFirst || hasNeighbourInSecond)) {
			removeList.pushBack(n);
			if (first.size() < second.size()) {
				if (!keepBalance || first.size() + 1 < ratio) {
					first.pushBack(n);
					assignments[n] = 1;
				}
			} else {
				if (!keepBalance || second.size() + 1 < ratio) {
					second.pushBack(n);
					assignments[n] = 2;
				}
			}
			continue;
		}

		// if the node has no neighbours in first, put it into second
		if (!hasNeighbourInFirst && (!keepBalance || second.size() + 1 < ratio)) {
			removeList.pushBack(n);
			second.pushBack(n);
			assignments[n] = 2;
			continue;
		}

		if (!hasNeighbourInSecond && (!keepBalance || first.size() + 1 < ratio)) {
			removeList.pushBack(n);
			first.pushBack(n);
			assignments[n] = 1;
		}
	}

	for (node n : removeList) {
		separator.removeFirst(n);
	}

	return !removeList.empty(); // returns whether we did something or not
}

void DMDecomposer::setupAssignments(const Graph& G, const List<node>& separator,
		const List<node>& first, const List<node>& second) {
	const List<node>& shortList = first.size() > second.size() ? second : first;
	const List<node>& longList = first.size() > second.size() ? first : second;

	assignments.init(G, -1);

	for (const node s : separator) {
		assignments[s] = 0;
	}
	for (const node s : shortList) {
		assignments[s] = 1;
	}
	for (const node s : longList) {
		assignments[s] = 2;
	}
}

void DMDecomposer::reset() {
	assignments.init();

	bipartSep.clear();
	bipartB.clear();

	SI.clear();
	SX.clear();
	SR.clear();

	BI.clear();
	BX.clear();
	BR.clear();
}

bool DMDecomposer::apply(const Graph& G, List<node>& separator, List<node>& first,
		List<node>& second) {
	reset();

	int sizeBefore = separator.size();

	setupAssignments(G, separator, first, second);

	Graph graph; // avoiding GraphCopy because that would be overkill and expensive

	isInS.init(graph, false);
	clone.init(G, nullptr);
	unclone.init(graph, nullptr);

	// 1. build bipartite graph
	bipartiteGraph(graph, separator);

	// 3. find Dulmage Mendelsohn Decompositon
	decompose(graph, bipartSep, SI, BX);
	decompose(graph, bipartB, SX, BI);

	// 4. calculate remainder lists SR and BR
	calculateRemainders(graph);

	// Now that we have the full decomposition, two possible separators emerge:
	// 1. SX u SR u BX   => Sep = Sep - SI + BX,            A = A + SI,       B = B - BX
	// 2. SX u BR u BX   => Sep = Sep - SI - SR + BX + BR,  A = A + SI + SR,  B = B - BX - BR

	// figure out which one has the better ratio
	chooseBalancedDecomposition(G, separator, first, second);

	return sizeBefore - separator.size() > 0;
}

void DMDecomposer::bipartiteGraph(Graph& graph, const List<node>& separator) {
	for (node s : separator) {
		node c = graph.newNode();
		bipartSep.pushBack(c);
		isInS[c] = true;
		clone[s] = c;
		unclone[c] = s;

		for (adjEntry adj : s->adjEntries) {
			node v = adj->twinNode();
			if (assignments[v] == 2) { // i.e. if v is in the larger list
				node w;
				if (clone[v] == nullptr) {
					w = graph.newNode();
					bipartB.pushBack(w);
					clone[v] = w;
					unclone[w] = v;
				} else {
					w = clone[v];
				}
				graph.newEdge(c, w);
			}
		}
	}

	// 2. solve matching problem in bipartite graph
	flow.init(graph, false);
	Matching::findMaximumCardinalityMatching(graph, bipartB, bipartSep, flow);
}

void DMDecomposer::chooseBalancedDecomposition(const Graph& G, List<node>& separator,
		List<node>& first, List<node>& second) {
	// shorthands for length of shorter / longer list
	int shorter = min(first.size(), second.size());
	int longer = max(first.size(), second.size());

	int lengthA = shorter + SI.size();
	int lengthB = longer - BX.size();

	double ratioSimple = min(lengthA, lengthB) / (double)max(lengthA, lengthB);

	lengthA += SR.size();
	lengthB -= BR.size();

	double ratioComplex = min(lengthA, lengthB) / (double)max(lengthA, lengthB);

	OGDF_ASSERT(SR.size() == BR.size()); // either I'm stupid or this should come out to be the same...

	// we definitely move SI and BX so let's do that
	for (node n : SI) {
		assignments[unclone[n]] = 1;
	}
	for (node n : BX) {
		assignments[unclone[n]] = 0;
	}

	// if complex ratio was better, we also need to move the other
	if (ratioSimple < ratioComplex) {
		for (node n : SR) {
			assignments[unclone[n]] = 1;
		}
		for (node n : BR) {
			assignments[unclone[n]] = 0;
		}
	}

	// finally, translate assignments-array back into lists
	translateAssignments(G, separator, first, second);
}

void assign(NodeArray<short>& loc, const List<node>& list, short val) {
	for (node n : list) {
		loc[n] = val;
	}
}

void DMDecomposer::calculateRemainders(const Graph& graph) {
	// next, calculate SR and BR
	NodeArray<short> location(graph,
			1); // dumb name, but assignments was taken - default value means is in BR

	assign(location, bipartSep, -1);
	assign(location, SI, -2);
	assign(location, SX, -3);

	assign(location, BI, 2);
	assign(location, BX, 3);

	for (node n : graph.nodes) {
		if (location[n] == -1) {
			SR.pushBack(n);
		}
	}

	for (node n : graph.nodes) {
		if (location[n] == 1) {
			BR.pushBack(n);
		}
	}
}

void DMDecomposer::translateAssignments(const Graph& G, List<node>& separator, List<node>& first,
		List<node>& second) const {
	separator.clear();
	first.clear();
	second.clear();

	for (node n : G.nodes) {
		if (assignments[n] == 0) {
			separator.pushBack(n);
		} else if (assignments[n] == 1) {
			first.pushBack(n);
		} else if (assignments[n] == 2) {
			second.pushBack(n);
		} else {
			OGDF_ASSERT(false); // every node should have one of the tree possible assignments
		}
	}
}

void DMDecomposer::decompose(const Graph& graph, const List<node>& bipart, List<node>& reachableSep,
		List<node>& reachableB) {
	List<node> bfsStartNodes;
	for (node n : bipart) {
		bool isStartNode = true;
		for (adjEntry adj : n->adjEntries) {
			edge e = adj->theEdge();
			if (flow[e]) {
				isStartNode = false;
				break;
			}
		}
		if (isStartNode) {
			bfsStartNodes.pushBack(n);
		}
	}

	alternatingBFS(graph, bfsStartNodes, reachableSep, reachableB);
}

// starts an alternating BFS at startNodes
bool DMDecomposer::alternatingBFS(const Graph& G, const List<node>& startNodes,
		List<node>& reachableSep, List<node>& reachableB) {
	bool flowVal = false;
	NodeArray<bool> marked(G, false);

	SListPure<List<node>> stack;
	stack.pushBack(startNodes);

	while (!stack.empty()) {
		List<node> nextLevel;

		for (node n : stack.popFrontRet()) {
			for (adjEntry adj : n->adjEntries) {
				if (flow[adj->theEdge()] == flowVal) {
					node next = adj->twinNode();
					if (!marked[next]) {
						nextLevel.pushBack(next);
						marked[next] = true;
					}
				}
			}
			if (isInS[n]) {
				reachableSep.pushBack(n);
			} else {
				reachableB.pushBack(n);
			}
		}
		if (!nextLevel.empty()) {
			stack.pushBack(nextLevel);
		}

		flowVal = !flowVal;
	}

	return true;
}

bool PlanarSeparatorModule::separateComponents(GraphCopy& G, List<node>& separator,
		List<node>& first, List<node>& second, bool skip) const {
	int n = G.numberOfNodes();

	OGDF_ASSERT(n != 0); // should be checked in setup

	if (n == 1) { // special case
		separator.pushBack(G.original(G.firstNode()));
		return true;
	}

	if (n == 2) { // special case
		if (isConnected(G)) {
			for (node no : G.nodes) {
				separator.pushBack(G.original(no));
			}
		} else {
			first.pushBack(G.original(G.nodes.head()));
			second.pushBack(G.original(G.nodes.tail()));
		}
		return true;
	}

	// if we are using this to assign nodes to lists after the algorithm has run, skip this distinction
	// - for very small graphs, we might terminate with a connected remainder < 2/3 n
	if (!skip) {
		if (isConnected(G)) {
			return false;
		} // in this case, just run the normal algorithm
	}
	// otherwise, prepare to manage multiple components

	NodeArray<int> comps;
	comps.init(G);
	std::map<int, int> componentSizes; // maps index to number of nodes of component

	connectedComponents(G, comps, componentSizes);

	// we need the map componentSizes sorted by value, so copy to vector and sort the vector
	std::vector<std::pair<int, int>> vec;
	for (const auto& item : componentSizes) {
		vec.emplace_back(item);
	}
	auto cmp = [](const std::pair<int, int>& x, const std::pair<int, int>& y) {
		return x.second > y.second;
	};
	std::sort(vec.begin(), vec.end(), cmp);

	auto compIt = vec.begin();

	// if we are using this for the dual-algorithm, we have to skip this part because the larger chunk
	// might actually be too large, which would not matter though
	if (!skip) {
		// 1. if there is a component > 2/3, drop everything else into the second list and solve main component

		if ((*compIt).second > 2.0 / 3.0 * n) {
			int biggestIndex = (*compIt).first;
			SListPure<node> markedForDeletion;
			for (const node& no : G.nodes) {
				if (comps[no] != biggestIndex) {
					second.pushBack(G.original(no));
					markedForDeletion.pushBack(no);
				}
			}
			for (const node& no : markedForDeletion) {
				G.delNode(no);
			}
			return false;
		}
	}

	// 2. otherwise, just add each component to the currently smaller list
	int fLength = first.size();
	int sLength = second.size();
	std::map<int, bool> targetList; // holds for each component whether to put it in first or second list

	for (; compIt != vec.end(); ++compIt) {
		int idx = (*compIt).first;
		int size = (*compIt).second;
		targetList[idx] = fLength < sLength;

		if (fLength < sLength) {
			fLength += size;
		} else {
			sLength += size;
		}
	}

	for (const node& no : G.nodes) {
		List<node>& list = targetList[comps[no]] ? first : second;
		list.pushBack(G.original(no));
	}
	return true;
}

int PlanarSeparatorModule::connectedComponents(const Graph& G, NodeArray<int>& component,
		std::map<int, int>& compSizes) const {
	int nComponent = 0;
	component.fill(-1);

	ArrayBuffer<node> S;

	for (node v : G.nodes) {
		if (component[v] != -1) {
			continue;
		}

		S.push(v);
		component[v] = nComponent;
		compSizes[nComponent] = 1;

		while (!S.empty()) {
			node w = S.popRet();
			for (adjEntry adj : w->adjEntries) {
				node x = adj->twinNode();
				if (component[x] == -1) {
					component[x] = nComponent;
					compSizes[nComponent]++;
					S.push(x);
				}
			}
		}

		++nComponent;
	}

	return nComponent;
}

// private constructor, for internal use only (while expanding cycle)
Cycle::Cycle(BFSTree* bfsTree, List<node>& nodeList, List<adjEntry>& edgeList, node root,
		bool clockwise)
	: tree {bfsTree}, nodes {nodeList}, edges {edgeList}, cycleRoot {root}, isClockwise {clockwise} {
	isOnCycle.init(*tree->getGraph(), false);
	isEdgeOnCycle.init(*tree->getGraph(), false);

	for (node no : nodes) {
		isOnCycle[no] = true;
	}
	for (adjEntry adj : edges) {
		isEdgeOnCycle[adj->theEdge()] = true;
	}
}

Cycle::Cycle(BFSTree* bfsTree, bool clockwise) : tree {bfsTree}, isClockwise {clockwise} {
	isOnCycle.init(*tree->getGraph(), false);
	isEdgeOnCycle.init(*tree->getGraph(), false);
}

void Cycle::init(List<node>& nodeList, List<adjEntry>& edgeList, node root) {
	nodes = nodeList;
	edges = edgeList;
	cycleRoot = root;

	for (node no : nodes) {
		isOnCycle[no] = true;
	}
	for (adjEntry adj : edges) {
		isEdgeOnCycle[adj->theEdge()] = true;
	}
}

void Cycle::popBackNode() {
	OGDF_ASSERT(!nodes.empty());
	node no = nodes.back();
	nodes.popBack();
	isOnCycle[no] = false;
}

void Cycle::popFrontNode() {
	OGDF_ASSERT(!nodes.empty());
	node no = nodes.front();
	nodes.popFront();
	isOnCycle[no] = false;
}

void Cycle::popBackEdge() {
	OGDF_ASSERT(!edges.empty());
	adjEntry adj = edges.back();
	edges.popBack();
	isEdgeOnCycle[adj->theEdge()] = false;
}

void Cycle::popFrontEdge() {
	OGDF_ASSERT(!edges.empty());
	adjEntry adj = edges.front();
	edges.popFront();
	isEdgeOnCycle[adj->theEdge()] = false;
}

void Cycle::pushFrontEdge(adjEntry adj) {
	edges.pushFront(adj);
	isEdgeOnCycle[adj->theEdge()] = true;
}

// expanding a cycle whose triangle has at least one tree edge
void Cycle::expandWithTreeEdge(node y, node v, node w, adjEntry vy, adjEntry yw) {
	// remove the current expand edge
	popFrontEdge();

	node removedNode = nullptr;

	// handle special case: if y is on the cycle already, we just need to cut off that entire triangle
	if (isOnCycle[y]) {
		if (tree->isInTree(vy->theEdge())) {
			popBackNode(); // remove v
			popBackEdge(); // remove y->v
			pushFrontEdge(yw);
			removedNode = v;
		} else if (tree->isInTree(yw->theEdge())) {
			popFrontNode(); // remove w
			popFrontEdge(); // remove w->y
			pushFrontEdge(vy);
			removedNode = w;
		} else {
			// Trying to expand tree edge, but no tree edge was present!
			OGDF_ASSERT(false);
		}

		// if root was removed, new root has to be y.
		if (removedNode == cycleRoot) {
			cycleRoot = y;
		}

		// in this scenario, we actually took a node out of the cycle and added it to the smaller side
		if (costClock >= costCounter) {
			costCounter++;
		} else {
			costClock++;
		}

		return;
	}

	// y was not on the cycle

	// depending on which edge was the tree edge, add y to beginning or end of cycle and add edges correspondingly
	if (tree->isInTree(vy->theEdge())) {
		nodes.pushBack(y);
		edges.pushBack(vy);
		edges.pushFront(yw); // the new startEdge has to be in the beginning of the list
	} else if (tree->isInTree(yw->theEdge())) {
		nodes.pushFront(y);
		edges.pushFront(yw); // order matters
		edges.pushFront(vy);
	} else {
		// Trying to expand tree edge, but no tree edge was present!
		OGDF_ASSERT(false);
	}

	// add y and the connecting edges to cycle
	isOnCycle[y] = true;
	isEdgeOnCycle[vy->theEdge()] = true;
	isEdgeOnCycle[yw->theEdge()] = true;

	if (tree->getLevelOfNode(y) < tree->getLevelOfNode(cycleRoot)) {
		cycleRoot = y;
	}

	if (costClock >= costCounter) {
		costClock--;
	} else {
		costCounter--;
	}
}

std::pair<node, node> Cycle::findPathToCycle(node& y, List<node>& pathNodes,
		List<adjEntry>& pathAdjEntries) const {
	List<node> nodesToRoot; // list of nodes down to the cycleRoot (starting at z)
	List<adjEntry> edgesToRoot; // list of adjEntries down to the cycleRoot (starting at z)

	int cycleRootLevel = tree->getLevelOfNode(this->cycleRoot);
	int yLevel = tree->getLevelOfNode(y);
	int levelDiff = cycleRootLevel - yLevel;

	node rootPathNode = cycleRoot;

	// compensate level difference between cycleRoot and y
	if (levelDiff > 0) {
		for (int i = 0; i < levelDiff; i++) {
			nodesToRoot.pushFront(rootPathNode);
			adjEntry adjToParent = tree->getAdjToParent(rootPathNode);
			edgesToRoot.pushFront(adjToParent->twin());
			rootPathNode = adjToParent->twinNode();
		}

	} else if (levelDiff < 0) {
		levelDiff = std::abs(levelDiff);

		int i = 0;
		while (!isOnCycle[y] && i < levelDiff) {
			pathNodes.pushBack(y);
			pathAdjEntries.pushBack(tree->getAdjToParent(y));
			y = tree->getParentOfNode(y);
			i++;
		}
	}

	// At this point, y and rootPathNode are on the same level.
	// Now, keep adding nodes to both parts of the path as long as
	// y hasn't met the cycle and they are not pointing at the same node.

	while (!isOnCycle[y] && y != rootPathNode) {
		// y-Path
		pathNodes.pushBack(y);
		pathAdjEntries.pushBack(tree->getAdjToParent(y));
		y = tree->getParentOfNode(y);

		// cycleRoot-Path, backwards
		nodesToRoot.pushFront(rootPathNode);
		adjEntry adjToParent = tree->getAdjToParent(rootPathNode);
		edgesToRoot.pushFront(adjToParent->twin());
		rootPathNode = adjToParent->twinNode();
	}

	// if y has met the cycle, we're happy and return, discarding the other path
	if (isOnCycle[y]) {
		return std::make_pair(y, cycleRoot);
	}

	// otherwise, the cycles have met, so now we need to add the node in which they met...
	pathNodes.pushBack(y);
	nodesToRoot.popBack(); // we don't actually want the root node on the path to itself

	// ...and glue the lists together
	pathNodes.conc(nodesToRoot);
	pathAdjEntries.conc(edgesToRoot);

	// not a mistake: if they met, the cycleRoot is the node at which they meet, and y now contains the highest node
	return std::make_pair(cycleRoot, y);
}

bool Cycle::findAlphaCycle(Cycle& cyc, const List<node>& pathNodes,
		const List<adjEntry>& pathAdjEntries, const node z, const node propRoot, const adjEntry yw,
		List<node>& oldNodes, List<adjEntry>& oldEdges) const {
	List<node> nodesOnC2;
	List<adjEntry> edgesOnC2;
	node rootOfC2;
	bool foundRoot =
			false; // whether we found the root node of the original cycle on our way from w to z

	oldEdges.popFront();

	// collect all nodes on the original cycle, from w to z
	node no = oldNodes.front();
	while (no != z) {
		nodesOnC2.pushBack(no);
		// if propRoot was actually just cycleRoot, and we come across cycleRoot on our way, use that as new root
		if (no == propRoot) {
			rootOfC2 = no;
			foundRoot = true;
		}
		edgesOnC2.pushBack(oldEdges.popFrontRet());

		oldNodes.popFront();
		no = oldNodes.front();
	}

	nodesOnC2.pushBack(z); // adding z
	rootOfC2 = foundRoot ? rootOfC2 : z;

	// adding all nodes on path in reversed order
	auto pathNodeIt = pathNodes.crbegin();
	auto pathEdgeIt = pathAdjEntries.crbegin();
	for (; pathNodeIt != pathNodes.crend(); ++pathNodeIt, ++pathEdgeIt) {
		nodesOnC2.pushBack(*pathNodeIt);
		edgesOnC2.pushBack((*pathEdgeIt)->twin());
	}
	edgesOnC2.pushFront(yw); // adding final non-tree edge

	// finally, if the proposed root was really better than cycleRoot, it's clear that this should be the new root
	if (tree->getLevelOfNode(propRoot) < tree->getLevelOfNode(this->cycleRoot)) {
		rootOfC2 = propRoot;
	}

	cyc.init(nodesOnC2, edgesOnC2, rootOfC2);

	return foundRoot;
}

void Cycle::findBetaCycle(Cycle& cyc, const List<node>& pathNodes,
		const List<adjEntry>& pathAdjEntries, const node z, const node propRoot, const adjEntry vy,
		List<node>& oldNodes, List<adjEntry>& oldEdges, bool foundRootOnAlpha) const {
	List<node> nodesOnC1(pathNodes);
	List<adjEntry> edgesOnC1(pathAdjEntries);
	node rootOfC1;

	nodesOnC1.conc(oldNodes);
	edgesOnC1.conc(oldEdges);

	edgesOnC1.pushFront(vy); // adding final non-tree edge

	// figuring out what the cycle root is
	if (tree->getLevelOfNode(propRoot) < tree->getLevelOfNode(cycleRoot)) {
		// we walked up the tree, propRoot is on path from y to z and is definitely the new root, end of discussion
		rootOfC1 = propRoot;
	} else {
		if (foundRootOnAlpha) {
			// cycleRoot was on alpha cycle already, so root of c1 has to be z
			rootOfC1 = z;
		} else {
			// alpha cycle never found cycleRoot, so it's actually on this one, so it's our root
			rootOfC1 = cycleRoot;
		}
	}

	cyc.init(nodesOnC1, edgesOnC1, rootOfC1);
}

Cycle Cycle::expandWithoutTreeEdges(node y, const node v, const node w, const adjEntry vy,
		const adjEntry yw) {
	List<node> pathNodes; // nodes on path from y to z
	List<adjEntry> pathAdjEntries; // adjEntries on path from y to z

	// 1. find all nodes on the path from y to z
	auto res = findPathToCycle(y, pathNodes, pathAdjEntries);
	node z = std::get<0>(res);
	node propRoot = std::get<1>(res);

	// making sure that isClockwise is correct, just out of paranoia
	OGDF_ASSERT(isClockwise == (costClock >= costCounter));

	// 3. create empty Cycles
	Cycle cycle2(tree, isClockwise);
	Cycle cycle1(tree, isClockwise);

	// 4. fill alpha cycle
	bool foundRootOnAlpha =
			findAlphaCycle(cycle2, pathNodes, pathAdjEntries, z, propRoot, yw, nodes, edges);

	// 5. fill beta cycle
	findBetaCycle(cycle1, pathNodes, pathAdjEntries, z, propRoot, vy, nodes, edges, foundRootOnAlpha);


	// walk over the insides of the two cycles - we are actually getting the insides because their iterator walks in the right direction
	auto c1It = cycle1.begin(), c2It = cycle2.begin();
	for (; c1It != cycle1.end() && c2It != cycle2.end(); ++c1It, ++c2It) {
		cycle1.increaseCost(*c1It, isClockwise);
		cycle2.increaseCost(*c2It, isClockwise);
	}

	// at this point, we know the inner cost of one of the two cycles, we just have to figure out which one is ready
	Cycle& completedCycle = c1It == cycle1.end() ? cycle1 : cycle2;
	Cycle& otherCycle = c1It == cycle1.end() ? cycle2 : cycle1;

	int ownCost = getInsideCost();

	if (isClockwise) {
		completedCycle.costCounter =
				tree->getGraphSize() - completedCycle.nodes.size() - completedCycle.costClock;

		otherCycle.costClock = ownCost - pathNodes.size() - completedCycle.costClock;
		otherCycle.costCounter =
				tree->getGraphSize() - otherCycle.nodes.size() - otherCycle.costClock;

		OGDF_ASSERT(costClock == completedCycle.costClock + pathNodes.size() + otherCycle.costClock);

		completedCycle.isClockwise = completedCycle.costClock >= completedCycle.costCounter;
		otherCycle.isClockwise = otherCycle.costClock >= otherCycle.costCounter;
	} else {
		completedCycle.costClock =
				tree->getGraphSize() - completedCycle.nodes.size() - completedCycle.costCounter;

		otherCycle.costCounter = ownCost - pathNodes.size() - completedCycle.costCounter;
		otherCycle.costClock =
				tree->getGraphSize() - otherCycle.nodes.size() - otherCycle.costCounter;

		OGDF_ASSERT(costCounter
				== completedCycle.costCounter + pathNodes.size() + otherCycle.costCounter);

		completedCycle.isClockwise = completedCycle.costClock >= completedCycle.costCounter;
		otherCycle.isClockwise = otherCycle.costClock >= otherCycle.costCounter;
	}

	if (completedCycle.getInsideCost() > otherCycle.getInsideCost()) {
		return std::move(otherCycle);
	} else {
		return std::move(completedCycle);
	}
}

void Cycle::increaseCost(adjEntry adj, bool clockwise) {
	// check if the corresponding edge is a tree-edge, we only care about those
	if (tree->isInTree(adj->theEdge())) {
		OGDF_ASSERT(!isOnCycle[adj->twinNode()]);

		int addedCost = 0;

		// if the edge leads from the cycle Root further up into the tree, we need to calculate its weight differently
		if (adj->theNode() == cycleRoot
				&& adj
						== tree->getAdjToParent(
								cycleRoot)) { // first comparison is redundant, but maybe faster to check?
			addedCost = tree->getGraphSize() - tree->getDescendantsOfNode(adj->theNode());
		} else {
			addedCost = tree->getDescendantsOfNode(adj->twinNode());
		}

		if (clockwise) {
			costClock += addedCost;
		} else {
			costCounter += addedCost;
		}
	}
}

void Cycle::collectChildrenOfNode(const node no, NodeArray<bool>& marked, List<node>& list,
		bool useRoot) const {
	if ((no != tree->getRoot()) || useRoot) {
		list.pushBack(tree->getGraph()->original(no));
	}

	marked[no] = true;

	for (node child : tree->getChildrenOfNode(no)) {
		if (!marked[child]) {
			collectChildrenOfNode(child, marked, list);
		}
	}
}

void Cycle::computeCosts() {
	costClock = 0;
	costCounter = 0;

	// walk over all adjEntries pointing to the inside of the cycle
	for (auto it = this->begin(), end = this->end(); it != end; ++it) {
		increaseCost(*it, isClockwise); // and add them to the appropriate cost
	}

	if (isClockwise) {
		costCounter = tree->getGraphSize() - nodes.size() - costClock;
	} else {
		costClock = tree->getGraphSize() - nodes.size() - costCounter;
	}

	isClockwise = costClock
			>= costCounter; // cycle is clockwise = the clockwise iteration of edges yields the bigger side
}

Cycle::Cycle(BFSTree* bfsTree, edge startEdge) : tree {bfsTree}, isClockwise {true} {
	isOnCycle.init(*tree->getGraph(), false);
	isEdgeOnCycle.init(*tree->getGraph(), false);

	// 'first' is the node with the deeper level
	node n1 = startEdge->source();
	node n2 = startEdge->target();
	node first = tree->getLevelOfNode(n1) > tree->getLevelOfNode(n2) ? n1 : n2;
	node second = tree->getLevelOfNode(n1) > tree->getLevelOfNode(n2) ? n2 : n1;

	// two halves of the cycle: up to common root from first, down from common node to second
	List<node> nodesA;
	List<node> nodesB;

	int diff = tree->getLevelOfNode(first) - tree->getLevelOfNode(second); // difference in levels

	// compensate difference first
	for (int i = 0; i < diff; i++) {
		nodesA.pushBack(first);
		isOnCycle[first] = true;
		first = tree->getParentOfNode(first);
	}

	// first and second are on the same level now, add both while they haven't met
	while (first != second) {
		nodesA.pushBack(first);
		isOnCycle[first] = true;
		first = tree->getParentOfNode(first);

		nodesB.pushFront(second);
		isOnCycle[second] = true;
		second = tree->getParentOfNode(second);
	}
	nodesA.pushBack(first); // only add the root of the cycle once
	isOnCycle[first] = true;
	cycleRoot = first;

	// compile a list of all cycle edges
	adjEntry startAdjEnt =
			startEdge->target() == nodesA.front() ? startEdge->adjSource() : startEdge->adjTarget();
	edges.pushBack(startAdjEnt);
	for (node no : nodesA) {
		if (no == cycleRoot) {
			continue; // just ignore this one
		}
		edges.pushBack(tree->getAdjToParent(no));
	}
	for (node no : nodesB) {
		if (no == cycleRoot) {
			continue; // just ignore this one (can never happen, anyway)
		}
		edges.pushBack(tree->getAdjToParent(no)->twin());
	}
	for (adjEntry adj : edges) {
		isEdgeOnCycle[adj->theEdge()] = true;
	}

	// the two branches have met in their common root, compile final list
	nodesA.conc(nodesB);
	nodes = nodesA;

	// compute the cost on both sides of the cycle
	computeCosts();
}

int Cycle::getInsideCost() const { return isClockwise ? costClock : costCounter; }

int Cycle::getOutsideCost() const { return isClockwise ? costCounter : costClock; }

adjEntry Cycle::getCurrentExpandEdge() const {
	OGDF_ASSERT(!edges.empty());
	return edges.front();
}

Cycle Cycle::expandCycle() {
	// step 1: find other two edges of inner triangle, in correct direction
	adjEntry vy; // edge from start of currentExpandEdge to y = v->y
	adjEntry yw; // edge from y to end of currentExpandEdge = y->w

	if (isClockwise) {
		vy = getCurrentExpandEdge()->cyclicPred();
		yw = vy->twin()->cyclicPred();
	} else {
		vy = getCurrentExpandEdge()->cyclicSucc();
		yw = vy->twin()->cyclicSucc();
	}
	node y = yw->theNode();
	node v = getCurrentExpandEdge()->theNode();
	node w = getCurrentExpandEdge()->twinNode();

	OGDF_ASSERT(vy->theNode() == v && vy->twinNode() == y);
	OGDF_ASSERT(yw->theNode() == y && yw->twinNode() == w);
	OGDF_ASSERT(!(tree->isInTree(vy->theEdge())
			&& tree->isInTree(yw->theEdge()))); // this should be impossible, if this happens, I messed up


	// step 2: inspect triangle v-y-w on the bigger side of the currentExpandEdge
	// scenario 1: one of the two edges is a tree edge
	if (tree->isInTree(vy->theEdge()) || tree->isInTree(yw->theEdge())) {
		expandWithTreeEdge(y, v, w, vy, yw);
		return std::move(*this);
	}

	// scenario 2: neither of the two edges is a tree edge
	return expandWithoutTreeEdges(y, v, w, vy, yw);
}

void Cycle::print() const {
	std::cout << "Nodes in cycle (" << nodes.size() << ") : ";
	for (node no : nodes) {
		std::cout << no << " ";
	}
	std::cout << std::endl;

	std::cout << "isClockwise: " << isClockwise << std::endl;
	std::cout << "cycle Root: " << cycleRoot << std::endl;
	std::cout << "Clockwise Cost: " << costClock << ", counter-Clockwise Cost: " << costCounter
			  << std::endl;
	OGDF_ASSERT(isClockwise == (costClock >= costCounter));
}

// fills supplied lists with separator = nodes on the cycle, first = inside, second = outside
void Cycle::fillLists(List<node>& separator, List<node>& first, List<node>& second, bool useRoot) {
	NodeArray<bool> marked; // contains for each node whether we have put it into one of the lists
	marked.init(*tree->getGraph(), false);
	if (!useRoot) {
		marked[tree->getRoot()] = true;
	}
	// step 1: nodes on the cycle are part of the separator
	for (const node no : nodes) {
		if ((no != tree->getRoot()) || useRoot) {
			separator.pushBack(tree->getGraph()->original(no));
		}
		marked[no] = true;
	}

	// precaution: enforcing that the inside is actually the bigger side
	OGDF_ASSERT(isClockwise == (costClock >= costCounter));

	// step 2: follow all tree-edges on the inside of the cycle and collect the found nodes in first list
	for (auto it = this->begin(), end = this->end(); it != end; ++it) {
		if (tree->isInTree((*it)->theEdge())) {
			if (it.isOutEdge()) {
				collectChildrenOfNode(tree->getRoot(), marked, first, useRoot);
			} else {
				collectChildrenOfNode((*it)->twinNode(), marked, first, useRoot);
			}
		}
	}

	// step 3: all remaining unmarked nodes belong into the second list
	for (const node no : tree->getGraph()->nodes) {
		if (!marked[no]) {
			second.pushBack(tree->getGraph()->original(no));
		}
	}
}

}
