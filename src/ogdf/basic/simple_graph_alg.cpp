/** \file
 * \brief Implementation of simple graph algorithms
 *
 * \author Carsten Gutwenger, Sebastian Leipert
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


#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/Stack.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/tuples.h>
#include <ogdf/basic/BoundedStack.h>


namespace ogdf {


//---------------------------------------------------------
// isLoopFree(), makeLoopFree()
// testing for self-loops, removing self-loops
//---------------------------------------------------------
bool isLoopFree(const Graph &G)
{
	for(edge e : G.edges)
		if(e->isSelfLoop()) return false;

	return true;
}


void makeLoopFree(Graph &G)
{
	edge e, eNext;
	for (e = G.firstEdge(); e; e = eNext) {
		eNext = e->succ();
		if (e->isSelfLoop()) G.delEdge(e);
	}
}


//---------------------------------------------------------
// isParallelFree(), makeParallelFree()
// testing for multi-edges, removing multi-edges
//---------------------------------------------------------

void parallelFreeSort(const Graph &G, SListPure<edge> &edges)
{
	G.allEdges(edges);

	BucketSourceIndex bucketSrc;
	edges.bucketSort(0,G.maxNodeIndex(),bucketSrc);

	BucketTargetIndex bucketTgt;
	edges.bucketSort(0,G.maxNodeIndex(),bucketTgt);
}


bool isParallelFree(const Graph &G)
{
	if (G.numberOfEdges() <= 1) return true;

	SListPure<edge> edges;
	parallelFreeSort(G,edges);

	SListConstIterator<edge> it = edges.begin();
	edge ePrev = *it, e;
	for(it = ++it; it.valid(); ++it, ePrev = e) {
		e = *it;
		if (ePrev->source() == e->source() && ePrev->target() == e->target())
			return false;
	}

	return true;
}


int numParallelEdges(const Graph &G)
{
	if (G.numberOfEdges() <= 1) return 0;

	SListPure<edge> edges;
	parallelFreeSort(G,edges);

	int num = 0;
	SListConstIterator<edge> it = edges.begin();
	edge ePrev = *it, e;
	for(it = ++it; it.valid(); ++it, ePrev = e) {
		e = *it;
		if (ePrev->source() == e->source() && ePrev->target() == e->target())
			++num;
	}

	return num;
}



//---------------------------------------------------------
// isParallelFreeUndirected(), makeParallelFreeUndirected()
// testing for (undirected) multi-edges, removing (undirected) multi-edges
//---------------------------------------------------------

void parallelFreeSortUndirected(const Graph &G,
	SListPure<edge> &edges,
	EdgeArray<int> &minIndex,
	EdgeArray<int> &maxIndex)
{
	G.allEdges(edges);

	for(edge e : G.edges) {
		int srcIndex = e->source()->index(), tgtIndex = e->target()->index();
		if (srcIndex <= tgtIndex) {
			minIndex[e] = srcIndex; maxIndex[e] = tgtIndex;
		} else {
			minIndex[e] = tgtIndex; maxIndex[e] = srcIndex;
		}
	}

	BucketEdgeArray bucketMin(minIndex), bucketMax(maxIndex);
	edges.bucketSort(0,G.maxNodeIndex(),bucketMin);
	edges.bucketSort(0,G.maxNodeIndex(),bucketMax);
}


bool isParallelFreeUndirected(const Graph &G)
{
	if (G.numberOfEdges() <= 1) return true;

	SListPure<edge> edges;
	EdgeArray<int> minIndex(G), maxIndex(G);
	parallelFreeSortUndirected(G,edges,minIndex,maxIndex);

	SListConstIterator<edge> it = edges.begin();
	edge ePrev = *it, e;
	for(it = ++it; it.valid(); ++it, ePrev = e) {
		e = *it;
		if (minIndex[ePrev] == minIndex[e] && maxIndex[ePrev] == maxIndex[e])
			return false;
	}

	return true;
}


int numParallelEdgesUndirected(const Graph &G)
{
	if (G.numberOfEdges() <= 1) return 0;

	SListPure<edge> edges;
	EdgeArray<int> minIndex(G), maxIndex(G);
	parallelFreeSortUndirected(G,edges,minIndex,maxIndex);

	int num = 0;
	SListConstIterator<edge> it = edges.begin();
	edge ePrev = *it, e;
	for(it = ++it; it.valid(); ++it, ePrev = e) {
		e = *it;
		if (minIndex[ePrev] == minIndex[e] && maxIndex[ePrev] == maxIndex[e])
			++num;
	}

	return num;
}



//---------------------------------------------------------
// isConnected(), makeConnected()
// testing connectivity, establishing connectivity
//---------------------------------------------------------

bool isConnected(const Graph &G)
{
	node v = G.firstNode();
	if (v == nullptr) return true;

	int count = 0;
	NodeArray<bool> visited(G,false);
	BoundedStack<node> S(G.numberOfNodes());

	S.push(v);
	visited[v] = true;
	while(!S.empty()) {
		v = S.pop();
		++count;

		for(adjEntry adj : v->adjEntries) {
			node w = adj->twinNode();
			if(!visited[w]) {
				visited[w] = true;
				S.push(w);
			}
		}
	}

	return count == G.numberOfNodes();
}


void makeConnected(Graph &G, List<edge> &added)
{
	added.clear();
	if (G.numberOfNodes() == 0) return;
	NodeArray<bool> visited(G,false);
	BoundedStack<node> S(G.numberOfNodes());

	node pred = nullptr;
	for(node u : G.nodes)
	{
		if (visited[u]) continue;

		node vMinDeg = u;
		int  minDeg  = u->degree();

		S.push(u);
		visited[u] = true;

		while(!S.empty())
		{
			node v = S.pop();

			for(adjEntry adj : v->adjEntries) {
				node w = adj->twinNode();
				if(!visited[w]) {
					visited[w] = true;
					S.push(w);

					int wDeg = w->degree();
					if (wDeg < minDeg) {
						vMinDeg = w;
						minDeg  = wDeg;
					}
				}
			}
		}

		if (pred)
			added.pushBack(G.newEdge(pred,vMinDeg));
		pred = vMinDeg;
	}
}


int connectedComponents(const Graph &G, NodeArray<int> &component)
{
	int nComponent = 0;
	component.fill(-1);

	StackPure<node> S;

	for(node v : G.nodes) {
		if (component[v] != -1) continue;

		S.push(v);
		component[v] = nComponent;

		while(!S.empty()) {
			node w = S.pop();
			for(adjEntry adj : w->adjEntries) {
				node x = adj->twinNode();
				if (component[x] == -1) {
					component[x] = nComponent;
					S.push(x);
				}
			}
		}

		++nComponent;
	}

	return nComponent;
}

//return the isolated nodes too, is used in incremental layout
int connectedIsolatedComponents(const Graph &G, List<node> &isolated,
								NodeArray<int> &component)
{
	int nComponent = 0;
	component.fill(-1);

	StackPure<node> S;

	for(node v : G.nodes) {
		if (component[v] != -1) continue;

		S.push(v);
		component[v] = nComponent;

		while(!S.empty()) {
			node w = S.pop();
			if (w->degree() == 0) isolated.pushBack(w);
			for(adjEntry adj : w->adjEntries) {
				node x = adj->twinNode();
				if (component[x] == -1) {
					component[x] = nComponent;
					S.push(x);
				}
			}
		}

		++nComponent;
	}

	return nComponent;
}//connectedIsolated


//---------------------------------------------------------
// isBiconnected(), makeBiconnected()
// testing biconnectivity, establishing biconnectivity
//---------------------------------------------------------
//! Build up a dfs-tree starting from the node root by assigning each reachable
//! node in the graph a discovery time (number) and a parent.
/**
 * @param root is the node which should be the root of the dfs-tree.
 * @param number is assigned the number (discovery time) for each node.
 *        The number of root is firstNr, the number of unvisited nodes is 0.
 * @param parent is assigned the parent in the dfs-tree for each node.
 * @param childNr is assigned the number of children for each node.
 * @param revS is assigned all visited nodes such that the top element of revS
 *        is the node that was visited last.
 * @param directed should be set to true if the directionality of edges should
 *        be respected.
 * @param firstNr is the index > 0 at which the numbering of the nodes starts.
 * @return the number of visited nodes, i.e., nodes in the dfs-tree.
 */
int buildDfsTree(const node &root,
		NodeArray<int> &number,
		NodeArray<node> &parent,
		NodeArray<int> &childNr,
		ArrayBuffer<node> &revS,
		bool directed = false,
		int firstNr = 1)
{
	OGDF_ASSERT(firstNr > 0);

	ArrayBuffer<node> S;
	S.push(root);

	int numCount = firstNr;
	childNr.fill(0);

	// Build up search tree and note discovery time and parent for each node.
	while (!S.empty()) {
		node v = S.popRet();

		// Ignore nodes that were already visited.
		if (number[v] != 0) {
			continue;
		}

		revS.push(v);

		// Set correct discovery time for v.
		number[v] = numCount++;

		// For all adjacent nodes w of v:
		for (adjEntry adj : v->adjEntries) {
			if (directed && adj->theEdge()->source() != v) {
				continue;
			}

			node w = adj->twinNode();

			// If w has not been visited yet:
			// Push it on the stack, remember its parent and number of children.
			if (number[w] == 0) {
				S.push(w);

				// If a parent was determined previously, revert that.
				if (parent[w] != nullptr) {
					childNr[parent[w]]--;
				}

				parent[w] = v;
				childNr[v]++;
			}
		}
	}

	return numCount - firstNr;
}

//! Find cut vertices and potential edges that could be added to turn the cut
//! vertices into non-cut vertices.
/**
 * The algorithm is applied to the graph whose nodes were pushed to the
 * ArrayBuffer revS. number, parent and revS can be obtained with buildDfsTree.
 *
 * @param number contains the number (discovery time) for each node.
 *        The number of root is 1, the number of unvisited nodes is 0.
 * @param parent contains the parent in the dfs-tree for each node.
 * @param revS contains the nodes of a graph such that the node that was visited
 *        last during the dfs-traversal is its top element.
 * @param cutVertices is assigned the cut vertices of the graph.
 * @param addEdges is assigned the tuples of nodes which have to be connected in
 *        order to turn each cut vertex into a non-cut vertex.
 * @param only_one should be set to true if the search should stop after finding
 *        one cut vertex, to false if all cut vertices should be found.
 * @return true if the graph contains at least one cut vertex, false otherwise.
 */
bool findCutVertices(NodeArray<int> &number,
		NodeArray<node> &parent,
		ArrayBuffer<node> &revS,
		ArrayBuffer<node> &cutVertices,
		ArrayBuffer<Tuple2<node,node>> &addEdges,
		bool only_one)
{
	NodeArray<int> lowpt(number);

	// Go backwards through the dfs-tree:
	// Calculate the lowpoint for each node and test for cut vertices.
	while (!revS.empty()) {
		node v = revS.popRet();
		node firstChild = nullptr;

		// For all adjacent nodes w of v:
		for (adjEntry adj : v->adjEntries) {
			node w = adj->twinNode();

			// Ignore self-loops and the parent of v.
			if (v == w || parent[v] == w) {
				continue;
			}

			// If v->w is a backedge in the dfs-tree, update v's lowpoint.
			if (number[v] > number[w] ) {
				if (lowpt[v] > number[w]) {
					lowpt[v] = number[w];
				}
			} else {
				// If w is v's child in the dfs-tree, update v's lowpoint.
				if (parent[w] == v) {
					if (lowpt[v] > lowpt[w]) {
						lowpt[v] = lowpt[w];
					}

					// See whether w is v's first son.
					if (firstChild == nullptr) {
						firstChild = w;
					}

					// Non-root v is a cut vertex if lowpt[w] >= number[v].
					if (parent[v] != nullptr && lowpt[w] >= number[v]) {
						// Suggest to add an edge between w and v's parent.
						cutVertices.push(v);
						addEdges.push(Tuple2<node,node>(w, parent[v]));

						if (only_one) {
							return true;
						}
					}

					// Root v is a cut vertex if v has two or more children.
					if  (parent[v] == nullptr && w != firstChild) {
						// Suggest to add an edge between those children.
						cutVertices.push(v);
						addEdges.push(Tuple2<node,node>(w, firstChild));

						if (only_one) {
							return true;
						}
					}
				}
			}
		}
	}

	return !cutVertices.empty();
}

bool isBiconnected(const Graph &G, node &cutVertex)
{
	cutVertex = nullptr;

	if (G.empty()) {
		return true;
	}

	NodeArray<int> number(G,0);        // discovery times
	NodeArray<node> parent(G,nullptr); // parents in the dfs tree
	ArrayBuffer<node> revS;            // nodes of the dfs tree in reverse order

	// Build the dfs-tree and get the number of visited nodes.
	NodeArray<int> childNr(G);
	int numCount = buildDfsTree(G.firstNode(), number, parent, childNr, revS);

	// If the graph is not connected, return false.
	if (numCount != G.numberOfNodes()) {
		return false;
	}

	// If there are cut vertices in the graph, return false, else true.
	ArrayBuffer<node> cutVertices;
	ArrayBuffer<Tuple2<node,node>> addEdges;
	if (findCutVertices(number, parent, revS, cutVertices, addEdges, true)) {
		cutVertex = cutVertices.top();
		return false;
	} else {
		return true;
	}
}

void makeBiconnected(Graph &G, List<edge> &added)
{
	if (G.empty()) {
		return;
	}

	makeConnected(G, added);

	NodeArray<int> number(G,0);        // discovery times
	NodeArray<node> parent(G,nullptr); // parents in the dfs tree
	ArrayBuffer<node> revS;            // nodes of the dfs tree in reverse order

	// Build the dfs-tree.
	NodeArray<int> childNr(G);
	buildDfsTree(G.firstNode(), number, parent, childNr, revS);

	// Find all cut vertices.
	ArrayBuffer<node> cutVertices;
	ArrayBuffer<Tuple2<node,node>> addEdges;
	findCutVertices(number, parent, revS, cutVertices, addEdges, false);

	// Add a new edge for each cut vertex to make the graph biconnected.
	for (Tuple2<node,node> nodes : addEdges) {
		added.pushBack(G.newEdge(nodes.x1(), nodes.x2()));
	}
}


//---------------------------------------------------------
// biconnectedComponents()
// computing biconnected components
//---------------------------------------------------------
static void dfsBiconComp (const Graph &G,
	node v,
	node father,
	NodeArray<int> &number,
	NodeArray<int> &lowpt,
	StackPure<node> &called,
	EdgeArray<int> &component,
	int &nNumber,
	int &nComponent)
{
	lowpt[v] = number[v] = ++nNumber;
	called.push(v);

	for(adjEntry adj : v->adjEntries) {
		node w = adj->twinNode();
		if (v == w) continue; // ignore self-loops

		if (number[w] == 0) {

			dfsBiconComp(G,w,v,number,lowpt,called,component,
				nNumber,nComponent);

			if (lowpt[w] < lowpt[v]) lowpt[v] = lowpt[w];

		} else {

			if (number[w] < lowpt[v]) lowpt[v] = number[w];
		}
	}

	if (father && (lowpt[v] == number[father])) {
		node w;
		do {
			w = called.top(); called.pop();

			for(adjEntry adj : w->adjEntries) {
				if (number[w] > number[adj->twinNode()])
					component[adj->theEdge()] = nComponent;
			}
		} while (w != v);

		++nComponent;
	}
}


int biconnectedComponents(const Graph &G, EdgeArray<int> &component)
{
	if (G.empty()) return 0;

	StackPure<node> called;
	NodeArray<int> number(G,0);
	NodeArray<int> lowpt(G);
	int nNumber = 0, nComponent = 0, nIsolated = 0;

	for(node v : G.nodes) {
		if (number[v] == 0) {
			bool isolated = true;
			for(adjEntry adj : v->adjEntries) {
				if (adj->twinNode() != v) {
					isolated = false; break;
				}
			}
			if (isolated)
				++nIsolated;
			else
				dfsBiconComp(G,v,nullptr,number,lowpt,called,component,
					nNumber,nComponent);
		}
	}

	return nComponent + nIsolated;
}


//---------------------------------------------------------
// isTriconnected()
// testing triconnectivity
//---------------------------------------------------------
bool isTriconnectedPrimitive(const Graph &G, node &s1, node &s2)
{
	s1 = s2 = nullptr;

	if (isConnected(G) == false)
		return false;

	if (isBiconnected(G,s1) == false)
		return false;

	if (G.numberOfNodes() <= 3)
		return true;

	// make a copy of G
	GraphCopySimple GC(G);

	// for each node v in G, we test if G \ v is biconnected
	for(node v : G.nodes)
	{
		node vC = GC.copy(v), wC;

		// store adjacent nodes
		SListPure<node> adjacentNodes;
		for(adjEntry adj : vC->adjEntries) {
			wC = adj->twinNode();
			// forget self-loops (vC would no longer be in GC!)
			if (wC != vC)
				adjacentNodes.pushBack(wC);
		}

		GC.delNode(vC);

		// test for biconnectivity
		if(isBiconnected(GC,wC) == false) {
			s1 = v; s2 = GC.original(wC);
			return false;
		}

		// restore deleted node with adjacent edges
		vC = GC.newNode(v);
		SListConstIterator<node> it;
		for(node uC : adjacentNodes)
			GC.newEdge(vC,uC);
	}

	return true;
}


//--------------------------------------------------------------------------
// triangulate()
//--------------------------------------------------------------------------
void triangulate(Graph &G)
{
	OGDF_ASSERT(isSimple(G));

	CombinatorialEmbedding E(G);

	OGDF_ASSERT(E.consistencyCheck());

	adjEntry succ, succ2, succ3;
	NodeArray<int> marked(E.getGraph(), 0);

	for(node v : E.getGraph().nodes) {
		marked.init(E.getGraph(), 0);

		for(adjEntry adj : v->adjEntries) {
			marked[adj->twinNode()] = 1;
		}

		// forall faces adj to v
		for(adjEntry adj : v->adjEntries) {
			succ = adj->faceCycleSucc();
			succ2 = succ->faceCycleSucc();

			if (succ->twinNode() != v && adj->twinNode() != v) {
				while (succ2->twinNode() != v) {
					if (marked[succ2->theNode()] == 1) {
						// edge e=(x2,x4)
						succ3 = succ2->faceCycleSucc();
						E.splitFace(succ, succ3);
					}
					else {
						// edge e=(v=x1,x3)
						edge e = E.splitFace(adj, succ2);
						marked[succ2->theNode()] = 1;

						// old adj is in wrong face
						adj = e->adjSource();
					}
					succ = adj->faceCycleSucc();
					succ2 = succ->faceCycleSucc();
				}
			}
		}
	}
}


//--------------------------------------------------------------------------
// isAcyclic(), isAcyclicUndirected(), makeAcyclic(), makeAcyclicByReverse()
// testing acyclicity, establishing acyclicity
//--------------------------------------------------------------------------
bool isAcyclic(const Graph &G, List<edge> &backedges)
{
	backedges.clear();

	NodeArray<int> number(G,0);        // discovery times
	NodeArray<node> parent(G,nullptr); // parents in the dfs tree
	NodeArray<int> childNr(G);         // number of children in the dfs tree
	ArrayBuffer<node> revS;

	ArrayBuffer<node> leaves;          // leaves of the dfs tree
	NodeArray<int> completion(G,0);    // completion times
	int complCount = 0;
	int numCount = 0;

	// For all unvisited nodes:
	for (node v : G.nodes) {
		if (number[v] == 0) {
			// Build the dfs-tree starting at v.
			numCount += buildDfsTree(v, number, parent, childNr, revS, true, numCount+1);

			// Get all leaves of the dfs-tree.
			while (!revS.empty()) {
				node w = revS.popRet();
				if (childNr[w] == 0) {
					leaves.push(w);
				}
			}

			node lastParent = parent[leaves.top()];

			// Go through leaves of the dfs-tree.
			while (!leaves.empty()) {
				node w = leaves.top();

				// If the new leaf is a child of the same parent as before,
				// assign it a completion time and pop it from the stack.
				if (parent[w] == lastParent) {
					completion[w] = complCount++;
					leaves.pop();

					// The last parent has now one child less. If it has no
					// children anymore, push it as a new leaf on the stack.
					if (lastParent != nullptr) {
						childNr[lastParent]--;
						if (childNr[lastParent] == 0) {
							leaves.push(lastParent);
							lastParent = parent[lastParent];
						}
					}
				} else {
					// Else just continue with the next leaves and their parent.
					lastParent = parent[w];
				}
			}
		}
	}

	// Remember backedges.
	for(edge e : G.edges) {
		node src = e->source();
		node tgt = e->target();

		if (number[src] >= number[tgt] && completion[src] <= completion[tgt]) {
			backedges.pushBack(e);
		}
	}

	return backedges.empty();
}


bool isAcyclicUndirected(const Graph &G, List<edge> &backedges)
{
	backedges.clear();

	NodeArray<int> number(G,0);        // discovery times
	NodeArray<node> parent(G,nullptr); // parents in the dfs tree
	ArrayBuffer<node> S;
	int numCount = 0;

	// For all unvisited nodes:
	for (node v : G.nodes) {
		if (number[v] == 0) {
			// Start depth first search at v.
			S.push(v);
			while (!S.empty()) {
				node w = S.popRet();

				// Ignore nodes that were already visited.
				if (number[w] != 0) {
					continue;
				}

				// Set correct discovery time for w.
				number[w] = ++numCount;
				bool parentSeen = false;

				// For all adjacent nodes u of w:
				for (adjEntry adj : w->adjEntries) {
					node u = adj->twinNode();

					// If u has not been visited yet,
					// push it on the stack and remember its parent.
					if (number[u] == 0) {
						S.push(u);
						parent[u] = w;
					} else if (parent[w] == u && !parentSeen) {
						// The first time you see w's parent, it is no backedge.
						parentSeen = true;
					} else if (w != u || adj->isSource()) {
						// Collect backedges (self-loops only in one direction).
						backedges.pushBack(adj->theEdge());
					}
				}
			}
		}
	}

	return backedges.empty();
}


void makeAcyclic(Graph &G)
{
	List<edge> backedges;
	isAcyclic(G,backedges);

	for(edge e : backedges)
		G.delEdge(e);
}


void makeAcyclicByReverse(Graph &G)
{
	List<edge> backedges;
	isAcyclic(G,backedges);

	for(edge e : backedges)
		if (!e->isSelfLoop()) G.reverseEdge(e);
}


//---------------------------------------------------------
// hasSingleSource(), hasSingleSink()
// testing for single source/sink
//---------------------------------------------------------
bool hasSingleSource(const Graph& G, node &s)
{
	s = nullptr;

	for(node v : G.nodes) {
		if (v->indeg() == 0) {
			if (s != nullptr) {
				s = nullptr;
				return false;
			} else s = v;
		}
	}
	return (G.empty() || s != nullptr);
}


bool hasSingleSink(const Graph& G, node &t)
{
	t = nullptr;

	for(node v : G.nodes) {
		if (v->outdeg() == 0) {
			if (t != nullptr) {
				t = nullptr;
				return false;
			} else t = v;
		}
	}
	return (G.empty() || t != nullptr);
}


//---------------------------------------------------------
// isStGraph()
// true <=> G is st-graph, i.e., is acyclic, contains exactly one source s
//   and one sink t, and the edge (s,t); returns single source s and single
//   sink t if contained (otherwise they are set to 0), and edge st if
//   contained (otherwise 0)
//---------------------------------------------------------
bool isStGraph(const Graph &G, node &s, node &t, edge &st)
{
	st = nullptr;

	hasSingleSource(G,s);
	hasSingleSink  (G,t);

	if (s == nullptr || t == nullptr || isAcyclic(G) == false) {
		s = t = nullptr;
		return false;
	}

	for(adjEntry adj : s->adjEntries) {
		edge e = adj->theEdge();

		if (e->target() == t) {
			st = e;
			break;
		}
	}

	return st != nullptr;
}


//---------------------------------------------------------
// topologicalNumbering()
// computes a topological numbering of an acyclic graph
//---------------------------------------------------------

void topologicalNumbering(const Graph &G, NodeArray<int> &num)
{
	BoundedStack<node> S(G.numberOfNodes());
	NodeArray<int> indeg(G);

	for(node v : G.nodes)
		if((indeg[v] = v->indeg()) == 0)
			S.push(v);

	int count = 0;
	while(!S.empty()) {
		node v = S.pop();
		num[v] = count++;

		for(adjEntry adj : v->adjEntries) {
			node u = adj->theEdge()->target();
			if(u != v) {
				if(--indeg[u] == 0)
					S.push(u);
			}
		}
	}
}

/**
 * Computes the strongly connected componets using tarjans algorithm.
 *
 * @param graph the graph to work on
 * @param v the node to start the traversal with
 * @param components maps nodes to resulting components
 * @param marked whether this node has been visited yet
 * @param lowLinks smallest node reachable from any node
 * @param set stack for maintaining the current components node
 * @param pNextIndex for numbering each node uniquely
 * @param pNextComponent for number each component uniquely
 */
void computeStrongComponents(
  const Graph &graph,
  const node v,
  NodeArray<int> &components,
  NodeArray<bool> &marked,
  NodeArray<int> &lowLinks,
  BoundedStack<node> &set,
  int *pNextIndex,
  int *pNextComponent)
{
	marked[v] = true;
	int min = lowLinks[v] = (*pNextIndex)++;
	set.push(v);

	for(adjEntry adj : v->adjEntries) {
		edge e = adj->theEdge();
		if(v == e->source()) {
			node w = e->target();
			if(!marked[w]) {
				computeStrongComponents(graph, w, components, marked, lowLinks, set, pNextIndex, pNextComponent);
			}
			if(lowLinks[w] < min) {
				min = lowLinks[w];
			}
		}
	}
	if(min < lowLinks[v]) {
		lowLinks[v] = min;
	} else {
		node w;
		do {
			w = set.pop();
			components[w] = *pNextComponent;
			lowLinks[w] = graph.numberOfNodes();
		} while(w != v);
		(*pNextComponent)++;
	}
}

int strongComponents(const Graph &graph, NodeArray<int> &components)
{
	int nNodes = graph.numberOfNodes();
	int result = 0;
	if(nNodes > 0) {
		if(nNodes == 1) {
			result = 1;
			components[graph.firstNode()] = 0;
		} else {
			NodeArray<bool> marked(graph, false);
			NodeArray<int> lowLinks(graph, -1);
			BoundedStack<node> set(nNodes);
			int nextIndex = 0;

			for(node v : graph.nodes) {
				if(!marked[v]) {
					computeStrongComponents(graph, v, components, marked, lowLinks, set, &nextIndex, &result);
				}
			}
		}
	}
	return result;
}

//---------------------------------------------------------
// makeBimodal()
// makes the DiGraph bimodal such that all embeddings of the
// graph are bimodal embeddings!
//---------------------------------------------------------

void makeBimodal(Graph &G, List<edge> &newEdge)
{
	List<node> nodes;
	G.allNodes(nodes);

	ListIterator<node> it_n = nodes.begin();
	while (it_n.valid()) {
		node v = *it_n;
		if (v->indeg() < 2 || v->outdeg() < 2) {
			++it_n; continue;
		}
		List<adjEntry> newOrder;
		for (adjEntry adj : v->adjEntries) {
			if (adj->theEdge()->target() == v)
				newOrder.pushFront(adj);
			else
				newOrder.pushBack(adj);
		}
		G.sort(v, newOrder);

		ListIterator<adjEntry> it = newOrder.begin();
		while ((*it)->theEdge()->target() == v)
			++it;
		node newNode = G.splitNode(newOrder.front(), *it);
		for (adjEntry adj : newNode->adjEntries) {
			if (adj->theEdge()->target() == newNode) {
				newEdge.pushBack(adj->theEdge());
				break;
			}
		}
		++it_n;
	}
}

//---------------------------------------------------------
// isFreeForest()
// testing if graph represents a free forest
//---------------------------------------------------------

bool isFreeForest(const Graph &G)
{
	NodeArray<bool> visited(G,false);

	for(node vFirst : G.nodes) {
		if (visited[vFirst]) continue;

		StackPure<Tuple2<node,node> > S;
		S.push(Tuple2<node,node>(vFirst,nullptr));

		while (!S.empty()) {
			Tuple2<node,node> t = S.pop();
			node v      = t.x1();
			node parent = t.x2();

			visited[v] = true;

			for(adjEntry adj : v->adjEntries) {
				node w = adj->twinNode();

				// skip edge to parent, but only once!
				if(w == parent) {
					parent = nullptr;
					continue;
				}

				if(visited[w] == true)
					return false;

				S.push(Tuple2<node,node>(w,v));
			}
		}
	}

	return true;
}


//---------------------------------------------------------
// isForest(), isArborescence()
// testing if graph represents a forest / an arborescence
//---------------------------------------------------------
static bool dfsIsForest (node v,
	NodeArray<bool> &visited,
	NodeArray<bool> &mark)
{
	SListPure<node> sons;

	visited[v] = true;

	for(adjEntry adj : v->adjEntries) {
		node w = adj->theEdge()->target();
		if (w != v && !mark[w]) {
			mark[w] = true;
			sons.pushBack(w);
		}
	}

	for(node w : sons)
		mark[w] = false;

	while(!sons.empty()) {
		node w = sons.front();
		sons.popFront();

		if (visited [w] || dfsIsForest(w,visited,mark) == false)
			return false;
	}

	return true;
}

bool isForest(const Graph& G, List<node> &roots)
{
	roots.clear();
	if (G.empty()) return true;

	NodeArray<bool> visited(G,false), mark(G,false);

	for(node v : G.nodes)
		if (v->indeg() == 0) {
			roots.pushBack(v);
			if (dfsIsForest(v,visited,mark) == false)
				return false;
		}

	for(node v : G.nodes)
		if (!visited[v]) return false;

	return true;
}


bool isArborescence (const Graph& G, node &root)
{
	List<node> roots;

	if (isForest(G,roots) && roots.size() == 1) {
		root = roots.front(); return true;
	}
	return false;
}

bool isRegular(const Graph& G) {
	if (G.numberOfEdges() == 0) {
		return true;
	}
	return isRegular(G, G.firstNode()->degree());
}

bool isRegular(const Graph& G, int d) {
	for (auto n: G.nodes) {
		if (n->degree() != d) {
			return false;
		}
	}
	return true;
}


} // end namespace ogdf
