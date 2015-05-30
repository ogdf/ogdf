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

		for(adjEntry adj : v->adjEdges) {
			node w = adj->twinNode();
			if(!visited[w]) {
				visited[w] = true;
				S.push(w);
			}
		}
	}

	return (count == G.numberOfNodes());
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

			for(adjEntry adj : v->adjEdges) {
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
			edge e;
			forall_adj_edges(e,w) {
				node x = e->opposite(w);
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
			edge e;
			forall_adj_edges(e,w) {
				node x = e->opposite(w);
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
static node dfsIsBicon (const Graph &G, node v, node father,
	NodeArray<int> &number, NodeArray<int> &lowpt, int &numCount)
{
	node first_son = nullptr;

	lowpt[v] = number[v] = ++numCount;

	edge e;
	forall_adj_edges(e,v) {
		node w = e->opposite(v);
		if (v == w) continue; // ignore self-loops

		if (number[w] == 0) {
			if (first_son == nullptr) first_son = w;

			node cutVertex = dfsIsBicon(G,w,v,number,lowpt,numCount);
			if (cutVertex) return cutVertex;

			// is v cut vertex ?
			if (lowpt[w] >= number[v] && (w != first_son || father != nullptr))
				return v;

			if (lowpt[w] < lowpt[v]) lowpt[v] = lowpt[w];

		} else {

			if (number[w] < lowpt[v]) lowpt[v] = number[w];
		}
	}

	return nullptr;
}


bool isBiconnected(const Graph &G, node &cutVertex)
{
	if (G.empty()) return true;

	NodeArray<int> number(G,0);
	NodeArray<int> lowpt(G);
	int numCount = 0;

	cutVertex = dfsIsBicon(G,G.firstNode(),nullptr,number,lowpt,numCount);

	return (numCount == G.numberOfNodes() && cutVertex == nullptr);
}


static void dfsMakeBicon (Graph &G,
	node v, node father,
	NodeArray<int> &number,
	NodeArray<int> &lowpt,
	int &numCount,
	List<edge> &added)
{
	node predSon = nullptr;

	lowpt[v] = number[v] = ++numCount;

	edge e;
	forall_adj_edges(e,v) {
		node w = e->opposite(v);
		if (v == w) continue; // ignore self-loops

		if (number[w] == 0) {

			dfsMakeBicon(G,w,v,number,lowpt,numCount,added);

			// is v cut vertex ?
			if (lowpt[w] >= number[v]) {
				if (predSon == nullptr && father != nullptr)
					added .pushBack(G.newEdge(w,father));

				else if (predSon != nullptr)
					added.pushBack(G.newEdge(w,predSon));
			}

			if (lowpt[w] < lowpt[v]) lowpt[v] = lowpt[w];
			predSon = w;

		} else {

			if (number[w] < lowpt[v]) lowpt[v] = number[w];
		}
	}
}


void makeBiconnected(Graph &G, List<edge> &added)
{
	if (G.empty()) return;

	makeConnected(G,added);

	NodeArray<int> number(G,0);
	NodeArray<int> lowpt(G);
	int numCount = 0;

	dfsMakeBicon(G,G.firstNode(),nullptr,number,lowpt,numCount,added);
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

	edge e;
	forall_adj_edges(e,v) {
		node w = e->opposite(v);
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

			forall_adj_edges(e,w) {
				if (number[w] > number[e->opposite(w)])
					component[e] = nComponent;
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
			edge e;
			forall_adj_edges(e,v)
				if (!e->isSelfLoop()) {
					isolated = false; break;
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
		edge eC;
		forall_adj_edges(eC,vC) {
			wC = eC->opposite(vC);
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

		for(adjEntry adj : v->adjEdges) {
			marked[adj->twinNode()] = 1;
		}

		// forall faces adj to v
		for(adjEntry adj : v->adjEdges) {
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
void dfsIsAcyclic(const Graph &G,
	node v,
	NodeArray<int> &number,
	NodeArray<int> &completion,
	int &nNumber,
	int &nCompletion)
{
	number[v] = ++nNumber;

	edge e;
	forall_adj_edges(e,v) {
		node w = e->target();

		if (number[w] == 0)
			dfsIsAcyclic(G,w,number,completion,nNumber,nCompletion);
	}

	completion[v] = ++nCompletion;
}


void dfsIsAcyclicUndirected(const Graph &G,
	node v,
	NodeArray<int> &number,
	int &nNumber,
	List<edge> &backedges)
{
	number[v] = ++nNumber;

	for(adjEntry adj : v->adjEdges) {
		node w = adj->twinNode();
		if (number[w] == 0) {
			dfsIsAcyclicUndirected(G,w,number,nNumber,backedges);
		} else {
			if (number[w] > number[v]) {
				backedges.pushBack(adj->theEdge());
			}
		}
	}
}


bool isAcyclic(const Graph &G, List<edge> &backedges)
{
	backedges.clear();

	NodeArray<int> number(G,0), completion(G);
	int nNumber = 0, nCompletion = 0;

	for(node v : G.nodes)
		if (number[v] == 0)
			dfsIsAcyclic(G,v,number,completion,nNumber,nCompletion);

	for(edge e : G.edges) {
		node src = e->source(), tgt = e->target();

		if (number[src] >= number[tgt] && completion[src] <= completion[tgt])
			backedges.pushBack(e);
	}

	return backedges.empty();
}


bool isAcyclicUndirected(const Graph &G, List<edge> &backedges)
{
	backedges.clear();
	int nNumber = 0;
	NodeArray<int> number(G,0);

	for(node v : G.nodes) {
		if (number[v] == 0) {
			dfsIsAcyclicUndirected(G,v,number,nNumber,backedges);
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

	edge e;
	forall_adj_edges(e,s) {
		if (e->target() == t) {
			st = e;
			break;
		}
	}

	return (st != nullptr);
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

		edge e;
		forall_adj_edges(e,v) {
			node u = e->target();
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

	for(adjEntry adj : v->adjEdges) {
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
		for (adjEntry adj : v->adjEdges) {
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
		for (adjEntry adj : newNode->adjEdges) {
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

			for(adjEntry adj : v->adjEdges) {
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

	edge e;
	forall_adj_edges(e,v) {
		node w = e->target();
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



} // end namespace ogdf
