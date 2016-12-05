/** \file
 * \brief Implements the class NonPlanarCore.
 *
 * \author Carsten Gutwenger, Mirko Wagner
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

#include <ogdf/planarity/NonPlanarCore.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Queue.h>
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/FaceArray.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/DualGraph.h>
#include <ogdf/graphalg/MinSTCut.h>
#include <ogdf/graphalg/MaxFlowSTPlanarDigraph.h>

namespace ogdf {

NonPlanarCore::NonPlanarCore(const Graph &G, bool nonPlanarityGuaranteed) :
  m_pOriginal(&G),
  m_T(G),
  m_orig(m_graph),
  m_real(m_graph, nullptr),
  m_cost(m_graph),
  m_mincut(m_graph),
  m_mapV(m_graph),
  m_mapE(m_graph),
  m_underlyingGraphs(m_graph),
  m_sNode(m_graph),
  m_tNode(m_graph) {
	call(G, nullptr, nullptr, nonPlanarityGuaranteed);
}

NonPlanarCore::NonPlanarCore(const Graph &G, const EdgeArray<int> &weight, MaxFlowModule<int> &maxFlowModule,
                             bool nonPlanarityGuaranteed) :
  m_pOriginal(&G),
  m_T(G),
  m_orig(m_graph),
  m_real(m_graph, nullptr),
  m_cost(m_graph),
  m_mincut(m_graph),
  m_mapV(m_graph),
  m_mapE(m_graph),
  m_underlyingGraphs(m_graph),
  m_sNode(m_graph),
  m_tNode(m_graph) {
	call(G, &weight, &maxFlowModule, nonPlanarityGuaranteed);
}

NonPlanarCore::NonPlanarCore(const Graph &G, const EdgeArray<int> &weight, bool nonPlanarityGuaranteed) :
  m_pOriginal(&G),
  m_T(G),
  m_orig(m_graph),
  m_real(m_graph, nullptr),
  m_cost(m_graph),
  m_mincut(m_graph),
  m_mapV(m_graph),
  m_mapE(m_graph),
  m_underlyingGraphs(m_graph),
  m_sNode(m_graph),
  m_tNode(m_graph) {
	call(G, &weight, nullptr, nonPlanarityGuaranteed);
}

void NonPlanarCore::call(const Graph &G, const EdgeArray<int> *weight, MaxFlowModule<int> *maxFlowModule,
                         bool nonPlanarityGuaranteed) {
	if(!nonPlanarityGuaranteed){
		if(isPlanar(G)){
			return;
		}
	}
	OGDF_ASSERT(!isPlanar(G));
	OGDF_ASSERT(isBiconnected(G));

	// mark tree nodes in the core
	NodeArray<bool> mark;
	markCore(mark);

	NodeArray<node> map(G, nullptr);
	NodeArray<node> mapAux(G, nullptr);
	const Graph &tree = m_T.tree();

	for(node v : tree.nodes){
		if(mark[v] == false){
			continue;
		}

		Skeleton &S = m_T.skeleton(v);
		edge lastCreatedEdge = nullptr;

		for(edge e : S.getGraph().edges){
			node src = S.original(e->source());
			node tgt = S.original(e->target());

			if(tgt == src){
				continue;
			}

			if(map[src] == nullptr){
				m_orig[map[src] = m_graph.newNode()] = S.original(e->source());
			}

			if(map[tgt] == nullptr){
				m_orig[map[tgt] = m_graph.newNode()] = S.original(e->target());
			}

			if(S.isVirtual(e)){
				node w = S.twinTreeNode(e);

				if(mark[w] == false){
					// new virtual edge in core graph
					lastCreatedEdge = m_graph.newEdge(map[src], map[tgt]);
					m_real[lastCreatedEdge] = nullptr;
					traversingPath(S, e, m_mincut[lastCreatedEdge], mapAux,
					               lastCreatedEdge, weight, maxFlowModule);
				}
			} else {
				// new real edge in core graph
				lastCreatedEdge = m_graph.newEdge(map[src], map[tgt]);
				m_real[lastCreatedEdge] = S.realEdge(e);
				traversingPath(S, e, m_mincut[lastCreatedEdge], mapAux,
				               lastCreatedEdge, weight, maxFlowModule);
			}
		}
	}

	if(weight != nullptr)
	{
		for(edge e : m_graph.edges)
		{
			int cost(0);
			for(auto cutEdge : m_mincut[e])
			{
				cost += (*weight)[cutEdge.e];
			}
			m_cost[e] = cost;
		}
	} else {
		for(edge e : m_graph.edges)
		{
			m_cost[e] = m_mincut[e].size();
		}
	}

	List<node> allNodes;
	m_graph.allNodes(allNodes);

	// The while loop is used to eliminate multiedges from the core. We're pruning P-Nodes.
	List<edge> winningMultiEdges;
	List<edge> losingMultiEdges;
	getAllMultiedges(winningMultiEdges, losingMultiEdges);
	while(!winningMultiEdges.empty() && !losingMultiEdges.empty()){
		edge winningMultiEdge = winningMultiEdges.popFrontRet();
		edge losingMultiEdge = losingMultiEdges.popFrontRet();
#ifdef OGDF_DEBUG
		int sizeOfWinCut = m_mincut[winningMultiEdge].size();
		int sizeOfLosCut = m_mincut[losingMultiEdge].size();
#endif

		glue(winningMultiEdge, losingMultiEdge);
		glueMincuts(winningMultiEdge, losingMultiEdge);

		OGDF_ASSERT(m_mincut[winningMultiEdge].size() == sizeOfWinCut + sizeOfLosCut);
		delete m_underlyingGraphs[losingMultiEdge];
		delete m_mapV[losingMultiEdge];
		delete m_mapE[losingMultiEdge];
		m_real[winningMultiEdge] = nullptr;
		m_real[losingMultiEdge] = nullptr;
		m_graph.delEdge(losingMultiEdge);
	}
	// The for loop is used to eliminate deg 2 nodes from the core. We're pruning S-Nodes.
	for(node v : allNodes){
		if(v->degree() != 2){
			continue;
		}
		edge outEdge = v->firstAdj()->theEdge();
		edge inEdge = v->lastAdj()->theEdge();

		if(m_cost[inEdge] > m_cost[outEdge])
		{
			swap(inEdge, outEdge);
		}
		// We join the embeddings of the underlying embeddings/graphs of both edges
		// so that outEdge gets integrated into inEdge
		glue(inEdge, outEdge);

		m_real[inEdge] = nullptr;
		m_real[outEdge] = nullptr;

		adjEntry adjSource = inEdge->adjSource()->cyclicSucc();
		adjEntry adjTarget = (outEdge->target() == v ? outEdge->adjSource()->cyclicSucc()
		                                             : outEdge->adjTarget()->cyclicSucc());
		if(inEdge->target() != v){
			adjSource = adjTarget;
			adjTarget = inEdge->adjTarget()->cyclicSucc();
		}
		m_graph.move(inEdge, adjSource, ogdf::Direction::before, adjTarget, ogdf::Direction::before);
		delete m_underlyingGraphs[outEdge];
		delete m_mapV[outEdge];
		delete m_mapE[outEdge];
		m_graph.delNode(v);
	}


	if(nonPlanarityGuaranteed){
		OGDF_ASSERT(!isPlanar(m_graph));
	}
}

void NonPlanarCore::markCore(NodeArray<bool> &mark){
	const Graph &tree = m_T.tree();

	// We mark every tree node that belongs to the core
	mark.init(tree, true);  // start with all nodes and unmark planar leaves
	NodeArray<int> degree(tree);

	Queue<node> Q;

	for(node v : tree.nodes){
		degree[v] = v->degree();
		if(degree[v] <= 1){ // also append deg-0 node (T has only one node)
			Q.append(v);
		}
	}

	while(!Q.empty()){
		node v = Q.pop();

		// if v has a planar skeleton
		if(m_T.typeOf(v) != SPQRTree::RNode ||
		   isPlanar(m_T.skeleton(v).getGraph()) == true){
			mark[v] = false; // unmark this leaf

			node w = nullptr;
			for(adjEntry adj : v->adjEntries){
				node x = adj->twinNode();
				if(mark[x] == true){
					w = x;
					break;
				}
			}

			if(w != nullptr){
				--degree[w];
				if(degree[w] == 1){
					Q.append(w);
				}
			}
		}
	}
}

struct OGDF_EXPORT QueueEntry {
	QueueEntry(node p, node v) : m_parent(p), m_current(v){ }

	node m_parent;
	node m_current;
};


void NonPlanarCore::traversingPath(const Skeleton &Sv, edge eS, List<CutEdge> &path, NodeArray<node> &mapV,
                                   edge coreEdge, const EdgeArray<int> *weight_src, MaxFlowModule<int> *maxFlowModule)
{
	List<CutEdge> &currPath = path;

	// Build the graph representing the planar st-component
	Graph *h_pointer = new Graph();
	Graph &H = *h_pointer;

	EdgeArray<edge> *mapE_src_pointer = new EdgeArray<edge>(H, nullptr);
	EdgeArray<edge> &mapE_src = *mapE_src_pointer;
	NodeArray<node> *mapV_src_pointer = new NodeArray<node>(H, nullptr);
	NodeArray<node> &mapV_src = *mapV_src_pointer;
	SListPure<node> nodes;

	if(Sv.isVirtual(eS)){
		Queue<QueueEntry> Q;
		Q.append(QueueEntry(Sv.treeNode(), Sv.twinTreeNode(eS)));

		while(!Q.empty()){
			QueueEntry x = Q.pop();
			node parent = x.m_parent;
			node current = x.m_current;

			const Skeleton &S = m_T.skeleton(current);
			for(edge e : S.getGraph().edges){
				if(S.isVirtual(e)){
					continue;
				}

				node src = S.original(e->source());
				node tgt = S.original(e->target());

				if(mapV[src] == nullptr){
					nodes.pushBack(src);
					mapV[src] = H.newNode();
					mapV_src[mapV[src]] = src;
				}
				if(mapV[tgt] == nullptr){
					nodes.pushBack(tgt);
					mapV[tgt] = H.newNode();
					mapV_src[mapV[tgt]] = tgt;
				}

				edge e_new = H.newEdge(mapV[src], mapV[tgt]);
				mapE_src[e_new] = S.realEdge(e);
				OGDF_ASSERT(mapE_src[e_new]->source() != nullptr);
			}

			for(adjEntry adj : current->adjEntries){
				node w = adj->twinNode();
				if(w != parent){
					Q.append(QueueEntry(current, w));
				}
			}
		}
	} else {
		nodes.pushBack(Sv.original(eS->source()));
		nodes.pushBack(Sv.original(eS->target()));
		mapV[Sv.original(eS->source())] = H.newNode();
		mapV_src[mapV[Sv.original(eS->source())]] = Sv.original(eS->source());
		mapV[Sv.original(eS->target())] = H.newNode();
		mapV_src[mapV[Sv.original(eS->target())]] = Sv.original(eS->target());
		mapE_src[H.newEdge(mapV[Sv.original(eS->source())], mapV[Sv.original(eS->target())])] = Sv.realEdge(eS);
	}
	// add st-edge
	edge e_st = H.newEdge(mapV[Sv.original(eS->source())], mapV[Sv.original(eS->target())]);
	m_sNode[coreEdge] = mapV[Sv.original(eS->source())];
	m_tNode[coreEdge] = mapV[Sv.original(eS->target())];

	// Compute planar embedding of H
#ifdef OGDF_DEBUG
	bool ok =
#endif
			planarEmbed(H);
	OGDF_ASSERT(ok);
	CombinatorialEmbedding E(H);

	// we rearange the adj Lists of s and t, so that adj(e_st) is the first adj
	List<adjEntry> adjListFront;
	List<adjEntry> adjListBack;
	e_st->source()->allAdjEntries(adjListFront);
	if(adjListFront.front() != e_st->adjSource()){
		adjListFront.split(adjListFront.search(e_st->adjSource()), adjListFront, adjListBack);
		adjListFront.concFront(adjListBack);
		H.sort(e_st->source(), adjListFront);
	}

	e_st->target()->allAdjEntries(adjListFront);
	if(adjListFront.front() != e_st->adjTarget()){
		adjListFront.split(adjListFront.search(e_st->adjTarget()), adjListFront, adjListBack, ogdf::Direction::before);
		adjListFront.concFront(adjListBack);
		H.sort(e_st->target(), adjListFront);
	}

	if(Sv.isVirtual(eS)){
		if(weight_src != nullptr)
		{
			GraphCopy doubleEdgedGraphCopy(H);
			node s = doubleEdgedGraphCopy.copy(e_st->source());
			node t = doubleEdgedGraphCopy.copy(e_st->target());
			doubleEdgedGraphCopy.delEdge(doubleEdgedGraphCopy.copy(e_st));
			List<edge> edges;
			doubleEdgedGraphCopy.allEdges(edges);
			EdgeArray<edge> originalEdge(doubleEdgedGraphCopy, nullptr);
			for(edge e : edges){
				edge revEdge = doubleEdgedGraphCopy.newEdge(e->target(), e->source());
				doubleEdgedGraphCopy.move(revEdge, e->adjTarget(), ogdf::before, e->adjSource(), ogdf::after);
				originalEdge[e] = originalEdge[revEdge] = doubleEdgedGraphCopy.original(e);
			}
			EdgeArray<int> flow(doubleEdgedGraphCopy, 0);
			EdgeArray<int> weight(doubleEdgedGraphCopy, 1);
			for(edge e : doubleEdgedGraphCopy.edges){
				weight[e] = (*weight_src)[mapE_src[originalEdge[e]]];
				OGDF_ASSERT(weight[e] > 0);
			}

			if(maxFlowModule){
				maxFlowModule->init(doubleEdgedGraphCopy);
				maxFlowModule->computeFlow(weight, s, t, flow);
			} else {
				MaxFlowSTPlanarDigraph<int> maxFlow(doubleEdgedGraphCopy);
				maxFlow.computeFlow(weight, s, t, flow);
			}

			MinSTCut<int> minSTCut;
			minSTCut.call(weight, flow, s, t);

			StackPure<edge> stack;
			EdgeArray<bool> visited(H, false);
			for(adjEntry adj : s->adjEntries){
				if(adj->theEdge()->adjTarget() != adj){
					stack.push(adj->theEdge());
				}
			}
			while(!stack.empty()){
				edge e = stack.pop();
				if(visited[originalEdge[e]]){
					continue;
				}
				node v = e->target();
				visited[originalEdge[e]] = true;
				if(minSTCut.isFrontCutEdge(e)){
					currPath.pushBack(CutEdge(mapE_src[originalEdge[e]], doubleEdgedGraphCopy.copy(originalEdge[e]) == e));
				} else {
					for(adjEntry adj : v->adjEntries){
						if(adj->theEdge()->adjTarget() != adj){
							stack.push(adj->theEdge());
						}
					}
				}
			}
		} else {
			DualGraph dual(E);
			node s = dual.dualNode(E.rightFace(e_st->adjSource()));
			node t = dual.dualNode(E.leftFace(e_st->adjSource()));

			NodeArray<edge> spPred(dual, nullptr);
			EdgeArray<node> prev(dual, nullptr);
			EdgeArray<bool> direction(dual, true);
			QueuePure<edge> queue;
			for(adjEntry adj : s->adjEntries) {
				if(dual.primalEdge(adj->theEdge()) != e_st) {
					queue.append(adj->theEdge());
					prev[adj->theEdge()] = s;
				}
			}
			// actual search (using bfs on directed dual)
			for(; ;){
				// next candidate edge
				edge eCand = queue.pop();
				bool dir = (eCand->source() == prev[eCand]);
				node v = (dir ? eCand->target() : eCand->source());

				if(eCand->source() == t || eCand->target() == t){
				}

				// leads to an unvisited node?
				if(spPred[v] == nullptr){
					// yes, then we set v's predecessor in search tree
					spPred[v] = eCand;
					direction[eCand] = dir;

					// have we reached t ...
					if(v == t){
						// ... then search is done.
						// constructed list of used edges (translated to crossed
						// edges entries in G) from t back to s (including first
						// and last!)

						do {
							edge eDual = spPred[v];
							OGDF_ASSERT(eDual != nullptr);
							// this should be the right way round
							edge eG = mapE_src[dual.primalEdge(eDual)];
							OGDF_ASSERT(eG != nullptr);
							currPath.pushBack(CutEdge(eG, !direction[eDual]));
							v = prev[eDual];
						} while(v != s);

						break;
					}

					// append next candidate edges to queue
					// (all edges leaving v)
					for(adjEntry adj : v->adjEntries){
						if(prev[adj->theEdge()] == nullptr){
							queue.append(adj->theEdge());
							prev[adj->theEdge()] = v;
						}
					}
				}
			}
		}
	}
	else {
		OGDF_ASSERT(Sv.realEdge(eS) != nullptr);
		currPath.pushFront(CutEdge(Sv.realEdge(eS), true));
	}
	H.delEdge(e_st);

	m_underlyingGraphs[coreEdge] = h_pointer;
	m_mapE[coreEdge] = mapE_src_pointer;
	m_mapV[coreEdge] = mapV_src_pointer;
#ifdef OGDF_DEBUG
	for(node v : H.nodes){
		OGDF_ASSERT(mapV_src[v] != nullptr);
	}
	for(edge e : H.edges){
		OGDF_ASSERT(mapE_src[e] != nullptr);
	}
#endif

	for(node v : nodes){
		mapV[v] = nullptr;
	}
}

void NonPlanarCore::getAllMultiedges(List<edge> &winningEdges, List<edge> &losingEdges){
	winningEdges.clear();
	losingEdges.clear();
	SListPure<edge> edges;
	EdgeArray<int> minIndex(m_graph), maxIndex(m_graph);
	parallelFreeSortUndirected(m_graph, edges, minIndex, maxIndex);

	SListConstIterator<edge> it = edges.begin();
	edge ePrev = *it, e;
	for(it = ++it; it.valid(); ++it, ePrev = e){
		e = *it;
		if(minIndex[ePrev] == minIndex[e] && maxIndex[ePrev] == maxIndex[e]){
			winningEdges.pushFront(ePrev);
			losingEdges.pushFront(e);
		}
	}
}

void NonPlanarCore::glue(edge eWinner, edge eLoser){
	GlueMap map(eWinner, eLoser, *this);

	// true iff this glueing is about a PNode (so a glueing at two common nodes)
	bool thisIsAboutAPNode = false;
	if((eLoser->source() == eWinner->source() && eLoser->target() == eWinner->target()) ||
	   (eLoser->target() == eWinner->source() && eLoser->source() == eWinner->target())){
		thisIsAboutAPNode = true;
	}

	// find the s- and t-nodes in their skeleton for both of the edges
	node sWinner = m_sNode[eWinner];
	node tWinner = m_tNode[eWinner];
	node sLoser = m_sNode[eLoser];
	node tLoser = m_tNode[eLoser];

	bool sameDirection = !(eWinner->source() == eLoser->target() && eWinner->target() == eLoser->source());

	// we get a list of all nodes of the loser graph
	List<node> allNodesButSt;
	map.getLoserGraph().allNodes(allNodesButSt);

#ifdef OGDF_DEBUG
	bool ok = true;
#endif

	// for both s and t of eLoser we check if it's either the s or the t node of eWinner
	// if one of it is, we delete it from the list of nodes, that are to be added
	// otherwise it stays in 'allNodesButSt' to be added later
	if(eLoser->source() == eWinner->source() || eLoser->source() == eWinner->target()){
#ifdef OGDF_DEBUG
		ok =
#endif
				allNodesButSt.removeFirst(sLoser);
		OGDF_ASSERT(ok);
		if(eLoser->source() == eWinner->source()){
			map.mapLoserToWinnerNode(sLoser, sWinner);
		} else {
			map.mapLoserToWinnerNode(sLoser, tWinner);
		}
		OGDF_ASSERT(!allNodesButSt.removeFirst(sLoser));
	}
	if(eLoser->target() == eWinner->source() || eLoser->target() == eWinner->target()){
#ifdef OGDF_DEBUG
		ok =
#endif
				allNodesButSt.removeFirst(tLoser);
		OGDF_ASSERT(ok);
		if(eLoser->target() == eWinner->source()){
			map.mapLoserToWinnerNode(tLoser, sWinner);
		} else {
			map.mapLoserToWinnerNode(tLoser, tWinner);
		}
		OGDF_ASSERT(!allNodesButSt.removeFirst(tLoser));
	}


	// insert the remaining nodes of the loser graph into the winner graph
	for(node v : allNodesButSt){
		map.mapLoserToNewWinnerNode(v);
	}

	// insert all edges of the loser graph into the the winner graph
	for(edge e : map.getLoserGraph().edges){
		map.mapLoserToNewWinnerEdge(e);
	}

	// reorder the adjEntries of every node of the loser graph in the winner graph,
	// to match the embedding in the loser graph
	List<node> allNodes = allNodesButSt;
	allNodes.pushBack(sLoser);
	allNodes.pushBack(tLoser);
	for(node v : allNodes){
		map.reorder(v, sameDirection, (tLoser == v && thisIsAboutAPNode));
	}
	if(!thisIsAboutAPNode){
		if(eWinner->source() == eLoser->source()){
			m_sNode[eWinner] = map.getWinnerNodeOfLoserNode(tLoser);
		}
		if(eWinner->target() == eLoser->source()){
			m_tNode[eWinner] = map.getWinnerNodeOfLoserNode(tLoser);
		}
		if(eWinner->source() == eLoser->target()){
			m_sNode[eWinner] = map.getWinnerNodeOfLoserNode(sLoser);
		}
		if(eWinner->target() == eLoser->target()){
			m_tNode[eWinner] = map.getWinnerNodeOfLoserNode(sLoser);
		}
	}
}

void NonPlanarCore::retransform(const GraphCopy &planarCore, GraphCopy &planarGraph){
#ifdef OGDF_DEBUG
	GraphCopy copyCore(planarCore);
	copyCore.removePseudoCrossings();
	OGDF_ASSERT(copyCore.numberOfNodes() == planarCore.numberOfNodes());
#endif
	m_endGraph = &planarGraph;
	m_planarCore = &planarCore;
	OGDF_ASSERT(m_planarCore->genus() == 0);
	m_endGraph->clear();
	m_endGraph->createEmpty(*m_pOriginal);
	List<node> allNodes;
	m_pOriginal->allNodes(allNodes);
	EdgeArray<edge> eCopy(*m_pOriginal, nullptr);
	m_endGraph->initByNodes(allNodes, eCopy);

#ifdef OGDF_DEBUG
	for(node v : m_endGraph->nodes){
		List<adjEntry> adjEntries;
		v->allAdjEntries(adjEntries);
		OGDF_ASSERT(v->degree() == adjEntries.size());
	}
#endif

	// For every node of the core we rearrange the adjacency order of the corresponding endGraph node
	// according to the planarized core.
	for(node v : m_planarCore->nodes){
		if(m_planarCore->isDummy(v)){
			continue;
		}
		List<adjEntry> pcOrder;
		v->allAdjEntries(pcOrder);
		List<adjEntry> newOrder;
		node coreNode = m_planarCore->original(v);
		OGDF_ASSERT(coreNode != nullptr);
		for(adjEntry adjPC : v->adjEntries) {
			edge coreEdge = m_planarCore->original(adjPC->theEdge());
			EdgeArray<edge> &mapE = *m_mapE[coreEdge];
			NodeArray<node> &mapV = *m_mapV[coreEdge];
			node stNode = (mapV[m_sNode[coreEdge]] == original(coreNode) ? m_sNode[coreEdge] : m_tNode[coreEdge]);
			// find the node of emb which represents the same node v represents
			for(adjEntry adjEmb : stNode->adjEntries){
				if(adjEmb->theEdge()->source() == adjEmb->theNode()){
					newOrder.pushBack(m_endGraph->copy(mapE[adjEmb->theEdge()])->adjSource());
				} else {
					newOrder.pushBack(m_endGraph->copy(mapE[adjEmb->theEdge()])->adjTarget());
				}
			}
		}
		m_endGraph->sort(m_endGraph->copy(original(coreNode)), newOrder);
	}
	List<node> splitdummies;
	for(edge e : m_graph.edges){
		List<edge> chain;
		// for every edge from the core we ensure, that the embedding of the subgraph it describes
		// is the same in both the original and the end graph
		importEmbedding(e);
		// reverse all cut edges, which are not the same direction as e
		normalizeCutEdgeDirection(e);
		// to ensure the right order of the inserted crossings, we insert dummy nodes
		// to split the edge in sections, each of which only has one crossing
		splitEdgeIntoSections(e, splitdummies);
	}

	// now we can start and insert the crossings of the planar core into the end graph.
	// a node represents a crossing if it's a dummy
	for(node v : m_planarCore->nodes){
		if(m_planarCore->isDummy(v)){
			inflateCrossing(v);
		}
	}
	OGDF_ASSERT(m_endGraph->genus() == 0);

	removeSplitdummies(splitdummies);
	for(edge e : m_graph.edges){
		normalizeCutEdgeDirection(e);
	}
}

void NonPlanarCore::normalizeCutEdgeDirection(edge coreEdge){
	for(auto cutE : m_mincut[coreEdge]) {
		if(!cutE.dir) {
			for(edge e : m_endGraph->chain(cutE.e)) {
				m_endGraph->reverseEdge(e);
			}
		}
	}
}

void NonPlanarCore::removeSplitdummies(List<node> &splitdummies){
	for(node v : splitdummies){
		edge eIn = v->firstAdj()->theEdge();
		edge eOut = v->lastAdj()->theEdge();
		if(eIn->target() == v){
			m_endGraph->unsplit(eIn, eOut);
		} else {
			m_endGraph->unsplit(eOut, eIn);
		}
	}
}

void NonPlanarCore::splitEdgeIntoSections(edge e, List<node> &splitdummies){
	List<edge> chain = m_planarCore->chain(e);
	int chainSize = chain.size();
	while(chainSize > 2){
		for(auto cutEdge : mincut(e)){
			splitdummies.pushBack(m_endGraph->split(m_endGraph->copy(cutEdge.e))->source());
		}
		chainSize--;
	}
#ifdef OGDF_DEBUG
	for(auto cutEdge : mincut(e)){
		if(chain.size() < 3){
			OGDF_ASSERT(m_endGraph->chain(cutEdge.e).size() == 1);
		} else {
			OGDF_ASSERT(m_endGraph->chain(cutEdge.e).size() == chain.size() - 1);
		}
		OGDF_ASSERT(m_endGraph->original(m_endGraph->chain(cutEdge.e).front()) == cutEdge.e);
	}
#endif
}

void NonPlanarCore::importEmbedding(edge e){
	const Graph &embG = *m_underlyingGraphs[e];
	// a map from the nodes of the emb to those in the end graph
	const EdgeArray<edge> &mapE_toOrig = *m_mapE[e];
	// bc the edges of the end graph are split for the crossing insertion,
	// a map of the emb might have more than one edge in the endgraph, we just
	// map the AdjEntries of both source and target of each edge of emb
	// to AdjEntries in the end graph
	const NodeArray<node> &mapV_toOrig = *m_mapV[e];
	AdjEntryArray<adjEntry> mapA_toFinal(embG, nullptr);
	for(auto it = mapE_toOrig.begin(); it != mapE_toOrig.end(); it++){
		OGDF_ASSERT(it.key() != nullptr);
		OGDF_ASSERT((*it) != nullptr);
		OGDF_ASSERT((*it)->graphOf() == m_pOriginal);
		mapA_toFinal[it.key()->adjSource()] = m_endGraph->chain(*it).front()->adjSource();
		mapA_toFinal[it.key()->adjTarget()] = m_endGraph->chain(*it).back()->adjTarget();
	}
	node s(m_sNode[e]), t(m_tNode[e]);
	List<node> nodesOfEmb;
	embG.allNodes(nodesOfEmb);
	// for every node of emb we order the adjEntries of the corresponding node
	// in the end graph, so that both match
	for(node v : nodesOfEmb){
		if(v == s || v == t){
			continue;
		}
		List<adjEntry> rightAdjOrder;
		for(adjEntry adj = v->firstAdj(); adj; adj = adj->succ()){
			rightAdjOrder.pushBack(mapA_toFinal[adj]);
		}
		m_endGraph->sort(m_endGraph->copy(mapV_toOrig[v]), rightAdjOrder);
	}
}

void NonPlanarCore::inflateCrossing(node v){
	// we want e1 and e2 to be these two edges
	//      ^
	//      |
	// e2-->v--->
	//      ^
	//      |
	//      e1
	edge e1 = v->firstAdj()->theEdge();
	while(e1->target() != v){
		e1 = e1->adjSource()->succ()->theEdge();
	}
	edge e2 = e1->adjTarget()->succ()->theEdge();
	while(e2->target() != v){
		e2 = e2->adjSource()->cyclicSucc()->theEdge();
	}
	if(e1 == e2->adjTarget()->cyclicSucc()->theEdge()){
		edge help = e1;
		e1 = e2;
		e2 = help;
	}
	OGDF_ASSERT(e2 == e1->adjTarget()->cyclicSucc()->theEdge());
	List<edge> e1cut;
	getMincut(e1, e1cut);
	List<edge> e2cut;
	getMincut(e2, e2cut);
	OGDF_ASSERT(e1 != e2);
	OGDF_ASSERT(e1cut.size() > 0);
	OGDF_ASSERT(e2cut.size() > 0);
	// the actual crossing insertion
	// for (auto it1 = e1cut.begin(); it1.valid(); it1++)
	for(int i = 0; i < e1cut.size(); i++){
		auto it1 = e1cut.get(i);
		edge crossingEdge = *it1;
		for(int j = 0; j < e2cut.size(); j++){
			auto it2 = e2cut.get(j);
			edge crossedEdge = *it2;
			m_endGraph->insertCrossing(*it1, crossedEdge, true);
			OGDF_ASSERT(crossedEdge == *it2);
			e2cut.insertAfter(crossedEdge, it2);
			e2cut.del(it2);
		}
		OGDF_ASSERT(crossingEdge != *it1);
		e1cut.insertAfter(crossingEdge, it1);
		e1cut.del(it1);
	}
}

void NonPlanarCore::getMincut(edge e, List<edge> &cut){
	OGDF_ASSERT(e->graphOf() == m_planarCore);

	cut.clear();
	// chain is a list of the edges of the planar core, that represent e
	List<edge> chain = m_planarCore->chain(m_planarCore->original(e));
	// this is the main part of this function:
	// as we know, the cut edges are split by splitdummies to partition the edge,
	// such that every crossing on the edge has its own section to be inserted into (denoted by pos)
	// cut_pre stores the first section for every cut edge
	List<CutEdge> cut_pre = mincut(m_planarCore->original(e));
	for(CutEdge eCut : cut_pre){
		OGDF_ASSERT(m_endGraph->chain(eCut.e).size() + 1 >= chain.size());
		// while iterating we have to differentiate between already inserted crossings and splitdummies
		// we can do that, by only counting the deg 2 nodes we pass while iterating through the chain of the cut edge
		auto it = m_endGraph->chain(eCut.e).begin();
		for(int i = 0; i < chain.pos(chain.search(e)); i++){
			it++;
			while((*it)->source()->degree() == 4){
				it++;
				OGDF_ASSERT(it.valid());
			}
		}
		cut.pushBack(*(it));
	}
	// cut is the result of this function
}

void NonPlanarCore::glueMincuts(edge eWinner, edge eLoser){
#ifdef OGDF_DEBUG
	if(eWinner->adjSource()->theNode() == eLoser->adjSource()->theNode()){
		OGDF_ASSERT(eWinner->adjSource()->theNode() == eLoser->adjSource()->theNode());
		OGDF_ASSERT(eWinner->adjTarget()->theNode() == eLoser->adjTarget()->theNode());
	}
	else {
		OGDF_ASSERT(eWinner->adjSource()->theNode() == eLoser->adjTarget()->theNode());
		OGDF_ASSERT(eWinner->adjTarget()->theNode() == eLoser->adjSource()->theNode());
	}
#endif
	List<CutEdge> wincut = m_mincut[eWinner];

	List<CutEdge> losecut = m_mincut[eLoser];

	if(eWinner->source() == eLoser->target()) {
		List<CutEdge> newLosecut;
		for(auto cutEit = losecut.begin(); cutEit != losecut.end(); cutEit++) {
			newLosecut.pushBack(CutEdge((*cutEit).e, !(*cutEit).dir));
		}
		losecut = newLosecut;
	}

	wincut.conc(losecut);
	m_mincut[eWinner] = wincut;
	m_cost[eWinner] += m_cost[eLoser];
}

NonPlanarCore::~NonPlanarCore(){
	for(auto pointer : m_mapE) {
		delete pointer;
	}
	for(auto pointer : m_mapV) {
		delete pointer;
	}
	for(auto pointer : m_underlyingGraphs) {
		delete pointer;
	}
}

GlueMap::GlueMap(edge eWinner, edge eLoser, NonPlanarCore &npc) : m_npc(npc), m_eWinner(eWinner),
                                                                  m_eLoser(eLoser){
	OGDF_ASSERT(m_eWinner != m_eLoser);

	OGDF_ASSERT(m_npc.m_underlyingGraphs[m_eLoser] != nullptr);
	OGDF_ASSERT(m_npc.m_underlyingGraphs[m_eWinner] != nullptr);

	m_gLoser = m_npc.m_underlyingGraphs[m_eLoser];
	m_gWinner = m_npc.m_underlyingGraphs[m_eWinner];

	OGDF_ASSERT(m_gWinner != m_gLoser);

	OGDF_ASSERT(m_npc.m_mapV[m_eWinner] != nullptr);
	OGDF_ASSERT(m_npc.m_mapV[m_eLoser] != nullptr);

	m_mapVwinner = m_npc.m_mapV[m_eWinner];
	m_mapVloser = m_npc.m_mapV[m_eLoser];

	OGDF_ASSERT(m_npc.m_mapE[m_eWinner] != nullptr);
	OGDF_ASSERT(m_npc.m_mapE[m_eLoser] != nullptr);

	m_mapEwinner = m_npc.m_mapE[m_eWinner];
	m_mapEloser = m_npc.m_mapE[m_eLoser];

	OGDF_ASSERT(m_mapEloser->graphOf() == m_gLoser);
	OGDF_ASSERT(m_mapVloser->graphOf() == m_gLoser);

	OGDF_ASSERT(m_mapEwinner->graphOf() == m_gWinner);
	OGDF_ASSERT(m_mapVwinner->graphOf() == m_gWinner);

	m_mapE_l2w = EdgeArray<edge>(*m_gLoser, nullptr);
	m_mapV_l2w = NodeArray<node>(*m_gLoser, nullptr);
}

void GlueMap::mapLoserToNewWinnerEdge(edge loser){
	edge newEdge = m_gWinner->newEdge(m_mapV_l2w[loser->source()], m_mapV_l2w[loser->target()]);
	m_mapE_l2w[loser] = newEdge;
	(*m_mapEwinner)[newEdge] = (*m_mapEloser)[loser];
}

void GlueMap::mapLoserToWinnerNode(node loser, node winner){
	m_mapV_l2w[loser] = winner;
	(*m_mapVwinner)[winner] = (*m_mapVloser)[loser];
}

void GlueMap::mapLoserToNewWinnerNode(node loser){
	node newNode = m_gWinner->newNode();
	m_mapV_l2w[loser] = newNode;
	(*m_mapVwinner)[newNode] = (*m_mapVloser)[loser];
}

void GlueMap::reorder(node vLoser, bool sameDirection, bool isTNodeOfPNode){
	node vWinner = m_mapV_l2w[vLoser];
	List<adjEntry> rightAdjOrder;
	List<adjEntry> wrongAdjOrder;
	vWinner->allAdjEntries(wrongAdjOrder);
	OGDF_ASSERT(wrongAdjOrder.size() == vWinner->degree());

	OGDF_ASSERT(vLoser->degree() <= vWinner->degree());
	// for every adjEntry of v in the "right" graph (the embedding which we want to get into the "wrong" graph)
	// we search for the corresponding adjEntry in the list of adjEntries of the "wrong" v
	for(adjEntry adj : vLoser->adjEntries){
		OGDF_ASSERT(m_mapE_l2w[adj->theEdge()] != nullptr);
		edge edgeInWinner = m_mapE_l2w[adj->theEdge()];
		adjEntry adj_in = (adj->theEdge()->adjSource() == adj ? edgeInWinner->adjSource() : edgeInWinner->adjTarget());
		rightAdjOrder.pushBack(adj_in);
	}
	List<adjEntry> adjOrder;
	vWinner->allAdjEntries(adjOrder);
	OGDF_ASSERT(vLoser->degree() <= adjOrder.size());
	if(!sameDirection){
		rightAdjOrder.reverse();
	}
	if(adjOrder.size() == rightAdjOrder.size()){
		adjOrder = rightAdjOrder;
	} else {
		List<adjEntry> helpList;
		adjOrder.split(adjOrder.get(adjOrder.size() - rightAdjOrder.size()), adjOrder, helpList);
		if(isTNodeOfPNode){
			adjOrder.concFront(rightAdjOrder);
		} else {
			adjOrder.conc(rightAdjOrder);
		}
	}
	m_gWinner->sort(vWinner, adjOrder);
}
} // end namespace ogdf
