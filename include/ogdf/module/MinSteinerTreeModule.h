/** \file
 * \brief Interface of Minimum Steiner Tree Algorithms
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

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/PriorityQueue.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/internal/steinertree/Triple.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/graphalg/AStarSearch.h>
#include <sstream>

namespace ogdf {

/*!
 * \brief This class serves as an interface for various methods to compute minimum
 * or approximations of minimum Steiner trees.
 *
 * Furthermore it supplies some helping methods.
 */
template<typename T>
class MinSteinerTreeModule {
public:
	MinSteinerTreeModule() { }

	virtual ~MinSteinerTreeModule() { }

	/*!
	 * \brief Builds a minimum Steiner tree given a weighted graph and a list of terminals
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param finalSteinerTree The final Steiner tree
	 * @return The objective value (sum of edge costs) of the final Steiner tree
	 */
	virtual T call(const EdgeWeightedGraph<T> &G,
			const List<node> &terminals,
			const NodeArray<bool> &isTerminal,
			EdgeWeightedGraphCopy<T> *&finalSteinerTree
			);

	/*!
	 * \brief Prune Steiner nodes with degree 1 and their paths to terminal or branching nodes
	 * @param steinerTree The given Steiner tree
	 * @param isTerminal Incidence vector indicating terminal nodes
	 * @return The edge weights of the removed edges (achieved improvement)
	 */
	static T pruneAllDanglingSteinerPaths(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal);

	/*!
	 * \brief Prune dangling Steiner paths beginning at given nodes only
	 */
	static T pruneDanglingSteinerPathsFrom(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal, const List<node> &start);

	/*!
	 * \brief Prune the dangling Steiner path beginning at a given node only
	 */
	static T pruneDanglingSteinerPathFrom(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal, node start);

	/*!
	 * \brief Remove remaining cycles from an Steiner "almost" tree
	 * @return The edge weights of the removed edges (achieved improvement)
	 */
	static T removeCyclesFrom(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal);

	/*!
	 * \brief Undirected Single-source-shortest-paths (Dijkstra) over non-terminal nodes of G.
	 * @param G The graph
	 * @param source The start terminal
	 * @param isTerminal Incidence vector indicating terminal nodes
	 * @param distance The distance matrix result
	 * @param pred The resulting shortest-path last edge of a path
	 */
	static void singleSourceShortestPaths(const EdgeWeightedGraph<T> &G, node source, const NodeArray<bool> &isTerminal, NodeArray<T> &distance, NodeArray<edge> &pred);
	static void singleSourceShortestPathsStrict(const EdgeWeightedGraph<T> &G, node source, const NodeArray<bool> &isTerminal, NodeArray<T> &distance, NodeArray<edge> &pred);
	static void singleSourceShortestPathsDetour(const EdgeWeightedGraph<T> &G, node source, const NodeArray<bool> &isTerminal, NodeArray<T> &distance, NodeArray<edge> &pred);

	/*!
	 * \brief Undirected All-pair-shortest-paths (Floyd-Warshall) over non-terminal nodes of G.
	 * @param G The graph
	 * @param nonterminals List of non-terminals (APSP is only computed over non-terminal nodes)
	 * @param distance The distance matrix result
	 * @param pred The resulting shortest-path last edge of a path
	 */
	static void allPairShortestPaths(const EdgeWeightedGraph<T> &G, const List<node> &nonterminals, NodeArray< NodeArray<T> > &distance, NodeArray< NodeArray<edge> > &pred);
	static void allPairShortestPathsStrict(const EdgeWeightedGraph<T> &G, const NodeArray<bool> &isTerminal, NodeArray< NodeArray<T> > &distance, NodeArray< NodeArray<edge> > &pred);
	static void allPairShortestPathsDetour(const EdgeWeightedGraph<T> &G, const List<node> &nonterminals, NodeArray< NodeArray<T> > &distance, NodeArray< NodeArray<edge> > &pred);

	/*!
	 * \brief Writes an SVG file of a minimum Steiner tree in the original graph
	 * @param G The original weighted graph
	 * @param isTerminal Incidence vector indicating terminal nodes
	 * @param steinerTree The Steiner tree of the given graph
	 * @param filename The name of the output file
	 */
	static void drawSVG(const EdgeWeightedGraph<T> &G, const NodeArray<bool> &isTerminal, const EdgeWeightedGraphCopy<T> &steinerTree, const char *filename);

	/*!
	 * \brief Writes a SVG that shows only the given Steiner tree
	 * @param steinerTree The Steiner tree to be drawn
	 * @param isTerminal Incidence vector indicating terminal nodes
	 * @param filename The name of the output file
	 */
	static void drawSteinerTreeSVG(const EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal, const char *filename);

	/*!
	 * \brief Checks in O(n) time if a given tree is acually a Steiner Tree
	 * @param G The original graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param steinerTree The Steiner tree to be checked
	 * @return true iff the given Steiner tree is actually one, false otherwise
	 */
	static bool isSteinerTree(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, const EdgeWeightedGraphCopy<T> &steinerTree);

	/*!
	 * \brief Checks in O(n + m) time if a given Steiner tree problem instance is quasi-bipartite
	 * @param G The original graph
	 * @param isTerminal A bool array of terminals
	 * @return true iff the given Steiner tree problem instance is quasi-bipartite
	 */
	static bool isQuasiBipartite(const EdgeWeightedGraph<T> &G, const NodeArray<bool> &isTerminal);

protected:
	virtual T computeSteinerTree(const EdgeWeightedGraph<T> &G,
	    const List<node> &terminals,
	    const NodeArray<bool> &isTerminal,
	    EdgeWeightedGraphCopy<T> *&finalSteinerTree
	  ) = 0;
};

template<typename T>
T MinSteinerTreeModule<T>::call(const EdgeWeightedGraph<T> &G,
			const List<node> &terminals,
			const NodeArray<bool> &isTerminal,
			EdgeWeightedGraphCopy<T> *&finalSteinerTree
			)
{
	if (terminals.size() > 2) {
		return this->computeSteinerTree(G, terminals, isTerminal, finalSteinerTree);
	}

	finalSteinerTree = new EdgeWeightedGraphCopy<T>();
	finalSteinerTree->createEmpty(G);
	switch (terminals.size()) {
	case 1:
		finalSteinerTree->newNode(terminals.front());
	case 0:
		return 0;
	case 2:
		T cost(0);
		AStarSearch<T> astar;
		NodeArray<edge> pred;
		NodeArray<T> dist;
		astar.call(G, G.edgeWeights(), terminals.front(), terminals.back(), pred);
		if(pred[terminals.back()] != nullptr) {
			finalSteinerTree->newNode(terminals.back());
			for (node t = terminals.back(); t != terminals.front(); t = pred[t]->opposite(t)) {
				const edge e = pred[t];
				finalSteinerTree->newNode(e->opposite(t));
				finalSteinerTree->newEdge(e, G.weight(e));
				cost += G.weight(e);
			}
			return cost;
		}
	}
	return -1; // unreachable
}

template<typename T>
bool MinSteinerTreeModule<T>::isSteinerTree(
	const EdgeWeightedGraph<T> &G,
	const List<node> &terminals,
	const NodeArray<bool> &isTerminal,
	const EdgeWeightedGraphCopy<T> &steinerTree)
{
	// the Steiner tree is actually a tree
	if (!isTree(steinerTree)) {
		return false;
	}

	// all terminal nodes are in the graph and have degree >= 1
	for(node v : terminals) {
		const node u = steinerTree.copy(v);
		if (!u || (terminals.size() > 1 && u->degree() < 1)) {
			return false;
		}
	}

	// all Steiner nodes are inner nodes
	for(node u : steinerTree.nodes) {
		if (!isTerminal[steinerTree.original(u)]) {
			if (u->degree() <= 1) {
				return false;
			}
		}
	}

	return true;
}

template<typename T>
bool MinSteinerTreeModule<T>::isQuasiBipartite(const EdgeWeightedGraph<T> &G, const NodeArray<bool> &isTerminal)
{
	for (node v = G.firstNode(); v; v = v->succ()) {
		if (!isTerminal[v]) {
			for (adjEntry adj = v->firstAdj(); adj; adj = adj->succ()) {
				if (!isTerminal[adj->twinNode()]) {
					return false;
				}
			}
		}
	}
	return true;
}

template<typename T>
T MinSteinerTreeModule<T>::pruneDanglingSteinerPathFrom(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal, const node start)
{
	OGDF_ASSERT(isConnected(steinerTree));
	T delWeights(0);
	node u = start;
	while (u->degree() == 1
	    && !isTerminal[steinerTree.original(u)]) {
		const adjEntry adj = u->firstAdj();
		const node v = adj->twinNode();
		delWeights += steinerTree.weight(adj->theEdge());
		steinerTree.delNode(u);
		u = v;
	}
	return delWeights;
}

template<typename T>
T MinSteinerTreeModule<T>::pruneDanglingSteinerPathsFrom(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal, const List<node> &start)
{
	T delWeights(0);
	for (node v : start) {
		delWeights += pruneDanglingSteinerPathFrom(steinerTree, isTerminal, v);
	}
	return delWeights;
}

template<typename T>
T MinSteinerTreeModule<T>::pruneAllDanglingSteinerPaths(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal)
{
	List<node> start;
	for (node u : steinerTree.nodes) {
		if (u->degree() == 1
		 && !isTerminal[steinerTree.original(u)]) {
			start.pushBack(u);
		}
	}

	return pruneDanglingSteinerPathsFrom(steinerTree, isTerminal, start);
}

template<typename T>
T MinSteinerTreeModule<T>::removeCyclesFrom(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal)
{
	if (steinerTree.numberOfEdges() > steinerTree.numberOfNodes() - 1) {
		EdgeArray<bool> isInTree(steinerTree);
		T oldCost(0);
		T newCost(computeMinST(steinerTree, steinerTree.edgeWeights(), isInTree));

		List<node> pendant; // collect resulting pendant edges
		for (edge nextEdge, e = steinerTree.firstEdge(); e; e = nextEdge) {
			oldCost += steinerTree.weight(e);
			nextEdge = e->succ();
			if (!isInTree[e]) {
				if (e->source()->degree() == 2) {
					pendant.pushBack(e->source());
				}
				if (e->target()->degree() == 2) {
					pendant.pushBack(e->target());
				}
				steinerTree.delEdge(e);
			}
		}
		newCost -= MinSteinerTreeModule<T>::pruneDanglingSteinerPathsFrom(steinerTree, isTerminal, pendant);
		return oldCost - newCost;
	}
	return 0;
}

template<typename T>
void MinSteinerTreeModule<T>::singleSourceShortestPaths(const EdgeWeightedGraph<T> &G, node source, const NodeArray<bool> &isTerminal, NodeArray<T> &distance, NodeArray<edge> &pred)
{
	PrioritizedMapQueue<node, T> queue(G);
	distance.init(G, numeric_limits<T>::max());
	pred.init(G, nullptr);

	// initialization
	distance[source] = 0;
	for (node v : G.nodes) {
		queue.push(v, distance[v]);
	}

	while (!queue.empty()) {
		node v = queue.topElement();
		queue.pop();

		if (isTerminal[v] && v != source) { // v is a terminal, ignore
			continue;
		}
		if (distance[v] == numeric_limits<T>::max()) { // min is unreachable, finished
			break;
		}
		for (adjEntry adj : v->adjEntries) {
			edge e = adj->theEdge();
			node w = adj->twinNode();
			if (distance[w] > distance[v] + G.weight(e)) {
				queue.decrease(w, (distance[w] = distance[v] + G.weight(e)));
				pred[w] = e;
			}
		}
	}
}

template<typename T>
void MinSteinerTreeModule<T>::singleSourceShortestPathsStrict(const EdgeWeightedGraph<T> &G, node source, const NodeArray<bool> &isTerminal, NodeArray<T> &distance, NodeArray<edge> &pred)
{
	PrioritizedMapQueue<node, T> queue(G);
	distance.init(G, numeric_limits<T>::max());
	pred.init(G, nullptr);

	// initialization
	distance[source] = 0;
	for (node v : G.nodes) {
		queue.push(v, distance[v]);
	}

	node v = queue.topElement();
	queue.pop();

	OGDF_ASSERT(v == source);
	for (adjEntry adj : v->adjEntries) {
		edge e = adj->theEdge();
		node w = adj->twinNode();
		if (distance[w] > G.weight(e)) { // this check is only necessary for multigraphs, otherwise this is always true
			queue.decrease(w, (distance[w] = G.weight(e)));
			pred[w] = e;
		}
	}

	while (!queue.empty()) {
		v = queue.topElement();
		queue.pop();

		if (distance[v] == numeric_limits<T>::max()) { // min is unreachable, finished
			break;
		}
		for (adjEntry adj : v->adjEntries) {
			edge e = adj->theEdge();
			node w = adj->twinNode();
			T dist = distance[v] + G.weight(e);
			if (distance[w] > dist) {
				queue.decrease(w, (distance[w] = dist));
				if (isTerminal[v] || pred[v] == nullptr) { // w is on a path with terminal
					pred[w] = nullptr;
				} else {
					pred[w] = e;
				}
			} else
			if (distance[w] == dist
			 && pred[w] != nullptr) { // tie
				if (isTerminal[v] || pred[v] == nullptr) { // w is on a path with terminal
					pred[w] = nullptr;
				} else {
					pred[w] = e;
				}
			}
		}
	}
}

template<typename T>
void MinSteinerTreeModule<T>::singleSourceShortestPathsDetour(const EdgeWeightedGraph<T> &G, node source, const NodeArray<bool> &isTerminal, NodeArray<T> &distance, NodeArray<edge> &pred)
{
	PrioritizedMapQueue<node, T> queue(G);
	distance.init(G, numeric_limits<T>::max());
	pred.init(G, nullptr);

	// initialization
	distance[source] = 0;
	for (node v : G.nodes) {
		queue.push(v, distance[v]);
	}

	while (!queue.empty()) {
		node v = queue.topElement();
		queue.pop();

		if (isTerminal[v] && v != source) { // v is a terminal, ignore
			continue;
		}
		if (distance[v] == numeric_limits<T>::max()) { // min is unreachable, finished
			break;
		}
		for (adjEntry adj : v->adjEntries) {
			edge e = adj->theEdge();
			node w = adj->twinNode();
			if (distance[w] > distance[v] + G.weight(e)) {
				queue.decrease(w, (distance[w] = distance[v] + G.weight(e)));
				pred[w] = e;
			}
		}
	}
}

template<typename T>
void MinSteinerTreeModule<T>::allPairShortestPaths(const EdgeWeightedGraph<T> &G, const List<node> &nonterminals, NodeArray< NodeArray<T> > &distance, NodeArray< NodeArray<edge> > &pred)
{
	// initialization
	for (node u = G.firstNode(); u; u = u->succ()) {
		distance[u].init(G, numeric_limits<T>::max());
		pred[u].init(G, nullptr);
	}
	for (edge e = G.firstEdge(); e; e = e->succ()) {
		const node u = e->source(), v = e->target();
		distance[u][v] = G.weight(e);
		distance[v][u] = G.weight(e);
		pred[u][v] = e;
		pred[v][u] = e;
	}

	// main loop
	for(node v : nonterminals) {
		for (node u : G.nodes) {
			const T duv = distance[u][v];
			if (duv < numeric_limits<T>::max()) {
				for (node w = u->succ(); w; w = w->succ()) {
					const T dvw = distance[v][w];
					if (dvw < numeric_limits<T>::max()
					 && duv + dvw < distance[u][w]) {
						distance[w][u] = distance[u][w] = duv + dvw;
						pred[u][w] = pred[v][w];
						pred[w][u] = pred[v][u];
					}
				}
			}
		}
	}
}

template<typename T>
void MinSteinerTreeModule<T>::allPairShortestPathsStrict(const EdgeWeightedGraph<T> &G, const NodeArray<bool> &isTerminal, NodeArray< NodeArray<T> > &distance, NodeArray< NodeArray<edge> > &pred)
{
	// initialization
	for (node u : G.nodes) {
		distance[u].init(G, numeric_limits<T>::max());
		pred[u].init(G, nullptr);
	}
	for (edge e : G.edges) {
		const node u = e->source(), v = e->target();
		distance[u][v] = distance[v][u] = G.weight(e);
		pred[u][v] = pred[v][u] = e;
	}

	// main loop
	for (node v : G.nodes) {
		if (isTerminal[v]) { // v is a terminal
			for (node u : G.nodes) {
				const T duv = distance[u][v];
				if (duv < numeric_limits<T>::max()) {
					for (node w = u->succ(); w; w = w->succ()) {
						const T dvw = distance[v][w];
						const T duvw = duv + dvw;
						if (dvw < numeric_limits<T>::max()
						 && duvw <= distance[u][w]) { // prefer terminals
							distance[w][u] = distance[u][w] = duvw;
							pred[w][u] = pred[u][w] = nullptr;
						}
					}
				}
			}
		} else { // v is not a terminal
			for (node u : G.nodes) {
				const T duv = distance[u][v];
				if (duv < numeric_limits<T>::max()) {
					for (node w = u->succ(); w; w = w->succ()) {
						const T dvw = distance[v][w];
						const T duvw = duv + dvw;
						if (dvw < numeric_limits<T>::max()
						 && duvw < distance[u][w]) { // do not prefer nonterminals
							distance[w][u] = distance[u][w] = duvw;
							pred[u][w] = (pred[u][v] ? pred[v][w] : nullptr);
							pred[w][u] = (pred[w][v] ? pred[v][u] : nullptr);
						}
					}
				}
			}
		}
	}
}

template<typename T>
void MinSteinerTreeModule<T>::allPairShortestPathsDetour(const EdgeWeightedGraph<T> &G, const List<node> &nonterminals, NodeArray< NodeArray<T> > &distance, NodeArray< NodeArray<edge> > &pred)
{
	// initialization
	for (node u = G.firstNode(); u; u = u->succ()) {
		distance[u].init(G, numeric_limits<T>::max());
		pred[u].init(G, nullptr);
	}
	for (edge e = G.firstEdge(); e; e = e->succ()) {
		const node u = e->source(), v = e->target();
		distance[u][v] = G.weight(e);
		distance[v][u] = G.weight(e);
		pred[u][v] = e;
		pred[v][u] = e;
	}

	// main loop
	for (node v : nonterminals) {
		for (node u : G.nodes) {
			const T duv = distance[u][v];
			if (duv < numeric_limits<T>::max()) {
				for (node w = u->succ(); w; w = w->succ()) {
					const T dvw = distance[v][w];
					if (dvw < numeric_limits<T>::max()
					 && duv + dvw < distance[u][w]) {
						distance[w][u] = distance[u][w] = duv + dvw;
						pred[u][w] = pred[v][w];
						pred[w][u] = pred[v][u];
					}
				}
			}
		}
	}
}

template<typename T>
void MinSteinerTreeModule<T>::drawSteinerTreeSVG(const EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal, const char *filename)
{
	GraphAttributes GA(steinerTree,
	  GraphAttributes::nodeGraphics |
	  GraphAttributes::nodeStyle |
	  GraphAttributes::nodeLabel |
	  GraphAttributes::edgeGraphics |
	  GraphAttributes::edgeStyle |
	  GraphAttributes::edgeLabel);

	GA.setDirected(false);

	string s;

	for (node v : steinerTree.nodes) {
		std::stringstream out;
		GA.width(v) = GA.height(v) = 25.0;
		if (isTerminal[steinerTree.original(v)]) {
			out << "T";
			GA.shape(v) = shRect;
			GA.fillColor(v) = Color::Red;
		} else {
			out << "S";
			GA.shape(v) = shEllipse;
			GA.fillColor(v) = Color::Gray;
		}
		out << steinerTree.original(v);
		GA.label(v) = out.str();
	}

	FMMMLayout fmmm;

	fmmm.useHighLevelOptions(true);
	fmmm.unitEdgeLength(44.0);
	fmmm.newInitialPlacement(true);
	fmmm.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);

	fmmm.call(GA);
	ofstream writeStream(filename, std::ofstream::out);
	GraphIO::drawSVG(GA, writeStream);
}

template<typename T>
void MinSteinerTreeModule<T>::drawSVG(const EdgeWeightedGraph<T> &G, const NodeArray<bool> &isTerminal, const EdgeWeightedGraphCopy<T> &steinerTree, const char *filename)
{
	GraphAttributes GA(G,
	  GraphAttributes::nodeGraphics |
	  GraphAttributes::nodeStyle |
	  GraphAttributes::nodeLabel |
	  GraphAttributes::edgeGraphics |
	  GraphAttributes::edgeStyle |
	  GraphAttributes::edgeLabel);

	GA.setDirected(false);

	for (edge e : G.edges) {
		GA.strokeColor(e) = Color::Black;
		GA.label(e) = to_string(G.weight(e));
		GA.strokeWidth(e) = 1;
	}
	for (edge e : steinerTree.edges) {
		GA.strokeColor(steinerTree.original(e)) = Color::Red;
		GA.strokeWidth(steinerTree.original(e)) = 2;
	}

	for (node v : G.nodes) {
		std::stringstream out;
		GA.width(v) = GA.height(v) = 25.0;
		GA.strokeColor(v) = Color::Black;
		if (isTerminal[v]) {
			out << "T" << v;
			GA.shape(v) = shRect;
			GA.fillColor(v) = Color::Red;
			GA.strokeWidth(v) = 2;
		} else {
			out << "S" << v;
			GA.shape(v) = shEllipse;
			if (steinerTree.copy(v)) {
				GA.fillColor(v) = Color::Gray;
				GA.strokeWidth(v) = 2;
			} else {
				GA.fillColor(v) = Color::White;
				GA.strokeWidth(v) = 1;
			}
		}
		GA.label(v) = out.str();
	}

	FMMMLayout fmmm;

	fmmm.useHighLevelOptions(true);
	fmmm.unitEdgeLength(44.0);
	fmmm.newInitialPlacement(true);
	fmmm.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);

	fmmm.call(GA);

	ofstream writeStream(filename, std::ofstream::out);
	GraphIO::drawSVG(GA, writeStream);
}

} // end namespace ogdf
