/** \file
 * \brief Declaration of class NonPlanarCore which represents the
 *        non-planar core reduction for biconnected graphs.
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

#pragma once

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/decomposition/StaticSPQRTree.h>
#include <ogdf/basic/HashArray.h>
#include <ogdf/module/MaxFlowModule.h>


namespace ogdf {

//! Non-planar core reduction.
/**
 * @ingroup ga-planarity
 *
 * The class ogdf::NonPlanarCore implements a reduction method that reduces a graph to a
 * smaller core graph which behaves invariant with respect to non-planarity measures
 * like crossing number, skewness, coarseness, and thickness. The core reduction is
 * based on the decomposition of a graph into its triconnected components and can
 * be computed in linear time.
 *
 * The implementation is based on the following publication:
 *
 * Markus Chimani, Carsten Gutwenger: <i>Non-planar core reduction of graphs</i>.
 * Discrete Mathematics 309(7) (2009) 1838-1855
 *
 */
class OGDF_EXPORT NonPlanarCore {

	friend class GlueMap;

public:
	/**
	 * Struct to represent an edge that needs to be crossed in order to cross an st-component.
	 *
	 */
	struct CutEdge {
		const edge e; //!< the edge
		const bool dir; //!< true, iff the edge is directed from the \a s partition to the \a t partion

		CutEdge(edge e, bool dir) : e(e), dir(dir) {};
	};
	//! Algorithm call and constructor
	/**
	 * This constructor computes the non-planar core of the graph \a G.
	 *
	 * \param G the graph of which the npc is to be made
	 * \param weight if the graph is weighted, the weights otherwise a nullptr
	 * \param nonPlanarityGuaranteed iff set to true the algorithm runs a bit faster for nonplanar graphs
	 * \param maxFlowModule the MaxFlowModule that should be used for calculating the traversing path
	 * for weighted graphs. default : MaxFlowSTPlanarDigraph
	 */
	NonPlanarCore(const Graph &G, const EdgeArray<int> &weight, MaxFlowModule<int> &maxFlowModule,
	              bool nonPlanarityGuaranteed = false);

	//! An slimmed version of the Algorithm call and constructor
	NonPlanarCore(const Graph &G, const EdgeArray<int> &weight, bool nonPlanarityGuaranteed = false);

	//! An slimmed version of the Algorithm call and constructor
	NonPlanarCore(const Graph &G, bool nonPlanarityGuaranteed = false);

	//!  Returns the non-planar core
	const Graph &core() const {
		return m_graph;
	}

	//!  Returns the original graph
	const Graph &originalGraph() const {
		return *m_pOriginal;
	}

	//!  Returns the node of the original graph, which is represented by \a v in the core
	node original(node v) const {
		return m_orig[v];
	}

	//!  Returns the edges of the original graph, which are represented by \a e in the core
	List<edge> original(edge e) const {
		List<edge> res;
		if (isVirtual(e)) {
			EdgeArray<edge> origEdgesOfThisSTComponent(*mapE(e));
			for (edge eInCop: origEdgesOfThisSTComponent.graphOf()->edges) {
				if (origEdgesOfThisSTComponent[eInCop] != nullptr) {
					res.pushBack(origEdgesOfThisSTComponent[eInCop]);
				}
			}
		} else {
			res.pushBack(realEdge(e));
		}
		return res;
	}

	//!  True iff the edge \a e in the core represents more than one orginal edge and therefore is virtual
	bool isVirtual(edge e) const {
		return m_real[e] == nullptr;
	}

	//!  Returns the edge of the orginal graph, which is represented by \a e or nullptr iff \a e is virtual
	edge realEdge(edge e) const {
		return m_real[e];
	}

	/**
	 * Returns the costs of the edges in the core, which is the number of original edges crossed,
	 * if \a e is crossed, i.e. one if an edge is real and `|`mincut(edge)`|` if an edge is virtual
	 */
	const EdgeArray<int> &cost() const {
		return m_cost;
	}

	/**
	 * Returns the t node of the skeleton of the st-component represented by the core edge /a e = (s,t)
	 * Note that this node is not contained in the input graph, but an internal auxiliary graph.
	 */
	node tNode(edge e) const {
		return m_tNode[e];
	}

	/**
	 * Returns the s node of the skeleton of the st-component represented by the core edge /a e = (s,t)
	 * Note that this node is not contained in the input graph, but an internal auxiliary graph.
	 */
	node sNode(edge e) const {
		return m_sNode[e];
	}

	/**
	 * Returns a map from the edges of the st-component represented by the core edge e to the original graph
	 */
	EdgeArray<edge>* mapE(edge e) const {
		return m_mapE[e];
	}

	/**
	 * Returns the cost of \a e, which is the number of original edges crossed, if \a e is crossed,
	 * i.e. 1 if \a e is real and `|`mincut(e)`|` if e is virtual
	 */
	int cost(edge e) const {
		return m_cost[e];
	}

	//!  Returns the mincut of the st-component represented by \a e
	const List<CutEdge> &mincut(edge e) const {
		return m_mincut[e];
	}

	// destructor
	virtual ~NonPlanarCore();

	//! Inserts the crossings from a copy of the core into a copy of the original graph.
	/**
	 * This method expects \a planarCore to be planarly embedded without pseudo-crossings.
	 * \param planarCore a GraphCopy of the core, in which dummy nodes are inserted to represent crossings
	 * \param planarGraph a GraphCopy which is replaced by a GraphCopy of the original graph
	 */
	void retransform(const GraphCopy &planarCore, GraphCopy &planarGraph);

protected:
	//! The private method behind the constructors.
	void call(const Graph &G, const EdgeArray<int> *weight, MaxFlowModule<int> *maxFlowModule,
              bool nonPlanarityGuaranteed);

	/**
	 * Checks for multiedges in the core.
	 * This method is a slightly modified version of ogdf::IsParallelFreeUndirected(),
	 * that adds the functionality that the found multiedges are stored in lists.
	 *
	 * \param winningEdges The list of edges that will survive the glue.
	 * \param losingEdges The list of edges that won't survive the glue.
	 */
	void getAllMultiedges(List<edge> &winningEdges, List<edge> &losingEdges);

	/**
	 * Computes the traversing path for a given edge and the unmarked tree rooted in the node of \a eS
	 * and saves the combinatorial embedding of the st-component which \a eS represents,
	 * i.e. a list of edges that are to be crossed, when the given edge is crossed in the core.
	 * This list is minimal.
	 *
	 * \param Sv the Skeleton of one of the marked nodes of \a m_T.
	 * \param eS an edge in \a Sv.
	 * \param path a container to write the traversing path to
	 * `true` iff the source of the edge is on the s side of the cut.
	 * \param mapV a NodeArray of the original graph to map original nodes to nodes created in this method.
	 * \param coreEdge the edge in the core that represents the st-component of which the
	 * traversing path is computed.
	 * \param weight_src the weight of the edges of the original graph
	 * \param maxFlowModule same as in the constructor
	 */
	void traversingPath(const Skeleton &Sv, edge eS, List<CutEdge> &path, NodeArray<node> &mapV,
	                    edge coreEdge, const EdgeArray<int> *weight_src, MaxFlowModule<int> *maxFlowModule);

	/**
	 * To be able to insert crossings correctly, an end graph edge ought to be split into n-1 sections
	 * if n is the number of crossings on the edge.
	 * This method does exactly that.
	 *
	 * \param e The edge to be split
	 * \param splitdummies To delete the inserted dummynodes later, we store all of them in here
	 */
	void splitEdgeIntoSections(edge e, List<node> &splitdummies);

	/**
	 * After inserting the crossings, the end graph edges don't need to be partitioned anymore
	 * so the \a splitdummies get removed.
	 */
	void removeSplitdummies(List<node> &splitdummies);

	/**
	 * Every edge of \a e 's cut that doesn't go the same direction as \a e gets reversed.
	 * This method is used to both normalize the cutedges and denormalize them
	 * after the retransformation.
	 *
	 * \param coreEdge the core edge
	 */
	void normalizeCutEdgeDirection(edge coreEdge);

	/**
	 * This method asserts that all parts of the end graph that are represented by edge \a e
	 * internally have the same embedding every time retransform is called,
	 * regardless of which planarization of the core is given.
	 * Only nodes that are present in the core can have a different embedding
	 * for a different planarization of the core. They are infact reassembling the planarization of the core.
	 */
	void importEmbedding(edge e);

	/**
	 * Marks all nodes of the underlying SPQRTree and prunes planar leaves until the marked nodes span a tree,
	 * which has only non-planar leaves, i.e. non-planar R-nodes.
	 *
	 * \param mark The array where the marking is done
	 */
	void markCore(NodeArray<bool> &mark);

	/**
	 * The crossing denoted by dummy node \a v from the planarized copy of the core
	 * get inserted into the end graph.
	 */
	void inflateCrossing(node v);

	/**
	 * Get the mincut of \a e with respect to its position in the chain of its original edge.
	 *
	 * \param e The edge that we want to know the position of in the
	 * graphcopy representing the planarized version of the original graph.
	 * \param cut A list to write the mincut to.
	 */
	void getMincut(edge e, List<edge> &cut);

	/**
	 * Glues together the skeletons of \a eWinner and \a eLoser for pruned P- and S-nodes.
	 * This is done, by inserting all nodes and edges of \a eLoser's skeleton into \a eWinner's skeleton,
	 * while preserving the embedding of both skeletons.
	 *
	 * \param eWinner the core edge that will represent the newly formed skeleton/embedding
	 * \param eLoser the core edge that is copied from
	 */
	void glue(edge eWinner, edge eLoser);

	/**
	 * Glues together the mincuts of the winner and the loser edge.
	 *
	 * \param eWinner the edge whose mincut gets augmented
	 * \param eLoser the edge whose mincut gets glued to the winner mincut
	 */
	void glueMincuts(edge eWinner, edge eLoser);

	//!  The core
	Graph m_graph;

	//!  The original graph.
	const Graph *m_pOriginal;

	/**
	 * A pointer to a copy of the core, in which crossings are replaced by dummy nodes.
	 * It's a `nullptr` unless NonPlanarCore::retransform was called.
	 */
	GraphCopy const *m_planarCore;

	/**
	 * A pointer to a copy of the original graph, in which crossings are replaced by dummy nodes.
	 * It's a `nullptr` unless NonPlanarCore::retransform was called.
	 */
	GraphCopy *m_endGraph;

	//! The SPQRTree that represents the original graph.
	StaticSPQRTree m_T;

	//! Corresp. original node
	NodeArray<node> m_orig;

	//! Corresp. original edge (0 if virtual)
	EdgeArray<edge> m_real;

	//! Costs to cross each edge of the core
	EdgeArray<int> m_cost;

	//! Traversing path for an edge in the core
	EdgeArray<List<CutEdge>> m_mincut;

	//! The mapping between the nodes of each embedding and their original
	EdgeArray<NodeArray<node> *> m_mapV;

	//! The mapping between the edges of each embedding and their original
	EdgeArray<EdgeArray<edge> *> m_mapE;

	//! The graph for the underlying skeleton of a virtual edge in the core
	EdgeArray<Graph *> m_underlyingGraphs;

	//! The s node of the st-component of a core edge
	EdgeArray<node> m_sNode;

	//! The t node of the st-component of a core edge
	EdgeArray<node> m_tNode;
};

/**
 * This is a helper class to make the glueing of two edges simpler.
 */
class OGDF_EXPORT GlueMap {
public:
	/**
	 * A GlueMap is created from an NonPlanarCore and two core edges that ought to be glued together.
	 * It holds many mappings, mostly to the original graph of the core.
	 *
	 * \param eWinner This edge gets extended.
	 * \param eLoser This edge gets deleted in the end and everything it represents is transferred
	 * to \a eWinner.
	 * \param npc The NonPlanarCore all of this exists in.
	 */
	GlueMap(edge eWinner, edge eLoser, NonPlanarCore &npc);

	/**
	 * A mapping from the \a eInLoser graph to a new edge in the winner graph is created.
	 */
	void mapLoserToNewWinnerEdge(edge eInLoser);

	/**
	 * A mapping from the \a vInLoser to the \a vInWinner is created.
	 */
	void mapLoserToWinnerNode(node vInLoser, node vInWinner);

	/**
	 * A mapping from the \a vInLoser to a new node in the winner graph is created.
	 */
	void mapLoserToNewWinnerNode(node vInLoser);

	/**
	 * This method reorders the adjacency order of \a vLoser's counterpart in the winner graph
	 * according to the AdjOrder of \a vLoser in the loser graph.
	 *
	 * \param vLoser the node in question
	 * \param sameDirection false iff this method is called while handling a P Node,
	 * for which the edges are not in the same direction.
	 * \param isTNodeOfPNode true iff \a vLoser is the target node of the loser graph and the glueing
	 * process is done on a P Node.
	 */
	void reorder(node vLoser, bool sameDirection, bool isTNodeOfPNode);

	/**
	 * Getter for #m_mapV_l2w
	 * @param v the loser node
	 * @return the winner node
	 */
	node getWinnerNodeOfLoserNode(node v) const {
		return m_mapV_l2w[v];
	}

	/**
	 * Getter for #m_gLoser
	 * @return the graph that loses this glueing
	 */
	const Graph &getLoserGraph() const {
		return *m_gLoser;
	}

protected:
	//! The NonPlanarCore on which this instance operates
	NonPlanarCore &m_npc;
	//! The core edge that will survive
	const edge m_eWinner;
	//! The core edge that will be deleted
	const edge m_eLoser;
	//! The graph that eWinner represents.
	Graph *m_gWinner;
	//! The graph that eLoser represents.
	const Graph *m_gLoser;
	//! A map from the edges of the winner graph to the original graph, to denote the original of each edge.
	EdgeArray<edge> *m_mapEwinner;
	//! A map from the edges of the loser graph to the original graph, to denote the original of each node.
	const EdgeArray<edge> *m_mapEloser;
	//! A map from the nodes of the winner graph to the original graph, to denote the original of each edge.
	NodeArray<node> *m_mapVwinner;
	//! A map from the nodes of the loser graph to the original graph, to denote the original of each node.
	const NodeArray<node> *m_mapVloser;
	//! A map from the nodes of the loser graph to their new home in the winner graph.
	NodeArray<node> m_mapV_l2w;
	//! A map from the edges of the loser graph to their new home in the winner graph.
	EdgeArray<edge> m_mapE_l2w;
};
} // end namespace ogdf
