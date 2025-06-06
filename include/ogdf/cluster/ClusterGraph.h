/** \file
 * \brief Derived class of GraphObserver providing additional functionality
 * to handle clustered graphs.
 *
 * \author Sebastian Leipert, Karsten Klein
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
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Observer.h>
#include <ogdf/basic/RegisteredArray.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/comparer.h>
#include <ogdf/basic/memory.h>

#include <functional>
#include <iosfwd>
#include <memory>
#include <utility>

namespace ogdf {

class ClusterElement; // IWYU pragma: keep
class ClusterGraph; // IWYU pragma: keep
class ClusterGraphObserver; // IWYU pragma: keep

using cluster = ClusterElement*; //!< The type of clusters.

//! Representation of clusters in a clustered graph.
/**
 * \see ClusterGraph
 */
class OGDF_EXPORT ClusterElement : private internal::GraphElement {
	friend class ClusterGraph;
	friend class internal::GraphList<ClusterElement>;

	int m_id; //!< The index of this cluster.
	int m_depth; //!< The depth of this cluster in the cluster tree.

public:
	/**
	 * @name Graph object containers
	 * These containers maintain the nodes and child clusters of the cluster, and provide iterators.
	 * If computed they also provide access to the sorted list of adjacency entries of edges leaving the cluster.
	 */
	//! @{

	//! The container containing the nodes lying (directly) in this cluster.
	ListContainer<node, ClusterElement> nodes;

	//! The container containing the child clusters (children in the cluster tree) of this cluster.
	ListContainer<cluster, ClusterElement> children;

	//! The container containing the sorted list of adjacency entries of edges leaving this cluster.
	/**
	 * This list is only available (i.e., non-empty) if explicitly computed (e.g., by a cluster-planar embedding algorithm).
	 */
	ListContainer<adjEntry, ClusterElement> adjEntries;

	//! @}

private:
	cluster m_parent; //!< The parent of this cluster.
	cluster m_pPrev; //!< The postorder predecessor of this cluster.
	cluster m_pNext; //!< The postorder successor of this cluster.
	ListIterator<cluster> m_it; //!< The position of this cluster within the children list of its parent.

#ifdef OGDF_DEBUG
	// we store the graph containing this cluster for debugging purposes only
	const ClusterGraph* m_pClusterGraph;
#endif


	//! Provides access to the encapsulated list of children (for friends of ClusterElement).
	List<cluster>& getChildren() { return children; }

	//! Provides access to the encapsulated list of nodes (for friends of ClusterElement).
	List<node>& getNodes() { return nodes; }

	//! Provides access to the encapsulated list of adjacency entries (for friends of ClusterElement).
	List<adjEntry>& getAdjEntries() { return adjEntries; }

	//! Traverses the inclusion tree and adds nodes to \p clusterNodes.
	/**
	 * Invoked by public function getClusterNodes(List<node> &clusterNodes).
	 */
	void getClusterInducedNodes(List<node>& clusterNodes);

	void getClusterInducedNodes(NodeArray<bool>& clusterNode, int& num);

public:
	//! Creates a new cluster element.
#ifdef OGDF_DEBUG
	ClusterElement(const ClusterGraph* pClusterGraph, int id)
		: m_id(id)
		, m_depth(0)
		, m_parent(nullptr)
		, m_pPrev(nullptr)
		, m_pNext(nullptr)
		, m_pClusterGraph(pClusterGraph) { }
#else
	explicit ClusterElement(int id)
		: m_id(id), m_depth(0), m_parent(nullptr), m_pPrev(nullptr), m_pNext(nullptr) { }
#endif


	/**
	 * @name Access methods
	 */
	//! @{

#ifdef OGDF_DEBUG
	const ClusterGraph* graphOf() const { return m_pClusterGraph; }
#endif


	//! Returns the (unique) index of the cluster.
	int index() const { return m_id; }

	//! Returns the depth of the cluster in the cluster tree.
	int depth() const { return m_depth; }

	//! Returns the parent of the cluster.
	cluster parent() { return m_parent; }

	//! Returns the successor of the cluster in the list of all clusters.
	cluster succ() const { return static_cast<cluster>(m_next); }

	//! Returns the predecessor of the cluster in the list of all clusters.
	cluster pred() const { return static_cast<cluster>(m_prev); }

	//! Returns the postorder successor of the cluster in the list of all clusters.
	cluster pSucc() const { return m_pNext; }

	//! Returns the postorder predecessor of the cluster in the list of all clusters.
	cluster pPred() const { return m_pPrev; }

	//! Returns the list of nodes in the cluster, i.e., all nodes in the subtree rooted at this cluster.
	/**
	 * Recursively traverses the cluster tree starting at this cluster.
	 */
	void getClusterNodes(List<node>& clusterNodes) {
		clusterNodes.clear();
		getClusterInducedNodes(clusterNodes);
	}

	//! Sets the entry for each node v to true if v is a member of the subgraph induced by the ClusterElement.
	/**
	 * All other entries remain unchanged!
	 * \pre \p clusterNode is a NodeArray initialized on the clustergraph this ClusterElement belongs to.
	 * \return the number of entries set to true.
	 */
	int getClusterNodes(NodeArray<bool>& clusterNode) {
		int num = 0;
		getClusterInducedNodes(clusterNode, num);
		return num;
	}

	//! Checks whether a cluster \p child is a descendant (i.e. child, child of a child, ...) of this cluster.
	/**
	 * @param child the cluster that might be a descendant of this
	 * @param allow_equal when not given or false, return false in the case that this == child
	 * \return whether child is a descendant of this
	 */
	bool isDescendant(ClusterElement* child, bool allow_equal = false) {
		OGDF_ASSERT(child != nullptr);
		if (!allow_equal) {
			child = child->parent();
		}
		while (child != this) {
			if (child == nullptr) {
				return false;
			}
			child = child->parent();
		}
		return true;
	}

	//! @}
	/**
	 * @name Iteration over tree structure
	 * Alternatively you can use the containers ClusterElement::nodes, ClusterElement::children and ClusterElement::adjEntries directly.
	 */
	//! @{

	//! Returns the first element in the list of child clusters.
	ListConstIterator<ClusterElement*> cBegin() const { return children.begin(); }

	//! Returns the last element in the list of child clusters.
	ListConstIterator<ClusterElement*> crBegin() const { return children.rbegin(); }

	//! Returns the number of child clusters.
	int cCount() { return children.size(); }

	//! Returns the first element in list of child nodes.
	ListConstIterator<node> nBegin() const { return nodes.begin(); }

	//! Returns the number of child nodes.
	int nCount() { return nodes.size(); }

	//! Returns the first adjacency entry in the list of outgoing edges.
	ListConstIterator<adjEntry> firstAdj() const { return adjEntries.begin(); }

	//! Returns the last adjacency entry in the list of outgoing edges.
	ListConstIterator<adjEntry> lastAdj() const { return adjEntries.rbegin(); }

	//! @}

	//! Standard Comparer (uses cluster indices).
	static int compare(const ClusterElement& x, const ClusterElement& y) { return x.m_id - y.m_id; }
	OGDF_AUGMENT_COMPARER(ClusterElement)

	OGDF_NEW_DELETE
};

//! \name Iteration macros
//! @{

//! Iterates over all outgoing edges (given by the outgoing adjacency entries).
//! @ingroup graphs
#define forall_cluster_adj(adj, c)                                          \
	for (ogdf::ListConstIterator<adjEntry> ogdf_loop_var = (c)->firstAdj(); \
			ogdf::test_forall_adj_entries_of_cluster(ogdf_loop_var, (adj)); \
			ogdf_loop_var = ogdf_loop_var.succ())

//! Iterates over all outgoing edges (given by the outgoing adjacency entries).
//! @ingroup graphs
#define forall_cluster_rev_adj(adj, c)                                      \
	for (ogdf::ListConstIterator<adjEntry> ogdf_loop_var = (c)->lastAdj();  \
			ogdf::test_forall_adj_entries_of_cluster(ogdf_loop_var, (adj)); \
			ogdf_loop_var = ogdf_loop_var.pred())

//! Iterates over all outgoing edges.
//! @ingroup graphs
#define forall_cluster_adj_edges(e, c)                                      \
	for (ogdf::ListConstIterator<adjEntry> ogdf_loop_var = (c)->firstAdj(); \
			ogdf::test_forall_adj_edges_of_cluster(ogdf_loop_var, (e));     \
			ogdf_loop_var = ogdf_loop_var.succ())

inline bool test_forall_adj_entries_of_cluster(ListConstIterator<adjEntry>& it, adjEntry& adj) {
	if (it.valid()) {
		adj = (*it);
		return true;
	} else {
		return false;
	}
}

inline bool test_forall_adj_edges_of_cluster(ListConstIterator<adjEntry>& it, edge& e) {
	adjEntry adj = (*it);
	if (adj) {
		e = adj->theEdge();
		return true;
	} else {
		return false;
	}
}

inline bool test_forall_adj_edges_of_cluster(adjEntry& adj, edge& e) {
	if (adj) {
		e = adj->theEdge();
		return true;
	} else {
		return false;
	}
}

//! Iteration over all clusters \p c of cluster graph \p C.
//! @ingroup graphs
#define forall_clusters(c, C) for ((c) = (C).firstCluster(); (c); (c) = (c)->succ())

//! Iteration over all clusters \p c of cluster graph \p C (in postorder).
//! @ingroup graphs
#define forall_postOrderClusters(c, C) \
	for ((c) = (C).firstPostOrderCluster(); (c); (c) = (c)->pSucc())

//! @}

using ClusterGraphRegistry = RegistryBase<cluster, ClusterGraph, internal::GraphIterator<cluster>>;

template<typename Value, bool WithDefault>
class ClusterArrayBase; // IWYU pragma: keep

#define OGDF_DECL_REG_ARRAY_TYPE(v, c) ClusterArrayBase<v, c>
OGDF_DECL_REG_ARRAY(ClusterArray)
#undef OGDF_DECL_REG_ARRAY_TYPE

/**
 * Abstract base class for cluster graph observers.
 *
 * @ingroup graphs
 *
 * If a class needs to keep track of changes in a clustered graph like addition or deletion
 * of clusters, you can derive it from ClusterGraphObserver and override the
 * notification methods clusterDeleted, clusterAdded.
 */
class OGDF_EXPORT ClusterGraphObserver : public Observer<ClusterGraph, ClusterGraphObserver> {
public:
	ClusterGraphObserver() = default;

	OGDF_DEPRECATED("calls registrationChanged with only partially-constructed child classes, "
					"see copy constructor of Observer for fix")

	explicit ClusterGraphObserver(const ClusterGraph* CG) { reregister(CG); }

	virtual void clusterDeleted(cluster v) = 0;
	virtual void clusterAdded(cluster v) = 0;
	virtual void clustersCleared() = 0;

	const ClusterGraph* getGraph() const { return getObserved(); }
};

//! Representation of clustered graphs.
/**
 * @ingroup graphs
 *
 * This class is derived from GraphObserver and handles hierarchical
 * clustering of the nodes in a graph, providing additional functionality.
 */
class OGDF_EXPORT ClusterGraph
	: private GraphObserver, // to get updates when graph nodes/edges added/removed
	  public Observable<ClusterGraphObserver, ClusterGraph>, // to publish updates when clusters added/removed
	  public ClusterGraphRegistry // for ClusterArrays
{
	using Obs = Observable<ClusterGraphObserver, ClusterGraph>;

	int m_clusterIdCount = 0; //!< The index assigned to the next created cluster.

	mutable cluster m_postOrderStart = nullptr; //!< The first cluster in postorder.
	cluster m_rootCluster = nullptr; //!< The root cluster.

	bool m_adjAvailable = false; //! True if the adjacency list for each cluster is available.
	bool m_allowEmptyClusters =
			true; //! Defines if empty clusters are deleted immediately if generated by operations.

	NodeArray<cluster> m_nodeMap; //!< Stores the cluster of each node.
	//! Stores for every node its position within the children list of its cluster.
	NodeArray<ListIterator<node>> m_itMap;

public:
	/**
	 * @name Iterators
	 * These types are used for graph object iterators, which are returned by graph object containers
	 * like nodes and edges.
	 */
	//! @{

	//! Provides a bidirectional iterator to a cluster in a clustered graph.
	using cluster_iterator = internal::GraphIterator<cluster>;

	//! @}

	/**
	 * @name Graph object containers
	 * These containers maintain the nodes and edges of the graph, and provide node and edge iterators.
	 */
	//! @{

	//! The container containing all cluster objects.
	internal::GraphObjectContainer<ClusterElement> clusters;

	//! @}


	//! Creates a cluster graph associated with no graph.
	/**
	 * After using this constructor, you will have to initialize it with a graph before you can use it, see init().
	 */
	ClusterGraph();

	//! Creates a cluster graph associated with graph \p G.
	/**
	 * All nodes in \p G are assigned to the root cluster.
	 */
	explicit ClusterGraph(const Graph& G);

	//! Copy constructor.
	/**
	 * The copy constructor only creates a copy of the cluster tree structure, not of the underlying graph.
	 * Consequently, this cluster graph and \p C will be associated with the same graph instance.
	 */
	ClusterGraph(const ClusterGraph& C);

	//! Copies the underlying graph of \p C into \p G and constructs a copy of \p C associated with \p G.
	/**
	 * The created cluster tree is a copy of \p C except for the associated nodes, which are the newly created copies in \p G.
	 */
	ClusterGraph(const ClusterGraph& C, Graph& G);

	//! Copies the underlying graph of \p C into \p G and constructs a copy of \p C associated with \p G.
	/**
	 * The created cluster tree is a copy of \p C except for the associated nodes, which are the newly created copies in \p G.
	 * Stores the new copies of the original nodes and clusters in the arrays \p originalNodeTable and \p originalClusterTable.
	 */
	ClusterGraph(const ClusterGraph& C, Graph& G, ClusterArray<cluster>& originalClusterTable,
			NodeArray<node>& originalNodeTable);

	//! Copies the underlying graph of \p C into \p G and constructs a copy of \p C associated with \p G.
	/**
	 * The created cluster tree is a copy of \p C except for the associated nodes, which are the newly created copies in \p G.
	 * Stores the new copies of the original nodes, edges, and clusters in the arrays
	 * \p originalNodeTable, \p edgeCopy, and \p originalClusterTable.
	 */
	ClusterGraph(const ClusterGraph& C, Graph& G, ClusterArray<cluster>& originalClusterTable,
			NodeArray<node>& originalNodeTable, EdgeArray<edge>& edgeCopy);

	//! Destructor.
	virtual ~ClusterGraph();

	/**
	 * @name Access methods
	 */
	//! @{

	using GraphObserver::getGraph;

	//! Returns the root cluster.
	cluster rootCluster() const { return m_rootCluster; }

	//! Returns the number of clusters.
	int numberOfClusters() const { return clusters.size(); }

	//! Returns the maximal used cluster index.
	int maxClusterIndex() const { return m_clusterIdCount - 1; }

	//! Returns the cluster to which a node belongs.
	inline cluster clusterOf(node v) const { return m_nodeMap[v]; }

	//! Returns depth of cluster c in cluster tree, starting with root depth 1.
	//should be called instead of direct c->depth call to assure
	//valid depth
	int& clusterDepth(cluster c) const {
		// updateDepth must be set to true if depth info shall be used
		OGDF_ASSERT(m_updateDepth);

		//initialize depth at first call
		if (!m_depthUpToDate) {
			computeSubTreeDepth(rootCluster());
		}
		OGDF_ASSERT(c->depth() != 0);
		return c->m_depth;
	}

	//! Returns the first cluster in the list of all clusters.
	cluster firstCluster() const { return clusters.head(); }

	//! Returns the last cluster in the list of all cluster.
	cluster lastCluster() const { return clusters.tail(); }

	//! Returns the first cluster in the list of post ordered clusters.
	cluster firstPostOrderCluster() const {
		if (!m_postOrderStart) {
			postOrder();
		}
		return m_postOrderStart;
	}

	//! Returns the list of all clusters in \p clusterList.
	template<class CLUSTERLIST>
	void allClusters(CLUSTERLIST& clusterList) const {
		clusterList.clear();
		for (cluster c : clusters) {
			clusterList.pushBack(c);
		}
	}

	//! @}
	/**
	 * @name Modification methods
	 */
	//! @{

	//! Removes all clusters except for the root cluster.
	void clear();

	//! Clears all cluster data and then reinitializes the instance with underlying graph \p G.
	void init(const Graph& G);

	//! Removes all clusters from the cluster subtree rooted at cluster C except for cluster C itself.
	void clearClusterTree(cluster C);

	//! Inserts a new cluster; makes it a child of the cluster \p parent.
	cluster newCluster(cluster parent, int id = -1);

	//! Creates an empty cluster with index \p clusterId and parent \p parent.
	cluster createEmptyCluster(const cluster parent = nullptr, int clusterId = -1);

	//! Creates a new cluster containing the nodes given by \p nodes; makes it a child of the cluster \p parent.
	/**
	 * The nodes are reassigned to the new cluster. If you turn off
	 * #m_allowEmptyClusters, an emptied cluster is deleted except if all
	 * nodes are put into the same cluster.
	 * @param nodes are the nodes that will be reassigned to the new cluster.
	 * @param parent is the parent of the new cluster.
	 * \return the created cluster.
	 */
	cluster createCluster(const SList<node>& nodes, const cluster parent = nullptr);

	//! Deletes cluster \p c.
	/**
	 * All subclusters become children of parent cluster of \p c.
	 * \pre \p c is not the root cluster.
	 */
	void delCluster(cluster c);

	//! Moves cluster \p c to a new parent \p newParent.
	void moveCluster(cluster c, cluster newParent);


	//! Reassigns node \p v to cluster \p c.
	void reassignNode(node v, cluster c);

	//! Clear cluster info structure, reinitializes with underlying graph \p G.
	//inserted mainly for use in gmlparser.
	void reInit(Graph& G) { reinitGraph(G); }

	//! Constructs a cluster tree
	void copyClusterTree(
			const ClusterGraph& C, const Graph& G, ClusterArray<cluster>& originalClusterTable,
			std::function<node(node)> nodeMap = [](node v) { return v; });

	//! Collapses all nodes in the list \p nodes to the first node; multi-edges are removed.
	template<class NODELIST>
	void collapse(NODELIST& nodes, Graph& G) {
		OGDF_ASSERT(&G == &constGraph());
		m_adjAvailable = false;

		m_postOrderStart = nullptr;
		node v = nodes.popFrontRet();
		while (!nodes.empty()) {
			node w = nodes.popFrontRet();
			adjEntry adj = w->firstAdj();
			while (adj != nullptr) {
				adjEntry succ = adj->succ();
				edge e = adj->theEdge();
				if (e->source() == v || e->target() == v) {
					G.delEdge(e);
				} else if (e->source() == w) {
					G.moveSource(e, v);
				} else {
					G.moveTarget(e, v);
				}
				adj = succ;
			}
			//because nodes can already be unassigned (they are always
			//unassigned if deleted), we have to check this
#if 0
			if (m_nodeMap[w])
			{
				cluster c = m_nodeMap[w];
				c->m_entries.del(m_itMap[w]);
			}
			removeNodeAssignment(w);
#endif
			G.delNode(w);
		}
	}

	//! @}


	/**
	 * @name Cluster tree queries
	 */
	//! @{

	/**
	 * Returns a random cluster.
	 *
	 * \c nullptr is returned if no feasible cluster exists.
	 *
	 * @see chooseIteratorFrom
	 */
	cluster chooseCluster(
			std::function<bool(cluster)> includeCluster = [](cluster) { return true; },
			bool isFastTest = true) const;

	//! Turns automatic update of node depth values on or off.
	void setUpdateDepth(bool b) const {
		m_updateDepth = b;
		//make sure that depth cant be used anymore
		//(even if it may still be valid a little while)
		if (!b) {
			m_depthUpToDate = false;
		}
	}

	//! Updates depth information in subtree after delCluster.
	void pullUpSubTree(cluster c);

	//! Computes depth of cluster tree, running time O(C).
	//maybe later we should provide a permanent depth member update
	int treeDepth() const;

	//! Computes depth of cluster tree hanging at \p c.
	void computeSubTreeDepth(cluster c) const;

	//! Returns lowest common cluster of nodes in list \p nodes.
	cluster commonCluster(SList<node>& nodes);

	//! Returns the lowest common cluster of \p v and \p w in the cluster tree
	/**
	 * \pre \p v and \p w are nodes in the graph.
	 */
	cluster commonCluster(node v, node w) const {
		cluster c1, c2;
		return commonClusterLastAncestors(v, w, c1, c2);
	}

	//! Returns the lowest common cluster lca and the highest ancestors on the path to lca.
	cluster commonClusterLastAncestors(node v, node w, cluster& c1, cluster& c2) const {
		List<cluster> eL;
		return commonClusterAncestorsPath(v, w, c1, c2, eL);
	}

	//! Returns lca of \p v and \p w and stores corresponding path in \p eL.
	/**
	 * The list \p eL is directed from \p v to \p w.
	 */
	cluster commonClusterPath(node v, node w, List<cluster>& eL) const {
		cluster c1, c2;
		return commonClusterAncestorsPath(v, w, c1, c2, eL);
	}

	//! Returns lca of \p v and \p w, stores corresponding path in \p eL and ancestors in \p c1, \p c2.
	cluster commonClusterAncestorsPath(node v, node w, cluster& c1, cluster& c2,
			List<cluster>& eL) const;

	//! Returns the list of clusters that are empty or only contain empty clusters.
	/**
	 * The list is constructed in an order that allows deletion and reinsertion.
	 * We never set rootcluster to be one of the empty clusters!!
	 * if checkClusters is given, only list elements are checked
	 * to allow efficient checking in the case
	 * that you know which clusters were recently changed (e.g. node reass.)
	 */
	void emptyClusters(SList<cluster>& emptyCluster, SList<cluster>* checkCluster = nullptr);

	//! Returns true if cluster \p c has only one node and no children.
	inline bool emptyOnNodeDelete(cluster c) //virtual?
	{
#if 0
		if (!c) {
			return false; //Allows easy use in loops
		}
#endif
		return c->nCount() == 1 && c->cCount() == 0;
	}

	//! Returns true if cluster \p c has only one child and no nodes.
	inline bool emptyOnClusterDelete(cluster c) //virtual?
	{
#if 0
		if (!c) {
			return false; //Allows easy use in loops
		}
#endif
		return c->nCount() == 0 && c->cCount() == 1;
	}

	//! @}
	/**
	 * @name Adjacent edges
	 */
	//! @{

	//! Returns the list of all edges adjacent to cluster \p c in \p edges.
	template<class EDGELIST>
	void adjEdges(cluster c, EDGELIST& edges) const {
		edges.clear();

		if (m_adjAvailable) {
			edge e;
			forall_cluster_adj_edges(e, c) {
				edges.pushBack(e);
			}
		}
	}

	//! Returns the list of all adjacency entries adjacent to cluster \p c in \p entries.
	template<class ADJLIST>
	void adjEntries(cluster c, ADJLIST& entries) const {
		entries.clear();
		if (m_adjAvailable) {
			for (adjEntry adj : c->adjEntries) {
				entries.pushBack(adj);
			}
		}
	}

	//! Computes the adjacency entry list for cluster \p c.
	template<class LISTITERATOR>
	void makeAdjEntries(cluster c, LISTITERATOR start) {
		c->adjEntries.clear();
		LISTITERATOR its;
		for (its = start; its.valid(); its++) {
			adjEntry adj = (*its);
			c->adjEntries.pushBack(adj);
		}
	}

	//! Gets the availability status of the adjacency entries.
	bool adjAvailable() const { return m_adjAvailable; }

	//! Sets the availability status of the adjacency entries.
	void adjAvailable(bool val) { m_adjAvailable = val; }

	//! @}
	/**
	 * @name Miscellaneous
	 */
	//! @{

	//! Checks the combinatorial cluster planar embedding.
	/**
	 * @return true if the current embedding (given by the adjacency lists of the clusters)
	 *         together with the combinatorial embedding of the underlying constGraph()
	 *         represents a cluster planar combinatorial embedding
	 */
	bool representsCombEmbedding() const;

	//! Checks the combinatorial cluster planar embedding.
	/**
	 * This only works when the underlying Graph is connected and represents a planar embedding,
	 * so check constGraph().representsCombEmbedding() first. Otherwise,
	 * use the slower representsCombEmbedding().
	 *
	 * @return true if the current embedding (given by the adjacency lists of the clusters)
	 *         represents a cluster planar combinatorial embedding
	 */
	bool representsConnectedCombEmbedding() const;

#ifdef OGDF_DEBUG
	//! Asserts consistency of this cluster graph.
	void consistencyCheck() const;
#endif

	//! @}
	/**
	 * @name Registering arrays and observers
	 * These methods are used by ClusterArray or ClusterGraphObserver.
	 * There should be no need to use them directly in user code.
	 */
	//! @{

	static inline int keyToIndex(cluster key) { return key->index(); }

	bool isKeyAssociated(cluster key) const {
		if (key == nullptr) {
			return false;
		}
#ifdef OGDF_DEBUG
		if (key->graphOf() == this) {
			OGDF_ASSERT(keyToIndex(key) < this->getArraySize());
			return true;
		} else {
			return false;
		}
#else
		return true;
#endif
	}

	int calculateArraySize(int add) const { return calculateTableSize(m_clusterIdCount + add); }

	int maxKeyIndex() const { return (m_clusterIdCount)-1; }

	cluster_iterator begin() const { return clusters.begin(); }

	cluster_iterator end() const { return clusters.end(); }

	//! @}
	/**
	 * @name Operators and conversion
	 */
	//! @{

	//! Conversion to const Graph reference (to underlying graph).
	operator const Graph&() const { return *getGraph(); }

	//! Returns a reference to the underlying graph.
	const Graph& constGraph() const { return *getGraph(); }

	//! Assignment operator.
	ClusterGraph& operator=(const ClusterGraph& C);

	//! @}

protected:
	mutable std::unique_ptr<ClusterArray<int>> m_lcaSearch; //!< Used to save last search run number for commoncluster.
	mutable int m_lcaNumber = 0; //!< Used to save last search run number for commoncluster.
	mutable std::unique_ptr<ClusterArray<cluster>> m_vAncestor; //!< Used to save last search run number for commoncluster.
	mutable std::unique_ptr<ClusterArray<cluster>> m_wAncestor; //!< Used to save last search run number for commoncluster.

	mutable bool m_updateDepth = false; //!< Depth of clusters is always updated if set to true.
	mutable bool m_depthUpToDate = false; //!< Status of cluster depth information.

	//! Creates new cluster containing nodes in parameter list
	//! with index \p clusterId.
	cluster doCreateCluster(const SList<node>& nodes, const cluster parent, int clusterId = -1);

	//! Creates new cluster containing nodes in parameter list and
	//! stores resulting empty clusters in list, cluster has index \p clusterId.
	cluster doCreateCluster(const SList<node>& nodes, SList<cluster>& emptyCluster,
			const cluster parent, int clusterId = -1);

	//! Clears all cluster data.
	void doClear();

	//! Copies lowest common ancestor info to copy of clustered graph.
	void copyLCA(const ClusterGraph& C);
	//int m_treeDepth; //should be implemented and updated in operations?

	//! Adjusts the post order structure for moved clusters.
	//we assume that c is inserted via pushback in newparent->children
	void updatePostOrder(cluster c, cluster oldParent, cluster newParent);

	//! Computes new predecessor for subtree at moved cluster \p c (nullptr if \p c is the root).
	cluster postOrderPredecessor(cluster c) const;

	//! Leftmost cluster in subtree rooted at c, gets predecessor of subtree.
	cluster leftMostCluster(cluster c) const;

	//! \name Functions inherited from GraphObserver (define how to cope with graph changes)
	//! @{

	//! Implementation of inherited method: Updates data if node deleted.
	void nodeDeleted(node v) override;

	//! Implementation of inherited method: Updates data if node added.
	void nodeAdded(node v) override { assignNode(v, rootCluster()); }

	//! Implementation of inherited method: Updates data if edge deleted.
	void edgeDeleted(edge /* e */) override { }

	//! Implementation of inherited method: Updates data if edge added.
	void edgeAdded(edge /* e */) override { }

	//! Clears cluster data without deleting root when underlying graphs' clear method is called.
	void cleared() override {
		//we don't want a complete clear, as the graph still exists
		//and can be updated from input stream
		clear();
	}

	void registrationChanged(const Graph* newG) override;

	// need to re-export both to allow parameter-type-based overload in the same namespace
	friend Observer<ClusterGraph, ClusterGraphObserver>;
	friend Observer<ClusterGraph, RegisteredObserver<ClusterGraph>>;
	using ClusterGraphRegistry::registerObserver;
	using ClusterGraphRegistry::unregisterObserver;
	using Obs::registerObserver;
	using Obs::unregisterObserver;

	//! @}

private:
	//! Fills \p emptyCluster with empty, non-root clusters from \p clusterList
	template<typename T>
	inline void fillEmptyClusters(SList<cluster>& emptyCluster, const T& clusterList) const {
		emptyCluster.clear();

		for (cluster cc : clusterList) {
			if (cc->cCount() + cc->nCount() == 0 && cc != rootCluster()) { // we dont add rootcluster
				emptyCluster.pushBack(cc);
			}
		}
	}

	//! Clears children of a cluster, used by recursive #clearClusterTree
	void recurseClearClusterTreeOnChildren(cluster c, List<node>& attached) {
		m_adjAvailable = false;
		for (auto child : c->getChildren()) {
			clearClusterTree(child, attached);
		}
	}

	//! Assigns node \p v to cluster \p C (\p v not yet assigned!).
	void assignNode(node v, cluster C);

	//! Unassigns node \p v from its cluster.
	void unassignNode(node v);

	//! Remove the assignment entries for nodes.
	//! Checks if node \p v is currently not assigned.
	void removeNodeAssignment(node v) {
		if (m_nodeMap[v]) //iff == 0, itmap == 0 !!?
		{
			cluster c2 = m_nodeMap[v];
			c2->nodes.del(m_itMap[v]);
			m_nodeMap[v] = nullptr;
			m_itMap[v] = ListIterator<node>();
		}
	}

	//! Performs a copy of the cluster structure of \p C,
	//! the underlying graph stays the same.
	void shallowCopy(const ClusterGraph& C);

	//! Perform a deep copy on \p C, \p C's underlying
	//! graph is copied into \p G.
	void deepCopy(const ClusterGraph& C, Graph& G);

	//! Perform a deep copy on \p C, \p C's underlying
	//! graph is copied into \p G. Stores associated nodes in \p originalNodeTable.
	void deepCopy(const ClusterGraph& C, Graph& G, ClusterArray<cluster>& originalClusterTable,
			NodeArray<node>& originalNodeTable);

	//! Perform a deep copy on \p C, \p C's underlying
	//! graph is copied into \p G. Stores associated nodes in \p originalNodeTable
	//! and edges in \p edgeCopy.
	void deepCopy(const ClusterGraph& C, Graph& G, ClusterArray<cluster>& originalClusterTable,
			NodeArray<node>& originalNodeTable, EdgeArray<edge>& edgeCopy);


	void clearClusterTree(cluster c, List<node>& attached);

	void initGraph(const Graph& G);

	//! Reinitializes instance with graph \p G.
	void reinitGraph(const Graph& G);

	//! Creates new cluster with given id, precondition: id not used
	cluster newCluster(int id);

	//! Creates new cluster.
	cluster newCluster();

	//! Create postorder information in cluster tree.
	void postOrder() const;

#ifdef OGDF_DEBUG
	//! Check postorder information in cluster tree.
	void checkPostOrder() const;
#endif

	void postOrder(cluster c, SListPure<cluster>& S) const;
};

//! RegisteredArray for labeling the \ref cluster "clusters" of a ClusterGraph.
template<class Value, bool WithDefault>
class ClusterArrayBase : public RegisteredArray<ClusterGraph, Value, WithDefault> {
	using RA = RegisteredArray<ClusterGraph, Value, WithDefault>;

public:
	using RA::RA;

	//! Default Constructor.
	ClusterArrayBase() = default;

	//! Creates a new cluster array with initial capacity \p size.
	ClusterArrayBase(const ClusterGraph& C, const Value& def, int size) : RA(C, def) {
		RA::resize(size, true);
	};

	//! Returns a pointer to the associated cluster graph.
	ClusterGraph* graphOf() const { return RA::registeredAt(); }
};

OGDF_EXPORT std::ostream& operator<<(std::ostream& os, cluster c);

//! Turn cluster borders into cycles of edges and cluster-border-edge-crossings into vertices.
/**
 * This subdivides each edge once for each cluster border it crosses and then cyclically connects
 * the subdivision vertices according to their order in adjEntries(), which usually corresponds
 * to a cluster-planar embedding (in which case the resulting graph will also be planarly embedded).
 * Graph /p G may be a copy of the constGraph() of this ClusterGraph, in which case the \p translate
 * function should translate nodes of constGraph() into nodes of \p G.
 * @param CG A cluster-planar embedded ClusterGraph.
 * @param G (A copy of) The graph underlying \p CG, into which we insert the new edges. Must have a corresponding planar embedding.
 * @param subdivisions If non-null, will be assigned information about the subdivisions that were created.
 * @param translate A mapping from the nodes of constGraph() to those of \p G.
 * @pre \p CG.adjAvailable() is true
 */
OGDF_EXPORT void planarizeClusterBorderCrossings(const ClusterGraph& CG, Graph& G,
		EdgeArray<List<std::pair<adjEntry, cluster>>>* subdivisions,
		const std::function<edge(edge)>& translate);

}
