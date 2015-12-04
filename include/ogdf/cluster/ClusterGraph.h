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

#pragma once

#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/GraphObserver.h>


namespace ogdf {

class OGDF_EXPORT ClusterGraph;
class OGDF_EXPORT ClusterGraphObserver;
class OGDF_EXPORT ClusterElement;

typedef ClusterElement *cluster; //!< The type of clusters.


//! Representation of clusters in a clustered graph.
/**
 * \see ClusterGraph
 */
class OGDF_EXPORT ClusterElement : private internal::GraphElement {

	friend class ClusterGraph;
	friend class internal::GraphList<ClusterElement>;

	int	m_id;         //!< The index of this cluster.
	int	m_depth;      //!< The depth of this cluster in the cluster tree.

public:

	/**
	 * @name Graph object containers
	 * These containers maintain the nodes and child clusters of the cluster, and provide iterators.
	 * If computed they also provide access to the sorted list of adjacency entries of edges leaving the cluster.
	 */
	//@{

	//! The container containing the nodes lying (directly) in this cluster.
	ListContainer<node, ClusterElement> nodes;

	//! The container containing the child clusters (children in the cluster tree) of this cluster.
	ListContainer<cluster, ClusterElement> children;

	//! The container containing the sorted list of adjacency entries of edges leaving this cluster.
	/**
	 * This list is only available (i.e., non-empty) if explicitly computed (e.g., by a cluster-planar embedding algorithm).
	 */
	ListContainer<adjEntry, ClusterElement>	adjEntries;

	//@}

private:

	cluster			m_parent;		//!< The parent of this cluster.
	cluster			m_pPrev;		//!< The postorder predecessor of this cluster.
	cluster			m_pNext;		//!< The postorder successor of this cluster.
	ListIterator<cluster> m_it;		//!< The position of this cluster within the children list of its parent.

#ifdef OGDF_DEBUG
	// we store the graph containing this cluster for debugging purposes only
	const ClusterGraph *m_pClusterGraph;
#endif


	//! Provides access to the encapsulated list of children (for friends of ClusterElement).
	List<cluster> &getChildren() {
		return children;
	}

	//! Provides access to the encapsulated list of nodes (for friends of ClusterElement).
	List<node> &getNodes() {
		return nodes;
	}

	//! Provides access to the encapsulated list of adjacency entries (for friends of ClusterElement).
	List<adjEntry> &getAdjEntries() {
		return adjEntries;
	}

	//! Traverses the inclusion tree and adds nodes to \a clusterNodes.
	/**
	 * Invoked by public function getClusterNodes(List<node> &clusterNodes).
	 */
	void getClusterInducedNodes(List<node> &clusterNodes);

	void getClusterInducedNodes(NodeArray<bool> &clusterNode, int& num);

public:

	//! Creates a new cluster element.
#ifdef OGDF_DEBUG
	ClusterElement(const ClusterGraph *pClusterGraph, int id) :
		m_id(id), m_depth(0), m_parent(nullptr), m_pPrev(nullptr), m_pNext(nullptr),
		m_pClusterGraph(pClusterGraph) { }
#else
	ClusterElement(int id) :
		m_id(id), m_depth(0), m_parent(nullptr), m_pPrev(nullptr), m_pNext(nullptr) { }
#endif


	/**
	* @name Access methods
	*/
	//@{

#ifdef OGDF_DEBUG
	const ClusterGraph *graphOf() const { return m_pClusterGraph; }
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
	void getClusterNodes(List<node> &clusterNodes) {
		clusterNodes.clear();
		getClusterInducedNodes(clusterNodes);
	}

	//! Sets the entry for each node v to true if v is a member of the subgraph induced by the ClusterElement.
	/**
	 * All other entries remain unchanged!
	 * \pre \a clusterNode is a NodeArray initialized on the clustergraph this ClusterElement belongs to.
	 * \return the number of entries set to true.
	 */
	int getClusterNodes(NodeArray<bool> &clusterNode) {
		int num = 0;
		getClusterInducedNodes(clusterNode, num);
		return num;
	}

	//@{
	/**
	* @name Iteration over tree structure
	* Alternatively you can use the containers ClusterElement::nodes, ClusterElement::children and ClusterElement::adjEntries directly.
	*/
	//@{

	//! Returns the first element in the list of child clusters.
	ListConstIterator<ClusterElement*> cBegin() const { return children.begin(); }

	//! Returns the last element in the list of child clusters.
	ListConstIterator<ClusterElement*> crBegin() const { return children.rbegin(); }

	//! Returns the number of child clusters.
	int cCount() { return children.size(); }

	//! Returns the first element in list of child nodes.
	ListConstIterator<node> nBegin() const{ return nodes.begin(); }

	//! Returns the number of child nodes.
	int nCount() { return nodes.size(); }


	//! Returns the first adjacency entry in the list of outgoing edges.
	ListConstIterator<adjEntry> firstAdj() const { return adjEntries.begin(); }

	//! Returns the last adjacency entry in the list of outgoing edges.
	ListConstIterator<adjEntry> lastAdj () const { return adjEntries.rbegin(); }

	//@}

	//! Standard Comparer (uses cluster indices).
	static int compare(const ClusterElement& x,const ClusterElement& y) { return x.m_id - y.m_id; }
	OGDF_AUGMENT_COMPARER(ClusterElement)

	OGDF_NEW_DELETE

};// class ClusterElement



//---------------------------------------------------------
// iteration macros
//---------------------------------------------------------

//! Iterates over all outgoing edges (given by the outgoing adjacency entries).
//! @ingroup graphs
#define forall_cluster_adj(adj,c)\
for(ogdf::ListConstIterator<adjEntry> ogdf_loop_var=(c)->firstAdj();\
	ogdf::test_forall_adj_entries_of_cluster(ogdf_loop_var,(adj));\
	ogdf_loop_var=ogdf_loop_var.succ())

//! Iterates over all outgoing edges (given by the outgoing adjacency entries).
//! @ingroup graphs
#define forall_cluster_rev_adj(adj,c)\
for(ogdf::ListConstIterator<adjEntry> ogdf_loop_var=(c)->lastAdj();\
	ogdf::test_forall_adj_entries_of_cluster(ogdf_loop_var,(adj));\
	ogdf_loop_var=ogdf_loop_var.pred())

//! Iterates over all outgoing edges.
//! @ingroup graphs
#define forall_cluster_adj_edges(e,c)\
for(ogdf::ListConstIterator<adjEntry> ogdf_loop_var=(c)->firstAdj();\
	ogdf::test_forall_adj_edges_of_cluster(ogdf_loop_var,(e));\
	ogdf_loop_var=ogdf_loop_var.succ())



inline bool test_forall_adj_entries_of_cluster(ListConstIterator<adjEntry> &it, adjEntry &adj)
{
	if (it.valid()) { adj = (*it); return true; }
	else return false;
}

inline bool test_forall_adj_edges_of_cluster(ListConstIterator<adjEntry> &it, edge &e)
{
	adjEntry adj = (*it);
	if (adj) { e = adj->theEdge(); return true; }
	else return false;
}

inline bool test_forall_adj_edges_of_cluster(adjEntry &adj, edge &e)
{
	if (adj) { e = adj->theEdge(); return true; }
	else return false;
}


//! Iteration over all clusters \a c of cluster graph \a C.
//! @ingroup graphs
#define forall_clusters(c,C) for((c)=(C).firstCluster(); (c); (c)=(c)->succ())

//! Iteration over all clusters \a c of cluster graph \a C (in postorder).
//! @ingroup graphs
#define forall_postOrderClusters(c,C)\
for((c)=(C).firstPostOrderCluster(); (c); (c)=(c)->pSucc())


class ClusterArrayBase;
template<class T>class ClusterArray;



//! Representation of clustered graphs.
/**
 * @ingroup graphs
 *
 * This class is derived from GraphObserver and handles hierarchical
 * clustering of the nodes in a graph, providing additional functionality.
 */
class OGDF_EXPORT ClusterGraph : public GraphObserver
{
	const Graph *m_pGraph;			  //!< The associated graph.

	int		m_clusterIdCount;		  //!< The index assigned to the next created cluster.
	int		m_clusterArrayTableSize;  //!< The current table size of cluster arrays.

	mutable cluster m_postOrderStart; //!< The first cluster in postorder.
	cluster	m_rootCluster;			  //!< The root cluster.

	bool    m_adjAvailable;		  //! True if the adjacency list for each cluster is available.
	bool    m_allowEmptyClusters; //! Defines if empty clusters are deleted immediately if generated by operations.

	NodeArray<cluster> m_nodeMap; //!< Stores the cluster of each node.
	//! Stores for every node its position within the children list of its cluster.
	NodeArray<ListIterator<node> >  m_itMap;

	mutable ListPure<ClusterArrayBase*> m_regClusterArrays; //!< The registered cluster arrays.
	mutable ListPure<ClusterGraphObserver*> m_regObservers; //!< The registered graph observers.

#ifndef OGDF_MEMORY_POOL_NTS
	mutable std::mutex m_mutexRegArrays; //!< The critical section for protecting shared acces to register/unregister methods.
#endif

public:

	/**
	* @name Iterators
	* These types are used for graph object iterators, which are returned by graph object containers
	* like nodes and edges.
	*/
	//@{

	//! Provides a bidirectional iterator to a cluster in a clustered graph.
	typedef internal::GraphIterator<cluster> cluster_iterator;

	//@}

	/**
	* @name Graph object containers
	* These containers maintain the nodes and edges of the graph, and provide node and edge iterators.
	*/
	//@{

	//! The container containing all cluster objects.
	internal::GraphObjectContainer<ClusterElement> clusters;

	//@}


	//! Creates a cluster graph associated with no graph.
	/**
	 * After using this constructor, you will have to initialize it with a graph before you can use it, see init().
	 */
	ClusterGraph();

	//! Creates a cluster graph associated with graph \a G.
	/**
	 * All nodes in \a G are assigned to the root cluster.
	 */
	ClusterGraph(const Graph &G);

	//! Copy constructor.
	/**
	 * The copy constructor only creates a copy of the cluster tree structure, not of the underlying graph.
	 * Consequently, this cluster graph and \a C will be associated with the same graph instance.
	 */
	ClusterGraph(const ClusterGraph &C);

	//! Copies the underlying graph of \a C into \a G and constructs a copy of \a C associated with \a G.
	/**
	 * The created cluster tree is a copy of \a C except for the associated nodes, which are the newly created copies in \a G.
	 */
	ClusterGraph(const ClusterGraph &C, Graph &G);

	//! Copies the underlying graph of \a C into \a G and constructs a copy of \a C associated with \a G.
	/**
	 * The created cluster tree is a copy of \a C except for the associated nodes, which are the newly created copies in \a G.
	 * Stores the new copies of the original nodes and clusters in the arrays \a originalNodeTable and \a originalClusterTable.
	 */
	ClusterGraph(
		const ClusterGraph &C,
		Graph &G,
		ClusterArray<cluster> &originalClusterTable,
		NodeArray<node> &originalNodeTable);

	//! Copies the underlying graph of \a C into \a G and constructs a copy of \a C associated with \a G.
	/**
	 * The created cluster tree is a copy of \a C except for the associated nodes, which are the newly created copies in \a G.
	 * Stores the new copies of the original nodes, edges, and clusters in the arrays
	 * \a originalNodeTable, \a edgeCopy, and \a originalClusterTable.
	 */
	ClusterGraph(
		const ClusterGraph &C,
		Graph &G,
		ClusterArray<cluster> &originalClusterTable,
		NodeArray<node> &originalNodeTable,
		EdgeArray<edge> &edgeCopy);

	//! Destructor.
	virtual ~ClusterGraph();


	/**
	* @name Access methods
	*/
	//@{

	//! Returns the root cluster.
	cluster rootCluster() const { return m_rootCluster; }

	//! Returns the number of clusters.
	int numberOfClusters() const { return clusters.size(); }

	//! Returns the maximal used cluster index.
	int maxClusterIndex() const { return m_clusterIdCount-1; }

	//! Returns the table size of cluster arrays associated with this graph.
	int clusterArrayTableSize() const { return m_clusterArrayTableSize; }

	//! Returns the cluster to which a node belongs.
	inline cluster clusterOf(node v) const {
		return m_nodeMap[v];
	}

	//! Returns depth of cluster c in cluster tree, starting with root depth 1.
	//should be called instead of direct c->depth call to assure
	//valid depth
	int& clusterDepth(cluster c) const
	{
		// updateDepth must be set to true if depth info shall be used
		OGDF_ASSERT(m_updateDepth);

		//initialize depth at first call
		if (!m_depthUpToDate)
			computeSubTreeDepth(rootCluster());
		OGDF_ASSERT(c->depth() != 0)
		return c->m_depth;
	}

	//! Returns the first cluster in the list of all clusters.
	cluster firstCluster() const { return clusters.head(); }

	//! Returns the last cluster in the list of all cluster.
	cluster lastCluster() const { return clusters.tail(); }

	//! Returns the first cluster in the list of post ordered clusters.
	cluster firstPostOrderCluster() const {
		if (!m_postOrderStart) postOrder();
		return m_postOrderStart;
	}

	//! Returns the list of all clusters in \a clusters.
	template<class CLUSTERLIST>
	void allClusters(CLUSTERLIST &clusterList) const {
		clusterList.clear();
		for (cluster c : clusters)
			clusterList.pushBack(c);
	}

	//@}
	/**
	* @name Modification methods
	*/
	//@{

	//! Removes all clusters except for the root cluster.
	void clear();

	//! Clears all cluster data and then reinitializes the instance with underlying graph \a G.
	void init(const Graph &G);

	//! Removes all clusters from the cluster subtree rooted at cluster C except for cluster C itself.
	void clearClusterTree(cluster C);

	//! Inserts a new cluster; makes it a child of the cluster \a parent.
	cluster newCluster(cluster parent, int id = -1);

	//! Creates an empty cluster with index \a clusterId and parent \a parent.
	cluster createEmptyCluster(const cluster parent = nullptr, int clusterId = -1);

	//! Creates a new cluster containing the nodes given by \a nodes; makes it a child of the cluster \a parent.
	/**
	 * The nodes are reassigned to the new cluster. If you turn off
	 * \a m_allowEmptyclusters, an emptied cluster is deleted except if all
	 * nodes are put into the same cluster.
	 * @param nodes are the nodes that will be reassigned to the new cluster.
	 * @param parent is the parent of the new cluster.
	 * \return the created cluster.
	 */
	cluster createCluster(SList<node>& nodes, const cluster parent = nullptr);

	//! Deletes cluster \a c.
	/**
	 * All subclusters become children of parent cluster of \a c.
	 * \pre \a c is not the root cluster.
	 */
	void delCluster(cluster c);

	//! Moves cluster \a c to a new parent \a newParent.
	void moveCluster(cluster c, cluster newParent);


	//! Reassigns node \a v to cluster \ c.
	void reassignNode(node v, cluster c);

	//! Clear cluster info structure, reinitializes with underlying graph \a G.
	//inserted mainly for use in gmlparser.
	void reInit(Graph& G) {
		reinitGraph(G);
	}

	//! Collapses all nodes in the list \a nodes to the first node; multi-edges are removed.
	template<class NODELIST>
	void collaps(NODELIST &nodes, Graph &G) {
		OGDF_ASSERT(&G == m_pGraph);
		m_adjAvailable = false;

		m_postOrderStart = nullptr;
		node v = nodes.popFrontRet();
		while (!nodes.empty())
		{
			node w = nodes.popFrontRet();
			adjEntry adj = w->firstAdj();
			while (adj != nullptr)
			{
				adjEntry succ = adj->succ();
				edge e = adj->theEdge();
				if (e->source() == v || e->target() == v)
					G.delEdge(e);
				else if (e->source() == w)
					G.moveSource(e, v);
				else
					G.moveTarget(e, v);
				adj = succ;
			}
			//because nodes can already be unassigned (they are always
			//unassigned if deleted), we have to check this
			/*
			if (m_nodeMap[w])
			{
			cluster c = m_nodeMap[w];
			c->m_entries.del(m_itMap[w]);
			}
			*/
			//removeNodeAssignment(w);
			G.delNode(w);
		}
	}

	//@}


	/**
	* @name Cluster tree queries
	*/
	//@{

	//! Turns automatic update of node depth values on or off.
	void setUpdateDepth(bool b) const {
		m_updateDepth = b;
		//make sure that depth cant be used anymore
		//(even if it may still be valid a little while)
		if (!b) m_depthUpToDate = false;
	}

	//! Updates depth information in subtree after delCluster.
	void pullUpSubTree(cluster c);

	//! Computes depth of cluster tree, running time O(C).
	//maybe later we should provide a permanent depth member update
	int treeDepth() const;

	//! Computes depth of cluster tree hanging at \a c.
	void computeSubTreeDepth(cluster c) const;

	//! Returns lowest common cluster of nodes in list \a nodes.
	cluster commonCluster(SList<node>& nodes);

	//! Returns the lowest common cluster of \a v and \a w in the cluster tree
	/**
	 * \pre \a v and \a w are nodes in the graph.
	 */
	cluster commonCluster(node v, node w) const {
		cluster c1, c2;
		return commonClusterLastAncestors(v, w, c1, c2);
	}

	//! Returns the lowest common cluster lca and the highest ancestors on the path to lca.
	cluster commonClusterLastAncestors(
		node v,
		node w,
		cluster& c1,
		cluster& c2) const
	{
		List<cluster> eL;
		return commonClusterAncestorsPath(v, w, c1, c2, eL);
	}

	//! Returns lca of \a v and \a w and stores corresponding path in \a eL.
	/**
	 * The list \a eL is directed from \a v to \a w.
	 */
	cluster commonClusterPath(
		node v,
		node w,
		List<cluster>& eL) const
	{
		cluster c1, c2;
		return commonClusterAncestorsPath(v, w, c1, c2, eL);
	}

	//! Returns lca of \a v and \a w, stores corresponding path in \a eL and ancestors in \a c1, \a c2.
	cluster commonClusterAncestorsPath(
		node v,
		node w,
		cluster& c1,
		cluster& c2,
		List<cluster>& eL) const;

	//! Returns the list of clusters that are empty or only contain empty clusters.
	/**
	 * The list is constructed in an order that allows deletion and reinsertion.
	 * We never set rootcluster to be one of the empty clusters!!
	 * if checkClusters is given, only list elements are checked
	 * to allow efficient checking in the case
	 * that you know which clusters were recently changed (e.g. node reass.)
	 */
	void emptyClusters(SList<cluster>& emptyCluster, SList<cluster>* checkCluster = 0);

	//! Returns true if cluster \a c has only one node and no children.
	inline bool emptyOnNodeDelete(cluster c) //virtual?
	{
		//if (!c) return false; //Allows easy use in loops
		return (c->nCount() == 1) && (c->cCount() == 0);
	}

	//! Returns true if cluster \a c has only one child and no nodes.
	inline bool emptyOnClusterDelete(cluster c) //virtual?
	{
		//if (!c) return false; //Allows easy use in loops
		return (c->nCount() == 0) && (c->cCount() == 1);
	}

	//@}
	/**
	* @name Adjacent edges
	*/
	//@{

	//! Returns the list of all edges adjacent to cluster \a c in \a edges.
	template<class EDGELIST>
	void adjEdges(cluster c, EDGELIST &edges) const {
		edges.clear();

		if (m_adjAvailable)
		{
			edge e;
			forall_cluster_adj_edges(e,c)
				edges.pushBack(e);
		}
	}

	//! Returns the list of all adjacency entries adjacent to cluster \a c in \a entries.
	template<class ADJLIST>
	void adjEntries(cluster c, ADJLIST &entries) const {
		entries.clear();
		if (m_adjAvailable)
		{
			for(adjEntry adj : c->adjEntries)
				entries.pushBack(adj);
		}
	}

	//! Computes the adjacency entry list for cluster \a c.
	template<class LISTITERATOR>
	void makeAdjEntries(cluster c,LISTITERATOR start) {

		c->adjEntries.clear();
		LISTITERATOR its;
		for (its = start; its.valid(); its++)
		{
			adjEntry adj = (*its);
			c->adjEntries.pushBack(adj);
		}
	}

	//! Sets the availability status of the adjacency entries.
	void adjAvailable(bool val) { m_adjAvailable = val; }


	//@}
	/**
	* @name Miscellaneous
	*/
	//@{

	//! Checks the combinatorial cluster planar embedding.
	bool representsCombEmbedding() const;

	//! Checks the consistency of the data structure.
	/**
	* \remark This method is meant for debugging purposes only.
	*
	* @return true if everything is ok, false if the data structure is inconsistent.
	*/
	bool consistencyCheck() const;

	//@}
	/**
	* @name Registering arrays and observers
	* These methods are used by ClusterArray or ClusterGraphObserver.
	* There should be no need to use them directly in user code.
	*/
	//@{

	//! Registers a cluster array.
	ListIterator<ClusterArrayBase*> registerArray(ClusterArrayBase *pClusterArray) const;

	//! Unregisters a cluster array.
	void unregisterArray(ListIterator<ClusterArrayBase*> it) const;

	//! Move the registration \a it of a cluster array to \a pClusterArray (used with move semantics for cluster arrays).
	void moveRegisterArray(ListIterator<ClusterArrayBase*> it, ClusterArrayBase *pClusterArray) const;

	//! Registers a cluster graph observer.
	ListIterator<ClusterGraphObserver*> registerObserver(ClusterGraphObserver *pObserver) const;

	//! Unregisters a cluster graph observer.
	void unregisterObserver(ListIterator<ClusterGraphObserver*> it) const;

	//@}
	/**
	* @name Operators and conversion
	*/
	//@{

	//! Conversion to const Graph reference (to underlying graph).
	operator const Graph &() const { return *m_pGraph; }

	//! Returns a reference to the underlying graph.
	const Graph &constGraph() const { return *m_pGraph; }

	//! Assignment operator.
	ClusterGraph &operator=(const ClusterGraph &C);

	//@}

protected:
	mutable ClusterArray<int>* m_lcaSearch; //!< Used to save last search run number for commoncluster.
	mutable int m_lcaNumber;//!< Used to save last search run number for commoncluster.
	mutable ClusterArray<cluster>* m_vAncestor;//!< Used to save last search run number for commoncluster.
	mutable ClusterArray<cluster>* m_wAncestor;//!< Used to save last search run number for commoncluster.

	mutable bool m_updateDepth; //!< Depth of clusters is always updated if set to true.
	mutable bool m_depthUpToDate; //!< Status of cluster depth information.

	//! Creates new cluster containing nodes in parameter list
	//! with index \a clusterid.
	cluster doCreateCluster(SList<node>& nodes,
		const cluster parent, int clusterId = -1);

	//! Creates new cluster containing nodes in parameter list and
	//! stores resulting empty clusters in list, cluster has index \a clusterid.
	cluster doCreateCluster(SList<node>& nodes,
		SList<cluster>& emptyCluster,
		const cluster parent, int clusterId = -1);

	//! Clears all cluster data.
	void doClear();

	//! Copies lowest common ancestor info to copy of clustered graph.
	void copyLCA(const ClusterGraph &C, ClusterArray<cluster>* clusterCopy = nullptr);
	//int m_treeDepth; //should be implemented and updated in operations?

	//! Adjusts the post order structure for moved clusters.
	//we assume that c is inserted via pushback in newparent->children
	void updatePostOrder(cluster c, cluster oldParent, cluster newParent);

	//! Computes new predecessor for subtree at moved cluster \a c (nullptr if \a c is the root).
	cluster postOrderPredecessor(cluster c) const;

	//! Leftmost cluster in subtree rooted at c, gets predecessor of subtree.
	cluster leftMostCluster(cluster c) const;

	//---------------------------------------
	//functions inherited from GraphObserver:
	//define how to cope with graph changes

	//! Implementation of inherited method: Updates data if node deleted.
	virtual void nodeDeleted(node v);

	//! Implementation of inherited method: Updates data if node added.
	virtual void nodeAdded(node v) {
		assignNode(v, rootCluster());
	}

	//! Implementation of inherited method: Updates data if edge deleted.
	virtual void edgeDeleted(edge /* e */) { }

	//! Implementation of inherited method: Updates data if edge added.
	virtual void edgeAdded(edge /* e */) { }

	//! Currently does nothing.
	virtual void reInit() { }

	//! Clears cluster data without deleting root when underlying graphs' clear method is called.
	virtual void cleared() {
		//we don't want a complete clear, as the graph still exists
		//and can be updated from input stream
		clear();
	}//Graph cleared

private:
	//! Assigns node \a v to cluster \a c (\a v not yet assigned!).
	void assignNode(node v, cluster C);

	//! Unassigns node \a v from its cluster.
	void unassignNode(node v);

	//! Remove the assignment entries for nodes.
	//! Checks if node is currently not assigned.
	void removeNodeAssignment(node v) {
		if (m_nodeMap[v]) //iff == 0, itmap == 0 !!?
		{
			cluster c2 = m_nodeMap[v];
			c2->nodes.del(m_itMap[v]);
			m_nodeMap[v] = nullptr;
			m_itMap[v] = ListIterator<node>();
		}
	}

	//! Performs a copy of the cluster structure of C,
	//! the underlying graph stays the same.
	void shallowCopy(const ClusterGraph &C);

	//! Perform a deep copy on C, C's underlying
	//! graph is copied into G.
	void deepCopy(const ClusterGraph &C,Graph &G);

	//! Perform a deep copy on C, C's underlying
	//! graph is copied into G. Stores associated nodes in \a originalNodeTable.

	void deepCopy(
		const ClusterGraph &C,Graph &G,
		ClusterArray<cluster> &originalClusterTable,
		NodeArray<node> &originalNodeTable);

	//! Perform a deep copy on C, C's underlying
	//! graph is copied into G.  Stores associated nodes in \a originalNodeTable
	//! and edges in \a edgeCopy.
	void deepCopy(
		const ClusterGraph &C,Graph &G,
		ClusterArray<cluster> &originalClusterTable,
		NodeArray<node> &originalNodeTable,
		EdgeArray<edge> &edgeCopy);


	void clearClusterTree(cluster c,List<node> &attached);

	void initGraph(const Graph &G);

	//! Reinitializes instance with graph \a G.
	void reinitGraph(const Graph &G);

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

	void postOrder(cluster c,SListPure<cluster> &S) const;

	void reinitArrays();
};



ostream &operator<<(ostream &os, ogdf::cluster c);


} // end namespace ogdf
