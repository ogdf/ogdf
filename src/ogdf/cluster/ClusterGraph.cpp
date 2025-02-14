/** \file
 * \brief Implements the class ClusterGraph, providing
 * extra functionality for clustered graphs.
 * A clustered graph C=(G,T) consists of an undirected graph G
 * and a rooted tree T in which the leaves of T correspond
 * to the vertices of G=(V,E).
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/Reverse.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/exceptions.h>
#include <ogdf/cluster/ClusterGraph.h>

#include <functional>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>

#include "ogdf/basic/GraphCopy.h"

namespace ogdf {

using Math::nextPower2;

void ClusterElement::getClusterInducedNodes(List<node>& clusterNodes) {
	for (node v : nodes) {
		clusterNodes.pushBack(v);
	}
	for (cluster c : children) {
		c->getClusterInducedNodes(clusterNodes);
	}
}

void ClusterElement::getClusterInducedNodes(NodeArray<bool>& clusterNode, int& num) {
	for (node v : nodes) {
		clusterNode[v] = true;
	}
	num += nodes.size();

	for (cluster c : children) {
		c->getClusterInducedNodes(clusterNode, num);
	}
}

ClusterGraph::ClusterGraph() { resizeArrays(); }

// Construction of a new cluster graph. All nodes
// are children of the root cluster
ClusterGraph::ClusterGraph(const Graph& G) { initGraph(G); }

ClusterGraph::ClusterGraph(const ClusterGraph& C)
	// need to explicitly call default parent class constructors in our copy constructor
	: GraphObserver(), Obs(), ClusterGraphRegistry() {
	shallowCopy(C);
}

ClusterGraph::ClusterGraph(const ClusterGraph& C, Graph& G,
		ClusterArray<cluster>& originalClusterTable, NodeArray<node>& originalNodeTable) {
	deepCopy(C, G, originalClusterTable, originalNodeTable);
}

ClusterGraph::ClusterGraph(const ClusterGraph& C, Graph& G,
		ClusterArray<cluster>& originalClusterTable, NodeArray<node>& originalNodeTable,
		EdgeArray<edge>& edgeCopy) {
	deepCopy(C, G, originalClusterTable, originalNodeTable, edgeCopy);
}

ClusterGraph::ClusterGraph(const ClusterGraph& C, Graph& G) { deepCopy(C, G); }

ClusterGraph::~ClusterGraph() {
	Obs::clearObservers();
	// this is only necessary because GraphObjectContainer simply deallocs its memory without calling destructors
	while (!clusters.empty()) {
		clusters.del(clusters.head());
	}
}

// Construction of a new cluster graph. All nodes
// are children of the root cluster
void ClusterGraph::init(const Graph& G) {
	doClear();
	m_clusterIdCount = 0;
	m_postOrderStart = nullptr;

	m_lcaNumber = 0;
	initGraph(G);
}

ClusterGraph& ClusterGraph::operator=(const ClusterGraph& C) {
	doClear();
	shallowCopy(C);

#ifdef OGDF_HEAVY_DEBUG
	consistencyCheck();
#endif
	return *this;
}

void ClusterGraph::copyClusterTree(const ClusterGraph& C, const Graph& G,
		ClusterArray<cluster>& originalClusterTable, std::function<node(node)> nodeMap) {
	for (cluster c : C.clusters) {
		if (c == C.m_rootCluster) {
			originalClusterTable[c] = m_rootCluster;
			// does not really need to be assigned HERE in for
			m_rootCluster->m_depth = 1;
			OGDF_ASSERT(c->depth() == 1);
		} else {
			originalClusterTable[c] = newCluster();
			originalClusterTable[c]->m_depth = c->depth();
		}
	}

	for (cluster c : C.clusters) {
		if (c != C.m_rootCluster) {
			originalClusterTable[c]->m_parent = originalClusterTable[c->m_parent];
			originalClusterTable[c->m_parent]->children.pushBack(originalClusterTable[c]);
			originalClusterTable[c]->m_it = originalClusterTable[c->m_parent]->getChildren().rbegin();
		}
	}

	for (node v : C.constGraph().nodes) {
		reassignNode(nodeMap(v), originalClusterTable[C.clusterOf(v)]);
	}

	copyLCA(C);
}

// Copy Function
void ClusterGraph::shallowCopy(const ClusterGraph& C) {
	resizeArrays(C.numberOfClusters() + 1);

	const Graph& G = C.constGraph();
	initGraph(G);

	m_updateDepth = C.m_updateDepth;
	m_depthUpToDate = C.m_depthUpToDate;

	ClusterArray<cluster> originalClusterTable(C);
	copyClusterTree(C, G, originalClusterTable);
}

// Initialize the graph
void ClusterGraph::initGraph(const Graph& G) {
	reregister(&G);

	m_adjAvailable = false;

	// root cluster must always get id 0
#ifdef OGDF_DEBUG
	m_rootCluster = new ClusterElement(this, 0);
#else
	m_rootCluster = new ClusterElement(0);
#endif

	OGDF_ASSERT(numberOfClusters() == 0);

	m_rootCluster->m_depth = 1;
	m_clusterIdCount++;
	m_nodeMap.init(G, m_rootCluster);
	m_itMap.init(G);
	// assign already existing nodes to root cluster (new nodes are assigned over nodeadded)
	for (node v : G.nodes) {
		m_itMap[v] = m_rootCluster->getNodes().pushBack(v);
	}

	clusters.pushBack(m_rootCluster);
	keyAdded(m_rootCluster);
	// we do not notify observers about the root cluster
}

void ClusterGraph::reinitGraph(const Graph& G) {
#ifdef OGDF_HEAVY_DEBUG
	G.consistencyCheck();
#endif

	if (numberOfClusters() != 0) {
		doClear();
	}

	initGraph(G); //already constructs root cluster, reassign
}

// Copy Function
void ClusterGraph::deepCopy(const ClusterGraph& C, Graph& G) {
	ClusterArray<cluster> originalClusterTable(C);
	NodeArray<node> originalNodeTable(C.constGraph());
	EdgeArray<edge> edgeCopy(C.constGraph());
	deepCopy(C, G, originalClusterTable, originalNodeTable, edgeCopy);
}

void ClusterGraph::deepCopy(const ClusterGraph& C, Graph& G,
		ClusterArray<cluster>& originalClusterTable, NodeArray<node>& originalNodeTable) {
	EdgeArray<edge> edgeCopy(C.constGraph());
	deepCopy(C, G, originalClusterTable, originalNodeTable, edgeCopy);
}

void ClusterGraph::deepCopy(const ClusterGraph& C, Graph& G,
		ClusterArray<cluster>& originalClusterTable, NodeArray<node>& originalNodeTable,
		EdgeArray<edge>& edgeCopy) {
	m_updateDepth = C.m_updateDepth;
	m_depthUpToDate = C.m_depthUpToDate;
	G.clear();
	resizeArrays(C.numberOfClusters() + 1);
	initGraph(G);
	G.insert(C.constGraph(), originalNodeTable, edgeCopy);
	if (!originalClusterTable.registeredAt()) {
		originalClusterTable.init(C, nullptr);
	}
	copyClusterTree(C, G, originalClusterTable, originalNodeTable);
}

//We search for the lowest common cluster of a set of nodes.
//We first compute the common path of two nodes, then update path if root
//path from other nodes hits it .
//We always stop if we encounter root cluster.
cluster ClusterGraph::commonCluster(SList<node>& nodes) {
	//worst case running time #nodes x clustertreeheight-1
	//always <= complete tree run
	//we could even use pathcompression...
	//at any time, we stop if root is encountered as lowest
	//common cluster of a node subset


	if (nodes.empty()) {
		return nullptr;
	}

	//For simplicity, we use cluster arrays
	ClusterArray<int> commonPathHit(*this, 0); //count for clusters path hits
	cluster pathCluster;
	SListIterator<node> sIt = nodes.begin();
	node v1 = *sIt;
	if (nodes.size() == 1) {
		return clusterOf(v1);
	}

	++sIt;
	cluster lowestCommon = commonCluster(v1, *sIt);
	commonPathHit[lowestCommon] = 2;
	pathCluster = lowestCommon;
	while (pathCluster->parent()) {
		pathCluster = pathCluster->parent();
		commonPathHit[pathCluster] = 2;
	}

	// we save direct lca access, it also lies on a runs hit path from root
	for (int runs = 2; runs < nodes.size() && lowestCommon != m_rootCluster; ++runs) {
		// runs is the number of nodes already considered
		++sIt;
		pathCluster = clusterOf(*sIt);
		while (commonPathHit[pathCluster] == 0) {
			OGDF_ASSERT(pathCluster->parent() != nullptr);
			pathCluster = pathCluster->parent();
		}
		// assign new (maybe same) lowest common
		if (commonPathHit[pathCluster] == runs) {
			lowestCommon = pathCluster;
		}
		++commonPathHit[pathCluster];
		if (pathCluster == m_rootCluster) {
			return m_rootCluster;
		}
		// update hits in path to root
		while (pathCluster->parent()) {
			pathCluster = pathCluster->parent();
			++commonPathHit[pathCluster];
		}
	}

	return lowestCommon;
}

//note that eL is directed from v to w
cluster ClusterGraph::commonClusterAncestorsPath(node v, node w, cluster& c1, cluster& c2,
		List<cluster>& eL) const {
	OGDF_ASSERT(v->graphOf() == &constGraph());
	OGDF_ASSERT(w->graphOf() == &constGraph());

	cluster cv = clusterOf(v);
	cluster cw = clusterOf(w);

	//clusters from v and w to common
	List<cluster> vList;
	List<cluster> wList;

	//CASE1 no search necessary
	//if both nodes are in the same cluster, we return this cluster
	//and have to check if c1 == c2 to have a (v,w) representation edge
	if (cv == cw) {
		c1 = c2 = cv;
		eL.pushBack(c1);
		return cv;
	}

	if (m_lcaNumber == std::numeric_limits<int>::max() - 1) {
		m_lcaNumber = 0;
	} else {
		m_lcaNumber++;
	}
	if (!m_lcaSearch) {
		m_lcaSearch.reset(new ClusterArray<int>(*this, -1));
		m_vAncestor.reset(new ClusterArray<cluster>(*this, nullptr));
		m_wAncestor.reset(new ClusterArray<cluster>(*this, nullptr));
	}

	//CASE2: one of the nodes hangs at root: save root as ancestor
	//any other case: save cluster of node as ancestor, too, to check this
	//case:: common = xCluster != yCluster
#if 0
	(*m_vAncestor)[rootCluster()] = rootCluster();
	(*m_wAncestor)[rootCluster()] = rootCluster();
#endif
	(*m_vAncestor)[cv] = nullptr;
	(*m_wAncestor)[cw] = nullptr;

	//we rely on the fact all nodes are in the rootcluster or
	//that parent is initialized to zero to terminate

	//we start with different clusters due to CASE1
	//save the ancestor information
	(*m_lcaSearch)[cw] = m_lcaNumber; //not really necessary, we won't return
	(*m_lcaSearch)[cv] = m_lcaNumber;
	vList.pushBack(cv);
	wList.pushBack(cw);

	// we break and return if we find a common node
	// before we reach the rootcluster
	do {
		if (cv->parent()) { // if root not reached on cv-path
			(*m_vAncestor)[cv->parent()] = cv;
			cv = cv->parent();
			// was cv visited on path from w
			if ((*m_lcaSearch)[cv] == m_lcaNumber) {
				c1 = (*m_vAncestor)[cv];
				c2 = (*m_wAncestor)[cv];

				// setup list
				for (auto c : vList) {
					eL.pushBack(c);
				}

				ListReverseIterator<cluster> itC;
				for (itC = wList.rbegin(); itC.valid() && *itC != cv; ++itC) {
					;
				}
				for (; itC.valid(); ++itC) {
					eL.pushBack(*itC);
				}

				return cv;
			}
			vList.pushBack(cv);
			(*m_lcaSearch)[cv] = m_lcaNumber;
		}

		if (cw->parent()) { // if root not reached on cw-path
			(*m_wAncestor)[cw->parent()] = cw;
			cw = cw->parent();
			// was cw visited on path from v
			if ((*m_lcaSearch)[cw] == m_lcaNumber) {
				c1 = (*m_vAncestor)[cw];
				c2 = (*m_wAncestor)[cw];

				// setup list
				for (auto itC = vList.begin(); itC.valid() && *itC != cw; ++itC) {
					eL.pushBack(*itC);
				}
				eL.pushBack(cw);

				for (cluster c : reverse(wList)) {
					eL.pushBack(c);
				}

				return cw;
			}
			wList.pushBack(cw);
			(*m_lcaSearch)[cw] = m_lcaNumber;
		}
	} while (cv->parent() || cw->parent());

	// v,w should be at least together in the rootcluster
	c1 = (*m_vAncestor)[rootCluster()];
	c2 = (*m_wAncestor)[rootCluster()];

	return rootCluster();
}

void ClusterGraph::copyLCA(const ClusterGraph& C) {
	if (C.m_lcaSearch) {
		//otherwise, initialization won't work
		m_lcaSearch.reset(new ClusterArray<int>(*this, -1)); //(*C.m_lcaSearch);
		m_vAncestor.reset(new ClusterArray<cluster>(*this, nullptr));
		m_wAncestor.reset(new ClusterArray<cluster>(*this, nullptr));
		//setting of clusters is not necessary!
	}
}

// check the graph for empty clusters
// we never set rootcluster to be one of the empty clusters!!
void ClusterGraph::emptyClusters(SList<cluster>& emptyCluster, SList<cluster>* checkCluster) {
	if (checkCluster) {
		fillEmptyClusters(emptyCluster, *checkCluster);
	} else {
		fillEmptyClusters(emptyCluster, clusters);
	}

	// other clusters can get empty, too, if we delete these
	ClusterArray<int> delCount(*this, 0);
	SList<cluster> emptyParent;
	for (cluster c : emptyCluster) {
		// count deleted children
		cluster runc = c->parent();
		if (runc) // is always the case as long as root was not inserted to list
		{
			delCount[runc]++;
			while ((runc->nCount() == 0) && (runc->cCount() == delCount[runc])) {
				if (runc == rootCluster()) {
					break;
				}
				emptyParent.pushBack(runc);
				runc = runc->parent();
				delCount[runc]++;
			}
		}
	}

	emptyCluster.conc(emptyParent);
	// for reinsertion, start at emptycluster's back
}

// Inserts a new cluster prescribing its parent
cluster ClusterGraph::newCluster(cluster parent, int id) {
	OGDF_ASSERT(parent);
	cluster c;
	if (id > 0) {
		c = newCluster(id);
	} else {
		c = newCluster();
	}
	parent->children.pushBack(c);
	c->m_it = parent->getChildren().rbegin();
	c->m_parent = parent;
	c->m_depth = parent->depth() + 1;

	return c;
}

//Insert a new cluster with given ID, precondition: id not used
//has to be updated in the same way as newcluster()
cluster ClusterGraph::newCluster(int id) {
	m_adjAvailable = false;
	m_postOrderStart = nullptr;
	if (id >= m_clusterIdCount) {
		m_clusterIdCount = id + 1;
	}
#ifdef OGDF_DEBUG
	cluster c = new ClusterElement(this, id);
#else
	cluster c = new ClusterElement(id);
#endif
	clusters.pushBack(c);
	keyAdded(c);
	for (ClusterGraphObserver* obs : Obs::getObservers()) {
		obs->clusterAdded(c);
	}

	return c;
}

// Inserts a new cluster
//has to be updated in the same way as newcluster(id)
cluster ClusterGraph::newCluster() {
	m_adjAvailable = false;
	m_postOrderStart = nullptr;
#ifdef OGDF_DEBUG
	cluster c = new ClusterElement(this, m_clusterIdCount++);
#else
	cluster c = new ClusterElement(m_clusterIdCount++);
#endif
	clusters.pushBack(c);
	keyAdded(c);
	for (ClusterGraphObserver* obs : Obs::getObservers()) {
		obs->clusterAdded(c);
	}

	return c;
}

cluster ClusterGraph::createEmptyCluster(const cluster parent, int clusterId) {
	//if no id given, use next free id
	if (clusterId < 0) {
		clusterId = m_clusterIdCount;
	}
	//create the new cluster
	cluster cnew;
	if (parent) {
		cnew = newCluster(parent, clusterId);
	} else {
		cnew = newCluster(m_rootCluster, clusterId);
	}
	return cnew;
}

cluster ClusterGraph::createCluster(const SList<node>& nodes, const cluster parent) {
	cluster c;
	if (m_allowEmptyClusters) {
		c = doCreateCluster(nodes, parent);
		return c;

	} else {
		SList<cluster> emptyCluster;

		c = doCreateCluster(nodes, emptyCluster, parent);

		for (cluster ec : emptyCluster) {
			delCluster(ec);
			//root cluster can never be empty, as we deleted a node
		}
	}

	return c;
}

cluster ClusterGraph::doCreateCluster(const SList<node>& nodes, const cluster parent, int clusterId) {
	if (nodes.empty()) {
		return nullptr;
	}

	//if no id given, use next free id
	if (clusterId < 0) {
		clusterId = m_clusterIdCount;
	}
	//create the new cluster
	cluster cnew;
	if (parent) {
		cnew = newCluster(parent, clusterId);
	} else {
		cnew = newCluster(m_rootCluster, clusterId);
	}

	//insert nodes in new cluster
	for (node v : nodes) {
		reassignNode(v, cnew);
	}

	return cnew;
}

cluster ClusterGraph::doCreateCluster(const SList<node>& nodes, SList<cluster>& emptyCluster,
		const cluster parent, int clusterId) {
	// Even if m_allowEmptyClusters is set we check if a cluster
	// looses all of its nodes and has
	// no more entries and childs. This can be used for special cluster
	// object handling or for deletion if m_allowEmptyClusters is not set
	// if it is not the new parent, it can be deleted
	// running time max(#cluster, length(nodelist))
	// TODO: Parameter, der dies auslaesst, da hohe Laufzeit
	// hier macht das nur Sinn, wenn es schneller ist als for all clusters,
	// sonst koennte man es ja auch aussen testen, aber bisher ist es nicht
	// schneller implementiert
	// Vorgehen: hash auf cluster index, falls nicht gesetzt, in liste einfuegen
	// und als checkcluster an emptycluster uebergeben

	if (nodes.empty()) {
		return nullptr;
	}

	//if no id given, use next free id
	if (clusterId < 0) {
		clusterId = m_clusterIdCount;
	}
	//create the new cluster
	cluster cnew;
	if (parent) {
		cnew = newCluster(parent, clusterId);
	} else {
		cnew = newCluster(m_rootCluster, clusterId);
	}

	//insert nodes in new cluster
	for (node v : nodes) {
		reassignNode(v, cnew);
	}

	//should be: only for changed clusters (see comment above)
	//it is important to save the cluster in an order
	//that allows deletion as well as reinsertion
	emptyClusters(emptyCluster);
	//for reinsertion, start at emptycluster's back

	return cnew;
}

// Deletes cluster c
// All subclusters become children of parent cluster
// Precondition: c is not the root cluster
// updating of cluster depth information pumps running time
// up to worst case O(#C)
void ClusterGraph::delCluster(cluster c) {
	OGDF_ASSERT(c != nullptr);
	OGDF_ASSERT(c->graphOf() == this);
	OGDF_ASSERT(c != m_rootCluster);

	// notify observers
	for (ClusterGraphObserver* obs : Obs::getObservers()) {
		obs->clusterDeleted(c);
	}
	keyRemoved(c);

	m_postOrderStart = nullptr;

	c->m_parent->children.del(c->m_it);
	c->m_it = ListIterator<cluster>();

	while (!c->children.empty()) {
		cluster trace = c->children.popFrontRet();
		trace->m_parent = c->m_parent;
		trace->m_parent->children.pushBack(trace);
		trace->m_it = trace->m_parent->getChildren().rbegin();

		//only recompute depth if option set and it makes sense
		if (m_updateDepth && m_depthUpToDate) {
			//update depth for all children in subtree
			OGDF_ASSERT(trace->depth() == trace->parent()->depth() + 2);
			pullUpSubTree(trace);
			//could just set depth-1 here
#if 0
			trace->depth() = trace->parent()->depth()+1;
#endif
		} else {
			m_depthUpToDate = false;
		}
	}
	while (!c->nodes.empty()) {
		node v = c->nodes.popFrontRet();
		m_nodeMap[v] = nullptr;
		reassignNode(v, c->m_parent);
	}

	clusters.del(c);
#ifdef OGDF_HEAVY_DEBUG
	consistencyCheck();
#endif
}

//pulls up depth of subtree located at c by one
//precondition: depth is consistent
//we dont ask for depthuptodate since the caller needs
//to know for himself if he wants the tree to be pulled
//for any special purpose
void ClusterGraph::pullUpSubTree(cluster c) {
	c->m_depth = c->depth() - 1;
	for (cluster ci : c->getChildren()) {
		pullUpSubTree(ci);
	}
}

void ClusterGraph::doClear() { // TODO merge with clear
	//split condition
	m_lcaSearch.reset();
	m_vAncestor.reset();
	m_wAncestor.reset();
	if (numberOfClusters() != 0) {
		clearClusterTree(m_rootCluster);
		clusters.del(m_rootCluster);
	}
	for (ClusterGraphObserver* obs : Obs::getObservers()) {
		obs->clustersCleared();
	}
	//no clusters, so we can restart at 0
	m_clusterIdCount = 0;
	keysCleared();
}

//don't delete root cluster
void ClusterGraph::clear() {
	m_lcaSearch.reset();
	m_vAncestor.reset();
	m_wAncestor.reset();
	if (numberOfClusters() != 0) {
		//clear the cluster structure under root cluster
		clearClusterTree(m_rootCluster);
		//now delete all rootcluster entries
		while (!m_rootCluster->nodes.empty()) {
			node v = m_rootCluster->nodes.popFrontRet();
			m_nodeMap[v] = nullptr;
		}
	}
	for (ClusterGraphObserver* obs : Obs::getObservers()) {
		obs->clustersCleared();
	}
	//no child clusters, so we can restart at 1
	m_clusterIdCount = 1;
	keysCleared();
	keyAdded(m_rootCluster);
}

// Removes the Clustering of a Tree and frees the allocated memory
void ClusterGraph::clearClusterTree(cluster c) {
	cluster parent = c->parent();
	m_postOrderStart = nullptr;

	List<node> attached;
	recurseClearClusterTreeOnChildren(c, attached);

	if (parent != nullptr) {
		for (ClusterGraphObserver* obs : Obs::getObservers()) {
			obs->clusterDeleted(c);
		}
		keyRemoved(c);
		for (node v : attached) {
			m_nodeMap[v] = parent;
			parent->nodes.pushBack(v);
			m_itMap[v] = parent->getNodes().rbegin();
		}
		clusters.del(c);
	} else if (c == m_rootCluster) {
		for (node v : attached) {
			m_nodeMap[v] = m_rootCluster;
			m_rootCluster->nodes.pushBack(v);
			m_itMap[v] = m_rootCluster->getNodes().rbegin();
		}
		m_rootCluster->children.clear();
	}
}

void ClusterGraph::clearClusterTree(cluster c, List<node>& attached) {
	for (ClusterGraphObserver* obs : Obs::getObservers()) {
		obs->clusterDeleted(c);
	}
	keyRemoved(c);
	attached.conc(c->nodes);
	recurseClearClusterTreeOnChildren(c, attached);
	clusters.del(c);
}

int ClusterGraph::treeDepth() const {
	//initialize depth at first call
	if (m_updateDepth && !m_depthUpToDate) {
		computeSubTreeDepth(rootCluster());
	}
	if (!m_updateDepth) {
		OGDF_THROW(AlgorithmFailureException);
	}
	int l_depth = 1;

	for (cluster c : clusters) {
		if (c->depth() > l_depth) {
			l_depth = c->depth();
		}
	}

	return l_depth;
}

//reassign cluster depth for clusters in subtree rooted at c
void ClusterGraph::computeSubTreeDepth(cluster c) const {
	if (c == rootCluster()) {
		m_depthUpToDate = true;
	}

	c->m_depth = (c->parent() == nullptr) ? 1 : c->parent()->depth() + 1;

	for (cluster child : c->children) {
		computeSubTreeDepth(child);
	}
}

//move cluster from old parent to an other
void ClusterGraph::moveCluster(cluster c, cluster newParent) {
	if (c == rootCluster()) {
		return;
	}
	if ((c == nullptr) || (newParent == nullptr)) {
		return; //no cheap tricks
	}
	if (c->parent() == newParent) {
		return; //no work to do
	}

	cluster oldParent = c->parent();
	//we dont move root
	OGDF_ASSERT(oldParent);

	//check if we move to a descendant
	cluster crun = newParent->parent();
	bool descendant = false;
	while (crun) {
		if (crun == c) {
			descendant = true;
			break;
		}
		crun = crun->parent();
	}

	//do not allow to move empty clusters to descendants
	if (descendant && (c->nCount() == 0)) {
		return;
	}

	//temporarily only recompute postorder for all clusters

	oldParent->children.del(c->m_it);
	newParent->children.pushBack(c);
	c->m_it = newParent->getChildren().rbegin();
	c->m_parent = newParent;

	//update the cluster depth information in the subtree
	//If moved to descendant, recompute
	//depth for parent (including all brother trees)
	if (descendant) {
		//how do we move:
		//only entries with c? => may be empty
		//we currently dont allow this, because it makes
		//no sense, you could just delete the cluster or move
		//the children
		//move all children to oldparent

		while (!c->children.empty()) {
			cluster child = c->children.popFrontRet();
			child->m_parent = oldParent;
			child->m_parent->children.pushBack(child);
			child->m_it = child->m_parent->getChildren().rbegin();
#if 0
			child++;
#endif
		}

		//recompute depth only if option set AND it makes sense at that point
		if (m_updateDepth && m_depthUpToDate) {
			computeSubTreeDepth(oldParent);
		} else {
			m_depthUpToDate = false;
		}
	} else {
		if (m_updateDepth && m_depthUpToDate) {
			computeSubTreeDepth(c);
		} else {
			m_depthUpToDate = false;
		}
	}

	// update postorder for new parent
	// we only recompute postorder for all clusters
	// because of special cases like move to descendant...
	postOrder();

	m_adjAvailable = false;

#if 0
	checkPostOrder();
#endif
}

//leftmostcluster in subtree rooted at c, has postorderpred for subtree
cluster ClusterGraph::leftMostCluster(cluster c) const {
	cluster result = c;
	if (!c) {
		return nullptr;
	}
	while (!result->children.empty()) {
		result = result->children.front();
	}
	return result;
}

//searches for predecessor of SUBTREE at c
cluster ClusterGraph::postOrderPredecessor(cluster c) const {
	//all clusters on a path from root to leftmost cluster in tree
	//have no predecessor for their subtree
	cluster run = c;
	do {
		//predecessor of clustertree is	0
		if (run == m_rootCluster) {
			return nullptr;
		}

		ListConstIterator<cluster> it = run->m_it;
		//a child to the left is the immediate predecessor,
		//otherwise we go one level up
		if (it == (run->m_parent)->children.begin()) {
			run = run->parent();
		} else {
			return *(it.pred());
		}

	} while (run);

	return nullptr;
}

void ClusterGraph::nodeDeleted(node v) {
	bool cRemove = false;
	cluster c = clusterOf(v);
	if (!c) {
		return;
	}
	//never allow totally empty cluster
#if 0
	if ((emptyOnNodeDelete(c)) &&
		(c != rootCluster()) ) {
		cRemove = true;
	}
#endif
	unassignNode(v);
	if (cRemove && !m_allowEmptyClusters) //parent exists
	{
		cluster nonEmpty = c->parent();
		cluster cRun = nonEmpty;
		delCluster(c);
		while ((cRun != rootCluster()) && (cRun->nCount() + cRun->cCount() == 0)) {
			nonEmpty = cRun->parent();
			delCluster(cRun);
			cRun = nonEmpty;
		}
	}
}

//node assignment
//Assigns a node to a new cluster
void ClusterGraph::assignNode(node v, cluster c) {
	m_adjAvailable = false;
	m_postOrderStart = nullptr;
	m_nodeMap[v] = c;
	c->nodes.pushBack(v);
	m_itMap[v] = c->getNodes().rbegin();
}

//Reassigns a node to a new cluster
void ClusterGraph::reassignNode(node v, cluster c) {
	OGDF_ASSERT(v->graphOf() == &constGraph());
	OGDF_ASSERT(c->graphOf() == this);

	unassignNode(v);
	m_nodeMap[v] = c;
	c->nodes.pushBack(v);
	m_itMap[v] = c->getNodes().rbegin();
}

//Unassigns a node of cluster
//Note: Nodes can already be unassigned by the nodeDeleted function.
void ClusterGraph::unassignNode(node v) {
	m_adjAvailable = false;
	m_postOrderStart = nullptr;

	removeNodeAssignment(v);
}

// Start function for post order
void ClusterGraph::postOrder() const {
	SListPure<cluster> L;
	postOrder(m_rootCluster, L);
	cluster c = nullptr;
	cluster prev = L.popFrontRet();
	prev->m_pPrev = nullptr;
	m_postOrderStart = prev;
	while (!L.empty()) {
		c = L.popFrontRet();
		prev->m_pNext = c;
		c->m_pPrev = prev;
		prev = c;
	}
	if (c != nullptr) {
		c->m_pNext = nullptr;
	} else {
		m_postOrderStart->m_pNext = nullptr;
	}

#ifdef OGDF_DEBUG
	for (cluster cl : clusters) {
		cluster cp = leftMostCluster(cl);
		OGDF_ASSERT(cp->pPred() == postOrderPredecessor(cl));
	}
#endif
}


#ifdef OGDF_DEBUG
void ClusterGraph::checkPostOrder() const {
	SListPure<cluster> L;
	postOrder(m_rootCluster, L);
	cluster c = nullptr;
	cluster prev = L.popFrontRet();
	OGDF_ASSERT(prev->m_pPrev == nullptr);

	while (!L.empty()) {
		c = L.popFrontRet();
		OGDF_ASSERT(prev->m_pNext == c);
		OGDF_ASSERT(c->m_pPrev == prev);
		prev = c;
	}
	if (c != nullptr) {
		OGDF_ASSERT(c->m_pNext == nullptr);
	} else {
		OGDF_ASSERT(m_postOrderStart->m_pNext == nullptr);
	}
}
#endif


// Recursive function for post order
void ClusterGraph::postOrder(cluster c, SListPure<cluster>& L) const {
	for (cluster ci : c->children) {
		postOrder(ci, L);
	}

	L.pushBack(c);
}


#ifdef OGDF_DEBUG
void ClusterGraph::consistencyCheck() const {
	ClusterArray<bool> visitedClusters((*this), false);
	NodeArray<bool> visitedNodes(constGraph(), false);
	int visitedClustersC = 0, visitedNodesC = 0;
	ClusterArray<int> clusterDepth;
	if (m_updateDepth && m_depthUpToDate) {
		clusterDepth.init((*this), -1);
	}

	for (cluster c = firstPostOrderCluster(); c != nullptr; c = c->pSucc()) {
		OGDF_ASSERT(!visitedClusters[c]);
		visitedClusters[c] = true;
		visitedClustersC++;
		if (m_updateDepth && m_depthUpToDate) {
			clusterDepth[c] = c->depth();
		}
		OGDF_ASSERT((c->parent() == nullptr) == (c == rootCluster()));

		for (node v : c->nodes) {
			OGDF_ASSERT(clusterOf(v) == c);
			OGDF_ASSERT(!visitedNodes[v]);
			visitedNodes[v] = true;
			visitedNodesC++;
		}
		for (cluster child : c->children) {
			OGDF_ASSERT(visitedClusters[child]);
			OGDF_ASSERT(child->parent() == c);
		}
	}
	OGDF_ASSERT(visitedClustersC == numberOfClusters());

	if (m_updateDepth && m_depthUpToDate) {
		computeSubTreeDepth(rootCluster());
	}
	for (cluster c : clusters) {
		OGDF_ASSERT(visitedClusters[c]);
		if (!m_allowEmptyClusters) {
			OGDF_ASSERT(c->nCount() + c->cCount() >= 1);
		}
		if (m_updateDepth && m_depthUpToDate) {
			OGDF_ASSERT(clusterDepth[c] == c->depth());
		}
	}

	OGDF_ASSERT(visitedNodesC == constGraph().numberOfNodes());
	for (node v : constGraph().nodes) {
		OGDF_ASSERT(visitedNodes[v]);
	}

	if (m_adjAvailable) {
		ClusterArray<List<edge>> boundaryCrossingEdges(*this);
		for (edge e : constGraph().edges) {
			List<cluster> path;
			cluster lca = commonClusterPath(e->source(), e->target(), path);
			for (cluster c : path) {
				boundaryCrossingEdges[c].pushBack(e);
			}
			OGDF_ASSERT(boundaryCrossingEdges[lca].back() == e);
			boundaryCrossingEdges[lca].popBack();
		}
		EdgeArray<cluster> clusterAssociation(constGraph(), nullptr);
		for (cluster c : clusters) {
			OGDF_ASSERT(c->adjEntries.size() == boundaryCrossingEdges[c].size());
			for (adjEntry adj : c->adjEntries) {
				OGDF_ASSERT(adj->graphOf() == getGraph());
				OGDF_ASSERT(!c->isDescendant(clusterOf(adj->twinNode()), true));
				OGDF_ASSERT(c->isDescendant(clusterOf(adj->theNode()), true));
				OGDF_ASSERT(clusterAssociation[adj] == nullptr);
				clusterAssociation[adj] = c;
			}
			for (edge e : boundaryCrossingEdges[c]) {
				OGDF_ASSERT(clusterAssociation[e] == c);
				clusterAssociation[e] = nullptr;
			}
		}
	}
}
#endif

bool ClusterGraph::representsCombEmbedding() const {
	if (!m_adjAvailable) {
		return false;
	}
#ifdef OGDF_DEBUG
	consistencyCheck();
#endif
	GraphCopySimple gcopy(constGraph());
	gcopy.setOriginalEmbedding();
	planarizeClusterBorderCrossings(*this, gcopy, nullptr,
			[&gcopy](edge e) -> edge { return gcopy.copy(e); });
	return gcopy.representsCombEmbedding();
}

bool ClusterGraph::representsConnectedCombEmbedding() const {
	if (!m_adjAvailable) {
		return false;
	}
#ifdef OGDF_DEBUG
	consistencyCheck();
#endif

	for (cluster c = firstPostOrderCluster(); c != m_rootCluster; c = c->pSucc()) {
		for (auto it = c->adjEntries.begin(); it.valid(); it++) {
			adjEntry adj = *it;
			adjEntry succAdj = *c->adjEntries.cyclicSucc(it);

			// run along the outer face of the cluster until you find the next outgoing edge
			adjEntry next = adj->cyclicSucc();
#ifdef OGDF_DEBUG
			int max = 2 * (constGraph().numberOfEdges() + 1);
			while (next != succAdj) {
				// this guards against an infinite loop if the underlying embedding is non-planar
				OGDF_ASSERT(max >= 0);
				max--;
#else
			while (next != succAdj) {
#endif
				if (next == adj->twin()) {
					// we searched the whole face but didn't find our successor along the cluster border
					return false;
				}
				next = next->twin()->cyclicSucc();
			}
		}
	}

	return true;
}

void ClusterGraph::registrationChanged(const ogdf::Graph* newG) {
	m_lcaSearch.reset();
	m_vAncestor.reset();
	m_wAncestor.reset();
	m_lcaNumber = 0;
	m_adjAvailable = false;
	for (cluster c : clusters) {
		c->nodes.clear();
		c->adjEntries.clear();
	}
}

cluster ClusterGraph::chooseCluster(std::function<bool(cluster)> includeCluster,
		bool isFastTest) const {
	return *chooseIteratorFrom<internal::GraphObjectContainer<ClusterElement>, cluster>(
			const_cast<internal::GraphObjectContainer<ClusterElement>&>(clusters),
			[&](const cluster& c) { return includeCluster(c); }, isFastTest);
}

std::ostream& operator<<(std::ostream& os, cluster c) {
	if (c != nullptr) {
		os << c->index();
	} else {
		os << "nil";
	}
	return os;
}

void planarizeClusterBorderCrossings(const ClusterGraph& CG, Graph& G,
		EdgeArray<List<std::pair<adjEntry, cluster>>>* subdivisions,
		const std::function<edge(edge)>& translate) {
	OGDF_ASSERT(CG.adjAvailable());
	for (cluster c = CG.firstPostOrderCluster(); c != CG.rootCluster(); c = c->pSucc()) {
		adjEntry prev_ray = nullptr, first_ray = nullptr;
		for (adjEntry adj : c->adjEntries) {
			bool reverse = adj->isSource();
			edge the_edge = translate(adj->theEdge());
			if (reverse) {
				G.reverseEdge(the_edge);
			}
			edge new_edge = G.split(the_edge);
			adjEntry spike_to_src =
					the_edge->adjTarget(); // adjacency of the new node toward the source of the_edge
			if (subdivisions != nullptr) {
				if (reverse) {
					(*subdivisions)[adj->theEdge()].emplaceBack(spike_to_src, c);
				} else {
					(*subdivisions)[adj->theEdge()].emplaceFront(spike_to_src, c);
				}
			}
			if (reverse) {
				G.reverseEdge(the_edge);
				G.reverseEdge(new_edge);
			}

			if (prev_ray != nullptr) {
				G.newEdge(prev_ray, spike_to_src->cyclicPred());
			}
			prev_ray = spike_to_src;
			if (first_ray == nullptr) {
				first_ray = spike_to_src;
			}
		}
		if (first_ray != nullptr) {
			G.newEdge(prev_ray, first_ray->cyclicPred());
		}
	}
}

}
