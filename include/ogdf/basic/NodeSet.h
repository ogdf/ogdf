/*
 * $Revision: 3951 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2014-03-03 13:57:46 +0100 (Mon, 03 Mar 2014) $
 ***************************************************************/

/** \file
 * \brief Declaration and implementation of class NodeSetSimple,
 *        NodeSetPure and NodeSet
 *
 * \author Carsten Gutwenger
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


#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_NODE_SET_H
#define OGDF_NODE_SET_H


#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/SList.h>



namespace ogdf {


//! Simple node sets.
/**
 * A node set maintains a subset \a S of the nodes contained in an associated
 * graph. This kind of node set only provides efficient operation for testing
 * membership, insertion, and clearing the set.
 *
 * \sa
 *   - NodeSet, NodeSetPure
 *   - FaceSet, FaceSetPure, FaceSetSimple
 */
class OGDF_EXPORT NodeSetSimple {
public:
	//! Creates an empty node set associated with graph \a G.
	NodeSetSimple(const Graph &G) : m_isContained(G,false) { }

	// destructor
	~NodeSetSimple() { }

	//! Inserts node \a v into \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	void insert(node v) {
		OGDF_ASSERT(v->graphOf() == m_isContained.graphOf());
		bool &isContained = m_isContained[v];
		if (isContained == false) {
			isContained = true;
			m_nodes.pushFront(v);
		}
	}

	//! Removes all nodes from \a S.
	/**
	 * After this operation, \a S is empty and still associated with the same graph.
	 * The runtime of this operations is O(k), where k is the number of nodes in \a S
	 * before this operation.
	 */
	void clear() {
		SListIterator<node> it;
		for(it = m_nodes.begin(); it.valid(); ++it) {
			m_isContained[*it] = false;
		}
		m_nodes.clear();
	}

	//! Returns true if node \a v is contained in \a S, false otherwise.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	bool isMember(node v) const {
		OGDF_ASSERT(v->graphOf() == m_isContained.graphOf());
		return m_isContained[v];
	}

	//! Returns a reference to the list of nodes contained in \a S.
	/**
	 * This list can be used for iterating over all nodes in \a S.
	 */
	const SListPure<node> &nodes() const {
		return m_nodes;
	}

private:
	//! m_isContained[v] is true iff \a v is contained in \a S.
	NodeArray<bool> m_isContained;

	//! The list of nodes contained in \a S.
	SListPure<node> m_nodes;
};



//! Node sets.
/**
 * A node set maintains a subset \a S of the nodes contained in an associated
 * graph. This kind of node set provides efficient operations for testing
 * membership, insertion and deletion of elements, and clearing the set.
 *
 * In contrast to NodeSet, a NodeSetPure does not provide efficient access
 * to the number of nodes stored in the set.
 *
 * \sa
 *   - NodeSet, NodeSetSimple
 *   - FaceSet, FaceSetPure, FaceSetSimple
 */
class OGDF_EXPORT NodeSetPure {
public:
	//! Creates an empty node set associated with graph \a G.
	NodeSetPure(const Graph &G) : m_it(G,ListIterator<node>()) { }

	// destructor
	~NodeSetPure() { }

	//! Inserts node \a v into \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	void insert(node v) {
		OGDF_ASSERT(v->graphOf() == m_it.graphOf());
		ListIterator<node> &itV = m_it[v];
		if (!itV.valid())
			itV = m_nodes.pushBack(v);
	}

	//! Removes node \a v from \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	void remove(node v) {
		OGDF_ASSERT(v->graphOf() == m_it.graphOf());
		ListIterator<node> &itV = m_it[v];
		if (itV.valid()) {
			m_nodes.del(itV);
			itV = ListIterator<node>();
		}
	}


	//! Removes all nodes from \a S.
	/**
	 * After this operation, \a S is empty and still associated with the same graph.
	 * The runtime of this operations is O(k), where k is the number of nodes in \a S
	 * before this operation.
	 */
	void clear() {
		ListIterator<node> it;
		for(it = m_nodes.begin(); it.valid(); ++it) {
			m_it[*it] = ListIterator<node>();
		}
		m_nodes.clear();
	}


	//! Returns true if node \a v is contained in \a S, false otherwise.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	bool isMember(node v) const {
		OGDF_ASSERT(v->graphOf() == m_it.graphOf());
		return m_it[v].valid();
	}

	//! Returns a reference to the list of nodes contained in \a S.
	/**
	 * This list can be used for iterating over all nodes in \a S.
	 */
	const ListPure<node> &nodes() const {
		return m_nodes;
	}

	//! Copy constructor.
	NodeSetPure(const NodeSetPure& V) : m_it(*V.m_it.graphOf(), ListIterator<node>()) {
		forall_listiterators(node, it, V.m_nodes) {
			insert(*it);
		}
	}

	//! Assignment operator.
	NodeSetPure &operator=(const NodeSetPure &V) {
		m_nodes.clear();
		m_it.init(*V.m_it.graphOf());
		forall_listiterators(node, it, V.m_nodes) {
			insert(*it);
		}
		return *this;
	}

private:
	//! m_it[v] contains the list iterator pointing to \a v if \a v is contained in S,
	//! an invalid list iterator otherwise.
	NodeArray<ListIterator<node> > m_it;

	//! The list of nodes contained in \a S.
	ListPure<node> m_nodes;
};



//! Node sets.
/**
 * A node set maintains a subset \a S of the nodes contained in an associated
 * graph. This kind of node set provides efficient operations for testing
 * membership, insertion and deletion of elements, and clearing the set.
 *
 * In contrast to NodeSetPure, a NodeSet provides efficient access
 * to the number of elements stored in the set.
 *
 * \sa
 *   - NodeSetPure, NodeSetSimple
 *   - FaceSet, FaceSetPure, FaceSetSimple
 */
class OGDF_EXPORT NodeSet {
public:
	//! Creates an empty node set associated with graph \a G.
	NodeSet(const Graph &G) : m_it(G,ListIterator<node>()) { }

	// destructor
	~NodeSet() { }

	//! Inserts node \a v into \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	void insert(node v) {
		OGDF_ASSERT(v->graphOf() == m_it.graphOf());
		ListIterator<node> &itV = m_it[v];
		if (!itV.valid())
			itV = m_nodes.pushBack(v);
	}

	//! Removes node \a v from \a S.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	void remove(node v) {
		OGDF_ASSERT(v->graphOf() == m_it.graphOf());
		ListIterator<node> &itV = m_it[v];
		if (itV.valid()) {
			m_nodes.del(itV);
			itV = ListIterator<node>();
		}
	}


	//! Removes all nodes from \a S.
	/**
	 * After this operation, \a S is empty and still associated with the same graph.
	 * The runtime of this operations is O(k), where k is the number of nodes in \a S
	 * before this operation.
	 */
	void clear() {
		ListIterator<node> it;
		for(it = m_nodes.begin(); it.valid(); ++it) {
			m_it[*it] = ListIterator<node>();
		}
		m_nodes.clear();
	}


	//! Returns true if node \a v is contained in \a S, false otherwise.
	/**
	 * This operation has constant runtime.
	 *
	 * \pre \a v is a node in the associated graph.
	 */
	bool isMember(node v) const {
		OGDF_ASSERT(v->graphOf() == m_it.graphOf());
		return m_it[v].valid();
	}

	//! Returns the size of \a S.
	/**
	 * This operation has constant runtime.
	 */
	int size() const {
		return m_nodes.size();
	}

	//! Returns a reference to the list of nodes contained in \a S.
	/**
	 * This list can be used for iterating over all nodes in \a S.
	 */
	const List<node> &nodes() const {
		return m_nodes;
	}

	//! Copy constructor.
	NodeSet(const NodeSet& V) : m_it(*V.m_it.graphOf(), ListIterator<node>()) {
		forall_listiterators(node, it, V.m_nodes) {
			insert(*it);
		}
	}

	//! Assignment operator.
	NodeSet &operator=(const NodeSet &V) {
		m_nodes.clear();
		m_it.init(*V.m_it.graphOf());
		forall_listiterators(node, it, V.m_nodes) {
			insert(*it);
		}
		return *this;
	}

private:
	//! m_it[v] contains the list iterator pointing to \a v if \a v is contained in S,
	//! an invalid list iterator otherwise.
	NodeArray<ListIterator<node> > m_it;

	//! The list of nodes contained in \a S.
	List<node> m_nodes;
};


} // end namespace ogdf


#endif
