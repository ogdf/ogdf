/*
 * $Revision: 3832 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 11:16:27 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration and implementation of Level class
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

#ifndef OGDF_LEVEL_H
#define OGDF_LEVEL_H



#include <ogdf/basic/Graph.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/tuples.h>
#include <ogdf/layered/CrossingMinInterfaces.h>

namespace ogdf {

class HierarchyLevels;
class LayerBasedUPRLayout;


template <class T = double>
class WeightComparer {
	const NodeArray<T> *m_pWeight;

public:
	WeightComparer(const NodeArray<T> *pWeight) : m_pWeight(pWeight) { }

	bool less(node v, node w) const { return (*m_pWeight)[v] < (*m_pWeight)[w]; }
	bool operator()(node v, node w) const { return (*m_pWeight)[v] < (*m_pWeight)[w]; }
};


//! Representation of levels in hierarchies.
/**
 * \see Hierarchy, SugiyamaLayout
 */
class OGDF_EXPORT Level : public LevelBase {

	friend class HierarchyLevels;
	friend class HierarchyLayoutModule;
	friend class LayerBasedUPRLayout;

	Array<node> m_nodes;     //!< The nodes on this level.
	HierarchyLevels *m_pLevels; //!< The hierarchy to which this level belongs.
	int m_index;             //!< The index of this level.

public:
	//! Creates a level with index \a index in hierarchy \a pHierarchy.
	/**
	 * @param pLevels is a pointer to the hierarchy to which the created level will belong.
	 * @param index   is the index of the level.
	 * @param num     is the number of nodes on this level.
	 */
	Level(HierarchyLevels *pLevels, int index, int num) :
		m_nodes(num), m_pLevels(pLevels), m_index(index) { }

	// destruction
	~Level() { }

	//! Returns the node at position \a i.
	const node &operator[](int i) const { return m_nodes[i]; }
	//! Returns the node at position \a i.
	node &operator[](int i) { return m_nodes[i]; }

	//! Returns the number of nodes on this level.
	int size() const { return m_nodes.size(); }

	//! Returns the maximal array index (= size()-1).
	int high() const { return m_nodes.high(); }

	//! Returns the array index of this level in the hierarchy.
	int index() const { return m_index; }

	//! Returns the (sorted) array of adjacent nodes of \a v (according to direction()).
	const Array<node> &adjNodes(node v) const;

	//! Returns the hierarchy to which this level belongs.
	const HierarchyLevels &levels() const { return *m_pLevels; }

	//! Exchanges nodes at position \a i and \a j.
	void swap(int i, int j);

	//! Sorts the nodes according to \a weight using quicksort.
	void sort(NodeArray<double> &weight);

	//! Sorts the nodes according to \a weight using bucket sort.
	void sort(NodeArray<int> &weight, int minBucket, int maxBucket);

	//! Sorts the nodes according to \a weight (without special placement for "isolated" nodes).
	void sortByWeightOnly(NodeArray<double> &weight);

	//!Sorts the nodes according to \a orderComparer
	template<class C>
	void sortOrder(C &orderComparer) {
		m_nodes.quicksort(orderComparer);
		recalcPos();
	}

	void recalcPos();

	friend ostream &operator<<(ostream &os, const Level &L) {
		os << L.m_nodes;
		return os;
	}

private:
	void getIsolatedNodes(SListPure<Tuple2<node,int> > &isolated) const;
	void setIsolatedNodes(SListPure<Tuple2<node,int> > &isolated);

	OGDF_MALLOC_NEW_DELETE
};

} // end namespace ogdf


#endif
