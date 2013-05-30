/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:39 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of Hierarchy class.
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

#ifndef OGDF_HIERARCHY_H
#define OGDF_HIERARCHY_H


#include <ogdf/basic/EdgeArray.h>
#include <ogdf/layered/Level.h>
#include <ogdf/basic/GraphCopy.h>


namespace ogdf {

	//! Representation of proper hierarchies used by Sugiyama-layout.
	/**
	 * \see Level, SugiyamaLayout
	 */
	class OGDF_EXPORT Hierarchy {

		friend class LayerBasedUPRLayout;

		GraphCopy m_GC; //!< The graph copy representing the topology of the proper hierarchy.
		NodeArray<int> m_rank; //!< The rank (level) of a node.
		Array<int> m_size;

	public:
		//! Creates an empty hierarchy.
		Hierarchy() { }
		//! Creates an hierarchy of graph \a G with node ranks \a rank.
		Hierarchy(const Graph &G, const NodeArray<int> &rank);

		// destruction
		~Hierarchy() { }

		void createEmpty(const Graph &G);
		void initByNodes(const List<node> &nodes,
			EdgeArray<edge> &eCopy,
			const NodeArray<int> &rank);

		//! Conversion to const GraphCopy reference.
		operator const GraphCopy &() const { return m_GC; }

		//! Returns the rank (level) of node \a v.
		int rank(node v) const { return m_rank[v]; }

		int maxRank() const { return m_size.high(); }

		int size(int i) const { return m_size[i]; }

		bool isLongEdgeDummy(node v) const {
			return (m_GC.isDummy(v) && v->outdeg() == 1);
		}

	private:
		void doInit(const NodeArray<int> &rank);
	};



	class OGDF_EXPORT HierarchyLevels {

	public:

		friend class Level;
		friend class LayerBasedUPRLayout;

		//! The traversing direction of the layer-by-layer sweep.
		enum TraversingDir { downward, upward };

	private:
		const Hierarchy &m_H;

		Array<Level *> m_pLevel; //!< The array of all levels.
		NodeArray<int> m_pos;    //!< The position of a node on its level.

		NodeArray<Array<node> > m_lowerAdjNodes; //!< (Sorted) adjacent nodes on lower level.
		NodeArray<Array<node> > m_upperAdjNodes; //!< (Sorted) adjacent nodes on upper level.

		NodeArray<int> m_nSet; //!< (Only used by buildAdjNodes().)

		TraversingDir m_direction; //!< The current direction of layer-by-layer sweep.

	public:
		explicit HierarchyLevels(const Hierarchy &H);
		~HierarchyLevels();

		const Hierarchy &hierarchy() const { return m_H; }

		//! Returns the current direction of layer-by-layer sweep.
		TraversingDir direction() const {
			return m_direction;
		}

		//! Sets the current direction of layer-by-layer sweep.
		void direction (TraversingDir dir) {
			m_direction = dir;
		}

		//! Returns the number of levels.
		int size() const { return m_pLevel.size(); }

		//! Returns the maximal array index of a level (= size()-1).
		int high() const { return m_pLevel.high(); }

		//! Returns the position of node \a v on its level.
		int pos(node v) const { return m_pos[v]; }

		//! Returns the adjacent nodes of \a v (according to direction()).
		const Array<node> &adjNodes(node v) const {
			return (m_direction == downward) ? m_lowerAdjNodes[v] :
			m_upperAdjNodes[v];
		}

		//! Returns the adjacent nodes of \a v.
		const Array<node> &adjNodes(node v, TraversingDir dir) const {
			return (dir == downward) ? m_lowerAdjNodes[v] :
			m_upperAdjNodes[v];
		}

		//! Returns the adjacent level of level \a i (according to direction()).
		const Level &adjLevel(int i) const {
			return (m_direction == downward) ? *m_pLevel[i-1] : *m_pLevel[i+1];
		}

		//! Returns the <i>i</i>-th level.
		const Level &operator[](int i) const { return *m_pLevel[i]; }

		//! Returns the <i>i</i>-th level.
		Level &operator[](int i) { return *m_pLevel[i]; }


		//! Computes the number of crossings between level \a i and \a i+1.
		int calculateCrossings(int i) const;
		//! Computes the total number of crossings.
		int calculateCrossings() const;

		//! Computes the number of crossings between level \a i and \a i+1 (for simultaneous drawing).
		int calculateCrossingsSimDraw(int i, const EdgeArray<__uint32> *edgeSubGraphs) const;
		//! Computes the total number of crossings (for simultaneous drawing).
		int calculateCrossingsSimDraw(const EdgeArray<__uint32> *edgeSubGraphs) const;

		//! Stores the position of nodes in \a oldPos.
		void storePos (NodeArray<int> &oldPos) const;
		//! Restores the position of nodes from \a newPos.
		void restorePos (const NodeArray<int> &newPos);

		//! Permutes the order of nodes on each level.
		void permute();

#ifdef OGDF_HAVE_CPP11
		template<class RNG>
		void permute(RNG &rng);
#endif

		//! Adjusts node positions such that nodes are ordered according to components numbers.
		void separateCCs(int numCC, const NodeArray<int> &component);

		bool transpose(node v);

		void print(ostream &os) const;

		void buildAdjNodes(int i);
		void buildAdjNodes();

		void check() const;

	private:
		int transposePart(const Array<node> &adjV, const Array<node> &adjW);

		OGDF_MALLOC_NEW_DELETE
	};


#ifdef OGDF_HAVE_CPP11
	template<class RNG>
	void HierarchyLevels::permute(RNG &rng)
	{
		for(int i = 0; i < m_pLevel.high(); ++i) {
			Level &level = *m_pLevel[i];
			level.m_nodes.permute(rng);
			for(int j = 0; j <= level.high(); ++j)
				m_pos[level[j]] = j;
		}

		buildAdjNodes();
	}
#endif


} // end namespace ogdf


#endif
