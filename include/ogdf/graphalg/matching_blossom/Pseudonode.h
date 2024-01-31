/** \file
 * \brief Helper structure representing a pseudonode in the Blossom algorithm.
 *
 * \author Joshua Sangmeister
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
#include <ogdf/graphalg/matching_blossom/Cycle.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <unordered_map>
#include <vector>

namespace ogdf {
namespace matching_blossom {

//! Helper class representing a pseudonode in the Blossom algorithm.
class OGDF_EXPORT Pseudonode {
	//! Helper class to store reference edges for all self loops.
	class ReferenceEdges {
	private:
		//! A mapping of all reference edges to all self loops that were removed because of them.
		std::unordered_map<edge, std::unordered_set<edge>> m_refToSelfLoops;

		//! A mapping of all self loops to their reference edges
		std::unordered_map<edge, edge> m_selfLoopToRef;

		//! A mapping of all self loops to the pseudonode they previously pointed to
		std::unordered_map<edge, Pseudonode*> m_selfLoopToOtherPseudonode;

	public:
		ValueIteratorContainer<edge, edge> refEdges;

		ReferenceEdges() : refEdges(m_selfLoopToRef) { }

		std::unordered_set<edge>& selfLoops(edge ref) { return m_refToSelfLoops[ref]; }

		//! Add a self loop \p selfLoop which was removed because of the reference edge \p ref and
		//! pointed previously to pseudonode \p other.
		void addReference(edge ref, edge selfLoop, Pseudonode* other) {
			m_refToSelfLoops[ref].insert(selfLoop);
			m_selfLoopToRef[selfLoop] = ref;
			m_selfLoopToOtherPseudonode[selfLoop] = other;
		}

		//! Remove the given \p selfLoop from the reference edges of the other pseudonode.
		void removeFromOther(edge selfLoop) {
			auto it = m_selfLoopToOtherPseudonode.find(selfLoop);
			if (it != m_selfLoopToOtherPseudonode.end() && it->second != nullptr) {
				auto& otherRefEdges = it->second->referenceEdges;
				auto edgeIt = otherRefEdges.m_selfLoopToRef.find(selfLoop);
				if (edgeIt != otherRefEdges.m_selfLoopToRef.end()) {
					otherRefEdges.m_refToSelfLoops[edgeIt->second].erase(selfLoop);
					otherRefEdges.m_selfLoopToRef.erase(edgeIt);
					otherRefEdges.m_selfLoopToOtherPseudonode.erase(selfLoop);
				}
			}
		}
	};

public:
	//! The node in the graph that this pseudonode is linked to.
	node graphNode;

	//! The odd-length cycle that this pseudonode represents
	Cycle* cycle;

	//! The ReferenceEdges for self loops of this pseudonode.
	ReferenceEdges referenceEdges;

	Pseudonode(node _graphNode, Cycle* _cycle);

	~Pseudonode();

	//! Add a self loop \p selfLoop which was removed because of the reference edge \p ref and
	//! pointed previously to pseudonode \p other.
	void addReference(edge ref, edge selfLoop, Pseudonode* other);
};

}
}
