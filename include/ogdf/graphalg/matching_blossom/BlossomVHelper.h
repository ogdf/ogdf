/** \file
 * \brief Utility class for the Blossom V algorithm.
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

#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/matching_blossom/AuxGraph.h>
#include <ogdf/graphalg/matching_blossom/BlossomHelper.h>
#include <ogdf/graphalg/matching_blossom/PQ.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

namespace ogdf {
namespace matching_blossom {

template<class TWeight>
class AuxGraph;

/**
 *  Helper class for the blossom matching algorithms.
 */
template<class TWeight>
class BlossomVHelper : public BlossomHelper<TWeight> {
	using BlossomHelper<TWeight>::m_y;
	using BlossomHelper<TWeight>::m_eps;

	//! Helper function to get the top element of any priority queue.
	template<class E>
	E getTopElement(BlossomPQ<E, TWeight>& pq) {
		if (!pq.empty()) {
			return pq.topElement();
		}
		return nullptr;
	}

protected:
	AuxGraph<TWeight>* m_auxGraph;

public:
	using BlossomHelper<TWeight>::c;
	using BlossomHelper<TWeight>::isZero;

	//! The current iteration of the algorithm.
	long currentIteration;

	/**
	 * @brief Construct a new Blossom V Helper object
	 *
	 * @param greedyInit whether or not to use the greedy initialization
	 */
	BlossomVHelper(bool greedyInit) : BlossomHelper<TWeight>(greedyInit) { }

	//! Initialize the helper with the given data.
	template<class WeightContainer>
	bool init(const Graph& graph, const WeightContainer& weights, AuxGraph<TWeight>* auxGraph) {
		m_auxGraph = auxGraph;
		if (!BlossomHelper<TWeight>::init(graph, weights)) {
			return false;
		}
		m_auxGraph->reset();
		currentIteration = 0;
		return true;
	}

	//! Gets the actual reduced weight of the edge, taking into account the delta of the
	//! corresponding trees.
	TWeight getRealReducedWeight(edge e) override {
		return c(e) - realY(e->source()) - realY(e->target());
	}

	//! Returns the actual y value of the node, taking into account the delta of the corresponding
	//! tree.
	TWeight realY(node v) {
		auto auxNode = m_auxGraph->treeAuxNode(v);
		return m_y[v] + (auxNode ? auxNode->delta(v) : 0);
	}

	//! Returns the realReducedWeight of \p e.
	TWeight realValue(edge e) { return getRealReducedWeight(e); }

	//! Returns the realY of \p v.
	TWeight realValue(node v) { return realY(v); }

	//! Checks if the realValue of \p b is zero.
	bool isZeroCostNode(node v) { return isZero(realY(v)); }

	//! Checks if the top element of the given \p pq has realValue == 0 and returns it, if
	//! possible.
	template<class E>
	E getTopEligibleElement(BlossomPQ<E, TWeight>& pq) {
		if (!pq.empty()) {
			E el = pq.topElement();
			if (isZero(realValue(el))) {
				return el;
			}
		}
		return nullptr;
	}

	//! Returns the real value (realReducedWeight or realY) of the top element in the priority
	//! queue, or infinity if the queue is empty.
	template<class E>
	TWeight getRealTopPriority(BlossomPQ<E, TWeight>& pq) {
		auto el = getTopElement(pq);
		if (el == nullptr) {
			return infinity<TWeight>();
		} else {
			return realValue(el);
		}
	}
};

}
}
