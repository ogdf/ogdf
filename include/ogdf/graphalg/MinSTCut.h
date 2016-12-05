/** \file
 * \brief Class for converting a maximum flow to a minimum cut.
 *
 * \author Tilo Wiedera
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

#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/EpsilonTest.h>
#include <memory>

namespace ogdf {

/**
 * Class for converting a maximum flow into a minimum cut.
 */
template<typename T>
class MinSTCut
{
public:
	/**
	 * The three types of cuts.
	 */
	enum cutType {
		FRONT_CUT,
		BACK_CUT,
		NO_CUT
	};

private:
	std::unique_ptr<EpsilonTest> m_et; // the module used for epsilon tests
	NodeArray<cutType> m_nodeSet; // holds the partition type for each node
	int m_frontCutCount; // the number of nodes in the front cut
	int m_backCutCount; // the number of nodes in the back cut
	int m_totalCount; // the total number of nodes in the graph

public:
	/**
	 * Creates a new MinSTCut instance.
	 */
	MinSTCut() : m_et(new EpsilonTest())
	{}

	/**
	 * Assigns a new epsilon test.
	 */
	void setEpsilonTest(EpsilonTest *et)
	{
		m_et.reset(et);
	}

	/**
	 * Partitions the nodes to front- and backcut.
	 *
	 * \param weights
	 *	the weights (aka capacity) of the edges
	 * \param flow
	 * 	the precomputed flow for each edge
	 * \param source
	 *  the source of the minimum cut
	 * \param target
	 * 	the target (aka sink) of the minimum cut
	 */
	void call(const EdgeArray<T> &weights, const EdgeArray<T> &flow, const node source, const node target)
	{
		m_frontCutCount = 0;
		m_backCutCount = 0;
		m_totalCount = weights.graphOf()->numberOfNodes();

		List<node> queue;
		m_nodeSet.init(*weights.graphOf(), NO_CUT);
		queue.pushBack(source);
		m_nodeSet[source] = FRONT_CUT;

		// front cut
		while(!queue.empty()) {
			const node v = queue.popFrontRet();
			for(adjEntry adj : v->adjEntries) {
				const node w = adj->twinNode();
				const edge e = adj->theEdge();
				if(m_nodeSet[w] == NO_CUT
				  && ((e->source() == v && m_et->less(flow[e], weights[e]))
				    || (e->target() == v && m_et->greater(flow[e], (T) 0)))) {
					queue.pushBack(w);
					m_nodeSet[w] = FRONT_CUT;
					m_frontCutCount++;
				}
			}
		}

		// back cut
		queue.pushBack(target);
		m_nodeSet[target] = BACK_CUT;

		while(!queue.empty()) {
			const node v = queue.popFrontRet();
			for(adjEntry adj : v->adjEntries) {
				const node w = adj->twinNode();
				const edge e = adj->theEdge();
				if(m_nodeSet[w] == NO_CUT
				  && ((e->target() == v && m_et->less(flow[e], weights[e]))
				    || (e->source() == v && m_et->greater(flow[e], (T) 0)))) {
					queue.pushBack(w);
					m_nodeSet[w] = BACK_CUT;
					m_backCutCount++;
				}
			}
		}
	}

	/**
	 * Returns whether the front cut is the complement of
	 * the backcut. i.e. there are no nodes not assigned to
	 * one of both cut types.
	 */
	bool frontCutIsComplementOfBackCut() const
	{
		return m_backCutCount + m_frontCutCount == m_totalCount;
	}

	/**
	 * Returns whether this edge is entering the back cut.
	 */
	bool isFrontCutEdge(const edge e) const
	{
		return m_nodeSet[e->source()] == FRONT_CUT
		  && m_nodeSet[e->target()] != FRONT_CUT;
	}

	/**
	 * Returns whether this edge is leaving the front cut.
	 */
	bool isBackCutEdge(const edge e) const
	{
		return m_nodeSet[e->target()] == BACK_CUT
		  && m_nodeSet[e->source()] != BACK_CUT;
	}

	/**
	 * Returns whether this node is part of the front cut.
	 * Meaning it is located in the same set as the source.
	 */
	bool isInFrontCut(const node v) const
	{
		return m_nodeSet[v] == FRONT_CUT;
	}

	/*
	 * Returns whether this node is part of the back cut.
	 * Meaning it is located in the same set as the target.
	 */
	bool isInBackCut(const node v)  const
	{
		return m_nodeSet[v] == BACK_CUT;
	}

	/**
	 * Return whether this node is of the specified type.
	 *
	 * \param v
	 * 	the node to be tested
	 * \param type
	 * 	the cut type to test for (see MinSTCut::CutType)
	 *
	 * \return
	 * 	true if the node is contained in the specified cut
	 */
	bool isOfType(const node v, cutType type) const
	{
		return m_nodeSet[v] == type;
	}
};

}
