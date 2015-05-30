/** \file
 * \brief Interface for Max Flow Algorithms
 *
 * \author Ivo Hedtke
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

#ifndef OGDF_MAX_FLOW_MODULE_H
#define OGDF_MAX_FLOW_MODULE_H

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/EpsilonTest.h>

namespace ogdf {

template <typename T>
class OGDF_EXPORT MaxFlowModule
{
protected:
	EpsilonTest *m_et; //!< Pointer to the used EpsilonTest.
	EdgeArray<T> *m_flow; //!< Pointer to (extern) flow array.
	const Graph *m_G; //!< Pointer to the given graph.
	const EdgeArray<T> *m_cap; //!< Pointer to the given capacity array.
	const node *m_s; //!< Pointer to the source node.
	const node *m_t; //!< Pointer to the sink node.

private:
	bool usingExternFlow = false; //!< Is an extern flow array given in the constructor?
	bool doingAReInit = false; //!< Is the next "init" call a re-init?

	void destroy()
	{
		if (!usingExternFlow) {
			delete m_flow;
		}
		delete m_et;
	}

public:
	//! Empty Constructor.
	MaxFlowModule()
		: m_et(nullptr), m_flow(nullptr), m_G(nullptr), m_s(nullptr), m_t(nullptr) { }

	//! Constructor that calls init.
	/**
	 * @param graph is the graph for the flow problem.
	 * @param flow is an optional argument that can be used to force the
	 * algorithm to work on an user given "external" EdgeArray \a flow
	 */
	MaxFlowModule(const Graph &graph, EdgeArray<T> *flow = nullptr)
		: m_s(nullptr), m_t(nullptr)
	{
		init(graph, flow);
	}

	//! Destructor that deletes m_flow if it is an internal flow array.
	virtual ~MaxFlowModule()
	{
		destroy();
	}

	//! Initialize the problem with a graph and optional flow array.
	//! If no \a flow array is given, a new ("internal") array will be created.
	//! If a \a flow array is giveb, the algorithm uses this "external" array.
	/**
	 * @param graph is the graph for the flow problem.
	 * @param flow is an optional argument that can be used to force the
	 * algorithm to work on an user given "external" EdgeArray \a flow
	 */
	virtual void init(const Graph &graph, EdgeArray<T> *flow = nullptr)
	{
		// if re-init after an init with internal flow:
		if (doingAReInit) {
			destroy();
		}
		m_G = &graph;
		if (flow) {
			m_flow = flow;
			usingExternFlow = true;
		}
		else {
			usingExternFlow = false; // in case of a re-init
			m_flow = new EdgeArray<T>(*m_G, 0);
		}
		m_et = new EpsilonTest();
		doingAReInit = true;
	}

	//! Change the used EpsilonTest from StandardEpsilonTest to a user given
	//! EpsilonTest.
	/**
	 * @param eps use an EpsilonTest with epsilon
	 */
	void useEpsilonTest(const double &eps)
	{
		delete m_et;
		m_et = new EpsilonTest(eps);
	}

	//! Compute only the value of the flow. There algorithms with 2 phases where
	//! where the value of the flow is known after the first phase, but the flow
	//! itself is not feasible or not known at this time.
	//! If source and target are the same node, the algorithm must return 0.
	/**
	 * @return The value of the flow.
	 * @param cap is the EdgeArray of capacities.
	 * @param s is the source.
	 * @param t is the sink.
	 */
	virtual T computeValue(const EdgeArray<T> &cap, const node &s, const node &t) = 0;

	//! Compute the flow itself after the flow value is already computed. Only
	//! used in algorithms with 2 phases. The flow is stored in the array that
	//! the user gave in the constructor.
	virtual void computeFlowAfterValue() = 0;

	//! Compute the flow itself after the flow value is already computed. Only
	//! used in algorithms with 2 phases. The flow is stored in the parameter
	//! \a flow.
	/**
	 * @param flow The "internal" flow array is given in \a flow.
	 */
	void computeFlowAfterValue(EdgeArray<T> &flow)
	{
		computeFlowAfterValue();
		flow = *m_flow;
	}

	//! Only a shortcut for computeValue and computeFlowAfterValue.
	/**
	 * @return The value of the flow.
	 * @param cap is the EdgeArray of capacities.
	 * @param s is the source.
	 * @param t is the sink.
	 * @param flow A copy of the "internal" flow array is given in \a flow.
	 */
	T computeFlow(EdgeArray<T> &cap, node &s, node &t, EdgeArray<T> &flow)
	{
		T val = computeValue(cap, s, t);
		computeFlowAfterValue(flow);
		return val;
	}

};

}// namespace ogdf

#endif //OGDF_MAX_FLOW_MODULE_H
