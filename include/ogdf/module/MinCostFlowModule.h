/** \file
 * \brief Template of base class of min-cost-flow algorithms.
 *
 * Includes some useful functions dealing with min-cost flow
 * (problem generator, problem checker).
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

#ifndef OGDF_MIN_COST_FLOW_MODULE_H
#define OGDF_MIN_COST_FLOW_MODULE_H


#include <ogdf/basic/Graph.h>


namespace ogdf {


/**
 * \brief Interface for min-cost flow algorithms.
 */
template<typename TCost>
class OGDF_EXPORT MinCostFlowModule
{
public:
	//! Initializes a min-cost flow module.
	MinCostFlowModule() { }

	// destruction
	virtual ~MinCostFlowModule() { }

	/**
	* \brief Computes a min-cost flow in the directed graph \a G.
	*
	* \pre \a G must be connected, \a lowerBound[\a e] \f$\leq\f$ \a upperBound[\a e]
	*      for all edges \a e, and the sum over all supplies must be zero.
	*
	* @param G is the directed input graph.
	* @param lowerBound gives the lower bound for the flow on each edge.
	* @param upperBound gives the upper bound for the flow on each edge.
	* @param cost gives the costs for each edge.
	* @param supply gives the supply (or demand if negative) of each node.
	* @param flow is assigned the computed flow on each edge.
	* \return true iff a feasible min-cost flow exists.
	*/
	virtual bool call(
		const Graph &G,                   // directed graph
		const EdgeArray<int> &lowerBound, // lower bound for flow
		const EdgeArray<int> &upperBound, // upper bound for flow
		const EdgeArray<TCost> &cost,     // cost of an edge
		const NodeArray<int> &supply,     // supply (if neg. demand) of a node
		EdgeArray<int> &flow)             // computed flow
	{
		NodeArray<TCost> dual(G);
		return call(G, lowerBound, upperBound, cost, supply, flow, dual);
	}

	/**
	 * \brief Computes a min-cost flow in the directed graph \a G.
	 *
	 * \pre \a G must be connected, \a lowerBound[\a e] \f$\leq\f$ \a upperBound[\a e]
	 *      for all edges \a e, and the sum over all supplies must be zero.
	 *
	 * @param G is the directed input graph.
	 * @param lowerBound gives the lower bound for the flow on each edge.
	 * @param upperBound gives the upper bound for the flow on each edge.
	 * @param cost gives the costs for each edge.
	 * @param supply gives the supply (or demand if negative) of each node.
	 * @param flow is assigned the computed flow on each edge.
	 * @param dual is assigned the computed dual variables.
	 * \return true iff a feasible min-cost flow exists.
	 */
	virtual bool call(
		const Graph &G,                   // directed graph
		const EdgeArray<int> &lowerBound, // lower bound for flow
		const EdgeArray<int> &upperBound, // upper bound for flow
		const EdgeArray<TCost> &cost,     // cost of an edge
		const NodeArray<int> &supply,     // supply (if neg. demand) of a node
		EdgeArray<int> &flow,			  // computed flow
		NodeArray<TCost> &dual            // computed dual variables
		) = 0;


	//
	// static functions
	//

	/**
	 * \brief Generates an instance of a min-cost flow problem with \a n nodes and
	 *        \a m+\a n edges.
	 */
	static void generateProblem(
		Graph &G,
		int n,
		int m,
		EdgeArray<int> &lowerBound,
		EdgeArray<int> &upperBound,
		EdgeArray<TCost> &cost,
		NodeArray<int> &supply);


	/**
	 * \brief Checks if a given min-cost flow problem instance satisfies
	 *        the preconditions.
	 * The following preconditions are checked:
	 *   - \a lowerBound[\a e] \f$\leq\f$ \a upperBound[\a e] for all edges \a e
	 *   - sum over all \a supply[\a v] = 0
	 *
	 * @param G is the input graph.
	 * @param lowerBound gives the lower bound for the flow on each edge.
	 * @param upperBound gives the upper bound for the flow on each edge.
	 * @param supply gives the supply (or demand if negative) of each node.
	 * \return true iff the problem satisfies the preconditions.
	 */
	static bool checkProblem(
		const Graph &G,
		const EdgeArray<int> &lowerBound,
		const EdgeArray<int> &upperBound,
		const NodeArray<int> &supply);



	/**
	 * \brief checks if a computed flow is a feasible solution to the given problem
	 *        instance.
	 *
	 * Checks in particular if:
	 *   - \a lowerBound[\a e] \f$\leq\f$ \a flow[\a e] \f$\leq\f$ \a upperBound[\a e]
	 *   - sum \a flow[\a e], \a e is outgoing edge of \a v minus
	 *     sum \a flow[\a e], \a e is incoming edge of \a v equals \a supply[\a v]
	 *     for each node \a v
	 *
	 * @param G is the input graph.
	 * @param lowerBound gives the lower bound for the flow on each edge.
	 * @param upperBound gives the upper bound for the flow on each edge.
	 * @param cost gives the costs for each edge.
	 * @param supply gives the supply (or demand if negative) of each node.
	 * @param flow is the flow on each edge.
	 * @param value is assigned the value of the flow.
	 * \return true iff the solution is feasible.
	 */
	static bool checkComputedFlow(
		const Graph &G,
		EdgeArray<int> &lowerBound,
		EdgeArray<int> &upperBound,
		EdgeArray<TCost> &cost,
		NodeArray<int> &supply,
		EdgeArray<int> &flow,
		TCost &value);

	/**
	 * \brief checks if a computed flow is a feasible solution to the given problem
	 *        instance.
	 *
	 * Checks in particular if:
	 *   - \a lowerBound[\a e] \f$\leq\f$ \a flow[\a e] \f$\leq\f$ \a upperBound[\a e]
	 *   - sum \a flow[\a e], \a e is outgoing edge of \a v minus
	 *     sum \a flow[\a e], \a e is incoming edge of \a v equals \a supply[\a v]
	 *     for each node \a v
	 *
	 * @param G is the input graph.
	 * @param lowerBound gives the lower bound for the flow on each edge.
	 * @param upperBound gives the upper bound for the flow on each edge.
	 * @param cost gives the costs for each edge.
	 * @param supply gives the supply (or demand if negative) of each node.
	 * @param flow is the flow on each edge.
	 * \return true iff the solution is feasible.
	 */
	static bool checkComputedFlow(
		const Graph &G,
		EdgeArray<int> &lowerBound,
		EdgeArray<int> &upperBound,
		EdgeArray<TCost> &cost,
		NodeArray<int> &supply,
		EdgeArray<int> &flow)
	{
		int value;
		return checkComputedFlow(
			G,lowerBound,upperBound,cost,supply,flow,value);
	}
};


} // end namespace ogdf


// Implementation

#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>

namespace ogdf {

template<typename TCost>
void MinCostFlowModule<TCost>::generateProblem(
	Graph &G,
	int n,
	int m,
	EdgeArray<int> &lowerBound,
	EdgeArray<int> &upperBound,
	EdgeArray<TCost> &cost,
	NodeArray<int> &supply)
{
	ogdf::randomGraph(G,n,m);

	node s = G.firstNode();
	node t = G.lastNode();

	for(node v : G.nodes) {
		G.newEdge(s,v);
		G.newEdge(v,t);
	}

	for(edge e : G.edges) {
		lowerBound[e] = 0;
		upperBound[e] = (e->source() != s) ? ogdf::randomNumber(1,10) : ogdf::randomNumber(2,13);
		cost[e] = static_cast<TCost>(ogdf::randomNumber(0,100));
	}


	for(node v = G.firstNode(), vl = G.lastNode(); true; v = v->succ(), vl = vl->pred()) {
		if (v == vl) {
			supply[v] = 0;
			break;
		}

		supply[v] = -(supply[vl] = ogdf::randomNumber(-1,1));

		if (vl == v->succ())
			break;
	}

}

template<typename TCost>
bool MinCostFlowModule<TCost>::checkProblem(
	const Graph &G,
	const EdgeArray<int> &lowerBound,
	const EdgeArray<int> &upperBound,
	const NodeArray<int> &supply)
{
	if(isConnected(G) == false)
		return false;

	for(edge e : G.edges) {
		if (lowerBound[e] > upperBound[e])
			return false;
	}

	int sum = 0;
	for(node v : G.nodes) {
		sum += supply[v];
	}

	return (sum == 0);
}


template<typename TCost>
bool MinCostFlowModule<TCost>::checkComputedFlow(
	const Graph &G,
	EdgeArray<int> &lowerBound,
	EdgeArray<int> &upperBound,
	EdgeArray<TCost> &cost,
	NodeArray<int> &supply,
	EdgeArray<int> &flow,
	TCost &value)
{
	value = 0;

	for (edge e : G.edges) {
		if (flow[e] < lowerBound[e] || upperBound[e] < flow[e]) {
			return false;
		}

		value += flow[e] * cost[e];
	}

	for (node v : G.nodes) {
		int sum = 0;
		edge e;
		forall_adj_edges(e,v) {
			if (e->isSelfLoop()) {
				continue;
			}

			if (e->source() == v) {
				sum += flow[e];
			} else {
				sum -= flow[e];
			}
		}
		if (sum != supply[v]) {
			return false;
		}
	}

	return true;
}


} // end namespace ogdf


#endif
