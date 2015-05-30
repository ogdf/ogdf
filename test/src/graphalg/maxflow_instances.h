/** \file
 * \brief Tests-Instances for Max Flow Algorithms
 *
 * \author Tilo Wiedera
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/EdgeArray.h>

#ifndef OGDF_MAX_FLOW_INSTANCES_H
#define OGDF_MAX_FLOW_INSTANCES_H

const int MAX_FLOW_INSTANCES_COUNT = 5;

using namespace ogdf;
using std::string;

/**
 * Defines which properties a graph fullfills
 * or an algorithm requires.
 */
enum MaxFlowRequirement {
	MFR_NONE = 0,
	MFR_PLANAR = 1,
	MFR_ST_PLANAR = 2,
	MFR_CONNECTED = 4,
	MFR_ST_NON_INCIDENT_FACE = 8,
};

MaxFlowRequirement operator|(MaxFlowRequirement a, MaxFlowRequirement b)
{
	return MaxFlowRequirement(int(a) | int(b));
}

/**
 * Returns an instance for testing maximum flow algorithms.
 *
 * @param index the index of this instance (between 0 and MAX_FLOW_INSTANCES_COUNT)
 * @param name the human readable description of the graph
 * @param graph the actual graph
 * @param caps the capacity of each edge in the graph
 * @param props the properties of this graph as specified by MaxFlowRequirements
 * @param s the source node for the flow problem
 * @param t the sink node for the flow problem
 * @param flow the optimal flow value from source to sink
 */
template<typename T>
bool getMaxFlowInstance(
  int index,
  std::string &name,
  Graph &graph,
  EdgeArray<T> &caps,
  MaxFlowRequirement &props,
  node &s,
  node &t,
  T &flow)
{
	graph.clear();
	caps.init(graph, 0);
	s = graph.newNode();
	t = graph.newNode();
	props = MFR_NONE;
	name = "### UNKNOWN ###";
	bool result = true;

	if(index == 0) {
		name = "an easy planar instance";
		props = MFR_CONNECTED | MFR_PLANAR | MFR_ST_NON_INCIDENT_FACE;
		flow = 4;

		node a = graph.newNode();
		node b = graph.newNode();

		edge sb = graph.newEdge(s, b);
		edge sa = graph.newEdge(s, a);
		edge st = graph.newEdge(s, t);
		edge ab = graph.newEdge(a, b);
		edge at = graph.newEdge(a, t);
		edge bt = graph.newEdge(b, t);

		caps[sb] = 1;
		caps[sa] = 2;
		caps[st] = 1;
		caps[ab] = 1;
		caps[at] = 1;
		caps[bt] = 2;

	} else if(index == 1) {
		name = "a not so easy planar instance";
		props = MFR_CONNECTED | MFR_PLANAR | MFR_ST_NON_INCIDENT_FACE;
		flow = 19;

		node v1 = graph.newNode();
		node v2 = graph.newNode();
		node v3 = graph.newNode();
		node v4 = graph.newNode();
		node v5 = graph.newNode();
		node v6 = graph.newNode();
		node v7 = graph.newNode();
		node v8 = graph.newNode();
		node v9 = graph.newNode();
		node v10 = graph.newNode();
		node v11 = graph.newNode();
		node v12= graph.newNode();
		node v13 = graph.newNode();
		node v14 = graph.newNode();

		edge e1 = graph.newEdge(s, v1);
		edge e2 = graph.newEdge(s, v2);
		edge e3 = graph.newEdge(s, v3);
		edge e4 = graph.newEdge(v1, v4);
		edge e5 = graph.newEdge(v1, v6);
		edge e6 = graph.newEdge(v2, v7);
		edge e7 = graph.newEdge(v2, v10);
		edge e8 = graph.newEdge(v3, v4);
		edge e9 = graph.newEdge(v3, v10);
		edge e10 = graph.newEdge(v3, v12);
		edge e11 = graph.newEdge(v4, v3);
		edge e12 = graph.newEdge(v4, v13);
		edge e13 = graph.newEdge(v4, v5);
		edge e14 = graph.newEdge(v5, t);
		edge e15 = graph.newEdge(v6, v5);
		edge e16 = graph.newEdge(v6, v14);
		edge e17 = graph.newEdge(v6, v7);
		edge e18 = graph.newEdge(v7, v6);
		edge e19 = graph.newEdge(v7, v14);
		edge e20 = graph.newEdge(v7, v8);
		edge e21 = graph.newEdge(v8, t);
		edge e22 = graph.newEdge(v9, v8);
		edge e23 = graph.newEdge(v9, v11);
		edge e24 = graph.newEdge(v10, v9);
		edge e25 = graph.newEdge(v10, v11);
		edge e26 = graph.newEdge(v11, t);
		edge e27 = graph.newEdge(v12, v11);
		edge e28 = graph.newEdge(v12, v13);
		edge e29 = graph.newEdge(v13, v12);
		edge e30 = graph.newEdge(v13, t);
		edge e31 = graph.newEdge(v14, t);

		caps[e1] = 10;
		caps[e2] = 6;
		caps[e3] = 9;
		caps[e4] = 5;
		caps[e5] = 4;
		caps[e6] = 4;
		caps[e7] = 4;
		caps[e8] = 3;
		caps[e9] = 3;
		caps[e10] = 3;
		caps[e11] = 3;
		caps[e12] = 2;
		caps[e13] = 5;
		caps[e14] = 3;
		caps[e15] = 2;
		caps[e16] = 3;
		caps[e17] = 1;
		caps[e18] = 3;
		caps[e19] = 3;
		caps[e20] = 8;
		caps[e21] = 3;
		caps[e22] = 2;
		caps[e23] = 6;
		caps[e24] = 6;
		caps[e25] = 1;
		caps[e26] = 6;
		caps[e27] = 2;
		caps[e28] = 2;
		caps[e29] = 3;
		caps[e30] = 6;
		caps[e31] = 3;

	} else if(index == 2) {
		name = "a little bit harder planar instance";
		props = MFR_CONNECTED | MFR_PLANAR | MFR_ST_NON_INCIDENT_FACE;
		flow = 6;

		// +------------+
		// | +------+   |
		// | |      |   |
		// | |  s - a - b
		// | |  | \ |   |
		// | +- c - d - e
		// |    |   | \ |
		// +--- f - g - t

		node a = graph.newNode();
		node b = graph.newNode();
		node c = graph.newNode();
		node d = graph.newNode();
		node e = graph.newNode();
		node f = graph.newNode();
		node g = graph.newNode();

		caps.init(graph, 1);
		graph.newEdge(s, a);
		caps[graph.newEdge(s, c)] = 4;
		graph.newEdge(s, d);
		caps[graph.newEdge(a, b)] = 2;
		graph.newEdge(a, d);
		graph.newEdge(b, e);
		graph.newEdge(b, f);
		caps[graph.newEdge(c, a)] = 2;
		graph.newEdge(c, d);
		graph.newEdge(c, f);
		graph.newEdge(d, e);
		graph.newEdge(d, g);
		graph.newEdge(d, t);
		caps[graph.newEdge(e, t)] = 2;
		caps[graph.newEdge(f, g)] = 2;
		caps[graph.newEdge(g, t)] = 3;

	} else if(index == 3) {
		name = "a tiny disconnected graph";
		props = MFR_PLANAR | MFR_ST_PLANAR;
		flow = 0;

		node v = graph.newNode();
		s = graph.newNode();
		t = graph.newNode();

		edge e_0_2 = graph.newEdge(v, t);
		caps[e_0_2] = 0;
		edge e_1_0 = graph.newEdge(s, v);
		caps[e_1_0] = 3;
		edge e_0_1 = graph.newEdge(v, s);
		caps[e_0_1] = 2;

	} else if(index == 4) {
		name = "an easy s-t-planar instance";
		props = MFR_CONNECTED | MFR_PLANAR | MFR_ST_PLANAR;
		flow = 3;

		node a = graph.newNode();
		node b = graph.newNode();

		edge sb = graph.newEdge(s, b);
		edge sa = graph.newEdge(s, a);
		edge st = graph.newEdge(s, t);
		edge at = graph.newEdge(a, t);
		edge bt = graph.newEdge(b, t);

		caps[sb] = 1;
		caps[sa] = 2;
		caps[st] = 1;
		caps[at] = 1;
		caps[bt] = 2;

	} else {
		result = false;
		cerr << "invalid call of getMaxFlowInstance: Index is " << index << std::endl;
	}
	return result;
}

#endif
