/** \file
 * \brief Tests for the basic dual graph class
 *
 * \author Mirko Wagner
 *
 * \par License:
 * This file is part of the Open myGraPath Halving, Drawing Framework (OGDF).
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
 * Software Foundation, Inc., 51 FRank,lin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <bandit/bandit.h>

#include <ogdf/basic/DualGraph.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <resources.h>
#include <ogdf/basic/Graph_d.h>

using namespace ogdf;
using namespace bandit;


go_bandit([](){
describe("DualGraph",[&](){
	Graph graph;
	CombinatorialEmbedding *combEmb;
	DualGraph *dualGraph;
	node u,v,w;
	edge e,f,g;

	before_each([&](){
		graph = Graph();
		u = graph.newNode();
		v = graph.newNode();
		w = graph.newNode();
		e = graph.newEdge(u, v);
		f = graph.newEdge(v, w);
		g = graph.newEdge(w, u);
		combEmb = new CombinatorialEmbedding(graph);
		dualGraph = new DualGraph(*combEmb);
	});

	after_each([&](){
		delete dualGraph;
		delete combEmb;
	});

	it("knows its primal embedding",[&](){
		AssertThat(&dualGraph->getPrimalEmbedding(), Equals(combEmb));
	});

	it("knows its primal graph",[&](){
		AssertThat(&dualGraph->getPrimalGraph(), Equals(&graph));
	});

	it("knows the dual edge to a primal one and the other way round",[&](){
		delete dualGraph;
		delete combEmb;
		graph = Graph();
		planarTriconnectedGraph(graph, 42, 42*3 - 6);
		combEmb = new CombinatorialEmbedding(graph);
		dualGraph = new DualGraph(*combEmb);
		e = graph.chooseEdge();
		AssertThat(dualGraph->primalEdge(dualGraph->dualEdge(e)), Equals(e));
		AssertThat(dualGraph->dualEdge(e), !Equals(e));
		AssertThat(dualGraph->primalFace(dualGraph->dualEdge(e)->source()), Equals(combEmb->rightFace(e->adjSource())) || Equals(combEmb->rightFace(e->adjTarget())));
	});

	it("knows the dual node to a primal face and the other way round",[&](){
		face a = combEmb->chooseFace();
		AssertThat(dualGraph->primalFace(dualGraph->dualNode(a)), Equals(a));
	});

	it("knows the dual face to a primal node and the other way round",[&](){
		AssertThat(dualGraph->primalNode(dualGraph->dualFace(v)), Equals(v));
	});
});
});
