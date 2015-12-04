/** \file
 * \brief Tests for the basic graph copy class
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

#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <resources.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Graph_d.h>

using namespace ogdf;
using namespace bandit;


go_bandit([](){
describe("Combinatorial Embedding", [&](){
	const int numberOfNodes(42);
	const int numberOfEdges(82);

	CombinatorialEmbedding *combEmb;
	const CombinatorialEmbedding *cCombEmb;
	Graph *graph;

	before_each([&](){
		graph = new Graph();
		planarTriconnectedGraph(*graph, numberOfNodes, numberOfEdges);
		combEmb = new CombinatorialEmbedding(*graph);
		cCombEmb = new CombinatorialEmbedding(*graph);
	});

	after_each([&](){
		if(combEmb) { delete combEmb; }
		if(cCombEmb) { delete cCombEmb; }
		delete graph;
	});

	it("knows its graph",[&](){
		AssertThat(&combEmb->getGraph(), Equals(graph));
		AssertThat(&cCombEmb->getGraph(), Equals(graph));
		AssertThat(&((Graph&)(*combEmb)), Equals(graph));
		AssertThat(&((const Graph&)(*combEmb)), Equals(graph));
	});

	describe("init",[&](){
		it("initializes w/o a graph",[&](){
			CombinatorialEmbedding *emptyCombEmb;
			emptyCombEmb = new CombinatorialEmbedding();
			AssertThat(&emptyCombEmb->getGraph(), IsNull());
			delete emptyCombEmb;
		});

		it("initializes w/o a graph const",[&](){
			ConstCombinatorialEmbedding *emptyCombEmb;
			emptyCombEmb = new ConstCombinatorialEmbedding();
			AssertThat(&emptyCombEmb->getGraph(), IsNull());
			delete emptyCombEmb;
		});

		it("initializes w a graph and the constructor",[&](){
			AssertThat(&combEmb->getGraph(), Equals(graph));
		});

		it("initializes w a graph and the constconstructor",[&](){
			ConstCombinatorialEmbedding *constCombEmb;
			constCombEmb = new ConstCombinatorialEmbedding(*graph);
			AssertThat(&constCombEmb->getGraph(), Equals(graph));
			delete constCombEmb;
		});

		it("won't initialize w a non-embedded graph and the constructor",[&](){
			delete graph;
			graph = new Graph();
			completeGraph(*graph, numberOfNodes);
			AssertThrows(PreconditionViolatedException, new CombinatorialEmbedding(*graph));
		});

		it("won't initialize w a non-embedded graph and the constconstructor",[&](){
			delete graph;
			graph = new Graph();
			completeGraph(*graph, numberOfNodes);
			AssertThrows(PreconditionViolatedException, new ConstCombinatorialEmbedding(*graph));
		});

		it("initializes w a graph and init()",[&](){
			combEmb->clear();
			combEmb->init(*graph);
			AssertThat(&combEmb->getGraph(), Equals(graph));
		});

		it("initializes w a graph and init() const",[&](){
			ConstCombinatorialEmbedding *constCombEmb;
			constCombEmb = new ConstCombinatorialEmbedding();
			constCombEmb->init(*graph);
			AssertThat(&constCombEmb->getGraph(), Equals(graph));
			delete constCombEmb;
		});

		it("won't initialize w a non-embedded graph and init()",[&](){
			delete graph;
			graph = new Graph();
			completeGraph(*graph, numberOfNodes);
			delete combEmb;
			combEmb = new CombinatorialEmbedding();
			AssertThrows(PreconditionViolatedException, combEmb->init(*graph));
		});

		it("won't initialize w a non-embedded graph and init() const",[&](){
			delete graph;
			graph = new Graph();
			completeGraph(*graph, numberOfNodes);
			ConstCombinatorialEmbedding *constCombEmb;
			constCombEmb = new ConstCombinatorialEmbedding();
			AssertThrows(PreconditionViolatedException, constCombEmb->init(*graph));
			delete constCombEmb;
		});
	});

	it("clears",[&](){
		combEmb->clear();
		AssertThat(graph->numberOfNodes(), Equals(0));
		AssertThat(graph->numberOfEdges(), Equals(0));
		AssertThat(combEmb->numberOfFaces(), Equals(0));
	});

	describe("update",[&](){
		it("splits an edge",[&](){
			int numberOfFaces = combEmb->numberOfFaces();
			edge splitEdgeBeginning = graph->chooseEdge();
			face leftFace = combEmb->leftFace(splitEdgeBeginning->adjSource());
			int leftFaceSize = leftFace->size();
			face rightFace = combEmb->rightFace(splitEdgeBeginning->adjSource());
			int rightFaceSize = rightFace->size();
			edge splitEdgeEnd = combEmb->split(splitEdgeBeginning);
			AssertThat(graph->numberOfNodes(), Equals(numberOfNodes+1));
			AssertThat(graph->numberOfEdges(), Equals(numberOfEdges+1));
			AssertThat(combEmb->numberOfFaces(), Equals(numberOfFaces));
			AssertThat(combEmb->leftFace(splitEdgeBeginning->adjSource()), Equals(leftFace));
			AssertThat(combEmb->rightFace(splitEdgeBeginning->adjSource()), Equals(rightFace));
			AssertThat(combEmb->leftFace(splitEdgeEnd->adjSource()), Equals(leftFace));
			AssertThat(combEmb->rightFace(splitEdgeEnd->adjSource()), Equals(rightFace));
			AssertThat(leftFace->size(), Equals(leftFaceSize+1));
			AssertThat(rightFace->size(), Equals(rightFaceSize+1));
		});

		it("unsplits an edge",[&](){
			int numberOfFaces = combEmb->numberOfFaces();
			edge splitEdgeBeginning = graph->chooseEdge();
			face leftFace = combEmb->leftFace(splitEdgeBeginning->adjSource());
			int leftFaceSize = leftFace->size();
			face rightFace = combEmb->rightFace(splitEdgeBeginning->adjSource());
			int rightFaceSize = rightFace->size();
			edge splitEdgeEnd = combEmb->split(splitEdgeBeginning);
			combEmb->unsplit(splitEdgeBeginning, splitEdgeEnd);
			AssertThat(graph->numberOfNodes(), Equals(numberOfNodes));
			AssertThat(graph->numberOfEdges(), Equals(numberOfEdges));
			AssertThat(combEmb->numberOfFaces(), Equals(numberOfFaces));
			AssertThat(leftFace->size(), Equals(leftFaceSize));
			AssertThat(rightFace->size(), Equals(rightFaceSize));
		});

		it("splits a node",[&](){
			node vl = graph->chooseNode();
			int degree = vl->degree();
			adjEntry adjStartLeft = vl->firstAdj();
			adjEntry adjStartRight =  vl->lastAdj();
			node vr = combEmb->splitNode(adjStartLeft, adjStartRight);
			AssertThat(vl->degree(), Equals(degree));
			AssertThat(vr->degree(), Equals(2));
			AssertThat(graph->searchEdge(vl, vr), !IsNull());
			AssertThat(vl->firstAdj()->theEdge(), Equals(vr->firstAdj()->theEdge()));
		});

		it("contracts an edge",[&](){
			node vl = graph->chooseNode();
			int degree = vl->degree();
			adjEntry adjStartLeft = vl->firstAdj();
			adjEntry adjStartRight =  vl->lastAdj();
			node vr = combEmb->splitNode(adjStartLeft, adjStartRight);
			node contractedNode = combEmb->contract(graph->searchEdge(vl, vr));
			AssertThat(contractedNode->degree(), Equals(degree));
		});

		it("reverses an edge",[&](){
			edge e = graph->chooseEdge();
			node eSrc = e->source();
			node eTgt = e->target();
			adjEntry eSrcAdj = e->adjSource();
			face rightFace = combEmb->rightFace(eSrcAdj);
			face leftFace = combEmb->leftFace(eSrcAdj);
			combEmb->reverseEdge(e);
			AssertThat(e->source(), Equals(eTgt));
			AssertThat(e->target(), Equals(eSrc));
			eSrcAdj = e->adjSource();
			AssertThat(combEmb->rightFace(eSrcAdj), Equals(leftFace));
			AssertThat(combEmb->leftFace(eSrcAdj), Equals(rightFace));
		});

		it("removes a degree-1 node",[&](){
			node deg1 = graph->newNode();
			node differentNode = graph->chooseNode();
			while(differentNode == deg1){
				differentNode = graph->chooseNode();
			}
			graph->newEdge(differentNode, deg1);
			delete combEmb;
			combEmb = new CombinatorialEmbedding(*graph);
			int numberOfFaces = combEmb->numberOfFaces();
			face faceDeg1 = combEmb->rightFace(deg1->firstAdj());
			int sizeOfFace = faceDeg1->size();
			AssertThat(combEmb->leftFace(deg1->firstAdj()), Equals(faceDeg1));
			combEmb->removeDeg1(deg1);
			AssertThat(combEmb->numberOfFaces(), Equals(numberOfFaces));
			AssertThat(faceDeg1->size(), Equals(sizeOfFace-2));
		});

		describe("face splitting",[&](){
			int sizeOfFace;
			face splitFace;
			adjEntry firstAdj;
			adjEntry secondAdj;

			before_each([&](){
				splitFace = combEmb->chooseFace();
				while(splitFace->size() <= 4){
					splitFace = combEmb->chooseFace();
				}
				firstAdj = splitFace->firstAdj();
				secondAdj = splitFace->nextFaceEdge(splitFace->nextFaceEdge(firstAdj));
				sizeOfFace = splitFace->size();
			});

			it("splits a face by two adjEntries of the same face",[&](){
				edge splitEdge = combEmb->splitFace(firstAdj, secondAdj);
				AssertThat(splitFace->size(), Equals(sizeOfFace-1));
				face newFace = combEmb->rightFace(firstAdj);
				AssertThat(newFace->size(), Equals(3));
				AssertThat(splitEdge, !IsNull());
				AssertThat(splitEdge->source(), Equals(firstAdj->theNode()));
				AssertThat(splitEdge->target(), Equals(secondAdj->theNode()));
			});

			it("splits a face by an adjEntry and a deg-0 node",[&](){
				node v = graph->newNode();
				edge splitEdge = combEmb->splitFace(v, secondAdj);
				AssertThat(splitFace->size(), Equals(sizeOfFace + 2));
				face newFace = combEmb->rightFace(secondAdj);
				AssertThat(newFace->size(), Equals(sizeOfFace + 2));
				AssertThat(splitEdge, !IsNull());
				AssertThat(splitEdge->source(), Equals(v));
				AssertThat(splitEdge->target(), Equals(secondAdj->theNode()));
			});

			it("splits a face by an adjEntry and a deg-0 node",[&](){
				node v = graph->newNode();
				edge splitEdge = combEmb->splitFace(firstAdj, v);
				AssertThat(splitFace->size(), Equals(sizeOfFace + 2));
				face newFace = combEmb->rightFace(firstAdj);
				AssertThat(newFace->size(), Equals(sizeOfFace + 2));
				AssertThat(splitEdge, !IsNull());
				AssertThat(splitEdge->target(), Equals(v));
				AssertThat(splitEdge->source(), Equals(firstAdj->theNode()));
			});

			it("splits a face by an adjEntry and a connected node",[&](){
				node v = graph->chooseNode();
				adjEntry adj = v->lastAdj()->faceCycleSucc()->faceCycleSucc();
				splitFace = combEmb->rightFace(adj);
				sizeOfFace = splitFace->size();
				edge splitEdge = combEmb->splitFace(adj, v);
				AssertThat(splitFace->size(), Equals(3));
				face newFace = combEmb->rightFace(adj);
				AssertThat(newFace->size(), Equals(sizeOfFace - 1));
				AssertThat(splitEdge, !IsNull());
				AssertThat(splitEdge->target(), Equals(v));
				AssertThat(splitEdge->source(), Equals(adj->theNode()));
			});

			it("splits a face by a connected node and an adjEntry",[&](){
				node v = graph->chooseNode();
				adjEntry adj = v->lastAdj()->faceCycleSucc()->faceCycleSucc();
				splitFace = combEmb->rightFace(adj);
				sizeOfFace = splitFace->size();
				edge splitEdge = combEmb->splitFace(v, adj);
				AssertThat(splitFace->size(), Equals(3));
				face newFace = combEmb->rightFace(adj);
				AssertThat(newFace->size(), Equals(sizeOfFace - 1));
				AssertThat(splitEdge, !IsNull());
				AssertThat(splitEdge->target(), Equals(v));
				AssertThat(splitEdge->source(), Equals(adj->theNode()));
			});

		});

		describe("face splitting donts",[&](){
			it("wont split a face by two adjEntries from different Faces",[&](){
				adjEntry v = graph->chooseEdge()->adjSource();
				adjEntry w = graph->chooseEdge()->adjSource();
				while(combEmb->rightFace(v) == combEmb->rightFace(w) ||
					combEmb->rightFace(v) == combEmb->leftFace(w) ||
					combEmb->leftFace(v) == combEmb->rightFace(w) ||
					combEmb->leftFace(v) == combEmb->leftFace(w)){
					w = graph->chooseEdge()->adjSource();
				}
				AssertThrows(PreconditionViolatedException, combEmb->splitFace(v, w));
			});
		});

		it("joins faces while removing the edge",[&](){
			int numberOfFaces = combEmb->numberOfFaces();
			edge delEdge = graph->chooseEdge();
			face faceLeft = combEmb->leftFace(delEdge->adjSource());
			face faceRight = combEmb->rightFace(delEdge->adjSource());
			int sizeLeft = faceLeft->size();
			int sizeRight = faceRight->size();
			face jointFace = combEmb->joinFaces(delEdge);
			AssertThat(jointFace->size(), Equals(sizeLeft + sizeRight - 2));
			AssertThat(combEmb->numberOfFaces(), Equals(numberOfFaces-1));
		});
	});

	describe("ConstCombinatorialEmbedding",[&](){
		it("knows its first face",[&](){
			face it = combEmb->firstFace();
			AssertThat(it->index(), Equals(0));
			AssertThat(it->pred(), IsNull());
			int counter = 0;
			for(; it != nullptr; it = it->succ()){
				counter++;
			}
			AssertThat(counter, Equals(combEmb->numberOfFaces()));
		});

		it("knows its last face",[&](){
			face it = combEmb->lastFace();
			AssertThat(it->index(), Equals(combEmb->maxFaceIndex()));
			AssertThat(it->succ(), IsNull());
			int counter = 0;
			for(; it != nullptr; it = it->pred()){
				counter++;
			}
			AssertThat(counter, Equals(combEmb->numberOfFaces()));
		});

		it("knows its number of faces",[&](){
			int counter = 0;
			int size = 0;
			for(face f : combEmb->faces){
				counter++;
				size += f->size();
			}
			AssertThat(size, Equals(graph->numberOfEdges()*2));
			AssertThat(counter, Equals(combEmb->numberOfFaces()));
		});

		it("knows a maximal face",[&](){
			int max = -1;
			face maxFace = nullptr;
			for(face it : combEmb->faces){
				if(it->size() > max){
					max = it->size();
					maxFace = it;
				}
			}
			AssertThat(maxFace->size(), Equals(combEmb->maximalFace()->size()));
		});

		it("chooses a face",[&](){
			for(int i=0; i<20; i++){
				AssertThat(combEmb->chooseFace(), !IsNull());
			}
		});

		it("knows which faces are incident to a node or edge",[&](){
			delete graph;
			delete combEmb;
			graph = new Graph();
			node u = graph->newNode();
			node v = graph->newNode();
			node w = graph->newNode();
			edge e = graph->newEdge(u, v);
			edge f = graph->newEdge(v, w);
			edge g = graph->newEdge(w, u);
			combEmb = new CombinatorialEmbedding(*graph);
			AssertThat(u->firstAdj()->theEdge(), Equals(e));
			face rightFace = combEmb->rightFace(e->adjSource());
			AssertThat(combEmb->rightFace(f->adjSource()), Equals(rightFace));
			AssertThat(combEmb->rightFace(g->adjSource()), Equals(rightFace));
			face leftFace = combEmb->leftFace(e->adjSource());
			AssertThat(combEmb->leftFace(f->adjSource()), Equals(leftFace));
			AssertThat(combEmb->leftFace(g->adjSource()), Equals(leftFace));
			AssertThat(combEmb->numberOfFaces(), Equals(2));
		});

		it("computes its faces",[&](){
			graph->newNode();
			graph->newEdge(graph->chooseNode(), graph->newNode());
			adjEntry adj;
			edge e = combEmb->splitFace(graph->newNode(), adj = graph->chooseNode()->firstAdj());
			int numberOfFaces = combEmb->numberOfFaces();
			int sizeOfFace = combEmb->rightFace(adj)->size();
			combEmb->computeFaces();
			AssertThat(combEmb->rightFace(e->adjSource()), Equals(combEmb->leftFace(e->adjSource())));
			AssertThat(combEmb->numberOfFaces(), Equals(numberOfFaces));
			AssertThat(combEmb->rightFace(adj)->size(), Equals(sizeOfFace));
		});

		it("has an arbitrary face set as the external face",[&](){
			face extFace = combEmb->chooseFace();
			combEmb->setExternalFace(extFace);
			AssertThat(combEmb->externalFace(), Equals(extFace));
		});

		it("recognizes bridges",[&](){
			edge bridgeEdge = graph->newEdge(graph->chooseNode(), graph->newNode());
			delete combEmb;
			combEmb = new CombinatorialEmbedding(*graph);
			AssertThat(combEmb->isBridge(bridgeEdge), IsTrue());

			delete graph;
			graph = new Graph();
			randomTree(*graph, numberOfNodes);
			delete combEmb;
			combEmb = new CombinatorialEmbedding(*graph);
			for(edge e : graph->edges){
				AssertThat(combEmb->isBridge(e), IsTrue());
			}
		});

		it("knows its faceArrayTableSize",[&](){
			delete graph;
			graph = new Graph();
			completeGraph(*graph, 1);
			delete combEmb;
			combEmb = new CombinatorialEmbedding(*graph);
			AssertThat(combEmb->faceArrayTableSize(), Equals(1<<4));
			delete graph;
			graph = new Graph();
			planarTriconnectedGraph(*graph, numberOfNodes*100, numberOfEdges*100);
			combEmb = new CombinatorialEmbedding(*graph);
			int numberOfFaces = combEmb->numberOfFaces();
			AssertThat(combEmb->faceArrayTableSize(), IsGreaterThan(numberOfFaces-1));
		});
	});
});
});
