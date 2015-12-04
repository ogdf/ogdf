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

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <resources.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/FaceSet.h>

using namespace ogdf;
using namespace bandit;


template<typename GCType>
EdgeArray<edge>* getCopyEdges(const Graph &graph, const GCType &graphCopy){
	EdgeArray<edge> *eCopy = new EdgeArray<edge>(graph);

	for(edge e : graph.edges){
		(*eCopy)[e] = graphCopy.copy(e);
	}

	return eCopy;
}

//! tests if a GraphCopy is initialized correctly
/**
 * \param vCopy a List of all the Nodes that should be active
 * \param eCopy an EdgeArray<edge> of the \a graph there every edge is assigned an edge
 * 		from GraphCopy if it is active or \a nullptr if it isn't active
 * \param graph a pointer to the graph of which graphCopy is a GraphCopy of
 * \param graphCopy a pointer to the GraphCopy
 * \param allAdjEdges if all or none of the nodes of a connected component should be active,
 *		every adjacent edge of an active node should be active too
 * \tparam the kind of graph-copy
 */
template<typename GCType>
void testInitGraph(const Graph &graph, const GCType &graphCopy,
					bool allAdjEdges,
					List<node> *vCopy = nullptr,
					EdgeArray<edge> *eCopy = nullptr)
{
	bool vCopySet = (vCopy == nullptr);
	bool eCopySet = (eCopy == nullptr);

	if(vCopySet){
		vCopy = new List<node>();
		graph.allNodes<List<node>>(*vCopy);
	}
	if(eCopySet){
		eCopy = getCopyEdges<GCType>(graph, graphCopy);
	}

	int numberOfActiveNodes = 0;
	int numberOfAdjActiveEdges = 0;
	for(node v : *vCopy){
		numberOfActiveNodes++;
		AssertThat(graphCopy.copy(v), !IsNull());
		for(adjEntry adj = v->firstAdj(); adj != nullptr; adj = adj->succ()){
			if(graphCopy.copy(adj->theEdge()) != nullptr){
				numberOfAdjActiveEdges++;
			} else {
				if(allAdjEdges){
					AssertThat(graphCopy.copy(adj->theEdge()), !IsNull());
				}
			}
		}
	}
	int nodeCounter = 0;
	for(node v : graph.nodes){
		if(graphCopy.copy(v) != nullptr){
			nodeCounter++;
		}
	}
	AssertThat(nodeCounter, Equals(numberOfActiveNodes));
	int edgeCounter = 0;
	for(edge e : graph.edges){
		if(graphCopy.copy(e) != nullptr){
			AssertThat((*eCopy)[e], !IsNull());
			edgeCounter++;
		}
	}
	AssertThat(numberOfAdjActiveEdges/2, Equals(edgeCounter));
	if(vCopySet){ delete vCopy; }
	if(eCopySet){ delete eCopy; }

}

/**
 * Tests common functionality of ogdf::GraphCopy and ogdf::GraphCopySimple.
 */
template<typename GCType>
void describeGraphCopySimple(int numberOfNodes)
{
	Graph graph;
	GCType *graphCopy;

	before_each([&](){
		randomGraph(graph,numberOfNodes,numberOfNodes*4);
		graphCopy = new GCType(graph);
	});

	after_each([&](){
#ifdef OGDF_DEBUG
		AssertThat(graphCopy->consistencyCheck(),IsTrue());
#endif
		delete graphCopy;
	});

	describe("simple initialization",[&](){
		it("is initialized with a given Graph",[&](){
			testInitGraph<GCType>(graph, *graphCopy, true);
		});

		it("is initialized with a given graph copy",[&](){
			GCType graphCopyCopy(*graphCopy);
			testInitGraph<GCType>(graph, graphCopyCopy, true);
		});

		it("is re-initialized with some other graph", [&](){
			Graph initialGraph;
			randomGraph(initialGraph, numberOfNodes/2, numberOfNodes*2);
			delete graphCopy;
			graphCopy = new GCType(initialGraph);
			graphCopy->init(graph);
			testInitGraph(graph, *graphCopy, true);
		});

		it("supports copy-construction", [&](){
			GraphCopy copy(*graphCopy);

			AssertThat(copy.numberOfNodes(), Equals(graphCopy->numberOfNodes()));
			AssertThat(copy.numberOfEdges(), Equals(graphCopy->numberOfEdges()));

			testInitGraph(graph, *graphCopy, true);
		});

		it("supports copy-construction on a modified copy", [&](){
			// slightly modify the copy to be copied
			for(int i = 0; i < numberOfNodes/4; i++) {
				graphCopy->delNode(graphCopy->chooseNode());
				graphCopy->delEdge(graphCopy->chooseEdge());
			}

			// create a single dummy (node, edge)
			node v = graphCopy->newNode();
			edge e = graphCopy->newEdge(graphCopy->chooseNode(), v);
			AssertThat(graphCopy->isDummy(v), IsTrue());
			AssertThat(graphCopy->isDummy(e), IsTrue());

			GCType copy(*graphCopy);

			AssertThat(copy.numberOfNodes(), Equals(graphCopy->numberOfNodes()));
			AssertThat(copy.numberOfEdges(), Equals(graphCopy->numberOfEdges()));

			// verify that there is exactly a single dummy (node, edge)
			bool foundDummy = false;
			for(edge e : copy.edges) {
				bool isDummy = copy.isDummy(e);
				AssertThat(foundDummy && isDummy, IsFalse());
				foundDummy |= isDummy;
			}
			AssertThat(foundDummy, IsTrue());
			foundDummy = false;
			for(node v : copy.nodes) {
				bool isDummy = copy.isDummy(v);
				AssertThat(foundDummy && isDummy, IsFalse());
				foundDummy |= isDummy;
			}
			AssertThat(foundDummy, IsTrue());
		});


		it("supports assignment", [&](){
			GraphCopy copy = *graphCopy;

			AssertThat(copy.numberOfNodes(), Equals(graphCopy->numberOfNodes()));
			AssertThat(copy.numberOfEdges(), Equals(graphCopy->numberOfEdges()));

			testInitGraph(graph, *graphCopy, true);
		});
	});

	it("manages copy and original",[&](){
		node originalNode=graph.chooseNode();
		AssertThat(graphCopy->original(graphCopy->copy(originalNode)), Equals(originalNode));
		edge originalEdge = graph.chooseEdge();
		AssertThat(graphCopy->original(graphCopy->copy(originalEdge)), Equals(originalEdge));
	});

	it("detects dummies",[&](){
		randomGraph(graph,numberOfNodes,0);
		delete graphCopy;
		graphCopy = new GCType(graph);
		AssertThat(graphCopy->isDummy(graphCopy->newEdge(graphCopy->chooseNode(),graphCopy->chooseNode())),IsTrue());
		AssertThat(graphCopy->isDummy(graphCopy->newNode()), IsTrue());
	});

	describe("edge adding",[&](){
		it("works using the original edge",[&](){
			edge origEdge = graph.chooseEdge();
			graphCopy->delEdge(graphCopy->copy(origEdge));
			edge copyEdge = graphCopy->newEdge(origEdge);
			AssertThat(copyEdge, !IsNull());
			AssertThat(graphCopy->copy(origEdge), Equals(copyEdge));
		});

		it("works using source and target",[&](){
			node firstNode = graphCopy->chooseNode();
			node secondNode = graphCopy->chooseNode();
			int degreeFirstNode = firstNode->degree();
			int degreeSecondNode = secondNode->degree();
			edge e = graphCopy->newEdge(firstNode, secondNode);
			AssertThat(e,!IsNull());
			AssertThat(e->source(), Equals(firstNode));
			AssertThat(e->target(), Equals(secondNode));
			AssertThat(firstNode->degree(), Equals(degreeFirstNode+1));
			AssertThat(secondNode->degree(), Equals(degreeSecondNode+1));
		});
	});

	it("adds new nodes",[&](){
		AssertThat(graphCopy->newNode(),!IsNull());
		AssertThat(graphCopy->numberOfNodes(),Equals(graph.numberOfNodes()+1));
	});

	it("un-splits edges",[&](){
		edge copyEdge = graphCopy->chooseEdge();
		edge copyCopyEdge = copyEdge;
		edge splitEdge = graphCopy->split(copyEdge);
		graphCopy->unsplit(copyEdge, splitEdge);
		AssertThat(graphCopy->original(copyEdge), Equals(graphCopy->original(copyCopyEdge)));
		AssertThat(copyEdge->source(), Equals(copyCopyEdge->source()));
		AssertThat(copyEdge->target(), Equals(copyCopyEdge->target()));
	});
}

go_bandit([](){
	const int numberOfNodes = 42;

	describe("GraphCopySimple", [&](){
		describeGraphCopySimple<GraphCopySimple>(numberOfNodes);
	});

	describe("GraphCopy", [&](){
		Graph graph;
		GraphCopy *graphCopy;

		before_each([&](){
			randomGraph(graph,numberOfNodes,numberOfNodes*4);
			graphCopy = new GraphCopy(graph);
		});

		after_each([&](){
			delete graphCopy;
		});

		describe("basic functionality",[&](){
			describeGraphCopySimple<GraphCopy>(numberOfNodes);
		});

		describe("initialization", [&](){
			EdgeArray<edge> eCopy;
			List<node> origNodes;

			it("can be assigned a given GraphCopy",[&](){
				GraphCopy graphCopyCopy;
				graphCopyCopy=*graphCopy;
				testInitGraph<GraphCopy>(graph, graphCopyCopy, true);
			});

			describe("creating empty copies",[&](){
				it("works with an empty graph",[&](){
					delete graphCopy;
					graphCopy = new GraphCopy();
					graph.clear();
					graphCopy->createEmpty(graph);
					AssertThat(graphCopy->numberOfNodes(),Equals(0));
					AssertThat(graphCopy->numberOfEdges(),Equals(0));
					AssertThat(&(graphCopy->original()),Equals(&graph));
				});

				it("works with a non-empty graph",[&](){
					delete graphCopy;
					graphCopy = new GraphCopy(graph);
					graphCopy->createEmpty(graph);
					AssertThat(graphCopy->numberOfNodes(),Equals(numberOfNodes));
					AssertThat(graphCopy->numberOfEdges(),Equals(numberOfNodes*4));
					AssertThat(&(graphCopy->original()),Equals(&graph));
					AssertThat(graphCopy->chooseNode(),!IsNull());
					AssertThat(graphCopy->chooseEdge(),!IsNull());
					AssertThat(graphCopy->copy(graph.chooseNode()), IsNull());
					AssertThat(graphCopy->copy(graph.chooseEdge()), IsNull());
					AssertThat(graphCopy->original(graphCopy->chooseNode()), IsNull());
					AssertThat(graphCopy->original(graphCopy->chooseEdge()), IsNull());
				});
			});

			it("is initialized by a given connected component",[&](){
				randomGraph(graph, numberOfNodes*2, numberOfNodes*3);
				Graph::CCsInfo ccs = Graph::CCsInfo(graph);
				delete graphCopy;
				graphCopy = new GraphCopy();
				int numberOfCC = ccs.numberOfCCs()-1;
				graphCopy->createEmpty(graph);
				graphCopy->initByCC(ccs, numberOfCC, eCopy);
				origNodes.clear();
				for(int i = ccs.startNode(numberOfCC); i< ccs.stopNode(numberOfCC); i++){
					origNodes.pushBack(ccs.v(i));
				}
				testInitGraph<GraphCopy>(graph, *graphCopy, false, &origNodes, &eCopy);
			});

			it("is initialized by either all or none of the nodes of a component",[&](){
				origNodes.clear();
				delete graphCopy;
				graphCopy = new GraphCopy();
				graphCopy->createEmpty(graph);
				graphCopy->initByNodes(origNodes, eCopy);
				testInitGraph<GraphCopy>(graph, *graphCopy, true, &origNodes, &eCopy);
				delete graphCopy;
				graphCopy = new GraphCopy();
				origNodes.clear();
				graph.allNodes<List<node>>(origNodes);
				eCopy = EdgeArray<edge>(graph);
				graphCopy->createEmpty(graph);
				graphCopy->initByNodes(origNodes, eCopy);
				testInitGraph<GraphCopy>(graph, *graphCopy, true, &origNodes, &eCopy);

#ifdef OGDF_DEBUG
				origNodes = List<node>();
				origNodes.pushBack(graph.firstNode());
				origNodes.pushBack(graph.lastNode());
				eCopy = EdgeArray<edge>(graph);
				delete graphCopy;
				graphCopy = new GraphCopy(graph);
				AssertThrows(PreconditionViolatedException, graphCopy->initByNodes(origNodes, eCopy));
#endif
			});

			it("is initialized by arbitrary nodes",[&](){
				eCopy = EdgeArray<edge>(graph);
				NodeArray<bool> activeNodes(graph, false);
				node actNode1 = graph.chooseNode();
				node actNode2 = actNode1->lastAdj()->twin()->theNode();
				activeNodes[actNode1] = true;
				activeNodes[actNode2] = true;
				origNodes.clear();
				origNodes.pushBack(actNode1);
				origNodes.pushBack(actNode2);
				graphCopy->createEmpty(graph);
				graphCopy->initByActiveNodes(origNodes, activeNodes, eCopy);
				List<node> asdf;
				graphCopy->allNodes<List<node>>(asdf);
				List<edge> asdfgh;
				graphCopy->allEdges(asdfgh);
				testInitGraph<GraphCopy>(graph, *graphCopy, false, &origNodes, &eCopy);
			});
		});
/*
#ifdef OGDF_DEBUG
		it("doesn't add a copied edge twice",[&](){
			AssertThrows(PreconditionViolatedException,graphCopy->newEdge(graph.chooseEdge()));
		});
#endif
*/
		it("adds copied nodes",[&](){
			int n = graph.numberOfNodes();
			AssertThat(graphCopy->newNode(graph.newNode()),!IsNull());
			AssertThat(graphCopy->numberOfNodes(),Equals(n+1));
		});

		it("returns the chain",[&](){
			edge originalEdge = graph.chooseEdge();
			List<edge> givenChain;
			givenChain.pushFront(graphCopy->split(graphCopy->copy(originalEdge)));
			givenChain.pushFront(graphCopy->split(graphCopy->copy(originalEdge)));
			givenChain.pushFront(graphCopy->copy(originalEdge));
			List<edge> returnedChain = graphCopy->chain(originalEdge);
			AssertThat(returnedChain.size(),Equals(3));
			AssertThat(returnedChain,Equals(givenChain));
		});

		it("detects reversed edges",[&](){
			edge reversedEdge;
			do {
				reversedEdge = graphCopy->chooseEdge();
			} while(reversedEdge->isSelfLoop());
			AssertThat(graphCopy->isReversed(graphCopy->original(reversedEdge)),IsFalse());
			graphCopy->reverseEdge(reversedEdge);
			AssertThat(graphCopy->isReversed(graphCopy->original(reversedEdge)),IsTrue());
		});

		describe("original embedding",[&](){
			before_each([&](){
				planarBiconnectedGraph(graph, numberOfNodes, static_cast<int>(min(numberOfNodes*2.5, numberOfNodes*3.0-6)));
				// shuffle adjacency order
				for(node v : graph.nodes) {
					for(adjEntry adj : v->adjEntries) {
						graph.swapAdjEdges(adj, randomNumber(0, 1) ? v->firstAdj() : v->lastAdj());
					}
				}
				delete graphCopy;
				graphCopy = new GraphCopy(graph);
			});

			it("works if the GraphCopy wasn't modified",[&](){
				planarEmbed(*graphCopy);
				AssertThat(graphCopy->representsCombEmbedding(), IsTrue());
				graphCopy->setOriginalEmbedding();
				AssertThat(graphCopy->genus(), Equals(graph.genus()));
			});

			it("doesn't embed split edges",[&](){
				graphCopy->split(graphCopy->chooseEdge());
				AssertThrows(PreconditionViolatedException, graphCopy->setOriginalEmbedding());
			});

			it("doesn't embed dummies",[&](){
				graphCopy->newNode();
				AssertThrows(PreconditionViolatedException, graphCopy->setOriginalEmbedding());
			});

			it("doesn't embed added edges",[&](){
				graphCopy->newEdge(graphCopy->chooseNode(),graphCopy->chooseNode());
				AssertThrows(PreconditionViolatedException, graphCopy->setOriginalEmbedding());
			});
		});

		describe("edge path",[&](){
			node t, u, v, w;
			edge tu, uv, vw, tw;

			before_each([&](){
				graph.clear();
				t = graph.newNode();
				u = graph.newNode();
				v = graph.newNode();
				w = graph.newNode();
				graph.newEdge(v, t);
				tu = graph.newEdge(t, u);
				uv = graph.newEdge(u, v);
				graph.newEdge(u, w);
				vw = graph.newEdge(v, w);
				tw = graph.newEdge(t, w);
				planarEmbed(graph);

				delete graphCopy;
				graphCopy = new GraphCopy(graph);
			});

			describe("non-embedded variant",[&](){
				before_each([&](){
					SList<adjEntry> crossedEdges;
					crossedEdges.pushBack(graphCopy->copy(uv)->adjTarget());
					graphCopy->insertEdgePath(tw, crossedEdges);
				});

				it("inserts a path",[&](){
					AssertThat(graphCopy->chain(tw).size(), Equals(2));
					AssertThat(graphCopy->chain(uv).size(), Equals(2));
					node newNode = graphCopy->lastNode();
					AssertThat(newNode->degree(), Equals(4));
					adjEntry adj = newNode->firstAdj();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(u));
					adj = adj->succ();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(v));
					adj = adj->succ();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(t));
					adj = adj->succ();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(w));
					adj = adj->succ();
				});

				it("removes a path",[&](){
					graphCopy->removeEdgePath(tw);
					AssertThat(graphCopy->chain(tw).size(), Equals(0));
					AssertThat(graphCopy->chain(uv).size(), Equals(1));
					AssertThat(graphCopy->copy(uv)->target(), Equals(graphCopy->copy(v)));
					AssertThat(graphCopy->copy(uv)->source(), Equals(graphCopy->copy(u)));
					AssertThat(graphCopy->numberOfNodes(), Equals(graph.numberOfNodes()));
					AssertThat(graphCopy->numberOfEdges(), Equals(graph.numberOfEdges() - 1));
				});
			});

			describe("embedded variant",[&](){
				CombinatorialEmbedding combEmb;
				SList<adjEntry> crossedEdges;

				before_each([&](){
					combEmb.init(*graphCopy);
					crossedEdges.clear();
					crossedEdges.pushBack(graphCopy->copy(tu)->adjSource());
					crossedEdges.pushBack(graphCopy->copy(uv)->adjTarget());
					crossedEdges.pushBack(graphCopy->copy(vw)->adjTarget());
					graphCopy->insertEdgePathEmbedded(tw, combEmb, crossedEdges);
				});

				it("inserts a path",[&](){
					AssertThat(graphCopy->chain(tw).size(), Equals(2));
					AssertThat(graphCopy->chain(uv).size(), Equals(2));
					AssertThat(graphCopy->numberOfEdges(), Equals(8));
					AssertThat(graphCopy->numberOfNodes(), Equals(5));
					node newNode = graphCopy->lastNode();
					AssertThat(newNode->degree(), Equals(4));
					adjEntry adj = newNode->firstAdj();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(u));
					adj = adj->succ();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(w));
					adj = adj->succ();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(v));
					adj = adj->succ();
					AssertThat(graphCopy->original(adj->twin()->theNode()), Equals(t));
					adj = adj->succ();
					AssertThat(combEmb.numberOfFaces(), Equals(5));
				});

				it("removes a path",[&](){
					FaceSetPure newFaces(combEmb);
					graphCopy->removeEdgePathEmbedded(combEmb, tw, newFaces);
					AssertThat(graphCopy->chain(tw).size(), Equals(0));
					edge newOldEdge = graphCopy->copy(tw);
					AssertThat(newOldEdge, IsNull());
					AssertThat(graphCopy->chain(uv).size(), Equals(1));
					AssertThat(graphCopy->numberOfEdges(), Equals(5));
					AssertThat(combEmb.rightFace(graphCopy->copy(tu)->adjSource())->size(), Equals(3));
					AssertThat(combEmb.leftFace(graphCopy->copy(tu)->adjSource())->size(), Equals(4));
					AssertThat(graphCopy->numberOfNodes(), Equals(4));
					AssertThat(combEmb.numberOfFaces(), Equals(3));
				});
			});
		});

		it("sets a copy edge and an original edge to be corresponding",[&](){
			completeGraph(graph, 2);
			delete graphCopy;
			graphCopy = new GraphCopy(graph);
			edge copyEdge = graphCopy->chooseEdge();
			edge origEdge = graphCopy->original(copyEdge);
			graphCopy->delEdge(copyEdge);
			graphCopy->setEdge(origEdge, copyEdge);
			AssertThat(graphCopy->original(copyEdge), Equals(origEdge));
			AssertThat(graphCopy->copy(origEdge), Equals(copyEdge));
		});

		it("inserts crossings",[&](){
			edge crossingEdge = graphCopy->chooseEdge();
			edge crossedEdge = graphCopy->chooseEdge();
			while(crossedEdge == crossingEdge){
				crossedEdge = graphCopy->chooseEdge();
			}
			graphCopy->insertCrossing(crossingEdge, crossedEdge, true);
			List<edge> crossingEdgeList = graphCopy->chain(graphCopy->original(crossingEdge));
			List<edge> crossedEdgeList = graphCopy->chain(graphCopy->original(crossedEdge));
			AssertThat(crossingEdgeList.size(), Equals(2));
			AssertThat(crossedEdgeList.size(), Equals(2));
			AssertThat(crossedEdgeList.front()->source(), Equals(crossedEdgeList.front()->source()));
		});

		it("removes pseudo crossings, where two edges merely touch",[&](){
			graph.clear();
			graph.newEdge(graph.newNode(), graph.newNode());
			graph.newEdge(graph.newNode(), graph.newNode());
			delete graphCopy;
			graphCopy = new GraphCopy(graph);
			edge eCopy = graphCopy->firstEdge();
			edge fCopy = graphCopy->lastEdge();
			edge eSplit = graphCopy->split(eCopy);
			edge fSplit = graphCopy->split(fCopy);
			graphCopy->contract(graphCopy->newEdge(eSplit->source(), fSplit->source()));
			graphCopy->removePseudoCrossings();
			AssertThat(graphCopy->chain(graphCopy->original(eCopy)).size(), Equals(1));
			AssertThat(graphCopy->chain(graphCopy->original(fCopy)).size(), Equals(1));
		});

#ifdef OGDF_DEBUG
		it("won't delete a split edge",[&](){
			edge splittedEdge=graphCopy->split(graphCopy->chooseEdge());
			AssertThrows(PreconditionViolatedException,graphCopy->delEdge(splittedEdge));
		});
#endif
	});
});
