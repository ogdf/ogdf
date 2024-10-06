/** \file
 * \brief Static tests for PCTrees
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/pctree/NodePCRotation.h>
#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>
#include <ogdf/basic/pctree/PCRegistry.h>
#include <ogdf/basic/pctree/PCTree.h>
#include <ogdf/basic/pctree/PCTreeIterators.h>
#include <ogdf/basic/pctree/util/IntrusiveList.h>
#include <ogdf/planarity/BoothLueker.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <list>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <testing.h>

using namespace ogdf;
using namespace pc_tree;
using namespace snowhouse;
using namespace bandit;

std::string treeToString(const PCTree& tree) {
	std::stringstream ss;
	ss << tree;
	return ss.str();
}

template<class Iterable>
void checkOrder(std::vector<PCNode*>& nodes, Iterable iter, std::vector<int> ids) {
	size_t index = 0;
	for (PCNode* node : iter) {
		AssertThat(node, Equals(nodes.at(ids.at(index))));
		index++;
	}
	AssertThat(index, Equals(ids.size()));
}

bool makeConsecutive(PCTree& tree, std::initializer_list<int> listIndexes) {
	std::vector<PCNode*> restriction;
	restriction.reserve(listIndexes.size());
	for (int i : listIndexes) {
		auto it = tree.getLeaves().begin();
		std::advance(it, i);
		restriction.push_back(*it);
	}
	return tree.makeConsecutive(restriction);
}

void testFromString(PCTree& tree) {
	it("produces equivalent trees when constructed with fromString", [&]() {
		std::stringstream ss;
		ss << tree;
		PCTree copy(ss.str(), true);
		AssertThat(tree.uniqueID(uid_utils::leafToID), Equals(copy.uniqueID(uid_utils::leafToID)));
	});
}

void testCopyCtor(PCTree& tree) {
	it("copies the tree correctly", [&]() {
		PCTreeNodeArray<PCNode*> map(tree);
		PCTree copy(tree, map, true);
		AssertThat(tree.uniqueID(uid_utils::leafToID), Equals(copy.uniqueID(uid_utils::leafToID)));
		FilteringPCTreeDFS dfs(tree, tree.getRootNode());
		for (PCNode* n : dfs) {
			AssertThat(map[n]->index(), Equals(n->index()));
		}
	});
}

void testLeafOrder(PCTree& tree) {
	it("considers the current leaf order a valid permutation", [&]() {
		std::vector<PCNode*> leafOrder;
		tree.currentLeafOrder(leafOrder);
		AssertThat(tree.isValidOrder(leafOrder), IsTrue());
	});
}

template<typename Iterable>
void testIterator(PCTree& tree, Iterable iterable) {
	PCTreeNodeArray<bool> visited(tree, false);

	size_t visitedLeafs = 0;
	size_t visitedPNodes = 0;
	size_t visitedCNodes = 0;
	for (PCNode* n : iterable) {
		AssertThat(visited[n], IsFalse());
		visited[n] = true;
		if (n->isLeaf()) {
			++visitedLeafs;
		} else if (n->getNodeType() == PCNodeType::PNode) {
			++visitedPNodes;
		} else {
			OGDF_ASSERT(n->getNodeType() == PCNodeType::CNode);
			++visitedCNodes;
		}
	}
	AssertThat(visitedLeafs, Equals(tree.getLeafCount()));
	AssertThat(visitedPNodes, Equals(tree.getPNodeCount()));
	AssertThat(visitedCNodes, Equals(tree.getCNodeCount()));
}

void testIterators(PCTree& tree) {
	testIterator(tree, FilteringPCTreeDFS(tree, tree.getRootNode()));
	testIterator(tree, FilteringPCTreeBFS(tree, tree.getRootNode()));
}

void testGeneric(PCTree& tree) {
	testFromString(tree);
	testCopyCtor(tree);
	testLeafOrder(tree);
	testIterators(tree);
}

void testPlanarity(int nodes, int edges, int seed, bool forcePlanar) {
	std::stringstream ss;
	ss << "correctly determines planarity of a graph (nodes: " << nodes << ", edges: " << edges
	   << ", seed: " << seed << ", forcePlanar: " << forcePlanar << ")";
	it(ss.str(), [&]() {
		Graph G;
		setSeed(seed);
		if (forcePlanar) {
			randomPlanarBiconnectedGraph(G, nodes, edges);
		} else {
			randomBiconnectedGraph(G, nodes, edges);
		}

		bool success = true;
		try {
			NodePCRotation N(G, G.lastNode());
			if (N.getLeafCount() > 2) {
				testGeneric(N);
			}
		} catch (GraphNotPlanarException& e) {
			success = false;
		}

		node last =
				G.chooseNode([&G](node n) -> bool { return n != G.lastNode() && n->degree() > 2; });
		if (last != nullptr) {
			bool success2 = true;
			try {
				NodePCRotation N(G, last);
				testGeneric(N);
			} catch (GraphNotPlanarException& e) {
				success2 = false;
			}
			AssertThat(success2, Equals(success));
		}

		if (forcePlanar) {
			AssertThat(success, IsTrue());
		} else {
			BoothLueker BL;
			AssertThat(success, Equals(BL.isPlanar(G)));
		}
	});
}

void testPlanarity() {
	std::list<std::function<int(int)>> edgeFuncs {[](int nodes) { return nodes; },
			[](int nodes) { return 2 * nodes; }, [](int nodes) { return 3 * nodes - 6; }};
	std::vector<int> seeds;
	for (int i = 0; i < 5; i++) {
		seeds.push_back(i);
		seeds.push_back((int)randomSeed());
	}

	for (int nodes = 10; nodes <= 100; nodes += 10) {
		for (auto& func : edgeFuncs) {
			for (int seed : seeds) {
				testPlanarity(nodes, func(nodes), seed, true);
				testPlanarity(nodes, func(nodes), seed, false);
			}
		}
	}
	testPlanarity(20, 40, 1705935965, true);
}

bool applyRestrictions(PCTree& t, std::initializer_list<std::initializer_list<int>> restrictions) {
	std::vector<PCNode*> leaves(t.getLeaves().begin(), t.getLeaves().end());
	for (auto restriction : restrictions) {
		std::vector<PCNode*> restrictionLeaves;
		for (int i : restriction) {
			restrictionLeaves.push_back(leaves.at(i));
		}
		if (!t.makeConsecutive(restrictionLeaves)) {
			return false;
		}
	}
	return true;
}

void testIntersection(int numLeaves, std::initializer_list<std::initializer_list<int>> restrictions1,
		std::initializer_list<std::initializer_list<int>> restrictions2) {
	PCTree t1(numLeaves);
	PCTree t2(numLeaves);
	PCTreeNodeArray<PCNode*> mapLeaves(t2);
	auto it1 = t1.getLeaves().begin();
	auto it2 = t2.getLeaves().begin();
	for (int i = 0; i < numLeaves; i++) {
		AssertThat(it1, !Equals(t1.getLeaves().end()));
		AssertThat(it2, !Equals(t2.getLeaves().end()));
		mapLeaves[*it2] = *it1;
		++it1;
		++it2;
	}
	AssertThat(it1, Equals(t1.getLeaves().end()));
	AssertThat(it2, Equals(t2.getLeaves().end()));

	AssertThat(applyRestrictions(t1, restrictions1), IsTrue());
	AssertThat(applyRestrictions(t2, restrictions2), IsTrue());

	PCTree check(numLeaves);
	AssertThat(applyRestrictions(check, restrictions1), IsTrue());
	bool possibleCheck = applyRestrictions(check, restrictions2);

	bool possibleIntersection = t1.intersect(t2, mapLeaves);
	AssertThat(possibleIntersection, Equals(possibleCheck));
	if (possibleCheck) {
		AssertThat(t1.uniqueID(uid_utils::leafToID), Equals(check.uniqueID(uid_utils::leafToID)));
	}
}

go_bandit([]() {
#ifdef OGDF_DEBUG
	PCTREE_DEBUG_CHECK_FREQ = 1;
#endif
	describe("PCTree", []() {
		it("allows creating a trivial instance", []() {
			std::vector<PCNode*> leaves;
			PCTree tree(5, &leaves);
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.isTrivial(), IsTrue());
			AssertThat(tree.getLeafCount(), Equals((size_t)5));
			AssertThat(tree.getPNodeCount(), Equals((size_t)1));
			AssertThat(tree.getCNodeCount(), Equals((size_t)0));
			AssertThat(std::equal(tree.getLeaves().begin(), tree.getLeaves().end(), leaves.begin(),
							   leaves.end()),
					IsTrue());
			AssertThat(tree.possibleOrders<size_t>(), Equals((size_t)24));
			AssertThat(treeToString(tree), Equals("0:(5, 4, 3, 2, 1)"));
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition), Equals("5:(4, 3, 2, 1, 0)"));
			PCNode* root = tree.getRootNode();
			AssertThat(root->getNodeType(), Equals(PCNodeType::PNode));
			AssertThat(root->getChildCount(), Equals((size_t)5));

			tree.makeConsecutive({leaves.at(1), leaves.at(2)});
			AssertThat(treeToString(tree), Equals("0:(6:(3, 2), 5, 4, 1)"));
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition), Equals("6:(5:[2, 1], 4, 3, 0)"));
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.isTrivial(), IsFalse());
		});

		it("correctly handles the first case where JPPZanetti fails", []() {
			std::vector<PCNode*> leaves;
			PCTree tree(7, &leaves);
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition), Equals("7:(6, 5, 4, 3, 2, 1, 0)"));
			tree.makeConsecutive({leaves.at(4), leaves.at(5)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("8:(7:[5, 4], 6, 3, 2, 1, 0)"));
			tree.makeConsecutive({leaves.at(3), leaves.at(4), leaves.at(5)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("9:(8:[7:[5, 4], 3], 6, 2, 1, 0)"));
			tree.makeConsecutive({leaves.at(0), leaves.at(1)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("10:(9:[8:[5, 4], 3], 7:[1, 0], 6, 2)"));
			tree.makeConsecutive({leaves.at(1), leaves.at(2)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("10:[9:[8:[5, 4], 3], 7:[2, 1, 0], 6]"));
			tree.makeConsecutive({leaves.at(2), leaves.at(3), leaves.at(4), leaves.at(5)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("9:[6, 8:[7:[5, 4], 3], 2, 1, 0]"));
			tree.makeConsecutive({leaves.at(3), leaves.at(4), leaves.at(5), leaves.at(6)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("9:[6, 8:[7:[5, 4], 3], 2, 1, 0]"));
			tree.makeConsecutive({leaves.at(3), leaves.at(4)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("8:[6, 7:[5, 4, 3], 2, 1, 0]"));
			tree.makeConsecutive({leaves.at(2), leaves.at(3), leaves.at(4), leaves.at(5)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(tree.uniqueID(uid_utils::nodeToPosition),
					Equals("8:[6, 7:[5, 4, 3], 2, 1, 0]"));
		});

		it("is iterated correctly", []() {
			std::vector<PCNode*> leaves;
			PCTree tree(7, &leaves);
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(treeToString(tree), Equals("0:(7, 6, 5, 4, 3, 2, 1)"));
			tree.makeConsecutive({leaves.at(0), leaves.at(1), leaves.at(2), leaves.at(3)});
			AssertThat(tree.checkValid(), IsTrue());
			AssertThat(treeToString(tree), Equals("0:(8:(4, 3, 2, 1), 7, 6, 5)"));

			FilteringPCTreeDFS walk = tree.allNodes();
			std::vector<PCNode*> nodes {walk.begin(), walk.end()};
			std::sort(nodes.begin(), nodes.end(),
					[](PCNode* a, PCNode* b) { return a->index() < b->index(); });
			AssertThat(nodes.size(), Equals((size_t)9));
			AssertThat(nodes.front()->index(), Equals((size_t)0));
			AssertThat(nodes.back()->index(), Equals((size_t)8));

			checkOrder(nodes, tree.getLeaves(), {1, 2, 3, 4, 5, 6, 7});
			checkOrder(nodes, tree.allNodes(), {0, 8, 4, 3, 2, 1, 7, 6, 5});
			checkOrder(nodes, tree.innerNodes(), {0, 8});

			FilteringPCTreeBFS bfs(tree, tree.getRootNode());
			checkOrder(nodes, bfs, {0, 5, 6, 7, 8, 1, 2, 3, 4});

			PCNode* node = tree.getRootNode()->getChild2();
			AssertThat(node->index(), Equals((size_t)8));
			checkOrder(nodes, node->children(), {1, 2, 3, 4});
			checkOrder(nodes, node->neighbors(), {1, 2, 3, 4, 0});
		});

		it("applies small restrictions correctly", []() {
			std::vector<PCNode*> added;
			PCTree T(10, &added);
			AssertThat(T.isTrivial(), IsTrue());
			AssertThat(makeConsecutive(T, {0, 1}), IsTrue());
			AssertThat(makeConsecutive(T, {2, 3}), IsTrue());
			AssertThat(makeConsecutive(T, {1, 2}), IsTrue());
			AssertThat(makeConsecutive(T, {3, 4, 5}), IsTrue());
			AssertThat(makeConsecutive(T, {1, 3}), IsFalse());
			AssertThat(T.isTrivial(), IsFalse());
			std::vector<PCNode*> expected(T.getLeaves().begin(), T.getLeaves().end());
			AssertThat(added, Equals(expected));

			testGeneric(T);
		});

		it("applies bigger restrictions correctly", []() {
			PCTree T(50);
			AssertThat(
					makeConsecutive(T,
							{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}),
					IsTrue());
			AssertThat(makeConsecutive(T, {25, 26, 27, 28, 29, 30, 31, 32, 33, 34}), IsTrue());
			AssertThat(makeConsecutive(T, {36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49}),
					IsTrue());
			AssertThat(makeConsecutive(T, {17, 42}), IsTrue());
			AssertThat(makeConsecutive(T, {5, 6, 7, 8, 9}), IsTrue());
			AssertThat(makeConsecutive(T, {11, 12, 13, 14}), IsTrue());
			AssertThat(makeConsecutive(T, {5, 6, 7, 8, 9, 10, 11, 12, 13, 14}), IsTrue());
			AssertThat(makeConsecutive(T, {8, 9}), IsTrue());
			AssertThat(makeConsecutive(T, {9, 10}), IsTrue());
			AssertThat(makeConsecutive(T, {10, 11}), IsTrue());
			AssertThat(makeConsecutive(T, {47, 48, 49, 36, 37, 38}), IsTrue());
			AssertThat(makeConsecutive(T, {47, 48}), IsTrue());
			AssertThat(makeConsecutive(T, {37, 38}), IsTrue());
			AssertThat(makeConsecutive(T, {36, 37, 38}), IsTrue());
			AssertThat(makeConsecutive(T, {47, 48, 49}), IsTrue());
			AssertThat(makeConsecutive(T, {49, 36}), IsTrue());
			AssertThat(makeConsecutive(T, {47, 32}), IsTrue());
			AssertThat(makeConsecutive(T, {34, 33, 25}), IsTrue());
			AssertThat(makeConsecutive(T, {34, 33, 32, 25}), IsTrue());
			AssertThat(makeConsecutive(T, {27, 33}), IsTrue());
			AssertThat(makeConsecutive(T, {12, 11, 10, 9, 8, 7, 6, 5, 4}), IsTrue());

			testGeneric(T);
		});

		it("detects validity of orders", []() {
			PCTree T(9);

			AssertThat(makeConsecutive(T, {0, 1, 2}), IsTrue());
			AssertThat(makeConsecutive(T, {3, 4, 5}), IsTrue());
			AssertThat(makeConsecutive(T, {0, 1}), IsTrue());
			AssertThat(makeConsecutive(T, {1, 2}), IsTrue());
			AssertThat(makeConsecutive(T, {0, 1, 2, 3}), IsTrue());

			std::vector<PCNode*> leaves(T.getLeaves().begin(), T.getLeaves().end());
			AssertThat(T.isValidOrder(leaves), IsTrue());
			std::reverse(leaves.begin(), leaves.end());
			AssertThat(T.isValidOrder(leaves), IsTrue());
			leaves.assign(T.getLeaves().begin(), T.getLeaves().end());
			std::swap(leaves.at(2), leaves.at(3));
			AssertThat(T.isValidOrder(leaves), IsFalse());

			testGeneric(T);
		});

		it("correctly applies restrictions on manually constructed trees", []() {
			PCTree T;
			auto root = T.newNode(PCNodeType::CNode);
			auto n1 = T.newNode(PCNodeType::PNode, root);
			T.insertLeaves(5, root);
			auto n2 = T.newNode(PCNodeType::PNode, root);
			T.insertLeaves(5, n2);
			T.insertLeaves(5, root);
			T.insertLeaves(5, n1);

			AssertThat(T.checkValid(), IsTrue());
			testGeneric(T);

			AssertThat(makeConsecutive(T, {7, 10, 15}), IsFalse());
			AssertThat(makeConsecutive(T, {6, 10, 11, 12, 13, 14, 17}), IsTrue());
			testGeneric(T);
		});

		it("correctly performs a simple intersection", []() {
			PCTree T(10);
			PCTree T2("0:(14:[15:(6, 5), 4, 3, 2, 1], 10, 9, 8, 7)");
			PCTreeNodeArray<PCNode*> leafMap(T2);
			auto it1 = T.getLeaves().begin();
			auto it2 = T2.getLeaves().begin();
			for (size_t i = 0; i < T.getLeafCount(); ++i) {
				AssertThat(it1, !Equals(T.getLeaves().end()));
				AssertThat(it2, !Equals(T2.getLeaves().end()));
				leafMap[*it2] = *it1;
				++it1;
				++it2;
			}
			AssertThat(it1, Equals(T.getLeaves().end()));
			AssertThat(it2, Equals(T2.getLeaves().end()));
			AssertThat(T.intersect(T2, leafMap), IsTrue());
			testGeneric(T);
		});

		describe("intersection", []() {
			it("correctly handles the trivial case", []() {
				testIntersection(10, {{0, 1, 2}}, {});
			});
			it("correctly handles another trivial case", []() {
				testIntersection(10, {}, {{0, 1, 2}});
			});
			it("correctly handles a tree with only P-Nodes", []() {
				testIntersection(10, {{3, 4, 5}}, {{0, 1, 2}, {6, 7, 8}});
			});
			it("correctly handles a simple intersection with a C-Node", []() {
				testIntersection(10, {{2, 3, 4}}, {{0, 1, 2}, {5, 6, 7}, {7, 8, 9}});
			});
			it("correctly handles a single C-Node", []() {
				testIntersection(5, {{1, 2, 3}}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});
			});
			it("correctly handles a more complicated intersection", []() {
				testIntersection(20, {{11, 12, 13, 14}, {0, 8}, {14, 9}},
						{{0, 1}, {1, 2}, {2, 3}, {6, 7, 8, 9, 10}, {11, 12}, {12, 13}, {13, 14},
								{15, 16}, {16, 17}, {17, 18}});
			});
			it("correctly handles an impossible intersection", []() {
				testIntersection(10, {{0, 1}, {1, 2}, {2, 3}}, {{0, 2}});
			});
		});
	});

	describe("NodePCRotation", []() {
		testPlanarity();

		it("computes bundle edges correctly", []() {
			Graph G12;
			auto a1 = G12.newNode(1);
			auto a2 = G12.newNode(2);
			auto a3 = G12.newNode(3);
			auto a4 = G12.newNode(4);
			auto a5 = G12.newNode(5);
			auto a6 = G12.newNode(6);
			auto a7 = G12.newNode(7);
			auto a8 = G12.newNode(8);
			auto a9 = G12.newNode(9);
			G12.newEdge(a1, a9);
			G12.newEdge(a1, a2);
			G12.newEdge(a2, a9);
			G12.newEdge(a2, a3);
			G12.newEdge(a3, a4);
			G12.newEdge(a4, a5);
			G12.newEdge(a5, a6);
			G12.newEdge(a5, a7);
			G12.newEdge(a6, a8);
			G12.newEdge(a7, a8);
			G12.newEdge(a8, a9);

			NodePCRotation test12(G12, a9, true);
			node partner = test12.getTrivialPartnerPole();
			AssertThat(partner != nullptr, IsTrue());
			EdgeSet adjEntries(G12);
			for (auto leaf : test12.getLeaves()) {
				AssertThat(test12.getPartnerEdgesForLeaf(leaf).empty(), IsFalse());
				for (edge twinEdge : test12.getPartnerEdgesForLeaf(leaf)) {
					adjEntries.insert(twinEdge);
				}
			}

			AssertThat(adjEntries.size(), Equals(partner->degree()));
			for (edge e : adjEntries.elements()) {
				AssertThat(e->isIncident(partner), IsTrue());
			}
		});
	});
});
