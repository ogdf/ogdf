#include <ogdf/basic/AdjEntryArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/decomposition/FourBlockTree.h>

#include <algorithm>
#include <array>
#include <random>
#include <unordered_set>
#include <vector>

#include <testing.h>

constexpr size_t N = 16;

/**
 * Generate a random 4-connected plane triangulation.
 *
 * Generates a random (3,1)-ordering, construct the corresponding irreducible
 * triangulation, add one edge.
 *
 * We take care not to add any paths of length 2 between east and west (other
 * than those via north and south), such that adding an edge {east,west} in the
 * external face creates no separating triangle.
 *
 * Invariant: nodes incident to the external face have the external face between
 * their firstAdj() and their lastAdj().
 */
template<typename RANDOM_ENGINE>
void generateFourConnectedPlaneTriangulation(Graph& g, adjEntry& externalFace, int n,
		RANDOM_ENGINE& eng) {
	const node south = g.newNode();
	const node east = g.newNode();
	const node west = g.newNode();
	g.newEdge(east, south);
	g.newEdge(south->lastAdj(), Direction::after, west, Direction::after, -1);
	externalFace =
			g.newEdge(west->lastAdj(), Direction::after, east->firstAdj(), Direction::before, -1)
					->adjTarget();
	const node left = g.newNode();
	const node right = g.newNode();
	g.newEdge(east->lastAdj(), Direction::after, right, Direction::after, -1);
	g.newEdge(south->firstAdj(), Direction::before, right->lastAdj(), Direction::after, -1);
	g.newEdge(right->lastAdj(), Direction::after, left, Direction::after, -1);
	g.newEdge(south->firstAdj(), Direction::before, left->lastAdj(), Direction::after, -1);
	g.newEdge(west->firstAdj(), Direction::before, left->lastAdj(), Direction::after, -1);
	size_t maxCoveredNodes = 0; // the number of nodes we may cover in the next step
	n -= 6; // we already count the node north here
	std::normal_distribution<double> normDist;
	using uniformDist = std::uniform_int_distribution<size_t>;
	while (n > 0) {
		const double randomValue = normDist(eng);
		if (randomValue < 0) {
			// we decrease the number of nodes on the external face by
			const size_t numCoveredNodes = double(maxCoveredNodes) * std::min(1.0, -randomValue / 2);
			// we pick a random spot on the external face
			node rightmostAdjacent = east;
			for (size_t i = uniformDist(0, maxCoveredNodes + 1 - numCoveredNodes)(eng); i > 0; --i) {
				rightmostAdjacent = rightmostAdjacent->lastAdj()->twinNode();
			}
			node nextAdjacent = rightmostAdjacent->lastAdj()->twinNode();
			// we add a singleton with numCoveredNodes + 3 neighbors
			const node v = g.newNode();
			g.newEdge(rightmostAdjacent->lastAdj(), Direction::after, v, Direction::after, -1);
			for (size_t i = 0; i < numCoveredNodes + 2; ++i) {
				g.newEdge(nextAdjacent->firstAdj(), Direction::before, v->lastAdj(),
						Direction::after, -1);
				nextAdjacent = nextAdjacent->lastAdj()->twinNode();
			}
			// we update the number of nodes on the external face
			maxCoveredNodes -= numCoveredNodes;
			// and the number of nodes we still have to create
			--n;
		} else {
			// we increase the number of nodes on the external face by
			const size_t numNewNodes = double(n) * std::min(1.0, randomValue / 3);
			// we pick a random spot on the external face where the singleton fits
			node rightAdjacent = east;
			for (size_t i = uniformDist(0, maxCoveredNodes + 1)(eng); i > 0; --i) {
				rightAdjacent = rightAdjacent->lastAdj()->twinNode();
			}
			const node middleAdjacent = rightAdjacent->lastAdj()->twinNode();
			const node leftAdjacent = middleAdjacent->lastAdj()->twinNode();
			// we add a fan with numNewNodes + 1 nodes
			node v = g.newNode();
			g.newEdge(rightAdjacent->lastAdj(), Direction::after, v, Direction::after, -1);
			g.newEdge(middleAdjacent->firstAdj(), Direction::before, v->lastAdj(), Direction::after,
					-1);
			for (size_t i = 0; i < numNewNodes; ++i) {
				const node w = g.newNode();
				g.newEdge(v->lastAdj(), Direction::after, w, Direction::after, -1);
				g.newEdge(middleAdjacent->firstAdj(), Direction::before, w->lastAdj(),
						Direction::after, -1);
				v = w;
			}
			g.newEdge(leftAdjacent->firstAdj(), Direction::before, v->lastAdj(), Direction::after,
					-1);
			// we update the number of nodes on the external face
			maxCoveredNodes += numNewNodes;
			// and the number of nodes we still have to create
			n -= numNewNodes + 1;
		}
	}
	const node north = g.newNode();
	g.newEdge(east->firstAdj(), Direction::before, north, Direction::after, -1);
	for (node v = east->lastAdj()->twinNode(); v != east; v = v->lastAdj()->twinNode()) {
		g.newEdge(v->firstAdj(), Direction::before, north->lastAdj(), Direction::after, -1);
	}
}

/**
 * Contains a decent amount of info about how the 4-block tree should look.
 */
struct FourBlockTreeStructure {
	Array<int> degDist; // as a surrogate for testing isomorphism
	std::vector<FourBlockTreeStructure> children;
};

template<typename RANDOM_ENGINE>
void generatePlaneTriangulation(Graph& g, adjEntry& externalFace, FourBlockTreeStructure& ref,
		int n, int d, RANDOM_ENGINE& eng) {
	n = std::max(n, d / 2 + 3); // so that there are enough faces for the number of children we want
	g.clear();

	std::vector<adjEntry> externalFaces;
	externalFaces.reserve(1 << d);
	std::vector<std::unordered_set<adjEntry>> originalFaces;
	originalFaces.reserve(1 << d);
	std::vector<FourBlockTreeStructure> nodes;
	nodes.reserve(1 << d);

	// generate 4-connected components
	const node dummySrc = g.newNode();
	const node dummyTgt = g.newNode();
	for (int i = 0; i < 1 << d; ++i) {
		Graph tempGraph;
		adjEntry tempExtFace;
		generateFourConnectedPlaneTriangulation(tempGraph, tempExtFace, n, eng);

		originalFaces.emplace_back();
		EdgeArray<edge> copies(tempGraph, nullptr);
		for (const edge e : tempGraph.edges) {
			const edge c = g.newEdge(dummySrc, dummyTgt);
			copies[e] = c;
			originalFaces.back().insert({c->adjSource(), c->adjTarget()});
		}
		for (const node v : tempGraph.nodes) {
			const node nodeCopy = g.newNode();
			for (const adjEntry a : v->adjEntries) {
				adjEntry adjEntryCopy;
				const edge edgeCopy = copies[a->theEdge()];
				if (a->isSource()) {
					adjEntryCopy = edgeCopy->adjSource();
					g.moveSource(edgeCopy, nodeCopy);
				} else {
					adjEntryCopy = edgeCopy->adjTarget();
					g.moveTarget(edgeCopy, nodeCopy);
				}
				if (a == tempExtFace) {
					externalFaces.push_back(adjEntryCopy);
				}
			}
		}
		const adjEntry e = externalFaces.back();
		originalFaces.back().erase(e);
		originalFaces.back().erase(e->faceCyclePred());
		originalFaces.back().erase(e->faceCycleSucc());

		nodes.emplace_back();
		degreeDistribution(tempGraph, nodes.back().degDist);
	}
	g.delNode(dummySrc);
	g.delNode(dummyTgt);

	// assemble final graph
	while (d-- > 0) {
		for (int i = 0; i < 1 << d; ++i) {
			nodes[i].children.push_back(std::move(nodes.back()));

			const auto it = originalFaces[i].begin();
			const adjEntry targetFace = *it;
			originalFaces[i].erase(it);
			originalFaces[i].erase(targetFace->faceCyclePred());
			originalFaces[i].erase(targetFace->faceCycleSucc());

			const adjEntry ext = externalFaces.back();
			std::array<adjEntry, 3> outerAdjEntries = {targetFace, targetFace->faceCycleSucc(),
											targetFace->faceCyclePred()},
									innerAdjEntries = {
											ext, ext->faceCyclePred(), ext->faceCycleSucc()};
			std::array<node, 3> innerNodes = {innerAdjEntries[0]->theNode(),
					innerAdjEntries[1]->theNode(), innerAdjEntries[2]->theNode()};
			for (const size_t j : {0, 1, 2}) {
				const node v = innerNodes[j];
				const adjEntry ao = outerAdjEntries[j];
				adjEntry ai = innerAdjEntries[j];
				while (v->degree() > 0) {
					const adjEntry tmp = ai->cyclicPred();
					if (ai->isSource()) {
						g.moveSource(ai->theEdge(), ao, Direction::after);
					} else {
						g.moveTarget(ai->theEdge(), ao, Direction::after);
					}
					ai = tmp;
				}
				g.delNode(v);
			}
			for (const adjEntry a : innerAdjEntries) {
				g.delEdge(a->theEdge());
			}

			externalFaces.pop_back();
			originalFaces.pop_back();
			nodes.pop_back();
		}
	}

	externalFace = externalFaces.front();
	ref = std::move(nodes.front());
}

void orderChildren(FourBlockTree& fbt) {
	fbt.preorder([](FourBlockTree& n) -> void {
		std::sort(n.children.begin(), n.children.end(), [](const auto& lhs, const auto& rhs) -> bool {
			return lhs->children.size() < rhs->children.size();
		});
	});
}

void validateTreeStructure(const FourBlockTree& fbt, const FourBlockTreeStructure& ref) {
	Array<int> degDist;
	degreeDistribution(*fbt.g, degDist);
	AssertThat(degDist, Equals(ref.degDist));
	const size_t numChildren = ref.children.size();
	AssertThat(fbt.children.size(), Equals(numChildren));
	for (size_t i = 0; i < numChildren; ++i) {
		validateTreeStructure(*fbt.children[i], ref.children[i]);
	}
}

void validateParentPointers(const FourBlockTree& fbt) {
	AssertThat(fbt.parent, Equals(nullptr));
	AssertThat(fbt.parentFace, Equals(nullptr));
	fbt.preorder([](const FourBlockTree& n) -> void {
		for (const auto& child : n.children) {
			AssertThat(child->parent, Equals(&n));
			AssertThat(child->parentFace->graphOf(), Equals(n.g.get()));
		}
	});
}

void validateFourConnectedComponents(const FourBlockTree& fbt, const Graph* g) {
	fbt.preorder([g](const FourBlockTree& n) -> void {
		AssertThat(isSimpleUndirected(*n.g), Equals(true));
		AssertThat(3 * n.g->numberOfNodes(), Equals(n.g->numberOfEdges() + 6));
		AssertThat(n.g->representsCombEmbedding(), Equals(true));
		AssertThat(n.externalFace->graphOf(), Equals(n.g.get()));
		for (const node v : n.g->nodes) {
			AssertThat(n.originalNodes[v]->graphOf(), Equals(g));
		}
	});
}

go_bandit([]() {
	describe("FourBlockTree", []() {
		std::mt19937_64 eng(std::random_device {}());
		it("works for 4-connected graphs", [&eng]() {
			for (int n : {10, 32, 100, 320}) {
				for (size_t i = 0; i < N; ++i) {
					Graph g;
					adjEntry externalFace = nullptr;
					generateFourConnectedPlaneTriangulation(g, externalFace, n, eng);
					const auto fbt = FourBlockTree::construct(g, externalFace);
					validateFourConnectedComponents(*fbt, &g);
					AssertThat(fbt->children, IsEmpty());
				}
			}
		});
		it("works for graphs with separating triangles", [&eng]() {
			std::uniform_int_distribution<int> nDist(6, 50);
			for (int d : {1, 3, 5}) {
				for (size_t i = 0; i < N; ++i) {
					Graph g;
					adjEntry externalFace = nullptr;
					FourBlockTreeStructure ref;
					const int n = nDist(eng);
					generatePlaneTriangulation(g, externalFace, ref, n, d, eng);
					auto fbt = FourBlockTree::construct(g, externalFace);
					validateFourConnectedComponents(*fbt, &g);
					validateParentPointers(*fbt);
					orderChildren(*fbt);
					validateTreeStructure(*fbt, ref);
				}
			}
		});
	});
});
