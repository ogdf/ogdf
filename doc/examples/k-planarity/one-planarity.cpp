#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GraphSets.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarityBacktracking.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarization.h>
#include <ogdf/k-planarity/1-planarity_backtracking/PartialSolutionFilter.h>

#include <memory>
#include <set>
#include <vector>

using namespace ogdf;
using namespace oneplan_backtracking;
using ReturnType = Module::ReturnType;

int main() {
	// Creating a solver that uses at most 1000 threads and extracts up to 1000 Kuratowski-subdivisions in each step
	OnePlanarityBacktracking solver(1000, 1000);

	Graph K5;
	completeGraph(K5, 5);
	// Testing 1-Planarity:
	OGDF_ASSERT(solver.testOnePlanarity(K5) == ReturnType::Feasible);

	// A corresponding planarization can be obtained as follows:
	OnePlanarization planarization;
	solver.testOnePlanarity(K5, &planarization);
	OGDF_ASSERT(isPlanar(planarization));
	OGDF_ASSERT(planarization.numberOfNodes() == 6); // 1 additional crossing vertex
	OGDF_ASSERT(planarization.kiteEdges().size() == 0); // No kite edges in a complete graph
	OGDF_ASSERT(planarization.numberOfEdges()
			== 12); // K5 minus 2 edges plus 4 edges incident to crossing vertex

	// kite edges are mapped to nullptr, edges incident to crossing vertices are mapped to the corresponding crossing edge:
	for (node crossingVertex : planarization.crossingVertices()) {
		OGDF_ASSERT(crossingVertex->degree() == 4);
		std::set<edge> originalEdges;
		for (adjEntry adj : crossingVertex->adjEntries) {
			originalEdges.insert(planarization.original(adj->theEdge()));
		}
		OGDF_ASSERT(originalEdges.size() == 2);
	}

	planarEmbed(planarization); // Represents a 1-planar embedding of K5

	// Testing IC- or NIC-planarity:
	OGDF_ASSERT(solver.testICPlanarity(K5) == ReturnType::Feasible);
	OGDF_ASSERT(solver.testNICPlanarity(K5) == ReturnType::Feasible);

	// No-instances:
	Graph K7;
	completeGraph(K7, 7);
	OGDF_ASSERT(solver.testOnePlanarity(K7) == ReturnType::NoFeasibleSolution);

	// Implementing custom filters:
	class BipartiteDensityFilter : public PartialSolutionFilter {
		// Bipartite 1-planar graphs have at most 3n-8 edges. This must also hold for the planarization of a partial solution.
		bool canCut(OnePlanarization& pl) override {
			return isBipartite(pl) && pl.numberOfEdges() > 3 * pl.numberOfNodes() - 8;
		}
	};

	class MyOneplanSolver : public OnePlanarityBacktracking {
	public:
		MyOneplanSolver() : OnePlanarityBacktracking() {
			m_filters.emplace_back(new BipartiteDensityFilter);
		}
	};

	Graph K45;
	completeBipartiteGraph(K45, 4, 5);
	MyOneplanSolver mySolver;
	OGDF_ASSERT(mySolver.testOnePlanarity(K45) == ReturnType::NoFeasibleSolution);
	OGDF_ASSERT(mySolver.processedNodes() <= 1);
}
