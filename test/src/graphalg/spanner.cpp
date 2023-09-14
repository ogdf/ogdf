/** \file
 * \brief Tests for spanner algorithms.
 *
 * \author Finn Stutzenstein
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
#include <ogdf/basic/graph_generators/deterministic.h>
#include <ogdf/graphalg/SpannerBasicGreedy.h>
#include <ogdf/graphalg/SpannerBaswanaSen.h>
#include <ogdf/graphalg/SpannerBerman.h>
#include <ogdf/graphalg/SpannerBermanDisconnected.h>
#include <ogdf/graphalg/SpannerElkinNeiman.h>
#include <ogdf/graphalg/SpannerKortsarzPeleg.h>

#include <graphs.h>
#include <testing.h>


template<typename T>
void makeWeighted(GraphAttributes& GA);

template<>
void makeWeighted<int>(GraphAttributes& GA) {
	GA.addAttributes(GraphAttributes::edgeIntWeight);
	int max = GA.constGraph().numberOfEdges() / 2;
	// max is m/2 so some edges have equal weight
	for (edge e : GA.constGraph().edges) {
		GA.intWeight(e) = randomNumber(0, max);
	}
};

template<>
void makeWeighted<double>(GraphAttributes& GA) {
	GA.addAttributes(GraphAttributes::edgeDoubleWeight);
	double max = static_cast<double>(GA.constGraph().numberOfEdges()) / 2.0;
	for (edge e : GA.constGraph().edges) {
		GA.doubleWeight(e) = randomDouble(0, max);
	}
};

template<template<typename> class TModule, typename TWeight>
void test(const std::set<GraphProperty>& requirements, bool directed, bool weighted,
		const std::vector<double>& stretches, GraphSizes sizes = GraphSizes()) {
	for (double stretch : stretches) {
		describe("stretch=" + to_string(stretch), [&] {
			forEachGraphItWorks(
					requirements,
					[&](const Graph& G) {
						GraphAttributes GA(G, 0);
						GA.directed() = directed;
						if (weighted) {
							makeWeighted<TWeight>(GA);
						}

						GraphCopySimple spanner;
						EdgeArray<bool> inSpanner;

						TModule<TWeight> sm;
						AssertThat(sm.call(GA, stretch, spanner, inSpanner),
								Equals(SpannerModule<TWeight>::ReturnType::Feasible));

						AssertThat(SpannerModule<TWeight>::isMultiplicativeSpanner(GA, spanner,
										   stretch),
								IsTrue());

						// check, that spanner and inSpanner are correctly set
						for (edge e : G.edges) {
							AssertThat(inSpanner[e], Equals(spanner.copy(e) != nullptr));
						}
					},
					sizes);
		});
	}
}

template<class TModule>
void testInvalidPreconditions(const GraphAttributes& GA, double stretch,
		const std::string& mustContain) {
	it("should fail precondition check", [&] {
		TModule sm;
		std::string error;
		AssertThat(sm.preconditionsOk(GA, stretch, error), IsFalse());
		AssertThat(error, Is().Containing(mustContain));
	});
}

template<class TModule>
void testInvalidPreconditions(bool directed, double stretch, const std::string& mustContain) {
	Graph G;
	GraphAttributes GA(G, 0);
	GA.directed() = directed;
	emptyGraph(G, 0);

	testInvalidPreconditions<TModule>(GA, stretch, mustContain);
}

template<class TModule>
void testTimelimit(bool directed, double stretch) {
	it("should timeout", [&] {
		Graph G;
		GraphAttributes GA(G, 0);
		GA.directed() = directed;
		completeGraph(G, 4); // Use an actual graph to avoid simple cases

		GraphCopySimple spanner;
		EdgeArray<bool> inSpanner;

		TModule sm;
		sm.setTimelimit(0);
		AssertThat(sm.call(GA, stretch, spanner, inSpanner),
				Equals(TModule::ReturnType::TimeoutInfeasible));
	});
}

template<template<typename> class TModule>
void testSpannerBerman(const std::set<GraphProperty>& requirements,
		const std::vector<double>& stretches) {
	//SpannerBerman<int>::logger.localLogLevel(Logger::Level::Default);
	auto testIntDouble = [&](bool directed, bool weighted) {
		std::string directedStr = directed ? "directed" : "undirected";
		std::string weightStr = weighted ? "weighted" : "unweighted";
		describe(directedStr + " int " + weightStr + " general graphs", [&] {
			test<TModule, int>(requirements, directed, weighted, stretches, GraphSizes(10, 20, 10));
		});
		describe(directedStr + " double " + weightStr + " general graphs", [&] {
			test<TModule, double>(requirements, directed, weighted, stretches,
					GraphSizes(10, 20, 10));
		});
	};
	testIntDouble(true, true);
	testIntDouble(false, true);
	testIntDouble(true, false);
	testIntDouble(false, false);
	describe("Graph with selfloops", [] {
		Graph G;
		GraphAttributes GA(G);
		GA.directed() = false;
		customGraph(G, 1, {{0, 0}});
		testInvalidPreconditions<TModule<int>>(GA, 2.0, "The graph is not simple");
	});
	describe("Graph with parallel edges", [] {
		Graph G;
		GraphAttributes GA(G);
		GA.directed() = false;
		customGraph(G, 2, {{0, 1}, {0, 1}});
		testInvalidPreconditions<TModule<int>>(GA, 2.0, "The graph is not simple");
	});
	testTimelimit<TModule<int>>(false, 1.0);
	testTimelimit<TModule<double>>(false, 1.0);
	testTimelimit<TModule<int>>(true, 1.0);
	testTimelimit<TModule<double>>(true, 1.0);
}

/**
 * Used for mocking a spanner algorithm used by the SpannerIteratedWrapper.
 * It allows for setting the amount of infeasible iterations, the return type,
 * a timeout in the nth iteration.
 *
 * The returned solution has m-iterations+1 amount of edges. It can be used
 * to check, if the wrapper returns the solution with the least amount of edges.
 */
template<typename TWeight>
class IteratedSpannerModule : public SpannerModule<TWeight> {
	bool m_preconditionsCalled = false;
	typename SpannerModule<TWeight>::ReturnType m_returnType;
	int m_failedIterations = 0;
	int m_timeoutAfter = -1;
	int m_iterations = 0;

public:
	IteratedSpannerModule(typename SpannerModule<TWeight>::ReturnType returnType =
								  SpannerModule<TWeight>::ReturnType::Feasible,
			int failedIterations = 0, int timeoutAfter = -1)
		: m_returnType(returnType)
		, m_failedIterations(failedIterations)
		, m_timeoutAfter(timeoutAfter) { }

	virtual bool preconditionsOk(const GraphAttributes& GA, double stretch,
			std::string& error) override {
		m_preconditionsCalled = true;
		return true;
	}

	bool preconditionsCalled() { return m_preconditionsCalled; }

private:
	virtual typename SpannerModule<TWeight>::ReturnType execute() override {
		m_iterations++;
		if (m_iterations <= m_failedIterations) {
			return SpannerModule<TWeight>::ReturnType::NoFeasibleSolution;
		}
		if (m_timeoutAfter != -1 && m_iterations >= m_timeoutAfter) {
			return SpannerModule<TWeight>::ReturnType::TimeoutInfeasible;
		}

		// setup spanner with m-iterations+1 edges
		const Graph& G = m_GA->constGraph();
		int i = G.numberOfEdges() - m_iterations + 1;
		for (edge e : G.edges) {
			if (i-- <= 0) {
				break;
			}
			m_spanner->newEdge(e);
			(*m_inSpanner)[e] = true;
		}

		return m_returnType;
	}

	using SpannerModule<TWeight>::m_GA;
	using SpannerModule<TWeight>::m_spanner;
	using SpannerModule<TWeight>::m_inSpanner;
};

void testSpannerIteratedWrapper(const std::string& should,
		std::function<void(GraphAttributes&, GraphCopySimple&, EdgeArray<bool>&)> test) {
	it(should, [&] {
		Graph G;
		randomGraph(G, 5, 10);
		GraphAttributes GA(G, 0);
		GraphCopySimple spanner;
		EdgeArray<bool> inSpanner;
		test(GA, spanner, inSpanner);
	});
}

go_bandit([] {
	describe("IteratedWrapper", [] {
		testSpannerIteratedWrapper("propagates errors",
				[](GraphAttributes& GA, GraphCopySimple& spanner, EdgeArray<bool>& inSpanner) {
					SpannerIteratedWrapper<int> wrapper(
							new IteratedSpannerModule<int>(SpannerModule<int>::ReturnType::Error, 1),
							10);
					AssertThat(wrapper.call(GA, 3.0, spanner, inSpanner),
							Equals(SpannerModule<int>::ReturnType::Error));
					AssertThat(wrapper.getExecutedIterations(), Equals(2));
				});
		testSpannerIteratedWrapper("stops after timeout",
				[](GraphAttributes& GA, GraphCopySimple& spanner, EdgeArray<bool>& inSpanner) {
					SpannerIteratedWrapper<int> wrapper(
							new IteratedSpannerModule<int>(SpannerModule<int>::ReturnType::Feasible,
									0, 0),
							10);
					AssertThat(wrapper.call(GA, 3.0, spanner, inSpanner),
							Equals(SpannerModule<int>::ReturnType::TimeoutInfeasible));
					AssertThat(wrapper.getExecutedIterations(), Equals(1));
				});
		testSpannerIteratedWrapper("returns at least one solution",
				[](GraphAttributes& GA, GraphCopySimple& spanner, EdgeArray<bool>& inSpanner) {
					SpannerIteratedWrapper<int> wrapper(
							new IteratedSpannerModule<int>(SpannerModule<int>::ReturnType::Feasible,
									9), // 9 failed iterations, the 10th one succeeds
							10);
					AssertThat(wrapper.call(GA, 3.0, spanner, inSpanner),
							Equals(SpannerModule<int>::ReturnType::Feasible));
					AssertThat(wrapper.getExecutedIterations(), Equals(10));
				});
		testSpannerIteratedWrapper("returns no feasible solution if every iteration is unfeasible",
				[](GraphAttributes& GA, GraphCopySimple& spanner, EdgeArray<bool>& inSpanner) {
					SpannerIteratedWrapper<int> wrapper(
							new IteratedSpannerModule<int>(SpannerModule<int>::ReturnType::Feasible,
									10),
							10);
					AssertThat(wrapper.call(GA, 3.0, spanner, inSpanner),
							Equals(SpannerModule<int>::ReturnType::NoFeasibleSolution));
					AssertThat(wrapper.getExecutedIterations(), Equals(10));
				});
		testSpannerIteratedWrapper("returns the best found solution",
				[](GraphAttributes& GA, GraphCopySimple& spanner, EdgeArray<bool>& inSpanner) {
					SpannerIteratedWrapper<int> wrapper(new IteratedSpannerModule<int>(), 6);
					AssertThat(wrapper.call(GA, 3.0, spanner, inSpanner),
							Equals(SpannerModule<int>::ReturnType::Feasible));
					AssertThat(wrapper.getExecutedIterations(), Equals(6));
					AssertThat(spanner.numberOfEdges(), Equals(5)); // 10-6+1 = 5
				});
		testSpannerIteratedWrapper("returns the best found solution with timeout",
				[](GraphAttributes& GA, GraphCopySimple& spanner, EdgeArray<bool>& inSpanner) {
					SpannerIteratedWrapper<int> wrapper(
							new IteratedSpannerModule<int>(SpannerModule<int>::ReturnType::Feasible,
									0, 4),
							6);
					AssertThat(wrapper.call(GA, 3.0, spanner, inSpanner),
							Equals(SpannerModule<int>::ReturnType::Feasible));
					AssertThat(wrapper.getExecutedIterations(), Equals(4));
					AssertThat(spanner.numberOfEdges(),
							Equals(8)); // 10-3+1 = 8, timeout is in 4th iteration
				});
	});
	describe("SpannerBaswanaSen", [] {
		describe("unweighted general graphs", [] {
			test<SpannerBaswanaSen, int>({GraphProperty::simple}, false, false, {1.0, 3.0, 5.0});
		});
		describe("int weighted general graphs", [] {
			test<SpannerBaswanaSen, int>({GraphProperty::simple}, false, true, {1.0, 3.0, 5.0});
		});
		describe("double weighted general graphs", [] {
			test<SpannerBaswanaSen, double>({GraphProperty::simple}, false, true, {1.0, 3.0, 5.0});
		});
		describe("Even stretch", [] {
			testInvalidPreconditions<SpannerBaswanaSen<int>>(false, 2.0, "The stretch must be odd");
		});
		describe("Non-integer stretch", [] {
			testInvalidPreconditions<SpannerBaswanaSen<int>>(false, 1.2,
					"The stretch is required to be an integer");
		});
		describe("Directed graphs", [] {
			testInvalidPreconditions<SpannerBaswanaSen<int>>(true, 1, "The graph must be undirected");
		});
		testTimelimit<SpannerBaswanaSen<int>>(false, 1.0);
		testTimelimit<SpannerBaswanaSen<double>>(false, 1.0);
	});
	describe("SpannerKortsarzPeleg", [] {
		describe("general graphs", [] {
			test<SpannerKortsarzPeleg, int>({GraphProperty::simple}, false, false, {2},
					GraphSizes(10, 50, 20));
		});
		describe("Stretch!=2", [] {
			testInvalidPreconditions<SpannerKortsarzPeleg<int>>(false, 3, "The stretch must be 2");
		});
		describe("Directed graphs", [] {
			testInvalidPreconditions<SpannerKortsarzPeleg<int>>(true, 2,
					"The graph must be undirected");
		});
		describe("Graph with selfloops", [] {
			Graph G;
			GraphAttributes GA(G);
			GA.directed() = false;
			customGraph(G, 1, {{0, 0}});
			testInvalidPreconditions<SpannerKortsarzPeleg<int>>(GA, 2.0, "The graph is not simple");
		});
		describe("Graph with parallel edges", [] {
			Graph G;
			GraphAttributes GA(G);
			GA.directed() = false;
			customGraph(G, 2, {{0, 1}, {0, 1}});
			testInvalidPreconditions<SpannerKortsarzPeleg<int>>(GA, 2.0, "The graph is not simple");
		});
		describe("Graph with int weighted edges", [] {
			Graph G;
			GraphAttributes GA(G, GraphAttributes::edgeIntWeight);
			GA.directed() = false;
			emptyGraph(G, 1);
			testInvalidPreconditions<SpannerKortsarzPeleg<int>>(GA, 2.0,
					"The graph must be unweighted");
		});
		describe("Graph with double weighted edges", [] {
			Graph G;
			GraphAttributes GA(G, GraphAttributes::edgeDoubleWeight);
			GA.directed() = false;
			emptyGraph(G, 1);
			testInvalidPreconditions<SpannerKortsarzPeleg<int>>(GA, 2.0,
					"The graph must be unweighted");
		});
		testTimelimit<SpannerKortsarzPeleg<int>>(false, 2.0);
		testTimelimit<SpannerKortsarzPeleg<double>>(false, 2.0);
	});
	describe("SpannerBasicGreedy", [] {
		describe("unweighted general graphs", [] {
			test<SpannerBasicGreedy, int>({GraphProperty::simple}, false, false,
					{2.0, 2.5, 3.0, 3.5, 4.0});
		});
		describe("int weighted general graphs", [] {
			test<SpannerBasicGreedy, int>({GraphProperty::simple}, false, true,
					{2.0, 2.5, 3.0, 3.5, 4.0});
		});
		describe("double weighted general graphs", [] {
			test<SpannerBasicGreedy, double>({GraphProperty::simple}, false, true,
					{2.0, 2.5, 3.0, 3.5, 4.0});
		});
		describe("Directed graphs", [] {
			testInvalidPreconditions<SpannerBasicGreedy<int>>(true, 2,
					"The graph must be undirected");
		});
		testTimelimit<SpannerBasicGreedy<int>>(false, 1.0);
		testTimelimit<SpannerBasicGreedy<double>>(false, 1.0);
	});
	describe("SpannerBerman", [] {
		testSpannerBerman<SpannerBerman>({GraphProperty::simple, GraphProperty::connected},
				{2.0, 2.5, 3.0, 3.5, 4.0});
		describe("Graph with multiple components", [] {
			Graph G;
			GraphAttributes GA(G);
			GA.directed() = false;
			emptyGraph(G, 2);
			testInvalidPreconditions<SpannerBerman<int>>(GA, 2.0, "The graph is not connected");
		});
	});
	describe("SpannerBermanDisconnected",
			[] { testSpannerBerman<SpannerBermanDisconnected>({GraphProperty::simple}, {2.0}); });
	describe("SpannerElkinNeiman", [] {
		describe("unweighted general graphs", [] {
			setSeed(42);
			for (double stretch : {1.0, 3.0, 5.0}) {
				describe("stretch=" + to_string(stretch), [&] {
					forEachGraphItWorks(
							{GraphProperty::simple},
							[&](const Graph& G) {
								GraphAttributes GA(G, 0);
								GA.directed() = false;

								GraphCopySimple spanner;
								EdgeArray<bool> inSpanner;

								SpannerElkinNeiman<int> sm;
								sm.setEpsilon(0.1);
								SpannerModule<int>::ReturnType result;

								int tries = 0;
								const int maxTries =
										3; // This works in combination with the seed above
								// if this raises errors in the future, you can just increse the amount of
								// max tries.
								do {
									result = sm.call(GA, stretch, spanner, inSpanner);
									AssertThat(result,
											Is().EqualTo(SpannerModule<int>::ReturnType::Feasible)
													.Or()
													.EqualTo(SpannerModule<
															int>::ReturnType::NoFeasibleSolution));
									tries++;
								} while (result != SpannerModule<int>::ReturnType::Feasible
										&& tries < maxTries);
								AssertThat(tries < maxTries, IsTrue());

								AssertThat(SpannerModule<int>::isMultiplicativeSpanner(GA, spanner,
												   stretch),
										IsTrue());

								// check, that spanner and inSpanner are correctly set
								for (edge e : G.edges) {
									AssertThat(inSpanner[e], Equals(spanner.copy(e) != nullptr));
								}
							},
							GraphSizes());
				});
			}
		});
		setSeed(42); // used for all tests below because mac/linux/windows may skip some tests so
		// the seed has to be equalized here.
		describe("Even stretch", [] {
			testInvalidPreconditions<SpannerElkinNeiman<int>>(false, 2.0, "The stretch must be odd");
		});
		describe("Non-integer stretch", [] {
			testInvalidPreconditions<SpannerElkinNeiman<int>>(false, 1.2,
					"The stretch is required to be an integer");
		});
		describe("Directed graphs", [] {
			testInvalidPreconditions<SpannerElkinNeiman<int>>(true, 1,
					"The graph must be undirected");
		});
		describe("Graph with int weighted edges", [] {
			Graph G;
			GraphAttributes GA(G, GraphAttributes::edgeIntWeight);
			GA.directed() = false;
			emptyGraph(G, 1);
			testInvalidPreconditions<SpannerElkinNeiman<int>>(GA, 3.0,
					"The graph must be unweighted");
		});
		describe("Graph with double weighted edges", [] {
			Graph G;
			GraphAttributes GA(G, GraphAttributes::edgeDoubleWeight);
			GA.directed() = false;
			emptyGraph(G, 1);
			testInvalidPreconditions<SpannerElkinNeiman<int>>(GA, 3.0,
					"The graph must be unweighted");
		});
		testTimelimit<SpannerElkinNeiman<int>>(false, 1.0);
		testTimelimit<SpannerElkinNeiman<double>>(false, 1.0);
	});
});
