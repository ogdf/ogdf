/** \file
 * \brief Implementation of the Blossom V algorithm by Kolmogorov (2009).
 *
 * \author Joshua Sangmeister
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

#pragma once

#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/graphalg/MatchingModule.h>
#include <ogdf/graphalg/matching_blossom/AlternatingTree.h>
#include <ogdf/graphalg/matching_blossom/AuxGraph.h>
#include <ogdf/graphalg/matching_blossom/BlossomVHelper.h>
#include <ogdf/graphalg/matching_blossom/Cycle.h>
#include <ogdf/graphalg/matching_blossom/Pseudonode.h>
#include <ogdf/graphalg/matching_blossom/utils.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Uncomment to print statistics
#define OGDF_BLOSSOMV_PRINT_STATS

// Helper macros for statistics
#ifdef OGDF_BLOSSOMV_PRINT_STATS
#	define OGDF_BLOSSOMV_START_TIMER() OGDF_BLOSSOMV_START_NAMED_TIMER(__timestamp)
#	define OGDF_BLOSSOMV_START_NAMED_TIMER(timer) auto timer = now()
#	define OGDF_BLOSSOMV_END_TIMER(stat) OGDF_BLOSSOMV_END_NAMED_TIMER(__timestamp, stat)
#	define OGDF_BLOSSOMV_END_NAMED_TIMER(timer, stat) m_stats[stat].add(end(timer))
#	define OGDF_BLOSSOMV_ADD_STAT(stat) m_stats[stat].add(0)
#else
#	define OGDF_BLOSSOMV_START_TIMER()
#	define OGDF_BLOSSOMV_START_NAMED_TIMER(timer)
#	define OGDF_BLOSSOMV_END_TIMER(name)
#	define OGDF_BLOSSOMV_END_NAMED_TIMER(timer, stat)
#	define OGDF_BLOSSOMV_ADD_STAT(stat)
#endif

namespace ogdf {

/**
 * Implementation of the Blossom-V algorithm by Kolmogorov (2009) for Minimum
 * Weight Perfect Matching.
 */
template<typename TWeight>
class MatchingBlossomV : public MatchingModule<TWeight> {
	BlossomVHelper<TWeight> m_helper;
	AuxGraph<TWeight> m_auxGraph;

	//! The current aux node in the main loop.
	AuxNode<TWeight>* m_currentAuxNode = nullptr;

	//! The epsilon test for floating point comparisons.
	EpsilonTest m_eps;

#ifdef OGDF_BLOSSOMV_PRINT_STATS
	//! Structure to store statistics.
	struct stats {
		int count = 0;
		long long time = 0;

		void add(long long t) {
			count++;
			time += t;
		}

		long long ms() { return time / 1000000; }
	};

	//! A mapping of all statistic names to their values.
	std::unordered_map<std::string, stats> m_stats;

	//! Get the current time point.
	std::chrono::high_resolution_clock::time_point now() {
		return std::chrono::high_resolution_clock::now();
	}

	//! Get the time difference between \p start and now in nanoseconds.
	long long end(std::chrono::high_resolution_clock::time_point start) {
		auto end = std::chrono::high_resolution_clock::now();
		return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	}
#endif

	using Logger::lout;

	//! Helper function to log with high priority.
	std::ostream& louth() { return lout(Logger::Level::High); }

	//! Helper function to advance the current aux node to the next one.
	void advanceCurrentAuxNode() {
		auto succ = m_currentAuxNode->graphNode()->succ();
		if (!succ) {
			succ = m_auxGraph.graph().firstNode();
		}
		m_currentAuxNode = m_auxGraph.auxNode(succ);
	}

#ifdef OGDF_HEAVY_DEBUG
	//! Debug function to assert that all helper variables are in a consistent state.
	void assertConsistency() {
		// aux graph & trees
		for (auto auxNode : m_auxGraph.nodes()) {
			OGDF_ASSERT(m_auxGraph.auxNode(auxNode->graphNode()) == auxNode);
			auto tree = auxNode->tree();
			for (node v : tree.evenNodes) {
				OGDF_ASSERT(m_auxGraph.treeAuxNode(v) == auxNode);
				OGDF_ASSERT(!tree.isOdd(v));
			}
			for (node v : tree.oddNodes) {
				OGDF_ASSERT(m_auxGraph.treeAuxNode(v) == auxNode);
				OGDF_ASSERT(!tree.isEven(v));
			}
			// tree structure
			for (node v : tree.evenNodes) {
				// all even nodes except the root have a corresponding matching edge
				if (v != tree.root()) {
					OGDF_ASSERT(m_helper.matching(v) == tree.evenBackEdge(v));
				}
			}
			for (node v : tree.oddNodes) {
				// all odd nodes have both a back edge in the tree and a forward matching
				// edge, since the matching edges are always defined for both ends
				node w = m_helper.matching(v)->opposite(v);
				OGDF_ASSERT(tree.isEven(w));
			}
			for (auto adj : auxNode->graphNode()->adjEntries) {
				OGDF_ASSERT(m_auxGraph.auxEdge(adj->theEdge()) != nullptr);
			}
		}
		// assert that all matching edges are correctly defined on both ends
		int matchedNodes = 0;
		for (node v : m_helper.graph().nodes) {
			if (edge matchingEdge = m_helper.matching(v)) {
				matchedNodes++;
				OGDF_ASSERT(matchingEdge->isIncident(v));
				OGDF_ASSERT(m_helper.matching(matchingEdge->opposite(v)) == matchingEdge);
			}
		}
		// matching
		std::unordered_set<edge> hiddenEdges;
		for (auto entry : m_helper.pseudonodes()) {
			auto pseudonode = entry.second;
			for (edge e : pseudonode->cycle->edgeOrder()) {
				hiddenEdges.insert(e);
			}
		}
		OGDF_ASSERT(matchedNodes + m_auxGraph.graph().numberOfNodes() + hiddenEdges.size()
				== (size_t)m_helper.graph().numberOfNodes());
		// priority queues
		auto checkAllPQs = [&](edge e, BlossomPQ<edge, TWeight>* shouldContain = nullptr) {
			for (auto auxNode : m_auxGraph.nodes()) {
				for (auto* pq : {&auxNode->evenFreeEdges(), &auxNode->evenEvenEdges()}) {
					OGDF_ASSERT(pq->contains(e) == (pq == shouldContain));
				}
			}
			for (auto auxEdge : m_auxGraph.edges()) {
				for (auto* pq : {&auxEdge->evenEvenEdges(), &auxEdge->evenOddEdges(),
							 &auxEdge->oddEvenEdges()}) {
					OGDF_ASSERT(pq->contains(e) == (pq == shouldContain));
				}
			}
			if (shouldContain != nullptr) {
				OGDF_ASSERT(shouldContain->priority(e) == m_helper.getReducedWeight(e));
			}
		};
		for (edge e : m_helper.graph().edges) {
			bool sourceInGraph = m_helper.repr(e->source()) == e->source();
			bool targetInGraph = m_helper.repr(e->target()) == e->target();
			OGDF_ASSERT(sourceInGraph == targetInGraph);
			if (!sourceInGraph) {
				checkAllPQs(e);
				continue;
			}
			OGDF_ASSERT(m_helper.getRealReducedWeight(e) >= 0);
			auto sourceNode = m_auxGraph.treeAuxNode(e->source());
			auto targetNode = m_auxGraph.treeAuxNode(e->target());
			bool sourceIsEven = sourceNode && sourceNode->tree().isEven(e->source());
			bool targetIsEven = targetNode && targetNode->tree().isEven(e->target());
			if (hiddenEdges.find(e) != hiddenEdges.end()) {
				OGDF_ASSERT(sourceNode == nullptr && targetNode == nullptr);
			}
			if (sourceNode == nullptr && targetNode == nullptr) {
				checkAllPQs(e);
				OGDF_ASSERT(m_helper.getRealReducedWeight(e) == m_helper.getReducedWeight(e));
			} else if (sourceNode == nullptr) {
				if (targetIsEven) {
					checkAllPQs(e, &targetNode->evenFreeEdges());
				}
			} else if (targetNode == nullptr) {
				if (sourceIsEven) {
					checkAllPQs(e, &sourceNode->evenFreeEdges());
				}
			} else if (sourceNode == targetNode) {
				if (sourceIsEven && targetIsEven) {
					checkAllPQs(e, &sourceNode->evenEvenEdges());
				} else {
					checkAllPQs(e);
				}
			} else {
				if (!sourceIsEven && !targetIsEven) {
					checkAllPQs(e);
				} else {
					edge graphAuxEdge = sourceNode->graphNode()->graphOf()->searchEdge(
							sourceNode->graphNode(), targetNode->graphNode());
					OGDF_ASSERT(graphAuxEdge != nullptr);
					auto auxEdge = m_auxGraph.auxEdge(graphAuxEdge);
					if (sourceIsEven && targetIsEven) {
						checkAllPQs(e, &auxEdge->evenEvenEdges());
					} else if (sourceIsEven) {
						checkAllPQs(e, &auxEdge->evenOddEdgesFromPerspective(sourceNode));
					} else {
						checkAllPQs(e, &auxEdge->evenOddEdgesFromPerspective(targetNode));
					}
				}
			}
		}
		for (auto entry : m_helper.pseudonodes()) {
			if (m_helper.repr(entry.first) == entry.first) {
				OGDF_ASSERT(m_helper.realY(entry.first) >= 0);
				auto auxNode = m_auxGraph.treeAuxNode(entry.first);
				for (auto other : m_auxGraph.nodes()) {
					OGDF_ASSERT(other->oddPseudonodes().contains(entry.first)
							== (other == auxNode && auxNode->tree().isOdd(entry.first)));
				}
			}
		}
	}
#endif

public:
	/**
	 * @brief Construct a MatchingBlossomV instance.
	 *
	 * @param greedyInit whether or not to use the greedy initialization
	 */
	MatchingBlossomV(bool greedyInit = true) : m_helper(greedyInit), m_auxGraph(m_helper) { }

private:
	bool doCall(const Graph& G, const EdgeArray<TWeight>& weights,
			std::unordered_set<edge>& matching) {
		return _doCall(G, weights, matching);
	}

	bool doCall(const GraphAttributes& GA, std::unordered_set<edge>& matching) {
		return _doCall(GA.constGraph(), GA, matching);
	}

	//! Helper for the main call function since abstract functions cannot be templated.
	template<class WeightContainer>
	bool _doCall(const Graph& G, const WeightContainer& weights, std::unordered_set<edge>& matching) {
		// init
		lout() << "start init" << std::endl;
#ifdef OGDF_BLOSSOMV_PRINT_STATS
		m_stats.clear();
#endif
		OGDF_BLOSSOMV_START_TIMER();
		if (!m_helper.init(G, weights, &m_auxGraph)) {
			return false;
		};
		OGDF_BLOSSOMV_END_TIMER("initialize");
		lout() << "finish init" << std::endl;
#ifdef OGDF_BLOSSOMV_PRINT_STATS
		louth() << "Trees: " << m_auxGraph.graph().numberOfNodes() << std::endl;
#endif
		return findMatching(matching);
	}

	/**
	 * @brief Main function of the algorithm. Finds a minimum weight perfect matching in the graph
	 * if one exists.
	 *
	 * @returns whether a matching was found
	 */
	bool findMatching(std::unordered_set<edge>& matching) {
#ifdef OGDF_BLOSSOMV_PRINT_STATS
		std::chrono::high_resolution_clock::time_point last_path, last_four_paths;
		// Calculate the next power of two for the number of aux nodes for logging purposes
		int p = 1;
		while (p < m_auxGraph.graph().numberOfNodes()) {
			p <<= 1;
		}
#endif
		// main loop
		while (m_auxGraph.graph().numberOfNodes() > 0) {
#ifdef OGDF_HEAVY_DEBUG
			assertConsistency();
#endif
#ifdef OGDF_BLOSSOMV_PRINT_STATS
			if (m_auxGraph.graph().numberOfNodes() == p / 2) {
				p /= 2;
				lout() << "." << p << std::flush;
				if (p == 2) {
					last_path = now();
				} else if (p == 8) {
					last_four_paths = now();
				}
			}
#endif
			if (!primalChange() && !dualChange()) {
				return false;
			}
		}
		lout() << "found matching" << std::endl;
#ifdef OGDF_BLOSSOMV_PRINT_STATS
		printParallelEdgesStats();
#endif
		OGDF_BLOSSOMV_START_TIMER();
		m_helper.getOriginalMatching(matching);
		OGDF_BLOSSOMV_END_TIMER("getOriginalMatching");
#ifdef OGDF_BLOSSOMV_PRINT_STATS
		m_stats["extraLastPathTime"].add(end(last_path));
		m_stats["extraLastFourPathsTime"].add(end(last_four_paths));
		printStatistics();
#endif
		return true;
	}

	//! Executes a primal change step.
	bool primalChange() {
#ifdef OGDF_HEAVY_DEBUG
		std::unordered_map<edge, TWeight> edgeWeights;
		m_helper.getReducedEdgeWeights(edgeWeights);
#endif
		if (!m_currentAuxNode) {
			m_currentAuxNode = m_auxGraph.auxNode(m_auxGraph.graph().firstNode());
		}
		auto start = m_currentAuxNode;
		AuxNode<TWeight>* auxNode;
		// iterate all aux nodes in cyclic order
		do {
			auxNode = m_currentAuxNode;
			advanceCurrentAuxNode();
			// invalidate currentEdges of all aux nodes
			m_helper.currentIteration++;
			if (findMatchingAugmentation(auxNode) || findTreeAugmentation(auxNode)
					|| findShrinkableCycle(auxNode) || findExpandablePseudonode(auxNode)) {
#ifdef OGDF_HEAVY_DEBUG
				m_helper.assertConsistentEdgeWeights(edgeWeights);
#endif
				return true;
			}
		} while (m_currentAuxNode != start);
		m_currentAuxNode = nullptr;
		return false;
	}

	/**
	 * @brief Finds and executes a matching augmentation in \p auxNode, if possible.
	 *
	 * @return whether the matching was augmented
	 */
	bool findMatchingAugmentation(AuxNode<TWeight>* auxNode) {
		OGDF_BLOSSOMV_START_TIMER();
		AuxEdge<TWeight>* auxEdge;
		AuxNode<TWeight>* other;
		edge minEdge;
		for (auto adj : auxNode->graphNode()->adjEntries) {
			auxEdge = m_auxGraph.auxEdge(adj->theEdge());
			other = m_auxGraph.auxNode(adj->twinNode());

			// set the currentEdge for all future operations in this iteration
			other->setCurrentEdge(auxEdge);
			// check if augmentation is possible
			minEdge = m_helper.getTopEligibleElement(auxEdge->evenEvenEdges());
			if (minEdge != nullptr) {
				OGDF_BLOSSOMV_END_TIMER("findAugment");
				augment(minEdge);
				return true;
			}
		}
		OGDF_BLOSSOMV_END_TIMER("findAugment");
		return false;
	}

	/**
	 * @brief Finds and executes a tree augmentation in \p auxNode, if possible.
	 *
	 * @return whether the tree was augmented
	 */
	bool findTreeAugmentation(AuxNode<TWeight>* auxNode) {
		OGDF_BLOSSOMV_START_TIMER();
		edge minEdge = m_helper.getTopEligibleElement(auxNode->evenFreeEdges());
		if (minEdge != nullptr) {
			OGDF_BLOSSOMV_END_TIMER("findGrow");
			grow(minEdge);
			return true;
		}
		OGDF_BLOSSOMV_END_TIMER("findGrow");
		return false;
	}

	/**
	 * @brief Finds and shrinks an odd cycle in \p auxNode, if possible.
	 *
	 * @return whether a cycle was shrunken
	 */
	bool findShrinkableCycle(AuxNode<TWeight>* auxNode) {
		OGDF_BLOSSOMV_START_TIMER();
		edge minEdge = m_helper.getTopEligibleElement(auxNode->evenEvenEdges());
		if (minEdge != nullptr) {
			OGDF_BLOSSOMV_END_TIMER("findShrink");
			shrink(minEdge);
			return true;
		}
		OGDF_BLOSSOMV_END_TIMER("findShrink");
		return false;
	}

	/**
	 * @brief Finds and expands an odd pseudonode in \p auxNode, if possible.
	 *
	 * @return whether a pseudonode was expanded
	 */
	bool findExpandablePseudonode(AuxNode<TWeight>* auxNode) {
		OGDF_BLOSSOMV_START_TIMER();
		node minNode = m_helper.getTopEligibleElement(auxNode->oddPseudonodes());
		if (minNode != nullptr) {
			auxNode->oddPseudonodes().pop();
			Pseudonode* pseudonode = m_helper.pseudonode(minNode);
			OGDF_BLOSSOMV_END_TIMER("findExpand");
			expand(pseudonode);
			return true;
		}
		OGDF_BLOSSOMV_END_TIMER("findExpand");
		return false;
	}

	//! Executes a dual change step.
	bool dualChange() {
		OGDF_BLOSSOMV_START_TIMER();
		NodeArray<TWeight> deltas(m_auxGraph.graph(), 0);
		TWeight delta, otherDelta;
		AuxEdge<TWeight>* auxEdge;
		node w;

		// Calculates the connected components of the aux graph and sets the same delta for nodes in
		// the same component.
		std::vector<std::unordered_set<node>> components;
		m_auxGraph.connectedComponents(components);
		for (auto& component : components) {
			delta = infinity<TWeight>();
			AuxNode<TWeight>* auxNode;
			for (node v : component) {
				auxNode = m_auxGraph.auxNode(v);
				if (!auxNode->evenEvenEdges().empty()) {
					delta = std::min(delta,
							m_helper.getRealTopPriority(auxNode->evenEvenEdges()) / 2);
				}
				delta = std::min({delta, m_helper.getRealTopPriority(auxNode->evenFreeEdges()),
						m_helper.getRealTopPriority(auxNode->oddPseudonodes())});
				for (auto& adj : v->adjEntries) {
					w = adj->twinNode();
					auxEdge = m_auxGraph.auxEdge(adj->theEdge());
					if (component.find(w) != component.end()) {
						// same component = same delta
						if (!auxEdge->evenEvenEdges().empty()) {
							delta = std::min(delta,
									m_helper.getRealTopPriority(auxEdge->evenEvenEdges()) / 2);
						}
					} else {
						otherDelta = deltas[w];
						if (!auxEdge->evenEvenEdges().empty()) {
							delta = std::min(delta,
									m_helper.getRealTopPriority(auxEdge->evenEvenEdges())
											- otherDelta);
						}
						auto& evenOddPQ = auxEdge->evenOddEdgesFromPerspective(auxNode);
						if (!evenOddPQ.empty()) {
							delta = std::min(delta,
									m_helper.getRealTopPriority(evenOddPQ) + otherDelta);
						}
					}
					// if already at zero, no dual change can be made in this tree
					if (m_helper.isZero(delta)) {
						break;
					}
				}
			}
			// apply delta to all nodes of the component
			if (delta > 0 && delta < infinity<TWeight>()) {
				for (auto& v : component) {
					deltas[v] = delta;
				}
			}
		}

		bool dualChange = false;
		for (auto& auxNode : m_auxGraph.nodes()) {
			delta = deltas[auxNode->graphNode()];
			if (delta > 0) {
				dualChange = true;
				auxNode->addDelta(delta);
			}
			lout() << "Delta for " << auxNode->tree().root() << ": " << delta << std::endl;
		}
		OGDF_BLOSSOMV_END_TIMER("dualChange");
		return dualChange;
	}

	//! Augment the matching with \p augmentationEdge.
	void augment(edge augmentationEdge) {
		OGDF_BLOSSOMV_START_TIMER();
		OGDF_ASSERT(augmentationEdge->graphOf() == &m_helper.graph());
		for (auto augmentationEdgeNode : augmentationEdge->nodes()) {
			auto auxNode = m_auxGraph.treeAuxNode(augmentationEdgeNode);
			auto oppositeAuxNode =
					m_auxGraph.treeAuxNode(augmentationEdge->opposite(augmentationEdgeNode));
			// merge PQs
			std::vector<edge> edgesToUpdate;
			for (auto adj : auxNode->graphNode()->adjEntries) {
				auto other = m_auxGraph.auxNode(adj->twinNode());
				if (other != auxNode && other != oppositeAuxNode) {
					auto otherEdge = m_auxGraph.auxEdge(adj->theEdge());
					for (auto pq : {&otherEdge->evenEvenEdges(),
								 &otherEdge->evenOddEdgesFromPerspective(other)}) {
						for (edge e : *pq) {
							// delay insertion into PQs until the weights have been updated
							edgesToUpdate.push_back(e);
						}
					}
				}
			}
			auxNode->tree().augmentMatching(augmentationEdge);
			if (m_currentAuxNode == auxNode) {
				advanceCurrentAuxNode();
			}
			m_auxGraph.deleteNode(auxNode);
			for (edge e : edgesToUpdate) {
				auto an = m_auxGraph.auxNodeForEdge(e);
				an->addEvenFreeEdge(e);
			}
		}
		OGDF_BLOSSOMV_END_TIMER("augment");
		lout() << "Matching augmented with " << augmentationEdge << std::endl;
	}

	//! Augment the corresponding tree with \p augmentationEdge.
	void grow(edge newEdge) {
		OGDF_BLOSSOMV_START_TIMER();
		auto auxNode = m_auxGraph.auxNodeForEdge(newEdge);
		auto& tree = auxNode->tree();
		// (tree - u) - v -- w, where v -- w is the matching edge
		node u = tree.commonNode(newEdge);
		node v = newEdge->opposite(u);
		edge matchingEdge = m_helper.matching(v);
		node w = matchingEdge->opposite(v);
		tree.grow(u, newEdge, matchingEdge);

		// adjust y values
		for (node x : matchingEdge->nodes()) {
			m_auxGraph.setAuxNode(x, auxNode);
			m_helper.y(x) -= auxNode->delta(x);
		}

		// update data priority queues
		if (m_helper.isPseudonode(v)) {
			auxNode->addOddPseudonode(v);
		}
		for (auto adj : v->adjEntries) {
			node x = adj->twinNode();
			edge e = adj->theEdge();
			auto xAuxNode = m_auxGraph.treeAuxNode(x);
			if (xAuxNode != nullptr && xAuxNode->tree().isEven(x)) {
				if (x != w) {
					xAuxNode->evenFreeEdges().remove(e);
				}
				if (xAuxNode != auxNode) {
					m_auxGraph.assertCurrentEdge(xAuxNode, auxNode)
							->addEvenOddEdgeFromPerspective(e, xAuxNode);
				}
			}
		}
		edge augmentationEdge = nullptr;
		for (auto adj : w->adjEntries) {
			node x = adj->twinNode();
			edge e = adj->theEdge();
			auto xAuxNode = m_auxGraph.treeAuxNode(x);
			if (xAuxNode == nullptr) {
				// x is free
				auxNode->addEvenFreeEdge(e);
			} else if (xAuxNode == auxNode) {
				// x belongs to the same tree
				if (tree.isEven(x)) {
					auxNode->evenFreeEdges().remove(e);
					auxNode->addEvenEvenEdge(e);
				}
			} else {
				auto auxEdge = m_auxGraph.assertCurrentEdge(xAuxNode, auxNode);
				auto xTree = xAuxNode->tree();
				if (xTree.isOdd(x)) {
					auxEdge->addEvenOddEdgeFromPerspective(e, auxNode);
				} else {
					xAuxNode->evenFreeEdges().remove(e);
					auxEdge->addEvenEvenEdge(e);
					if (augmentationEdge == nullptr && m_helper.isEqualityEdge(e)) {
						augmentationEdge = e;
					}
				}
			}
		}
		OGDF_BLOSSOMV_END_TIMER("grow");
		lout() << "Tree augmented with " << newEdge << std::endl;
		if (augmentationEdge != nullptr) {
			OGDF_BLOSSOMV_ADD_STAT("extraAugmentShortcut");
			augment(augmentationEdge);
		}
	}

	//! Shrink the odd cycle induced by \p cycleEdge and the tree determined by its endpoints.
	void shrink(edge cycleEdge) {
		OGDF_BLOSSOMV_START_TIMER();
		auto auxNode = m_auxGraph.auxNodeForEdge(cycleEdge);
		auto& tree = auxNode->tree();
		Cycle* cycle = tree.getCycle(cycleEdge);
		OGDF_BLOSSOMV_END_TIMER("extraGetCycle");
		OGDF_BLOSSOMV_START_NAMED_TIMER(beforeShrinkTimer);
		// update priority queues
		std::vector<edge> oddEdges;
		std::unordered_map<node, edge> bestEdges;
		edge augmentationEdge = nullptr;
		for (node u : cycle->nodes()) {
			bool isOdd = tree.isOdd(u);
			for (auto adj : u->adjEntries) {
				edge e = adj->theEdge();
				node v = adj->twinNode();
				if (isOdd) {
					if (!cycle->contains(v)) {
						oddEdges.push_back(e);
						OGDF_BLOSSOMV_ADD_STAT("extraShrinkFoundOddEdge");
					}
				} else if (e->source() == v && tree.isEven(v) && cycle->contains(v)) {
					// only remove if v is the source of the edge to prevent duplicate removal
					auxNode->evenEvenEdges().remove(e);
					OGDF_BLOSSOMV_ADD_STAT("extraShrinkRemoveEvenEvenEdge");
				}
			}
			if (isOdd && m_helper.isPseudonode(u)) {
				auxNode->oddPseudonodes().remove(u);
				OGDF_BLOSSOMV_ADD_STAT("extraShrinkOddPseudonodeChild");
			}
		}

		// update y values and tree associations
		for (node u : cycle->nodes()) {
			m_helper.y(u) = m_helper.realY(u);
			m_auxGraph.setAuxNode(u, nullptr);
		}
		OGDF_BLOSSOMV_END_NAMED_TIMER(beforeShrinkTimer, "extraShrinkUpdatePQsBefore");

		// do the actual shrinking
		std::vector<std::tuple<edge, bool>> selfLoops;
		OGDF_BLOSSOMV_START_NAMED_TIMER(actualShrinkTimer);
		Pseudonode* pseudonode = tree.shrink(cycle, selfLoops);
		OGDF_BLOSSOMV_END_NAMED_TIMER(actualShrinkTimer, "extraActualShrink");
		OGDF_BLOSSOMV_START_NAMED_TIMER(afterShrinkTimer);
		node newNode = pseudonode->graphNode;
		m_auxGraph.setAuxNode(newNode, auxNode);
		// set y to -delta so that realY is 0 (newNode is always even)
		m_helper.y(newNode) = -auxNode->delta();

		// remove self loops from pqs
		std::unordered_set<edge> hiddenEdges;
		for (auto& entry : selfLoops) {
			edge e = std::get<0>(entry);
			bool isEven = std::get<1>(entry);
			node vInner = m_helper.getOppositeBaseNode(e, newNode);
			node v = m_helper.repr(vInner);
			hiddenEdges.insert(e);
			auto other = m_auxGraph.treeAuxNode(v);
			if (other == nullptr) {
				if (isEven) {
					auxNode->evenFreeEdges().remove(e);
				}
			} else if (other == auxNode) {
				if (isEven && tree.isEven(v)) {
					auxNode->evenEvenEdges().remove(e);
				}
			} else {
				if (isEven && other->tree().isEven(v)) {
					other->currentEdge()->evenEvenEdges().remove(e);
				}
				if (isEven && other->tree().isOdd(v)) {
					other->currentEdge()->evenOddEdgesFromPerspective(auxNode).remove(e);
				}
				if (!isEven && other->tree().isEven(v)) {
					other->currentEdge()->evenOddEdgesFromPerspective(other).remove(e);
				}
			}
		}

		// add edges back to correct priority queues after the costs have been updated
		for (edge e : oddEdges) {
			if (hiddenEdges.find(e) != hiddenEdges.end()) {
				continue;
			}
			node v = e->opposite(newNode);
			auto other = m_auxGraph.treeAuxNode(v);
			if (other == nullptr) {
				// v is free -> add edge to evenFreeEdges
				auxNode->addEvenFreeEdge(e);
			} else if (other->tree().isEven(v)) {
				if (other == auxNode) {
					auxNode->addEvenEvenEdge(e);
				} else {
					// u-v is in the evenOddEdges/oddEvenEdges of other->currentEdge and has to
					// be switched to the evenEvenEdges, since the new pseudonode is even
					auto auxEdge = m_auxGraph.assertCurrentEdge(other, auxNode);
					auxEdge->evenOddEdgesFromPerspective(other).remove(e);
					auxEdge->addEvenEvenEdge(e);
					if (augmentationEdge == nullptr && m_helper.isEqualityEdge(e)) {
						augmentationEdge = e;
					}
				}
			} else if (other != auxNode) {
				m_auxGraph.assertCurrentEdge(other, auxNode)->addEvenOddEdgeFromPerspective(e, auxNode);
			}
		}
		OGDF_BLOSSOMV_END_NAMED_TIMER(afterShrinkTimer, "extraShrinkUpdatePQsAfter");
		OGDF_BLOSSOMV_END_TIMER("shrink");
		OGDF_BLOSSOMV_ADD_STAT("extraShrinkSize" + std::to_string(cycle->nodes().size()));
		lout() << std::endl << "Shrank at " << cycleEdge << " into " << newNode << std::endl;
		if (augmentationEdge != nullptr) {
			OGDF_BLOSSOMV_ADD_STAT("extraAugmentShortcut");
			augment(augmentationEdge);
		}
	}

	//! Expand the given \p pseudonode.
	void expand(Pseudonode* pseudonode) {
		OGDF_BLOSSOMV_START_TIMER();
		auto auxNode = m_auxGraph.treeAuxNode(pseudonode->graphNode);
		auto& tree = auxNode->tree();
		auto cycle = pseudonode->cycle;
		int nodeIndex = pseudonode->graphNode->index();
		tree.expand(pseudonode);

		// update y values and tree associations
		for (node u : cycle->nodes()) {
			if (tree.contains(u)) {
				m_auxGraph.setAuxNode(u, auxNode);
				m_helper.y(u) -= auxNode->delta(u);
			}
		}

		// update priority queues
		edge augmentationEdge = nullptr;
		for (node u : cycle->nodes()) {
			if (tree.isEven(u)) {
				for (auto adj : u->adjEntries) {
					edge e = adj->theEdge();
					node v = adj->twinNode();
					auto other = m_auxGraph.treeAuxNode(v);
					if (other == nullptr) {
						// v is free
						auxNode->addEvenFreeEdge(e);
					} else if (other == auxNode) {
						// v belongs to the same tree
						// prevent double adding for cycle edges by only adding if either v is not
						// part of the cycle or v is the source of the edge
						if (tree.isEven(v) && (!cycle->contains(v) || e->source() == v)) {
							auxNode->addEvenEvenEdge(e);
						}
					} else {
						auto auxEdge = m_auxGraph.assertCurrentEdge(other, auxNode);
						if (other->tree().isEven(v)) {
							if (auxEdge->evenOddEdgesFromPerspective(other).contains(e)) {
								auxEdge->evenOddEdgesFromPerspective(other).remove(e);
							}
							auxEdge->addEvenEvenEdge(e);
							if (augmentationEdge == nullptr && m_helper.isEqualityEdge(e)) {
								augmentationEdge = e;
							}
						} else {
							auxEdge->addEvenOddEdgeFromPerspective(e, auxNode);
						}
					}
				}
			} else if (tree.isOdd(u)) {
				if (m_helper.isPseudonode(u)) {
					auxNode->addOddPseudonode(u);
				}
				for (auto adj : u->adjEntries) {
					edge e = adj->theEdge();
					node v = adj->twinNode();
					auto other = m_auxGraph.treeAuxNode(v);
					if (other && other != auxNode && other->tree().isEven(v)) {
						auto auxEdge = m_auxGraph.assertCurrentEdge(other, auxNode);
						if (!auxEdge->evenOddEdgesFromPerspective(other).contains(e)) {
							auxEdge->addEvenOddEdgeFromPerspective(e, other);
						}
					}
				}
			} else {
				// u is free
				for (auto adj : u->adjEntries) {
					edge e = adj->theEdge();
					node v = adj->twinNode();
					auto other = m_auxGraph.treeAuxNode(v);
					if (other && other->tree().isEven(v)) {
						if (other != auxNode) {
							if (other->currentEdge()->evenOddEdgesFromPerspective(other).contains(e)) {
								other->currentEdge()->evenOddEdgesFromPerspective(other).remove(e);
							}
						}
						// only add if v is not part of the cycle, otherwise it is already added in
						// the if-clause above
						if (!cycle->contains(v)) {
							other->addEvenFreeEdge(e);
						}
					}
				}
			}
		}
		delete pseudonode;
		OGDF_BLOSSOMV_END_TIMER("expand");
		lout() << std::endl << "Expanded " << nodeIndex << std::endl;
		if (augmentationEdge != nullptr) {
			OGDF_BLOSSOMV_ADD_STAT("extraAugmentShortcut");
			augment(augmentationEdge);
		}
	}

#ifdef OGDF_BLOSSOMV_PRINT_STATS
	//! Print all statistics.
	void printStatistics() {
		long long total = 0;
		for (std::string key : {"initialize", "augment", "grow", "shrink", "expand", "dualChange",
					 "findAugment", "findGrow", "findShrink", "findExpand"}) {
			total += processStatisticEntry(key);
		}

		std::vector<std::string> extraKeys;
		std::vector<std::string> otherKeys;
		for (auto entry : m_stats) {
			if (entry.first.substr(0, 5) == "extra") {
				extraKeys.push_back(entry.first);
			} else {
				otherKeys.push_back(entry.first);
			}
		}
		std::sort(extraKeys.begin(), extraKeys.end());
		std::sort(otherKeys.begin(), otherKeys.end());

		for (auto key : otherKeys) {
			total += processStatisticEntry(key);
		}
		louth() << "Total tracked time: " << total << " ms" << std::endl;
		for (auto key : extraKeys) {
			processStatisticEntry(key);
		}
	}

	//! Print additional statistics about parallel edges.
	void printParallelEdgesStats() {
		int max = 0, sum = 0, total = 0;
		std::unordered_map<int, std::unordered_map<int, int>> edgeCount;
		for (edge e : m_helper.graph().edges) {
			if (m_helper.repr(e->source()) == e->source()
					&& m_helper.repr(e->target()) == e->target()) {
				OGDF_ASSERT(!e->isSelfLoop());
				total++;
				int i = std::min(e->source()->index(), e->target()->index());
				int j = std::max(e->source()->index(), e->target()->index());
				edgeCount[i][j]++;
			}
		}
		for (auto iEntry : edgeCount) {
			for (auto jEntry : iEntry.second) {
				int count = jEntry.second;
				if (count > max) {
					max = count;
				}
				sum += count - 1;
			}
		}
		louth() << "Max number of parallel edges: " << max << std::endl;
		louth() << "Total number of parallel edges: " << sum << " (" << ((double)sum / total * 100)
				<< "%)" << std::endl;
	}

	//! Print a single statistics entry.
	long long processStatisticEntry(const std::string& key) {
		louth() << key << ": " << m_stats[key].count << " (" << m_stats[key].ms() << " ms)"
				<< std::endl;
		auto time = m_stats[key].ms();
		m_stats.erase(key);
		return time;
	}
#endif
};

}
