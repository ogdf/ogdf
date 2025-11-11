/** \file
 * \brief Implementation of OnePlanarityBacktracking methods.
 *
 * \author Matthias Pfretzschner
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
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/Module.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/Stopwatch.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/k-planarity/1-planarity_backtracking/EdgePairPartition.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarityBacktracking.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarization.h>
#include <ogdf/k-planarity/1-planarity_backtracking/PartialSolutionFilter.h>
#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/planarity/ExtractKuratowskis.h>

#include <algorithm>
#include <climits>
#include <deque>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <vector>

using namespace std;
using namespace ogdf;
using namespace ogdf::oneplan_backtracking;

bool OnePlanarityBacktracking::DFSThread::nextStep() {
	bool toBranch = m_branchStack.top();
	m_branchStack.pop();
	if (toBranch) {
		return branch();
	} else {
		cleanUp(*this);
	}

	return false;
}

void OnePlanarityBacktracking::DFSThread::cleanUp(DFSThread& thread) {
	EdgePairPartition& epp = thread.m_epp;

	if (epp.hasTodoCrossing()) {
		EdgePair ep = epp.getNextTodoCrossing();
		epp.startTransaction();
		epp.crossEdgePair(ep);
		thread.m_branchStack.push(false);
		m_solver->push(&thread, &epp);
	} else {
		epp.undoTransaction();
	}
}

//! Computes all crossable edge pairs of a given Kuratowski-subdivision.
void extractCrossableKuratowskiEdgePairs(OnePlanarization& p, EdgePairPartition& epp,
		KuratowskiWrapper& kuratowski, vector<EdgePair>& out) {
	Graph::HiddenEdgeSet hiddenEdges(p);
	vector<edge> edgeListCopy;
	for (edge e : p.edges) {
		edgeListCopy.push_back(e);
	}
	for (edge e : edgeListCopy) {
		hiddenEdges.hide(e);
	}
	for (edge e : kuratowski.edgeList) {
		hiddenEdges.restore(e);
	}

	GraphCopy minor(*reinterpret_cast<Graph*>(&p));

	EdgeArray<edge> rep(p, nullptr);
	rep[p.firstEdge()] = minor.copy(p.firstEdge());
	stack<adjEntry> stk {{p.firstEdge()->adjSource()}};

	while (!stk.empty()) {
		adjEntry adj = stk.top();
		stk.pop();
		node n = adj->twinNode();

		if (n->degree() == 2) {
			adjEntry nextAdj = n->firstAdj()->twin() == adj ? n->lastAdj() : n->firstAdj();
			rep[nextAdj->theEdge()] = rep[adj->theEdge()];
			stk.push(nextAdj);
		} else {
			for (adjEntry nextAdj : n->adjEntries) {
				if (nextAdj->twin() == adj || rep[nextAdj->theEdge()] != nullptr) {
					continue;
				}
				rep[nextAdj->theEdge()] = minor.copy(nextAdj->theEdge());
				stk.push(nextAdj);
			}
		}
	}
	for (edge e : p.edges) {
		if (rep[e] != minor.copy(e)) {
			minor.contract(minor.copy(e));
		}
	}
	OGDF_ASSERT(minor.numberOfEdges() == 10 || minor.numberOfEdges() == 9);
	for (edge e : kuratowski.edgeList) {
		OGDF_ASSERT(rep[e]);
	}

	hiddenEdges.restore();

	out.clear();
	for (edge e1 : kuratowski.edgeList) {
		for (edge e2 : kuratowski.edgeList) {
			if (e1 == e2) {
				break;
			}
			if (rep[e1] == rep[e2] || rep[e1]->isAdjacent(rep[e2])) {
				continue;
			}
			edge m1 = p.original(e1);
			edge m2 = p.original(e2);
			if (m1 == m2 || !m1 || !m2) {
				continue;
			}
			EdgePair pa(m1, m2);
			if (epp.isFree(pa)) {
				out.push_back(pa);
			}
		}
	}
}

/**
 * Either writes a set of edge pairs that need to be crossed to \p forcedCrossings or outputs an
 * exhaustive set of edge pairs that can be branched over to \p childPairs.
 */
bool computeBranches(EdgePairPartition* epp, int nKuratowskis, vector<EdgePair>& childPairs,
		set<EdgePair>& forcedCrossings) {
	OnePlanarization p(epp);
	SList<KuratowskiWrapper> output;
	BoyerMyrvold bm;
	bm.planarEmbed(p, output, nKuratowskis);

	map<EdgePair, int> pairs;
	vector<EdgePair> minKuratowskiPairs;
	int minKuratowskiPairSize = INT_MAX;

	forcedCrossings.clear();
	childPairs.clear();
	for (auto& kuratowski : output) {
		vector<EdgePair> out;
		extractCrossableKuratowskiEdgePairs(p, *epp, kuratowski, out);
		if (out.empty()) {
			return false;
		}
		if (out.size() == 1) { // only one crossable edge pair -> must be crossed
			forcedCrossings.insert(out.front());
		} else {
			for (auto& ep : out) {
				pairs[ep]++;
			}
			if ((int)out.size() < minKuratowskiPairSize) {
				minKuratowskiPairSize = out.size();
				std::swap(out, minKuratowskiPairs);
			}
		}
	}
	if (forcedCrossings.empty()) {
		OGDF_ASSERT(minKuratowskiPairs.size() > 1);
		std::swap(minKuratowskiPairs, childPairs);
		std::sort(childPairs.begin(), // Order pairs by the number of their occurrences in subdivisions.
				childPairs.end(), [&](auto& a, auto& b) { return pairs[a] > pairs[b]; });
	}
	return true;
}

bool OnePlanarityBacktracking::DFSThread::branch() {
	auto res = m_solver->verifyNode(&m_epp);
	if (res == NodeStatus::SOL) {
		return true;
	} else if (res == NodeStatus::CUT) {
		return false;
	}

	vector<EdgePair> childPairs;
	set<EdgePair> forcedCrossings;
	if (!computeBranches(&m_epp, m_solver->m_maxExtractedKuratowskis, childPairs, forcedCrossings)) {
		return false;
	}

	if (!forcedCrossings.empty()) { // apply forced crossings and push
		m_epp.startTransaction();
		for (auto forcedCrossing : forcedCrossings) {
			OGDF_ASSERT(!OnePlanarization(&m_epp).isPlanar());
			if (!m_epp.isFree(forcedCrossing)) {
				m_epp.undoTransaction();
				return false;
			}
			m_epp.crossEdgePair(forcedCrossing);
		}
		m_solver->push(this, &m_epp);
	} else {
		m_epp.startTransaction(childPairs); // store edge pairs corresponding to children
		m_branchStack.push(false); // mark on stack
		EdgePair ep = m_epp.getNextTodoCrossing();
		m_epp.startTransaction();
		m_epp.crossEdgePair(ep); // cross first edge pair and push
		m_solver->push(this, &m_epp);
	}
	return false;
}

Module::ReturnType OnePlanarityBacktracking::test(oneplan_backtracking::OneplanMode mode,
		const Graph& G, OnePlanarization* out) {
	OGDF_ASSERT(isSimple(G));
	m_processedNodes = 0;
	m_threads.clear();
	if (G.numberOfNodes() > 6 && G.numberOfEdges() > 4 * G.numberOfNodes() - 8) {
		return Module::ReturnType::NoFeasibleSolution;
	}

	StopwatchWallClock timer;
	timer.start();
	m_threads.emplace_back(G, mode, this);

	while (!m_threads.empty()) {
		for (auto it = m_threads.begin(); it != m_threads.end();) {
			if (m_timeLimit >= 0 && timer.milliSeconds() / 1000.0 > m_timeLimit) {
				return Module::ReturnType::TimeoutInfeasible;
			}

			DFSThread& t = *it;
			bool solutionFound = t.nextStep();
			if (solutionFound) {
				if (out != nullptr) {
					out->init(&t.m_epp);
					OGDF_ASSERT(out->isPlanar());
				}
				return Module::ReturnType::Feasible;
			}

			if (t.finished()) {
				it = m_threads.erase(it);
			} else {
				++it;
			}
		}
	}
	return Module::ReturnType::NoFeasibleSolution;
}

void OnePlanarityBacktracking::push(OnePlanarityBacktracking::DFSThread* t, EdgePairPartition* epp) {
	t->m_branchStack.push(false);
	if ((int)m_threads.size() < m_maxThreads) {
		m_threads.emplace_back(*epp, this);
	} else {
		t->m_branchStack.push(true);
	}
}

OnePlanarityBacktracking::NodeStatus OnePlanarityBacktracking::verifyNode(EdgePairPartition* epp) {
	m_processedNodes++;
	OnePlanarization pl(epp);
	if (pl.isPlanar()) {
		return NodeStatus::SOL;
	}

	for (auto& f : m_filters) {
		if (f->canCut(pl)) {
			return NodeStatus::CUT;
		}
	}
	OGDF_ASSERT(!OnePlanarization(epp).isPlanar());

	return NodeStatus::CNT;
}

OnePlanarityBacktracking::OnePlanarityBacktracking(int maxThreads, int maxKuratowskis)
	: m_maxThreads(maxThreads), m_maxExtractedKuratowskis(maxKuratowskis) {
	m_filters.emplace_back(new PlanarizationDensityFilter);
	m_filters.emplace_back(new SaturatedSubgraphPlanarityFilter);
	m_filters.emplace_back(new FreeEdgesDensityFilter);
	m_filters.emplace_back(new SeparatingCycleFilter);
}

OnePlanarityBacktracking::~OnePlanarityBacktracking() = default;
