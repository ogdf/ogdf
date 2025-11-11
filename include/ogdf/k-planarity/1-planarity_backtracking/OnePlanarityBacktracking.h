/** \file
 * \brief The main class of the 1-Planarity backtracker.
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

#pragma once

#include <ogdf/basic/Module.h>
#include <ogdf/basic/Timeouter.h>
#include <ogdf/basic/basic.h>
#include <ogdf/k-planarity/1-planarity_backtracking/EdgePairPartition.h>
#include <ogdf/k-planarity/1-planarity_backtracking/OnePlanarization.h> // IWYU pragma: keep

#include <list>
#include <memory>
#include <stack>
#include <vector>

namespace ogdf {
class Graph;
} // namespace ogdf

namespace ogdf::oneplan_backtracking {
class PartialSolutionFilter;

//! A backtracking implementation for the 1-Planarity problem.
/**
 * Implements the best configuration of the backtracking procedure for 1-Planarity evaluated in:
 * \remark Simon D. Fink, Miriam Münch, Matthias Pfretzschner and Ignaz Rutter. Heuristics for Exact
 * 1-Planarity Testing. To appear in the Proc. of the 33rd International Symposium on Graph Drawing
 * and Network Visualization, GD 2025, LIPIcs, Volume 357, 2025.
 *
 * This solver finds edge pairs to branch over by extracting multiple Kuratowski-Subdivisions of the
 * graph and uses filters to reject non-realizable partial solutions early to shrink the search space.
 * Moreover, it explores the search space using multiple depth-first searches running alternatingly
 * to increase the likelihood of finding solutions with few crossings early.
 */
class OGDF_EXPORT OnePlanarityBacktracking : public Timeouter {
private:
	//! A class representing one depth-first search in the backtracking tree.
	class DFSThread {
	public:
		//! The current partial solution considered by the DFSThread.
		EdgePairPartition m_epp;

		//! The associated solver.
		OnePlanarityBacktracking* m_solver;

		//! The DFS-stack. True indicates that new branches should be computed, false means that
		//! the changes in the partial solution should be reverted.
		std::stack<bool> m_branchStack;

		//! Creates a new DFS-thread associated with \p s that starts with a copy of \p x.
		DFSThread(EdgePairPartition& x, OnePlanarityBacktracking* s) : m_epp(x), m_solver(s) {
			m_branchStack.push(true);
		}

		//! Creates a new DFS-thread associated with \p s and initializes it with a new
		//! \ref EdgePairPartition for \p G with mode \p mode.
		DFSThread(const Graph& G, OneplanMode mode, OnePlanarityBacktracking* s)
			: m_epp(G, mode), m_solver(s) {
			m_branchStack.push(true);
		}

		//! Executes one step of the DFS and returns true if a solution was found.
		bool nextStep();

		//! Returns whether the DFS-thread has terminated.
		bool finished() const { return m_branchStack.empty(); }

		//! Computes the children of the current node in the search tree using Kuratowski-subdivisions.
		/**
		 *  If the maximum number of DFS-threads has not been reached, the children will be moved to
		 *  new DFS-threads.
		 *  @return true if a solution was found, false otherwise.
		 */
		bool branch();

		//! Proceeds to the next child if there are any left, otherwise reverts the changes to the partial solution.
		void cleanUp(DFSThread& thread);
	};

	//! Classification of partial solutions in the search tree.
	enum class NodeStatus {
		CNT, //!< Unknown whether partial solution is realizable or not.
		CUT, //!< Partial solution is not realizable and can be rejected.
		SOL //!< Partial solution is a solution.
	};

	std::list<DFSThread> m_threads;
	std::vector<std::unique_ptr<PartialSolutionFilter>> m_filters;

	int m_maxThreads;
	int m_maxExtractedKuratowskis;
	int m_processedNodes = 0;

public:
	//! Creates a new solver.
	/**
	 * @param maxThreads The maximum number of DFS-threads that are executed alternatingly.
	 * A higher number may help find solutions with few crossings faster, but comes at a cost of
	 * increased memory consumption.
	 * @param maxKuratowskis The maximum number of Kuratowski-Subdivisions that are taken into
	 * consideration when branching. A higher number may shrink the search space, but increases the
	 * time spent at each node in the search tree.
	 */
	explicit OnePlanarityBacktracking(int maxThreads = 1000, int maxKuratowskis = 1000);

	~OnePlanarityBacktracking();

	OnePlanarityBacktracking(const OnePlanarityBacktracking&) = delete;

	OnePlanarityBacktracking& operator=(const OnePlanarityBacktracking&) = delete;

	//! Tests whether \p G is 1-planar.
	/**
	 * A graph is 1-planar if it admits a drawing with at most one crossing per edge.
	 *
	 * @param G is the input graph.
	 * @param out is assigned a solution, if one is found.
	 * @return \ref Module::ReturnType::TimeoutInfeasible if a \ref timeLimit() was specified and
	 * exceeded. Otherwise returns \ref Module::ReturnType::Feasible if \p G is 1-planar and
	 * \ref Module::ReturnType::NoFeasibleSolution if it is not.
	 *
	 * \remark For better performance, first decompose the graph into its biconnected components.
	 */
	Module::ReturnType testOnePlanarity(const Graph& G, OnePlanarization* out = nullptr) {
		return test(OneplanMode::Normal, G, out);
	}

	//! Tests whether \p G is IC-Planar.
	/**
	 * A graph is IC-Planar if it admits a drawing with at most one crossing per edge and no two
	 * crossed edges may share an endpoint.
	 *
	 *
	 * @param G is the input graph.
	 * @param out is assigned a solution, if one is found.
	 * @return \ref Module::ReturnType::TimeoutInfeasible if a \ref timeLimit() was specified and
	 * exceeded. Otherwise returns \ref Module::ReturnType::Feasible if \p G is IC-planar and
	 * \ref Module::ReturnType::NoFeasibleSolution if it is not.
	 *
	 * \remark In contrast to 1-Planarity and NIC-Planarity, a graph is not necessarily IC-planar if
	 * and only if its biconnected components are.
	 */
	Module::ReturnType testICPlanarity(const Graph& G, OnePlanarization* out = nullptr) {
		return test(OneplanMode::IC, G, out);
	}

	//! Tests whether \p G is NIC-Planar.
	/**
	 * A graph is NIC-Planar if it admits a drawing with at most one crossing per edge and two
	 * pairs of crossed edges may have at most one common endpoint.
	 *
	 * @param G is the input graph.
	 * @param out is assigned a solution, if one is found.
	 * @return \ref Module::ReturnType::TimeoutInfeasible if a \ref timeLimit() was specified and
	 * exceeded. Otherwise returns \ref Module::ReturnType::Feasible if \p G is NIC-planar and
	 * \ref Module::ReturnType::NoFeasibleSolution if it is not.
	 *
	 * \remark For better performance, first decompose the graph into its biconnected components.
	 */
	Module::ReturnType testNICPlanarity(const Graph& G, OnePlanarization* out = nullptr) {
		return test(OneplanMode::NIC, G, out);
	}

	//! Returns the number of search tree nodes that were processed in the previous run.
	int processedNodes() const { return m_processedNodes; }

private:
	//! Runs the backtracking on \p G and writes the solution to \p out.
	Module::ReturnType test(OneplanMode mode, const Graph& G, OnePlanarization* out = nullptr);

	//! Determines whether a partial solution is a solution, non-realizable, or the search must continue.
	virtual NodeStatus verifyNode(EdgePairPartition* epp);

	//! Creates a new DFS-thread for \p epp if there is space, otherwise pushes it to \p t.
	virtual void push(DFSThread* t, EdgePairPartition* epp);
};
}
