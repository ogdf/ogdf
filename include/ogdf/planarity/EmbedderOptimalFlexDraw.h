/*
 * $Revision: 3835 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 13:18:01 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief The algorithm computes a planar embedding with minimum cost.
 *
 * \author Dzmitry Sledneu
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 * \par
 * Copyright (C)\n
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation
 * and appearing in the files LICENSE_GPL_v2.txt and
 * LICENSE_GPL_v3.txt included in the packaging of this file.
 *
 * \par
 * In addition, as a special exception, you have permission to link
 * this software with the libraries of the COIN-OR Osi project
 * (http://www.coin-or.org/projects/Osi.xml), all libraries required
 * by Osi, and all LP-solver libraries directly supported by the
 * COIN-OR Osi project, and distribute executables, as long as
 * you follow the requirements of the GNU General Public License
 * in regard to all of the software in the executable aside from these
 * third-party libraries.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_EMBEDDER_OPTIMAL_FLEX_DRAW_H
#define OGDF_EMBEDDER_OPTIMAL_FLEX_DRAW_H

#include <ogdf/module/EmbedderModule.h>
#include <ogdf/basic/ModuleOption.h>
#include <ogdf/module/MinCostFlowModule.h>
#include <ogdf/decomposition/StaticPlanarSPQRTree.h>

namespace ogdf {

/// The algorithm computes a planar embedding with minimum cost.
///
/// See paper "Optimal Orthogonal Graph Drawing with Convex Bend Costs"
/// by Thomas Blasius, Ignaz Rutter, Dorothea Wagner (2012) for details.
class OGDF_EXPORT EmbedderOptimalFlexDraw : public EmbedderModule
{
public:

	EmbedderOptimalFlexDraw();

	~EmbedderOptimalFlexDraw() { }

	void call(Graph &G, adjEntry &adjExternal);

	/// Sets the module option to compute min-cost flow.
	void setMinCostFlowComputer(MinCostFlowModule *pMinCostFlowComputer) {
		m_minCostFlowComputer.set(pMinCostFlowComputer);
	}

	/// Sets bend costs for each edge.
	void cost(EdgeArray<int> *cost) {
		m_cost = cost;
	}

private:

	ModuleOption<MinCostFlowModule> m_minCostFlowComputer;

	EdgeArray<int> *m_cost;

	void createNetwork(
			node parent,
		node mu,
		int bends,
		NodeArray<int> cost[],
		NodeArray<long long> embedding[],
		Skeleton &skeleton,
		EdgeArray<node> &edgeNode,
		Graph &N,
		EdgeArray<int> &upper,
		EdgeArray<int> &perUnitCost,
		NodeArray<int> &supply);

	void optimizeOverEmbeddings(
		StaticPlanarSPQRTree &T,
		node parent,
		node mu,
		int bends,
		NodeArray<int> cost[],
		NodeArray<long long> embedding[]);

	void computePrincipalSplitComponentCost(
		StaticPlanarSPQRTree &T,
		NodeArray<int> cost[],
		NodeArray<long long> embedding[],
		node parent,
		node mu);
};

} // end namespace ogdf

#endif
