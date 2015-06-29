/** \file
 * \brief Declaration of the subgraph wrapper class of the Boyer-Myrvold planarity test
 *
 * \author Tilo Wiedera
 *
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
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
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#include <ogdf/planarity/BoyerMyrvoldSubgraph.h>

namespace ogdf {

Module::ReturnType BoyerMyrvoldSubgraph::doCall(
	const Graph &graph,
	const List<edge> &preferedEdges,
	List<edge> &delEdges,
	const EdgeArray<int>  *pCosts,
	bool preferedImplyPlanar)
{
	int bestCost = -1;

	for(int i = 0; i < m_runs; i++) {
		SListPure<KuratowskiStructure> tmp;
		GraphCopy copy(graph);
		EdgeArray<int> *costs = nullptr;

		if(pCosts != nullptr) {
			costs = new EdgeArray<int>(copy);

			for(edge e : copy.edges) {
				(*costs)[e] = (*pCosts)[copy.original(e)];
			}
		}

		BoyerMyrvoldPlanar bmp(copy, false, BoyerMyrvoldPlanar::doFindUnlimited, false, tmp, m_randomness, true, true, costs);
		std::minstd_rand rand(m_rand());
		bmp.seed(rand);
		bmp.start();

		OGDF_ASSERT(m_planModule.isPlanar(copy));
		OGDF_ASSERT(copy.numberOfEdges() == graph.numberOfEdges());

		if(m_doMaximize) {
			maximizeSubgraph(copy);
		}

		int totalCost = 0;
		if(i != 0) {
			for(edge e : graph.edges) {
				if(isRemoved(copy, e)) {
					totalCost += costs == nullptr ? 1 : (*pCosts)[e];
				}
			}
		}

		if(i == 0 || totalCost < bestCost) {
			bestCost = totalCost;
			delEdges.clear();
			for(edge e : graph.edges) {
				if(isRemoved(copy, e)) {
					delEdges.pushBack(e);
				}
			}
		}

		if(costs != nullptr) {
			delete costs;
		}
	}

	return Module::ReturnType::retFeasible;
}

void BoyerMyrvoldSubgraph::maximizeSubgraph(GraphCopy &copy) {
	const Graph &graph = copy.original();
	List<edge> removedEdges;

	for(edge e : graph.edges) {
		if(isRemoved(copy, e)) {
			removedEdges.pushBack(e);
			copy.delEdge(copy.copy(e));
		}
	}

	bool inserted = true;
	while(inserted) {
		inserted = false;
		removedEdges.permute();

		for(ListIterator<edge> it = removedEdges.begin(); it != removedEdges.end();) {
			ListIterator<edge> tmp = it;
			it = it.succ();
			edge e = copy.newEdge(*tmp);

			if(m_planModule.isPlanar(copy)) {
				inserted = true;
				removedEdges.del(tmp);
			} else {
				copy.delEdge(e);
			}
		}
	}
}

}
