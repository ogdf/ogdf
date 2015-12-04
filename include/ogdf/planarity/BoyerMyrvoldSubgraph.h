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

#pragma once

#include <random>
#include <ogdf/module/PlanarSubgraphModule.h>
#include <ogdf/internal/planarity/BoyerMyrvoldPlanar.h>
#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/planarity/BoothLueker.h>

namespace ogdf {

class OGDF_EXPORT BoyerMyrvoldSubgraph : public PlanarSubgraphModule
{
private:
	bool m_doMaximize;
	int m_runs;
	double m_randomness;
	BoothLueker m_planModule;
	std::minstd_rand m_rand;

public:
	/**
	 * Creates a new Boyer-Myrvold subgraph module.
	 *
	 * @param doMaximize iff true the returned subgraph will be maximal, but there will be no linear running time
	 * @param runs The number of times to start the algorithm, the greatest found subgraph is returned
	 * @param randomness If randomness equals 1, each edge is chosen with the same probability.
	 *                   A randomness of zero chooses each edge according to its cost.
	 *                   Any value between 0 and 1 is allowed and will result in a specific random influence.
	 *                   When performing multiple runs, a randomness greater zero should be chosen.
	 */
	BoyerMyrvoldSubgraph(bool doMaximize = true, int runs = 1, double randomness = 0)
	  : m_doMaximize(doMaximize),
	    m_runs(runs),
	    m_randomness(randomness),
	    m_rand(rand())
	{};

	~BoyerMyrvoldSubgraph() {};

	virtual BoyerMyrvoldSubgraph *clone() const override {
		return new BoyerMyrvoldSubgraph(m_doMaximize, m_runs);
	};

	//! Seeds the random generator for performing a random DFS.
	//! If this method is never called the random generator will be seeded by a value
	//! extracted from the global random generator.
	void seed(std::minstd_rand rand) {
		m_rand = rand;
	}

protected:
	/**
	 * Constructs a planar subgraph according to the options supplied to the constructor.
	 *
	 * @param graph the Graph to be planarized
	 * @param preferedEdges ignored
	 * @param delEdges will contain the list of edges that can be deleted to achieve planarity
	 * @param pCost the costs for removing each edge
	 * @param preferedImplyPlanar ignored
	 */
	virtual ReturnType doCall(
		const Graph &graph,
		const List<edge> &preferedEdges,
		List<edge> &delEdges,
		const EdgeArray<int> *pCost,
		bool preferedImplyPlanar) override;

	/**
	 * Inserts removed edges into the subgraph as long as it remains planar.
	 */
	void maximizeSubgraph(GraphCopy &graph);

	/**
	 * Returns true iff this edge could not be embedded.
	 */
	bool isRemoved(const GraphCopy &copy, const edge e) {
		return copy.copy(e) == nullptr || copy.copy(e)->source() != copy.copy(e->source()) || copy.copy(e)->target() != copy.copy(e->target());
	}
};

}
