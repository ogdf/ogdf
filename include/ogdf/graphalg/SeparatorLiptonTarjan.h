/** \file
 * \brief Declaration of class SeparatorLiptonTarjan.
 *
 * \author Thomas Klein
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

#include <ogdf/graphalg/PlanarSeparatorModule.h>

namespace ogdf {

//! Computes planar separators according to Lipton and Tarjan 1979
/**
 * Computes separators according to the paper
 * "A Separator Theorem for Planar Graphs" by Richard J. Lipton and Robert Endre Tarjan,
 * published in the SIAM Journal on Applied Mathematics, 1979.
 *
 * @ingroup ga-plansep
 */
class OGDF_EXPORT SeparatorLiptonTarjan : public PlanarSeparatorModule {
	friend class Cycle;

public:
	/**
	 * Constructor.
	 *
	 * @param useTriangulatingBFS whether to use triangulating BFS or not
	 * @param treeHeightIt how many iterations of tree height maximization to perform
	 */
	SeparatorLiptonTarjan(bool useTriangulatingBFS = false, unsigned int treeHeightIt = 0)
		: useTriBFS {useTriangulatingBFS}, treeHeightIterations(treeHeightIt + 1) { }

	virtual double getMaxSeparatorSize(int n) const override { return sqrt(8) * sqrt(n); }


protected:
	bool useTriBFS;
	unsigned int treeHeightIterations;

	virtual bool doSeparate(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second) override;

	std::shared_ptr<BFSTreeClassical> tree;

	virtual std::string getSpecificName() const override {
		std::string name = "LT";
		if (useTriBFS) {
			name += "-triBFS";
		}
		if (treeHeightIterations > 1) {
			name += "-THM-" + std::to_string(treeHeightIterations - 1);
		}

		return name;
	}

	/**
	 * Creates the BFS tree used by the algorithm.
	 */
	virtual void makeTree();

	/**
	 * Fills the lists with the cycle / inside / outside once the cycle is ready.
	 *
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 */
	void fillLists(List<node>& separator, List<node>& first, List<node>& second) const;

	/**
	 * Chooses the initial edge for the very first cycle.
	 *
	 * @return a random edge
	 */
	edge chooseEdge() const;
};


}
