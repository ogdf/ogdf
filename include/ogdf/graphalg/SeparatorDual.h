/** \file
 * \brief Declaration of class SeparatorDual.
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

#include <ogdf/basic/FaceArray.h>
#include <ogdf/graphalg/PlanarSeparatorModule.h>
#include <ogdf/graphalg/SeparatorLiptonTarjan.h>
#include <ogdf/graphalg/SeparatorLiptonTarjanFC.h>
#include <ogdf/graphalg/planar_separator/SeparatorDualHelper.h>

#include <unordered_set>

namespace ogdf {

//! Computes planar separators using the Dual of the graph.
/**
 * Computes planar separators using the dual of the graph as presented
 * in the textbook "The Design and Analysis of Algorithms" by D. Kozen, Springer, 1992.
 *
 * @ingroup ga-plansep
 */
class OGDF_EXPORT SeparatorDual : public SeparatorLiptonTarjan {
public:
	/**
	 * Constructor.
	 *
	 * @param useTriangulatingBFS whether to use triangulating BFS for the second phase or not
	 * @param treeHeightIt how many iterations of tree height maximization to perform in the first phase
	 */
	SeparatorDual(bool useTriangulatingBFS = false, unsigned int treeHeightIt = 0)
		: useTriBFS {useTriangulatingBFS}, treeHeightIterations(treeHeightIt + 1) { }

	virtual std::string getSpecificName() const override {
		std::string name = "Dual";
		if (useTriBFS) {
			name += "-triBFS";
		}
		if (treeHeightIterations > 1) {
			name += "-THM-" + std::to_string(treeHeightIterations - 1);
		}

		return name;
	}

	virtual double getMaxSeparatorSize(int n) const override { return 4 * sqrt(n); }

protected:
	virtual void makeTree() override;

	virtual bool doSeparate(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second) override;

	bool useTriBFS; // whether to use triangulating BFS for the second phase or not
	unsigned int treeHeightIterations; // how many iterations of tree height maximisation to run
};

}
