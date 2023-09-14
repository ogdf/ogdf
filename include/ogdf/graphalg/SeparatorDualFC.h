/** \file
 * \brief Declaration of class SeparatorDualFC, which applies the Fundamental Cycle Lemma directly to obtain a cycle.
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
#include <ogdf/graphalg/SeparatorDual.h>
#include <ogdf/graphalg/SeparatorLiptonTarjanFC.h>
#include <ogdf/graphalg/planar_separator/SeparatorDualHelper.h>

namespace ogdf {

//! Computes planar separators by applying the Fundamental Cycle Lemma directly, without trying tree levels first.
/**
 * @ingroup ga-plansep
 */
class OGDF_EXPORT SeparatorDualFC : public SeparatorLiptonTarjanFC {
public:
	/**
	 * Constructor.
	 *
	 * @param useTriBFS whether to use triangulating BFS or not
	 */
	SeparatorDualFC(bool useTriBFS = false) : useTriangulatingBFS {useTriBFS} { }

	/** Maximum separator size depends on diameter of graph, so returns -1. */
	virtual double getMaxSeparatorSize(int n) const override { return -1; }

protected:
	bool useTriangulatingBFS;
	std::shared_ptr<ArrayBFSTree> tree;

	/**
	 * Builds the BFS tree.
	 */
	void makeTree();

	virtual bool doSeparate(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second) override;

	/**
	 * Finds a suitable cycle by performing a DFS over the faces of the dual of the graph.
	 *
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @return true on success
	 */
	virtual bool findCycle(List<node>& separator, List<node>& first, List<node>& second) override;

	virtual std::string getSpecificName() const override {
		std::string name = "DualFC";
		if (useTriangulatingBFS) {
			name += "-triBFS";
		}
		return name;
	}
};

}
