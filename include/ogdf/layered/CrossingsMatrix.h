/** \file
 * \brief Declaration of class CrossingsMatrix.
 *
 * \author Andrea Wagner
 *         Michael Schulz
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Array2D.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>

#include <cstdint>

namespace ogdf {
class HierarchyLevels;
class Level;

//! Implements crossings matrix which is used by some
//! TwoLayerCrossingMinimization heuristics (e.g. split)
class OGDF_EXPORT CrossingsMatrix {
public:
	CrossingsMatrix() : matrix(0, 0, 0, 0) { m_bigM = 10000; }

	explicit CrossingsMatrix(const HierarchyLevels& levels);

	~CrossingsMatrix() { }

	int operator()(int i, int j) const { return matrix(map[i], map[j]); }

	void swap(int i, int j) { map.swap(i, j); }

	//! ordinary init
	void init(Level& L);

	//! SimDraw init
	void init(Level& L, const EdgeArray<uint32_t>* edgeSubGraphs);

private:
	Array<int> map;
	Array2D<int> matrix;
	//! need this for SimDraw to grant epsilon-crossings instead of zero-crossings
	int m_bigM; // is set to some big number in both constructors
};

}
