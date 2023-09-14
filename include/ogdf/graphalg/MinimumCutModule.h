/** \file
 * \brief Declaration of ogdf::MinimumCutModule
 *
 * \author Max Ilsen
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

#include <ogdf/basic/Graph_d.h>

namespace ogdf {

/**
 * Serves as an interface for various methods to compute minimum cuts with or
 * without edge weights.
 *
 * @tparam T The type of the edge weights of the mincut instance
 */
template<typename T>
class MinimumCutModule {
public:
	//! Do nothing on destruction
	virtual ~MinimumCutModule() { }

	/**
	 * Computes the minimum cut of \p G.
	 *
	 * @param G The input graph
	 * @return The minimum cut value
	 */
	virtual T call(const Graph& G) = 0;

	/**
	 * Computes the minimum cut of \p G with edge weights \p weights.
	 *
	 * @param G The input graph
	 * @param weights The edge weights
	 * @return The minimum cut value
	 */
	virtual T call(const Graph& G, const EdgeArray<T>& weights) = 0;

	//! Returns the edges defining the computed mincut.
	virtual const ArrayBuffer<edge>& edges() = 0;

	//! Returns a list of nodes belonging to one side of the bipartition.
	virtual const ArrayBuffer<node>& nodes() = 0;

	//! Returns the value of the last minimum cut computation.
	virtual T value() const = 0;
};

}
