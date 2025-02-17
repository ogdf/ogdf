/** \file
 * \brief Provides a simple, dfs-based algorithm for biconnectivity augmentation.
 *
 * \author Carsten Gutwenger
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>
#include <ogdf/augmentation/AugmentationModule.h>

namespace ogdf {
template<class E>
class List;

/**
 * \brief Implementation of a DFS-based algorithm for biconnectivity augmentation.
 *
 * @ingroup ga-augment
 *
 * The class DfsMakeBiconnected implements an augmentation algorithms
 * that augments a graph to a biconnected graph. In addition, if the graph was
 * planar before augmentation, the resulting graph will be biconnected and
 * planar.
 * The algorithm simply uses DFS and, whenever a cut vertex is discovered,
 * a new edge is added.
 */

class OGDF_EXPORT DfsMakeBiconnected : public AugmentationModule {
public:
	//! Creates an instance of DFS-based biconnectivity augmentation.
	DfsMakeBiconnected() { }

	//! Destruction
	~DfsMakeBiconnected() { }

protected:
	//! Implements the algorithm call.
	virtual void doCall(Graph& G, List<edge>& L) override;
};

}
