/** \file
 * \brief Applies the node coloring approximation specified by Boppana&Halldorsson.
 *
 * \author Jan-Niklas Buckow
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
#include <ogdf/graphalg/NodeColoringModule.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Boppana&Halldorsson
 * which colors the graph by finding independent sets with the Ramsey-algorithm.
 */
class OGDF_EXPORT NodeColoringBoppanaHalldorsson : public NodeColoringModule {
public:
	/**
	 * The constructor.
	 * Initializes the Ramsey-procedure with the smallest index procedure.
	 */
	NodeColoringBoppanaHalldorsson() {
		NodeColoringModule::m_ramseyProcedure = NodeColoringModule::RamseyProcedure::smallestDegree;
	}

	/**
	 * Sets the Ramsey-procedure of findings nodes to a specific value.
	 * @param ramseyProcedure The given Ramsey-procedure.
	 */
	inline void setRamseyProcedure(RamseyProcedure ramseyProcedure) {
		NodeColoringModule::m_ramseyProcedure = ramseyProcedure;
	}

	virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors,
			NodeColor start = 0) override;
};
}
