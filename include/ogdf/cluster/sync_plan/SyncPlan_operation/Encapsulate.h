/** \file
 * \brief Utilities for the SyncPlan::encapsulate() operation.
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <ostream>

namespace ogdf::sync_plan::internal {

//! Information on a single block adjacent to a cut-vertex that is about to be encapsulated.
struct EncapsulatedBlock {
	node bicon = nullptr;
	node bicon_rep = nullptr;
	node star_rep = nullptr;
	PipeBij bij;

	explicit EncapsulatedBlock(node _bicon) : bicon(_bicon) { }

	friend std::ostream& operator<<(std::ostream& os, const EncapsulatedBlock& block);
};

}
