/** \file
 * \brief Consistency checks for debugging the SyncPlan algorithm.
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

#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/SyncPlanDrawer.h>

#include <string>

namespace ogdf::sync_plan {
class SyncPlan;

//! Consistency checks for debugging the SyncPlan algorithm.
class OGDF_EXPORT SyncPlanConsistency {
	SyncPlan& pq;
	SyncPlanDrawer draw;
	int checkCounter = 0;

public:
	static bool doWriteOut;

	explicit SyncPlanConsistency(SyncPlan& _pq) : pq(_pq), draw(&_pq) {};

	bool consistencyCheck(bool force_check_components = false);

	void writeOut(std::string name = "", bool format = true, bool components = true);

	void checkComponentRegeneration();

	int getCheckCounter() const { return checkCounter; }
};
}
