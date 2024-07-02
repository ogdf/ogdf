/** \file
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
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/PipeOrder.h>
#include <ogdf/cluster/sync_plan/SyncPlan.h>
#include <ogdf/cluster/sync_plan/SyncPlanComponents.h>

using namespace ogdf::sync_plan::internal;

namespace ogdf::sync_plan {

bool PipeQueueByDegreePreferContract::comparePipes(const Pipe* x, const Pipe* y) const {
	if (x->heap_data != y->heap_data) {
		return x->heap_data < y->heap_data;
	}

	if (invert_degree) {
		return x->degree() < y->degree();
	} else {
		return x->degree() > y->degree();
	}
}

bool PipeQueueByDegreePreferContract::isQueue1(Pipe* p) const {
	if (p->pipe_priority >= 0) {
		return true;
	}
	bool ret = p->degree() <= 3 || PQ->canContract(p);
	if (invert_contract) {
		ret = !ret;
	}
	return ret;
}

}
