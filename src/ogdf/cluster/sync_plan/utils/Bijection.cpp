/** \file
 * \brief TODO Document
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
#include <ogdf/basic/AdjEntryArray.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

PipeBijIterator getPipeBijection(node u, node v) {
	OGDF_ASSERT(u->degree() == v->degree());
	return PipeBijIterator(u->adjEntries.begin(), u->adjEntries.end(), v->adjEntries.rbegin(),
			v->adjEntries.rend());
}

void getPipeBijection(node u, node v, PipeBij& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out.emplaceBack(u_adj, (*v_adj_it));
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getPipeBijection(node u, node v, AdjEntryArray<adjEntry>& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->graphOf() == out.graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out[u_adj] = *v_adj_it;
		out[*v_adj_it] = u_adj;
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getPipeBijection(node u, node v, EdgeArray<edge>& out) {
	OGDF_ASSERT(u->graphOf() == v->graphOf());
	OGDF_ASSERT(u->graphOf() == out.graphOf());
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		out[u_adj->theEdge()] = (*v_adj_it)->theEdge();
		v_adj_it++;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void getFrozenPipeBijection(node u, node v, FrozenPipeBij& bij) {
	OGDF_ASSERT(u->degree() == v->degree());
	auto v_adj_it = v->adjEntries.rbegin();
	for (adjEntry u_adj : u->adjEntries) {
		OGDF_ASSERT(v_adj_it != v->adjEntries.rend());
		bij.emplaceBack(u_adj->theEdge()->index(), (*v_adj_it)->theEdge()->index());
		++v_adj_it;
	}
	OGDF_ASSERT(v_adj_it == v->adjEntries.rend());
}

void freezePipeBijection(const PipeBij& in, FrozenPipeBij& out) {
	for (const PipeBijPair& pair : in) {
		out.emplaceBack(pair.first->theEdge()->index(), pair.second->theEdge()->index());
	}
}
