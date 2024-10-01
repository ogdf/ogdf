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
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>
#include <ogdf/cluster/sync_plan/utils/Logging.h>

#include <functional>
#include <memory>
#include <ostream>

namespace ogdf::sync_plan {
using namespace ogdf::sync_plan::internal;

Pipe::Pipe(node _node1, node _node2)
	: node1(_node1)
	, node2(_node2)
#ifdef OGDF_DEBUG
	, dbg_degree(_node1->degree())
#endif
{
}

int Pipe::degree() const {
	OGDF_ASSERT(dbg_degree == node1->degree());
	OGDF_ASSERT(dbg_degree == node2->degree());
	return node1->degree();
}

node Pipe::getTwin(node n) const {
	if (n == node1) {
		return node2;
	} else {
		OGDF_ASSERT(n == node2);
		return node1;
	}
}

std::ostream& operator<<(std::ostream& os, const Pipe& pipe) {
	os << "(°" << pipe.degree() << " #" << pipe.node1->index() << " = #" << pipe.node2->index()
	   << ")";
	return os;
}

bool PMatching::isMatchedPVertex(node n) const {
	OGDF_ASSERT(n != nullptr);
	return nodes[n] != nullptr;
}

const Pipe* PMatching::getPipe(node n) const {
	OGDF_ASSERT(n != nullptr);
	return nodes[n];
}

node PMatching::getTwin(node n) const {
	OGDF_ASSERT(n != nullptr);
	Pipe* pipe = nodes[n];
	OGDF_ASSERT(pipe != nullptr);
	node twin = pipe->getTwin(n);
	OGDF_ASSERT(nodes[twin] == pipe);
	OGDF_ASSERT(twin->degree() == n->degree());
	return twin;
}

node PMatching::getTwinOrNull(node n) const {
	OGDF_ASSERT(n != nullptr);
	Pipe* pipe = nodes[n];
	if (pipe == nullptr) {
		return nullptr;
	}
	node twin = pipe->getTwin(n);
	OGDF_ASSERT(nodes[twin] == pipe);
	OGDF_ASSERT(twin->degree() == n->degree());
	return twin;
}

const Pipe& PMatching::getTopPipe() const {
	OGDF_ASSERT(queue);
	return *(queue->getTop());
}

void PMatching::matchNodes(node f, node s) { // room for improvement: don't match deg-1 nodes
	OGDF_ASSERT(!isMatchedPVertex(f));
	OGDF_ASSERT(!isMatchedPVertex(s));
	OGDF_ASSERT(f->degree() == s->degree());
	// OGDF_ASSERT(f->degree() > 3); // this is violated by the RemCut before and Contract Rays after Propagate
	List<Pipe>::iterator it = pipes_list.emplaceBack(f, s);
	(*it).list_entry = it;
	if (queue) {
		queue->addPipe(&(*it));
	}
	nodes[f] = nodes[s] = &(*it);
}

node PMatching::removeMatching(node n, node t) {
	OGDF_ASSERT(n != nullptr);
	Pipe* pipe = nodes[n];
	OGDF_ASSERT(pipe != nullptr);
	if (pipe->pipe_priority >= 0) {
		priority_pipes--;
	}
	node twin = pipe->getTwin(n);
	OGDF_ASSERT(t == nullptr || t == twin);
	nodes[n] = nodes[twin] = nullptr;
	if (queue) {
		queue->removePipe(pipe);
	}
	pipes_list.del(pipe->list_entry);
	return twin;
}

void PMatching::makePriority(node n) {
	Pipe* p = nodes[n];
	bool incr = p->pipe_priority < 0;
	if (queue) {
		p->pipe_priority = queue->getTop()->pipe_priority + 1;
		queue->removePipe(p);
		queue->addPipe(p);
	} else {
		p->pipe_priority = priority_pipes;
	}
	if (incr) {
		priority_pipes++;
	}
}

void PMatching::rebuildHeap() {
	if (queue) {
		queue->rebuild(pipes_list);
	}
}

void PMatching::nodeDeleted(node v) { OGDF_ASSERT(!isMatchedPVertex(v)); }

List<Pipe>::const_iterator PMatching::begin() const { return pipes_list.begin(); }

List<Pipe>::const_iterator PMatching::end() const { return pipes_list.end(); }

const List<Pipe>& PMatching::getPipes() const { return pipes_list; }

int PMatching::getPipeCount() const { return pipes_list.size(); }

bool PMatching::isReduced() const { return pipes_list.empty(); }

PipeBijRange PMatching::getIncidentEdgeBijection(node u) const {
	return getPipeBijection(u, getTwin(u));
}

void PMatching::getIncidentEdgeBijection(node u, PipeBij& out) const {
	getPipeBijection(u, getTwin(u), out);
}

void PMatching::getIncidentEdgeBijection(node u, AdjEntryArray<adjEntry>& out) const {
	getPipeBijection(u, getTwin(u), out);
}

void PMatching::getIncidentEdgeBijection(node u, EdgeArray<edge>& out) const {
	getPipeBijection(u, getTwin(u), out);
}

adjEntry PMatching::translateIncidentEdge(adjEntry e) const {
	node u = e->theNode();
	node v = getTwin(u);
	auto ve = v->adjEntries.rbegin();
	for (auto ue : u->adjEntries) {
		if (ue == e) {
			return (*ve);
		} else {
			ve++;
		}
		OGDF_ASSERT(ve != v->adjEntries.rend());
	}
	OGDF_ASSERT(false);
	return nullptr;
}

std::function<std::ostream&(std::ostream&)> PMatching::printBijection(node u) const {
	node v = getTwin(u);
	const auto bij = getIncidentEdgeBijection(u);
	return [u, v, bij](std::ostream& ss) -> std::ostream& {
		return ss << "u" << u->index() << " = v" << v->index() << ": "
				  << sync_plan::internal::printBijection(bij);
	};
}

}
