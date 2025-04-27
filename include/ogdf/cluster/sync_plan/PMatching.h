/** \file
 * \brief Manages the matching of P-nodes via \ref ogdf::sync_plan::Pipe "pipes" in an instance of SyncPlan.
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
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <functional>
#include <memory>
#include <ostream>
#include <utility>

namespace ogdf::sync_plan {
enum class PipeType { BlockBlock, BlockCut, CutCut };

//! A pair of matched vertices of the same degree, whose rotation shall be synchronized.
struct OGDF_EXPORT Pipe {
	node node1, node2;
	int pipe_priority = -1;
	List<Pipe>::iterator list_entry;
	void* heap_entry = nullptr;
	int heap_data = 0;
#ifdef OGDF_DEBUG
	int dbg_degree;
#endif

	Pipe(node node1, node node2);

	int degree() const;

	node getTwin(node n) const;

	friend std::ostream& operator<<(std::ostream& os, const Pipe& pipe);
};

//! A queue of all pipes, ordered by an arbitrary comparator function.
struct OGDF_EXPORT PipeQueue {
	virtual ~PipeQueue() = default;

	virtual bool empty() = 0;

	virtual int size() = 0;

	virtual Pipe* getTop() = 0;

	virtual void addPipe(Pipe* p) = 0;

	virtual void removePipe(Pipe* p) = 0;

	virtual void rebuild(List<Pipe>& pipes_list) = 0;

	virtual void clear() = 0;
};

//! Manages the matching of P-nodes via \ref Pipe "pipes" in an instance of SyncPlan.
class OGDF_EXPORT PMatching : protected GraphObserver {
	friend class SyncPlanConsistency;

private:
	int priority_pipes = 0;
	List<Pipe> pipes_list;
	NodeArray<Pipe*> nodes;
	std::unique_ptr<PipeQueue> queue;

public:
	explicit PMatching(const Graph* G) : GraphObserver(), nodes(*G, nullptr) { reregister(G); }

	bool isMatchedPVertex(node n) const;

	//! For a matched vertex, return the vertex it is matched with.
	node getTwin(node n) const;

	node getTwinOrNull(node n) const;

	const Pipe* getPipe(node n) const;

	const Pipe& getTopPipe() const;

	List<Pipe>::const_iterator begin() const;

	List<Pipe>::const_iterator end() const;

	const List<Pipe>& getPipes() const;

	int getPipeCount() const;

	//! Whether there are no pipes left.
	bool isReduced() const;

	//! Get the bijection between the edges incident to the two endpoints of a pipe.
	PipeBijRange getIncidentEdgeBijection(node u) const;

	void getIncidentEdgeBijection(node u, PipeBij& out) const;

	void getIncidentEdgeBijection(node u, AdjEntryArray<adjEntry>& out) const;

	void getIncidentEdgeBijection(node u, EdgeArray<edge>& out) const;

	adjEntry translateIncidentEdge(adjEntry e) const;

	std::function<std::ostream&(std::ostream&)> printBijection(node u) const;

	void matchNodes(node f, node s);

	node removeMatching(node n, node t = nullptr);

	//! Mark the pipe of this ndoe to be processed before all other pipes, no matter the order in the PipeQueue.
	void makePriority(node n);

	int getPriorityPipeCount() const { return priority_pipes; }

	//! Rebuild the PipeQueue, e.g. after priorities were changed externally.
	void rebuildHeap();

	void setPipeQueue(PipeQueue* q) {
		queue.reset(q);
		rebuildHeap();
	}

	void setPipeQueue(std::unique_ptr<PipeQueue> q) {
		queue = std::move(q);
		rebuildHeap();
	}

	PipeQueue* getPipeQueue() const { return queue.get(); }

protected:
	void nodeDeleted(node v) override;

	void nodeAdded(node v) override {};

	void edgeDeleted(edge e) override {};

	void edgeAdded(edge e) override {};

	void cleared() override {
		pipes_list.clear();
		if (queue) {
			queue->clear();
		}
		nodes.fill(nullptr);
	};
};
}
