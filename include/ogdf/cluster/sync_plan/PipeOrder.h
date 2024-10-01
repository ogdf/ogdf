/** \file
 * \brief Different PipeQueue implementations with varying sorting functions.
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

#include <ogdf/basic/List.h>
#include <ogdf/basic/PriorityQueue.h>
#include <ogdf/basic/basic.h>
#include <ogdf/cluster/sync_plan/PMatching.h>

#include <algorithm>
#include <memory>
#include <random>

namespace ogdf::sync_plan {
// using PipeCmp = std::function<bool(const Pipe *, const Pipe *)>;

//! A null-safe and priority aware comparator (wrapper) for \ref Pipe "pipes".
template<typename PipeCmp>
struct PipeCmpPtr {
	const PipeCmp* cmp;

	PipeCmpPtr() = delete;

	PipeCmpPtr(const PipeCmp* queue) : cmp(queue) {};

	bool operator()(const Pipe* x, const Pipe* y) const {
		if (x == nullptr) {
			return y != nullptr; // true
		} else if (y == nullptr) {
			return false;
		} else if (x->pipe_priority != y->pipe_priority) {
			return x->pipe_priority > y->pipe_priority;
		} else {
			return cmp->comparePipes(x, y);
		}
	}

	bool checkOrder(const Pipe* first, const Pipe* second) const {
		return first == second || !(*this)(second, first);
	}
};

//! PipeQueue CRTP base class for ordering \ref Pipe "pipes" by some simple comparator function.
template<typename PipeCmp>
class SimplePipeQueue : public PipeQueue {
public:
	using PipesHeap = PriorityQueue<Pipe*, PipeCmpPtr<PipeCmp>>;
	using PipesHeapHandle = typename PipesHeap::handle;

protected:
	std::unique_ptr<PipesHeap> pipes_heap;

public:
	SimplePipeQueue() = default;

	SimplePipeQueue(const SimplePipeQueue& copy) = delete;

	SimplePipeQueue(SimplePipeQueue&& move) = delete;

	SimplePipeQueue& operator=(const SimplePipeQueue& copy) = delete;

	SimplePipeQueue& operator=(SimplePipeQueue&& move) = delete;

	bool empty() override { return pipes_heap->empty(); }

	int size() override { return pipes_heap->size(); }

	Pipe* getTop() override { return pipes_heap->top(); }

	void addPipe(Pipe* p) override {
#ifdef OGDF_DEBUG
		OGDF_ASSERT(p->node1->degree() == p->node2->degree());
		p->dbg_degree = p->node1->degree();
#endif
		p->heap_entry = pipes_heap->push(p);
	}

	void removePipe(Pipe* pipe) override {
		OGDF_ASSERT(pipes_heap->value((PipesHeapHandle)pipe->heap_entry) == pipe);
		pipes_heap->decrease((PipesHeapHandle)pipe->heap_entry, nullptr);
		OGDF_ASSERT(pipes_heap->top() == nullptr);
		pipes_heap->pop();
		OGDF_ASSERT(pipes_heap->empty() || pipes_heap->top() != nullptr);
	}

	void rebuild(List<Pipe>& pipes_list) override {
		clear();
		for (Pipe& p : pipes_list) {
			addPipe(&p);
		}
	}

	void clear() override {
#if 0
		Pipe* top = nullptr;
		while (!pipes_heap->empty()) {
			OGDF_ASSERT(checkOrder(top, pipes_heap->top()));
			top = pipes_heap->top();
			pipes_heap->pop();
		}
#else
		pipes_heap->clear();
#endif
	}
};

//! PipeQueue yielding pipes in order of descending or ascending degree.
struct OGDF_EXPORT PipeQueueByDegree : public SimplePipeQueue<PipeQueueByDegree> {
	bool invert_degree;

	/**
	 * @param invert if true, will order by ascending degree; if false (the default), will order by descending degree.
	 */
	explicit PipeQueueByDegree(bool invert = false) : invert_degree(invert) {
		pipes_heap = std::make_unique<PipesHeap>(this);
	};

	bool comparePipes(const Pipe* x, const Pipe* y) const {
		if (invert_degree) {
			return x->degree() < y->degree();
		} else {
			return x->degree() > y->degree();
		}
	}
};

//! PipeQueue yielding pipes in some random (but stable and deterministic) order.
struct OGDF_EXPORT PipeQueueRandom : public SimplePipeQueue<PipeQueueRandom> {
	using engine = std::minstd_rand;
	mutable engine gen;

	explicit PipeQueueRandom() { pipes_heap = std::make_unique<PipesHeap>(this); };

	engine::result_type hash(const Pipe* p) const {
		gen.seed(max(p->node1->index(), p->node2->index()));
		gen.discard(5 + p->degree() % 20);
		return gen();
	}

	bool comparePipes(const Pipe* x, const Pipe* y) const {
		if (x == nullptr) {
			return true;
		} else if (y == nullptr) {
			return false;
		}
		return hash(x) > hash(y);
	}
};

//! Base class for PipeQueues providing a "priority lane" for some pipes and sorting with different functions in both lanes/queues.
template<typename PipeCmp1, typename PipeCmp2>
class DoublePipeQueue : public SimplePipeQueue<PipeCmp1> {
public:
	using PipesHeap2 = PriorityQueue<Pipe*, PipeCmpPtr<PipeCmp2>>;
	using PipesHeapHandle2 = typename PipesHeap2::handle;
	using Base = SimplePipeQueue<PipeCmp1>;

protected:
	using Base::pipes_heap;
	std::unique_ptr<PipesHeap2> pipes_heap2;

	virtual bool isQueue1(Pipe* p) const = 0;

public:
	bool empty() override { return Base::empty() && pipes_heap2->empty(); }

	int size() override { return Base::size() + pipes_heap2->size(); }

	Pipe* getTop() override {
		while (!Base::empty()) {
			Pipe* p = Base::getTop();
			if (isQueue1(p)) {
				return p;
			} else {
				removePipe(p);
				addPipe(p);
			}
		}
		OGDF_ASSERT(!pipes_heap2->empty());
		return pipes_heap2->top();
	}

	void addPipe(Pipe* p) override {
#ifdef OGDF_DEBUG
		OGDF_ASSERT(p->node1->degree() == p->node2->degree());
		p->dbg_degree = p->node1->degree();
#endif
		if (isQueue1(p)) {
			p->heap_data = 1;
			p->heap_entry = pipes_heap->push(p);
		} else {
			p->heap_data = 2;
			p->heap_entry = pipes_heap2->push(p);
		}
	}

	void removePipe(Pipe* pipe) override {
		if (pipe->heap_data == 2) {
			OGDF_ASSERT(pipes_heap2->value((PipesHeapHandle2)pipe->heap_entry) == pipe);
			pipes_heap2->decrease((PipesHeapHandle2)pipe->heap_entry, nullptr);
			OGDF_ASSERT(pipes_heap2->top() == nullptr);
			pipes_heap2->pop();
			OGDF_ASSERT(pipes_heap2->empty() || pipes_heap2->top() != nullptr);
		} else {
			Base::removePipe(pipe);
		}
	}

	void clear() override {
#if 0
		Pipe* top = nullptr;
		while (!pipes_heap->empty()) {
			OGDF_ASSERT(checkOrder(top, pipes_heap->top()));
			top = pipes_heap->top();
			pipes_heap->pop();
		}
		while (!pipes_heap2->empty()) {
			OGDF_ASSERT(checkOrder(top, pipes_heap2->top()));
			top = pipes_heap2->top();
			pipes_heap2->pop();
		}
#else
		pipes_heap->clear();
		pipes_heap2->clear();
#endif
	}
};

class SyncPlan;

//! PipeQueue yielding contractable pipes first (or last), in order of descending (or ascending) degree.
struct OGDF_EXPORT PipeQueueByDegreePreferContract
	: public DoublePipeQueue<PipeQueueByDegreePreferContract, PipeQueueByDegreePreferContract> {
	SyncPlan* PQ;
	bool invert_degree, invert_contract;

	/**
	 * @param pq the SyncPlan instance we're working on.
	 * @param invertDegree if true, will order by ascending degree; if false (the default), will order by descending degree.
	 * @param invertContract if true, will yield non-contractable pipes first; if false (the default), will yield contractable pipes first.
	 */
	explicit PipeQueueByDegreePreferContract(SyncPlan* pq, bool invertDegree = false,
			bool invertContract = false)
		: PQ(pq), invert_degree(invertDegree), invert_contract(invertContract) {
		pipes_heap = std::make_unique<PipesHeap>(this);
		pipes_heap2 = std::make_unique<PipesHeap2>(this);
	}

	bool comparePipes(const Pipe* x, const Pipe* y) const;

	bool isQueue1(Pipe* p) const override;
};
}
