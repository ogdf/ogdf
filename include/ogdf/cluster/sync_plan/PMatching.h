#pragma once

#include <ogdf/basic/AdjEntryArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/GraphObserver.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/cluster/sync_plan/basic/Iterators.h>
#include <ogdf/cluster/sync_plan/utils/Bijection.h>

#include <memory>
#include <ostream>

class PMatching;

enum class PipeType { BlockBlock, BlockCut, CutCut };

struct Pipe {
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

struct PipeQueue {
	virtual ~PipeQueue() = default;

	virtual bool empty() = 0;

	virtual int size() = 0;

	virtual Pipe* getTop() = 0;

	virtual void addPipe(Pipe* p) = 0;

	virtual void removePipe(Pipe* p) = 0;

	virtual void rebuild(List<Pipe>& pipes_list) = 0;

	virtual void clear() = 0;
};

class PMatching : protected GraphObserver {
	friend class PQPlanarityConsistency;

private:
	int priority_pipes = 0;
	List<Pipe> pipes_list;
	NodeArray<Pipe*> nodes;
	std::unique_ptr<PipeQueue> queue;

public:
	explicit PMatching(const Graph* G) : nodes(*G, nullptr), GraphObserver(G) { }

	bool isMatchedPVertex(node n) const;

	node getTwin(node n) const;

	node getTwinOrNull(node n) const;

	const Pipe* getPipe(node n) const;

	const Pipe& getTopPipe() const;

	List<Pipe>::const_iterator begin() const;

	List<Pipe>::const_iterator end() const;

	const List<Pipe>& getPipes() const;

	int getPipeCount() const;

	bool isReduced() const;

	PipeBijIterator getIncidentEdgeBijection(node u) const;

	void getIncidentEdgeBijection(node u, PipeBij& out) const;

	void getIncidentEdgeBijection(node u, AdjEntryArray<adjEntry>& out) const;

	void getIncidentEdgeBijection(node u, EdgeArray<edge>& out) const;

	adjEntry translateIncidentEdge(adjEntry e) const;

	std::function<std::ostream&(std::ostream&)> printBijection(node u) const;

	void matchNodes(node f, node s);

	node removeMatching(node n, node t = nullptr);

	void makePriority(node n);

	int getPriorityPipeCount() const { return priority_pipes; }

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
