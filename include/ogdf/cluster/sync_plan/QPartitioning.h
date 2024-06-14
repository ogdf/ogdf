#pragma once

#include <ogdf/basic/GraphObserver.h>

#include <utils/RegisteredArray.h>

using namespace ogdf;

#define NO_PARTITION (-1)

class QPartitioning;

template<class Value>
using PartitionArray = RegisteredArray<QPartitioning, int, Value>;

class QPartitioning : protected GraphObserver, public RegistryBase<int, QPartitioning> {
private:
	PartitionArray<List<node>> partitioned_nodes;
	NodeArray<int> partitions;
	int q_vertex_count = 0;
	int partition_next_id = 0;
	int partition_table_size = MIN_TABLE_SIZE;

public:
	explicit QPartitioning(const Graph* G) : partitions(*G, NO_PARTITION), GraphObserver(G) {
		partitioned_nodes.init(*this);
	}

	bool isQVertex(node n) const;

	int getPartitionOf(node n) const;

	// XXX it is assumed that the graph structure ensures that only two rotations are possible, i.e. by calling makeWheel
	int makeQVertex(node n, int p = NO_PARTITION);

	void releaseQVertex(node n);

	int partitionCount() const { return partition_next_id; }

	int qVertexCount() const { return q_vertex_count; }

	const List<node>& nodesInPartition(int partition) const { return partitioned_nodes[partition]; }

	bool isKeyAssociated(int key) const override { return 0 <= key && key <= maxKeyIndex(); }

	int keyToIndex(int key) const override { return key; }

	int keyArrayTableSize() const override { return partition_table_size; }

	int maxKeyIndex() const override { return partition_next_id - 1; }

protected:
	void nodeDeleted(node v) override;

	void nodeAdded(node v) override {};

	void edgeDeleted(edge e) override {};

	void edgeAdded(edge e) override {};

	void cleared() override {
		partitions.fill(NO_PARTITION);
		partition_next_id = 1;
		partition_table_size = MIN_TABLE_SIZE;
		reinitArrays();
	};
};
