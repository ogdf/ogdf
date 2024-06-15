#include <ogdf/basic/Math.h>
#include <ogdf/cluster/sync_plan/QPartitioning.h>

using namespace ogdf;

bool QPartitioning::isQVertex(node n) const { return partitions[n] != NO_PARTITION; }

int QPartitioning::getPartitionOf(node n) const {
	OGDF_ASSERT(isQVertex(n));
	return partitions[n];
}

int QPartitioning::makeQVertex(node n, int p) {
	OGDF_ASSERT(!isQVertex(n));
	OGDF_ASSERT(n->degree() > 2);
	// this could actually be > 3, but then we have to manually preserve the edge bijection / order
	// OGDF_ASSERT(!pq.components.isCutVertex(n));
	if (p == NO_PARTITION) {
		p = partition_next_id;
		partition_next_id++;
	}
	int new_size = Math::nextPower2(MIN_TABLE_SIZE, partition_next_id);
	if (new_size != partition_table_size) {
		partition_table_size = new_size;
		enlargeArrayTables();
	}
	partitioned_nodes[p].emplaceBack(n);
	partitions[n] = p;
	q_vertex_count++;
	return p;
}

void QPartitioning::releaseQVertex(node n) {
	OGDF_ASSERT(isQVertex(n));
	int p = partitions[n];
	partitions[n] = NO_PARTITION;
	// FIXME the list is very short in our case, but do we really need to do this in our undo operations?
	bool has_removed = partitioned_nodes[p].removeFirst(n);
	OGDF_ASSERT(has_removed);
	q_vertex_count--;
}

void QPartitioning::nodeDeleted(node v) {
	if (isQVertex(v)) {
		partitioned_nodes[getPartitionOf(v)].removeFirst(v);
		// partition might have turned empty, but we don't care
	}
}
