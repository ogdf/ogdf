/** \file
 * \brief Manages the partitioning of Q-nodes in an instance of SyncPlan.
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

namespace ogdf::sync_plan {
class QPartitioning; // IWYU pragma: keep

#define OGDF_DECL_REG_ARRAY_TYPE(v, c) RegisteredArray<QPartitioning, v, c>
//! RegisteredArray for labeling the partitions in a QPartitioning with an arbitrary \p Value.
OGDF_DECL_REG_ARRAY(PartitionArray)
#undef OGDF_DECL_REG_ARRAY_TYPE

//! Manages the partitioning of Q-nodes in an instance of SyncPlan.
class OGDF_EXPORT QPartitioning : protected GraphObserver,
								  public RegistryBase<int, QPartitioning, int> {
private:
	PartitionArray<List<node>> partitioned_nodes;
	NodeArray<int> partitions;
	int q_vertex_count = 0;
	int partition_next_id = 0;

public:
	static inline int NO_PARTITION = -1;

	explicit QPartitioning(const Graph* G) : GraphObserver(), partitions(*G, NO_PARTITION) {
		partitioned_nodes.init(*this);
		reregister(G);
	}

	bool isQVertex(node n) const;

	int getPartitionOf(node n) const;

	//! Mark a node as Q-node.
	//! @note it is assumed that the graph structure ensures that only two rotations are possible, i.e. by calling SyncPlan::makeWheel()
	int makeQVertex(node n, int p = NO_PARTITION);

	void releaseQVertex(node n);

	int partitionCount() const { return partition_next_id; }

	int qVertexCount() const { return q_vertex_count; }

	const List<node>& nodesInPartition(int partition) const { return partitioned_nodes[partition]; }

	//! Returns the index of \p key.
	static inline int keyToIndex(int key) { return key; }

	//! Returns whether \p key is associated with this registry.
	bool isKeyAssociated(int key) const { return 0 <= key && key <= maxKeyIndex(); }

	//! Returns the maximum index of all keys managed by this registry.
	int maxKeyIndex() const { return partition_next_id - 1; }

	//! Returns the array size currently requested by this registry.
	int calculateArraySize(int add) const { return calculateTableSize(partition_next_id + add); }

	int begin() const { return 0; }

	int end() const { return partition_next_id; }

protected:
	void nodeDeleted(node v) override;

	void nodeAdded(node v) override {};

	void edgeDeleted(edge e) override {};

	void edgeAdded(edge e) override {};

	void cleared() override {
		partitioned_nodes.fillWithDefault();
		partitions.fill(NO_PARTITION);
		q_vertex_count = 0;
		partition_next_id = 0;
	};
};
}
