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
#pragma once

#include <ogdf/basic/GraphObserver.h>
#include <ogdf/cluster/sync_plan/utils/RegisteredArray.h>

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
