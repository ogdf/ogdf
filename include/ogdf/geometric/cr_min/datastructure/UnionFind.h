/** \file
 *
 * \author Marcel Radermacher
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

#include <ogdf/basic/basic.h>

#include <vector>

namespace ogdf {
namespace internal {
namespace gcm {
namespace datastructure {

/**
 * Implements the Union Find Datastructure to maintain disjoint sets efficiently.
 */
class UnionFind {
private:
	std::vector<int> data;

public:
	/**
	 * Create a new set representation with not more than @p max_element elements.
	 * Initialy every element is in its own set.
	 * @param max_element maximum number of elements
	 */
	UnionFind(unsigned int max_element) : data(max_element) { all_to_singletons(); }

	/**
	 * Assigns every element to a singleton set.
	 * Set id is equal to element id.
	 */


	void all_to_singletons() { data.assign(data.size(), -1); }

	/**
	 * Find the representative to element @p u
	 * @param u element
	 * @return representative of set containing @p u
	 */
	unsigned int find(unsigned int u) {
		OGDF_ASSERT(u < data.size());
		if (data[u] >= 0) {
			data[u] = find(data[u]);
			return data[u];
		} else {
			return u;
		}
	}

	unsigned int operator[](unsigned int u) { return find(u); }

	/**
	 *  Merge the two sets containing @p u and @p v
	 *  @param u element u
	 *  @param v element v
	 */
	void merge(unsigned int u, unsigned int v) {
		unsigned int set_u = find(u);
		unsigned int set_v = find(v);
		if (set_u == set_v) {
			return;
		}

		if (data[set_u] > data[set_v]) {
			data[set_u] = set_v;
		} else if (data[set_v] > data[set_u]) {
			data[set_v] = set_u;
		} else {
			data[set_u] = set_v;
			--data[set_v];
		}
	}
};
}
}
}
}
