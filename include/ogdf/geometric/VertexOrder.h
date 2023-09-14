/** \file
 * \brief Computes an vertex order based on the number of crossings in a given (straight-line) drawing.
 *
 * \author Marcel Radermacher
 *
 * \pre Requires CGAL! See README.md in this folder.
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

#include <ogdf/basic/GraphAttributes.h>

namespace ogdf {

enum class OrderEnum { asc, desc, rnd };
enum class MeasureEnum { zero, log, sum, squared };

/**
 * \pre Requires CGAL! See README.md in this folder.
 */
class OGDF_EXPORT CrossingVertexOrder {
private:
	GraphAttributes ga;

	using QElement = std::pair<ogdf::node, unsigned int>;

	std::vector<QElement> vertex_order;
	OrderEnum o;
	MeasureEnum m;

	void sort();

	double crossings(int c);

	void init_all();
	void init();


	void init_cr(edge e);

public:
	CrossingVertexOrder(GraphAttributes& _ga, OrderEnum _o, MeasureEnum _m)
		: ga(_ga), o(_o), m(_m) {
		//nothing to do
	}

	List<node> get_vertex_order() {
		init();
		List<node> order;
		for (auto v : vertex_order) {
			order.pushBack(v.first);
		}
		return order;
	}

	List<node> get_vertex_order_by_crossed_edges(edge e) {
		vertex_order.clear();
		init_cr(e);
		List<node> order;
		for (auto v : vertex_order) {
			order.pushBack(v.first);
		}
		return order;
	}
};
}
