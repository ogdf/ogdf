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

#ifdef OGDF_INCLUDE_CGAL

//CGAL example
#	include <ogdf/geometric/cr_min/geometry/objects/Polygon.h>

#	include <iostream>

#	include <CGAL/Constrained_Delaunay_triangulation_2.h>
#	include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#	include <CGAL/Point_2.h>
#	include <CGAL/Polygon_2.h>
#	include <CGAL/Triangulation_face_base_with_info_2.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace geometry {
/*! Wrapper for the CGAL restricted triangulation
 */
template<typename Kernel>
class RestrictedTriangulation {
private:
	struct FaceInfo2 {
		FaceInfo2() { }

		int nesting_level;

		bool in_domain() { return nesting_level % 2 == 1; }
	};

	using K = Kernel;
	using Vb = CGAL::Triangulation_vertex_base_2<K>;
	using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>;
	using Fb = CGAL::Constrained_triangulation_face_base_2<K, Fbb>;
	using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
	using Itag = CGAL::Exact_predicates_tag;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>;
	using Point = CGAL::Point_2<K>;
	using Polygon_2 = CGAL::Polygon_2<K>;

	void mark_domains(CDT& ct, typename CDT::Face_handle start, int index,
			std::list<typename CDT::Edge>& border) {
		if (start->info().nesting_level != -1) {
			return;
		}
		std::list<typename CDT::Face_handle> queue;
		queue.push_back(start);
		while (!queue.empty()) {
			typename CDT::Face_handle fh = queue.front();
			queue.pop_front();
			if (fh->info().nesting_level == -1) {
				fh->info().nesting_level = index;
				for (int i = 0; i < 3; i++) {
					typename CDT::Edge e(fh, i);
					typename CDT::Face_handle n = fh->neighbor(i);
					if (n->info().nesting_level == -1) {
						if (ct.is_constrained(e)) {
							border.push_back(e);
						} else {
							queue.push_back(n);
						}
					}
				}
			}
		}
	}

	//explore set of facets connected with non constrained edges,
	//and attribute to each such set a nesting level.
	//We start from facets incident to the infinite vertex, with a nesting
	//level of 0. Then we recursively consider the non-explored facets incident
	//to constrained edges bounding the former set and increase the nesting level by 1.
	//Facets in the domain are those with an odd nesting level.
	void mark_domains(CDT& cdt) {
		for (typename CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end();
				++it) {
			it->info().nesting_level = -1;
		}
		std::list<typename CDT::Edge> border;
		mark_domains(cdt, cdt.infinite_face(), 0, border);
		while (!border.empty()) {
			typename CDT::Edge e = border.front();
			border.pop_front();
			typename CDT::Face_handle n = e.first->neighbor(e.second);
			if (n->info().nesting_level == -1) {
				mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
			}
		}
	}

	void insert_polygon(CDT& cdt, const Polygon_2& polygon) {
		if (polygon.is_empty()) {
			return;
		}
		typename CDT::Vertex_handle v_prev = cdt.insert(*CGAL::cpp11::prev(polygon.vertices_end()));
		for (typename Polygon_2::Vertex_iterator vit = polygon.vertices_begin();
				vit != polygon.vertices_end(); ++vit) {
			typename CDT::Vertex_handle vh = cdt.insert(*vit);
			cdt.insert_constraint(vh, v_prev);
			v_prev = vh;
		}
	}

public:
	CDT run(const Polygon_t<Kernel>& polygon) {
		CDT cdt;
		insert_polygon(cdt, polygon);
		mark_domains(cdt);
		return std::move(cdt);
	}
};
}
}
}
}

#endif
