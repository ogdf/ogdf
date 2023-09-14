/** \file
 * \brief Computes a position that induces a minimal number of crossings for a given vertex and straight-line drawing.
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

#define OGDF_GEOMETRIC_INEXACT_NUMBER_TYPE true

#include <ogdf/geometric/CrossingMinimalPosition.h>

#ifdef OGDF_INCLUDE_CGAL

#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/graph/GeometricDrawing.h>
#	include <ogdf/geometric/cr_min/graph/OGDFGraphWrapper.h>
#	include <ogdf/geometric/cr_min/graph/geometric_crossing_min/CrossingMinimalPositionRnd.h>
#	include <ogdf/geometric/cr_min/tools/math.h>

#	include <limits>

#	include <CGAL/Gmpq.h>
#	include <CGAL/Simple_cartesian.h>

namespace ogdf {


template<typename FT>
DPoint CrossingMinimalPosition<FT>::call(GraphAttributes& GA, node v) {
	using namespace internal::gcm;
	using Kernel = CGAL::Simple_cartesian<FT>;

	using _Graph = graph::OGDFGraphWrapper;
	using Drawing = graph::GeometricDrawing<Kernel, _Graph>;

	_Graph _g(GA.constGraph());
	Drawing drawing(_g);

#	ifdef OGDF_GEOMETRIC_CR_MIN_DEBUG
	std::cout << "[MOVE VERTEX] " << v << "(" << drawing.get_point(v) << ", degree: " << v->degree()
			  << ")" << std::endl;
#	endif

	graph::ogdf_attributes_to_geometric_drawing(GA, drawing);

	rnd.seed(ogdf::randomSeed());

	std::vector<edge> samples;
	for (auto e : _g.edges()) {
		if (!e->isIncident(v)) {
			samples.push_back(e);
		}
	}

	geometry::Rectangle_t<Kernel> bbox(m_x_min, m_y_min, m_x_max, m_y_max);

	if (v->degree() > 1) {
		auto opt = graph::CrossingMinimalPositionRnd<Kernel, _Graph>::compute(drawing, v, samples,
				m_number_of_point_samples, m_neighborhood_threshold, bbox, m_within_region, rnd);

		drawing.set_point(v, opt);
	}

	std::uniform_real_distribution<double> dist(-1, 1);
	while (geometry::has_collinear_pair(drawing, v)) {
		double d_x = dist(rnd);
		double d_y = dist(rnd);
		auto p = drawing.get_point(v);
		drawing.set_point(v, {p.x() + d_x, p.y() + d_y});
	}

	DPoint p;
	p.m_x = CGAL::to_double(drawing.get_point(v).x());
	p.m_y = CGAL::to_double(drawing.get_point(v).y());
	return p;
}

template class CrossingMinimalPosition<double>;
template class CrossingMinimalPosition<CGAL::Gmpq>;

}
#else
namespace ogdf {

template<typename FT>
DPoint CrossingMinimalPosition<FT>::call(GraphAttributes& GA, node v) {
	OGDF_THROW_PARAM(LibraryNotSupportedException, LibraryNotSupportedCode::Cgal);
}

template class CrossingMinimalPosition<double>;

}


#endif
