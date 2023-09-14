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

#	include <ogdf/geometric/cr_min/geometry/objects/LineSegment.h>
#	include <ogdf/geometric/cr_min/geometry/objects/Point.h>

namespace ogdf {
namespace internal {
namespace gcm {

namespace tools {

template<typename stream = std::ostream>
class GeogebraExporter {
private:
	stream& outstream;
	unsigned int point_counter;

public:
	GeogebraExporter(stream& outstream_) : outstream(outstream_), point_counter(0) {
		//nothing to do
	}

	template<typename type>
	void add(const geometry::Point_t<type>& point) {
		outstream << "A" << point_counter << "=" << point << "\n";
		++point_counter;
	}

	template<typename Polygon>
	void add(const Polygon& polygon) {
		unsigned int old_counter = point_counter;
		for (unsigned int i = 0; i < polygon.size(); ++i) {
			add(polygon[i]);
		}

		outstream << "input = polygon[";

		for (unsigned int i = old_counter; i < point_counter; ++i) {
			outstream << "A" << i;
			if (i + 1 < point_counter) {
				outstream << ", ";
			}
		}

		outstream << "]\n";
	}

	template<typename type>
	void add(const geometry::LineSegment_t<type>& l) {
		unsigned int i = point_counter;
		add(l.source());
		add(l.target());
		outstream << "segment[A" << i << ", A" << i + 1 << "]\n";
	}
};
}
}
}
}

#endif
