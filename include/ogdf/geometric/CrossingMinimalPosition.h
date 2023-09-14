/** \file
 * \brief Computes a position that induces a minimal number of crossings for a given vertex and straight-line drawing.
 *
 * Based on the implementation and techniques of the following papers:
 *
 * Marcel Radermacher, Klara Reichard, Ignaz Rutter, and Dorothea Wagner.
 * Geometric Heuristics for Rectilinear Crossing Minimization.
 * Journal of Experimental Algorithmics 24:1, pages 1.12:1–1.12:21, 2019. doi:10.1145/3325861 .
 *
 * Marcel Radermacher and Ignaz Rutter.
 * Geometric Crossing Minimization - A Scalable Randomized Approach.
 * In Proceedings of the 27th Annual European Symposium on Algorithms (ESA’19).
 * Ed. by Michael A. Bender, Ola Svensson, and Grzegorz Herman. Leibniz International Proceedings in Informatics (LIPIcs), pages 76:1–76:16.
 * Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2019. doi: 10.4230/LIPIcs.ESA
 *
 * \pre Requires CGAL! See README.md in this folder.
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

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/geometric/VertexPositionModule.h>

#ifdef OGDF_INCLUDE_CGAL
#	include <CGAL/Gmpq.h>
#endif

namespace ogdf {

/**
 * \brief Compute a crossing minimal position for a vertex
 * @tparam FT internal floating point type. Selecting an exact floating point type ensures robust computations.
 *
 * \pre Requires CGAL! See README.md in this folder.
 */
template<typename FT>
class OGDF_EXPORT CrossingMinimalPosition : public VertexPositionModule {
public: // ~Initialize vertex position module
	CrossingMinimalPosition() { }

	~CrossingMinimalPosition() { }

	/**
	 */
	DPoint call(GraphAttributes& GA, node v);

	void setExactComputation() {
		m_number_of_edge_samples = -1;
		m_number_of_point_samples = 1;
		m_neighborhood_threshold = -1;
		m_within_region = true;
	}

	/**set the number of edges that are randomly selected to compute the new vertex postion and the number of points that are tested within the best region.
	* @param number_of_edge_samples  number of randomly selected edges
	* @param number_of_point_samples number of randomly selected point
	*/
	void setSampleSize(const unsigned int number_of_edge_samples,
			const unsigned int number_of_point_samples) {
		m_number_of_edge_samples = number_of_edge_samples;
		m_number_of_point_samples = number_of_point_samples;
	}

	/** If the value \p within_region is set to false, the algorithm samples points outside the optimal region.
	 * @param within_region indicates whether the new position has to be within the crossing minimal region
	 */
	void computePositionInOptimalRegion(const bool within_region) {
		m_within_region = within_region;
	}

	/** The algortihm randomly partitions the neighbor of the vertex into blocks of size \p threshold.
	 * A large vertex degree results in long running times and considerable memory consumption.
	 * @param threshold indicates the block size
	 */
	void setNeighboorhoodThreshold(const unsigned int threshold) {
		m_neighborhood_threshold = threshold;
	}

	unsigned int neighborThreshold() const { return m_neighborhood_threshold; }

protected:
	unsigned int m_number_of_edge_samples = -1;
	unsigned int m_number_of_point_samples = 1;
	unsigned int m_neighborhood_threshold = 100;
	bool m_within_region = true;

	std::mt19937_64 rnd;
};

using CrossingMinimalPositionFast = CrossingMinimalPosition<double>;

#ifdef OGDF_INCLUDE_CGAL
using CrossingMinimalPositionPrecise = CrossingMinimalPosition<CGAL::Gmpq>;
#endif


class OGDF_EXPORT CrossingMinimalPositionApx : public CrossingMinimalPositionFast {
public:
	CrossingMinimalPositionApx() {
		m_number_of_edge_samples = 512;
		m_number_of_point_samples = 1000;
		m_neighborhood_threshold = 100;
		m_within_region = true;
	}
};

class OGDF_EXPORT CrossingMinimalPositionApxWeighted : public CrossingMinimalPositionFast {
public:
	CrossingMinimalPositionApxWeighted() {
		m_number_of_edge_samples = 512;
		m_number_of_point_samples = 1000;
		m_neighborhood_threshold = 100;
		m_within_region = false;
	}
};


}
