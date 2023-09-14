/** \file
 * \brief Optimal Vertex Position interface
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
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/geometry.h>

namespace ogdf {

/**
 * \brief Interface for computing a good / optimal vertex position
 */
class OGDF_EXPORT VertexPositionModule {
public:
	// ~Initialize vertex position module
	VertexPositionModule() { }

	~VertexPositionModule() { }

	/** Vertex has to be moved within the given bound
	 * \param x_min minimum x-coordinate
	 * \param y_min minimum y-coordinate
	 * \param x_max maximum x-coordinate
	 * \param y_max maximum y-coordinate
	 */
	void setBoundingBox(double x_min, double y_min, double x_max, double y_max) {
		m_x_min = x_min;
		m_y_min = y_min;
		m_x_max = x_max;
		m_y_max = y_max;
	}

	//! computes a good position for the vertex \p v with respect to \p GA
	virtual DPoint call(GraphAttributes& GA, node v) = 0;

	//! computes a good position for the vertex \p v with respect to \p GA
	DPoint operator()(GraphAttributes& GA, node v) { return call(GA, v); }

protected:
	double m_x_min = 0;
	double m_y_min = 0;
	double m_x_max = 1;
	double m_y_max = 1;
};

}
