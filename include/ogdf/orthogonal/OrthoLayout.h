/*
 * $Revision: 3188 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-01-10 09:53:32 +0100 (Do, 10. Jan 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class OrthoLayout which represents an
 *        orthogonal planar drawing algorithm.
 *
 * \author Carsten Gutwenger, Sebastian Leipert, Karsten Klein
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/


#ifdef _MSC_VER
#pragma once
#endif


#ifndef OGDF_ORTHO_LAYOUT_H
#define OGDF_ORTHO_LAYOUT_H


#include <ogdf/module/LayoutPlanRepModule.h>
#include <ogdf/orthogonal/OrthoRep.h>


namespace ogdf {


//! The Orthogonal layout algorithm for planar graphs.
class OGDF_EXPORT OrthoLayout : public LayoutPlanRepModule
{
public:
	//! Creates an instance of Orthogonal layout and sets options to default values.
	OrthoLayout();


	// calls planar UML layout algorithm. Input is a planarized representation
	// PG of a connected component of the graph, output is a layout of the
	// (modified) planarized representation in drawing
	//! Calls the layout algorithm for planarized representation \a PG.
	/**
	 * \pre \a PG is embedded and \a adjExternal is an adjecenty entry in \a PG.
	 *
	 * @param PG          is the planarized representation for which a layout shall be computed.
	 *                    \a PG may be modified during the call.
	 * @param adjExternal is an adjaceny entry in \a PG that shall be on the external
	 *                    face of the drawing.
	 * @param drawing     is assigned the final layout.
	 */
	void call(PlanRep &PG, adjEntry adjExternal, Layout &drawing);

	/** @}
	 *  @name Optional parameters
	 *  @{
	 */

	//! Returns the minimum distance between edges and vertices.
	double separation() const {
		return m_separation;
	}

	//! Sets the minimum distance between vertices.
	void separation(double sep) {
		m_separation = sep;
	}

	//! Returns the option \a cOverhang, which specifies the minimal distance of incident edges to the corner of a vertex.
	/**
	 * \a cOverhang * \a separation is the minimum distance between the glue point of an edge and a corner of the vertex boundary.
	 */
	double cOverhang() const {
		return m_cOverhang;
	}

	//! Sets the option \a cOverhang, which specifies the minimal distance of incident edges to the corner of a vertex.
	void cOverhang(double c) {
		m_cOverhang = c;
	}

	//! Returns the desired margin around the drawing.
	/**
	 * This is the distance between the tight bounding box and the boundary of the drawing.
	 */
	double margin() const {
		return m_margin;
	}

	//! Sets the desired margin around the drawing.
	void margin(double m) {
		m_margin = m;
	}

	//! Returns whether the currently selected orthogonaliaztion model is \e progressive.
	bool progressive() const { return m_progressive; }

	//! Selects if the progressive (true) or traditional (false) orthogonalization model is used.
	/**
	 * In the progressive model, 180 degree angles at degree-2 nodes are preferred.
	 */
	void progressive(bool b) { m_progressive = b; }

	//! Returns whether scaling is used in the compaction phase.
	bool scaling() const { return m_useScalingCompaction; }

	//! Selects if scaling is used in the compaction phase.
	void scaling(bool b) { m_useScalingCompaction = b; }

	//! Set bound on the number of bends
	void bendBound(int i) {
		if(i >= 0) m_bendBound = i;
	}

	//! @}

private:
	// compute bounding box and move final drawing such that it is 0 aligned
	// respecting margins
	void computeBoundingBox(const PlanRep &PG, Layout &drawing);

	// options

	double m_separation;	//!< minimum distance between obkects
	double m_cOverhang;		//!< distance to corner (relative to node size)
	double m_margin;		//!< margin around drawing

	bool m_progressive; //!< use progressive ortho style (prefer 180 degree angles on deg-2 vertices).
	int m_bendBound;	//!< bounds the number of bends per edge in ortho shaper

	bool m_useScalingCompaction;	//!< use scaling for compaction
	int m_scalingSteps;				//!< number of scaling steps (NOT REALLY USED!)
};


} // end namespace ogdf


#endif
