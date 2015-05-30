/** \file
 * \brief Implementations of force-models for Spring-Embedder algorithm
 *
 * \author Carsten Gutwenger
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

#include <ogdf/internal/energybased/SEGV_ForceModel.h>
#include <ogdf/basic/Math.h>


namespace ogdf {

	//-----------------------------------------------------
	// Fruchterman / Reingold
	//-----------------------------------------------------

	DPoint SpringEmbedderGridVariant::ForceModelFR::computeDisplacement(int j, double boxLength) const
	{
		const double cEps = eps();

		const NodeInfo &vj = m_vInfo[j];
		int grid_x = vj.m_gridX;
		int grid_y = vj.m_gridY;

		// repulsive forces on node j: F_rep(d) = iel^2 / d
		DPoint force(0,0);
		for(int gi = -1; gi <= 1; gi++) {
			for(int gj = -1; gj <= 1; gj++) {
				for(int u : m_gridCell(grid_x+gi,grid_y+gj)) {

					if(u == j) continue;
					DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
					double d = dist.norm();

					if(d < boxLength) {
						dist /= d * d + cEps;
						force += dist;
					}
				}
			}
		}

		force *= m_idealEdgeLength * m_idealEdgeLength;
		DPoint disp(force);

		// attractive forces on j: F_attr(d) = -d^2 / iel
		force = DPoint(0,0);
		for(int l = vj.m_adjBegin; l != vj.m_adjStop; ++l) {
			int u = m_adjLists[l];

			DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
			double d = dist.norm();

			dist *= d;
			force -= dist;
		}

		force /= m_idealEdgeLength;
		disp += force;

		return disp;
	}


	//-----------------------------------------------------
	// Fruchterman / Reingold
	//   with modified attractive forces
	//-----------------------------------------------------

	DPoint SpringEmbedderGridVariant::ForceModelFRModAttr::computeDisplacement(int j, double boxLength) const
	{
		const double cEps = eps();

		const NodeInfo &vj = m_vInfo[j];
		int grid_x = vj.m_gridX;
		int grid_y = vj.m_gridY;

		// repulsive forces on node j: F_rep(d) = iel^3 / d
		DPoint force(0,0);
		for(int gi = -1; gi <= 1; gi++) {
			for(int gj = -1; gj <= 1; gj++) {
				for(int u : m_gridCell(grid_x+gi,grid_y+gj)) {

					if(u == j) continue;
					DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
					double d = dist.norm();

					if(d < boxLength) {
						dist /= d * d + cEps;
						force += dist;
					}
				}
			}
		}

		force *= m_idealEdgeLength * m_idealEdgeLength * m_idealEdgeLength;
		DPoint disp(force);

		// attractive forces on j: F_attr(d) = -d^3 / iel
		force = DPoint(0,0);
		for(int l = vj.m_adjBegin; l != vj.m_adjStop; ++l) {
			int u = m_adjLists[l];

			DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
			double d = dist.norm();

			dist *= d*d;
			force -= dist;
		}

		force /= m_idealEdgeLength;
		disp += force;

		return disp;
	}


	//-----------------------------------------------------
	// Fruchterman / Reingold
	//   with modified repulsive repulsive
	//-----------------------------------------------------

	DPoint SpringEmbedderGridVariant::ForceModelFRModRep::computeDisplacement(int j, double boxLength) const
	{
		const double cEps = eps();

		const NodeInfo &vj = m_vInfo[j];
		int grid_x = vj.m_gridX;
		int grid_y = vj.m_gridY;

		// repulsive forces on node j: F_rep(d) = iel^2 / d^2
		DPoint force(0,0);
		for(int gi = -1; gi <= 1; gi++) {
			for(int gj = -1; gj <= 1; gj++) {
				for(int u : m_gridCell(grid_x+gi,grid_y+gj)) {

					if(u == j) continue;
					DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
					double d = dist.norm();

					if(d < boxLength) {
						dist /= d * d * d + cEps;
						force += dist;
					}
				}
			}
		}

		force *= m_idealEdgeLength * m_idealEdgeLength;
		DPoint disp(force);

		// attractive forces on j: F_attr(d) = -d^2 / iel^2
		force = DPoint(0,0);
		for(int l = vj.m_adjBegin; l != vj.m_adjStop; ++l) {
			int u = m_adjLists[l];

			DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
			double d = dist.norm();

			dist *= d;
			force -= dist;
		}

		force /= m_idealEdgeLength * m_idealEdgeLength;
		disp += force;

		return disp;
	}


	//-----------------------------------------------------
	// Eades
	//-----------------------------------------------------

	DPoint SpringEmbedderGridVariant::ForceModelEades::computeDisplacement(int j, double boxLength) const
	{
		const double cEps = eps();
		const double cIELEps = m_idealEdgeLength + cEps;

		const NodeInfo &vj = m_vInfo[j];
		int grid_x = vj.m_gridX;
		int grid_y = vj.m_gridY;

		// repulsive forces on node j: F_rep(d) = iel^2 / d
		DPoint force(0,0);
		for(int gi = -1; gi <= 1; gi++) {
			for(int gj = -1; gj <= 1; gj++) {
				for(int u : m_gridCell(grid_x+gi,grid_y+gj)) {

					if(u == j) continue;
					DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
					double d = dist.norm();

					if(d < boxLength) {
						dist /= d * d + cEps;
						force += dist;
					}
				}
			}
		}

		force *= m_idealEdgeLength * m_idealEdgeLength;
		DPoint disp(force);

		// attractive forces on j: F_attr(d) = -c iel log_2(d/iel)
		const double c = 0.1;

		DPoint forceRepSub(0,0); // subtract rep. force on adjacent vertices
		force = DPoint(0,0);
		for(int l = vj.m_adjBegin; l != vj.m_adjStop; ++l) {
			int u = m_adjLists[l];

			DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
			double d = dist.norm();

			force -= Math::log2((d+cEps) / cIELEps) * dist;
			if(d < boxLength) {
				double f = 1.0 / (d * d + cEps);
				forceRepSub += f * dist;
			}
		}

		force       *= c * m_idealEdgeLength;
		forceRepSub *= m_idealEdgeLength * m_idealEdgeLength;
		disp += force - forceRepSub;

		return disp;
	}


	//-----------------------------------------------------
	// Hachul (new method)
	//-----------------------------------------------------

	DPoint SpringEmbedderGridVariant::ForceModelHachul::computeDisplacement(int j, double boxLength) const
	{
		const double cEps = eps();
		const double cIELEps = m_idealEdgeLength + cEps;

		const NodeInfo &vj = m_vInfo[j];
		int grid_x = vj.m_gridX;
		int grid_y = vj.m_gridY;

		// repulsive forces on node j: F_rep(d) = iel^2 / d
		DPoint force(0,0);
		for(int gi = -1; gi <= 1; gi++) {
			for(int gj = -1; gj <= 1; gj++) {
				for(int u : m_gridCell(grid_x+gi,grid_y+gj)) {

					if(u == j) continue;
					DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
					double d = dist.norm();

					if(d < boxLength) {
						dist /= d * d + cEps;
						force += dist;
					}
				}
			}
		}

		force *= m_idealEdgeLength * m_idealEdgeLength;
		DPoint disp(force);

		// attractive forces on j: F_attr(d) = -d^2 log_2(d/iel) / iel
		DPoint forceRepSub(0,0); // subtract rep. force on adjacent vertices
		force = DPoint(0,0);
		for(int l = vj.m_adjBegin; l != vj.m_adjStop; ++l) {
			int u = m_adjLists[l];

			DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
			double d = dist.norm();

			force -= d * Math::log2((d+cEps) / cIELEps) * dist;
			if(d < boxLength) {
				double f = 1.0 / (d * d + cEps);
				forceRepSub += f * dist;
			}
		}

		force       /= m_idealEdgeLength;
		forceRepSub *= m_idealEdgeLength * m_idealEdgeLength;
		disp += force - forceRepSub;

		return disp;
	}


	//-----------------------------------------------------
	// Gronemann
	//-----------------------------------------------------

	DPoint SpringEmbedderGridVariant::ForceModelGronemann::computeDisplacement(int j, double boxLength) const
	{
		const double cEps = eps();
		const double cIELEps = m_idealEdgeLength + cEps;

		const NodeInfo &vj = m_vInfo[j];
		int grid_x = vj.m_gridX;
		int grid_y = vj.m_gridY;

		// repulsive forces on node j: F_rep(d) = iel^2 / d
		DPoint force(0,0);
		for(int gi = -1; gi <= 1; gi++) {
			for(int gj = -1; gj <= 1; gj++) {
				for(int u : m_gridCell(grid_x+gi,grid_y+gj)) {

					if(u == j) continue;
					DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
					double d = dist.norm();

					if(d < boxLength) {
						dist /= d * d + cEps;
						force += dist;
					}
				}
			}
		}

		force *= m_idealEdgeLength * m_idealEdgeLength;
		DPoint disp(force);

		// attractive forces on j: F_attr(d) = c / deg(v) * ln(d/iel)
		const double c = 0.5;

		DPoint forceRepSub(0,0); // subtract rep. force on adjacent vertices
		force = DPoint(0,0);
		for(int l = vj.m_adjBegin; l != vj.m_adjStop; ++l) {
			int u = m_adjLists[l];

			DPoint dist = vj.m_pos - m_vInfo[u].m_pos;
			double d = dist.norm();

			force -= log((d+cEps) / cIELEps) * dist;
			if(d < boxLength) {
				double f = 1.0 / (d * d + cEps);
				forceRepSub += f * dist;
			}
		}

		force       *= c * (vj.m_adjStop-vj.m_adjBegin);
		forceRepSub *= m_idealEdgeLength * m_idealEdgeLength;
		disp += force - forceRepSub;

		return disp;
	}


}
