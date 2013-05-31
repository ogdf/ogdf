/*
 * $Revision: 3521 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-05-31 14:52:33 +0200 (Fr, 31. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Implement class ClusterGraphAttributes
 *
 * \author Karsten Klein, Carsten Gutwenger
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


#include <ogdf/cluster/ClusterGraphAttributes.h>


namespace ogdf {


ClusterGraphAttributes::ClusterGraphAttributes(ClusterGraph& cg, long initAttributes)
	: GraphAttributes(cg.constGraph(), initAttributes | edgeType | nodeType | nodeGraphics | edgeGraphics),
	m_pClusterGraph(&cg), m_clusterInfo(cg), m_clusterTemplate(cg)
{
	//should we always fill the cluster infos here?
}//constructor


//reinitialize graph
void ClusterGraphAttributes::init(ClusterGraph &cg, long initAttributes)
{
	GraphAttributes::init(cg.constGraph(), initAttributes);

	m_pClusterGraph = &cg;
	m_clusterInfo.init(cg);
	m_clusterTemplate.init(cg);
}


//
// calculates the bounding box of the graph including clusters
const DRect ClusterGraphAttributes::boundingBox() const
{
	DRect bb = GraphAttributes::boundingBox();
	double minx = bb.p1().m_x;
	double miny = bb.p1().m_y;
	double maxx = bb.p2().m_x;
	double maxy = bb.p2().m_y;

	cluster c;
	forall_clusters(c,*m_pClusterGraph)
	{
		if(c == m_pClusterGraph->rootCluster())
			continue;

		double x1 = x(c);
		double y1 = y(c);
		double x2 = x1 + width(c);
		double y2 = y1 + height(c);

		if (x1 < minx) minx = x1;
		if (x2 > maxx) maxx = x2;
		if (y1 < miny) miny = y1;
		if (y2 > maxy) maxy = y2;
	}

	return DRect(minx, miny, maxx, maxy);
}


void ClusterGraphAttributes::updateClusterPositions(double boundaryDist)
{
	cluster c;
	//run through children and nodes and update size accordingly
	//we use width, height temporarily to store max values
	forall_postOrderClusters(c,*m_pClusterGraph)
	{
		ListIterator<node> nit = c->nBegin();
		ListConstIterator<ClusterElement*> cit = c->cBegin();
		//Initialize with first element
		if (nit.valid())
		{
			x(c) = m_x[*nit] - m_width[*nit]/2;
			y(c) = m_y[*nit] - m_height[*nit]/2;
			width(c) = m_x[*nit] + m_width[*nit]/2;
			height(c) = m_y[*nit] + m_height[*nit]/2;
			nit++;
		}
		else
		{
			if (cit.valid())
			{
				x(c) = x(*cit);
				y(c) = y(*cit);
				width(c) = x(*cit) + width(*cit);
				height(c) = y(*cit) + height(*cit);
				cit++;
			}
			else
			{
				x(c) = 0.0;
				y(c) = 0.0;
				width(c) = 1.0;
				height(c) = 1.0;
			}
		}
		//run through elements and update
		while (nit.valid())
		{
			if (x(c) > m_x[*nit] - m_width[*nit]/2)
				x(c) = m_x[*nit] - m_width[*nit]/2;
			if (y(c) > m_y[*nit] - m_height[*nit]/2)
				y(c) = m_y[*nit] - m_height[*nit]/2;
			if (width(c) < m_x[*nit] + m_width[*nit]/2)
				width(c) = m_x[*nit] + m_width[*nit]/2;
			if (height(c) < m_y[*nit] + m_height[*nit]/2)
				height(c) = m_y[*nit] + m_height[*nit]/2;
			nit++;
		}
		while (cit.valid())
		{
			if (x(c) > x(*cit))
				x(c) = x(*cit);
			if (y(c) > y(*cit))
				y(c) = y(*cit);
			if (width(c) < x(*cit) + width(*cit))
				width(c) = x(*cit) + width(*cit);
			if (height(c) < y(*cit) + height(*cit))
				height(c) = y(*cit) + height(*cit);
			cit++;
		}
		x(c) -= boundaryDist;
		y(c) -= boundaryDist;
		width(c) = width(c) - x(c) + boundaryDist;
		height(c) = height(c) - y(c) + boundaryDist;
	}
}



ostream &operator<<(ostream &os, ogdf::cluster c)
{
	if (c) os << c->index(); else os << "nil";
	return os;
}


} // end namespace ogdf

