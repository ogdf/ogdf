/** \file
 * \brief Implementation of EdgeComparer.
 *
 * \author Karsten Klein
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

#include <ogdf/basic/EdgeComparer.h>
#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/planarity/PlanRep.h>

namespace ogdf {

//compares outgoing adjEntries
int EdgeComparer::compare(const adjEntry& e1, const adjEntry& e2) const {
	//identify if the edges meet at a common point, otherwise
	//sort by index (stable)
	//are the adjentries sourceentries (edges outgoing) ?
	bool sAdj1 = (e1 == e1->theEdge()->adjSource());
	bool sAdj2 = (e2 == e2->theEdge()->adjSource());
	//generic nodes/edges
	node s1 = (m_PR ? m_PR->original(e1->theNode()) : e1->theNode());
	node s2 = (m_PR ? m_PR->original(e2->theNode()) : e2->theNode());
	edge ed1 = (m_PR ? m_PR->original(e1->theEdge()) : e1->theEdge());
	edge ed2 = (m_PR ? m_PR->original(e2->theEdge()) : e2->theEdge());
	node t1 = (m_PR ? m_PR->original(e1->twinNode()) : e1->twinNode());
	node t2 = (m_PR ? m_PR->original(e2->twinNode()) : e2->twinNode());

	//a: adjentry node(s)
	//b: twin nodes

	double x1a, x1b, x2a, x2b, y1a, y1b, y2a, y2b;

#if 0
	//debug stuff
	bool output = false;
	if ( ( (s1->degree() == 5) || (s2->degree() == 5)) &&
		( (m_AG->type(ed1)== Graph::generalization) ||
		(m_AG->type(ed2)== Graph::generalization)))
	{
		output = true;
	}
#endif

	x1a = m_AG->x(s1);
	x2a = m_AG->x(s2);
	y1a = m_AG->y(s1);
	y2a = m_AG->y(s2);


	//meet check not yet implemented, assume same start point
	OGDF_ASSERT(
			!((!OGDF_GEOM_ET.equal(x1a, x2a)) && (!(OGDF_GEOM_ET.equal(y1a, y2a)))) || (s1 != s2));

	//check if we have bends without representation node
	//use them as second end point
	const DPolyline& bends1 = m_AG->bends(ed1);
	const DPolyline& bends2 = m_AG->bends(ed2);

	ListConstIterator<DPoint> it1;
	ListConstIterator<DPoint> it2;

	bool hasBends1 = bends1.size() > 0;
	bool hasBends2 = bends2.size() > 0;

	if (hasBends1) {
		if (sAdj1) {
			it1 = bends1.begin();
		} else {
			it1 = bends1.rbegin();
		}

		x1b = (*it1).m_x;
		y1b = (*it1).m_y;
	} else {
		x1b = m_AG->x(t1);
		y1b = m_AG->y(t1);
	}
	if (hasBends2) {
		if (sAdj2) {
			it2 = bends2.begin();
		} else {
			it2 = bends2.rbegin();
		}

		x2b = (*it2).m_x;
		y2b = (*it2).m_y;
	} else {
		x2b = m_AG->x(t2);
		y2b = m_AG->y(t2);
	}
	//special condition: if we have parallel segments, e.g. merger edges,
	//we don't use the bend position
	if (hasBends1 && hasBends2) {
		bool end1 = false;
		//try to run over the bends / edge endpoint to differing points
		while ((x1b == x2b) && (y1b == y2b) && (x1a == x2a) && (y1a == y2a)) {
			//move one step further
			if (it1.valid()) {
				if (sAdj1) {
					++it1;
				} else {
					--it1;
				}
			}
			if (it2.valid()) {
				if (sAdj2) {
					++it2;
				} else {
					--it2;
				}
			}

			if (it1.valid()) {
				x1b = (*it1).m_x;
				y1b = (*it1).m_y;
			} else {
				x1b = m_AG->x(t1);
				y1b = m_AG->y(t1);
				end1 = true;
			}
			if (it2.valid()) {
				x2b = (*it2).m_x;
				y2b = (*it2).m_y;
			} else {
				x2b = m_AG->x(t2);
				y2b = m_AG->y(t2);
				//stop searching at both endpoints
				if (end1) {
					break;
				}
			}
		}
	}

	//now we have all the points necessary to sort
#if 0
	double dx1, dx2, dy1, dy2;
	dx1 = x1b - x1a;
	dx2 = x2b - x2a;
	dy1 = y1b - y1a;
	dy2 = y2b - y2a;
#endif

#if 0
	//debug
	std::ofstream f("c:\\Temp\\Karsten\\ASorting.txt", std::ios::app);

	f << "\nEntries at node: " << s1->index() <<"\n";
	f << "Compare " << s1->index() <<"->"<<t1->index() << " " <<m_AG->type(ed1)<<" , \n"
		<< s2->index() <<"->"<<t2->index() << " " <<m_AG->type(ed1)<< "\n";

	f << "Resultat: "<<-orientation(DPoint(x1a, y1a),
		DPoint(x1b, y1b),
		DPoint(x2b, y2b)) << "\n";
#endif

	//hier rueckgabe Vergleich e1->id, e2->id, falls nicht selber knoten
	//...
	//clockwise ist -, counterclockwise ist +
	//in AGD alle +!
#if 0
	return -compareVectors(dx1, dy1, dx2, dy2);
#endif
	//HIER fuer debuggen -davor

	if (s1 == s2) {
		node uvw = s1->firstAdj()->twinNode();
		edge ed3 = s1->firstAdj()->theEdge();
		const DPolyline& bends3 = m_AG->bends(ed3);
		DPoint dp;
		//as we use different comparison points for edges
		//at different compares, we have to assure that
		//we always have the same position of a edge's comparison
		//point compared to our reference vector
		//
		//
		//           xxxxxxx
		//           |  |  |l           l = Mess-Vektor
		//           x  |  x l
		//              |     l
		//              |-------o
		//--
		//workaround: wir kreisen, bis wir nicht auf einer Generalisierung
		//liegen. Das muss nicht immer moeglich sein und selbst dann koennen
		//Assoziationen uebereinander liegen.
		//Oder: Wir gehen immer bis zum zweiten Knick (immer diff bei merger)
		if (bends3.size() > 0) {
			if (bends3.size() > 1) {
				ListConstIterator<DPoint> itb;
				if (s1->firstAdj() == ed3->adjSource()) {
					itb = bends3.begin();
					++itb;
					dp = (*itb);
				} else {
					itb = bends3.rbegin();
					--itb;
					dp = (*itb);
				}
			} else {
				if (s1->firstAdj() == ed3->adjSource()) {
					dp = bends3.front();
				} else {
					dp = bends3.back();
				}
			}
		}
		//--
		else {
			dp = DPoint(m_AG->x(uvw) - 2, m_AG->y(uvw) + 1);
		}
		DPoint dp1a(x1a, y1a);
		double w1 = dp1a.angle(dp, DPoint(x1b, y1b));
		double w2 = dp1a.angle(dp, DPoint(x2b, y2b));
		OGDF_ASSERT(w1 >= 0);
		OGDF_ASSERT(w2 >= 0);

		//workaround shortcut (can be inserted above)
		if (ed1 == ed3) {
			return 1; //Reference edge is first edge
		}
		if (ed2 == ed3) {
			return -1;
		}

#if 0
		// debug stuff
		if (output)
		{
			std::ofstream fout("c:\\temp\\ogdl\\EComp.txt", std::ios::app);
			fout << "Knoten an Position "<<x1a<<"/"<<y1a<<"\n\n";
			fout << "Positionen: \n" << "DPReferenz: "<<dp.m_x<<"/"<<dp.m_y<<"\n";
			fout << "Punkt 1: "<<x1b<<"/"<<y1b<<" Typ "<<m_AG->type(ed1)<<
				" Winkel: "<<w1<< "\n";
			fout << "Punkt 2: "<<x2b<<"/"<<y2b<<" Typ "<<m_AG->type(ed2)<<
				" Winkel: "<<w2<<"\n";
			fout << "1 vor 2? " << (w1<w2? "JA" :(w1>w2?"NEIN":"Vielleicht")) <<"\n";
		}
#endif

		if (w1 < w2) {
			return 1;
		} else if (w1 > w2) {
			return -1;
		} else {
			return 0;
		}
	} else {
		return orientation(DPoint(x1a, y1a), DPoint(x1b, y1b), DPoint(x2b, y2b));
	}
}

bool EdgeComparer::before(const DPoint& u, const DPoint& v, const DPoint& w) const {
#if 0
	double dx1 = v.m_x - u.m_x;
	double dx2 = w.m_x - u.m_x;
	double dy1 = v.m_y - u.m_y;
	double dy2 = w.m_y - u.m_y;

	return (compareVectors(dx1, dy1, dx2, dy2) < 0);
#else
	return orientation(u, v, w) > 0;
#endif
}

//Is based on the sign of the determinant of a matrix
//defined by the point coordinates
//respects the flipping of y axis!!
//TODO: shift into geometric
int EdgeComparer::orientation(const DPoint& u, const DPoint& v, const DPoint& w) const {
	double plus1 = v.m_x * w.m_y;
	double plus2 = w.m_x * u.m_y;
	double plus3 = u.m_x * v.m_y;
	double minus1 = v.m_x * u.m_y;
	double minus2 = w.m_x * v.m_y;
	double minus3 = u.m_x * w.m_y;

	double E = plus1 + plus2 + plus3 - minus1 - minus2 - minus3;

	if (E > 0) {
		return 1;
	}
	if (E < 0) {
		return -1;
	}
	return 0;
}

}
