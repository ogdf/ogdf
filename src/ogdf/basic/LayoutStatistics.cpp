/** \file
 * \brief Implements class LayoutStatistics which provides various
 *        functions for computing statistical measures of a layout.
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

#include <ogdf/basic/LayoutStatistics.h>
#include <ogdf/basic/Math.h>


namespace ogdf {

	double LayoutStatistics::edgeLengths(
		const GraphAttributes &ga,
		double *pMinLength,
		double *pMaxLength,
		double *pAvgLength,
		double *pStdDeviation,
		bool    considerSelfLoops)
	{
		const Graph &G = ga.constGraph();
		int m = G.numberOfEdges();

		double totalLength = 0, minLength = numeric_limits<double>::max(), maxLength = -numeric_limits<double>::max();

		EdgeArray<double> len(G);
		int nSelfLoops = 0;

		for(edge e : G.edges) {
			if(!considerSelfLoops && e->isSelfLoop()) {
				nSelfLoops++;
				continue;
			}

			const DPolyline &dpl = ga.bends(e);

			if(!dpl.empty()) {
				len[e] = dpl.length();

			} else {
				DPoint pv = DPoint(ga.x(e->source()),ga.y(e->source()));
				DPoint pw = DPoint(ga.x(e->target()),ga.y(e->target()));
				len[e] = pv.distance(pw);
			}

			totalLength += len[e];
			minLength = min(minLength, len[e]);
			maxLength = max(maxLength, len[e]);
		}

		m -= nSelfLoops;

		double avgEdgeLength = totalLength / m;
		if(pAvgLength) *pAvgLength = avgEdgeLength;
		if(pMinLength) *pMinLength = minLength;
		if(pMaxLength) *pMaxLength = maxLength;

		if(pStdDeviation) {
			double sum = 0;
			for(edge e : G.edges) {
				if(!considerSelfLoops && e->isSelfLoop())
					continue;
				double d = len[e] - avgEdgeLength;
				sum += d*d;
			}

			*pStdDeviation = sqrt(sum / m);
		}

		return totalLength;
	}


	int LayoutStatistics::numberOfBends(
		const GraphAttributes &ga,
		int *pMinBendsPerEdge,
		int *pMaxBendsPerEdge,
		double *pAvgBendsPerEdge,
		double *pStdDeviation,
		bool    considerSelfLoops)
	{
		const Graph &G = ga.constGraph();
		int m = G.numberOfEdges();

		int totalBends = 0, minBends = numeric_limits<int>::max(), maxBends = 0;

		EdgeArray<int> bends(G);
		int nSelfLoops = 0;

		for(edge e : G.edges) {
			if(!considerSelfLoops && e->isSelfLoop()) {
				nSelfLoops++;
				continue;
			}

			const DPolyline &dpl = ga.bends(e);

			bends[e] = max(0, dpl.size() - 2);

			totalBends += bends[e];
			minBends = min(minBends, bends[e]);
			maxBends = max(maxBends, bends[e]);
		}

		m -= nSelfLoops;

		double avgBends = double(totalBends) / m;
		if(pAvgBendsPerEdge) *pAvgBendsPerEdge = avgBends;
		if(pMinBendsPerEdge) *pMinBendsPerEdge = minBends;
		if(pMaxBendsPerEdge) *pMaxBendsPerEdge = maxBends;

		if(pStdDeviation) {
			double sum = 0;
			for(edge e : G.edges) {
				if(!considerSelfLoops && e->isSelfLoop())
					continue;
				double d = bends[e] - avgBends;
				sum += d*d;
			}

			*pStdDeviation = sqrt(sum / m);
		}

		return totalBends;
	}


	double LayoutStatistics::angularResolution(
		const GraphAttributes &ga,
		double *pMaxAngle,
		double *pAvgAngle,
		double *pStdDeviation,
		bool    considerBends)
	{
		const Graph &G = ga.constGraph();

		double minAngle = 2*Math::pi, maxAngle = 0, sumAngles = 0;

		int numAngles = 0;
		ListPure<double> allAngles;

		for (node v : G.nodes) {

			double vx = ga.x(v), vy = ga.y(v);

			List<double> angles;
			for (adjEntry adj : v->adjEdges) {
				const DPolyline &dpl = ga.bends(adj->theEdge());
				double ex, ey;
				if (dpl.empty()) {
					ex = ga.x(adj->twinNode());
					ey = ga.y(adj->twinNode());
				}
				else {
					ex = dpl.front().m_x;
					ey = dpl.front().m_y;
				}

				angles.pushBack(atan2(ex-vx, ey-vy));
			}

			if (angles.size() < 2)
				continue;

			numAngles += angles.size();
			angles.quicksort();

			double lastAngle = angles.back();
			for (double psi : angles) {
				double alpha = psi - lastAngle;
				if (pStdDeviation)
					allAngles.pushBack(alpha);

				sumAngles += alpha;
				minAngle = min(minAngle, alpha);
				maxAngle = max(maxAngle, alpha);

				lastAngle = psi;
			}
		}

		if (considerBends) {
			for (edge e : G.edges) {
				DPolyline dpl = ga.bends(e);

				dpl.pushFront( DPoint(ga.x(e->source()), ga.y(e->source())) );
				dpl.pushBack ( DPoint(ga.x(e->target()), ga.y(e->target())) );
				dpl.normalize();

				if (dpl.size() < 3)
					continue;

				for (ListConstIterator<DPoint> it = dpl.begin().succ(); it != dpl.rbegin(); ++it) {
					double bx = (*it).m_x, by = (*it).m_y;

					const DPoint &p1 = *it.pred();
					double psi1 = atan2(p1.m_x-bx, p1.m_y-by);

					const DPoint &p2 = *it.succ();
					double psi2 = atan2(p2.m_x - bx, p2.m_y - by);

					double alpha = fabs(psi1 - psi2);
					if (alpha > Math::pi)
						alpha -= Math::pi;

					sumAngles += 2 * Math::pi;
					minAngle = min(minAngle, alpha);
					maxAngle = max(maxAngle, alpha + Math::pi);

					if (pStdDeviation) {
						numAngles += 2;
						allAngles.pushBack(alpha);
						allAngles.pushBack(alpha*Math::pi);
					}
				}
			}
		}

		double avgAngle = sumAngles / numAngles;
		if (pAvgAngle) *pAvgAngle = avgAngle;
		if (pMaxAngle) *pMaxAngle = maxAngle;

		if (pStdDeviation) {
			double sum = 0;
			for (double alpha : allAngles) {
				double d = alpha - avgAngle;
				sum += d*d;
			}

			*pStdDeviation = sqrt(sum / numAngles);
		}

		return minAngle;
	}


	int LayoutStatistics::numberOfCrossings(const GraphAttributes &ga)
	{
		Graph H;
		NodeArray<DPoint> points;

		NodeArray<node> origNode;
		EdgeArray<edge> origEdge;
		intersectionGraph(ga, H, points, origNode, origEdge);

		int ncr = 0;
		for(node v : H.nodes) {
			node vOrig = origNode[v];
			int d = (vOrig != nullptr) ? vOrig->degree() : 0;
			int k = (v->degree() - d) / 2;
			if(k >= 2)
				ncr += (k * (k-1)) / 2;
		}

		return ncr;
	}


} // end namespace ogdf
