/*
 * $Revision: 3481 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-05-02 17:00:38 +0200 (Do, 02. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Implements GML write functionality of class GraphIO.
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

#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/Queue.h>


namespace ogdf {


//---------------------------------------------------------
// class GraphIO::SVGSettings
//---------------------------------------------------------

GraphIO::SVGSettings GraphIO::svgSettings;


GraphIO::SVGSettings::SVGSettings()
{
	m_margin = 1;

	m_fontSize = 10;
	m_fontColor = "#000000";
	m_fontFamily = "Arial";
}


//---------------------------------------------------------
// GraphIO::drawSVG
//---------------------------------------------------------

static void write_svg_header(ostream &os, double xmin, double ymin, double xmax, double ymax)
{
	os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
	os << "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" version=\"1.1\" baseProfile=\"full\" ";

	os << "width=\"" << (xmax - xmin) << "px\" ";
	os << "height=\"" << (ymax - ymin) << "px\" ";
	os << "viewBox=\"" << 0 << " " << 0 << " " << (xmax - xmin) << " " << (ymax - ymin) << "\">\n";
}

static void write_svg_footer(ostream &os)
{
	os << "</svg>\n";
}


static void write_dasharray(StrokeType lineStyle, double lineWidth, ostream &os)
{
	if(lineStyle == stNone || lineStyle == stSolid)
		return;

	os << "stroke-dasharray=\"";
	switch(lineStyle) {
	case stDash:
		os << 4*lineWidth << "," << 2*lineWidth;
		break;
	case stDot:
		os << 1*lineWidth << "," << 2*lineWidth;
		break;
	case stDashdot:
		os << 4*lineWidth << "," << 2*lineWidth << "," << 1*lineWidth << "," << 2*lineWidth;
		break;
	case stDashdotdot:
		os << 4*lineWidth << "," << 2*lineWidth << "," << 1*lineWidth << "," << 2*lineWidth << "," << 1*lineWidth << "," << 2*lineWidth;
		break;
	default:
		; // nothing to do
	}
	os << "\" ";
}


static void compute_bounding_box(const GraphAttributes &A, double &xmin, double &ymin, double &xmax, double &ymax)
{
	const Graph &G = A.constGraph();
	if(G.numberOfNodes() == 0) {
		xmin = xmax = ymin = ymax = 0;
		return;
	}

	node v = G.firstNode();
	xmin = xmax = A.x(v),
	ymin = ymax = A.y(v);

	forall_nodes(v, G) {
		double lw = (A.attributes() & GraphAttributes::nodeStyle) ? 0.5*A.strokeWidth(v) : 0.5;

		xmax = max(xmax, A.x(v) + A.width (v)/2 + lw);
		ymax = max(ymax, A.y(v) + A.height(v)/2 + lw);
		xmin = min(xmin, A.x(v) - A.width (v)/2 - lw);
		ymin = min(ymin, A.y(v) - A.height(v)/2 - lw);
	}

	if (A.attributes() & GraphAttributes::edgeGraphics) {
		edge e;
		forall_edges(e, G) {
			double lw = (A.attributes() & GraphAttributes::edgeStyle) ? 0.5*A.strokeWidth(e) : 0.5;

			const DPolyline &dpl = A.bends(e);
			if (!dpl.empty()) {
				for(ListConstIterator<DPoint> it = dpl.begin(); it.valid(); ++it) {
					xmax = max(xmax, (*it).m_x + lw);
					ymax = max(ymax, (*it).m_y + lw);
					xmin = min(xmin, (*it).m_x - lw);
					ymin = min(ymin, (*it).m_y - lw);
				}
			}
		}
	}
}

static void compute_bounding_box(const ClusterGraphAttributes &A, double &xmin, double &ymin, double &xmax, double &ymax)
{
	compute_bounding_box((const GraphAttributes &)A, xmin, ymin, xmax, ymax);

	const ClusterGraph &C = A.constClusterGraph();

	cluster c;
	forall_clusters(c,C) {
		if (c == C.rootCluster())
			continue;

		double lw = (A.attributes() & GraphAttributes::nodeStyle) ? 0.5*A.strokeWidth(c) : 0.5;

		xmax = max(xmax, A.x(c) + A.width (c) + lw);
		ymax = max(ymax, A.y(c) + A.height(c) + lw);
		xmin = min(xmin, A.x(c) - lw);
		ymin = min(ymin, A.y(c) - lw);
	}

}


// compute the list of clusters in level order, i.e., clusters with increasing depth
static void getLevelOrderClusters(const ClusterGraph &cg, SListPure<cluster> &clusters)
{
	Queue<cluster> Q;
	Q.append(cg.rootCluster());

	while(!Q.empty()) {
		cluster c = Q.pop();
		clusters.pushBack(c);

		ListConstIterator<cluster> it;
		for(it = c->cBegin(); it.valid(); ++it)
			Q.append(*it);
	}
}

// the arrow must not point to the center of the node but to its border
// XXX: we always assume rectangular nodes here
static void
fix_arrow_target_coordinates(DPoint const &source, DPoint &target, double w, double h)
{
	const double dx = target.m_x - source.m_x;
	const double dy = target.m_y - source.m_y;
	// check if coordinates are useful at all; if not: give up
	//if (w > dx*2 || h > dy*2) {
	//	return;
	//}
	double lx, ly;
	if (dx < 0) {
		lx = (dx + w) / dx;
	} else
	if (dx > 0) {
		lx = (dx - w) / dx;
	} else {
		lx = 0;
	}
	if (dy < 0) {
		ly = (dy + h) / dy;
	} else
	if (dy > 0) {
		ly = (dy - h) / dy;
	} else {
		ly = 0;
	}
	double l = (lx > ly ? lx : ly);
	// l should be in the interval [0, 1]
	target.m_x = source.m_x + l * dx;
	target.m_y = source.m_y + l * dy;
}

static void
arrow_head_coordinates(DPoint const &source, DPoint const &target, DPoint &arrow1, DPoint& arrow2, double length)
{
	const double w = 0.309016994374 * length; // arrowhead width/2
	// arrow head with 20 operations (1 sqrt, 8 mult/div, 11 add/sub):
	const double dx = target.m_x - source.m_x;
	const double dy = target.m_y - source.m_y;
	const double dist = sqrt(dx * dx + dy * dy);
	const double perc = 1 - length / dist;
	const double xC = source.m_x + perc * dx;
	const double yC = source.m_y + perc * dy;
	const double wd = w / dist;
	const double x = dy * wd;
	const double y = -dx * wd;
	arrow1.m_x = xC + x;
	arrow1.m_y = yC + y;
	arrow2.m_x = xC - x;
	arrow2.m_y = yC - y;
}

static void
draw_arrow_head(const GraphAttributes &A, double xmin, double ymin, ostream &os, edge e, bool reverse = false)
{
	StrokeType lineStyle = (A.attributes() & GraphAttributes::edgeStyle) ? A.strokeType(e) : stSolid;
	const DPolyline &dpl = A.bends(e);

	DPoint source, target, arrow1, arrow2;
	node u, v;
	if (reverse) {
		u = e->target();
		v = e->source();
	} else {
		u = e->source();
		v = e->target();
	}
	const double arrowLength = 0.15450849718747*(A.width(v) + A.height(v));
	bool fixme = true;
	if (dpl.empty()) { // single-line
		source.m_x = A.x(u);
		source.m_y = A.y(u);
		target.m_x = A.x(v);
		target.m_y = A.y(v);
	} else { // multi-line
		ListConstIterator<DPoint> it = (reverse ? dpl.begin() : dpl.rbegin());
		if ((*it).m_x < A.x(v) - A.width(v)/2
		 || (*it).m_x > A.x(v) + A.width(v)/2
		 || (*it).m_y < A.y(v) - A.height(v)/2
		 || (*it).m_y > A.y(v) + A.height(v)/2) {
			target.m_x = A.x(v);
			target.m_y = A.y(v);
		} else {
			target = (*it);
			fixme = false;
			if (reverse) {
				++it;
			} else {
				--it;
			}
		}
		if (it.valid()) {
			source = (*it);
		} else {
			source.m_x = A.x(u);
			source.m_y = A.y(u);
		}
	}
	if (fixme) {
		fix_arrow_target_coordinates(source, target, A.width(v)/2, A.height(v)/2);
	}
	arrow_head_coordinates(source, target, arrow1, arrow2, arrowLength);
	GraphIO::indent(os,1)
	  << "<polyline fill=\"none\" points=\""
	  << arrow1.m_x - xmin << "," << arrow1.m_y - ymin << " "
	  << target.m_x - xmin << "," << target.m_y - ymin << " "
	  << arrow2.m_x - xmin << "," << arrow2.m_y - ymin << "\" ";
	if (lineStyle != stNone) {
		if (A.attributes() & GraphAttributes::edgeStyle) {
			os << "stroke=\"" << A.strokeColor(e) << "\" ";
			write_dasharray(lineStyle, A.strokeWidth(e), os);
			os << "stroke-width=\"" << A.strokeWidth(e) << "px\" ";
		} else {
			os << "stroke=\"#000000\" ";
		}
	}
	os << " />\n";
}

static void write_svg_node_edges(
	const GraphAttributes &A,
	double xmin, double ymin,
	ostream &os,
	const GraphIO::SVGSettings &settings)
{
	const Graph &G = A.constGraph();

	edge e;
	forall_edges(e, G)
	{
		const DPolyline &dpl = A.bends(e);
		if (A.attributes() & GraphAttributes::edgeGraphics)
		{
			if(dpl.empty())
				GraphIO::indent(os,1) << "<line ";
			else
				GraphIO::indent(os,1) << "<polyline fill=\"none\" ";

			StrokeType lineStyle = (A.attributes() & GraphAttributes::edgeStyle) ? A.strokeType(e) : stSolid;

			if(lineStyle != stNone) {
				if (A.attributes() & GraphAttributes::edgeStyle) {
					os << "stroke=\"" << A.strokeColor(e) << "\" ";
					write_dasharray(lineStyle, A.strokeWidth(e), os);
					os << "stroke-width=\"" << A.strokeWidth(e) << "px\" ";
				} else {
					os << "stroke=\"#000000\" ";
				}
			}

			if (!dpl.empty()) { //polyline
				os << "points=\"";
				node v = e->source();
				// under certain circumstances start at the center of the source
				if (dpl.front().m_x < A.x(v) - A.width(v)/2
				 || dpl.front().m_x > A.x(v) + A.width(v)/2
				 || dpl.front().m_y < A.y(v) - A.height(v)/2
				 || dpl.front().m_y > A.y(v) + A.height(v)/2) {
					os << (A.x(e->source()) - xmin) << "," << (A.y(e->source()) - ymin) << " ";
				}

				// connect points
				for (ListConstIterator<DPoint> it = dpl.begin(); it.valid(); ++it) {
					os << ((*it).m_x - xmin) << "," << ((*it).m_y - ymin) << " ";
				}

				v = e->target();
				// under certain circumstances finish at the center of the target
				if (dpl.back().m_x < A.x(v) - A.width(v)/2
				 || dpl.back().m_x > A.x(v) + A.width(v)/2
				 || dpl.back().m_y < A.y(v) - A.height(v)/2
				 || dpl.back().m_y > A.y(v) + A.height(v)/2) {
					os << (A.x(e->target()) - xmin) << "," << (A.y(e->target()) - ymin) << " ";
				}
				os << "\" ";

			} else { // single line
				os << "x1=\"" << A.x(e->source()) - xmin << "\" ";
				os << "y1=\"" << A.y(e->source()) - ymin << "\" ";
				os << "x2=\"" << A.x(e->target()) - xmin << "\" ";
				os << "y2=\"" << A.y(e->target()) - ymin<< "\" ";
			}
			os << "/>\n";

			if (A.attributes() & GraphAttributes::edgeLabel
			 && !A.label(e).empty()) {
				double x, y;
				if (dpl.empty()) { // single-line
					x = (A.x(e->source()) + A.x(e->target())) * .5;
					y = (A.y(e->source()) + A.y(e->target())) * .5;
				} else { // poly-line
					ListConstIterator<DPoint> it = dpl.begin();
					const DPoint *lastPoint = &(*it);
					const DPoint *curPoint;
					double step = 0;
					double distance = 0;
					for (++it; it.valid(); ++it) {
						curPoint = &(*it);
						step += curPoint->distance(*lastPoint);
						lastPoint = curPoint;
					}
					step *= .5;

					it = dpl.begin();
					lastPoint = &(*it);
					for (++it; it.valid(); ++it) {
						curPoint = &(*it);
						distance = curPoint->distance(*lastPoint);
						if (distance <= step) {
							step -= distance;
						} else {
							break;
						}
						lastPoint = curPoint;
					}
					x = lastPoint->m_x;
					y = lastPoint->m_y;
					if (distance != 0) {
						x += (curPoint->m_x - lastPoint->m_x) * step / distance;
						y += (curPoint->m_y - lastPoint->m_y) * step / distance;
					}
				}
				GraphIO::indent(os, 1)
				  << "<text x=\"" << x - xmin
				  << "\" y=\"" << y - ymin
				  << "\" text-anchor=\"middle\" dominant-baseline=\"middle"
				  << "\" font-family=\"" << settings.fontFamily()
				  << "\" font-size=\"" << settings.fontSize()
				  << "\" fill=\"" << settings.fontColor() << "\">"
				  << A.label(e) << "</text>\n";
			}

			// draw arrows if G is directed or if arrow types are defined for the edge
			bool drawFirst = false;
			bool drawLast = false;
			if (A.attributes() & GraphAttributes::edgeArrow) {
				switch (A.arrowType(e)) {
				case eaUndefined:
					if (A.directed()) {
						drawLast = true;
					}
					break;
				case eaFirst:
					drawFirst = true;
					break;
				case eaBoth:
					drawFirst = true;
				case eaLast:
					drawLast = true;
				case eaNone:
					break;
				}
			} else
			if (A.directed()) {
				drawLast = true;
			}
			if (drawFirst) {
				draw_arrow_head(A, xmin, ymin, os, e, true);
			}
			if (drawLast) {
				draw_arrow_head(A, xmin, ymin, os, e);
			}
		}
	}

	node v;
	forall_nodes(v,G) {
		if (A.attributes() & GraphAttributes::nodeGraphics) {
			const double
			  x = A.x(v) - xmin,
			  y = A.y(v) - ymin,
			  hw1 = 0.5 * A.width(v),
			  hh1 = 0.5 * A.height(v),
			  qw1 = 0.5 * hw1,
			  qh1 = 0.43301270189222 * A.height(v),
			  pw1 = 0.475528258147577 * A.width(v),
			  ph1 = 0.154508497187474 * A.height(v),
			  pw2 = 0.293892626146236 * A.width(v),
			  ph2 = 0.404508497187474 * A.height(v),
			  ow1 = 0.461939766255643 * A.width(v),
			  oh1 = 0.191341716182545 * A.height(v),
			  ow2 = 0.191341716182545 * A.width(v),
			  oh2 = 0.461939766255643 * A.height(v);
			// values are precomputed to save expensive sin/cos calls
			switch (A.shape(v))
			{
			case shEllipse:
				GraphIO::indent(os,1) << "<ellipse ";
				os << "cx=\"" << x << "\" ";
				os << "cy=\"" << y << "\" ";
				os << "rx=\"" << hw1 << "\" ";
				os << "ry=\"" << hh1 << "\" ";
				break;
			case shTriangle:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x << ","
				  << y - hh1 << " "
				  << x - hw1 << ","
				  << y + hh1 << " "
				  << x + hw1 << ","
				  << y + hh1 << "\" ";
				break;
			case shInvTriangle:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x << ","
				  << y + hh1 << " "
				  << x - hw1 << ","
				  << y - hh1 << " "
				  << x + hw1 << ","
				  << y - hh1 << "\" ";
				break;
			case shPentagon:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x << "," << y - hh1 << " "
				  << x + pw1 << "," << y - ph1 << " "
				  << x + pw2 << "," << y + ph2 << " "
				  << x - pw2 << "," << y + ph2 << " "
				  << x - pw1 << "," << y - ph1 << "\" ";
				break;
			case shHexagon:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x + qw1 << "," << y + qh1 << " "
				  << x - qw1 << "," << y + qh1 << " "
				  << x - hw1 << "," << y << " "
				  << x - qw1 << "," << y - qh1 << " "
				  << x + qw1 << "," << y - qh1 << " "
				  << x + hw1 << "," << y << "\" ";
				break;
			case shOctagon:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x + ow1 << "," << y + oh1 << " "
				  << x + ow2 << "," << y + oh2 << " "
				  << x - ow2 << "," << y + oh2 << " "
				  << x - ow1 << "," << y + oh1 << " "
				  << x - ow1 << "," << y - oh1 << " "
				  << x - ow2 << "," << y - oh2 << " "
				  << x + ow2 << "," << y - oh2 << " "
				  << x + ow1 << "," << y - oh1 << "\" ";
				break;
			case shRhomb:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x + hw1 << "," << y << " "
				  << x << "," << y + hh1 << " "
				  << x - hw1 << "," << y << " "
				  << x << "," << y - hh1 << "\" ";
				break;
			case shTrapeze:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x - hw1 << "," << y + hh1 << " "
				  << x + hw1 << "," << y + hh1 << " "
				  << x + qw1 << "," << y - hh1 << " "
				  << x - qw1 << "," << y - hh1 << "\" ";
				break;
			case shInvTrapeze:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x - hw1 << "," << y - hh1 << " "
				  << x + hw1 << "," << y - hh1 << " "
				  << x + qw1 << "," << y + hh1 << " "
				  << x - qw1 << "," << y + hh1 << "\" ";
				break;
			case shParallelogram:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x - hw1 << "," << y + hh1 << " "
				  << x + qw1 << "," << y + hh1 << " "
				  << x + hw1 << "," << y - hh1 << " "
				  << x - qw1 << "," << y - hh1 << "\" ";
				break;
			case shInvParallelogram:
				GraphIO::indent(os,1)
				  << "<polygon points=\""
				  << x - hw1 << "," << y - hh1 << " "
				  << x + qw1 << "," << y - hh1 << " "
				  << x + hw1 << "," << y + hh1 << " "
				  << x - qw1 << "," << y + hh1 << "\" ";
				break;
			case shRect:
			case shRoundedRect:
			default: // unsupported: shImage
				GraphIO::indent(os,1) << "<rect ";
				os << "x=\"" << x - hw1 << "\" ";
				os << "y=\"" << y - hh1 << "\" ";
				if (A.shape(v) == shRoundedRect) {
					os << "rx=\"" << A.width(v) / 10 << "\" ";
					os << "ry=\"" << A.height(v) / 10 << "\" ";
				}
				os << "width=\"" << A.width(v) << "\" ";
				os << "height=\"" << A.height(v) << "\" ";
				break;
			}

			StrokeType lineStyle = (A.attributes() & GraphAttributes::nodeStyle) ? A.strokeType(v) : stSolid;

			if (A.attributes() & GraphAttributes::nodeStyle)
			{
				os << "fill=\"" << A.fillColor(v) << "\" ";
				if(lineStyle == stNone)
					os << "stroke=\"none\" ";
				else
					os << "stroke=\"" << A.strokeColor(v) << "\" ";
				write_dasharray(lineStyle, A.strokeWidth(v), os);
				os << "stroke-width=\"" << A.strokeWidth(v) << "px\" ";
			}

			os << "/>\n";

			if(A.attributes() & GraphAttributes::nodeLabel){
				GraphIO::indent(os,1) << "<text x=\"" << A.x(v) - xmin << "\" y=\"" << A.y(v) - ymin
					<< "\" text-anchor=\"middle\" dominant-baseline=\"middle"
					<< "\" font-family=\"" << settings.fontFamily()
					<< "\" font-size=\"" << settings.fontSize()
					<< "\" fill=\"" << settings.fontColor() << "\">"
					<< A.label(v) << "</text>\n";
			}
		}
	}
}


static void write_svg_clusters(
	const ClusterGraphAttributes &A,
	double xmin, double ymin,
	ostream &os,
	const GraphIO::SVGSettings &settings)
{
	const ClusterGraph &C = A.constClusterGraph();

	SListPure<cluster> clusters;
	getLevelOrderClusters(C, clusters);

	SListConstIterator<cluster> itC;
	for(itC = clusters.begin(); itC.valid(); ++itC)
	{
		cluster c = *itC;
		if(c == C.rootCluster())
			continue;

		double x = A.x(c);
		double y = A.y(c);
		double w = A.width(c);
		double h = A.height(c);

		GraphIO::indent(os,1) << "<rect ";
		os << "x=\"" << x - xmin << "\" ";
		os << "y=\"" << y - ymin << "\" ";
		os << "width=\"" << w << "\" ";
		os << "height=\"" << h << "\" ";

		if(A.fillPattern(c) == fpNone)
			os << "fill=\"none\" ";
		else
			os << "fill=\"" << A.fillColor(c) << "\" ";

		StrokeType lineStyle = A.strokeType(c);

		if(lineStyle == stNone)
			os << "stroke=\"none\" ";
		else
			os << "stroke=\"" << A.strokeColor(c) << "\" ";

		write_dasharray(lineStyle, A.strokeWidth(c), os);
		os << "stroke-width=\"" << A.strokeWidth(c) << "px\" ";
		os << "/>\n";
	}
}


bool GraphIO::drawSVG(const GraphAttributes &A, ostream &os, const SVGSettings &settings)
{
	double xmin, ymin, xmax, ymax;
	compute_bounding_box(A, xmin, ymin, xmax, ymax);

	double m = settings.margin();
	xmin -= m;
	ymin -= m;
	write_svg_header(os, xmin, ymin, xmax+m, ymax+m);

	write_svg_node_edges(A, xmin, ymin, os, settings);
	write_svg_footer(os);

	return true;
}

bool GraphIO::drawSVG(const ClusterGraphAttributes &A, ostream &os, const SVGSettings &settings)
{
	double xmin, ymin, xmax, ymax;
	compute_bounding_box(A, xmin, ymin, xmax, ymax);

	double m = settings.margin();
	xmin -= m;
	ymin -= m;
	write_svg_header(os, xmin, ymin, xmax+m, ymax+m);

	write_svg_clusters(A, xmin, ymin, os, settings);
	write_svg_node_edges(A, xmin, ymin, os, settings);

	write_svg_footer(os);

	return true;
}

}
