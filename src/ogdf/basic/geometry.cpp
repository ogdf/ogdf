/** \file
 * \brief Implementation of Geometry classes like ogdf::DPoint, ogdf::DRect,
 * ogdf::DIntersectableRect, ogdf::DPolygon.
 *
 * \author Joachim Kupke
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


#include <ogdf/basic/EpsilonTest.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/basic/graphics.h>

#include <algorithm>
#include <cmath>
#include <ostream>
#include <utility>

namespace ogdf {

const EpsilonTest OGDF_GEOM_ET(1.0e-6);

// output the rect
std::ostream& operator<<(std::ostream& os, const DRect& dr) {
	os << "\nLower left corner: " << dr.m_p1;
	os << "\nUpper right corner: " << dr.m_p2;
	os << "\nWidth: " << dr.width();
	os << "\nHeight: " << dr.height();
	return os;
}

double DRect::parallelDist(const DSegment& d1, const DSegment& d2) const {
	OGDF_ASSERT((d1.isHorizontal() && d2.isHorizontal()) || (d1.isVertical() && d2.isVertical()));
	double d1min, d1max, d2min, d2max, paraDist, dist;
	if (d1.isVertical()) {
		d1min = d1.start().m_y;
		d1max = d1.end().m_y;
		d2min = d2.start().m_y;
		d2max = d2.end().m_y;
		paraDist = fabs(d1.start().m_x - d2.start().m_x);
	} else {
		d1min = d1.start().m_x;
		d1max = d1.end().m_x;
		d2min = d2.start().m_x;
		d2max = d2.end().m_x;
		paraDist = fabs(d1.start().m_y - d2.start().m_y);
	}
	if (d1min > d1max) {
		std::swap(d1min, d1max);
	}
	if (d2min > d2max) {
		std::swap(d2min, d2max);
	}
	if (d1min > d2max || d2min > d1max) { // no overlap
		dist = pointDist(d1.start(), d2.start());
		Math::updateMin(dist, pointDist(d1.start(), d2.end()));
		Math::updateMin(dist, pointDist(d1.end(), d2.start()));
		Math::updateMin(dist, pointDist(d1.end(), d2.end()));
	} else {
		dist = paraDist; // segments overlap
	}
	return dist;
}

// output the rect
std::ostream& operator<<(std::ostream& os, const DIntersectableRect& dr) {
	os << static_cast<DRect>(dr);
	os << "\nCenter: " << dr.center();
	os << "\nArea: " << dr.area();
	return os;
}

void DIntersectableRect::initAreaAndCenter() {
	m_area = (m_p2.m_x - m_p1.m_x) * (m_p2.m_y - m_p1.m_y);
	m_center.m_x = m_p1.m_x + 0.5 * (m_p2.m_x - m_p1.m_x);
	m_center.m_y = m_p1.m_y + 0.5 * (m_p2.m_y - m_p1.m_y);
}

void DIntersectableRect::move(const DPoint& point) {
	double dX = point.m_x - m_center.m_x;
	double dY = point.m_y - m_center.m_y;
	m_center = point;
	m_p1.m_x += dX;
	m_p1.m_y += dY;
	m_p2.m_x += dX;
	m_p2.m_y += dY;
}

double DIntersectableRect::distance(const DIntersectableRect& other) const {
	double dist = 0.0;
	if (!intersects(other)) {
		dist = parallelDist(top(), other.bottom());
		Math::updateMin(dist, parallelDist(left(), other.right()));
		Math::updateMin(dist, parallelDist(right(), other.left()));
		Math::updateMin(dist, parallelDist(bottom(), other.top()));
	}
	return dist;
}

bool DIntersectableRect::intersects(const DIntersectableRect& rectangle) const {
	bool intersect = false;
	if (contains(rectangle.m_center) || rectangle.contains(m_center)) {
		intersect = true;
	} else {
		DPoint p1(rectangle.m_p1.m_x, rectangle.m_p2.m_y);
		DPoint p2(rectangle.m_p2.m_x, rectangle.m_p1.m_y);
		intersect = contains(p1) || contains(p2) || contains(rectangle.m_p1)
				|| contains(rectangle.m_p2);
	}
	return intersect;
}

DIntersectableRect DIntersectableRect::intersection(const DIntersectableRect& other) const {
	double top1 = m_p2.m_y;
	double bottom1 = m_p1.m_y;
	double left1 = m_p1.m_x;
	double right1 = m_p2.m_x;

	double top2 = other.m_p2.m_y;
	double bottom2 = other.m_p1.m_y;
	double left2 = other.m_p1.m_x;
	double right2 = other.m_p2.m_x;

	OGDF_ASSERT(top1 >= bottom1);
	OGDF_ASSERT(left1 <= right1);
	OGDF_ASSERT(top2 >= bottom2);
	OGDF_ASSERT(left2 <= right2);

	double bottomInter = max(bottom1, bottom2);
	double topInter = min(top1, top2);
	double leftInter = max(left1, left2);
	double rightInter = min(right1, right2);

	if (bottomInter > topInter) {
		return DIntersectableRect();
	}
	if (leftInter > rightInter) {
		return DIntersectableRect();
	}

	return DIntersectableRect(DPoint(leftInter, bottomInter), DPoint(rightInter, topInter));
}

// gives the segment starting at point 'it'
DSegment DPolygon::segment(ListConstIterator<DPoint> it) const {
	OGDF_ASSERT(!empty());
	OGDF_ASSERT(size() != 1);
	return DSegment(*it, *cyclicSucc(it));
}

// Assignment operator (for assigning from a rectangle).
DPolygon& DPolygon::operator=(const DRect& rect) {
	clear();
	DRect r1(rect);
	DRect r2(rect);
	if (m_counterclock) {
		r2.xInvert();
	} else {
		r2.yInvert();
	}

	pushBack(r1.p1());
	pushBack(r2.p1());
	pushBack(r1.p2());
	pushBack(r2.p2());

	unify();
	return *this;
}

// inserts the point p, which must ly on the boarder of the polygon, between the two points p1 and p2
// returns the index to that point, which is inserted only once
ListIterator<DPoint> DPolygon::insertPoint(const DPoint& p, ListIterator<DPoint> p1,
		ListIterator<DPoint> p2) {
	ListIterator<DPoint> i = p1;

	do {
		DSegment seg = segment(i);
		if (seg.contains(p)) {
			if (seg.start() == p) {
				return i;
			} else if (seg.end() == p) {
				i = cyclicSucc(i);
				return i;
			} else {
				return insertAfter(p, i);
			}
		}

		i = cyclicSucc(i);
	} while (i != p2);

	OGDF_ASSERT(false); // Point not in polygon, should not be reached!
	return i;
}

// inserts 'p' on every segment (a,b) with p in the open range ]a, b[
void DPolygon::insertCrossPoint(const DPoint& p) {
	ListIterator<DPoint> i = begin();

	do {
		DSegment seg = segment(i);
		if (seg.contains(p) && seg.start() != p && seg.end() != p) {
			i = insertAfter(p, i);
		}

		i = cyclicSucc(i);
	} while (i != begin());
}

//
int DPolygon::getCrossPoints(const DPolygon& p, List<DPoint>& crossPoints) const {
	crossPoints.clear();

	for (auto i = begin(); i.valid(); ++i) {
		DSegment s1 = segment(i);
		for (auto j = p.begin(); j.valid(); ++j) {
			DSegment s2 = p.segment(j);
			DPoint intersec;

			// TODO: What to do when IntersectionType::Overlapping is returned?
			if (s1.intersection(s2, intersec) == IntersectionType::SinglePoint) {
				crossPoints.pushBack(intersec);
			}
		}
	}
	// unify the list
	for (auto i = crossPoints.begin(); i.valid(); ++i) {
		for (auto j = i.succ(); j.valid(); ++j) {
			if (*i == *j) {
				--j;
				crossPoints.del(crossPoints.cyclicSucc(j));
			}
		}
	}

	return crossPoints.size();
}

// delete all consecutive double-points
void DPolygon::unify() {
	ListIterator<DPoint> iter, next;
	for (iter = begin(); iter.valid(); ++iter) {
		next = cyclicSucc(iter);
		while (*iter == *next) {
			del(next);
			next = cyclicSucc(iter);
			if (iter == next) {
				break;
			}
		}
	}
}

// deletes all points, which are not facets
void DPolygon::normalize() {
	unify();

	ListIterator<DPoint> iter, next;
	for (iter = begin(); iter.valid(); ++iter) {
		for (;;) {
			next = cyclicSucc(iter);
			DSegment s1 = segment(iter);
			DSegment s2 = segment(next);
			DRect r(*iter, *cyclicSucc(next));
			if (s1.slope() == s2.slope() && r.contains(*next)) {
				del(next);
			} else {
				break; // while
			}
		}
	}
}

// Checks wether a Point /a p is inside the Poylgon or not.
bool DPolygon::containsPoint(DPoint& p) const {
	if (size() < 3) {
		return false;
	}

	double angle = 0.0;
	DPolygon::const_iterator i = cyclicPred(begin());
	double lastangle = atan2((*i).m_y - p.m_y, (*i).m_x - p.m_x);
	for (const DPoint& q : *this) {
		double tempangle = atan2(q.m_y - p.m_y, q.m_x - p.m_x);
		double step = lastangle - tempangle;
		while (step > Math::pi) {
			step -= 2.0 * Math::pi;
		}
		while (step < -Math::pi) {
			step += 2.0 * Math::pi;
		}
		angle += step;
		lastangle = tempangle;
	}

	double d = angle / (2.0 * Math::pi);
	int rounds = static_cast<int>(d < 0 ? d - .5 : d + .5);

	return (rounds % 2) != 0;
}

// outputs the polygon
std::ostream& operator<<(std::ostream& os, const DPolygon& dop) {
	print(os, dop, ' ');
	return os;
}

int orientation(const DPoint& p, const DPoint& q, const DPoint& r) {
	double d1 = (p.m_x - q.m_x) * (p.m_y - r.m_y);
	double d2 = (p.m_y - q.m_y) * (p.m_x - r.m_x);

	if (d1 == d2) {
		return 0;
	} else {
		return (d1 > d2) ? +1 : -1;
	}
}

OGDF_EXPORT bool isPointCoveredByNode(const DPoint& point, const DPoint& v, const DPoint& vSize,
		const Shape& shape) {
	const double epsilon = 1e-6;
	const double trapeziumWidthOffset = vSize.m_x * 0.275;
	DPolyline polygon;

	auto isInConvexCCWPolygon = [&] {
		for (int i = 0; i < polygon.size(); i++) {
			DPoint edgePt1 = v + *polygon.get(i);
			DPoint edgePt2 = v + *polygon.get((i + 1) % polygon.size());

			if ((edgePt2.m_x - edgePt1.m_x) * (point.m_y - edgePt1.m_y)
							- (edgePt2.m_y - edgePt1.m_y) * (point.m_x - edgePt1.m_x)
					< -epsilon) {
				return false;
			}
		}
		return true;
	};

	auto isInRegularPolygon = [&](unsigned int sides) {
		polygon.clear();
		double radius = (max(vSize.m_x, vSize.m_y) / 2.0);
		for (unsigned int i = 0; i < sides; ++i) {
			double angle = -(Math::pi / 2) + Math::pi / sides + i * (2.0 * Math::pi / sides);
			polygon.pushBack(DPoint(radius * cos(angle), radius * sin(angle)));
		}
		return isInConvexCCWPolygon();
	};

	switch (shape) {
	// currently these tikz polygons are only supported as regular polygons, i.e. width=height
	case Shape::Pentagon:
		return isInRegularPolygon(5);
	case Shape::Hexagon:
		return isInRegularPolygon(6);
	case Shape::Octagon:
		return isInRegularPolygon(8);
	case Shape::Triangle:
		return isInRegularPolygon(3);
	// Non-regular polygons
	case Shape::InvTriangle:
		polygon.pushBack(DPoint(0, -vSize.m_y * 2.0 / 3.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0, vSize.m_y * 1.0 / 3.0));
		polygon.pushBack(DPoint(-vSize.m_x / 2.0, vSize.m_y * 1.0 / 3.0));
		return isInConvexCCWPolygon();
	case Shape::Rhomb:
		polygon.pushBack(DPoint(vSize.m_x / 2.0, 0));
		polygon.pushBack(DPoint(0, vSize.m_y / 2.0));
		polygon.pushBack(DPoint(-vSize.m_x / 2.0, 0));
		polygon.pushBack(DPoint(0, -vSize.m_y / 2.0));
		return isInConvexCCWPolygon();
	case Shape::Trapeze:
		polygon.pushBack(DPoint(-vSize.m_x / 2.0, -vSize.m_y / 2.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0, -vSize.m_y / 2.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0 - trapeziumWidthOffset, +vSize.m_y / 2.0));
		polygon.pushBack(DPoint(-vSize.m_x / 2.0 + trapeziumWidthOffset, +vSize.m_y / 2.0));
		return isInConvexCCWPolygon();
	case Shape::InvTrapeze:
		polygon.pushBack(DPoint(vSize.m_x / 2.0, vSize.m_y / 2.0));
		polygon.pushBack(DPoint(-vSize.m_x / 2.0, vSize.m_y / 2.0));
		polygon.pushBack(DPoint(-vSize.m_x / 2.0 + trapeziumWidthOffset, -vSize.m_y / 2.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0 - trapeziumWidthOffset, -vSize.m_y / 2.0));
		return isInConvexCCWPolygon();
	case Shape::Parallelogram:
		polygon.pushBack(DPoint(-vSize.m_x / 2.0, -vSize.m_y / 2.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0 - trapeziumWidthOffset, -vSize.m_y / 2.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0, +vSize.m_y / 2.0));
		polygon.pushBack(DPoint(-vSize.m_x / 2.0 + trapeziumWidthOffset, +vSize.m_y / 2.0));
		return isInConvexCCWPolygon();
	case Shape::InvParallelogram:
		polygon.pushBack(DPoint(-vSize.m_x / 2.0 + trapeziumWidthOffset, -vSize.m_y / 2.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0, -vSize.m_y / 2.0));
		polygon.pushBack(DPoint(vSize.m_x / 2.0 - trapeziumWidthOffset, vSize.m_y / 2.0));
		polygon.pushBack(DPoint(-vSize.m_x / 2.0, vSize.m_y / 2.0));
		return isInConvexCCWPolygon();
	// Ellipse
	case Shape::Ellipse:
		return pow((point.m_x - v.m_x) / (vSize.m_x * 0.5), 2)
				+ pow((point.m_y - v.m_y) / (vSize.m_y * 0.5), 2)
				< 1;
	// Simple x y comparison
	case Shape::Rect:
	case Shape::RoundedRect:
	default:
		return point.m_x + epsilon >= v.m_x - vSize.m_x / 2.0
				&& point.m_x - epsilon <= v.m_x + vSize.m_x / 2.0
				&& point.m_y + epsilon >= v.m_y - vSize.m_y / 2.0
				&& point.m_y - epsilon <= v.m_y + vSize.m_y / 2.0;
	}
}

DPoint contourPointFromAngle(double angle, int n, double rotationOffset, const DPoint& center,
		const DPoint& vSize) {
	// math visualised: https://www.desmos.com/calculator/j6iktd7fs4
	double nOffset = floor((angle - rotationOffset) / (2 * Math::pi / n)) * 2 * Math::pi / n;
	double polyLineStartAngle = rotationOffset + nOffset;
	double polyLineEndAngle = polyLineStartAngle + 2 * Math::pi / n;
	DLine polyLine = DLine(-cos(polyLineStartAngle), -sin(polyLineStartAngle),
			-cos(polyLineEndAngle), -sin(polyLineEndAngle));

	DLine originLine = DLine(0, 0, cos(angle), sin(angle));

	DPoint intersectionPoint;
	originLine.intersection(polyLine, intersectionPoint);
	intersectionPoint = DPoint(intersectionPoint.m_x * vSize.m_x, intersectionPoint.m_y * vSize.m_y);
	return intersectionPoint + center;
}

DPoint contourPointFromAngle(double angle, Shape shape, const DPoint& center, const DPoint& vSize) {
	angle = std::fmod(angle, 2 * Math::pi);
	if (angle < 0) {
		angle += Math::pi * 2;
	}

	switch (shape) {
	case Shape::Triangle:
		return contourPointFromAngle(angle, 3, Math::pi / 2, center, vSize * .5);
	case Shape::InvTriangle:
		return center - contourPointFromAngle(angle + Math::pi, Shape::Triangle, DPoint(), vSize);
	case Shape::Image:
	case Shape::RoundedRect:
	case Shape::Rect:
		return contourPointFromAngle(angle, 4, Math::pi / 4, center, vSize / sqrt(2));
	case Shape::Pentagon:
		return contourPointFromAngle(angle, 5, Math::pi / 2, center, vSize / 2);
	case Shape::Hexagon:
		return contourPointFromAngle(angle, 6, 0, center, vSize / 2);
	case Shape::Octagon:
		return contourPointFromAngle(angle, 8, Math::pi / 8, center, vSize / 2);
	case Shape::Rhomb:
		return contourPointFromAngle(angle, 4, Math::pi / 2, center, vSize / 2);
	case Shape::Trapeze:
		if (angle < atan(2) || angle >= Math::pi * 7 / 4) {
			DPoint other = contourPointFromAngle(Math::pi - angle, Shape::Trapeze, DPoint(), vSize);
			other.m_x *= -1;
			return other + center;
		} else if (angle < Math::pi - atan(2)) {
			return contourPointFromAngle(angle, Shape::Rect, center, vSize);
		} else if (angle < Math::pi * 5 / 4) {
			DLine tLine = DLine(.5, -1, 1, 1);
			DLine eLine = DLine(0, 0, 2 * cos(angle), 2 * sin(angle));
			DPoint iPoint;
			tLine.intersection(eLine, iPoint);
			iPoint = DPoint(iPoint.m_x * vSize.m_x * .5, iPoint.m_y * vSize.m_y * .5);
			return iPoint + center;
		} else { // angle < Math::pi * 7 / 4
			return contourPointFromAngle(angle, Shape::Rect, center, vSize);
		}
	case Shape::InvTrapeze:
		return center - contourPointFromAngle(angle + Math::pi, Shape::Trapeze, DPoint(), vSize);
	case Shape::Parallelogram:
		if (angle < atan(2) || angle > Math::pi * 7 / 4) {
			DLine tLine = DLine(-.5, -1, -1, 1);
			DLine eLine = DLine(0, 0, 2 * cos(angle), 2 * sin(angle));
			DPoint iPoint;
			tLine.intersection(eLine, iPoint);
			iPoint = DPoint(iPoint.m_x * vSize.m_x * .5, iPoint.m_y * vSize.m_y * .5);
			return iPoint + center;
		} else if (angle < Math::pi * 3 / 4) {
			return contourPointFromAngle(angle, Shape::Rect, center, vSize);
		} else if (angle < Math::pi + atan(2)) {
			DLine tLine = DLine(.5, 1, 1, -1);
			DLine eLine = DLine(0, 0, 2 * cos(angle), 2 * sin(angle));
			DPoint iPoint;
			tLine.intersection(eLine, iPoint);
			iPoint = DPoint(iPoint.m_x * vSize.m_x * .5, iPoint.m_y * vSize.m_y * .5);
			return iPoint + center;
		} else { // angle < Math::pi * 7 / 4
			return contourPointFromAngle(angle, Shape::Rect, center, vSize);
		}
	case Shape::InvParallelogram: {
		DPoint p = contourPointFromAngle(Math::pi - angle, Shape::Parallelogram, DPoint(), vSize);
		p.m_x *= -1;
		return p + center;
	}
	case Shape::Ellipse:
	default:
		return DPoint(-vSize.m_x * .5 * cos(angle), -vSize.m_y * .5 * sin(angle)) + center;
	}
}

}
