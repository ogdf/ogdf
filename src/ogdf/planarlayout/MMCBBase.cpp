/** \file
 * \brief Implementation of Mixed-Model crossings beautifiers.
 *
 * \author Carsten Gutwenger
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

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/GridLayout.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Math.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/planarity/PlanRep.h>
#include <ogdf/planarlayout/MMCBBase.h>
#include <ogdf/planarlayout/MMCBDoubleGrid.h>
#include <ogdf/planarlayout/MMCBLocalStretch.h>

#include <functional>

namespace ogdf {

void MMCBBase::insertBend(GridLayout& gl, edge e, node v, int x, int y) {
	if (v == e->target()) {
		gl.bends(e).pushBack(IPoint(x, y));
	} else {
		gl.bends(e).pushFront(IPoint(x, y));
	}
}

void MMCBBase::copyOn(int old_a[], int new_a[]) {
	for (int i = 0; i < 3; ++i) {
		new_a[i] = old_a[i];
	}
}

const static int CHANGE_NONE = 0;
const static int CHANGE_X = 1;
const static int CHANGE_Y = 2;

int MMCBBase::workOn(GridLayout& gl, node v) {
	edge L[4];
	int ev[4][3];
	int count = 0;

	for (adjEntry adj : v->adjEntries) {
		edge e = adj->theEdge();

		IPolyline& ip = gl.bends(e);
		int xc, yc;
		int xxc = gl.x(v), yyc = gl.y(v);
		int add_l = 0;
		do {
			if (ip.size() < (1 + add_l)) {
				if (v == e->target()) {
					xc = gl.x(e->source());
					yc = gl.y(e->source());
				} else {
					xc = gl.x(e->target());
					yc = gl.y(e->target());
				}
			} else {
				IPoint p;
				if (v == e->target()) {
					p = *ip.get(ip.size() - add_l - 1);
				} else {
					p = *ip.get(1 + add_l - 1);
				}
				xc = p.m_x;
				yc = p.m_y;
			}
			++add_l;
		} while ((xc == xxc) && (yc == yyc));

		if (xc < gl.x(v)) {
			ev[count][0] = -1;
		} else if (xc == gl.x(v)) {
			ev[count][0] = 0;
		} else {
			ev[count][0] = 1;
		}

		if (yc < gl.y(v)) {
			ev[count][1] = -1;
		} else if (yc == gl.y(v)) {
			ev[count][1] = 0;
		} else {
			ev[count][1] = 1;
		}

		L[count] = e;
		++count;
		if (3 < count) {
			break;
		}
	}
	OGDF_ASSERT(count == 4);

	for (auto& elem : ev) {
		if (elem[0] > 0) {
			elem[2] = 4 + elem[1] + 2;

		} else if (elem[0] == 0) {
			if (elem[1] > 0) {
				elem[2] = 0;
			} else {
				elem[2] = 4;
			}

		} else {
			elem[2] = -elem[1] + 2;
		}
	}

	int ew[4][3];
	edge Lw[4];

	int k = 4;
	for (int i = 0; (i < 8) && (k >= 0); i++) {
		for (int j = 0; j < 4; j++) {
			if (i == ev[j][2]) {
				copyOn(ev[j], ew[4 - k]);
				edge edge_k = L[j];
				Lw[4 - k] = edge_k;
				--k;
			}
		}
	}

	//so now in ew the edges are in an order against the clock starting at position 0,1
	int ea, eb;
	int crossingCase = 0;

	ea = ew[2][2] - ew[0][2];
	if (ea > 4) {
		ea = 8 - ea;
	}
	eb = ew[3][2] - ew[1][2];
	if (eb > 4) {
		eb = 8 - eb;
	}

	// first case: there is an edge with angle of pi/2
	// in this case i'm going to, well, cut this edge
	//   e4 O      _/e3             Ooo     _/
	//      O    _/                  Ooo _/
	//   \  O  _/                \     OXo
	//      O /                      _/  Ooo
	//  --  *ooooooo^   ==>    --   /      Ooo
	//             e2
	//   /  |  \                 /  |   \_
	//  e1

	int lw_add = 0;
	if (ea <= 2) {
		if (ew[2][2] - ew[0][2] > 4) {
			lw_add = 2;
			crossingCase = 1;
		} else {
			crossingCase = 1;
		}
	} else if (eb <= 2) {
		if (ew[3][2] - ew[1][2] > 4) {
			lw_add = 1;
			crossingCase = 1;
		} else {
			lw_add = 3;
			crossingCase = 1;
		}

		// second case: both angles are 3pi/4, in this case
		// i'm going to take both of them back with a row
		//   ^     _/                 ooO      _/
		//   O   _/                 ooO      _/
		//   O _/                  Ooo     _/
		//   O/                      Ooo _/
		//   *o          ==>           OXo
		//   |Ooo                     _/ Ooo
		//   |  Oo                  _/     Ooo
		//   |   Ooo                \_       Ooo
		//   |     Oo                 \_       Oo
		//
	} else if (ea == 3) {
		if (eb == 3) {
			if (ew[1][2] - ew[0][2] >= 3) {
				if (ew[1][2] - ew[0][2] == 3) {
					crossingCase = 2;
				} else {
					crossingCase = 3;
				}
			} else if (ew[2][2] - ew[1][2] >= 3) {
				lw_add = 3;
				if (ew[2][2] - ew[1][2] == 3) {
					crossingCase = 2;
				} else {
					crossingCase = 3;
				}
			} else if (ew[3][2] - ew[2][2] >= 3) {
				lw_add = 2;
				if (ew[3][2] - ew[2][2] == 3) {
					crossingCase = 2;
				} else {
					crossingCase = 3;
				}
			} else {
				lw_add = 1;
				if (ew[3][2] - ew[0][2] == 3) {
					crossingCase = 2;
				} else {
					crossingCase = 3;
				}
			}

			// third case: one of the angles is 3pi/4 and the other one is pi
			//          ^                            O
			//          O                            O
			//          O                            O
			//          O                            O
			//  --------*o-------    ==>    ---------X--------
			//           Ooo                         O
			//             Ooo                       O
			//               Ooo                     O
			//                 Oo                    Ooooooooo
			//
		} else {
			if (ew[1][2] - ew[0][2] == 2) {
				crossingCase = 4;
			} else {
				lw_add = 2;
				crossingCase = 4;
			}
		}
	} else if (eb == 3) {
		if (ew[2][2] - ew[1][2] == 2) {
			lw_add = 3;
			crossingCase = 4;
		} else {
			lw_add = 1;
			crossingCase = 4;
		}
	}

	// v, Lw, ev
	copyOn(ew[(0 + lw_add) % 4], ev[0]);
	copyOn(ew[(1 + lw_add) % 4], ev[1]);
	copyOn(ew[(2 + lw_add) % 4], ev[2]);
	copyOn(ew[(3 + lw_add) % 4], ev[3]);

	edge e0 = Lw[(0 + lw_add) % 4];
	edge e1 = Lw[(1 + lw_add) % 4];
	edge e2 = Lw[(2 + lw_add) % 4];
	edge e3 = Lw[(3 + lw_add) % 4];

	int retVal = CHANGE_NONE;

	switch (crossingCase) {
	case 4:
		if ((ev[0][0] * ev[0][1] == 0) || (gl.bends(e3).size() > 0)) {
			insertBend(gl, e3, v, gl.x(v) - ev[1][0], gl.y(v) - ev[1][1]);
		} else {
			insertBend(gl, e0, v, gl.x(v) + ev[3][0], gl.y(v) + ev[3][1]);
			insertBend(gl, e2, v, gl.x(v) - ev[3][0], gl.y(v) - ev[3][1]);
			insertBend(gl, e1, v, gl.x(v) + ev[0][0] + ev[3][0], gl.y(v) + ev[0][1] + ev[3][1]);
		}
		break;

	case 3:
		if (ev[0][0] * ev[0][1] != 0) {
			insertBend(gl, e0, v, gl.x(v) - ev[2][0], gl.y(v) - ev[2][1]);
			insertBend(gl, e1, v, gl.x(v) - ev[3][0], gl.y(v) - ev[3][1]);
		} else {
			gl.x(v) = gl.x(v) + ev[2][0] - ev[1][0];
			gl.y(v) = gl.y(v) + ev[2][1] - ev[1][1];
			if (ev[2][0] - ev[1][0] == 0) {
				retVal = CHANGE_Y;
			} else {
				retVal = CHANGE_X;
			}
		}
		break;

	case 2:
		insertBend(gl, e1, v, gl.x(v) - ev[3][0], gl.y(v) - ev[3][1]);
		insertBend(gl, e2, v, gl.x(v) - ev[0][0], gl.y(v) - ev[0][0]);
		break;

	case 1:
		if (ev[0][0] * ev[0][1] == 0) {
			int old_x = gl.x(v);
			int old_y = gl.y(v);
			int x_plus = (ev[0][0] + ev[2][0]);
			int y_plus = (ev[0][1] + ev[2][1]);
			insertBend(gl, e0, v, old_x + 2 * ev[0][0], old_y + 2 * ev[0][1]);
			insertBend(gl, e2, v, old_x + 2 * ev[2][0], old_y + 2 * ev[2][1]);
			insertBend(gl, e1, v, old_x, old_y);
			insertBend(gl, e3, v, old_x, old_y);
			gl.x(v) = old_x + x_plus;
			gl.y(v) = old_y + y_plus;
			retVal = CHANGE_X | CHANGE_Y;
		} else {
			int old_x = gl.x(v);
			int old_y = gl.y(v);

			gl.x(v) = gl.x(v) + ev[1][0];
			gl.y(v) = gl.y(v) + ev[1][1];
			insertBend(gl, e3, v, old_x, old_y);
			insertBend(gl, e1, v, gl.x(v) + ev[1][0], gl.y(v) + ev[1][1]);
			if (ev[1][0] != 0) {
				retVal = CHANGE_X;
			} else {
				retVal = CHANGE_Y;
			}
		}
		break;

	case 0:
		break;

	default:
		OGDF_ASSERT(false);
	}

	return retVal;
}

static void doForEachCoordinate(const PlanRep& PG, GridLayout& gl,
		std::function<void(int&, int&)> func) {
	for (edge e : PG.edges) {
		for (IPoint& p : gl.bends(e)) {
			func(p.m_x, p.m_y);
		}
	}

	for (node v : PG.nodes) {
		func(gl.x(v), gl.y(v));
	}
}

void MMCBDoubleGrid::doCall(const PlanRep& PG, GridLayout& gl, const List<node>& L) {
	doForEachCoordinate(PG, gl, [](int& x, int& y) {
		x *= 2;
		y *= 2;
	});

	for (node v : L) {
		workOn(gl, v);
	}
}

static void fillSignum(Array<int>& array) {
	for (int i = array.low(); i <= array.high(); i++) {
		array[i] = Math::sgn(i);
	}
}

void MMCBLocalStretch::doCall(const PlanRep& PG, GridLayout& gl, const List<node>& L) {
	int max_x = 0, max_y = 0;
	int min_x = 0, min_y = 0;

	doForEachCoordinate(PG, gl, [&](int& x, int& y) {
		Math::updateMax(max_x, x);
		Math::updateMin(min_x, x);
		Math::updateMax(max_y, y);
		Math::updateMin(min_y, y);
		x *= 2;
		y *= 2;
	});

	Array<int> change_x(min_x, max_x);
	fillSignum(change_x);
	Array<int> change_y(min_y, max_y);
	fillSignum(change_y);

	auto index = [](int pos) { return pos >= 0 ? (pos + 1) / 2 : (pos - 1) / 2; };

	for (node v : L) {
		int val = workOn(gl, v);
		if (val & CHANGE_X) {
			change_x[index(gl.x(v))] = 0;
		}
		if (val & CHANGE_Y) {
			change_y[index(gl.y(v))] = 0;
		}
	}

	for (int i = -1; i >= min_x; i--) {
		change_x[i] += change_x[i + 1];
	}
	for (int i = 1; i <= max_x; i++) {
		change_x[i] += change_x[i - 1];
	}
	for (int i = -1; i >= min_y; i--) {
		change_y[i] += change_y[i + 1];
	}
	for (int i = 1; i <= max_y; i++) {
		change_y[i] += change_y[i - 1];
	}

	doForEachCoordinate(PG, gl, [&](int& x, int& y) {
		x -= change_x[index(x)];
		y -= change_y[index(y)];
	});
}

}
