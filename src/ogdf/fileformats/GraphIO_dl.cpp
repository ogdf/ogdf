/*
 * $Revision: 3837 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-11-13 15:19:30 +0100 (Mi, 13. Nov 2013) $
 ***************************************************************/

/** \file
 * \brief Implements UCINET DL write functionality of class GraphIO.
 *
 * \author Łukasz Hanuszczak
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

#include <vector>

namespace ogdf {


static void writeMatrix(
	std::ostream &os,
	const Graph &G, const GraphAttributes *GA)
{
	os << "DATA:\n";
	const long attrs = GA ? GA->attributes() : 0;
	const int n = G.numberOfNodes();
	std::vector<double> matrix(n * n, 0);

	edge e;
	forall_edges(e, G) {
		const int vs = e->source()->index();
		const int vt = e->target()->index();

		if(attrs & GraphAttributes::edgeDoubleWeight) {
			matrix[vs * n + vt] = GA->doubleWeight(e);
		} else if(attrs & GraphAttributes::edgeIntWeight) {
			matrix[vs * n + vt] = GA->intWeight(e);
		} else {
			matrix[vs * n + vt] = 1;
		}
	}

	node v, u;
	forall_nodes(v, G) {
		bool space = false;
		forall_nodes(u, G) {
			if(space) {
				os << " ";
			}
			space = true;

			const int vs = v->index(), vt = u->index();
			os << matrix[vs * n + vt];
		}
		os << "\n";
	}
}


static void writeEdges(
	std::ostream &os,
	const Graph &G, const GraphAttributes *GA)
{
	os << "DATA:\n";
	const long attrs = GA ? GA->attributes() : 0;

	edge e;
	forall_edges(e, G) {
		os << (e->source()->index() + 1) << " " << (e->target()->index() + 1);

		if(attrs & GraphAttributes::edgeDoubleWeight) {
			os << " " << GA->doubleWeight(e);
		} else if(attrs & GraphAttributes::edgeIntWeight) {
			os << " " << GA->intWeight(e);
		}

		os << "\n";
	}
}


static void writeGraph(
	std::ostream &os,
	const Graph &G, const GraphAttributes *GA)
{
	const long long n = G.numberOfNodes(), m = G.numberOfEdges();

	os << "DL N = " << n << "\n";

	// We pick output format basing on edge density.
	enum { matrix, edges } format = (m > (n * n / 2)) ? matrix : edges;

	// Specify output format.
	os << "FORMAT = ";
	if(format == matrix) {
		os << "fullmatrix\n";
		writeMatrix(os, G, GA);
	} else if(format == edges) {
		os << "edgelist1\n";
		writeEdges(os, G, GA);
	}
}


bool GraphIO::writeDL(const Graph &G, std::ostream &os)
{
	writeGraph(os, G, NULL);
	return true;
}


bool GraphIO::writeDL(const GraphAttributes &GA, std::ostream &os)
{
	writeGraph(os, GA.constGraph(), &GA);
	return true;
}


} // end namespace ogdf

