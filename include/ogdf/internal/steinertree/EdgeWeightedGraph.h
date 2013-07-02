/*
 * $Revision: 3609 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-07-01 20:33:08 +0200 (Mo, 01. Jul 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class WeightedGraph
 *
 * \author Matthias Woste
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

#ifndef OGDF_EDGE_WEIGHTED_GRAPH_H_
#define OGDF_EDGE_WEIGHTED_GRAPH_H_

#include <ogdf/basic/GraphCopy.h>

namespace ogdf {

template<typename T>
class EdgeWeightedGraph: public Graph {
public:
	EdgeWeightedGraph()
	  : Graph()
	  , m_edgeWeight(*this)
	{
	}

	EdgeWeightedGraph(GraphCopy &gC)
	{
	}

	virtual ~EdgeWeightedGraph()
	{
	}

	edge newEdge(const node v, const node w, T weight)
	{
		edge e = Graph::newEdge(v, w);
		m_edgeWeight[e] = weight;
		return e;
	}

	node newNode()
	{
		node u = Graph::newNode();
		return u;
	}

	T weight(const edge e) const
	{
		return m_edgeWeight[e];
	}

	EdgeArray<T> edgeWeights() const
	{
		return m_edgeWeight;
	}

protected:
	EdgeArray<T> m_edgeWeight;
};

}

#endif /* OGDF_EDGE_WEIGHTED_GRAPH_H_ */
