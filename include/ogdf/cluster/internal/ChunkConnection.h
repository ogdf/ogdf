/** \file
 * \brief implementation of initial cut-constraint class for the Branch&Cut algorithm
 * for the Maximum C-Planar SubGraph problem.
 *
 * A feasible ILP solution has to imply a completely connected, planar Sub-Clustergraph.
 * For each cluster that is not connected, additional connection edges have to be inserted
 * between the chunks of the cluster, to obtain c-connectivity.
 * Thus, initial constraints are added that guarantee initial c-connectivity, if the number of chunks
 * is at most 3. If some cluster consists of more than 3 chunks, additional constraints
 * have to be separated during the optimization.
 *
 * \author Mathias Jansen
 *
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

#pragma once

#include <ogdf/basic/Array.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/cluster/internal/EdgeVar.h>
#include <ogdf/cluster/internal/basics.h>

#include <ostream>

namespace abacus {
class Master;
class Variable;
} // namespace abacus

namespace ogdf {
template<class E, class INDEX>
class ArrayBuffer;

namespace cluster_planarity {

class ChunkConnection : public BaseConstraint {
#ifdef OGDF_DEBUG
	//Mainly for debugging output purposes
	friend class MaxCPlanarMaster;
	friend class MaxCPlanarSub;
	friend class CPlanarMaster;
	friend class CPlanarSub;
#endif
public:
	ChunkConnection(abacus::Master* master, const ArrayBuffer<node>& chunk,
			const ArrayBuffer<node>& cochunk);

	virtual ~ChunkConnection();

	// Computes and returns the coefficient for the given variable
	virtual double coeff(const abacus::Variable* v) const override {
		const EdgeVar* ev = static_cast<const EdgeVar*>(v);
		//Safe for both clustered planarity testing and maximum c-planar subgraph
		return (ev->theEdgeType() != EdgeVar::EdgeType::Connect)
				? 0.0
				: (double)coeff(ev->sourceNode(), ev->targetNode());
	}

	inline int coeff(const NodePair& n) const override { return coeff(n.source, n.target); }

	int coeff(node v1, node v2) const;

	void printMe(std::ostream& out) const {
		out << "[ChunkCon: (";
		for (node v : m_chunk) {
			Logger::slout() << v << ",";
		}
		out << "|";
		for (node v : m_cochunk) {
			Logger::slout() << v << ",";
		}
		out << ")]";
	}

private:
	// The nodePairs corresponding to the constraint
	Array<node> m_chunk;
	Array<node> m_cochunk;
};

}
}
