/*
 * $Revision: 3005 $
 *
 * last checkin:
 *   $Author: chimani $
 *   $Date: 2012-11-12 14:19:48 +0100 (Mo, 12. Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Declaration of a constraint class for the Branch&Cut algorithm
 * for the Maximum C-Planar SubGraph problem.
 *
 * These constraints represent the planarity-constraints belonging to the
 * ILP formulation. These constraints are dynamically separated.
 * For the separation the planarity test algorithm by Boyer and Myrvold is used.
 *
 * \author Mathias Jansen
 *
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

#ifndef OGDF_KURATOWSKI_CONSTRAINT_H
#define OGDF_KURATOWSKI_CONSTRAINT_H

#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/internal/cluster/Cluster_EdgeVar.h>
#include <ogdf/internal/cluster/basics.h>
//#include <ogdf/internal/cluster/MaxCPlanar_Master.h>
//#include <ogdf/abacus/master.h>

#include <ogdf/abacus/constraint.h>

namespace ogdf {


class ClusterKuratowskiConstraint : public abacus::Constraint {

public:

	ClusterKuratowskiConstraint(abacus::Master *master, int nEdges, SListPure<nodePair> &ks);

	virtual ~ClusterKuratowskiConstraint();

	// Computes and returns the coefficient for the given variable
	virtual double coeff(const abacus::Variable *v) const;

	void printMe(ostream& out) const {
		out << "[KuraCon: ";
		forall_listiterators(nodePair, it, m_subdivision) {
			(*it).printMe(out);
			out << ",";
		}
		out << "]";
	}

private:

	// The subdivision containing edges forming a SubGraph that is not planar
	List<nodePair> m_subdivision;

};

}

#endif
