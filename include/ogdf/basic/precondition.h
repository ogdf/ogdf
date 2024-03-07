/** \file
 * \brief Declaration of functions for drawing module precondition
 *        handling.
 *
 * \author Karsten Klein
 *
 * \attention This is legacy code from UML class diagram handling,
 * and it should be checked if it is still required.
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

#include <ogdf/orthogonal/EdgeRouter.h>
#include <ogdf/uml/UMLGraph.h>

namespace ogdf {

//descent the hierarchy tree at "sink" v recursively
bool dfsGenTreeRec(UMLGraph& UG, EdgeArray<bool>& used,
		NodeArray<int>& hierNumber, //number of hierarchy tree
		// A node is visited if its hierNumber != 0
		int hierNum, node v,
		List<edge>& fakedGens, //temporary
		bool fakeTree);

edge firstOutGen(UMLGraph& UG, node v, EdgeArray<bool>& /* used */);

bool dfsGenTree(UMLGraph& UG, List<edge>& fakedGens, bool fakeTree);

}
