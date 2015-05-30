/** \file
 * \brief Tests for Minisat wrapper
 *
 * \author Robert Zeranski, Stephan Beyer
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

#include <bandit/bandit.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/minisat/Minisat.h>
#include <iostream>

using namespace bandit;
using namespace ogdf;

static void satisfiableTest()
{
	Minisat::Formula F;
	Minisat::Model model;
	Minisat::clause cl = F.newClause();

	cl->addMultiple(4,-1,-2,-3,4);
	//cl->writeToConsole();
	F.finalizeClause(cl);

	F.newVars(11);

	for (int i = 1; i < 10; i++) {
			Minisat::clause c = F.newClause();
			if (i % 2 == 0) {
				c->add(i);
			}
			else {
				c->add(-i);
			}
			c->add(i+1);
			F.finalizeClause(c);
	}

	bool satisfiable = F.solve(model);

	AssertThat(satisfiable, IsTrue());
#if 0
	cout << "#vars = " << F.getVariableCount() << endl;
	cout << "#clauses = " << F.getClauseCount() << endl;
	if (val) {
		model.printModel();
	}
	F.reset();
#endif
#if 0
	F.writeFormulaToDimacs("example.cnf");
#endif
}

static void nonsatisfiableTest()
{
	Minisat::Formula F;
	Minisat::Model model;
	Minisat::clause cl;
	cl = F.newClause();
	cl->addMultiple(2, 1, 2);
	F.finalizeClause(cl);
	cl = F.newClause();
	cl->addMultiple(2, 1, -2, 3);
	F.finalizeClause(cl);
	cl = F.newClause();
	cl->addMultiple(2, -1, 2);
	F.finalizeClause(cl);
	cl = F.newClause();
	cl->addMultiple(2, -1, -2);
	F.finalizeClause(cl);
	cl = F.newClause();
	cl->addMultiple(1, -3);
	F.finalizeClause(cl);

	bool satisfiable = F.solve(model);

	AssertThat(satisfiable, IsFalse());
}

go_bandit([]() {
	describe("Minisat wrapper", []() {
		it("solves a satisfiable formula", []() {
			satisfiableTest();
		});
		it("solves a non-satisfiable formula", []() {
			nonsatisfiableTest();
		});
	});
});
