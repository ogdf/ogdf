/** \file
 * \brief Tests for SugiyamaLayout
 *
 * \author Carsten Gutwenger, Tilo Wiedera
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

#include <bandit/bandit.h>

#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/MedianHeuristic.h>
#include <ogdf/layered/OptimalHierarchyLayout.h>
#include <ogdf/layered/OptimalRanking.h>
#include <ogdf/layered/GreedyCycleRemoval.h>

#include "layout_helpers.h"

using namespace ogdf;

go_bandit([](){ bandit::describe("Sugiyama layouts", [](){
	SugiyamaLayout sugi, sugiOpt, sugiTrans, sugiRuns;

	sugi.setLayout(new FastHierarchyLayout);
	describeLayoutModule("Sugiyama with fast hierarchy", sugi, 0, {}, 100);

	std::string desc = "Sugiyama with optimal ranking, median";
	sugiOpt.setLayout(new OptimalHierarchyLayout);
	desc += " and optimal hierarchy";
	OptimalRanking *optr = new OptimalRanking;
	optr->setSubgraph(new GreedyCycleRemoval);
	sugiOpt.setRanking(optr);
	sugiOpt.setCrossMin(new MedianHeuristic);
	sugiOpt.transpose(false);
	describeLayoutModule(desc.c_str(), sugiOpt, 0, {}, 50);

	sugiTrans.transpose(true);
	describeLayoutModule("Sugiyama with transpositions", sugiTrans, 0, {}, 100);

	sugiRuns.runs(40);
	describeLayoutModule("Sugiyama with 40 runs", sugiRuns, 0, {}, 50);
}); });
