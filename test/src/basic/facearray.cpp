/** \file
 * \brief Tests for ogdf::FaceArray
 *
 * \author Matthias Pfretzschner
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
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/FaceSet.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/basic.h>
#include <ogdf/basic/graph_generators/randomized.h>

#include <functional>
#include <list>
#include <string>

#include "array_helper.h"
#include <testing.h>

go_bandit([]() {
	auto chooseFace = [](const CombinatorialEmbedding& C) { return C.chooseFace(); };

	auto allFaces = [](const CombinatorialEmbedding& C, List<face>& list) {
		list.clear();
		for (face f : C.faces) {
			list.pushBack(f);
		}
	};

	auto createFace = [](CombinatorialEmbedding& C) {
		adjEntry a1 = C.lastFace()->firstAdj();
		adjEntry a2 = C.lastFace()->nextFaceEdge(a1);
		C.splitFace(a1, a2);
		return C.lastFace();
	};
	auto deleteFace = [](CombinatorialEmbedding& C, face f) {
		face o = C.joinFaces(f->firstAdj()->twin());
		OGDF_ASSERT(o != f);
	};
	auto clearFaces = [](CombinatorialEmbedding& C) { C.clear(); };

	std::list<Graph> graphs;
	auto init = [&](CombinatorialEmbedding& C) {
		graphs.emplace_back();
		randomPlanarConnectedGraph(graphs.back(), 42, 168);
		C.init(graphs.back());
	};

	runBasicArrayTests<CombinatorialEmbedding, FaceArray, face>( //
			"FaceArray", init, chooseFace, allFaces, createFace);

	graphs.clear();
	runBasicSetTests<CombinatorialEmbedding, FaceSet, face>("FaceSet", init, chooseFace, allFaces,
			createFace, deleteFace, clearFaces);
});
