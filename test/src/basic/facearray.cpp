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

#include "array_helper.h"

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

	Graph G;

	auto init = [&](CombinatorialEmbedding& C) {
		randomPlanarConnectedGraph(G, 42, 168);
		C.init(G);
	};

	describeArray<CombinatorialEmbedding, FaceArray, face, int>("FaceArray filled with ints", 42,
			43, init, chooseFace, allFaces, createFace);
	describeArray<CombinatorialEmbedding, FaceArray, face, List<int>>(
			"FaceArray filled with lists of ints", {1, 2, 3}, {42}, init, chooseFace, allFaces,
			createFace);
});
