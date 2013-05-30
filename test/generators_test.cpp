/*
 * $Revision$
 *
 * last checkin:
 *   $Author$
 *   $Date$
 ***************************************************************/

/** \file
 * \brief Implementation of a command line based tool to run
 * tests.
 *
 * \author Christoph Schulz
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

#include "gtest/gtest.h"
#include "ogdf/basic/Graph.h"
#include "ogdf/basic/graph_generators.h"
#include <cstdlib>

double randomValue()
{
	return static_cast<double>(rand()) / RAND_MAX;
}

void randomGraph(ogdf::Graph &G, size_t index)
{
    int n = static_cast<int>(randomValue() * 10 + 5);
    int m = std::max(n + 1, static_cast<int>((randomValue() * n * (n - 1) / 2)));
	double p = 0.618;
	switch (index) {
	case 0:
		G.clear();
		ogdf::randomGraph(G, n, m);
		break;
	case 1:
		G.clear();
		ogdf::randomSimpleGraph(G, n, m);
		break;
	case 2:
		G.clear();
		ogdf::randomBiconnectedGraph(G, n, m);
		break;
	case 3:
		G.clear();
		ogdf::randomTriconnectedGraph(G, n, p, 1.0 - p);
		break;
	case 4:
		G.clear();
		ogdf::randomTree(G, n);
		break;
	case 5:
		G.clear();
		ogdf::randomTree(G, n, n / 2, 1024);
		break;
	case 6:
		G.clear();
		ogdf::randomHierarchy(G, n, m, true, false, true);
		break;
	case 7:
		G.clear();
		ogdf::randomDiGraph(G, n, p);
		break;
	default:
		return;
	}
}

TEST(GeneratorsTest, RandomGraph)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 0);
}

TEST(GeneratorsTest, RandomSimpleGraph)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 1);
}

TEST(GeneratorsTest, RandomBiconnectedGraph)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 2);
}

TEST(GeneratorsTest, RandomTriconnectedGraph)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 3);
}

TEST(GeneratorsTest, RandomTree1)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 4);
}

TEST(GeneratorsTest, RandomTree2)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 5);
}

TEST(GeneratorsTest, RandomHierarchy)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 6);
}

TEST(GeneratorsTest, RandomDiGraph)
{
	srand(8609);
	ogdf::Graph G;
	randomGraph(G, 7);
}
