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

static unsigned int randomSeed = 8609;

//void randomize(unsigned int seed)
//{
//	randomSeed = seed + 123456;
//}

double randomNumber()
{
	const unsigned int a = 16807;
	const unsigned int c = 0;
	const unsigned int m = 0x7fffffff;

	randomSeed = (a * randomSeed + c) & m;
	return static_cast<double>(randomSeed) / m;
}

TEST(GeneratorsTest, RandomGraph)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	const int m = std::max(n + 1, static_cast<int>((randomNumber() * n * (n - 1) / 2)));
	ogdf::Graph G;
	ogdf::randomGraph(G, n, m);
}

TEST(GeneratorsTest, RandomSimpleGraph)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	const int m = std::max(n + 1, static_cast<int>((randomNumber() * n * (n - 1) / 2)));
	ogdf::Graph G;
	ogdf::randomSimpleGraph(G, n, m);
}

TEST(GeneratorsTest, RandomBiconnectedGraph)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	const int m = std::max(n + 1, static_cast<int>((randomNumber() * n * (n - 1) / 2)));
	ogdf::Graph G;
	ogdf::randomBiconnectedGraph(G, n, m);
}

TEST(GeneratorsTest, RandomTriconnectedGraph)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	double p = 0.618;
	ogdf::Graph G;
	ogdf::randomTriconnectedGraph(G, n, p, 1.0 - p);
}

TEST(GeneratorsTest, RandomTree1)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	ogdf::Graph G;
	ogdf::randomTree(G, n);
}

TEST(GeneratorsTest, RandomTree2)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	ogdf::Graph G;
	ogdf::randomTree(G, n, n / 2, 1024);
}

TEST(GeneratorsTest, RandomHierarchy)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	const int m = std::max(n + 1, static_cast<int>((randomNumber() * n * (n - 1) / 2)));
	ogdf::Graph G;
	ogdf::randomHierarchy(G, n, m, true, false, true);
}

TEST(GeneratorsTest, RandomDiGraph)
{
	const int n = static_cast<int>(randomNumber() * 10 + 5);
	double p = 0.618;
	ogdf::Graph G;
	ogdf::randomDiGraph(G, n, p);
}
