/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 14:49:39 +0200 (Do, 16. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the class UpwardPlanarity.
 *
 * \author Robert Zeranski
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

#include <ogdf/upward/UpwardPlanarity.h>

#include <ogdf/internal/upward/UpwardPlanarityEmbeddedDigraph.h>
#include <ogdf/internal/upward/UpwardPlanaritySingleSource.h>

#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/upward/FaceSinkGraph.h>


namespace ogdf {

	//
	// General digraphs
	//

	bool UpwardPlanarity::isUpwardPlanar_embedded(const Graph &G)
	{
			if (G.representsCombEmbedding() && isAcyclic(G)) {
				UpwardPlanarityEmbeddedDigraph p(G);
				return p.isUpwardPlanarEmbedded();
			}
			return false;
	}


	bool UpwardPlanarity::isUpwardPlanar_embedded(const Graph &G, List<adjEntry> &possibleExternalFaces)
	{
			if (G.representsCombEmbedding() && isAcyclic(G)) {
				UpwardPlanarityEmbeddedDigraph p(G);
				return p.isUpwardPlanarEmbedded(possibleExternalFaces);
			}
			return false;
	}


	//
	// Triconnected digraphs
	//


	bool UpwardPlanarity::isUpwardPlanar_triconnected(const Graph &G)
	{
		if (isTriconnected(G) && isAcyclic(G)) {
			Graph H(G);
			BoyerMyrvold p;
			if (!p.planarEmbed(H)) return false;
			return isUpwardPlanar_embedded(H);
		}
		return false;
	}


	bool UpwardPlanarity::upwardPlanarEmbed_triconnected(Graph &G)
	{
		if (isTriconnected(G) && isAcyclic(G)) {
			BoyerMyrvold p;
			if (!p.planarEmbed(G)) return false;
			return isUpwardPlanar_embedded(G);
		}
		return false;
	}

	//
	// Single-source digraphs
	//

	bool UpwardPlanarity::isUpwardPlanar_singleSource(const Graph &G)
	{
		NodeArray<SListPure<adjEntry> > adjacentEdges;
		return UpwardPlanaritySingleSource::testAndFindEmbedding(G, false, adjacentEdges);
	}


	bool UpwardPlanarity::upwardPlanarEmbed_singleSource(Graph &G)
	{
		NodeArray<SListPure<adjEntry> > adjacentEdges(G);
		if(UpwardPlanaritySingleSource::testAndFindEmbedding(G, true, adjacentEdges) == false)
			return false;

		node superSink;
		SList<edge> augmentedEdges;
		UpwardPlanaritySingleSource::embedAndAugment(G, adjacentEdges, false, superSink, augmentedEdges);

		return true;
	}


	bool UpwardPlanarity::upwardPlanarAugment_singleSource(Graph &G)
	{
		node superSink;
		SList<edge> augmentedEdges;

		return upwardPlanarAugment_singleSource(G, superSink, augmentedEdges);
	}


	bool UpwardPlanarity::upwardPlanarAugment_singleSource(
		Graph &G,
		node &superSink,
		SList<edge> &augmentedEdges)
	{
		NodeArray<SListPure<adjEntry> > adjacentEdges(G);
		if(UpwardPlanaritySingleSource::testAndFindEmbedding(G, true, adjacentEdges) == false)
			return false;

		UpwardPlanaritySingleSource::embedAndAugment(G, adjacentEdges, true, superSink, augmentedEdges);
		return true;
	}


	bool UpwardPlanarity::isUpwardPlanar_singleSource_embedded(
		const ConstCombinatorialEmbedding &E,
		SList<face> &externalFaces)
	{
		const Graph &G = E;
		OGDF_ASSERT(G.representsCombEmbedding());

		externalFaces.clear();

		// trivial cases
		if(G.empty())
			return true;

		if(isAcyclic(G) == false)
			return false;

		// determine the single source in G
		node s;
		if(!hasSingleSource(G,s))
			return false;

		// construct face-sink graph anf find possible external faces
		FaceSinkGraph F(E,s);
		F.possibleExternalFaces(externalFaces);

		return !externalFaces.empty();
	}


	bool UpwardPlanarity::upwardPlanarAugment_singleSource_embedded(
		Graph &G,
		node  &superSink,
		SList<edge> &augmentedEdges)
	{
		OGDF_ASSERT(G.representsCombEmbedding());

		// trivial cases
		if(G.empty())
			return true;

		if(isAcyclic(G) == false)
			return false;

		// determine the single source in G
		node s;
		if(!hasSingleSource(G,s))
			return false;

		// construct embedding represented by G and face-sink graph
		ConstCombinatorialEmbedding E(G);
		FaceSinkGraph F(E,s);

		// find possible external faces
		SList<face> externalFaces;
		F.possibleExternalFaces(externalFaces);

		if (externalFaces.empty())
			return false;

		else {
			F.stAugmentation(F.faceNodeOf(externalFaces.front()), G, superSink, augmentedEdges);
			return true;
		}
	}


 } // end namespace ogdf
