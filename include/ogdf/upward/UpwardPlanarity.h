/*
 * $Revision: 3475 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-05-02 10:43:51 +0200 (Do, 02. Mai 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class UpwardPlanarity, which implements different types
 * 			of algorithms testing upward planarity of graphs with different restrictions.
 * 			Actually, restrictions are:
 * 			- fixed embedding
 *          - single source
 *          - triconnected
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_UPWARDPLANARITY_H
#define OGDF_UPWARDPLANARITY_H


#include <ogdf/basic/Graph.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/CombinatorialEmbedding.h>


namespace ogdf {


//! Upward planarity testing and embedding.
/**
 * This class provides various static functions for upward planarity testing
 * and embedding. These functions perform different tasks (testing, embedding,
 * augmenting) and pose different restrictions on the input graph (general,
 * biconnected, single source).
 *
 * We use a strict naming scheme to make it clear what the functions are doing
 * and which restrictions they have. The prefix of the function name denotes
 * the <em>task</em>:
 *   - <tt>isUpwardPlanar</tt>: Tests if the input graph is upward planar.
 *   - <tt>upwardPlanarEmbed</tt>: First tests if the input graph is upward planar, and if yes
 *     upward planarly embeds it.
 *   - <tt>upwardPlanarAugment</tt>: Adds additional edges to the input graph such that the graph
 *     remains upward planar and fulfills a special property like single source.
 *
 * The suffix of the function name (if present) describes <em>additional restrictions</em>:
 *   - <tt>singleSource</tt>: The input graph has exactly one source (but possibly several sinks).
 *   - <tt>seriesParallel</tt>: The input graph is a series-parallel graph.
 *
 * Some of the functions take a combinatorial embedding (e.g., given by the order in the adjacency
 * lists) as input and test this given embedding. These functions are appended by <tt>_embedded</tt>.
 */
class OGDF_EXPORT UpwardPlanarity
{
public:
	/**
	 * @name General digraphs
	 */
	//@{

	//! Tests whether graph \a G is upward planarly embedded.
	/**
	 * The fixed embedding of \a G is given by the order of \a G's adjacency lists.
	 *
	 * \param G is the input graph to be tested.
	 * \return true if \a G is upward planarly embedded, false otherwise.
	 */
	static bool isUpwardPlanar_embedded(const Graph &G);

	//! Tests whether graph \a G is upward planarly embedded and computes the set of possible external faces.
	static bool isUpwardPlanar_embedded(const Graph &G, List<adjEntry> &possibleExternalFaces);


	//@}
	/**
	 * @name Triconnected digraphs
	 */
	//@{

	//! Tests whether the triconnected digraph \a G is upward planar.
	/**
	 * \remark If \a G is not triconnected the function returns false.
	 *
	 * \param G is the (triconnected) input digraph.
	 * \return true if \a G was triconnected and upward planar, false otherwise.
	 */
	static bool isUpwardPlanar_triconnected(const Graph &G);

	//! Upward planarly embeds the triconnected digraph \a G.
	/**
	 * \remark If \a G is not triconnected the function returns false.
	 *
	 * \param G is the (triconnected) input digraph.
	 * \return true if \a G was triconnected and upward planar, false otherwise.
	 */
	static bool upwardPlanarEmbed_triconnected(Graph &G);


	//@}
	/**
	 * @name Single-source digraphs
	 */
	//@{

	//! Tests whether the single-source digraph \a G is upward planar.
	/**
	 * \remark If \a G is not single-source the function returns false.
	 *
	 * \param G is the (single-source) input digraph.
	 * \return true if \a G was single-source and upward planar, false otherwise.
	 */
	static bool isUpwardPlanar_singleSource(const Graph &G);

	//! Upward planarly embeds the single-source digraph \a G.
	/**
	 * \remark If \a G is not single-source the function returns false.
	 *
	 * \param G is the (single-source) input digraph.
	 * \return true if \a G was single-source and upward planar, false otherwise.
	 */
	static bool upwardPlanarEmbed_singleSource(Graph &G);

	//! Tests whether single-source digraph \a G is upward planar, and if yes augments it to a planar st-digraph.
	/**
	 * \remark If \a G is not single-source the function returns false.
	 *
	 * If \a G is upward planar, this method adds a super sink node \a t and adds further edges such that the
	 * resulting digraph is a planar st-digraph.
	 *
	 * \param G is the input digraph which gets augmented.
	 * \return true if \a G is single-source and upward planar, false otherwise.
	 */
	static bool upwardPlanarAugment_singleSource(Graph &G);

	//! Tests whether single-source digraph \a G is upward planar, and if yes augments it to a planar st-digraph.
	/**
	 * \remark If \a G is not single-source the function returns false.
	 *
	 * If \a G is upward planar, this method adds a super sink node \a t and adds further edges such that the
	 * resulting digraph is a planar st-digraph.
	 *
	 * \param G              is the input digraph which gets augmented.
	 * \param superSink      is assigned the inserted super sink node.
	 * \param augmentedEdges is assigned the list of inserted edges.
	 * \return true if \a G is single-source and upward planar, false otherwise.
	 */
	static bool upwardPlanarAugment_singleSource(
		Graph &G,
		node &superSink,
		SList<edge> &augmentedEdges);


	//! Tests whether the embedding \a E of a single-source digraph is upward planar.
	/**
	 * \param E             is the given combinatorial embedding to be tested.
	 * \param externalFaces is assigned the list of possible external faces such that \a E is upward planar.
	 * \return true if \a E is upward planar, false otherwise.
	 */
	static bool isUpwardPlanar_singleSource_embedded(
		const ConstCombinatorialEmbedding &E,
		SList<face> &externalFaces);

	//! Tests if \a single-source digraph G is upward planarly embedded and augments it to a planar st-digraph.
	/**
	 * \param G              is the embedded input graph.
	 * \param superSink      is assigned the added super sink.
	 * \param augmentedEdges is assigned the list of added edges.
	 * \return true if \a G is upward planarly embedded (in this case \a G is also augmented by adding a
	 *         super sink and additional edges), false otherwise.
	 */
	static bool upwardPlanarAugment_singleSource_embedded(
		Graph &G,
		node  &superSink,
		SList<edge> &augmentedEdges);

	//@}
};


} // end namespace ogdf


#endif
