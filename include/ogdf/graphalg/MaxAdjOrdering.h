/** \file
 * \brief Calculate one or all Maximum Adjacency Ordering(s) of a given simple undirected graph.
 *
 * \author Sebastian Semper
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


#ifndef OGDF_MAXADJORDERING_H
#define OGDF_MAXADJORDERING_H

#include <ogdf/basic/basic.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/basic/graphics.h>
#include <ogdf/misclayout/LinearLayout.h>

#include <vector>
#include <cmath>
#include <math.h>


namespace ogdf {
	/**
		* \brief Calculate one or all Maximum Adjacency Ordering(s) of a given simple undirected graph.
		*
		*
		* The class \a MaxAdjOrdering provides one algorithm to calculate a MAO or all MAOs of a given graph.
		* It returns a ListPure of NodeElements or a list of ListPures that contain the ordering.
		*
		*/
	class OGDF_EXPORT MaxAdjOrdering{
	private:
		/**
		* \brief This method gets called recursively to generate all MAOs
		* via backtracking
		* @param n Number of nodes
		* @param currentOrder The partial MAO to this point
		* @param currentUnsorted Nodes that are still left to sort
		* @param r Values on the nodes that count the edges to the partial MAO
		* @param MAOs The list of all MAOs to this point
		*/
		void m_calcAllMAOs_recursion(
				int n,
				ListPure<NodeElement *> currentOrder,
				ListPure<NodeElement *> currentUnsorted,
				std::vector<int> r,
				ListPure< ListPure<NodeElement *> > *MAOs
				);

		/**
		* \brief This method gets called recursively to generate all MAOs and their
		* induced forest decompositions of the edges
		* via backtracking
		* @param n Number of nodes
		* @param currentOrder The partial MAO to this point
		* @param currentForest The partial forest decomposition to this point
		* @param currentUnsorted Nodes that are still left to sort
		* @param r Values on the nodes that count the edges to the partial MAO
		* @param MAOs The list of all MAOs to this point
		* @param Fs The list of all forests to this point
		*/
		void m_calcAllMAOs_recursion(
				int n,
				ListPure<NodeElement *> currentOrder,
				ListPure< ListPure <EdgeElement *> > currentForest,
				ListPure<NodeElement *> currentUnsorted,
				std::vector<int> r,
				ListPure< ListPure<NodeElement *> > *MAOs,
				ListPure< ListPure< ListPure <EdgeElement *> > > *Fs
				);
	public:
		/**
		* \brief Standard constructor
		*/
		MaxAdjOrdering();
		/**
		* \brief Standard destructor
		*/
		~MaxAdjOrdering();

		/**
		* \brief Calculates one MAO starting with the node with index 0
		* @param G is the Graph to work on
		* @param MAO is the pointer to a list that stores the MAO
		*/
		void calc(
				const Graph *G,
				ListPure<NodeElement *> *MAO
				);
                /**
                * \brief Calculates one MAO starting with the node with index 0
                * and lex-bfs tie breaking.
                * @param G is the Graph to work on
                * @param MAO is the pointer to a list that stores the MAO
                */
                void calcBfs(
                                const Graph *G,
                                ListPure<NodeElement *> *MAO
                                );

		/**
		* \brief Calculates one MAO starting with a given node
		* @param G is the Graph to work on
		* @param s Node to start the MAO with
		* @param MAO is the pointer to a list that stores the MAO
		*/
		void calc(
				const Graph *G,
				NodeElement *s,
				ListPure<NodeElement *> *MAO
				);

		/**
		* \brief Calculates one MAO starting with the node with index 0
		* @param G is the Graph to work on
		* @param MAO is the pointer to a list that stores the MAO
		* @param Forests is the pointer to a list that stores the forest decomposition associated with the MAO
		*/
		void calc(
				const Graph *G,
				ListPure<NodeElement *> *MAO,
				ListPure< ListPure<EdgeElement *> > *Forests
				);

		/**
		* \brief Calculates one MAO starting with a given node
		* @param G is the Graph to work on
		* @param s Node to start the MAO with
		* @param MAO is the pointer to a list that stores the MAO
		* @param Forests is the pointer to a list that stores the forest decomposition associated with the MAO
		*/
		void calc(
				const Graph *G,
				NodeElement *s,
				ListPure<NodeElement *> *MAO,
				ListPure< ListPure<EdgeElement *> > *Forests
				);

		/**
		* \brief Calculates all MAOs of a given graph
		* @param G is the Graph to work on
		* @param MAOs List of all MAOs. So it must be a list of lists
		*/
		void calcAll(
				const Graph *G,
				ListPure< ListPure<NodeElement *> > *MAOs
				);

		/**
		* \brief Calculates all MAOs including their associated forest decompositions of a given graph
		* @param G is the graph to work on
		* @param MAOs List of all MAOs. So it must be a list of lists of vertices.
		* @param Fs List of all forest decompositions. So it must be a list of lists of lists of edges.
		*/
		void calcAll(
				const Graph *G,
				ListPure< ListPure< NodeElement *> > *MAOs,
				ListPure< ListPure< ListPure <EdgeElement *> > > *Fs
				);
		/**
		 * \brief Test if a given ordering is a MAO
		 * @param G is the graph to work on
		 * @param Ordering is a \a ListPure that contains a permutation of the nodes
		*/
		bool testIfMAO(
				const Graph *G,
				ListPure< NodeElement *> *Ordering
				);
        /**
         * \brief Test if a given ordering is a MAO that follows lex-bfs tie breaking
         * @param G is the graph to work on
         * @param Ordering is a \a ListPure that contains a permutation of the nodes
        */
        bool testIfMAOBfs(
                const Graph *G,
                ListPure< NodeElement *> *Ordering
                );


		/**
		 * @brief testIfAllMAOs checks all permutations (must be provided) if they are a MAO and if
		 * yes searches this ordering in the provided list. If one permutation is no MAO
                 * it still gets searched to rule out any false elements in \a Perms. So we make sure we generated
		 * all MAOs of \a G.
		 * @param G is the graph to work on
		 * @param Orderings contains Lists that are supposedly MAOs
		 * @param Perms contains all permutations of a graph with the same size as \a G
		 * @param verbose determines whether to generate output while testing
		 * @return
		 */
		bool testIfAllMAOs(
				const Graph *G,
				ListPure< ListPure< NodeElement *> > *Orderings,
				ListPure< ListPure< NodeElement *> > *Perms,
				bool verbose
				);
		/**
		* \brief Convenient way to visualize a MAO with the \a LinearLayout class.
		* @param GA Graphattributes of the Graph to paint
		* @param MAO Mao to use for ordering the nodes
		*/
		void visualize(
				GraphAttributes *GA,
				ListPure<NodeElement *> *MAO
				);

		/**
		* \brief Convenient way to visualize a MAO with the \a LinearLayout class.
		* @param GA Graphattributes of the Graph to paint
		* @param MAO Mao to use for ordering the nodes
		* @param F A Forest can also be included
		*/
		void visualize(
				GraphAttributes *GA,
				ListPure<NodeElement *> *MAO,
				ListPure< ListPure <EdgeElement *> > *F
				);

		};

	} // end namespace ogdf


#endif
