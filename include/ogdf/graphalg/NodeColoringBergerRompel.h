/** \file
 * \brief Applies the node coloring approximation specified by Berger&Rompel.
 *
 * \author Jan-Niklas Buckow
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

#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/NodeColoringJohnson.h>
#include <ogdf/graphalg/NodeColoringModule.h>
#include <ogdf/graphalg/NodeColoringSequential.h>
#include <ogdf/graphalg/NodeColoringSimple.h>
#include <ogdf/graphalg/NodeColoringWigderson.h>

namespace ogdf {

/**
 * Approximation algorithms for the node coloring problem in graphs.
 * This class implements the approximation given by Berger&Rompel which
 * colors the graph by finding independent sets.
 */
class NodeColoringBergerRompel : public NodeColoringModule {
public:
    /**
     * Declares the procedure dealing with the brute force coloring.
     */
    enum class BruteForceProcedure {
        trivialColoring,          ///< Uses trivial coloring, i.e. colors each node with a different color
        sequentialColoringSimple, ///< Uses sequential coloring with the trivial node permutation
        sequentialColoringDegree, ///< Uses sequential coloring with the nodes sorted decreasing by their degree
        johnsonColoring,          ///< Uses Johnson's coloring algorithm
        wigdersonColoring,        ///< Uses Wigderson's coloring algorithm
        pooled,                   ///< Pools all nodes from the recursion anchor and then applies sequential coloring
    };

    /**
     * The constructor.
     * Initializes the search procedure with the binary-search by default.
     * The parameter alpha is initialized with 1.0 by default.
     */
    NodeColoringBergerRompel()
      : m_simpleColoring()
      , m_sequentialColoring()
      , m_johnsonColoring()
      , m_wigdersonColoring()
      , m_searchProcedure(SearchProcedure::wigdersonSearch)
      , m_bruteForceProcedure(BruteForceProcedure::pooled)
      , m_alpha(0.4) {}

    /**
     * Sets the search procedure to find the smallest possible parameter k such as
     * the graph is k-colorable.
     * @param searchProcedure The desired search procedure
     */
    void setSearchProcedure(SearchProcedure searchProcedure) { m_searchProcedure = searchProcedure; }

    /**
     * Sets the brute force procedure.
     * @param bruteForceProcedure The desired brute force coloring procedure
     */
    void setBruteForceProcedure(BruteForceProcedure bruteForceProcedure) {
        m_bruteForceProcedure = bruteForceProcedure;
    }

    /**
     * Sets the parameter alpha which controls the accuracy of the algorithm.
     * A higher value of alpha results in a better approximation ratio.
     * Also, the running time increases with alpha.
     * @param alpha The control parameter alpha
     */
    void setAlpha(double alpha) {
        OGDF_ASSERT(alpha > 0.1);
        m_alpha = alpha;
    }

    virtual NodeColor call(const Graph& graph, NodeArray<NodeColor>& colors, NodeColor start = 0) override {
        // Copy the input graph
        GraphCopy graphMain = GraphCopy(graph);
        preprocessGraph(graphMain);
        NodeArray<NodeColor> colorsMain(graphMain);

        // Perform a search to find the smallest possible parameter k
        auto searchWrapper = SearchWrapperBergerRompel(*this, graphMain, colorsMain, m_alpha);
        int k;
        switch (m_searchProcedure) {
            case SearchProcedure::linearSearch:
                k = searchLinear(&searchWrapper, 2, graphMain.numberOfNodes());
                break;
            case SearchProcedure::binarySearch:
                k = searchBinary(&searchWrapper, 2, graphMain.numberOfNodes());
                break;

            case SearchProcedure::wigdersonSearch:
                k = searchWigderson(&searchWrapper);
                break;

            default:
                k = searchWigderson(&searchWrapper);
                break;
        }

        // Perform Berger&Rompel with the smallest possible k
        NodeColor numberColorsBegin = start;
        bool bergerRompelResult = bergerRompelParameterized(graphMain, colorsMain, start, k, m_alpha);
        OGDF_ASSERT(bergerRompelResult);
        for (node v : graph.nodes) {
            colors[v] = colorsMain[graphMain.copy(v)];
        }
        return start - numberColorsBegin;
    }

private:
    NodeColoringSimple m_simpleColoring;
    NodeColoringSequential m_sequentialColoring;
    NodeColoringJohnson m_johnsonColoring;
    NodeColoringWigderson m_wigdersonColoring;
    SearchProcedure m_searchProcedure;
    BruteForceProcedure m_bruteForceProcedure;
    double m_alpha;

    bool bergerRompelParameterized(const Graph& graph,
      NodeArray<NodeColor>& colors,
      NodeColor& color,
      int k,
      double alpha) {
        // 1. Check special easy cases
        OGDF_ASSERT(k >= 1);
        auto startColor = color;
        if (k <= 2) {
            NodeArray<bool> partitions(graph);
            bool bipartite = isBipartite(graph, partitions);
            if (!bipartite) {
                return false;
            }
            bool partitionUsedA = false;
            bool partitionUsedB = false;
            for (node v : graph.nodes) {
                if (partitions[v]) {
                    colors[v] = color;
                    partitionUsedA = true;
                } else {
                    colors[v] = color + 1;
                    partitionUsedB = true;
                }
            }
            color += partitionUsedA + partitionUsedB;
            OGDF_ASSERT(checkColoring(graph, colors));
            return true;
        }

        // 2. Check if Johnson algorithm is sufficient
        const int n = graph.numberOfNodes();
        if (k > min(std::sqrt(n), std::pow(n, alpha))) {
            color += m_johnsonColoring.call(graph, colors, color);
            return true;
        }

        // 3. Copy the input graph
        GraphCopy graphMain = GraphCopy(graph);
        preprocessGraph(graphMain);

        // 4. Prepare the main loop
        const int m = std::floor(alpha * (std::log(n) / std::log(k)));

        // 5. Perform the coloring until the number of nodes is small enough
        while (graphMain.numberOfNodes() >= k * m) {
            // Search for big independent sets
            GraphCopy gSubgraph = graphMain;

            while (gSubgraph.numberOfNodes() >= k * m) {
                Array<Array<node>> buckets;
                createBuckets(gSubgraph, k * m, buckets);
                auto finished = false;
                for (Array<node>& bucket : buckets) {
                    for (NodeColoringBergerRompel::SubsetIterator subsetIterator(bucket, m); subsetIterator.isOk();
                         subsetIterator.advance()) {
                        auto subset = subsetIterator.currentSubset();
                        List<node> neighbors;
                        getNeighbors<ListIterator<node>>(gSubgraph, subset.begin(), neighbors);
                        auto numberNodesSubgraph = gSubgraph.numberOfNodes();
                        if (checkIndependentSet<ListIterator<node>>(gSubgraph, subset.begin()) &&
                            neighbors.size() <= (numberNodesSubgraph - numberNodesSubgraph / k)) {

                            // Color the nodes in the subset
                            for (node v : subset) {
                                colors[gSubgraph.original(v)] = color;
                            }

                            // Delete nodes in the subgraph
                            List<node> nodesToDelete;
                            mergeNodeLists<ListIterator<node>>(
                              gSubgraph, subset.begin(), neighbors.begin(), nodesToDelete);
                            for (node v : nodesToDelete) {
                                gSubgraph.delNode(v);
                            }

                            // Delete the already colored nodes from the original graph
                            for (node v : subset) {
                                graphMain.delNode(graphMain.copy(gSubgraph.original(v)));
                            }

                            finished = true;
                            break;
                        }
                    }
                    if (finished) {
                        break;
                    }
                }
                if (!finished) {
                    return false;
                }
            }
            // Increment the color the next independent set
            color++;
        }

        // 6. Brute force coloring stage
        NodeArray<NodeColor> colorsCopy(graphMain);
        // Case distinction
        switch (m_bruteForceProcedure) {
            case BruteForceProcedure::trivialColoring:
                color += m_simpleColoring.call(graphMain, colorsCopy, color);
                break;
            case BruteForceProcedure::sequentialColoringSimple:
                color += m_sequentialColoring.colorByIndex(graphMain, colorsCopy, color);
                break;
            case BruteForceProcedure::sequentialColoringDegree:
                color += m_sequentialColoring.colorByDegree(graphMain, colorsCopy, color);
                break;
            case BruteForceProcedure::johnsonColoring:
                color += m_johnsonColoring.call(graphMain, colorsCopy, color);
                break;
            case BruteForceProcedure::wigdersonColoring:
                color += m_wigdersonColoring.call(graphMain, colorsCopy, color);
                break;
            case BruteForceProcedure::pooled:
                List<node> nodesToBeColored;
                for (node v : graphMain.nodes) {
                    nodesToBeColored.emplaceBack(graphMain.original(v));
                }
                m_sequentialColoring.sortByDegree(nodesToBeColored);
                m_sequentialColoring.fromPermutation(graph, colors, nodesToBeColored, startColor, true);
                auto maxColorIndex = getMaximumNodeColor(colors);
                color = maxColorIndex + 1 - startColor;
                break;
        }
        // Copy the coloring result to the original graph
        if (m_bruteForceProcedure != BruteForceProcedure::pooled) {
            for (node v : graphMain.nodes) {
                colors[graphMain.original(v)] = colorsCopy[v];
            }
        }

        // Check the coloring
        OGDF_ASSERT(checkColoring(graph, colors));
        return true;
    }

    /**
     * Wraps the parameterized Berger&Rompel algorithm
     */
    struct SearchWrapperBergerRompel : public SearchWrapper {

        /**
         * Creates the wrapper.
         * @param coloringBergerRompel Reference to the NodeColoringBergerRompel
         * @param graph The graph to color
         * @param colors The array of colors to be assigned
         * @param alpha The alpha control parameter of the Berger&Rompel algorithm
         */
        SearchWrapperBergerRompel(NodeColoringBergerRompel& coloringBergerRompel,
          const Graph& graph,
          NodeArray<NodeColor>& colors,
          double alpha)
          : m_coloring(coloringBergerRompel)
          , m_graph(graph)
          , m_colors(colors)
          , m_alpha(alpha) {}

        bool step(int k) override {
            NodeColor start = NodeColor(1);
            return m_coloring.bergerRompelParameterized(m_graph, m_colors, start, k, m_alpha);
        }

        NodeColoringBergerRompel& m_coloring;
        const Graph& m_graph;
        NodeArray<NodeColor>& m_colors;
        double m_alpha;
    };
};
}
