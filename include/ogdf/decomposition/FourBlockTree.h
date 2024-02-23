/** \file
 * \brief Declaration of FourBlockTree.
 *
 * \author Gregor Diatzko
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

#include <ogdf/basic/AdjEntryArray.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/NodeArray.h>

#include <memory>
#include <vector>

namespace ogdf {

/**
 * A node in a 4-block tree.
 *
 * Since each node contains its children, the root is the entire tree.
 */
struct FourBlockTree {
    /**
     * The 4-connected component.
     */
    Graph g;

    /**
     * The nodes in the original graph corresponding to the nodes in g.
     *
     * Since nodes may appear in multiple 4-connected components, these
     * need not be unique across nodes of the 4-block tree.
     */
    NodeArray<node> originalNodes;

    /**
     * A half-edge in g such that the external face of g is to its right.
     */
    adjEntry externalFace;

    /**
     * The parent node of this node in the 4-block tree.
     *
     * If this node is the root node, parent is nullptr.
     */
    FourBlockTree* parent;

    /**
     * The half-edge in parent->g corresponding to externalFace.
     *
     * If this node is the root node, parentFace is nullptr.
     */
    adjEntry parentFace;

    /**
     * The child nodes of this nodes.
     */
    std::vector<std::unique_ptr<FourBlockTree>> children;

    /**
     * Construct a 4-block tree of the given graph.
     *
     * @param g The plane triangulated graph whose 4-block tree shall be constructed.
     *          This graph will be used destructively.
     *          Edge directions in g are not respected.
     *          The order of edges at each node is used as the combinatorial
     *          embedding.
     * @param externalFace A half-edge in g such that the external face of g
     *                     lies to its right.
     */
    static FourBlockTree construct(const Graph& g, adjEntry externalFace);
};

} // namespace ogdf
