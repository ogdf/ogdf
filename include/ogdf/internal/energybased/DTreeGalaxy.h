/** \file
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

#ifndef OGDF_DTREE_MULTILEVEL_H__
#define OGDF_DTREE_MULTILEVEL_H__
#include <ogdf/basic/Graph.h>

namespace ogdf {

// simple implementation of the slightly modified version of hachul by gronemann
class GalaxyLevel
{
public:
    // constructor for the finest level i.e. the original graph
    // the graph has to be connected!
    GalaxyLevel(const Graph& graph);

    // destructor. Note this deletes this level and all subsequent i.e coarser ones
    ~GalaxyLevel();

    // returns the graph
    const Graph& graph() const;

    // returns the parent node of a node on the coarser level
    node parent(node v) const;

    // returns the weight of a node
    double weight(node v) const;

    // returns the edge weight of e
    double edgeWeight(edge e) const;

    // returns the weight of a node
    void setWeight(node v, double weight);

    // returns the edge weight of e
    void setEdgeWeight(edge e, double weight);

    // returns true if this is the level of the original graph
    bool isFinestLevel() const;

    // returns true if this is the coarsest level
    bool isCoarsestLevel() const;

    // return the next coarser one
    GalaxyLevel* nextCoarser();

    // return the next finer one
    GalaxyLevel* nextFiner();

    // builds all levels until the graph has less then maxNumNodes
    // in case there are already coarser levels attached to this it will
    // extend the coarsest if necessary. if not i.e. the coarsest has < maxNumNodes
    // it will return the coarsest
    GalaxyLevel* buildLevelsUntil(int maxNumNodes);

private:
    // creates a new coarser version of this graph
    GalaxyLevel* buildNextCoarserLevel(int numLabels = 3);

    // private constructor for creating a coarser level
    GalaxyLevel(GalaxyLevel* pNextFiner);

    // remove par edges with weight
    void removeParEdgesWithWeight();

    // pointer to the next finer level
    GalaxyLevel* m_pNextFiner;

    // pointer to the next coarser
    GalaxyLevel* m_pNextCoarser;

    // the graph
    Graph* m_pGraph;

    // the weight of the node is the sum of weights of the children
    NodeArray<double> m_nodeWeight;

    // pointer to the parent node on the coarser level
    NodeArray<node> m_parent;

    // edge weight
    EdgeArray<double> m_edgeWeight;
};

}

#endif
