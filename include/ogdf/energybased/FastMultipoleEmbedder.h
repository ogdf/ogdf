/** \file
 * \brief Declaration of Fast-Multipole-Embedder layout algorithm.
 *
 * \author Martin Gronemann
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

#include <ogdf/basic/Graph.h>
#include <ogdf/module/LayoutModule.h>
#include <ogdf/internal/energybased/MultilevelGraph.h>
#include <ogdf/internal/energybased/FMEThread.h>
#include <ogdf/internal/energybased/FMEFunc.h>
#include <ogdf/internal/energybased/GalaxyMultilevel.h>

namespace ogdf {

//! The fast multipole embedder approach for force-directed layout.
/**
 * @ingroup gd-energy
 */
class OGDF_EXPORT FastMultipoleEmbedder : public LayoutModule
{
public:
	//! constructor
	FastMultipoleEmbedder();

	//! destructor
	~FastMultipoleEmbedder();

#if 0
	//! Calls the algorithm for graph \a MLG.
	//Does not do anything smart, can be safely removed
	void call(MultilevelGraph &MLG);
#endif

	//! Calls the algorithm for graph \a G with the given edgelength and returns the layout information in \a nodeXPosition, nodeYPosition.
	void call(
		const Graph& G,
		NodeArray<float>& nodeXPosition,
		NodeArray<float>& nodeYPosition,
		const EdgeArray<float>& edgeLength,
		const NodeArray<float>& nodeSize);

	//! Calls the algorithm for graph \a GA with the given edgelength and returns the layout information in \a GA.
	void call(GraphAttributes &GA, const EdgeArray<float>& edgeLength, const NodeArray<float>& nodeSize);

	//! Calls the algorithm for graph \a GA and returns the layout information in \a GA.
	virtual void call(GraphAttributes &GA) override;

	//! sets the maximum number of iterations
	void setNumIterations(uint32_t numIterations) { m_numIterations = numIterations; }

	//! sets the number of coefficients for the expansions. default = 4
	void setMultipolePrec(uint32_t precision) { m_precisionParameter = precision; }

	//! if true, layout algorithm will randomize the layout in the beginning
	void setRandomize(bool b) { m_randomize = b; }

	//!
	void setDefaultEdgeLength(float edgeLength) { m_defaultEdgeLength = edgeLength; }

	//!
	void setDefaultNodeSize(float nodeSize) { m_defaultNodeSize = nodeSize; }

	//!
	void setNumberOfThreads(uint32_t numThreads) {
#ifndef OGDF_MEMORY_POOL_NTS
		m_maxNumberOfThreads = numThreads;
#endif
	}

#if 0
	void setEnablePostProcessing(bool b) { m_doPostProcessing = b; }
#endif

private:
	void initOptions();

	void runMultipole();

	void runSingle();

	//! runs the simulation with the given number of iterations
	void run(uint32_t numIterations);

	//! allocates the memory
	void allocate(uint32_t numNodes, uint32_t numEdges);

	//! frees the memory
	void deallocate();

	uint32_t m_numIterations;

	ArrayGraph* m_pGraph;

	FMEThreadPool* m_threadPool;

	FMEGlobalOptions* m_pOptions;

	uint32_t m_precisionParameter;

	bool m_randomize;

	float m_defaultEdgeLength;

	float m_defaultNodeSize;

	uint32_t m_numberOfThreads;

	uint32_t m_maxNumberOfThreads;
};


//! The fast multipole multilevel embedder approach for force-directed multilevel layout.
/**
 * @ingroup gd-energy
 */
class OGDF_EXPORT FastMultipoleMultilevelEmbedder : public LayoutModule
{
public:
	//! Constructor, just sets number of maximum threads
	FastMultipoleMultilevelEmbedder() : m_iMaxNumThreads(1) {}
	//! Calls the algorithm for graph \a GA and returns the layout information in \a GA.
	void call(GraphAttributes &GA) override;

	//! sets the bound for the number of nodes for multilevel step
	void multilevelUntilNumNodesAreLess(int nodesBound) { m_multiLevelNumNodesBound = nodesBound; }

	void maxNumThreads(int numThreads) { m_iMaxNumThreads = numThreads; }
private:
	//! internal function to compute a good edgelength
	void computeAutoEdgeLength(const GraphAttributes& GA, EdgeArray<float>& edgeLength, float factor = 1.0f);

	//! internal main function for the multilevel layout
	void run(GraphAttributes& GA, const EdgeArray<float>& edgeLength);

	//! creates all multilevels
	void createMultiLevelGraphs(Graph* pGraph, GraphAttributes& GA, const EdgeArray<float>& edgeLength);

	//! init the original graphs multilevel
	void initFinestLevel(GraphAttributes &GA, const EdgeArray<float>& edgeLength);

	//! calls the fast multipole embedder on current level
	void layoutCurrentLevel();

	//! assigns the nodes in the current level coords by coarser level
	void assignPositionsFromPrevLevel();

	//! writes the current level to graph attributes. used for output
	void writeCurrentToGraphAttributes(GraphAttributes& GA);

	//! refine
	void nextLevel();

	//! initialize datastructure by current level
	void initCurrentLevel();

	//! clean up the multilevel graphs
	void deleteMultiLevelGraphs();

	//! for debugging only
	void dumpCurrentLevel(const char *filename);

	//! computes the maximum number of iterations by level nr
	uint32_t numberOfIterationsByLevelNr(uint32_t levelNr);

	int				  m_iMaxNumThreads;
	int				  m_iNumLevels;
	int				  m_multiLevelNumNodesBound;

	GalaxyMultilevel* m_pCurrentLevel;
	GalaxyMultilevel* m_pFinestLevel;
	GalaxyMultilevel* m_pCoarsestLevel;

	Graph*			  m_pCurrentGraph;
	NodeArray<float>* m_pCurrentNodeXPos;
	NodeArray<float>* m_pCurrentNodeYPos;
	EdgeArray<float>* m_pCurrentEdgeLength;
	NodeArray<float>* m_pCurrentNodeSize;
	NodeArray<float>  m_adjustedNodeSize;
	int				  m_iCurrentLevelNr;

	NodeArray<float>* m_pLastNodeXPos;
	NodeArray<float>* m_pLastNodeYPos;
};

} // end of namespace ogdf
