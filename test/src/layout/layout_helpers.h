//****************:*****************************************
//  Test helpers for layout algorithms
//
//  Author: Carsten Gutwenger, Tilo Wiedera
//*********************************************************
#ifndef LAYOUT_HELPERS_H
#define LAYOUT_HELPERS_H

#include <random>
#include <bandit/bandit.h>

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/module/LayoutModule.h>
#include <ogdf/module/GridLayoutModule.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>

#define MAX_NODES 200
// FIXME: Crashes upon setting MIN_NODES to 0
//         probably related to planarBiconnectedGraph
#define MIN_NODES 25
#define STEP_SIZE 25

using namespace ogdf;

enum GraphRequirementFlags {
	GR_ALL = 0,
	GR_PLANAR = 1,
	GR_TRIPLE_CONNECTED = 2,
	GR_CONNECTED = 4
};

extern void insertGraph(Graph &g, const Graph &g2);

extern void createDisconnectedGraph(
	Graph &G,
	int n_max,
	double density_min,
	double density_max,
	int cc);

extern void createAlmostPlanarGraph(Graph &G, int n, int m, int add_m);

extern void getRandomLayout(GraphAttributes &GA);

extern int64_t callLayout(const Graph &G, LayoutModule &L, bool isGridLayout);

/**
 * Runs several tests for a given layout module.
 * The layout algorithm is executed for different graphs.
 * There are no assertions yet.
 *
 * \param name
 * 	the name to be used for describing this module
 * \param L
 * 	the module to be tested
 * \param req
 * 	the requirements for graphs to be drawn, see GraphRequirementFlags for details
 * \param maxNodes
 * 	the maximum number of nodes the algorithm should be tested on
 * \param isGridLayout
 * 	set this to true if the layout module is a grid layout module
 **/
extern void describeLayoutModule(
  const std::string name,
  LayoutModule &L,
  long extraAttributes = 0,
  int req = GR_ALL,
  int maxNodes = MAX_NODES,
  bool isGridLayout = false);

extern void describeGridLayoutModule(
  const std::string name,
  GridLayoutModule &L,
  int req = GR_ALL,
  int maxNodes = MAX_NODES);

#endif
