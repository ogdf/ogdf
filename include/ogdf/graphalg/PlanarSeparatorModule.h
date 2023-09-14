/** \file
 * \brief Declaration of base class of all planar separator algorithms.
 *
 * \author Thomas Klein
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
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/simple_graph_alg.h>

#include <map>
#include <memory>
#include <set>

namespace ogdf {
namespace planar_separators {

/**
 * Abstract description of a Breadth First Search tree.
 */
class OGDF_EXPORT BFSTree {
public:
	virtual ~BFSTree() = default;

	/**
	 * Allows access to a copy of the graph.
	 * @return the GraphCopy
	 */
	virtual GraphCopy* getGraph() const = 0;

	/**
	 * Gets the current root node of the tree.
	 * @return the root node
	 */
	virtual node getRoot() const = 0;

	/**
	 * Gets the number of nodes of the graph.
	 * @return the graph size
	 */
	virtual int getGraphSize() const = 0;

	/**
	 * Checks if an edge is a tree-edge.
	 *
	 * @param e the edge to be checked
	 * @return if the edge is in the tree
	 */
	virtual bool isInTree(edge e) const = 0;

	/**
	 * Returns the level (=depth in the tree) for a node.
	 *
	 * @param n the node
	 * @return the depth of \p n in the tree
	 */
	virtual int getLevelOfNode(node n) const = 0;

	/**
	 * Returns the node that is the parent of \p n in the tree.
	 *
	 * @param n the child node
	 * @return the parent node
	 */
	virtual node getParentOfNode(node n) const = 0;

	/**
	 * Returns the total number of children, grandchildren etc. in the subtree rooted at \p n.
	 *
	 * @param n the node
	 * @return the total number of children
	 */
	virtual int getDescendantsOfNode(node n) const = 0;

	/**
	 * Returns all (immediate) children of a node.
	 *
	 * @param n the node
	 * @return a list of all children of \p n
	 */
	virtual List<node> getChildrenOfNode(node n) const = 0;

	/**
	 * Returns the adjEntry that leads up to the parent of \p n.
	 *
	 * @param n the node
	 * @return the adjEntry to parent
	 */
	virtual adjEntry getAdjToParent(node n) const = 0;
};

/**
 * Abstract BFSTree that is realized via NodeArrays.
 */
class OGDF_EXPORT ArrayBFSTree : public BFSTree {
public:
	/**
	 * Constructor.
	 *
	 * @param G the graph spanned by the tree
	 * @param rootNode the root node for the tree
	 */
	ArrayBFSTree(GraphCopy& G, node rootNode) : pGraph {&G}, root {rootNode} { init(); }

	/**
	 * Initializes all internal arrays.
	 */
	void init() {
		// init arrays
		levelOfNode.init(*pGraph, -1);
		parentOfNode.init(*pGraph, nullptr);
		childrenOfNode.init(*pGraph);
		edgeToParent.init(*pGraph, nullptr);
		descendantsOfNode.init(*pGraph,
				1); // number of descendants of a node, the node itself counts as its own descendant
		inTree.init(*pGraph, false);
		mark.init(*pGraph, false); // records if that node was visited by BFS

		// mark root node
		mark[root] = true;
	}

	GraphCopy* getGraph() const override { return pGraph; }

	node getRoot() const override { return root; }

	int getGraphSize() const override {
		OGDF_ASSERT(pGraph != nullptr);
		return pGraph->numberOfNodes();
	}

	bool isInTree(edge e) const override { return inTree[e]; }

	int getLevelOfNode(node n) const override { return levelOfNode[n]; }

	node getParentOfNode(node n) const override { return parentOfNode[n]; }

	int getDescendantsOfNode(node n) const override { return descendantsOfNode[n]; }

	List<node> getChildrenOfNode(node n) const override { return childrenOfNode[n]; }

	adjEntry getAdjToParent(node n) const override { return edgeToParent[n]; }

#ifdef OGDF_HEAVY_DEBUG
	// this is very expensive
	void assertTreeConsistency() const {
		NodeArray<int> marked;
		marked.init(*pGraph);
		for (node no : pGraph->nodes) {
			marked[no] = no->index();
		}

		for (edge e : pGraph->edges) {
			if (isInTree(e)) {
				node s = e->source();
				node t = e->target();
				OGDF_ASSERT(marked[s] != marked[t]);
				int newIdx = marked[s] < marked[t] ? marked[s] : marked[t];
				int otherIdx = marked[s] < marked[t] ? marked[t] : marked[s];

				for (node no : pGraph->nodes) {
					if (marked[no] == otherIdx) {
						marked[no] = newIdx;
					}
				}
			}
		}

		// afterwards, assert that all nodes are in fact marked
		int lowest = marked[pGraph->firstNode()];
		for (edge e : pGraph->edges) {
			if (!isInTree(e)) {
				OGDF_ASSERT(marked[e->source()] == lowest && marked[e->target()] == lowest);
			}
		}
	}
#endif


protected:
	GraphCopy* pGraph;
	node root; // root node

	NodeArray<int> levelOfNode; // holds for each node on which level it is
	NodeArray<node> parentOfNode; // holds for each node which node is his parent in the tree
	NodeArray<List<node>> childrenOfNode; // holds for each node a list of children in the tree
	NodeArray<adjEntry> edgeToParent; // holds for each node the edge which connects it to its parent
	NodeArray<int> descendantsOfNode; // holds for each node how many descendants it has in the tree
	NodeArray<bool> mark;
	EdgeArray<bool> inTree; // hold for each edge whether it is in the BFS tree
};

/**
 * BFS tree used by both classical algorithms (LT and Dual).
 */
class OGDF_EXPORT BFSTreeClassical : public ArrayBFSTree {
public:
	/**
	 * Constructor.
	 *
	 * @param G the graph
	 * @param rootNode node at which tree should be rooted
	 * @param heightMaxIterations how many iterations of tree height maximization to run
	 * @param findLevelsSimple whether the levels should be found using the simple method or not
	 */
	BFSTreeClassical(GraphCopy& G, node rootNode, unsigned int heightMaxIterations,
			bool findLevelsSimple = false);

	~BFSTreeClassical() { }

	/**
	 * Constructs the tree.
	 *
	 * @param rootNode the root node at which construction starts in this iteration
	 * @param numIterations how many iterations of tree height maximization to run
	 */
	void construct(node rootNode, unsigned int numIterations);

	/**
	 * Reconstructs the tree using triangulating bfs.
	 */
	void reconstruct();

	/**
	 * Creates a new root node for the graph, replacing all levels below t0.
	 *
	 * @param useTriBFS whether to use triangulating BFS or not
	 */
	void createNewRoot(bool useTriBFS = false);

	/**
	 * Gets the size of a specific level.
	 *
	 * @param level index of the desired level
	 * @return the number of nodes in that level
	 */
	int getSizeOfLevel(int level) const;

	/**
	 * Returns a level of the tree.
	 *
	 * @param level index of the desired level
	 * @return List of nodes in that level
	 */
	List<node> getLevel(int level) const;

	/**
	 * Returns all nodes between two levels.
	 *
	 * @param start index of the start level
	 * @param end index of the end level
	 * @return a list containing all nodes between levels start and end
	 */
	List<node> getNodesFromTo(int start, int end) const;

	/**
	 * Returns all nodes from a given level onwards.
	 *
	 * @param start index of the start level
	 * @return all nodes from this level on
	 */
	List<node> getNodesFrom(int start) const;

	/**
	 * Restructures the tree by adding a new root and deleting all nodes below t0 and above t2, adding all nodes
	 * that are removed to the second component and keeping levels t0 and t2 in the separator.
	 *
	 * @param separator the list of separator nodes
	 * @param second the second component of the separation
	 * @param useTriBFS whether to use triangulating BFS
	 */
	void restructure(List<node>& separator, List<node>& second, bool useTriBFS = false);

	/**
	 * Removes the two separator levels t0 and t2 from the tree.
	 *
	 * @param separator the list of separator nodes
	 * @param second the second component
	 */
	void removeSeparatorLevels(List<node>& separator, List<node>& second);

	/**
	 * Checks whether a node is visited by BFS.
	 *
	 * @param n the node
	 * @return true if this node was visited
	 */
	bool isVisited(node n) const { return mark[n]; }

	int getSeparatorLevel() const { return currentLevel; }

	int get_t0() const { return t0; }

	int get_t1() const { return t1; }

	int get_t2() const { return t2; }

protected:
	/**
	 * Finds the levels t0 and t2 of the tree that might serve as separators.
	 * This is the original version described by Lipton & Tarjan.
	 */
	void findLevels();

	/**
	 * Simplified version of findLevels that simply finds a level smaller than sqrt(n).
	 */
	void findLevelsSimple();


private:
	List<List<node>> levels; // contains all levels, each as a list of nodes

	unsigned int heightMaxIterations;
	bool simple; // stupid workaround for being able to change generation of level t0 and t2

	// construction
	int currentLevel = 0;
	bool belowMiddle = true;
	double m_ratio;

	// separator levels
	int t0;
	int t1;
	int t2;
	int k; // number of nodes in levels 0 through t1

	void visit(node v, node parent, adjEntry adj, SListPure<node>& bfs);
};

/**
 * Abstract description of postprocessors.
 */
class OGDF_EXPORT Postprocessor {
public:
	Postprocessor() {};

	virtual ~Postprocessor() {};

	/**
	 * Applies the postprocessor to a given separation.
	 *
	 * @param G the graph
	 * @param separator list of separator nodes
	 * @param first first half of nodes
	 * @param second second half of nodes
	 * @return true if everything worked
	 */
	virtual bool apply(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second) = 0;

	/**
	 * Returns the human-readable identifier of this postprocessor.
	 * @return the name
	 */
	virtual std::string getName() const = 0;
};

/**
 * NodeExpulsor: Remove all nodes that are not connected to both halves of the separation.
 */
class OGDF_EXPORT NodeExpulsor : public Postprocessor {
public:
	/**
	 * Constructor. NE can either aggressively reduce separator size, even at the expense of violating the
	 * maximal size of 2/3 for one component, or only reduce separator size until the ratio is met exactly.
	 *
	 * @param balance whether to keep the components balanced
	 */
	NodeExpulsor(bool balance = true) : keepBalance {balance} { }

	bool apply(const Graph& G, List<node>& separator, List<node>& first, List<node>& second) override;

	std::string getName() const override { return "NE"; }

private:
	bool keepBalance;
};

/**
 * Dulmage-Mendelsohn-Decomposition. Improves a separation by creating a bipartite graph out of the separator and
 * the larger of the two components, finding a maximal matching in this graph, and combining certain subsets of the
 * two components to form a new separator in a way that is provably optimal.
 */
class OGDF_EXPORT DMDecomposer : public Postprocessor {
public:
	bool apply(const Graph& G, List<node>& separator, List<node>& first, List<node>& second) override;

	std::string getName() const override { return "DMD"; }

private:
	NodeArray<short> assignments; // assigns numbers to nodes: 0 = n is in separator, 1 = n is in shorter list, 2 = n is in longer list

	List<node> bipartSep; // clones of separator nodes in the bipartite graph
	List<node> bipartB; // clones of B-nodes in the bipartite graph

	// arrays from original graph G
	NodeArray<bool> isInS; // contains whether a node is in the separator or not
	NodeArray<node> clone; // maps from G to bipartite graph

	// arrays from bipartite graph
	NodeArray<node> unclone; // maps from bipartite graph to G
	EdgeArray<bool> flow; // contains whether an edge is in the matching or not

	// Dulmage Mendelsohn Decomposition
	List<node> SI;
	List<node> SX;
	List<node> SR;

	List<node> BI;
	List<node> BX;
	List<node> BR;

	/**
	 * Resets the component so it can be reused.
	 */
	void reset();

	/**
	 * Fills a the NodeArray \p assignments with the values 0, 1 and 2 to represent the assignments of nodes to
	 * the separator / first / second component, respectively.
	 *
	 * @param G the graph
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 */
	void setupAssignments(const Graph& G, const List<node>& separator, const List<node>& first,
			const List<node>& second);

	/**
	 * Creates a bipartite graph from the nodes in the separator and those in the bigger component.
	 *
	 * @param graph the graph
	 * @param separator the nodes in the separator
	 */
	void bipartiteGraph(Graph& graph, const List<node>& separator);

	/**
	 * Translates the assignments stored in the NodeArray assignments back to the lists, which can then be returned.
	 *
	 * @param G the graph
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 */
	void translateAssignments(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second) const;

	/**
	 * Calculates the subset SR and BR once SI / SX and BI / BX have been found.
	 *
	 * @param graph the graph
	 */
	void calculateRemainders(const Graph& graph);

	/**
	 * Given the subsets SI / SX / SR and BI / BX / BR, creates a separator with minimal size so that the two
	 * components are as balanced as possible.
	 *
	 * @param G the graph
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 */
	void chooseBalancedDecomposition(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second);

	/**
	 * Performs an alternating breadth first search to find all nodes in the separator that are reachable, and all
	 * those in the larger component B that are reachable.
	 *
	 * @param G the graph
	 * @param startNodes starting points for the BFS
	 * @param reachableSep nodes in the separator that could be reached
	 * @param reachableB nodes in component B that could be reached
	 * @return true if everything worked
	 */
	bool alternatingBFS(const Graph& G, const List<node>& startNodes, List<node>& reachableSep,
			List<node>& reachableB);

	/**
	 * Given all nodes in the bipartite graph, selects starting points for the BFS and finds the reachable nodes in
	 * the separator and in B.
	 *
	 * @param graph the graph
	 * @param bipart the bipartite graph
	 * @param reachableSep the reachable nodes from the separator
	 * @param reachableB the reachable nodes in the bigger component B
	 */
	void decompose(const Graph& graph, const List<node>& bipart, List<node>& reachableSep,
			List<node>& reachableB);
};

/**
 * Auxiliary data structure to represent Cycles in planar graphs.
 */
class OGDF_EXPORT Cycle {
public:
	/**
	 * Constructor. Cycles are created from a tree and a non-tree edge.
	 *
	 * @param tree the tree
	 * @param startEdge the non-tree start edge
	 */
	Cycle(BFSTree* tree, edge startEdge);

	Cycle& operator=(Cycle&& other) {
		tree = other.tree;
		nodes = std::move(other.nodes);
		edges = std::move(other.edges);
		isOnCycle = std::move(other.isOnCycle);
		isEdgeOnCycle = std::move(other.isEdgeOnCycle);
		cycleRoot = other.cycleRoot;
		costClock = other.costClock;
		costCounter = other.costCounter;
		isClockwise = other.isClockwise;

		other.tree = nullptr;
		other.cycleRoot = nullptr;

		return *this;
	}

	Cycle(Cycle&& other)
		: tree {other.tree}
		, nodes(std::move(other.nodes))
		, edges(std::move(other.edges))
		, isOnCycle(std::move(other.isOnCycle))
		, isEdgeOnCycle(std::move(other.isEdgeOnCycle))
		, cycleRoot {other.cycleRoot}
		, costClock {other.costClock}
		, costCounter {other.costCounter}
		, isClockwise {other.isClockwise} {
		other.tree = nullptr;
		other.cycleRoot = nullptr;
	}

	/**
	 * Returns the size of the cycle = the number of nodes on the cycle.
	 *
	 * @return the size
	 */
	int getSize() const { return nodes.size(); }

	/**
	 * Checks if the cycle is clockwise, i.e. you get the inward-pointing edges of each node on the cycle by
	 * rotating clockwise from cycle-edge to cycle-edge.
	 *
	 * @return true if the cycle is clockwise
	 */
	bool getClockwise() const { return isClockwise; }

	/**
	 * Gets the nodes on the cycle.
	 *
	 * @return read-only list of nodes
	 */
	const List<node>& getNodes() const { return nodes; }

	/**
	 * Gets the adjEntries on the cycle
	 *
	 * @return read-only list of adjEntries
	 */
	const List<adjEntry>& getEdges() const { return edges; }

	/**
	 * Gives access to the root node of the cycle.
	 *
	 * @return read-only access to the root node
	 */
	const node& getRoot() const { return cycleRoot; }

	/**
	 * Gets the cost on the inside of the cycle. (=the number of nodes on the inside of the cycle, where the inside
	 * is defined as the larger side.)
	 *
	 * @return the inside cost
	 */
	int getInsideCost() const;

	/**
	 * Gets the cost on the outside of the cycle. (=the number of nodes on the outside of the cycle, where the inside
	 * is defined as the larger side.)
	 *
	 * @return the inside cost
	 */
	int getOutsideCost() const;

	/**
	 * Expands the cycle, as described by Lipton and Tarjan 1979, by examining the triangle adjacent to the current
	 * expand edge on the inside of the cycle and integrating it into the cycle.
	 *
	 * @return the new, expanded cycle
	 */
	Cycle expandCycle();

	/**
	 * Utility method for printing the cycle to the console. For debugging only.
	 */
	void print() const;

	/**
	 * Fills the lists of nodes, by putting all nodes on this cycle into the list of separator nodes and putting
	 * the nodes on the sides of the cycle into first and second, respectively. (Of course only once a suitable
	 * cycle is found).
	 *
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @param useRoot whether to add the root node as well - sometimes, the root node is artificial and should be
	 * 					ignored instead of being added
	 */
	void fillLists(List<node>& separator, List<node>& first, List<node>& second,
			bool useRoot = false);

private:
	/**
	 * Special iterator to walk over the inward-pointing edges of the cycle.
	 */
	class Iterator {
	private:
		const Cycle* cycle;
		adjEntry m_current;
		const bool isClockwise;

		/**
		 * Yields the next adjEntry, given the current one (internal use only).
		 *
		 * @param current the current AdjEntry
		 * @return the next one
		 */
		adjEntry next(adjEntry current) const {
			adjEntry next;

			if (isClockwise) {
				if (cycle->isEdgeOnCycle[m_current->theEdge()]) {
					next = current->twin()->cyclicSucc();
				} else {
					next = current->cyclicSucc();
				}
			} else {
				if (cycle->isEdgeOnCycle[m_current->theEdge()]) {
					next = current->twin()->cyclicPred();
				} else {
					next = current->cyclicPred();
				}
			}
			if (next->theEdge() == cycle->getCurrentExpandEdge()->theEdge()) {
				return nullptr;
			}
			return next;
		}

	public:
		/**
		 * Constructor. Standard constructor that starts at the first clockwise adjEntry within the cycle.
		 *
		 * @param cyc the cycle to which this iterator belongs
		 * @param clockwise whether the cycle is clockwise
		 */
		Iterator(const Cycle* cyc, bool clockwise)
			: cycle {cyc}, m_current {cyc->getCurrentExpandEdge()}, isClockwise {clockwise} {
			++(*this);
		}

		/**
		 * Constructor. Constructs a dead Iterator, as returned by Cycle.end().
		 *
		 * @param cyc the cycle this iterator belongs to
		 */
		Iterator(const Cycle* cyc) : cycle(cyc), m_current(nullptr), isClockwise {false} { }

		/**
		 * Checks whether the current adjEntry is the one that leads up to the root
		 * @return true if current adjEntry leads to root
		 */
		bool isOutEdge() {
			return m_current->theNode() == cycle->cycleRoot
					&& m_current == cycle->tree->getAdjToParent(cycle->cycleRoot);
		}

		adjEntry operator*() const { return m_current; }

		Iterator& operator++() {
			m_current = next(m_current);
			while (m_current != nullptr
					&& cycle->isEdgeOnCycle[m_current->theEdge()]) { //m_current == cycle->sep.edgeToParent[m_current->theNode()]) {
				m_current = next(m_current);
			}
			return *this;
		}

		Iterator operator++(int) {
			Iterator tmp = *this;
			++(*this);
			return tmp;
		}

		friend bool operator==(const Iterator& a, const Iterator& b) {
			return a.m_current == b.m_current;
		};

		friend bool operator!=(const Iterator& a, const Iterator& b) {
			return a.m_current != b.m_current;
		};
	};

	BFSTree* tree; // the BFS-Tree with which this cycle is associated

	List<node> nodes; // nodes on the cycle
	List<adjEntry> edges; // edges on the cycle
	NodeArray<bool> isOnCycle; // holds for each node whether it is on the cycle or not
	EdgeArray<bool> isEdgeOnCycle; // holds for each edge whether it is on the cycle or not

	node cycleRoot; // root node of the cycle, i.e. node in which its two arms meet

	int costClock {0}; // cost of clockwise nodes
	int costCounter {0}; // cost of counter-clockwise nodes
	bool isClockwise;

	/**
	 * Constructor. Creates a Cycle from a complete description (used by the expansion procedure)
	 *
	 * @param tree the tree to create this cycle from
	 * @param nodeList list of nodes that should be on the cycle
	 * @param edgeList list of adjEntries that should be on the cycle
	 * @param root the root node of the cycle (= node with lowest level)
	 * @param clockwise whether the cycle is clockwise or not
	 */
	Cycle(BFSTree* tree, List<node>& nodeList, List<adjEntry>& edgeList, node root, bool clockwise);

	/**
	 * Constructor. Creates an empty cycle (used only during the expansion procedure).
	 *
	 * @param tree the tree
	 * @param clockwise whether the cycle is clockwise
	 */
	Cycle(BFSTree* tree, bool clockwise);

	/**
	 * Initializes the cycle (used by constructors).
	 *
	 * @param nodeList the list of nodes on the cycle
	 * @param edgeList the list of adjEntries on the cycle
	 * @param root the root node of the cycle
	 */
	void init(List<node>& nodeList, List<adjEntry>& edgeList, node root);

	/** Iterators for clockwise iteration over edges on the "inside", starting and ending at currentExpandEdge */
	Iterator begin() const { return Iterator(this, isClockwise); }

	Iterator end() const { return Iterator(this); }

	/**
	 * Access methods for adding and removing nodes from the cycle.
	 */
	void popBackNode();

	void popBackEdge();

	void popFrontNode();

	void popFrontEdge();

	void pushFrontEdge(adjEntry adj);

	/**
	 * Gets the current non-tree edge on the cycle, which is the one that will be used to expand the cycle further.
	 *
	 * @return the current non-tree edge
	 */
	adjEntry getCurrentExpandEdge() const;

	/**
	 * Expand the Cycle when one of the inner edges of the triangle is a tree edge.
	 * This method changes this cycle to expand it.
	 *
	 * @param y the tip of the triangle, inside of the cycle
	 * @param v the starting point of the triangle-edge
	 * @param w the end point of the triangle-edge
	 * @param vy the edge that leads to the tip
	 * @param yw the edge that leads away from the tip
	 */
	void expandWithTreeEdge(node y, node v, node w, adjEntry vy, adjEntry yw);

	/**
	 * Expand the cycle if none of the inner edges of the triangle is a tree edge.
	 * This method does not change this cycle, but creates a new one (because during the process, two candidates
	 * are created and one of them is selected).
	 *
	 * @param y the tip of the triangle, inside of the cycle
	 * @param v the starting point of the triangle-edge
	 * @param w the end point of the triangle-edge
	 * @param vy the edge that leads to the tip
	 * @param yw the edge that leads away from the tip
	 * @return the next cycle
	 */
	Cycle expandWithoutTreeEdges(node y, const node v, const node w, const adjEntry vy,
			const adjEntry yw);

	/**
	 * Used during expansion without tree edges. Finds a path from y back to the cycle by following parent pointers.
	 * Returns a pair of nodes (z, r) where
	 * z = the node on the cycle in which the path beginning at y enters the cycle, and
	 * r = the node on the cycle with the shortest distance to overall root
	 *
	 * @param y the tip of the triangle
	 * @param pathNodes the nodes on the path
	 * @param pathAdjEntries the adjEntries on the path
	 * @return nodes (z, r)
	 */
	std::pair<node, node> findPathToCycle(node& y, List<node>& pathNodes,
			List<adjEntry>& pathAdjEntries) const;

	/**
	 * Used during expansion without tree edges. Finds the first cycle bounded by the path from y to z and the
	 * cycle itself.
	 *
	 * @param cyc the original cycle
	 * @param pathNodes the nodes on the path from y to z
	 * @param pathAdjEntries the adjEntries on the path from y to z
	 * @param z the node z where the path meets the cycle
	 * @param propRoot the proposed root node = the root node of the old cycle, which may be improved on
	 * @param yw the side of the triangle
	 * @param oldNodes nodes on old cycle
	 * @param oldEdges edges on old cycle
	 * @return whether we found the root node of the original cycle on our way from w to z
	 */
	bool findAlphaCycle(Cycle& cyc, const List<node>& pathNodes,
			const List<adjEntry>& pathAdjEntries, const node z, const node propRoot,
			const adjEntry yw, List<node>& oldNodes, List<adjEntry>& oldEdges) const;

	/**
	 * Used during expansion without tree edges. Finds the second cycle bounded by the path from y to z and the
	 * cycle itself.
	 *
	 * @param cyc the original cycle
	 * @param pathNodes the nodes on the path from y to z
	 * @param pathAdjEntries the adjEntries on the path from y to z
	 * @param z the node z where the path meets the cycle
	 * @param propRoot the proposed root node = the root node of the old cycle, which may be improved on
	 * @param vy the side of the triangle
	 * @param oldNodes nodes on old cycle
	 * @param oldEdges edges on old cycle
	 * @param foundRootOnAlpha whether the alpha cycle found the old root node
	 */
	void findBetaCycle(Cycle& cyc, const List<node>& pathNodes, const List<adjEntry>& pathAdjEntries,
			const node z, const node propRoot, const adjEntry vy, List<node>& oldNodes,
			List<adjEntry>& oldEdges, bool foundRootOnAlpha) const;

	/**
	 * Increases the cost on the inside of this cycle by following the given adjEntry and adding the cost of the
	 * node on the other end of the adjEntry.
	 *
	 * @param adj the adjEntry whose cost should be added
	 * @param clockwise whether the cost should be added to the clockwise cost
	 */
	void increaseCost(adjEntry adj, bool clockwise);

	/**
	 * Recursively creates a list containing all descendants of a node. Collects the nodes in the original graph,
	 * not those in the copy.
	 *
	 * @param no the node whose children should be collected
	 * @param marked holds which nodes were visited already
	 * @param list the list containing all children
	 * @param useRoot whether to treat the root as a normal node or not
	 */
	void collectChildrenOfNode(const node no, NodeArray<bool>& marked, List<node>& list,
			bool useRoot = false) const;

	/**
	 * Computes the costs on both sides of the cycle. Costs are stored in the respective variables.
	 */
	void computeCosts();
};

} // namespace planar_separators

using namespace planar_separators;

//! Abstract description of all planar separator algorithms.
/**
 * @ingroup ga-plansep
 */
class OGDF_EXPORT PlanarSeparatorModule {
public:
	PlanarSeparatorModule() { }

	virtual ~PlanarSeparatorModule() { }

	/**
	 * Adds a postprocessor to this separator, which will always be applied.
	 *
	 * @param post the postprocessor
	 */
	void addPostProcessor(Postprocessor& post) { postProcessors.push_back(&post); }

	/**
	 * Deletes all appended postprocessors from this separator.
	 */
	void clearPostProcessors() { postProcessors.clear(); }

	void setStartIndex(int index) { startNodeIndex = index; }

	/**
	 * Separates a planar graph.
	 * This method takes care of multiple components, makes sure that all preconditions are fulfilled and applies
	 * postprocessing, both general postprocessing (only useful for small graphs) and via the added postprocessors.
	 *
	 * @pre The input graph is planar, simple and undirected. If you know that all conditions hold and the graph is already
	 * 		planarly embedded, you can set checkPreconditions = false for a small speedup.
	 *
	 * @param G the graph to be separated
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @param checkPreconditions whether to ensure that the current embedding is planar and G is simple undirected
	 * @return true on success
	 */
	virtual bool separate(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second, bool checkPreconditions = true) final {
		if (setup(G, separator, first, second, checkPreconditions)) {
			return true;
		}

		reset(); // reset everything to ensure the module can be reused

		bool result = doSeparate(G, separator, first, second); // call core algorithm

		if (result) {
			cleanup(G, separator, first, second);

			postProcess(G, separator, first, second);
		}

		return result;
	}

	/**
	 * Separates a planar graph. Represents a solution as a NodeArray that assigns the component to each node.
	 * 0 = separator node, 1 = first component, 2 = second component
	 *
	 * @param G the graph to be separated
	 * @param assignments the NodeArray containing the assignments
	 * @param checkPreconditions whether to check if the current embedding is planar and G is simple undirected
	 * @return true on success
	 */
	virtual bool separate(const Graph& G, NodeArray<short>& assignments,
			bool checkPreconditions = true) final {
		OGDF_ASSERT(assignments.graphOf() == &G);

		List<node> separator;
		List<node> first;
		List<node> second;

		bool result = separate(G, separator, first, second, checkPreconditions);

		if (!result) {
			return false;
		}

		for (node n : separator) {
			assignments[n] = 0;
		}
		for (node n : first) {
			assignments[n] = 1;
		}
		for (node n : second) {
			assignments[n] = 2;
		}
		return true;
	}

	/**
	 * Returns the full name of this algorithm. Useful when running experiments. The name consists of the specific
	 * name of the algorithm and the names of the postprocessors.
	 *
	 * @return the full name
	 */
	virtual std::string getName() const {
		std::string name = getSpecificName();
		for (const auto& post : postProcessors) {
			name += "_" + post->getName();
		}
		return name;
	}

	/**
	 * Returns the exitPoint, i.e. a string describing the point at which the algorithm returned.
	 *
	 * @return the exitPoint
	 */
	std::string getExitPoint() const { return exitPoint; }

	/**
	 * Provides the maximal separator size that this algorithm guarantees as a function of the number of nodes of
	 * the graph, or a negative value if the guarantee cannot be expressed through such a function.
	 * See e.g. SeparatorHarPeled or SeparatorDualFC for examples.
	 *
	 * @param n the number of nodes of the graph
	 * @return the maximal separator size
	 */
	virtual double getMaxSeparatorSize(int n) const = 0;

protected:
	std::vector<Postprocessor*> postProcessors;

	std::shared_ptr<GraphCopy> graph;

	int startNodeIndex = -1;

	// the algorithms can return at different points - this variable stores where that happened, for analysis only
	std::string exitPoint;

	/**
	 * Selects the starting node for the BFS.
	 *
	 * @param G the graph to be solved
	 * @return random node if no desired index was set, node with that index otherwise
	 */
	node getStartNode(const Graph& G) const {
		if (startNodeIndex == -1) {
			return G.chooseNode();
		} else {
			return G.chooseNode([&](node n) { return (n->index() == startNodeIndex); });
		}
	}

	/**
	 * Core of the specific separation algorithm - override this in inheriting classes.
	 *
	 * @pre G is planar, simple-undirected, connected and represents a Combinatorial Embedding = is planarly embedded already.
	 * @param G the graph to be separated
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @return true on success
	 */
	virtual bool doSeparate(const Graph& G, List<node>& separator, List<node>& first,
			List<node>& second) = 0;

	/**
	 * Performs some initial setup to ensure that all preconditions hold and takes trivial steps to separate the
	 * graph. Asserts that the graph is planar, simple and undirected.
	 *
	 * @param G the graph to be separated
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @param checkPreconditions whether to check if the graph represents a comb-embedding and is simple and undirected
	 * @return whether the graph could already be solved by trivial operations (i.e. if true, we are done)
	 */
	bool setup(const Graph& G, List<node>& separator, List<node>& first, List<node>& second,
			bool checkPreconditions = true) {
		exitPoint = "graph_trivial";
		if (G.empty()) {
			return true; // an empty graph is separated without us doing anything
		}

		OGDF_ASSERT(isPlanar(G));

		graph = std::make_shared<GraphCopy>(G);

		// call separateComponents to check if we even need to do anything, return success if not
		if (separateComponents(*graph, separator, first, second)) {
			return true;
		}

		// even checking these conditions is expensive, so if you know that they all hold, you can skip this
		if (checkPreconditions) {
			if (!graph->representsCombEmbedding()) {
				planarEmbedPlanarGraph(*graph);
			}
			if (!isSimpleUndirected(*graph)) {
				makeSimpleUndirected(*graph);
			}
		}

		OGDF_ASSERT(graph->representsCombEmbedding());
		OGDF_ASSERT(isConnected(*graph));
		OGDF_ASSERT(isSimpleUndirected(*graph));

		return false;
	}

	/**
	 * Performs built-in post-processing: For small instances, it can happen that all nodes are assigned to the
	 * separator, while both components are empty, which can be fixed by moving half of the nodes to the first list.
	 *
	 * @param G the graph
	 * @param separator the separator
	 * @param first the first component
	 * @param second the second component
	 * @return whether the cleanup procedure was applied or not
	 */
	bool cleanup(const Graph& G, List<node>& separator, List<node>& first, List<node>& second) {
		if (first.empty() && second.empty()) {
			for (int i = 0; i < G.numberOfNodes() / 2.0; i++) {
				first.pushBack(separator.popFrontRet());
			}
			return true;
		}
		return false;
	}

	/**
	 * Checks if the graph consists of multiple connected components, takes necessary steps for fixing that,
	 * returns true if this already solved the graph, false if the core algorithm still needs to run.
	 *
	 * @param G the graph to be separated
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @param skip whether to skip the connectedness-step
	 * @return true if the graph could be separated already
	 */
	bool separateComponents(GraphCopy& G, List<node>& separator, List<node>& first,
			List<node>& second, bool skip = false) const;

	/**
	 * Apply all postprocessors.
	 *
	 * @param G the graph to be separated
	 * @param separator the separator nodes
	 * @param first the first component
	 * @param second the second component
	 * @return true on success
	 */
	bool postProcess(const Graph& G, List<node>& separator, List<node>& first, List<node>& second) {
		for (Postprocessor* post : postProcessors) {
			post->apply(G, separator, first, second);
		}
		return true;
	}

	/**
	 * Returns the unique name of the core algorithm, to be combined with postprocessors later.
	 * Override this in inheriting methods.
	 *
	 * @return the specific name as a string
	 */
	virtual std::string getSpecificName() const = 0;

	/**
	 * Reset everything to enable reuse of the module.
	 */
	virtual void reset() {};

private:
	/**
	 * Finds all connected components within the graph. Essentially a modified version of ogdf::connectedComponents
	 * that also fills a map to maintain component sizes.
	 *
	 * @param G the graph
	 * @param component assigns a component to each node
	 * @param compSizes maps component number to its size
	 * @return the number of connected components
	 */
	int connectedComponents(const Graph& G, NodeArray<int>& component,
			std::map<int, int>& compSizes) const;
};


}
