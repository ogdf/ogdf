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
     * The IDs of the vertices in the original graph as specified in the second
     * argument to the constructor of FourBlockTreeBuilder.
     *
     * Since vertices may appear in multiple 4-connected components, these IDs
     * need not be unique across nodes of the 4-block tree.
     */
    NodeArray<size_t> vertexIds;

    /**
     * A half-edge in g such that the external face of g is to its left.
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
};

/**
 * A class that constructs the 4-block tree of a given graph.

 * For details see https://doi.org/10.48550/arXiv.2308.16020
 */
class FourBlockTreeBuilder {
    using triangle_t = std::array<edge, 3>;
    using returnSide_t = unsigned char;
    static constexpr returnSide_t LEFT = 0b10, RIGHT = 0b01; // when root is drawn at the top

    Graph& m_g;
    NodeArray<size_t>& m_vertexIds;
    adjEntry m_externalFace;

    /**
     * A node in the external face.
     *
     * This will be used as the root for DFS.
     */
    node m_root;

    /**
     * Index of every adjEntry in its adjencency list.
     */
    AdjEntryArray<size_t> m_indices;

    /**
     * Populate m_indices.
     */
    void populateIndices();

    /**
     * All separating triangles.
     */
    std::vector<triangle_t> m_sepTriangles;

    /**
     * Populate m_sepTriangles.
     */
    void populateSepTriangles();

    /**
     * The length of the root-v-path along tree edges for each vertex v.
     */
    NodeArray<size_t> m_depth;

    /**
     * The half-edge along which each vertex was found during the first DFS.
     *
     * m_parentEdge[v]->theNode() = v
     */
    NodeArray<adjEntry> m_parentEdge; // m_parentEdge[v]->theNode() = v

    /**
     * Whether each edge is a tree edge or a back edge.
     */
    EdgeArray<bool> m_isTreeEdge;

    /**
     * The lowpoint of each edge.
     */
    EdgeArray<node> m_lowpoint;

    /**
     * The return side of each edge.
     */
    EdgeArray<returnSide_t> m_returnSide;

    /**
     * The angular distance of each edge.
     */
    EdgeArray<size_t> m_angularDistance;

    /**
     * Run the first DFS.
     *
     * Populate m_depth, m_parentEdge, m_isTreeEdge, m_lowpoint, m_returnSide,
     * m_angularDistance.
     */
    void firstDfs();

    /**
     * Order the elements of m_sepTriangles inside out.
     */
    void orderTriangles();

    /**
     * Split m_g along the elements of m_sepTriangles and build the resulting
     * 4-block tree.
     */
    FourBlockTree buildTree();

public:
    /**
     * Prepares all necessary data structures.
     *
     * @param g The plane triangulated graph whose 4-block tree shall be constructed.
     *          This graph will be used destructively.
     *          Edge directions in g are not respected.
     *          The order of edges at each vertex is used as the combinatorial
     *          embedding.
     * @param externalFace A half-edge in g such that the external face of g
     *                     lies to its left.
     */
    FourBlockTreeBuilder(Graph& g, NodeArray<size_t>& vertexIds, adjEntry externalFace);

    /**
     * Run the algorithm.
     */
    FourBlockTree call();
};

} // namespace ogdf
