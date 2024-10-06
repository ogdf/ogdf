/** \file
 * \brief Declaration of FourBlockTree.
 *
 * Based on the implementation and techiques of the following papers:
 *
 * Norishige Chiba, and Takao Nishizeki.
 * Arboricity and subgraph listing algorithms.
 * SIAM Journal on computing 14.1 (1985): 210-223.
 * [doi:10.1137/0214017](https://doi.org/10.1137/0214017).
 *
 * Sabine Cornelsen, and Gregor Diatzko.
 * Decomposing Triangulations into 4-Connected Components.
 * arXiv [cs.DS], 2023.
 * [doi:10.48550/arXiv.2308.16020](https://doi.org/10.48550/arXiv.2308.16020).
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/basic.h>

#include <memory>
#include <vector>

namespace ogdf {

/**
 * A node in a 4-block tree.
 *
 * Since each node contains its children, the root is the entire tree.
 */
struct OGDF_EXPORT FourBlockTree {
	FourBlockTree() = default;
	FourBlockTree(const FourBlockTree&) = delete;
	FourBlockTree(FourBlockTree&&) = delete;
	FourBlockTree& operator=(const FourBlockTree&) = delete;
	FourBlockTree& operator=(FourBlockTree&&) = delete;

	~FourBlockTree() {
		// free manually bottom-up to avoid stack overflow in case of deep tree
		postorder([](FourBlockTree& treeNode) -> void { treeNode.children.clear(); });
	}

	/**
	 * The 4-connected component.
	 */
	std::unique_ptr<Graph> g = std::make_unique<Graph>();

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
	static std::unique_ptr<FourBlockTree> construct(const Graph& g, adjEntry externalFace);

	/**
	 * Perform a pre-order traversal of the 4-block tree.
	 *
	 * Each child is processed after its parent.
	 *
	 * @tparam _F The type of callback, something like
	 *            `void (*)(const FourBlockTree&)`.
	 * @param callback The function to be called for each node of the tree.
	 */
	template<typename _F>
	void preorder(_F callback) const {
		struct stackEntry {
			const FourBlockTree* node;
			std::vector<std::unique_ptr<FourBlockTree>>::const_iterator nextChild;
		};

		std::vector<stackEntry> stack;
		stack.push_back({this, children.begin()});
		callback(*this);
		while (!stack.empty()) {
			auto& it = stack.back().nextChild;
			if (it != stack.back().node->children.end()) {
				const FourBlockTree* child = it->get();
				++it;
				stack.push_back({child, child->children.begin()});
				callback(*child);
			} else {
				stack.pop_back();
			}
		}
	}

	/**
	 * Perform a pre-order traversal of the 4-block tree.
	 *
	 * Each child is processed after its parent.
	 *
	 * @tparam _F The type of callback, something like
	 *            `void (*)(FourBlockTree&)`.
	 * @param callback The function to be called for each node of the tree.
	 */
	template<typename _F>
	void preorder(_F callback) {
		struct stackEntry {
			FourBlockTree* node;
			std::vector<std::unique_ptr<FourBlockTree>>::iterator nextChild;
		};

		std::vector<stackEntry> stack;
		stack.push_back({this, children.begin()});
		callback(*this);
		while (!stack.empty()) {
			auto& it = stack.back().nextChild;
			if (it != stack.back().node->children.end()) {
				FourBlockTree* child = it->get();
				++it;
				stack.push_back({child, child->children.begin()});
				callback(*child);
			} else {
				stack.pop_back();
			}
		}
	}

	/**
	 * Perform a post-order traversal of the 4-block tree.
	 *
	 * Each child is processed before its parent.
	 *
	 * @tparam _F The type of callback, something like
	 *            `void (*)(const FourBlockTree&)`.
	 * @param callback The function to be called for each node of the tree.
	 */
	template<typename _F>
	void postorder(_F callback) const {
		struct stackEntry {
			const FourBlockTree* node;
			std::vector<std::unique_ptr<FourBlockTree>>::const_iterator nextChild;
		};

		std::vector<stackEntry> stack;
		stack.push_back({this, children.begin()});
		while (!stack.empty()) {
			auto& it = stack.back().nextChild;
			if (it != stack.back().node->children.end()) {
				const FourBlockTree* child = it->get();
				++it;
				stack.push_back({child, child->children.begin()});
			} else {
				callback(*stack.back().node);
				stack.pop_back();
			}
		}
	}

	/**
	 * Perform a post-order traversal of the 4-block tree.
	 *
	 * Each child is processed before its parent.
	 *
	 * @tparam _F The type of callback, something like
	 *            `void (*)(FourBlockTree&)`.
	 * @param callback The function to be called for each node of the tree.
	 */
	template<typename _F>
	void postorder(_F callback) {
		struct stackEntry {
			FourBlockTree* node;
			std::vector<std::unique_ptr<FourBlockTree>>::iterator nextChild;
		};

		std::vector<stackEntry> stack;
		stack.push_back({this, children.begin()});
		while (!stack.empty()) {
			auto& it = stack.back().nextChild;
			if (it != stack.back().node->children.end()) {
				FourBlockTree* child = it->get();
				++it;
				stack.push_back({child, child->children.begin()});
			} else {
				callback(*stack.back().node);
				stack.pop_back();
			}
		}
	}
};

} // namespace ogdf
