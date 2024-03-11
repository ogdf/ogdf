/** \file
 * \brief Utils for PCTree::allNodes(), PCTree::innerNodes(), PCNode::children() and PCNode::neighbors().
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
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

#include <ogdf/basic/pctree/PCEnum.h>
#include <ogdf/basic/pctree/PCNode.h>

#include <deque>
#include <utility>

namespace pc_tree {
class OGDF_EXPORT PCNodeIterator {
	friend struct PCNodeChildrenIterable;
	friend struct PCNodeNeighborsIterable;

	PCNode* node = nullptr;
	PCNode* pred = nullptr;
	PCNode* curr = nullptr;

	PCNodeIterator(PCNode* p_node, PCNode* p_pred, PCNode* p_curr)
		: node(p_node), pred(p_pred), curr(p_curr) { }

public:
	using iterator_category = std::forward_iterator_tag;
	using value_type = PCNode;
	using pointer = PCNode*;
	using reference = PCNode&;
	using difference_type = std::ptrdiff_t;

	PCNodeIterator() = default;

	PCNode& operator->() const { return *curr; }

	PCNode* operator*() const { return curr; }

	//! Increment operator (prefix, returns result).
	PCNodeIterator& operator++();

	//! Increment operator (postfix, returns previous value).
	PCNodeIterator operator++(int);

	bool operator==(const PCNodeIterator& rhs) const {
		return node == rhs.node && pred == rhs.pred && curr == rhs.curr;
	}

	bool operator!=(const PCNodeIterator& rhs) const { return !(rhs == *this); }

	PCNode* nodeOf() const { return node; }

	bool isParent();
};

struct OGDF_EXPORT PCNodeChildrenIterable {
	PCNode* const node;

	explicit PCNodeChildrenIterable(PCNode* p_node) : node(p_node) { }

	PCNodeIterator begin() const noexcept;

	PCNodeIterator end() const noexcept;

	unsigned long count() const;
};

struct OGDF_EXPORT PCNodeNeighborsIterable {
	PCNode* const node;
	PCNode* const first;

	explicit PCNodeNeighborsIterable(PCNode* p_node, PCNode* p_first = nullptr)
		: node(p_node)
		, first(p_first != nullptr
						  ? p_first
						  : (p_node->child1 != nullptr ? p_node->child1 : p_node->getParent())) {
		if (this->first == nullptr) {
			OGDF_ASSERT(this->node->getDegree() == 0);
		} else {
			OGDF_ASSERT(this->node->isParentOf(this->first) || this->first->isParentOf(this->node));
		}
	}

	PCNodeIterator begin() const noexcept;

	PCNodeIterator end() const noexcept;

	unsigned long count() const;
};

/**
 * A DFS or BFS through a PCTree.
 * @sa FilteringBFS
 */
template<bool dfs, bool reverse = false>
class FilteringPCTreeWalk {
	using container_type =
			typename std::conditional<dfs, std::vector<PCNode*>, std::deque<PCNode*>>::type;

	container_type m_pending;
	std::function<bool(PCNode*)> m_visit;
	std::function<bool(PCNode*)> m_descend;

public:
	// iterator traits
	using iterator_category = std::input_iterator_tag;
	using value_type = PCNode*;
	using difference_type = std::ptrdiff_t;
	using pointer = PCNode**;
	using reference = PCNode*&;

	static bool return_true([[maybe_unused]] PCNode* n) { return true; }

	explicit FilteringPCTreeWalk() = default;

	explicit FilteringPCTreeWalk([[maybe_unused]] const PCTree& T, PCNode* start,
			std::function<bool(PCNode*)> visit = return_true,
			std::function<bool(PCNode*)> descend_from = return_true)
		: m_pending({start}), m_visit(std::move(visit)), m_descend(std::move(descend_from)) {
		if (!m_pending.empty() && !m_visit(top())) {
			next();
		}
	}

	bool operator==(const FilteringPCTreeWalk& rhs) const { return m_pending == rhs.m_pending; }

	bool operator!=(const FilteringPCTreeWalk& rhs) const { return m_pending != rhs.m_pending; }

	FilteringPCTreeWalk& begin() { return *this; }

	FilteringPCTreeWalk end() const { return FilteringPCTreeWalk(); }

	PCNode* top() {
		OGDF_ASSERT(!m_pending.empty());
		if constexpr (dfs) {
			return m_pending.back();
		} else {
			return m_pending.front();
		}
	}

	PCNode* operator*() { return top(); }

	//! Increment operator (prefix, returns result).
	FilteringPCTreeWalk& operator++() {
		next();
		return *this;
	}

	//! Increment operator (postfix, returns previous value).
	OGDF_DEPRECATED("Calling FilteringPCTreeWalk++ will copy the array of pending nodes")

	FilteringPCTreeWalk operator++(int) {
		FilteringPCTreeWalk before = *this;
		next();
		return before;
	}

	void next() {
		do {
			OGDF_ASSERT(!m_pending.empty());
			PCNode* node = top();
			if constexpr (dfs) {
				m_pending.pop_back();
			} else {
				m_pending.pop_front();
			}
			if (m_descend(node)) {
				std::copy(node->children().begin(), node->children().end(),
						std::back_inserter(m_pending));
				if constexpr (reverse) {
					std::reverse(m_pending.end() - node->getChildCount(), m_pending.end());
				}
			}
		} while (!m_pending.empty() && !m_visit(top()));
	}

	explicit operator bool() const { return valid(); }

	bool valid() const { return !m_pending.empty(); }

	void append(PCNode* a) { m_pending.push(a); }

	int pendingCount() const { return m_pending.size(); }
};

using FilteringPCTreeDFS = FilteringPCTreeWalk<true>;
using FilteringPCTreeBFS = FilteringPCTreeWalk<false>;
}
