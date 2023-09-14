/** \file
 *
 * \author Marcel Radermacher
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

#include <vector>

namespace ogdf {
namespace internal {
namespace gcm {

namespace datastructure {

template<typename T, typename Graph>
class NodeVector : public std::vector<T> {
private:
	using Node = typename Graph::Node;
	using reference = typename std::vector<T>::reference;
	using const_reference = typename std::vector<T>::const_reference;

	const Graph* g;
	T default_value;


public:
	NodeVector() : g(nullptr) { }

	NodeVector(const Graph& _g) : std::vector<T>(_g.max_node_index() + 1), g(&_g) { }

	NodeVector(const Graph& _g, T& v)
		: std::vector<T>(_g.max_node_index() + 1, v), g(&_g), default_value(v) { }

	NodeVector(const Graph& _g, T v)
		: std::vector<T>(_g.max_node_index() + 1, v), g(&_g), default_value(v) { }

	NodeVector(const NodeVector<T, Graph>& x)
		: std::vector<T>(x), g(x.g), default_value(x.default_value) { }

	inline void adapt() {
		if ((std::size_t)g->max_node_index() >= std::vector<T>::size()) {
			std::vector<T>::resize(g->max_node_index() + 1, default_value);
		}
	}

	reference operator[](const Node& v) {
		adapt();
		return std::vector<T>::operator[](v->index());
	}

	const_reference operator[](const Node& v) const {
		OGDF_ASSERT((std::size_t)v->index() < this->size());
		return std::vector<T>::operator[](v->index());
	}

	NodeVector<T, Graph>& operator=(const NodeVector<T, Graph>& x) {
		std::vector<T>::operator=(x);
		default_value = x.default_value;
		return *this;
	}

	NodeVector<T, Graph>& operator=(NodeVector<T, Graph>&& x) {
		std::vector<T>::operator=(std::move(x));
		default_value = x.default_value;
		return *this;
	}
};

template<typename T, typename Graph>
class EdgeVector : public std::vector<T> {
private:
	using Edge = typename Graph::Edge;
	using reference = typename std::vector<T>::reference;
	using const_reference = typename std::vector<T>::const_reference;
	const Graph* g;
	T default_value;

public:
	EdgeVector() : g(nullptr) { }

	EdgeVector(const Graph& _g) : std::vector<T>(_g.max_edge_index() + 1), g(&_g) { }

	EdgeVector(const Graph& _g, T& v)
		: std::vector<T>(_g.max_edge_index() + 1, v), g(&_g), default_value(v) { }

	EdgeVector(const Graph& _g, T v)
		: std::vector<T>(_g.max_edge_index() + 1, v), g(&_g), default_value(v) { }

	EdgeVector(const EdgeVector& x) : std::vector<T>(x), g(&x.g), default_value(x.default_value) { }

	EdgeVector(const EdgeVector&& x)
		: std::vector<T>(x), g(std::move(x.g)), default_value(x.default_value) { }

	inline void adapt() {
		if (g && (std::size_t)g->max_edge_index() >= std::vector<T>::size()) {
			std::vector<T>::resize(g->max_edge_index() + 1, default_value);
		}
	}

	reference operator[](const Edge& e) {
		adapt();
		return std::vector<T>::operator[](e->index());
	}

	const_reference operator[](const Edge& e) const {
		OGDF_ASSERT((std::size_t)e->index() < this->size());
		return std::vector<T>::operator[](e->index());
	}

	EdgeVector<T, Graph>& operator=(const EdgeVector<T, Graph>& x) {
		std::vector<T>::operator=(x);
		default_value = x.default_value;
		return *this;
	}

	EdgeVector<T, Graph>& operator=(EdgeVector<T, Graph>&& x) {
		std::vector<T>::operator=(std::move(x));
		default_value = x.default_value;
		return *this;
	}
};

}

}
}
}
