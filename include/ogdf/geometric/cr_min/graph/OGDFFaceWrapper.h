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

#ifdef OGDF_INCLUDE_CGAL

#	include <ogdf/basic/CombinatorialEmbedding.h>
#	include <ogdf/geometric/cr_min/graph/OGDFGraphWrapper.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

/*!
 *! Supplies easy access to a ogdf face element. Be carefull reserves space in the size of the number of nodes of the graph!
 */
class OGDFFaceWrapper {
private:
	using Node = node;
	face m_face;
	std::shared_ptr<CombinatorialEmbedding> m_ce;
	std::vector<unsigned int> is_on_flag;
	std::vector<Node> nodes; //TODO defince cyclic iterator?
public:
	OGDFFaceWrapper() : m_ce(new CombinatorialEmbedding()) { }

	CombinatorialEmbedding& get_combinatorial_embedding() { return *m_ce; }

	void set_face(adjEntry external_entry) {
		set_face(get_combinatorial_embedding().leftFace(external_entry));
	}

	void set_face(face _face) {
		is_on_flag.assign(m_ce->getGraph().maxNodeIndex() + 1, false);
		m_face = _face;
		unsigned int c = 1;
		nodes.clear();
		for_all_nodes([&](const Node v) {
			is_on_flag[v->index()] = c++;
			nodes.push_back(v);
		});
	}

	face ogdf_face() { return m_face; }

	adjEntry first_adj() { return m_face->firstAdj(); }

	std::vector<Node>::iterator begin() { return nodes.begin(); }

	std::vector<Node>::const_iterator begin() const { return nodes.cbegin(); }

	std::vector<Node>::iterator end() { return nodes.end(); }

	std::vector<Node>::const_iterator end() const { return nodes.cend(); }

	template<typename Handler>
	void for_all_nodes(Handler&& handler) const {
		adjEntry current = m_face->firstAdj();
		do {
			handler(current->theNode());
			current = current->faceCycleSucc();
			OGDF_ASSERT(current != nullptr);
		} while (current != m_face->firstAdj());
	}

	bool has_on(Node v) const { return is_on_flag[v->index()]; }

	unsigned int ordering_number(Node v) const { return is_on_flag[v->index()]; }

	size_t number_of_nodes() const { return nodes.size(); }
};

}
}
}
}

#endif
