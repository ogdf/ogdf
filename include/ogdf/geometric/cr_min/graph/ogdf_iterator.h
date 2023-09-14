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

#include <ogdf/basic/Graph_d.h>

namespace ogdf {
namespace internal {
namespace gcm {
namespace graph {

/*! iterate over entries of an ogdf graph
 */
template<typename _Entry>
class OGDFEntryIterator {
public:
	using Entry = _Entry;
	using Iterator = OGDFEntryIterator<Entry>;

protected:
	Entry m_cur;

public:
	OGDFEntryIterator(Entry _cur) : m_cur(_cur) {
		//nothing to do
	}

	Iterator& operator--() {
		m_cur = m_cur->pred();
		return *this;
	}

	Iterator& operator++() {
		m_cur = m_cur->succ();
		return *this;
	}

	bool operator==(const Iterator& b) const { return m_cur == b.m_cur; }

	bool operator!=(const Iterator& b) const { return m_cur != b.m_cur; }
};

template<typename Entry>
class EntryIterator : public OGDFEntryIterator<Entry> {
private:
	using parent = OGDFEntryIterator<Entry>;

public:
	using T = Entry;
	using parent::parent;

	Entry operator*() { return this->m_cur; }
};

/*! iterate over incident edges */
class IncidentEdgeIterator : public OGDFEntryIterator<adjEntry> {
private:
	using parent = OGDFEntryIterator<adjEntry>;

public:
	using T = edge;
	using parent::parent;

	edge operator*() { return this->m_cur->theEdge(); }
};

/* Iterate over adjacent nodes
 */
class AdjacentNodeIterator : public OGDFEntryIterator<adjEntry> {
private:
	using parent = OGDFEntryIterator<adjEntry>;

public:
	using T = node;
	using parent::parent;

	node operator*() { return this->m_cur->twin()->theNode(); }
};

using AdjEntryIterator = EntryIterator<adjEntry>;
using NodeIterator = EntryIterator<node>;
using EdgeIterator = EntryIterator<edge>;


}
}
}
}
