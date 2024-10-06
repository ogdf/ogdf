/** \file
 * \brief A registry that allows labelling the nodes of a PC-tree.
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

#include <ogdf/basic/basic.h>
#include <ogdf/basic/pctree/PCEnum.h>

namespace ogdf::pc_tree {
class PCNode;
class PCTreeForest;

/**
 * A registry that allows labelling the nodes of a PC-tree.
 */
class OGDF_EXPORT PCTreeRegistry : public ogdf::RegistryBase<PCNode*, PCTreeRegistry> {
	PCTreeForest* m_pForest;

public:
	PCTreeRegistry(PCTreeForest* pcTreeForest) : m_pForest(pcTreeForest) { }

	//! Returns the index of \p key.
	static inline int keyToIndex(PCNode* key);

	//! Returns whether \p key is associated with this registry.
	bool isKeyAssociated(PCNode* key) const;

	//! Returns the maximum index of all keys managed by this registry.
	int maxKeyIndex() const;

	//! Returns the array size currently requested by this registry.
	int calculateArraySize(int add) const;

	operator PCTreeForest&() const { return *m_pForest; }

	operator PCTreeForest*() const { return m_pForest; }

	PCTreeForest& getForest() const { return *m_pForest; }
};
}
