/** \file
 * \brief Predeclaration of various PC-tree related classes and enums.
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

#include <ogdf/basic/RegisteredArray.h>
#include <ogdf/basic/RegisteredSet.h>
#include <ogdf/basic/basic.h>

#include <ostream>

namespace ogdf::pc_tree {
enum class NodeLabel { Unknown, Partial, Full, Empty = Unknown };

enum class PCNodeType { PNode, CNode, Leaf };

class PCNode;
class PCTree;
class PCTreeRegistry;

#define OGDF_DECL_REG_ARRAY_TYPE(v, c) ogdf::RegisteredArray<PCTreeRegistry, v, c>
OGDF_DECL_REG_ARRAY(PCTreeNodeArray)
#undef OGDF_DECL_REG_ARRAY_TYPE

using PCTreeNodeSet = ogdf::RegisteredSet<PCTreeRegistry>;

OGDF_EXPORT std::ostream& operator<<(std::ostream&, ogdf::pc_tree::NodeLabel);

OGDF_EXPORT std::ostream& operator<<(std::ostream&, ogdf::pc_tree::PCNodeType);

OGDF_EXPORT std::ostream& operator<<(std::ostream&, const ogdf::pc_tree::PCTree*);

OGDF_EXPORT std::ostream& operator<<(std::ostream&, const ogdf::pc_tree::PCNode*);

OGDF_EXPORT std::ostream& operator<<(std::ostream&, const ogdf::pc_tree::PCTree&);

OGDF_EXPORT std::ostream& operator<<(std::ostream&, const ogdf::pc_tree::PCNode&);
}
