/** \file
 * \brief // TODO DESCRIBE WHAT IS IMPLEMENTED
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

#include <ostream>

namespace pc_tree {
enum class NodeLabel { Unknown, Partial, Full, Empty = Unknown };

enum class PCNodeType { PNode, CNode, Leaf };

class PCTree;

class PCTreeForest;

template<class Key>
class PCTreeRegistry;

class PCNode;

#define OGDF_DECL_REG_ARRAY_TYPE(v, c) ogdf::RegisteredArray<PCTreeRegistry<PCNode*>, v, c>
OGDF_DECL_REG_ARRAY(PCTreeNodeArray)
#undef OGDF_DECL_REG_ARRAY_TYPE

template<bool SupportFastSizeQuery = true>
using PCTreeNodeSet = ogdf::RegisteredSet<PCTreeRegistry<PCNode*>, SupportFastSizeQuery>;
}

std::ostream& operator<<(std::ostream&, pc_tree::NodeLabel);

std::ostream& operator<<(std::ostream&, pc_tree::PCNodeType);

std::ostream& operator<<(std::ostream&, const pc_tree::PCTree*);

std::ostream& operator<<(std::ostream&, const pc_tree::PCNode*);

std::ostream& operator<<(std::ostream&, const pc_tree::PCTree&);

std::ostream& operator<<(std::ostream&, const pc_tree::PCNode&);
