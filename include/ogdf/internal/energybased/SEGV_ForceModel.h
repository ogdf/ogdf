/** \file
 * \brief Declarations of force-models for Spring-Embedder algorithm
 *
 * \author Carsten Gutwenger
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
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
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#pragma once

#include <ogdf/energybased/SpringEmbedderGridVariant.h>

namespace ogdf {

class SpringEmbedderGridVariant::ForceModelBase
{
public:
	ForceModelBase(const Array<NodeInfo> &vInfo, const Array<int> &adjLists, const Array2D<ListPure<int>> &gridCell, double idealEdgeLength)
		: m_vInfo(vInfo), m_adjLists(adjLists), m_gridCell(gridCell), m_idealEdgeLength(idealEdgeLength) { }

	virtual ~ForceModelBase() { }

	virtual DPoint computeDisplacement(int j, double boxLength) const = 0;

	double eps() const { return 0.01*m_idealEdgeLength; }

protected:
	const Array<NodeInfo>        &m_vInfo;
	const Array<int>             &m_adjLists;
	const Array2D<ListPure<int>> &m_gridCell;

	double m_idealEdgeLength;
};


class SpringEmbedderGridVariant::ForceModelFR : public ForceModelBase
{
public:
	ForceModelFR(const Array<NodeInfo> &vInfo, const Array<int> &adjLists, const Array2D<ListPure<int>> &gridCell, double idealEdgeLength)
		: ForceModelBase(vInfo, adjLists, gridCell, idealEdgeLength) { }

	DPoint computeDisplacement(int j, double boxLength) const;
};

class SpringEmbedderGridVariant::ForceModelFRModAttr : public ForceModelBase
{
public:
	ForceModelFRModAttr(const Array<NodeInfo> &vInfo, const Array<int> &adjLists, const Array2D<ListPure<int>> &gridCell, double idealEdgeLength)
		: ForceModelBase(vInfo, adjLists, gridCell, idealEdgeLength) { }

	DPoint computeDisplacement(int j, double boxLength) const;
};

class SpringEmbedderGridVariant::ForceModelFRModRep : public ForceModelBase
{
public:
	ForceModelFRModRep(const Array<NodeInfo> &vInfo, const Array<int> &adjLists, const Array2D<ListPure<int>> &gridCell, double idealEdgeLength)
		: ForceModelBase(vInfo, adjLists, gridCell, idealEdgeLength) { }

	DPoint computeDisplacement(int j, double boxLength) const;
};


class SpringEmbedderGridVariant::ForceModelEades : public ForceModelBase
{
public:
	ForceModelEades(const Array<NodeInfo> &vInfo, const Array<int> &adjLists, const Array2D<ListPure<int>> &gridCell, double idealEdgeLength)
		: ForceModelBase(vInfo, adjLists, gridCell, idealEdgeLength) { }

	DPoint computeDisplacement(int j, double boxLength) const;
};


class SpringEmbedderGridVariant::ForceModelHachul : public ForceModelBase
{
public:
	ForceModelHachul(const Array<NodeInfo> &vInfo, const Array<int> &adjLists, const Array2D<ListPure<int>> &gridCell, double idealEdgeLength)
		: ForceModelBase(vInfo, adjLists, gridCell, idealEdgeLength) { }

	DPoint computeDisplacement(int j, double boxLength) const;
};


class SpringEmbedderGridVariant::ForceModelGronemann : public ForceModelBase
{
public:
	ForceModelGronemann(const Array<NodeInfo> &vInfo, const Array<int> &adjLists, const Array2D<ListPure<int>> &gridCell, double idealEdgeLength)
		: ForceModelBase(vInfo, adjLists, gridCell, idealEdgeLength) { }

	DPoint computeDisplacement(int j, double boxLength) const;
};

}
