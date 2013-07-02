/*
 * $Revision: 3553 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-06-07 19:36:02 +0200 (Fr, 07. Jun 2013) $
 ***************************************************************/

/** \file
 * \brief Merges nodes with neighbour to get a Multilevel Graph
 *
 * \author Gereon Bartel
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_RANDOM_MERGER_H
#define OGDF_RANDOM_MERGER_H

#include <ogdf/energybased/multilevelmixer/MultilevelBuilder.h>

namespace ogdf {

class OGDF_EXPORT RandomMerger : public MultilevelBuilder
{
private:
	double m_levelSizeFactor;

	bool buildOneLevel(MultilevelGraph &MLG);

public:
	RandomMerger();
	void setFactor(double factor);
};

} // namespace ogdf

#endif
