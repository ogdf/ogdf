/** \file
 * \brief Places Nodes at the Positio of the merge-partner
 *
 * \author Gereon Bartel
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

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphList.h>
#include <ogdf/basic/basic.h>
#include <ogdf/energybased/multilevel_mixer/MedianPlacer.h>
#include <ogdf/energybased/multilevel_mixer/MultilevelGraph.h>

#include <algorithm>
#include <vector>

namespace ogdf {

void MedianPlacer::placeOneLevel(MultilevelGraph& MLG) {
	int level = MLG.getLevel();
	while (MLG.getLevel() == level && MLG.getLastMerge() != nullptr) {
		placeOneNode(MLG);
	}
}

void MedianPlacer::placeOneNode(MultilevelGraph& MLG) {
	node merged = MLG.undoLastMerge();
	int i = 0;
	std::vector<double> xVector;
	std::vector<double> yVector;
	for (adjEntry adj : merged->adjEntries) {
		i++;
		xVector.push_back(MLG.x(adj->twinNode()));
		yVector.push_back(MLG.y(adj->twinNode()));
	}
	std::nth_element(xVector.begin(), xVector.begin() + (i / 2), xVector.end());
	std::nth_element(yVector.begin(), yVector.begin() + (i / 2), yVector.end());
	double x = xVector[i / 2];
	double y = yVector[i / 2];
	if (i % 2 == 0) {
		std::nth_element(xVector.begin(), xVector.begin() + (i / 2) - 1, xVector.end());
		std::nth_element(yVector.begin(), yVector.begin() + (i / 2) - 1, yVector.end());
		x += xVector[i / 2 - 1];
		y += yVector[i / 2 - 1];
		x /= 2.0;
		y /= 2.0;
	}
	MLG.x(merged, x + ((m_randomOffset) ? (float)randomDouble(-1.0, 1.0) : 0.f));
	MLG.y(merged, y + ((m_randomOffset) ? (float)randomDouble(-1.0, 1.0) : 0.f));
}

}
