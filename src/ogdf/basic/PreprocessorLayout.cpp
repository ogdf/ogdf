/** \file
 * \brief Preprocessor Layout simplifies Graphs for use in other Algorithms
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

#include <ogdf/basic/PreprocessorLayout.h>

namespace ogdf {

PreprocessorLayout::PreprocessorLayout()
: m_randomize(false)
{

}


void PreprocessorLayout::call(GraphAttributes &GA)
{
	if (m_secondaryLayout.valid()) {
		MultilevelGraph MLG(GA);
		call(MLG);
		MLG.exportAttributes(GA);
	}
}


void PreprocessorLayout::call(MultilevelGraph &MLG)
{
	m_deletedEdges.clear();
	Graph * G = &(MLG.getGraph());

	double sqrsize;
	if (m_randomize) sqrsize = 2.0*sqrt((double)G->numberOfNodes())*MLG.averageRadius();

	for(node v : G->nodes) {
		if (MLG.radius(v) <= 0) {
			MLG.radius(v, 1.0);
		}
		if (m_randomize) {
			MLG.x(v, randomDouble( -sqrsize, sqrsize ));//-5.0, 5.0));
			MLG.y(v, randomDouble( -sqrsize, sqrsize ));
		}
	}
	if (m_secondaryLayout.valid()) {

		call(*G, MLG);

		m_secondaryLayout.get().call(MLG.getGraphAttributes());
		MLG.updateReverseIndizes();

		for( const EdgeData &ed : m_deletedEdges ) {
			int index = ed.edgeIndex;
			edge temp = G->newEdge(MLG.getNode(ed.sourceIndex), MLG.getNode(ed.targetIndex), index);
			MLG.weight(temp, (float)ed.weight);
		}
	}
}


void PreprocessorLayout::call(Graph &G, MultilevelGraph &MLG)
{
	std::vector<edge> deletedEdges;

	for(edge e : G.edges) {
		int index = e->index();
		if (e->isSelfLoop()) {
			deletedEdges.push_back(e);
			m_deletedEdges.push_back(EdgeData(index, e->source()->index(), e->target()->index(), MLG.weight(e)));
		} else {
			for(adjEntry adj : e->source()->adjEdges) {
				if (adj->theEdge()->index() < index && adj->twinNode() == e->target()) {
					deletedEdges.push_back(e);
					m_deletedEdges.push_back(EdgeData(index, e->source()->index(), e->target()->index(), MLG.weight(e)));
					break;
				}
			}
		}
	}

	for (edge e : deletedEdges) {
		G.delEdge(e);
	}
}


} // namespace ogdf
