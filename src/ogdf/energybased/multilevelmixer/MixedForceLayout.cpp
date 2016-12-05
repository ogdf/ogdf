/** \file
 * \brief Uses Fruchtermann Rheingold and Fast Multipole Embedder for faster and better FR results.
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

#include <ogdf/energybased/multilevelmixer/MixedForceLayout.h>
#include <ogdf/energybased/SpringEmbedderGridVariant.h>
#include <ogdf/energybased/FastMultipoleEmbedder.h>

namespace ogdf {

MixedForceLayout::MixedForceLayout()
{
	SpringEmbedderGridVariant * FR = new SpringEmbedderGridVariant();
	FR->scaling(SpringEmbedderGridVariant::Scaling::input);
	m_FR = FR;

	FastMultipoleEmbedder * FME = new FastMultipoleEmbedder();
	FME->setNumIterations(1000);
	FME->setRandomize(false);
	FME->setNumberOfThreads(2);

	m_FME = FME;
}


void MixedForceLayout::call(GraphAttributes &GA)
{
	m_FME->call(GA);
	m_FR->call(GA);
}


void MixedForceLayout::call(MultilevelGraph &MLG)
{
	m_FME->call(MLG.getGraphAttributes());
	m_FR->call(MLG.getGraphAttributes());
}

} // namespace ogdf
