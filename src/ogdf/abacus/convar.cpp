/*
 * $Revision: 3386 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-10 14:00:02 +0200 (Mi, 10. Apr 2013) $
 ***************************************************************/

/*!\file
 * \author Matthias Elf
 *
 * \par License:
 * This file is part of ABACUS - A Branch And CUt System
 * Copyright (C) 1995 - 2003
 * University of Cologne, Germany
 *
 * \par
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * \par
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * \see http://www.gnu.org/copyleft/gpl.html
 */

#include <ogdf/abacus/convar.h>
#include <ogdf/abacus/master.h>

namespace abacus {


void ConVar::_expand() const
{
	if(expanded_) {
		Logger::ifout() << "WARNING: ConVar::_expand(): ";
		Logger::ifout() << "constraint already expanded" << endl;
		return;
	}
	expand();
	expanded_ = true;
}


void ConVar::_compress() const
{
	if(!expanded_) {
		Logger::ifout() << "WARNING: ConVar::_compress(): ";
		Logger::ifout() << "constraint already compressed" << endl;
		return;
	}
	compress();
	expanded_ = false;
}


void ConVar::print(ostream &out) const
{
	out << "ConVar::print() is only a dummy." << endl;
}


unsigned ConVar::hashKey() const
{
	Logger::ifout() << "ConVar::hashKey() must be defined in derived class.\n";
	OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcConvar);
}


const char *ConVar::name() const
{
	Logger::ifout() << "ConVar::name() must be defined in derived class.\n";
	OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcConvar);
}


bool ConVar::equal(const ConVar * /* cv */) const
{
	Logger::ifout() << "ConVar::equal() must be defined in derived class.\n";
	OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcConvar);
}
} //namespace abacus
