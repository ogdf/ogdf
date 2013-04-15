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
 *
 * $Id: tailoff.cc,v 2.6 2008-10-27 11:02:25 baumann Exp $
 */

#include <ogdf/abacus/tailoff.h>

namespace abacus {


ostream &operator<<(ostream &out, const TailOff &rhs)
{
	out << "LP-history:" << endl;
	if (rhs.lpHistory_)
		out << *(rhs.lpHistory_);
	else
		out << "no LP-history available" << endl;

	return out;
}


bool TailOff::tailOff() const
{
	if (!lpHistory_) return false;

	if (!lpHistory_->filled()) return false;  //!< not enough iterations

	//FIXME
	double den = fabs(lpHistory_->oldest()) < 1e-30 ? 1e-30 : lpHistory_->oldest();

	if (fabs((lpHistory_->oldest() - lpHistory_->newest())*100.0
		/den)
		< master_->tailOffPercent()) return true;
	else return false;
}


int TailOff::diff(int nLps, double &d) const
{
	double oldVal;
	if (lpHistory_->previous(nLps, oldVal))
		return 1;

	double lastVal = lpHistory_->newest();

	d = fabs((lastVal - oldVal)*100.0/oldVal);

	return 0;
}
} //namespace abacus
