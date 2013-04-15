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

#include <ogdf/abacus/setbranchrule.h>
#include <ogdf/abacus/master.h>
#include <ogdf/abacus/sub.h>

namespace abacus {


ostream &operator<<(ostream &out, const SetBranchRule &rhs)
{
	return out << "x" << rhs.variable_ << " = " << (int) rhs.status_;
}


int SetBranchRule::extract(Sub *sub)
{
	if (sub->fsVarStat(variable_)->contradiction(status_))
		return 1;

	sub->fsVarStat(variable_)->status(status_);
	return 0;
}


void SetBranchRule::extract(LpSub *lp)
{
	if (status_ == FSVarStat::SetToLowerBound) {
		oldLpBound_ = lp->uBound(variable_);
		lp->changeUBound(variable_, lp->lBound(variable_));
	}
	else {
		oldLpBound_ = lp->lBound(variable_);
		lp->changeLBound(variable_, lp->uBound(variable_));
	}
}


void SetBranchRule::unExtract(LpSub *lp)
{
	if (status_ == FSVarStat::SetToLowerBound)
		lp->changeUBound(variable_, oldLpBound_);
	else
		lp->changeLBound(variable_, oldLpBound_);
}
} //namespace abacus
