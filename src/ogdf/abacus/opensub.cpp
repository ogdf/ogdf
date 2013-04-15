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

#include <ogdf/abacus/opensub.h>
#include <ogdf/abacus/master.h>
#include <ogdf/abacus/optsense.h>
#include <ogdf/abacus/sub.h>

namespace abacus {


void OpenSub::insert(Sub *sub)
{
	// update the dual bound of all subproblems in the \a list_
	if (empty())
		dualBound_ = sub->dualBound();
	else {
		if (master_->optSense()->max()) {
			if (sub->dualBound() > dualBound_)
				dualBound_ = sub->dualBound();
		}
		else
			if (sub->dualBound() < dualBound_)
				dualBound_ = sub->dualBound();
	}

	list_.pushBack(sub);
}


Sub* OpenSub::select()
{
	if(list_.empty())
		return 0;

	ogdf::ListIterator<Sub*> itMin = list_.begin();

	for(ogdf::ListIterator<Sub*> it = list_.begin(); it.valid(); ++it) {
		Sub *s = *it;
		if (s->status() == Sub::Dormant) {
			s->newDormantRound();
			if (s->nDormantRounds() < master_->minDormantRounds())
				continue;
		}
		if (master_->enumerationStrategy(s, *itMin) > 0)
			itMin = it;
	}
	Sub* min = *itMin;
	list_.del(itMin);

	updateDualBound();

	return min;
}


double OpenSub::dualBound() const
{
	double ret;
	if (empty()) {
		if (master_->optSense()->max()) ret = -master_->infinity();
		else                            ret = master_->infinity();
	}
	else
		ret = dualBound_;
	return ret;
}


void OpenSub::updateDualBound()
{
	if (master_->optSense()->max()) {
		dualBound_ = -master_->infinity();

		for(ogdf::ListIterator<Sub*> it = list_.begin(); it.valid(); ++it) {
			Sub *s = *it;
			if (s->dualBound() > dualBound_)
				dualBound_ = s->dualBound();
		}
	}
	else {
		dualBound_ = master_->infinity();

		for(ogdf::ListIterator<Sub*> it = list_.begin(); it.valid(); ++it) {
			Sub *s = *it;
			if (s->dualBound() < dualBound_)
				dualBound_ = s->dualBound();
		}
	}
}
} //namespace abacus
