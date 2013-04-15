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

#ifndef ABA__POOLSLOT_INC
#define ABA__POOLSLOT_INC

#include <ogdf/abacus/poolslot.h>
#include <ogdf/abacus/pool.h>
#include <ogdf/abacus/convar.h>

namespace abacus {

template<class BaseType, class CoType>
PoolSlot<BaseType, CoType>::PoolSlot(
	Master *master,
	Pool<BaseType, CoType> *pool,
	BaseType *convar) : master_(master), conVar_(convar), pool_(pool)
{
	version_ = (convar) ? 1 : 0;
}


template<class BaseType, class CoType>
PoolSlot<BaseType, CoType>::~PoolSlot()
{
	if (conVar_ && conVar_->nReferences()) {
		Logger::ifout() << "~PoolSlot(): it is not allowed to destruct objects of class PoolSlot with a constraint/variable with positive reference counter = " << conVar_->nReferences() << ".\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcPoolslot);
	}

	delete conVar_;
}


template<class BaseType, class CoType>
void PoolSlot<BaseType, CoType>::insert(BaseType *convar)
{
	if (conVar_ != 0) {
		Logger::ifout() << "PoolSlot::insert(): insertion failed, the slot is not void\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcPoolslot);
	}

	if (version_ == numeric_limits<unsigned long>::max()) {
		Logger::ifout() << "PoolSlot::insert(): insertion failed, maximum version number ULONG_MAX reached\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcPoolslot);
	}

	conVar_ = convar;
	++version_;
}

} //namespace abacus

#endif   // PoolSlot_INC
