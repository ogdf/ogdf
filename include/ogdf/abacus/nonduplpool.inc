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

#ifndef ABA__NONDUPLPOOL_INC
#define ABA__NONDUPLPOOL_INC

#include <ogdf/abacus/nonduplpool.h>

namespace abacus {


template<class BaseType, class CoType>
PoolSlot<BaseType, CoType> * NonDuplPool<BaseType, CoType>::insert(
	BaseType *cv)
{

	PoolSlot<BaseType, CoType>* slot;

	slot = present(cv);
	if (slot == 0) {
		slot = StandardPool<BaseType, CoType>::insert(cv);
		if (slot)
			hash_.insert(cv->hashKey(), slot);
	}
	else {
		delete cv;
		nDuplications_++;
	}
	return slot;
}


template<class BaseType, class CoType>
PoolSlot<BaseType, CoType> *NonDuplPool<BaseType, CoType>::present(
	BaseType *cv)
{
	int key = cv->hashKey();

	PoolSlot<BaseType, CoType> **cand = hash_.initializeIteration(key);
	while(cand) {
		if (cv->equal((*cand)->conVar())) {
			return *cand;
		}
		cand = hash_.next(key);
	}
	return 0;
}

template<class BaseType, class CoType>
const PoolSlot<BaseType, CoType> *NonDuplPool<BaseType, CoType>::present(
	const BaseType *cv) const
{
	int key = cv->hashKey();

	PoolSlot<BaseType, CoType> *const *cand = hash_.initializeIteration(key);
	while(cand) {
		if (cv->equal((*cand)->conVar())) {
			return *cand;
		}
		cand = hash_.next(key);
	}
	return 0;
}

template<class BaseType, class CoType>
int NonDuplPool<BaseType, CoType>::softDeleteConVar(PoolSlot<BaseType, CoType> *slot)
{
	int key = slot->conVar()->hashKey();

	if (Pool<BaseType, CoType>::softDeleteConVar(slot) == 0) {
		if (hash_.remove(key, slot)) {
			// Commented out in ABACUS 3.2
			//Logger::ifout() << "NonDuplPool::softDeleteCon(): slot not found in hash table.\n";
			//OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcNonDuplPool);
		}
		return 0;
	}
	return 1;
}


template<class BaseType, class CoType>
void NonDuplPool<BaseType, CoType>::hardDeleteConVar(
	PoolSlot<BaseType, CoType> *slot)
{
	if (hash_.remove(slot->conVar()->hashKey(), slot)) {
		// Commented out in ABACUS 3.2
		//Logger::ifout() << "NonDuplPool::hardDeleteConVar(): constraint not found in hash table.\n";
		//OGDF_THROW_PARAM(AlgorithmFailureException, ogdf::afcNonDuplPool);
	}
	Pool<BaseType, CoType>::hardDeleteConVar(slot);
}

} // namespace abacus

#endif   // NonDuplPool_INC
