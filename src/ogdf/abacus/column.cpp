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

#include <ogdf/abacus/column.h>
#include <ogdf/abacus/global.h>

namespace abacus {


ostream &operator<<(ostream &out, const Column &rhs)
{
	out << "objective function coefficient: " << rhs.obj() << endl
	 << "bounds: " << rhs.lBound_ << " <= x <= " << rhs.uBound_ << endl
	 << "nonzero rows of column :" << endl;

	for (int i = 0; i < rhs.nnz_; i++)
		out << 'r' << rhs.support_[i] << ": " << rhs.coeff_[i] << endl;

	return out;
}


void Column::copy(const Column &col)
{
	SparVec::copy(col);

	obj_    = col.obj_;
	lBound_ = col.lBound_;
	uBound_ = col.uBound_;
}
} //namespace abacus
