/** \file
 * \brief Implementation of mathematical constants, functions.
 *
 * \author Carsten Gutwenger
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


#include <ogdf/basic/Math.h>


namespace ogdf {


	const double Math::pi     = 3.14159265358979323846;
	const double Math::pi_2   = 1.57079632679489661923;
	const double Math::log_of_4 = log(4.0);

	int Math::binomial(int n, int k)
	{
		if(k>n/2) k = n-k;
		if(k == 0) return 1;
		int r = n;
		for(int i = 2; i<=k; ++i)
			r = (r * (n+1-i))/i;
		return r;
	}

	double Math::binomial_d(int n, int k)
	{
		if(k>n/2) k = n-k;
		if(k == 0) return 1.0;
		double r = n;
		for(int i = 2; i<=k; ++i)
			r = (r * (n+1-i))/i;
		return r;
	}

	int Math::factorial(int n)
	{
		return (int) std::tgamma(n+1);
	}

	double Math::factorial_d(int n)
	{
		return std::tgamma(n+1);
	}

} // namespace ogdf
