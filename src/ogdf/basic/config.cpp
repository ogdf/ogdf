/*
 * $Revision: 2963 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-11-05 14:17:50 +0100 (Mo, 05. Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Implementation of basic configuration utilities.
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

#include <ogdf/internal/basic/config.h>


namespace ogdf {

	static string system_str[Configuration::sysSTOP+1] = {
		"unknown", "Windows", "Unix/linux", "Apple OSX", "STOP"
	};

	static string lpsolver_str[Configuration::lpsSTOP+1] = {
		"N/A", "COIN-OR LP (Clp)", "Symphony", "CPLEX", "Gurobi", "STOP"
	};

	static string mm_str[Configuration::mmSTOP+1] = {
		"pool allocator (thread-safe)", "pool allocator (not thread-safe)", "malloc", "STOP"
	};


	const string &Configuration::toString(System sys)
	{
		return system_str[(sys < sysSTOP) ? sys : sysSTOP];
	}

	const string &Configuration::toString(LPSolver lps)
	{
		return lpsolver_str[(lps < lpsSTOP) ? lps : lpsSTOP];
	}

	const string &Configuration::toString(MemoryManager mm)
	{
		return mm_str[(mm < mmSTOP) ? mm : mmSTOP];
	}


	Configuration::System Configuration::whichSystem()
	{
#ifdef OGDF_SYSTEM_WINDOWS
		return sysWindows;
#elif defined(OGDF_SYSTEM_OSX)
		return sysOSX;
#elif defined(OGDF_SYSTEM_UNIX)
		return sysUnix;
#else
		return sysUnknown
#endif
	}


	bool Configuration::haveLPSolver()
	{
#ifdef OGDF_LP_SOLVER
		return true;
#else
		return false;
#endif
	}


	Configuration::LPSolver Configuration::whichLPSolver()
	{
#if defined(COIN_OSI_CLP)
		return lpsClp;
#elif defined(COIN_OSI_SYM)
		return lpsSymphony;
#elif defined(COIN_OSI_CPX)
		return lpsCPLEX;
#elif defined(COIN_OSI_GRB)
		return lpsGurobi;
#else
		return lpsNone;
#endif
	}


	bool Configuration::haveCoin()
	{
#ifdef USE_COIN
		return true;
#else
		return false;
#endif
	}


	bool Configuration::haveAbacus()
	{
#ifdef USE_ABACUS
		return true;
#else
		return false;
#endif
	}


	Configuration::MemoryManager Configuration::whichMemoryManager()
	{
#ifdef OGDF_MEMORY_POOL_TS
		return mmPoolTS;
#elif defined(OGDF_MEMORY_POOL_NTS)
		return mmPoolNTS;
#elif defined(OGDF_MEMORY_MALLOC_TS)
		return mmMalloc;
#else
		return sysSTOP
#endif
	}

}
