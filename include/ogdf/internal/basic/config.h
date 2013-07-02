/*
 * $Revision: 3533 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-06-03 18:22:41 +0200 (Mo, 03. Jun 2013) $
 ***************************************************************/

/** \file
 * \brief Basic configuration file
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_CONFIG_H
#define OGDF_CONFIG_H

#include <ogdf/internal/version.h>

#include <ogdf/internal/config_autogen.h>

// define minimal MS runtime version for mingw32
#if defined(__MINGW32__) && !defined(__MINGW64__)
#ifndef __MSVCRT_VERSION__
#define __MSVCRT_VERSION__ 0x0700
#endif
#endif


// common stdlib includes
#include <iostream>
#include <string>

// generally used <iostream> members
using std::ios;
using std::istream;
using std::ostream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::swap;

// generally used <string> members
using std::string;


//---------------------------------------------------------
// detection of the system
//---------------------------------------------------------

#if defined(unix) || defined(__unix__) || defined(__unix) || defined(_AIX) || defined(__APPLE__)
#define OGDF_SYSTEM_UNIX
#endif

#if defined(__WIN32__) || defined(_WIN32) || defined(__NT__)
#define OGDF_SYSTEM_WINDOWS
#endif

// Note: Apple OS X machines will be both OGDF_SYSTEM_UNIX and OGDF_SYSTEM_OSX
#if defined(__APPLE__)
#define OGDF_SYSTEM_OSX
#endif


// COIN and ABACUS
#if defined(USE_COIN)

#define OGDF_LP_SOLVER
#define USE_ABACUS
#define ABACUS_LP_OSI
#define OSI_CLP
#define OSI_SYM

#if defined(COIN_OSI_CPX) && !defined(OSI_CPX)
#define OSI_CPX
#endif
#if defined(COIN_OSI_GRB) && !defined(OSI_GRB)
#define OSI_GRB
#endif

#if !defined(COIN_OSI_CPX) && !defined(COIN_OSI_GRB) && !defined(COIN_OSI_SYM) && !defined(COIN_OSI_CLP)
#error "Compiler-flag USE_COIN requires an additional COIN_OSI_xxx-flag to choose the default LP solver backend."
#endif

#endif



//---------------------------------------------------------
// C++ standard
//---------------------------------------------------------

#if __cplusplus >= 201103
#define OGDF_HAVE_CPP11

#elif defined(_MSC_VER)
#if _MSC_VER >= 1600
#define OGDF_HAVE_CPP11
#endif

#elif defined(__GNUC__)
#ifdef __GXX_EXPERIMENTAL_CXX0X__
#define OGDF_HAVE_CPP11
#endif

#endif


#if defined(__CYGWIN__) || defined(__APPLE__) || defined(__sparc__)
#define OGDF_NO_COMPILER_TLS
#elif defined(__GNUC__)
#if __GNUC__ < 4
#define OGDF_NO_COMPILER_TLS
#endif
#endif


//---------------------------------------------------------
// macros for compiling OGDF as DLL
//---------------------------------------------------------

#ifdef OGDF_SYSTEM_WINDOWS
#ifdef OGDF_DLL

#ifdef OGDF_INSTALL
#define OGDF_EXPORT __declspec(dllexport)
#else
#define OGDF_EXPORT __declspec(dllimport)
#endif

#else
#define OGDF_EXPORT
#endif

#else
#define OGDF_EXPORT
#endif


//---------------------------------------------------------
// compiler adaptions
//---------------------------------------------------------

#ifdef _MSC_VER

#ifdef OGDF_DLL
// disable useless warnings
// missing dll-interface

// warning C4251: 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'
#pragma warning(disable : 4251)
// warning C4275: non – DLL-interface classkey 'identifier' used as base for DLL-interface classkey 'identifier'
#pragma warning(disable : 4275)
#endif

// warning C4355: 'this' : used in base member initializer list
#pragma warning (disable : 4355)

#endif


//---------------------------------------------------------
// memory manager
//
// OGDF_MEMORY_POOL_TS:
//   buffered-pool allocator per thread pool (thread-safe)
//
// OGDF_MEMORY_POOL_NTS:
//   pool allocator (not thread-safe)
//
// OGDF_MEMORY_MALLOC_TS:
//   just using malloc/free (thread-safe)
//
// default (nothing defined): OGDF_MEMORY_POOL_TS
//---------------------------------------------------------

// By default, we use the thread-safe pool allocator
#if !defined(OGDF_MEMORY_POOL_NTS) && !defined(OGDF_MEMORY_MALLOC_TS) && !defined(OGDF_MEMORY_POOL_TS)
#define OGDF_MEMORY_POOL_TS
#endif


namespace ogdf {

	class OGDF_EXPORT Configuration {
	public:
		//! Specifies the operating system for which OGDF has been configured/built.
		enum System {
			sysUnknown,	//!< not known (inproper configuration)
			sysWindows,	//!< Windows
			sysUnix,	//!< Unix/Linux
			sysOSX,		//!< Apple OSX
			sysSTOP
		};

		//! Specifies the LP-solver used by OGDF.
		enum LPSolver {
			lpsNone,		//!< no LP-solver available
			lpsClp,			//!< COIN-OR LP-solver (Clp)
			lpsSymphony,	//!< Symphony
			lpsCPLEX,		//!< CPLEX (commercial)
			lpsGurobi,		//!< Gurobi (commercial)
			lpsSTOP
		};

		//! Specifies the memory-manager used by OGDF.
		enum MemoryManager {
			mmPoolTS,	//!< thread-safe pool allocator
			mmPoolNTS,	//!< non-thread-safe pool allocator
			mmMalloc,	//!< malloc/free allocator
			mmSTOP
		};

		//! Returns the operating system for which OGDF has been configured.
		static System whichSystem();

		//! Returns whether OGDF has been configured with LP-solver support.
		static bool haveLPSolver();

		//! Returns the LP-solver used by OGDF.
		static LPSolver whichLPSolver();


		//! Returns whether OGDF has been configured with COIN support.
		/**
		 * COIN is used as LP-solver by some OGDF algorithms. If OGDF is configured
		 * without COIN support, this functionality is not available.
		 */
		static bool haveCoin();

		//! Returns whether OGDF has been configured with ABACUS support.
		/**
		 * ABACUS is used as branch-and-cut-solver by some OGDF algorithms.
		 * If OGDF is configured without ABACUS support, this functionality is not available.
		 */
		static bool haveAbacus();

		//! Returns the memory-manager used by OGDF.
		static MemoryManager whichMemoryManager();


		//! Converts \a sys to a (readable) string.
		static const string &toString(System sys);

		//! Converts \a lps to a (readable) string.
		static const string &toString(LPSolver lps);

		//! Converts \a mm to a (readable) string.
		static const string &toString(MemoryManager mm);
	};


	//! Output operator for Configuration::System (uses Configuration::toString(Configuration::System)).
	inline ostream &operator<<(ostream &os, Configuration::System sys)
	{
		os << Configuration::toString(sys);
		return os;
	}

	//! Output operator for Configuration::LPSolver (uses Configuration::toString(Configuration::LPSolver)).
	inline ostream &operator<<(ostream &os, Configuration::LPSolver lps)
	{
		os << Configuration::toString(lps);
		return os;
	}

	//! Output operator for Configuration::MemoryManager (uses Configuration::toString(Configuration::MemoryManager)).
	inline ostream &operator<<(ostream &os, Configuration::MemoryManager mm)
	{
		os << Configuration::toString(mm);
		return os;
	}
}


#endif
