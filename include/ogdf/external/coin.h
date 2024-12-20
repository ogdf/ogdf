/** \file
 * \brief Definition of ogdf::CoinManager
 *
 * \author Markus Chimani, Stephan Beyer
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
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
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */

#pragma once

// IWYU pragma: always_keep

#include <ogdf/basic/Logger.h>

#include <ogdf/lib/abacus/osiinclude.h>
#include <coin/OsiSolverInterface.hpp> // IWYU pragma: export
// not used here but always necessary when using COIN
#include <coin/CoinPackedVector.hpp> // IWYU pragma: export

namespace ogdf {

//! If you use COIN-OR, you should use this class
class OGDF_EXPORT CoinManager {
public:
	//! The OGDF Logger which will determine the log level for a new instance returned by createCorrectOsiSolverInterface.
	/**
	 *  Use updateLogging to propagate changes made to this value to already existing solver instances.
	 */
	static Logger CoinLog;

	//! Get a new solver and set its initial log level according to the level of CoinLog.
	static OsiSolverInterface* createCorrectOsiSolverInterface() {
#ifdef COIN_OSI_CPX
		OsiCpxSolverInterface* ret = new OsiCpxSolverInterface(); // CPLEX
#elif defined(COIN_OSI_GRB)
		OsiGrbSolverInterface* ret = new OsiGrbSolverInterface(); // Gurobi
#elif defined(COIN_OSI_SYM)
		OsiSymSolverInterface* ret = new OsiSymSolverInterface(); // Symphony
		ret->setSymParam(OsiSymVerbosity, -2);
#else // COIN_OSI_CLP
		OsiClpSolverInterface* ret = new OsiClpSolverInterface(); // Coin-OR LP
#endif
		updateLogging(ret);
		return ret;
	}

	//! Update the log level of the CoinMessageHandler associated with \p osi to match the log level of the ogdf::Logger CoinLog.
	static void updateLogging(OsiSolverInterface* osi) {
		if (CoinLog.effectiveStatisticMode()) {
			osi->messageHandler()->setLogLevel(0);
		} else {
			switch (CoinLog.effectiveLogLevel()) {
			//- 0 - none
			//- 1 - minimal
			//- 2 - normal low
			//- 3 - normal high
			//- 4 - verbose
			case Logger::Level::Minor:
				osi->messageHandler()->setLogLevel(4);
				break;
			case Logger::Level::Medium:
				osi->messageHandler()->setLogLevel(3);
				break;
			case Logger::Level::Default:
				osi->messageHandler()->setLogLevel(2);
				break;
			case Logger::Level::High:
				osi->messageHandler()->setLogLevel(1);
				break;
			case Logger::Level::Alarm:
				osi->messageHandler()->setLogLevel(0);
				break;
			case Logger::Level::Force:
				osi->messageHandler()->setLogLevel(0);
				break;
			}
		}
	}
};

}
