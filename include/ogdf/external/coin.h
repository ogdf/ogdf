/** \file
 * \brief Handles connection to the COIN library, by offering
 * helper classes.
 *
 * If you use Coin, you need to include this file.
 *
 * \todo Currently, there is only a single implementation of the
 * CoinCallback-class declared herein (necc. for userdefined cuts).
 * This implementation is CPLEX specific.
 * -- with current coin, it might not even work anymore!
 *
 * \author Markus Chimani
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

#include <ogdf/basic/basic.h>

#include <coin/OsiSolverInterface.hpp>
#include <coin/CoinPackedVector.hpp>

namespace ogdf {

class OGDF_EXPORT CoinCallbacks {
	friend class OGDF_EXPORT CoinManager;
public:
	enum CallbackType { CT_Cut = 1, CT_Heuristic = 2, CT_Incumbent = 4, CT_Branch  = 8 };
	enum CutReturn { CR_Error, CR_SolutionValid, CR_AddCuts, CR_DontAddCuts, CR_NoCutsFound };
	enum HeuristicReturn { HR_Error, HR_Ignore, HR_Update };
	enum IncumbentReturn { IR_Error, IR_Ignore, IR_Update };
#if 0
	enum BranchReturn { BR_Error, ... };
#endif
	virtual CutReturn cutCallback(const double /* objValue */, const double* /* fracSolution */, OsiCuts* /* addThese */) { return CR_Error; }
	virtual HeuristicReturn heuristicCallback(double& /* objValue */, double* /* solution */) { return HR_Error; }
	virtual IncumbentReturn incumbentCallback(const double /* objValue */, const double* /* solution */) { return IR_Error; }
#if 0
	virtual BranchReturn branchCallback() { return BR_Error; };
#endif

	virtual ~CoinCallbacks() {}
private:
	bool registerCallbacks(OsiSolverInterface* _posi, int callbackTypes);
};

class OGDF_EXPORT CoinManager {
public:
	static OsiSolverInterface* createCorrectOsiSolverInterface();
	static OsiSolverInterface* createCorrectOsiSolverInterface(CoinCallbacks* ccc, int callbackTypes) {
		OsiSolverInterface* posi = createCorrectOsiSolverInterface();
		if(ccc->registerCallbacks(posi, callbackTypes))
			return posi;
		delete posi;
		return nullptr;
	}
	static void logging(OsiSolverInterface* osi, bool logMe);
};

}
