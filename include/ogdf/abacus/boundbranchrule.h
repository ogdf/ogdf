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
 * \brief This class implements a branching rule for modifying the lower and
 * the upper bound of a variable.
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

#ifndef ABA__BOUNDBRANCHRULE_H
#define ABA__BOUNDBRANCHRULE_H

#include <ogdf/abacus/branchrule.h>
#include <ogdf/abacus/sub.h>

namespace abacus {

//! Implements a branching rule for modifying the lower and the upper bound of a variable.
class  BoundBranchRule :  public BranchRule  {
public:

	//! Creates a bound branch rule for given branching \a variable, lower bound \a lBound and upper bound \a uBound.
	/**
	 * \param master   A pointer to the corresponding master of the optimization.
	 * \param variable The branching variable.
	 * \param lBound   The lower bound of the branching variable.
	 * \param uBound   The upper bound of the branching variable.
	 */
	BoundBranchRule(
		Master *master,
		int variable,
		double lBound,
		double uBound)
		: BranchRule(master), variable_(variable), lBound_(lBound), uBound_(uBound) { }


	virtual ~BoundBranchRule() { }

	//! Output operator for bound branch rules.
	/**
	 * Writes the branching variable together with its lower and upper bound to an output stream.
	 *
	 * \param out The output stream.
	 * \param rhs The branch rule being output.
	 *
	 * \return A reference to the output stream.
	 */
	friend ostream &operator<<(ostream &out, const BoundBranchRule &rhs);


	//! Modifies a subproblem by changing the lower and the upper bound of the branching variable.
	/**
	 * \param sub The subproblem being modified.
	 *
	 * \return 0 If the subproblem is successfully modified.
	 * \return 1 If the modification causes a contradiction.
	 */
	virtual int extract(Sub *sub) {
		if (sub->fsVarStat(variable_)->fixedOrSet())
			return 1;

		sub->lBound(variable_, lBound_);
		sub->uBound(variable_, uBound_);

		return 0;
	}


	//! Pverloaded to modify directly the linear programming relaxation.
	/**
	 * This required to evaluate the quality of a branching rule.
	 */
	virtual void extract(LpSub *lp) {
		oldLpLBound_ = lp->lBound(variable_);
		oldLpUBound_ = lp->uBound(variable_);

		lp->changeLBound(variable_, lBound_);
		lp->changeUBound(variable_, uBound_);
	}


	virtual void unExtract(LpSub *lp) {
		lp->changeLBound(variable_, oldLpLBound_);
		lp->changeUBound(variable_, oldLpUBound_);
	}


	//! Returns the number of the branching variable.
	int variable() const {
		return variable_;
	}

	//! Returns the lower bound of the branching variable.
	double lBound() const {
		return lBound_;
	}

	//! Returns the upper bound of the branching variable.
	double uBound() const {
		return uBound_;
	}

private:
	int    variable_;     //!< The branching variable.
	double lBound_;       //!< The lower bound of the branching variable.
	double uBound_;       //!< The upper bound of the branching variable.
	double oldLpLBound_;
	double oldLpUBound_;
};


inline ostream &operator<<(ostream &out, const BoundBranchRule &rhs)
{
	return out << rhs.lBound_ << " <= x" << rhs.variable_ << " <= " << rhs.uBound_;
}

} //namespace abacus

#endif  // !BoundBranchRule_H
