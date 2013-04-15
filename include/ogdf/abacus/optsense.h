/*
 * $Revision: 3386 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-10 14:00:02 +0200 (Mi, 10. Apr 2013) $
 ***************************************************************/

/*!\file
 * \author Matthias Elf
 * \brief sense of optimization.
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

#ifndef ABA__OPTSENSE_H
#define ABA__OPTSENSE_H

#include <ogdf/abacus/abacusroot.h>

namespace abacus {


//! Sense of optimization
/**
 * We can either minimize or maximize the objective function.
 * We encapsulate this information in a class since it is required
 * in various classes, e.g., in the master of the branch-and-bound algorithm and
 * in the linear program.
 */
class  OptSense :  public AbacusRoot  {
public:

	//! The enumeration defining the sense of optimization.
	enum SENSE {
		Min,	//!< Minimization problem.
		Max,	//!< Maximization problem.
		Unknown	//!< Unknown optimization sense, required to recognize uninitialized object.
	};

	//! Initializes the optimization sense to\a s.
	/**
	 * \param s The sense of the optimization. The default value is \a Unknown.
	 */
	OptSense(SENSE s = Unknown) : sense_(s) { }

	//! Output operator for optimization senses.
	/**
	 * The output operator writes the optimization sense on an output stream
	 * in the form <tt>maximize</tt>, <tt>minimize</tt>, or <tt>unknown</tt>.
	 *
	 * \param out The output stream.
	 * \param rhs The sense being output.
	 *
	 * \return The output stream.
	 */
	friend ostream &operator<<(ostream& out, const OptSense &rhs);

	//! Sets the optimization sense to \a s.
	/**
	 * \param s The new sense of the optimization.
	 */
	void sense(SENSE s) { sense_ = s; }


	//! Returns the sense of the optimization.
	SENSE sense() const { return sense_; }


	//! Returns true If it is minimization problem,, false otherwise.
	bool min() const {
		return (sense_ == Min);
	}


	//! Returns true if it is maximization problem,, false otherwise.
	bool max() const {
		return (sense_ == Max);
	}


	//! Returns true if the optimization sense  is unknown,, false otherwise.
	bool unknown() const {
		return (sense_ == Unknown);
	}

private:

	SENSE sense_; //!< The optimization sense.
};


} //namespace abacus

#endif  // OptSense_H
