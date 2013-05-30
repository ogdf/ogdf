/* $Id: CoinPresolveDupcol.hpp 1372 2011-01-03 23:31:00Z lou $ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CoinPresolveDupcol_H
#define CoinPresolveDupcol_H

#include "CoinPresolveMatrix.hpp"

/*!
  \file
*/

#define	DUPCOL	10

/*! \class dupcol_action
    \brief Detect and remove duplicate columns

    The general technique is to sum the coefficients a_(*,j) of each column.
    Columns with identical sums are duplicates. The obvious problem is that,
    <i>e.g.</i>, [1 0 1 0] and [0 1 0 1] both add to 2. To minimize the
    chances of false positives, the coefficients of each row are multipled by
    a random number r_i, so that we sum r_i*a_ij.

   Candidate columns are checked to confirm they are identical. Where the
   columns have the same objective coefficient, the two are combined. If the
   columns have different objective coefficients, complications ensue. In order
   to remove the duplicate, it must be possible to fix the variable at a bound.
*/

class dupcol_action : public CoinPresolveAction {
  dupcol_action();
  dupcol_action(const dupcol_action& rhs);
  dupcol_action& operator=(const dupcol_action& rhs);

  struct action {
    double thislo;
    double thisup;
    double lastlo;
    double lastup;
    int ithis;
    int ilast;

    double *colels;
    int nincol;
  };

  const int nactions_;
  // actions_ is owned by the class and must be deleted at destruction
  const action *const actions_;

  dupcol_action(int nactions, const action *actions,
		const CoinPresolveAction *next) :
      CoinPresolveAction(next),
      nactions_(nactions),
      actions_(actions) {}

 public:
  const char *name() const;

  static const CoinPresolveAction *presolve(CoinPresolveMatrix *prob,
					 const CoinPresolveAction *next);

  void postsolve(CoinPostsolveMatrix *prob) const;

  ~dupcol_action();

};


/*! \class duprow_action
    \brief Detect and remove duplicate rows

    The algorithm to detect duplicate rows is as outlined for dupcol_action.

    If the feasible interval for one constraint is strictly contained in the
    other, the tighter (contained) constraint is kept. If the feasible
    intervals are disjoint, the problem is infeasible. If the feasible
    intervals overlap, both constraints are kept.

    duprow_action is definitely a work in progress; #postsolve is
    unimplemented.
    This doesn't matter as it uses useless_constraint.
*/

class duprow_action : public CoinPresolveAction {
  struct action {
    int row;
    double lbound;
    double ubound;
  };

  const int nactions_;
  const action *const actions_;

  duprow_action():CoinPresolveAction(NULL),nactions_(0),actions_(NULL) {}
  duprow_action(int nactions,
		      const action *actions,
		      const CoinPresolveAction *next) :
    CoinPresolveAction(next),
    nactions_(nactions), actions_(actions) {}

 public:
  const char *name() const;

  static const CoinPresolveAction *presolve(CoinPresolveMatrix *prob,
					 const CoinPresolveAction *next);

  void postsolve(CoinPostsolveMatrix *prob) const;

  //~duprow_action() { delete[]actions_; }
};

/*! \class gubrow_action
    \brief Detect and remove entries whose sum is known

    If we have an equality row where all entries same then
    For other rows where all entries for that equality row are same
    then we can delete entries and modify rhs
    gubrow_action is definitely a work in progress; #postsolve is
    unimplemented.
*/

class gubrow_action : public CoinPresolveAction {
  struct action {
    int row;
    double lbound;
    double ubound;
  };

  const int nactions_;
  const action *const actions_;

  gubrow_action():CoinPresolveAction(NULL),nactions_(0),actions_(NULL) {}
  gubrow_action(int nactions,
		      const action *actions,
		      const CoinPresolveAction *next) :
    CoinPresolveAction(next),
    nactions_(nactions), actions_(actions) {}

 public:
  const char *name() const;

  static const CoinPresolveAction *presolve(CoinPresolveMatrix *prob,
					 const CoinPresolveAction *next);

  void postsolve(CoinPostsolveMatrix *prob) const;

  //~gubrow_action() { delete[]actions_; }
};

#endif

