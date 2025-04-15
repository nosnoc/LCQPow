/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2022 by Jonas Hall et al.
 *
 *	LCQPow is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPow is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef LCQPOW_SUBSOLVEROSQP_HPP
#define LCQPOW_SUBSOLVEROSQP_HPP

#include "SubsolverBase.hpp"

extern "C"
{
#include <osqp.h>
}

namespace LCQPow
{
class SubsolverOSQP : public SubsolverBase
{

public:
  /** Default constructor. */
  SubsolverOSQP();

  /** Constructor for sparse matrices.
   *
   * @param Q The Hessian matrix in sparse csc format.
   * @param A The linear constraint matrix in sparse csc format (should include the rows of the complementarity selector
   * matrices).
   */
  SubsolverOSQP(const csc* const _Q, const csc* const _A);

  /** Copy constructor. */
  SubsolverOSQP(const SubsolverOSQP& rhs);

  /** Destructor. */
  ~SubsolverOSQP();

  /** Clear memory. */
  void clear();

  /** Assignment operator (deep copy). */
  virtual SubsolverOSQP& operator=(const SubsolverOSQP& rhs);

  /** Set OSQP settings. */
  void setOptions(OSQPSettings* settings);

  /** Implementation for applying the subsolver to solve the QP.
   *
   * @param initialSolver A flag indicating whether the call should initialize the sequence.
   * @param iterations A reference to write the number of subsolver iterates to.
   * @param _g The (potentially) updated objective linear component.
   * @param _lbA The (potentially) updated lower bounds of the linear constraints.
   * @param _ubA The (potentially) updated upper bounds of the linear constraints.
   * @param _x0 The primal initial guess. NULL pointer can be passed.
   * @param _y0 The dual initial guess. NULL pointer can be passed.
   * @param _lb This entry is ignored in this solver (only required in declaration due to inflexibility of abstract
   * classes).
   * @param _ub This entry is ignored in this solver (only required in declaration due to inflexibility of abstract
   * classes).
   */
  ReturnValue solve(bool initialSolve, int& iterations, int& exit_flag, const double* const _g,
                    const double* const _lbA, const double* const _ubA, const double* const x0 = 0,
                    const double* const y0 = 0, const double* const _lb = 0, const double* const _ub = 0);

  /** Get the primal and dual solution.
   *
   * @param x Pointer to the (assumed to be allocated) primal solution vector.
   * @param y Pointer to the (assumed to be allocated) dual solution vector.
   */
  void getSolution(double* x, double* y);

protected:
  /** Copies all members from given rhs object. */
  void copy(const SubsolverOSQP& rhs);

private:
  int nV; /**< Number of optimization variables. */
  int nC; /**< Number of constraints. */

  OSQPWorkspace* work = NULL;    /**< OSQP workspace. */
  OSQPSettings* settings = NULL; /**< OSQP settings. */
  OSQPData* data = NULL;         /**< OSQP data. */

  csc* Q = NULL; /**< Hessian matrix in csc format (must be upper triagonal). */
  csc* A = NULL; /**< Constraint matrix in csc format (should contain rows of compl. sel. matrices). */
};
} // namespace LCQPow

#endif // LCQPOW_SUBSOLVEROSQP_HPP