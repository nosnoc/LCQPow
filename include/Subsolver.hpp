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

#ifndef LCQPOW_SUBSOLVER_HPP
#define LCQPOW_SUBSOLVER_HPP

#include "SubsolverOSQP.hpp"
#include "SubsolverQPOASES.hpp"

extern "C"
{
#include <osqp.h>
}

namespace LCQPow
{

class Subsolver
{

public:
  /** Default constructor. */
  Subsolver();

  /** Constructor for dense matrices (qpOASES).
   *
   * @param nV The number of optimization variables.
   * @param nC The number of linear constraints (should include the complementarity pairs).
   * @param Q The Hessian matrix in dense format.
   * @param A The linear constraint matrix (should include the rows of the complementarity selector matrices).
   */
  Subsolver(int nV, int nC, double* Q, double* A);

  /** Constructor for sparse matrices (qpOASES/OSQP).
   *
   * @param nV The number of optimization variables.
   * @param nC The number of linear constraints (should include the complementarity pairs).
   * @param Q The Hessian matrix in sparse csc format.
   * @param A The linear constraint matrix in sparse csc format (should include the rows of the complementarity selector
   * matrices).
   * @param qpSolver The QP subproblem solver to be used.
   */
  Subsolver(int nV, int nC, csc* Q, csc* A, QPSolver qpSolver);

  /** Copy constructor. */
  Subsolver(const Subsolver& rhs);

  /** Assignment operator (deep copy). */
  virtual Subsolver& operator=(const Subsolver& rhs);

  /** Write solution to x and y. */
  void getSolution(double* x, double* y);

  /** Abstract method for solving the QP. */
  ReturnValue solve(bool initialSolve, int& iterations, int& exit_flag, const double* g, const double* lbA,
                    const double* ubA, const double* x0 = 0, const double* y0 = 0, const double* lb = 0,
                    const double* ub = 0);

  /** Setting the user options. */
  void setOptions(qpOASES::Options& options);

  /** Setting the user options. */
  void setOptions(OSQPSettings* settings);

protected:
  /** Copies all members from given rhs object. */
  void copy(const Subsolver& rhs);

private:
  // The solver type
  QPSolver qpSolver; /**< Inidicating which qpSolver to use. */

  // The different solvers
  SubsolverQPOASES solverQPOASES; /**< When using qpOASES. */
  SubsolverOSQP solverOSQP;       /**< When using OSQP. */
};
} // namespace LCQPow

#endif // LCQPOW_SUBSOLVER_HPP