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

#ifndef LCQPOW_SUBSOLVERBASE_HPP
#define LCQPOW_SUBSOLVERBASE_HPP

#include "Utilities.hpp"

namespace LCQPow
{
class SubsolverBase
{

public:
  /** Get the primal and dual solution.
   *
   * @param x Pointer to the (assumed to be allocated) primal solution vector.
   * @param y Pointer to the (assumed to be allocated) dual solution vector.
   */
  virtual void getSolution(double* x, double* y) = 0;

  /** Abstract method for solving the QP.
   *
   * @param initialSolver A flag indicating whether the call should initialize the sequence.
   * @param iterations A reference to write the number of subsolver iterates to.
   * @param _g The (potentially) updated objective linear component.
   * @param _lbA The (potentially) updated lower bounds of the linear constraints.
   * @param _ubA The (potentially) updated upper bounds of the linear constraints.
   * @param _x0 The primal initial guess. NULL pointer can be passed.
   * @param _y0 The dual initial guess. NULL pointer can be passed.
   * @param _lb The (potentially) updated lower box constraints. NULL pointer can be passed.
   * @param _ub The (potentially) updated upper box constraints. NULL pointer can be passed.
   */
  virtual ReturnValue solve(bool initialSolve, int& iterations, int& exit_flag, const double* const _g,
                            const double* const _lbA, const double* const _ubA, const double* const x0 = 0,
                            const double* const y0 = 0, const double* const _lb = 0, const double* const _ub = 0)
      = 0;
};
} // namespace LCQPow

#endif // LCQPOW_SUBSOLVERBASE_HPP
