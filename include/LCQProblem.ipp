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

#include <cstring>

namespace LCQPow
{

inline ReturnValue LCQProblem::setQ(const double* const Q_new)
{
  if(nV <= 0)
    return LCQPOBJECT_NOT_SETUP;

  Q = new double[nV * nV];
  memcpy(Q, Q_new, (size_t) (nV * nV) * sizeof(double));

  return SUCCESSFUL_RETURN;
}

inline ReturnValue LCQProblem::setG(const double* const g_new)
{
  if(nV == 0)
    return LCQPOBJECT_NOT_SETUP;

  if(Utilities::isNullPtr(g_new))
    return INVALID_OBJECTIVE_LINEAR_TERM;

  g = new double[nV];
  memcpy(g, g_new, (size_t) nV * sizeof(double));

  return SUCCESSFUL_RETURN;
}

inline ReturnValue LCQProblem::setLB(const double* const lb_new)
{
  if(nV == 0)
    return LCQPOBJECT_NOT_SETUP;

  lb = new double[nV];

  if(Utilities::isNotNullPtr(lb_new))
  {
    memcpy(lb, lb_new, (size_t) nV * sizeof(double));
  }
  else
  {
    /* if no lower bounds are specified, set them to -infinity */
    for(int i = 0; i < nV; i++)
      lb[i] = -INFINITY;
  }

  return SUCCESSFUL_RETURN;
}

inline ReturnValue LCQProblem::setLB(int number, double value)
{
  if(nV == 0)
    return LCQPOBJECT_NOT_SETUP;

  if((number >= 0) && (number < nV))
  {
    lb[number] = value;
    return SUCCESSFUL_RETURN;
  }
  else
  {
    return INDEX_OUT_OF_BOUNDS;
  }
}

inline ReturnValue LCQProblem::setUB(const double* const ub_new)
{
  if(nV == 0)
    return LCQPOBJECT_NOT_SETUP;

  ub = new double[nV];

  if(Utilities::isNotNullPtr(ub_new))
  {
    memcpy(ub, ub_new, (size_t) nV * sizeof(double));
  }
  else
  {
    /* if no upper bounds are specified, set them to infinity */
    for(int i = 0; i < nV; ++i)
      ub[i] = INFINITY;
  }

  return SUCCESSFUL_RETURN;
}

inline ReturnValue LCQProblem::setUB(int number, double value)
{
  if(nV == 0)
    return LCQPOBJECT_NOT_SETUP;

  if((number >= 0) && (number < nV))
  {
    ub[number] = value;

    return SUCCESSFUL_RETURN;
  }
  else
  {
    return INDEX_OUT_OF_BOUNDS;
  }
}

inline ReturnValue LCQProblem::setInitialGuess(const double* const _x0, const double* const _y0)
{
  if(nV == 0 || nComp == 0)
    return LCQPOBJECT_NOT_SETUP;

  xk = new double[nV]();

  if(Utilities::isNotNullPtr(_x0))
  {
    memcpy(xk, _x0, (size_t) nV * sizeof(double));
  }

  if(Utilities::isNotNullPtr(_y0))
  {
    // If user passes dual constraints, let us for now assume they have all of the constraint guesses:
    //    1) box, 2) Linear, 3) Complementarity
    int dualGuessLength = nV + nC + 2 * nComp;

    yk = new double[dualGuessLength];
    for(int i = 0; i < dualGuessLength; i++)
      yk[i] = _y0[i];
  }
  else
  {
    yk = (double*) 0;
  }

  return SUCCESSFUL_RETURN;
}

inline void LCQProblem::setOptions(const Options& _options) { options = _options; }
} // namespace LCQPow
