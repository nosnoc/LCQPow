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

#include "SubsolverOSQP.hpp"

extern "C"
{
#include <osqp.h>
}

namespace LCQPow
{
SubsolverOSQP::SubsolverOSQP()
{
  data = NULL;
  settings = NULL;
  work = NULL;
}

SubsolverOSQP::SubsolverOSQP(const csc* const _Q, const csc* const _A)
{
  // Store dimensions
  nV = _Q->n;
  nC = _A->m;

  // Copy matrices
  Q = Utilities::copyCSC(_Q, true);
  A = Utilities::copyCSC(_A);
}

SubsolverOSQP::SubsolverOSQP(const SubsolverOSQP& rhs) { copy(rhs); }

SubsolverOSQP::~SubsolverOSQP() { clear(); }

void SubsolverOSQP::clear()
{
  osqp_cleanup(work);

  if(Utilities::isNotNullPtr(data))
  {
    if(Utilities::isNotNullPtr(data->l))
    {
      free(data->l);
      data->l = NULL;
    }

    if(Utilities::isNotNullPtr(data->u))
    {
      free(data->u);
      data->u = NULL;
    }

    if(Utilities::isNotNullPtr(data->q))
    {
      free(data->q);
      data->q = NULL;
    }

    c_free(data);
    data = NULL;
  }

  if(Utilities::isNotNullPtr(Q))
  {
    Utilities::ClearSparseMat(Q);
    Q = NULL;
  }

  if(Utilities::isNotNullPtr(A))
  {
    Utilities::ClearSparseMat(A);
    A = NULL;
  }

  if(Utilities::isNotNullPtr(settings))
  {
    c_free(settings);
    settings = NULL;
  }
}

SubsolverOSQP& SubsolverOSQP::operator=(const SubsolverOSQP& rhs)
{
  if(this != &rhs)
  {
    copy(rhs);
  }

  return *this;
}

void SubsolverOSQP::setOptions(OSQPSettings* _settings)
{
  if(Utilities::isNotNullPtr(settings))
  {
    c_free(settings);
    settings = NULL;
  }

  if(Utilities::isNotNullPtr(_settings))
    settings = copy_settings(_settings);
}

ReturnValue SubsolverOSQP::solve(bool initialSolve, int& iterations, int& exit_flag, const double* const _g,
                                 const double* const _lbA, const double* const _ubA, const double* const x0,
                                 const double* const y0, const double* const _lb, const double* const _ub)
{
  // Make sure that lb and ub are null pointers, as OSQP does not handle box constraints
  if(Utilities::isNotNullPtr(_lb) || Utilities::isNotNullPtr(_ub))
  {
    return ReturnValue::INVALID_OSQP_BOX_CONSTRAINTS;
  }

  // Setup workspace on initial solve
  if(initialSolve)
  {
    double* l = (double*) malloc((size_t) nC * sizeof(double));
    double* u = (double*) malloc((size_t) nC * sizeof(double));
    double* g = (double*) malloc((size_t) nV * sizeof(double));
    memcpy(l, _lbA, (size_t) nC * sizeof(double));
    memcpy(u, _ubA, (size_t) nC * sizeof(double));
    memcpy(g, _g, (size_t) nV * sizeof(c_float));

    data = (OSQPData*) c_malloc(sizeof(OSQPData));
    data->n = nV;
    data->m = nC;
    data->P = Q;
    data->A = A;
    data->q = g;
    data->l = l;
    data->u = u;
    osqp_setup(&work, data, settings);

    if(Utilities::isNotNullPtr(x0))
      if(osqp_warm_start_x(work, x0) != 0)
        return ReturnValue::OSQP_INITIAL_PRIMAL_GUESS_FAILED;

    if(Utilities::isNotNullPtr(y0))
      if(osqp_warm_start_y(work, y0) != 0)
        return ReturnValue::OSQP_INITIAL_DUAL_GUESS_FAILED;
  }
  else
  {
    // Update linear cost and bounds
    osqp_update_lin_cost(work, _g);
    osqp_update_bounds(work, _lbA, _ubA);
  }

  if(work == 0)
  {
    return ReturnValue::OSQP_WORKSPACE_NOT_SET_UP;
  }

  // Solve Problem
  int errorflag = osqp_solve(work);

  // Get number of iterations
  iterations = work->info->iter;
  exit_flag = work->info->status_val;

  // Either pass error
  if(errorflag != 0 || exit_flag <= 0)
    return ReturnValue::SUBPROBLEM_SOLVER_ERROR;

  // Or pass successful return
  return ReturnValue::SUCCESSFUL_RETURN;
}

void SubsolverOSQP::getSolution(double* x, double* y)
{
  OSQPSolution* sol(work->solution);

  if(Utilities::isNotNullPtr(sol->x))
  {
    memcpy(x, sol->x, (size_t) nV * (sizeof(double)));
  }

  // Copy duals with negative sign
  for(int i = 0; i < nC; i++)
  {
    y[i] = -sol->y[i];
  }
}

void SubsolverOSQP::copy(const SubsolverOSQP& rhs)
{
  clear();

  nV = rhs.nV;
  nC = rhs.nC;

  Q = copy_csc_mat(rhs.Q);
  A = copy_csc_mat(rhs.A);

  setOptions(rhs.settings);

  if(Utilities::isNotNullPtr(rhs.data))
  {
    double* l = (double*) malloc((size_t) nC * sizeof(double));
    double* u = (double*) malloc((size_t) nC * sizeof(double));
    double* g = (double*) malloc((size_t) nV * sizeof(double));
    memcpy(l, rhs.data->l, (size_t) nC * sizeof(double));
    memcpy(u, rhs.data->u, (size_t) nC * sizeof(double));
    memcpy(g, rhs.data->q, (size_t) nV * sizeof(c_float));

    data = (OSQPData*) c_malloc(sizeof(OSQPData));
    data->n = nV;
    data->m = nC;
    data->P = Q;
    data->A = A;
    data->q = g;
    data->l = l;
    data->u = u;

    osqp_setup(&work, data, settings);
  }
}
} // namespace LCQPow