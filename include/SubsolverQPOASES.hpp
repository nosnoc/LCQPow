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

#ifndef LCQPOW_SUBSOLVERQPOASES_HPP
#define LCQPOW_SUBSOLVERQPOASES_HPP

#include "SubsolverBase.hpp"
#include <qpOASES.hpp>

namespace LCQPow
{
class SubsolverQPOASES : public SubsolverBase
{

public:
  /** Default constructor. */
  SubsolverQPOASES();

  /** Constructor for dense matrices.
   *
   * @param nV Number of optimization variables.
   * @param nC Number of linear constraints (should include complementarity pairs).
   * @param Q The Hessian matrix in dense format.
   * @param A The linear constraint matrix in dense format (should include the rows of the complementarity selector
   * matrices).
   */
  SubsolverQPOASES(int nV, int nC, double* Q, double* A);

  /** Constructor for sparse matrices.
   *
   * @param nV Number of optimization variables.
   * @param nC Number of linear constraints (should include complementarity pairs).
   * @param Q The Hessian matrix in sparse csc format.
   * @param A The linear constraint matrix in sparse csc format (should include the rows of the complementarity selector
   * matrices).
   */
  SubsolverQPOASES(int nV, int nC, csc* Q, csc* A);

  /** Copy constructor. */
  SubsolverQPOASES(const SubsolverQPOASES& rhs);

  /** Destructor. */
  ~SubsolverQPOASES();

  /** Assignment operator (deep copy). */
  virtual SubsolverQPOASES& operator=(const SubsolverQPOASES& rhs);

  /** Setting the user options. */
  void setOptions(qpOASES::Options options);

  /** Implementation for applying the subsolver to solve the QP.
   *
   * @param initialSolver A flag indicating whether the call should initialize the sequence.
   * @param iterations A reference to write the number of subsolver iterates to.
   * @param _g The (potentially updated) objective linear component.
   * @param _lbA The (potentially updated) lower bounds of the linear constraints.
   * @param _ubA The (potentially updated) upper bounds of the linear constraints.
   * @param _x0 The primal initial guess. NULL pointer can be passed.
   * @param _y0 The dual initial guess. NULL pointer can be passed.
   * @param _lb The (potentially updated) lower box constraints. NULL pointer can be passed.
   * @param _ub The (potentially updated) upper box constraints. NULL pointer can be passed.
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
  void copy(const SubsolverQPOASES& rhs);

  /** Clear the memory. */
  void clear();

private:
  int nV; /**< Number of optimization variables. */
  int nC; /**< Total number of dual variables. */

  bool isSparse = false; /**< A flag storing whether data is given in sparse or dense format. */
  bool useSchur = false; /**< A flag indicating whether to use the Shur Complement method. */

  double* Q = NULL; /**< Hessian matrix in dense format. */
  double* A = NULL; /**< Constraint matrix in dense format (should contain rows of compl. sel. matrices). */

  qpOASES::SymSparseMat* Q_sparse = NULL; /**< Hessian matrix as qpOASES symmetric sparse matrix. */
  qpOASES::SparseMatrix* A_sparse
      = NULL; /**< Constraint matrix as qpOASES sparse matrix (should contain rows of compl. sel. matrices). */

  double* Q_x
      = NULL; /**< Hessian matrix sparse data (required because one cannot copy a symmetric(sprase) qpOASES matrix). */
  int* Q_i
      = NULL; /**< Hessian matrix sparse rows (required because one cannot copy a symmetric(sprase) qpOASES matrix). */
  int* Q_p = NULL; /**< Hessian matrix sparse col pointers (required because one cannot copy a symmetric(sprase) qpOASES
                      mat                           rix). */

  double* A_x = NULL; /**< Constraint matrix sparse data (required because one cannot copy a symmetric(sprase) qpOASES
                         matrix). */
  int* A_i = NULL;    /**< Constraint matrix sparse rows (required because one cannot copy a symmetric(sprase) qpOASES
                         matrix). */
  int* A_p = NULL;    /**< Constraint matrix sparse col pointers (required because one cannot copy a symmetric(sprase)
                         qpOASES matrix). */

  qpOASES::QProblem qp; /**< Store a QP class and call it sequentially (using its hotstart functionality). */
  qpOASES::SQProblemSchur
      qpSchur; /**< Store a Schur Complement QP class and call it sequentially (using its hotstart functionality). */
};
} // namespace LCQPow

#endif // LCQPOW_SUBSOLVERQPOASES_HPP
