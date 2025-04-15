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

#ifndef LCQPOW_OUTPUTSTATISTICS_HPP
#define LCQPOW_OUTPUTSTATISTICS_HPP

#include "Utilities.hpp"
#include <vector>

namespace LCQPow
{

class OutputStatistics
{
public:
  /** Default constructor. */
  OutputStatistics();

  /** Assignment operator.
   *
   * @param rhs The obejct from which to assign.
   */
  OutputStatistics& operator=(const OutputStatistics& rhs);

  /** Resets the statistics. */
  void reset();

  /** Update total iteration counter.
   *
   * @return Success or specifies the invalid argument.
   */
  ReturnValue updateIterTotal(int delta_iter);

  /** Update total outer iteration counter.
   *
   * @return Success or specifies the invalid argument.
   */
  ReturnValue updateIterOuter(int delta_iter);

  /** Update total number of working set changes counter.
   *
   * @return Success or specifies the invalid argument.
   */
  ReturnValue updateSubproblemIter(int delta_iter);

  /** Update rho at solution.
   *
   * @return Success or specifies the invalid argument.
   */
  ReturnValue updateRhoOpt(double _rho);

  /** Update the solution status.
   *
   * @return Success or specifies the invalid argument.
   */
  ReturnValue updateSolutionStatus(AlgorithmStatus _status);

  /** Update the QP solver exit flag.
   *
   * @return Success.
   */
  ReturnValue updateQPSolverExitFlag(int _flag);

  /** Update tracking vectors.
   *
   * @return Success or specifies the invalid argument.
   */
  ReturnValue updateTrackingVectors(double* xStep, int innerIters, int subproblemIters, double stepLength,
                                    double stepSize, double statVals, double objVals, double phiVals, double meritVals,
                                    int nV);

  /** Get the total number of iterations. */
  int getIterTotal() const;

  /** Get the total number of outer iterations. */
  int getIterOuter() const;

  /** Get the total number of subproblem iterations. */
  int getSubproblemIter() const;

  /** Get the penalty parameter at the optimal solution (if found). */
  double getRhoOpt() const;

  /** Get the solution status (if solved it will return the stationarity type). */
  AlgorithmStatus getSolutionStatus() const;

  /** Get the most recent exit flag of the QP solver. */
  int getQPSolverExitFlag() const;

  /** Get values of inner loop iterates.*/
  int* getInnerIters() const;

  /** Get values of inner loop iterates.*/
  std::vector<int> getInnerItersStdVec() const;

  /** Get values of subsolver iterates.*/
  int* getSubproblemIters() const;

  /** Get values of subsolver iterates.*/
  std::vector<int> getSubproblemItersStdVec() const;

  /** Get accumulated number of subsolver iterates.*/
  int* getAccuSubproblemIters() const;

  /** Get accumulated number of subsolver iterates.*/
  std::vector<int> getAccuSubproblemItersStdVec() const;

  /** Get values of alpha.*/
  double* getStepLength() const;

  /** Get values of alpha.*/
  std::vector<double> getStepLengthStdVec() const;

  /** Get values of norm of pk.*/
  double* getStepSize() const;

  /** Get values of norm of pk.*/
  std::vector<double> getStepSizeStdVec() const;

  /** Get Lagrangian's gradient violaton values.*/
  double* getStatVals() const;

  /** Get Lagrangian's gradient violaton values.*/
  std::vector<double> getStatValsStdVec() const;

  /** Get objective function values.*/
  double* getObjVals() const;

  /** Get objective function values.*/
  std::vector<double> getObjValsStdVec() const;

  /** Get penalty function values.*/
  double* getPhiVals() const;

  /** Get penalty function values.*/
  std::vector<double> getPhiValsStdVec() const;

  /** Get merit function values.*/
  double* getMeritVals() const;

  /** Get merit function values.*/
  std::vector<double> getMeritValsStdVec() const;

  /** Get primal iterates. */
  std::vector<std::vector<double>> getxStepsStdVec() const;

private:
  int iterTotal = 0;      /**< Total number of iterations, i.e., total number of inner iterations. */
  int iterOuter = 0;      /**< Total number of outer iterations, i.e., number of penalty updates. */
  int subproblemIter = 0; /**< Total number of subsolver iterations. */
  double rhoOpt = 0.0;    /**< Value of penalty parameter at the final iterate. */
  AlgorithmStatus status
      = PROBLEM_NOT_SOLVED; /**< Status of the solver. This is set to the solution type on success. */
  int qpSolver_exit_flag
      = 0; /**< The exit flag of the most recent QP solved (refer to the respective QP solver docs for meanings). */

  // Tracking vectors
  std::vector<int> innerIters;             /**< Number of inner iterations (accumulated per inner loop). */
  std::vector<int> subproblemIters;        /**< Number of subsolver iterations for each inner loop. */
  std::vector<int> accuSubproblemIters;    /**< Accumulated number of subsolver iterations. */
  std::vector<std::vector<double>> xSteps; /**< Primal iterates. */
  std::vector<double> stepLength;          /**< Track values of alpha. */
  std::vector<double> stepSize;            /**< Track norm of pk. */
  std::vector<double> statVals;            /**< Track values of Lagrangian's gradient violation. */
  std::vector<double> objVals;             /**< Track objective function values. */
  std::vector<double> phiVals;             /**< Track values of complementarity penalty function. */
  std::vector<double> meritVals;           /**< Track merit function values. */
};
} // namespace LCQPow

#endif // LCQPOW_OUTPUTSTATISTICS_HPP