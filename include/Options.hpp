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

#ifndef LCQPOW_OPTIONS_HPP
#define LCQPOW_OPTIONS_HPP

#include "Utilities.hpp"

namespace LCQPow
{

class Options
{

public:
  /** Default constructor. */
  Options();

  /** Copy constructor (deep copy).
   *
   * @param rhs The object to be copied.
   */
  Options(const Options& rhs);

  /** Destructor. */
  ~Options();

  /** Assignment operator.
   *
   * @param rhs The obejct from which to assign.
   */
  Options& operator=(const Options& rhs);

  /** Sets all options to default values. */
  void setToDefault();

  /** Get stationarity tolerance. */
  double getStationarityTolerance();

  /** Set stationarity tolerance. */
  ReturnValue setStationarityTolerance(double val);

  /** Get complementarity tolerance. */
  double getComplementarityTolerance();

  /** Set complementarity tolerance. */
  ReturnValue setComplementarityTolerance(double val);

  /** Get initial penalty parameter. */
  double getInitialPenaltyParameter();

  /** Set complementarity tolerance. */
  ReturnValue setInitialPenaltyParameter(double val);

  /** Get penalty parameter update factor. */
  double getPenaltyUpdateFactor();

  /** Set penalty parameter update factor. */
  ReturnValue setPenaltyUpdateFactor(double val);

  /** Get whether to solve for (complement.) unconstrained global minumum first. */
  bool getSolveZeroPenaltyFirst();

  /** Set whether to solve for (complement.) unconstrained global minumum first. */
  ReturnValue setSolveZeroPenaltyFirst(bool val);

  /** Get whether to perform step perturbation. */
  bool getPerturbStep();

  /** Set whether to perform step perturbation. */
  ReturnValue setPerturbStep(bool val);

  /** Get maximum number of iterations. */
  int getMaxIterations();

  /** Set maximum number of iterations. */
  ReturnValue setMaxIterations(int val);

  /** Get maximum penalty value. */
  double getMaxPenaltyParameter();

  /** Set maximum penalty value. */
  ReturnValue setMaxPenaltyParameter(double val);

  /** Get number of previous complementarity values to check against (Leyffer strategy). */
  int getNDynamicPenalty();

  /** Set number of previous complementarity values to check against (Leyffer strategy). */
  ReturnValue setNDynamicPenalty(int val);

  /** Get fraction required for complementarity loss (Leyffer strategy). */
  double getEtaDynamicPenalty();

  /** Get fraction required for complementarity loss (Leyffer strategy). */
  ReturnValue setEtaDynamicPenalty(double val);

  /** Get print level. */
  PrintLevel getPrintLevel();

  /** Set print level. */
  ReturnValue setPrintLevel(PrintLevel val);

  /** Set print level (using an integer). */
  ReturnValue setPrintLevel(int val);

  /** Get whether to store information per iterate. */
  bool getStoreSteps();

  /** Set whether to store information per iterate. */
  ReturnValue setStoreSteps(bool val);

  /** Get QP solver. */
  QPSolver getQPSolver();

  /** Set print level. */
  ReturnValue setQPSolver(QPSolver val);

  /** Set print level (using an integer). */
  ReturnValue setQPSolver(int val);

  /** Pass options for qpOASES. */
  ReturnValue setqpOASESOptions(const qpOASES::Options& _options);

  /** Get options for qpOASES. */
  qpOASES::Options& getqpOASESOptions();

  /** Get options for OSQP. */
  ReturnValue setOSQPOptions(OSQPSettings* _options);

  /** Pass options for OSQP. */
  OSQPSettings* getOSQPOptions();

protected:
  void copy(const Options& rhs); /**< Copy each property. */

  double stationarityTolerance;    /**< Tolerance for 1-Norm of stationarity violation. */
  double complementarityTolerance; /**< Complementarity tolerance. */
  double initialPenaltyParameter;  /**< Start value for complementarity penalty term. */
  double penaltyUpdateFactor;      /**< Factor for updating penaltised complementarity term. */

  bool solveZeroPenaltyFirst; /**< Flag indicating whether first QP should ignore penalization. */

  bool perturbStep; /**< Flag whether to perform step perturbation. */

  int maxIterations;          /**< Maximum number of iterations to be performed. */
  double maxPenaltyParameter; /**< Maximum penalty value. */

  int nDynamicPenalty; /**< Number of previous iterates to compare complementarity loss (only enabled if positive). */
  double etaDynamicPenalty; /**< Parameter describing fraction of required complementarity loss. */

  PrintLevel printLevel; /**< Print level. */

  bool storeSteps; /**< Whether to store detailed information for each iterate (time consuming). */

  QPSolver qpSolver;              /**< The QP solver to be used. */
  qpOASES::Options qpOASES_opts;  /**< qpOASES options. */
  OSQPSettings* OSQP_opts = NULL; /**< OSQP options. */
};
} // namespace LCQPow

#endif // LCQPOW_UTILITIES_HPP
