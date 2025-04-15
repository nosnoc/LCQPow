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

#include "MessageHandler.hpp"

namespace LCQPow
{

ReturnValue MessageHandler::PrintMessage(ReturnValue ret, MessageType type)
{

  // Print the message type
  switch(type)
  {
  case MESSAGE:
    break;

  case WARNING:
    printf("WARNING: ");
    break;

  case ERROR:
    printf("ERROR: ");
    break;

  default:
    break;
  }

  switch(ret)
  {
  case SUCCESSFUL_RETURN:
    break;

  case NOT_YET_IMPLEMENTED:
    printf("This method has not yet been implemented.\n");
    break;

  case LCQPOBJECT_NOT_SETUP:
    printf("The LCQP object has not been set up correctly.\n");
    break;

  case INDEX_OUT_OF_BOUNDS:
    printf("Index out of bounds.\n");
    break;

  case SUBPROBLEM_SOLVER_ERROR:
    printf("The subproblem solver produced an error.\n");
    break;

  case UNABLE_TO_READ_FILE:
    printf("Unable to read file.\n");
    break;

  case MAX_ITERATIONS_REACHED:
    printf("Maximum number of iterations reached.\n");
    break;

  case MAX_PENALTY_REACHED:
    printf("Maxium penalty value reached.\n");
    break;

  case INITIAL_SUBPROBLEM_FAILED:
    printf("Failed to solve initial QP.\n");
    break;

  case INVALID_ARGUMENT:
    printf("Invalid argument passed.\n");
    break;

  case INVALID_NUMBER_OF_OPTIM_VARS:
    printf("Invalid optimization variable dimension passed (required to be > 0).\n");
    break;

  case INVALID_NUMBER_OF_COMP_VARS:
    printf("Invalid complementarity dimension passed (required to be > 0).\n");
    break;

  case INVALID_NUMBER_OF_CONSTRAINT_VARS:
    printf("Invalid number of optimization variables passed (required to be >= 0).\n");
    break;

  case INVALID_QPSOLVER:
    printf("Invalid QPSolver passed.\n");
    break;

  case INVALID_COMPLEMENTARITY_TOLERANCE:
    printf("Ignoring invalid complementarity tolerance.\n");
    break;

  case INVALID_INITIAL_PENALTY_VALUE:
    printf("Invalid argument passed (initial penalty value).\n");
    break;

  case INVALID_PENALTY_UPDATE_VALUE:
    printf("Ignoring invalid penalty update value.\n");
    break;

  case INVALID_MAX_ITERATIONS_VALUE:
    printf("Ignoring invalid number of maximum iterations.\n");
    break;

  case INVALID_MAX_RHO_VALUE:
    printf("Ignoring invalid number of maximum penalty value.\n");
    break;

  case DENSE_SPARSE_MISSMATCH:
    printf("The solver was initialized with dense (sparse) matrices but a sparse (dense) method was chosen.\n");
    break;

  case INVALID_STATIONARITY_TOLERANCE:
    printf("Ignoring invalid stationarity tolerance.\n");
    break;

  case INVALID_INDEX_POINTER:
    printf("Invalid index pointer passed in csc format.\n");
    break;

  case INVALID_INDEX_ARRAY:
    printf("Invalid index array passed in csc format.\n");
    break;

  case INVALID_OSQP_BOX_CONSTRAINTS:
    printf("Invalid constraints passed to OSQP solver: This solver does not handle box constraints, please pass them "
           "through linear constraints.\n");
    break;

  case INVALID_TOTAL_ITER_COUNT:
    printf("Invalid total number of iterations delta passed to output statistics (must be non-negative integer).\n");
    break;

  case INVALID_TOTAL_OUTER_ITER:
    printf(
        "Invalid total number of outer iterations delta passed to output statistics (must be non-negative integer).\n");
    break;

  case IVALID_SUBPROBLEM_ITER:
    printf("Invalid total number of subproblem solver iterates delta passed to output statistics (must be non-negative "
           "integer).\n");
    break;

  case INVALID_RHO_OPT:
    printf("Invalid rho value at solution passed to output statistics (must be positive double).\n");
    break;

  case INVALID_PRINT_LEVEL_VALUE:
    printf("Ignoring invalid integer to be parsed to print level passed (must be in range of enum).\n");
    break;

  case INVALID_OBJECTIVE_LINEAR_TERM:
    printf("Invalid objective linear term passed (must be a double array of length n).\n");
    break;

  case INVALID_CONSTRAINT_MATRIX:
    printf("Invalid constraint matrix passed (matrix was null pointer but number of constraints is positive).\n");
    break;

  case INVALID_COMPLEMENTARITY_MATRIX:
    printf("Invalid complementarity matrix passed (can not be null pointer).\n");
    break;

  case INVALID_ETA_VALUE:
    printf("Invalid etaDynamicPenalty value, which describes the fraction of loss required for complementarity "
           "progress (must be in (0,1)).");
    break;

  case OSQP_INITIAL_PRIMAL_GUESS_FAILED:
    printf("OSQP failed to use the primal initial guess.\n");
    break;

  case OSQP_INITIAL_DUAL_GUESS_FAILED:
    printf("OSQP failed to use the dual initial guess.\n");
    break;

  case INVALID_LOWER_COMPLEMENTARITY_BOUND:
    printf("Lower complementarity bound must be bounded below.\n");
    break;

  case FAILED_SYM_COMPLEMENTARITY_MATRIX:
    printf("Failed to compute the symmetric complementarity matrix C.\n");
    break;

  case FAILED_SWITCH_TO_SPARSE:
    printf("Failed to switch to sparse mode (a to be created sparse matrix was nullpointer).\n");
    break;

  case FAILED_SWITCH_TO_DENSE:
    printf("Failed to switch to dense mode (an array to be created was nullpointer).\n");
    break;

  case OSQP_WORKSPACE_NOT_SET_UP:
    printf("OSQP Workspace is not set up (please check for OSQP errors).\n");
    break;
  }

  fflush(stdout);

  return ret;
}

AlgorithmStatus MessageHandler::PrintSolution(AlgorithmStatus algoStat)
{

  if(algoStat != PROBLEM_NOT_SOLVED)
    printf("\n\n#################################\n");

  switch(algoStat)
  {
  case PROBLEM_NOT_SOLVED:
    printf("The LCQP has not been solved.\n");
    break;

  case W_STATIONARY_SOLUTION:
    printf("## W-Stationary solution found ##\n");
    break;

  case C_STATIONARY_SOLUTION:
    printf("## C-Stationary solution found ##\n");
    break;

  case M_STATIONARY_SOLUTION:
    printf("## M-Stationary solution found ##\n");
    break;

  case S_STATIONARY_SOLUTION:
    printf("## S-Stationary solution found ##\n");
    break;
  }

  if(algoStat != PROBLEM_NOT_SOLVED)
    printf("#################################\n\n");

  return algoStat;
}
} // namespace LCQPow
