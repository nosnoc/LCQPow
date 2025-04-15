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

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "LCQProblem.hpp"
#include <cstdint>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
 public:
  void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;

    // Check inputs are valid
    checkArguments(outputs, inputs);
    // Get self value
    matlab::data::TypedArray<std::uintptr_t> self_array(std::move(matlab->getProperty(inputs[0], u"self")));

    if (self_array.isEmpty())
    {
      return;
    }

    std::uintptr_t self_ptr = reinterpret_cast<std::uintptr_t>((std::uintptr_t) (self_array[0]));

    // Use raw pointer because we will have to handle this as an implicit unique pointer stored in matlab.
    LCQPow::LCQProblem* problem = reinterpret_cast<LCQPow::LCQProblem*>(self_ptr);

    // Run the LCQPow algorithm on setup object.
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    LCQPow::ReturnValue ret = problem->runSolver();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsed_secs = (double)(end - begin).count()/1000.0/1000.0/1000.0;

    // Get dimensions and options for future use
    LCQPow::Options opts = problem->getOptions();
    int nV = problem->getNV();
    int nC = problem->getNC();
    int nComp = problem->getNComp();
    int nDuals = problem->getNumberOfDuals();

    // Setup output arrays
    double* x_opt = new double[nV];
    double* y_opt = new double[nDuals];
    LCQPow::OutputStatistics stats;

    // Get output data
    problem->getPrimalSolution(x_opt);
    problem->getDualSolution(y_opt);
    problem->getOutputStatistics(stats);

    // Generate output array. This does copy because mathworks changed how createArrayFromBuffer
    // works in 2024b but didn't document this change. Therefore we can't use it
    matlab::data::TypedArray<double> x_out = factory.createArray<double>({nV, 1}, x_opt, x_opt + nV);

    // Move (not copy) the primal output to the first output
    outputs[0] = std::move(x_out);

    // Repeat this for y
    if(outputs.size() >= 2)
    {
      matlab::data::TypedArray<double> y_out = factory.createArray<double>({nDuals, 1}, y_opt, y_opt + nV);
      outputs[1] = std::move(y_out);
    }

    // Process output statistics if necessary.
    if(outputs.size() >= 3)
    {
      matlab::data::Array stats_struct;
      if(opts.getStoreSteps())
      {
        stats_struct = std::move(factory.createStructArray({1,1}, {"iters_total", "iters_outer", "iters_subproblem", "rho_opt",
              "elapsed_time", "exit_flag", "solution_type", "qp_exit_flag",
              "innerIters", "xSteps", "accumulatedSubproblemIters", "stepLength", "stepSize",
              "statVals", "objVals", "phiVals", "meritVals", "subproblemIters"}));
      }
      else
      {
        stats_struct = std::move(factory.createStructArray({1,1}, {"iters_total", "iters_outer", "iters_subproblem", "rho_opt",
              "elapsed_time", "exit_flag", "solution_type", "qp_exit_flag"}));
      }

      // Add always present entries
      int iters_total = stats.getIterTotal();
      stats_struct[0]["iters_total"] = factory.createScalar(iters_total);
      stats_struct[0]["iters_outer"] = factory.createScalar(stats.getIterOuter());
      stats_struct[0]["iters_subproblem"] = factory.createScalar(stats.getSubproblemIter());
      stats_struct[0]["rho_opt"] = factory.createScalar(stats.getRhoOpt());
      stats_struct[0]["elapsed_time"] = factory.createScalar(elapsed_secs);
      stats_struct[0]["exit_flag"] = factory.createScalar((int) ret);
      stats_struct[0]["solution_type"] = factory.createScalar((int) stats.getSolutionStatus());
      stats_struct[0]["qp_exit_flag"] = factory.createScalar(stats.getQPSolverExitFlag());

      // Add iter values:
      if(opts.getStoreSteps() && stats.getIterTotal() > 0)
      {
        // xSteps
        // TODO(@anton) perhaps this is one too many copies? but annoyingly nothing much can be done I think.
        std::vector<std::vector<double>> xSteps_vec = stats.getxStepsStdVec();
        std::vector<double> xSteps;
        xSteps.reserve(iters_total*nV);
        for(auto step : xSteps_vec)
        {
          xSteps.insert(xSteps.end(), step.begin(), step.end());
        }
        stats_struct[0]["xSteps"] = factory.createArray({iters_total, nV}, xSteps.begin(), xSteps.end());
        // inner iters
        std::vector<int> innerIters = stats.getInnerItersStdVec();
        stats_struct[0]["innerIters"] = factory.createArray({iters_total, 1}, innerIters.begin(), innerIters.end());
        // accuSubproblemItters
        std::vector<int> accuSubproblemIters = stats.getInnerItersStdVec();
        stats_struct[0]["accumulatedSubproblemIters"] = factory.createArray({iters_total, 1}, accuSubproblemIters.begin(), accuSubproblemIters.end());
        // stepLength
        std::vector<double> stepLength = stats.getStepLengthStdVec();
        stats_struct[0]["stepLength"] = factory.createArray({iters_total, 1}, stepLength.begin(), stepLength.end());
        // stepSize
        std::vector<double> stepSize = stats.getStepSizeStdVec();
        stats_struct[0]["stepSize"] = factory.createArray({iters_total, 1}, stepSize.begin(), stepSize.end());
        // statVals
        std::vector<double> statVals = stats.getStatValsStdVec();
        stats_struct[0]["statVals"] = factory.createArray({iters_total, 1}, statVals.begin(), statVals.end());
        // objVals
        std::vector<double> objVals = stats.getObjValsStdVec();
        stats_struct[0]["objVals"] = factory.createArray({iters_total, 1}, objVals.begin(), objVals.end());
        // phiVals
        std::vector<double> phiVals = stats.getPhiValsStdVec();
        stats_struct[0]["phiVals"] = factory.createArray({iters_total, 1}, phiVals.begin(), phiVals.end());
        // meritVals
        std::vector<double> meritVals = stats.getMeritValsStdVec();
        stats_struct[0]["meritVals"] = factory.createArray({iters_total, 1}, meritVals.begin(), meritVals.end());
      }
      outputs[2] = std::move(stats_struct);
    }
  }

  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    if(outputs.size() > 3)
    {
      matlab->feval(u"error", 0,
                    std::vector<matlab::data::Array>({ factory.createScalar("A maximum of 3 outputs are provided") }));

    }
    if(inputs.size() > 1)
    {
      matlab->feval(u"error", 0,
                    std::vector<matlab::data::Array>({ factory.createScalar("Incorrect number of inputs") }));

    }


  }
};
