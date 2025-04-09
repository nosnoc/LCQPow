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

    // Get the uintptr_t array from self member property
    matlab::data::TypedArray<std::uintptr_t> self_array(std::move(matlab->getProperty(inputs[0], u"self")));

    // Check if we have no value, which _should_ never happen but fail out anyway.
    if(self_array.isEmpty())
    {
      // TODO(@anton) probably need to throw an exception here.
      return;
    }

    std::uintptr_t self_ptr = reinterpret_cast<std::uintptr_t>((std::uintptr_t) (self_array[0]));

    LCQPow::LCQProblem* problem = reinterpret_cast<LCQPow::LCQProblem*>(self_ptr);
    // Check if self ptr is null and fail out immediately.
    if(problem == nullptr)
    {
      // TODO(@anton) probably need to throw an exception here.
      return;
    }

    // TODO(@anton) support also sparse
    boolean dense = true;
    // Build problem:
    if dense
    {
      buildDense(problem, outputs, inputs);
    }
    else
    {
      buildSparse(problem, outputs, inputs);
    }

    // Parse options

    // call create problem

    // call set options
  }

  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    if(outputs.size() > 0)
    {
      matlab->feval(u"error", 0,
                    std::vector<matlab::data::Array>({ factory.createScalar("no outputs returned") }));

    }
    if(inputs.size() < 8 || inputs.size() > 14)
    {
      matlab->feval(u"error", 0,
                    std::vector<matlab::data::Array>({ factory.createScalar("Incorrect number of inputs") }));

    }
  }

  void buildDense(LCQPow::LCQProblem* problem, matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
  {
    // PARSE INPUTS
    // Get Q. Q is assumed symmetric PSD so we don't care that matlab uses column major order.
    matlab::data::TypedArray<double> Q_arr(std::move(inputs[1]));
    double* Q = Q_arr.release().release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER

    matlab::data::TypedArray<double> L_arr(std::move(inputs[2]));
    double* L = L_arr.release().release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER

    matlab::data::TypedArray<double> R_arr(std::move(inputs[3]));
    double* R = R_arr.release().release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER

    matlab::data::TypedArray<double> A_arr(std::move(inputs[8]));
    double* A = A_arr.release().release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER


  }

  void buildSparse(LCQPow::LCQProblem* problem, matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
  {
    // PARSE INPUTS
    // Get Q:

  }

};
