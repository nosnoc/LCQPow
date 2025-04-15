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

    // Assume the inputs are integral, this is generally not true but, shrug, it is fine.
    int nV = (int) inputs[1][0];
    int nC = (int) inputs[2][0];
    int nComp = (int) inputs[3][0];

    // Use raw pointer because we will have to handle this as an implicit unique pointer stored in matlab.
    LCQPow::LCQProblem* problem = new LCQPow::LCQProblem(nV,nC,nComp);

    // Reinterpret the pointer to an unsigned integer
    std::uintptr_t self = reinterpret_cast<std::uintptr_t>(problem);

    // Set the self property the current object
    matlab->setProperty(inputs[0], u"self", factory.createScalar(self));
  }

  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    if(outputs.size() > 0)
    {
      matlab->feval(u"error", 0,
                    std::vector<matlab::data::Array>({ factory.createScalar("no outputs returned") }));

    }
    if(inputs.size() != 4)
    {
      matlab->feval(u"error", 0,
                    std::vector<matlab::data::Array>({ factory.createScalar("Incorrect number of inputs") }));

    }
    if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE || inputs[1].getNumberOfElements() != 1)
    {
      matlab->feval(u"error", 0,
                       std::vector<matlab::data::Array>({ factory.createScalar("nV must be a scalar double") }));
    }
    if (inputs[2].getType() != matlab::data::ArrayType::DOUBLE || inputs[2].getNumberOfElements() != 1)
    {
      matlab->feval(u"error", 0,
                       std::vector<matlab::data::Array>({ factory.createScalar("nC must be a scalar double") }));
    }
    if (inputs[3].getType() != matlab::data::ArrayType::DOUBLE || inputs[3].getNumberOfElements() != 1)
    {
      matlab->feval(u"error", 0,
                       std::vector<matlab::data::Array>({ factory.createScalar("nComp must be a scalar double") }));
    }
  }
};
