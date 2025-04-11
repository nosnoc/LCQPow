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
#include <iostream>
#include <sstream>
#include <string>
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
    bool dense = true;
    // Build problem:
    if(dense)
    {
      buildDense(problem, outputs, inputs, nullptr, nullptr);
    }
    else
    {
      buildSparse(problem, outputs, inputs, nullptr, nullptr);
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
    if(inputs.size() != 15)
    {
      matlab->feval(u"error", 0,
                    std::vector<matlab::data::Array>({ factory.createScalar("Incorrect number of inputs") }));

    }
  }

  void buildDense(LCQPow::LCQProblem* problem, matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs, const double* const x0, const double* const y0)
  {
    double* Q = nullptr;
    double* g = nullptr;
    double* L = nullptr;
    double* R = nullptr;
    double* lbL = nullptr;
    double* ubL = nullptr;
    double* lbR = nullptr;
    double* ubR = nullptr;
    double* lb = nullptr;
    double* ub = nullptr;
    double* A = nullptr;
    double* lbA = nullptr;
    double* ubA = nullptr;

    // Temporary matrices to keep the column major formats
    double* L_col = nullptr;
    double* R_col = nullptr;
    double* A_col = nullptr;

    int nV = problem->getNV();
    int nC = problem->getNC();
    int nComp = problem->getNComp();

    // PARSE INPUTS
    // Get Q. Q is assumed symmetric PSD so we don't care that matlab uses column major order.
    if(!checkDimensionAndTypeDouble(inputs[1], nV, nV, "Q")) return;
    matlab::data::TypedArray<double> Q_arr(std::move(inputs[1]));
    Q = Q_arr.release().release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER

    // Get L
    if(!checkDimensionAndTypeDouble(inputs[3], nComp, nV, "L")) return;
    matlab::data::TypedArray<double> L_arr(std::move(inputs[3]));
    matlab::data::buffer_ptr_t<double> L_col_ptr = L_arr.release();
    auto del_L = L_col_ptr.get_deleter();
    L_col = L_arr.release().release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER
    if(L_col != nullptr)
    {
      L = new double[nComp*nV];
      colMajorToRowMajor(L_col, L, nComp, nV);
      del_L(L_col);
    }

    // Get R
    if(!checkDimensionAndTypeDouble(inputs[4], nComp, nV, "R")) return;
    matlab::data::TypedArray<double> R_arr(std::move(inputs[4]));
    matlab::data::buffer_ptr_t<double> R_col_ptr = R_arr.release();
    auto del_R = R_col_ptr.get_deleter();
    R_col = R_col_ptr.release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER
    if(R_col != nullptr)
    {
      R = new double[nComp*nV];
      colMajorToRowMajor(R_col, R, nComp, nV);
      del_R(R_col);
    }

    // Get A
    if(!checkDimensionAndTypeDouble(inputs[9], nC, nV, "A")) return;
    matlab::data::TypedArray<double> A_arr(std::move(inputs[9]));
    if(!A_arr.isEmpty())
    {
      matlab::data::buffer_ptr_t<double> A_col_ptr = A_arr.release();
      auto del_A = A_col_ptr.get_deleter();
      A_col = A_col_ptr.release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER
      A = new double[nC*nV];
      colMajorToRowMajor(A_col, A, nC, nV);
      del_A(A_col);
    }

    // Get vectors
    if(!checkDimensionAndTypeDouble(inputs[2], nV, 1, "g")) return;
    matlab::data::TypedArray<double> g_arr(std::move(inputs[2]));
    g = g_arr.release().release();

    if(!checkDimensionAndTypeDouble(inputs[5], nComp, 1, "lbL", true)) return;
    if(!inputs[5].isEmpty())
    {
      matlab::data::TypedArray<double> lbL_arr(std::move(inputs[5]));
      lbL = lbL_arr.release().release();
    }

    if(!checkDimensionAndTypeDouble(inputs[6], nComp, 1, "ubL", true)) return;
    if(!inputs[6].isEmpty())
    {
      matlab::data::TypedArray<double> ubL_arr(std::move(inputs[6]));
      ubL = ubL_arr.release().release();
    }

    if(!checkDimensionAndTypeDouble(inputs[7], nComp, 1, "lbR", true)) return;
    if(!inputs[7].isEmpty())
    {
      matlab::data::TypedArray<double> lbR_arr(std::move(inputs[7]));
      lbR = lbR_arr.release().release();
    }

    if(!checkDimensionAndTypeDouble(inputs[8], nComp, 1, "ubR", true)) return;
    if(!inputs[8].isEmpty())
    {
      matlab::data::TypedArray<double> ubR_arr(std::move(inputs[8]));
      ubR = ubR_arr.release().release();
    }

    if(!checkDimensionAndTypeDouble(inputs[10], nC, 1, "lbA", true)) return;
    if(!inputs[10].isEmpty())
    {
      matlab::data::TypedArray<double> lbA_arr(std::move(inputs[10]));
      lbA = lbA_arr.release().release();
    }

    if(!checkDimensionAndTypeDouble(inputs[11], nC, 1, "ubA", true)) return;
    if(!inputs[11].isEmpty())
    {
      matlab::data::TypedArray<double> ubA_arr(std::move(inputs[11]));
      ubA = ubA_arr.release().release();
    }

    if(!checkDimensionAndTypeDouble(inputs[12], nV, 1, "lb", true)) return;
    if(!inputs[12].isEmpty())
    {
      matlab::data::TypedArray<double> lb_arr(std::move(inputs[12]));
      lb = lb_arr.release().release();
    }

    if(!checkDimensionAndTypeDouble(inputs[13], nV, 1, "ub", true)) return;
    if(!inputs[13].isEmpty())
    {
      matlab::data::TypedArray<double> ub_arr(std::move(inputs[13]));
      ub = ub_arr.release().release();
    }

    LCQPow::ReturnValue ret = problem->loadLCQP(Q, g, L, R, lbL, ubL, lbR, ubR, A, lbA, ubA, lb, ub, x0, y0);
  }

  void buildSparse(LCQPow::LCQProblem* problem, matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs, const double* const x0, const double* const y0)
  {
    // PARSE INPUTS
    // Get Q:

  }

  void colMajorToRowMajor(double* col_maj, double* row_maj, int m, int n)
  {
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        row_maj[i*n + j] = col_maj[j*m + i];
  }

  bool checkDimensionAndTypeDouble(const matlab::data::Array& arr, int m, int n, std::string name, bool allowEmpty = false)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    if (allowEmpty && arr.isEmpty())
    {
      return true;
    }

    if (arr.getType() != matlab::data::ArrayType::DOUBLE)
    {
      std::ostringstream msg;
      msg << "Invalid type: " << name << " must be of type double.\n";
      matlab->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return false;
    }

    matlab::data::ArrayDimensions dims = arr.getDimensions();
    if (dims.size() != 2 || dims[0] != m || dims[1] != n) {
      std::ostringstream msg;
      msg << "Invalid dimension: " << name << " must be " << m << " x " << n << "\n";
      matlab->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return false;
    }

    return true;
  }

};
