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

/*
 * This is a macro which processes a matlab::data::Array into a double* pointer
 */
#define UNWRAP_ARRAY(name, idx, size) if(!checkDimensionAndTypeDouble(inputs[idx], size, 1, #name, true)) return; \
  bool name##_present = !inputs[idx].isEmpty();                         \
  matlab::data::TypedArray<double> name##_arr(std::move(inputs[idx]));  \
  matlab::data::buffer_ptr_t<double> name##_arr_ptr = name##_arr.release(); \
  auto del_##name = name##_arr_ptr.get_deleter();                       \
  if(name##_present)                                                    \
  {                                                                     \
    name = name##_arr_ptr.release();                                    \
    std::cout << "Getting " << #name << std::endl;                      \
  }

#define MAYBE_CLEANUP_VECTOR(name) if(name##_present) \
  {                                                   \
    del_##name(name);                                 \
    std::cout << "Deleteing " << #name << std::endl;  \
  }
    
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
    matlab::data::buffer_ptr_t<double> Q_arr_ptr = Q_arr.release();
    auto del_Q = Q_arr_ptr.get_deleter();
    Q = Q_arr_ptr.release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER

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
    matlab::data::buffer_ptr_t<double> g_arr_ptr = g_arr.release();
    auto del_g = g_arr_ptr.get_deleter();
    g = g_arr_ptr.release();

    // Get all bounding (optional) vectors
    UNWRAP_ARRAY(lbL,5,nComp);
    UNWRAP_ARRAY(ubL,6,nComp);
    UNWRAP_ARRAY(lbR,7,nComp);
    UNWRAP_ARRAY(ubR,8,nComp);
    UNWRAP_ARRAY(lbA,10,nC);
    UNWRAP_ARRAY(ubA,11,nC);
    UNWRAP_ARRAY(lb,12,nV);
    UNWRAP_ARRAY(ub,13,nV);

    LCQPow::ReturnValue ret = problem->loadLCQP(Q, g, L, R, lbL, ubL, lbR, ubR, A, lbA, ubA, lb, ub, x0, y0);

    // Cleanup arrays
    if(A != nullptr)
    {
      delete[] A;
      std::cout << "Deleting A" << std::endl;
    }
    if(L != nullptr)
    {
      delete[] L;
      std::cout << "Deleting L" << std::endl;
    }
    if(R != nullptr)
    {
      delete[] R;
      std::cout << "Deleting R" << std::endl;
    }
    if(Q != nullptr)
    {
      del_Q(Q);
      std::cout << "Deleting Q" << std::endl;
    }
    // Cleanup vectors
    del_g(g);
    std::cout << "Deleting g" << std::endl;
    MAYBE_CLEANUP_VECTOR(lbL);
    MAYBE_CLEANUP_VECTOR(ubL);
    MAYBE_CLEANUP_VECTOR(lbR);
    MAYBE_CLEANUP_VECTOR(ubR);
    MAYBE_CLEANUP_VECTOR(lbA);
    MAYBE_CLEANUP_VECTOR(ubA);
    MAYBE_CLEANUP_VECTOR(lb);
    MAYBE_CLEANUP_VECTOR(ub);
    
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
