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
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>

/*
 * This is a macro which processes a matlab::data::Array into a double* pointer
 */
#define UNWRAP_ARRAY(NAME, IDX, SIZE) if(!checkDimensionAndTypeDouble(inputs[IDX], SIZE, 1, #NAME, true)) return 1; \
  bool NAME##_present = !inputs[IDX].isEmpty();                         \
  matlab::data::TypedArray<double> NAME##_arr(std::move(inputs[IDX]));  \
  matlab::data::buffer_ptr_t<double> NAME##_arr_ptr = NAME##_arr.release(); \
  auto del_##NAME = NAME##_arr_ptr.get_deleter();                       \
  if(NAME##_present)                                                    \
  {                                                                     \
    NAME = NAME##_arr_ptr.release();                                    \
  }

#define MAYBE_CLEANUP_VECTOR(NAME) if(NAME##_present) \
  {                                                   \
    del_##NAME(NAME);                                 \
  } 

#define GET_SCALAR_DOUBLE_OPTION(NAME, SETTER) if(name == #NAME)        \
  {                                                                     \
    if(!checkDimensionAndTypeDouble(field_arr, 1, 1, "params." #NAME)) return 1; \
    matlab::data::TypedArray<double> field(std::move(field_arr));       \
    options.SETTER(field[0]);                                           \
    continue;                                                           \
  }

#define GET_SCALAR_INT_OPTION(NAME, SETTER) if(name == #NAME)           \
  {                                                                     \
    if(!checkDimensionAndTypeDouble(field_arr, 1, 1, "params." #NAME)) return 1; \
    matlab::data::TypedArray<double> field(std::move(field_arr));       \
    options.SETTER((int) (field[0]));                                   \
    continue;                                                           \
  }

#define GET_SCALAR_BOOL_OPTION(NAME, SETTER) if(name == #NAME)          \
  {                                                                     \
    if(!checkDimensionAndTypeBool(field_arr, 1, 1, "params." #NAME)) return 1; \
    matlab::data::TypedArray<bool> field(std::move(field_arr));         \
    options.SETTER(field[0]);                                           \
    continue;                                                           \
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

    // Initialize values needed for options handling
    LCQPow::Options options;
    double* x0 = nullptr;
    double* y0 = nullptr;

    if(!inputs[14].isEmpty())
    {
      if(!checkTypeStruct(inputs[14], "opts")) return;
      matlab::data::StructArray opts_struct_arr(std::move(inputs[14]));

      handleOptions(problem, opts_struct_arr[0], options, &x0, &y0);
    }

    // Check if we are using sparse matrices
    bool dense = inputs[1].getType() == matlab::data::ArrayType::DOUBLE;
    // Build problem:
    if(dense)
    {
      buildDense(problem, outputs, inputs, x0, y0);
    }
    else
    {
      buildSparse(problem, outputs, inputs, x0, y0);
    }
    // call set options
    // Set options and print them
    std::cout << "AAAAAAAAAAAAAAAAhhhhhh" << std::endl;
    problem->setOptions(options);
  }

  int handleOptions(LCQPow::LCQProblem* problem, matlab::data::Struct opts_struct, LCQPow::Options& options, double** x0, double** y0)
  {
    const std::string params_fieldnames[] = {
      "x0",
      "y0",
      "stationarityTolerance",
      "complementarityTolerance",
      "initialPenaltyParameter",
      "penaltyUpdateFactor",
      "solveZeroPenaltyFirst",
      "maxIterations",
      "maxPenaltyParameter",
      "nDynamicPenalty",
      "etaDynamicPenalty",
      "printLevel",
      "storeSteps",
      "qpSolver",
      "perturbStep",
      "qpOASES_options",
      "OSQP_options"
    };

    int nV = problem->getNV();
    int nC = problem->getNC();
    int nComp = problem->getNComp();
    bool dual_guess_passed = false;
    matlab::data::Array dual_guess_field_arr;
    for(auto name : params_fieldnames)
    {
      matlab::data::Array field_arr;
      try
      {
        field_arr = std::move(opts_struct[name]); // don't copy just move.
      }
      catch(matlab::data::InvalidFieldNameException err) // I hate there isn't a way to check without exception handling.
      {
        continue; // missing, do nothing.
      }

      GET_SCALAR_DOUBLE_OPTION(stationarityTolerance, setStationarityTolerance);
      GET_SCALAR_DOUBLE_OPTION(complementarityTolerance, setComplementarityTolerance);
      GET_SCALAR_DOUBLE_OPTION(initialPenaltyParameter, setInitialPenaltyParameter);
      GET_SCALAR_DOUBLE_OPTION(penaltyUpdateFactor, setPenaltyUpdateFactor);
      GET_SCALAR_DOUBLE_OPTION(etaDynamicPenalty, setEtaDynamicPenalty);
      GET_SCALAR_DOUBLE_OPTION(maxPenaltyParameter, setMaxPenaltyParameter);
      GET_SCALAR_BOOL_OPTION(solveZeroPenaltyFirst, setSolveZeroPenaltyFirst);
      GET_SCALAR_BOOL_OPTION(storeSteps, setStoreSteps);
      GET_SCALAR_INT_OPTION(maxIterations, setMaxIterations);
      GET_SCALAR_INT_OPTION(nDynamicPenatly, setNDynamicPenalty);
      GET_SCALAR_INT_OPTION(printLevel, setPrintLevel);
      GET_SCALAR_INT_OPTION(printLevel, setPrintLevel);
      GET_SCALAR_INT_OPTION(qpSolver, setQPSolver);

      if(name == "x0") {
        if(!checkDimensionAndTypeDouble(field_arr, nV, 1, "params.x0")) return 1;

        matlab::data::TypedArray<double> field(std::move(field_arr));
        *x0 = new double[nV]; // copy because otherwise we need to carry the deleter around :(
        std::copy(field.begin(), field.end(), *x0);
      }

      if(name == "y0") {
        dual_guess_passed = true;
        dual_guess_field_arr = std::move(field_arr);
      }

      if(name == "qpOASES_options") {
        if(!checkTypeStruct(field_arr, "params.qpOASES_options")) return 1;

        matlab::data::StructArray field(std::move(field_arr));
        qpOASES::Options opts;
        setupqpOASESOptions(&opts, field);
        options.setqpOASESOptions(opts);
      }

      if(name == "OSQP_options") {
        if(!checkTypeStruct(field_arr, "params.OSQP_options")) return 1;

        matlab::data::StructArray field(std::move(field_arr));
        OSQPSettings* opts = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
        osqp_set_default_settings(opts);
        setupOSQPOptions(opts, field);
        options.setOSQPOptions(opts);
        c_free(opts);
      }
    }

    if(dual_guess_passed) {            
      LCQPow::QPSolver sol = options.getQPSolver();
      int nDualsIn = 0;
      if(sol < LCQPow::QPSolver::OSQP_SPARSE) {
        nDualsIn = nV + nC + 2*nComp;
        if(!checkDimensionAndTypeDouble(dual_guess_field_arr, nDualsIn, 1, "params.y0")) return 1;
      } else {
        nDualsIn = nC + 2*nComp;
        if(!checkDimensionAndTypeDouble(dual_guess_field_arr, nDualsIn, 1, "params.y0")) return 1;
      }
      matlab::data::TypedArray<double> field(std::move(dual_guess_field_arr));
      *y0 = new double[nDualsIn]; // copy because otherwise we need to carry the deleter around :(
      std::copy(field.begin(), field.end(), *y0);
    }

    return 0;
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

  int buildDense(LCQPow::LCQProblem* problem, matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs, const double* const x0, const double* const y0)
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
    if(!checkDimensionAndTypeDouble(inputs[1], nV, nV, "Q")) return 1;
    matlab::data::TypedArray<double> Q_arr(std::move(inputs[1]));
    matlab::data::buffer_ptr_t<double> Q_arr_ptr = Q_arr.release();
    auto del_Q = Q_arr_ptr.get_deleter();
    Q = Q_arr_ptr.release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER

    // Get L
    if(!checkDimensionAndTypeDouble(inputs[3], nComp, nV, "L")) return 1;
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
    if(!checkDimensionAndTypeDouble(inputs[4], nComp, nV, "R")) return 1;
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
    if(nC>0 && !checkDimensionAndTypeDouble(inputs[9], nC, nV, "A")) return 1;
    matlab::data::TypedArray<double> A_arr(std::move(inputs[9]));
    if(nC>0 && !A_arr.isEmpty())
    {
      matlab::data::buffer_ptr_t<double> A_col_ptr = A_arr.release();
      auto del_A = A_col_ptr.get_deleter();
      A_col = A_col_ptr.release(); // NOTE: WE ARE NOW RESPONSIBLE FOR FREEING THIS POINTER
      A = new double[nC*nV];
      colMajorToRowMajor(A_col, A, nC, nV);
      del_A(A_col);
    }

    // Get vectors
    if(!checkDimensionAndTypeDouble(inputs[2], nV, 1, "g")) return 1;
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
    }
    if(L != nullptr)
    {
      delete[] L;
    }
    if(R != nullptr)
    {
      delete[] R;
    }
    if(Q != nullptr)
    {
      del_Q(Q);
    }
    // Cleanup vectors
    del_g(g);
    MAYBE_CLEANUP_VECTOR(lbL);
    MAYBE_CLEANUP_VECTOR(ubL);
    MAYBE_CLEANUP_VECTOR(lbR);
    MAYBE_CLEANUP_VECTOR(ubR);
    MAYBE_CLEANUP_VECTOR(lbA);
    MAYBE_CLEANUP_VECTOR(ubA);
    MAYBE_CLEANUP_VECTOR(lb);
    MAYBE_CLEANUP_VECTOR(ub);

    return ret;
  }

  int buildSparse(LCQPow::LCQProblem* problem, matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs, const double* const x0, const double* const y0)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    int nV = problem->getNV();
    int nC = problem->getNC();
    int nComp = problem->getNComp();
    if ( !(inputs[1].getType() == matlab::data::ArrayType::SPARSE_DOUBLE) ||
         !(inputs[3].getType() == matlab::data::ArrayType::SPARSE_DOUBLE) ||
         !(inputs[4].getType() == matlab::data::ArrayType::SPARSE_DOUBLE) ||
         (nC > 0 && !(inputs[9].getType() == matlab::data::ArrayType::SPARSE_DOUBLE)))
    {
      std::ostringstream msg;
      msg << "If using the sparse mode, please make sure to provide all matrices in sparse format!\n";
      matlab->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return LCQPow::ReturnValue::DENSE_SPARSE_MISSMATCH;
    }
    matlab::data::SparseArray<double> Q_arr(std::move(inputs[1]));
    matlab::data::SparseArray<double> L_arr(std::move(inputs[3]));
    matlab::data::SparseArray<double> R_arr(std::move(inputs[4]));
    //matlab::data::SparseArray<double> Q_arr(inputs[1]);
    //matlab::data::SparseArray<double> L_arr(inputs[3]);
    //matlab::data::SparseArray<double> R_arr(inputs[4]);
    // Read sparse matrices
    csc* Q = readSparseMatrix(Q_arr, nV, nV);
    csc* L = readSparseMatrix(L_arr, nComp, nV);
    csc* R = readSparseMatrix(R_arr, nComp, nV);
    csc* A = nullptr;
    if (nC > 0) {
      matlab::data::SparseArray<double> A_arr(std::move(inputs[9]));
      //matlab::data::SparseArray<double> A_arr(inputs[9]);
      A = readSparseMatrix(A_arr, nC, nV);
    }

    // Read vectors
    double* g;
    double* lbL = nullptr;
    double* ubL = nullptr;
    double* lbR = nullptr;
    double* ubR = nullptr;
    double* lbA = nullptr;
    double* ubA = nullptr;
    double* lb = nullptr;
    double* ub = nullptr;

    // Get vectors
    if(!checkDimensionAndTypeDouble(inputs[2], nV, 1, "g")) return 1;
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

    std::cout << "pre load Q: " << Q << std::endl;
    LCQPow::ReturnValue ret = problem->loadLCQP(Q, g, L, R, lbL, ubL, lbR, ubR, A, lbA, ubA, lb, ub, x0, y0);
    std::cout << Q << std::endl;
    
    // Cleanup vectors
    std::cout << "cleaning up vectors" << std::endl;
    del_g(g);
    MAYBE_CLEANUP_VECTOR(lbL);
    MAYBE_CLEANUP_VECTOR(ubL);
    MAYBE_CLEANUP_VECTOR(lbR);
    MAYBE_CLEANUP_VECTOR(ubR);
    MAYBE_CLEANUP_VECTOR(lbA);
    MAYBE_CLEANUP_VECTOR(ubA);
    MAYBE_CLEANUP_VECTOR(lb);
    MAYBE_CLEANUP_VECTOR(ub);

    // TODO(@anton): why does this double free?
    // TODO(@anton): this now leaks memory, fix it.
    std::cout << "clearing Q" << std::endl;
    // Clear sparse specific memory
    LCQPow::Utilities::ClearSparseMat(Q);
    std::cout << "clearing L" << std::endl;
    LCQPow::Utilities::ClearSparseMat(L);
    std::cout << "clearing R" << std::endl;
    LCQPow::Utilities::ClearSparseMat(R);
    std::cout << "clearing A" << std::endl;
    if (A != 0) LCQPow::Utilities::ClearSparseMat(A);
    std::cout << "returning buildSparse" << std::endl;
    std::cout << ret << std::endl;
    return ret;
  }
  
  csc* readSparseMatrix(matlab::data::SparseArray<double>& mat, int nRow, int nCol)
  {
    double* x = (double*) malloc(mat.getNumberOfNonZeroElements()*sizeof(double));
    int* i = (int*) malloc(mat.getNumberOfNonZeroElements()*sizeof(int));
    int n = nCol;
    int m = nRow;
    int* p = (int*) malloc((nCol+1)*sizeof(int));
    std::cout << p << std::endl;

    // Get necessary values
    auto mat_it = mat.begin();
    auto mat_it_end = mat.end();
    double* x_it = x;
    int* p_it = p;
    int curr_col = 0;
    p[0] = 0;
    int nc = 0;
    int idx = 0;
    // Iterate over all nonzeros;
    while(mat_it != mat_it_end)
    {
      // get index for current nonzero
      auto mat_idx = mat.getIndex(mat_it);
      int row = mat_idx.first;
      int col = mat_idx.second;
      // If we have moved on to the next column: populate pointers in p
      if(col > curr_col)
      {        
        for(int ii = curr_col+1; ii < col; ii++)
        {
          p[ii] = p[curr_col];
          std::cout << "ii " << ii << " nCoL " << nCol+1 << std::endl;
        }
        p[col] = idx;
        curr_col = col;
      }
      // update row data
      i[idx] = row;
      // copy nonzero into buffer
      *x_it = *mat_it;
      std::cout << "idx " << idx << ", nnz " << mat.getNumberOfNonZeroElements() << std::endl;
      std::cout << "x_it - x " << x_it - x << ", nnz " << mat.getNumberOfNonZeroElements() << std::endl;
      idx++;x_it++;mat_it++; // move pointers
    }
    p[nCol] = idx;

    csc* M = (csc *)malloc(sizeof(csc));

    if(M==nullptr) return nullptr;

		M->m = m;
		M->n = n;
		M->p = p;
		M->i = i;
		M->x = x;
		M->nz = -1;
		M->nzmax = mat.getNumberOfNonZeroElements();

    return M;
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

  bool checkDimensionAndTypeBool(const matlab::data::Array& arr, int m, int n, std::string name, bool allowEmpty = false)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    if (allowEmpty && arr.isEmpty())
    {
      return true;
    }

    if (arr.getType() != matlab::data::ArrayType::LOGICAL)
    {
      std::ostringstream msg;
      msg << "Invalid type: " << name << " must be of type logical.\n";
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

  bool checkTypeStruct(const matlab::data::Array& arr, std::string name)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    if(arr.isEmpty())
    {
      return true;
    }

    if(arr.getType() != matlab::data::ArrayType::STRUCT)
    {
      std::ostringstream msg;
      msg << "Invalid type: " << name << " must be of type struct.\n";
      matlab->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return false;
    }

    matlab::data::ArrayDimensions dims = arr.getDimensions();
    if(dims.size() != 2 || dims[0] != 1 || dims[1] != 1)
    {
      std::ostringstream msg;
      msg << "Invalid dimension: " << name << " must be a scalar struct\n";
      matlab->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return false;
    }

    return true;
  }

  /******************************/ 
  /** QPOASES OPTIONS HANDLING **/
  /******************************/

  /*
   *	has options value
   */
  bool hasOptionsValue(const matlab::data::Struct& options, const std::string& name, double* optionValue)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    matlab::data::Array option;
    try
    {
      option = std::move(options[name]);
    }
    catch(matlab::data::InvalidFieldNameException err) // I hate there isn't a way to check without exception handling.
    {
      std::ostringstream msg;
      msg << "Option struct does not contain entry " << name << ", using default value instead!\n";
      matlab->feval(u"warning", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return qpOASES::BT_FALSE;
    }

    if(option.getType()== matlab::data::ArrayType::DOUBLE && !option.isEmpty() && option.getNumberOfElements() == 1)
    {
      matlab::data::TypedArray<double> val_arr(std::move(option));
      *optionValue = val_arr[0];
      return qpOASES::BT_TRUE;
    }
    else
    {
      std::ostringstream msg;
      msg << "Option " << name << "is not scalar, using default value instead!\n";
      matlab->feval(u"warning", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return qpOASES::BT_FALSE;
    }
  }

  /*
   * TODO(@anton): This is currently broken.
   */
  bool hasOptionsValue(const matlab::data::Struct& options, const std::string& name, char** optionValue)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    matlab::data::Array option;
    try
    {
      option = std::move(options[name]);
    }
    catch(matlab::data::InvalidFieldNameException err) // I hate there isn't a way to check without exception handling.
    {
      std::ostringstream msg;
      msg << "Option struct does not contain entry " << name << ", using default value instead!\n";
      matlab->feval(u"warning", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return qpOASES::BT_FALSE;
    }

    if(option.getType()== matlab::data::ArrayType::MATLAB_STRING && !option.isEmpty() && option.getNumberOfElements() == 1)
    {
      matlab::data::TypedArray<matlab::data::MATLABString> val_arr(std::move(option));
      matlab::data::MATLABString str = val_arr[0];
      // *optionValue = val_arr[0];
      // fix this
      return qpOASES::BT_TRUE;
    }
    else
    {
      std::ostringstream msg;
      msg << "Option " << name << "is not scalar, using default value instead!\n";
      matlab->feval(u"warning", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
      return qpOASES::BT_FALSE;
    }
  }

/*
 *	setup qpOASES options
 */
  void setupqpOASESOptions(qpOASES::Options* options, const matlab::data::StructArray& optionsStructArr)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    double optionValue;
    qpOASES::int_t optionValueInt;
    matlab::data::Struct optionsStruct(std::move(optionsStructArr[0]));

    /* Check for correct number of option entries;
     * may occur, e.g., if user types options.<misspelledName> = <someValue>; */
    if(optionsStructArr.getNumberOfFields() != 31)
    {
      std::ostringstream msg;
      msg << "qpOASES options might be set incorrectly as struct has wrong number of entries!" << std::endl;
      matlab->feval(u"warning", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
    }

    if(hasOptionsValue(optionsStruct, "printLevel", &optionValue))
    {
#ifdef __SUPPRESSANYOUTPUT__
      options->printLevel = qpOASES::PL_NONE;
#else
      optionValueInt = (qpOASES::int_t)optionValue;
      options->printLevel = (REFER_NAMESPACE_QPOASES PrintLevel)optionValueInt;
      if(options->printLevel < qpOASES::PL_DEBUG_ITER)
        options->printLevel = qpOASES::PL_DEBUG_ITER;
      if(options->printLevel > qpOASES::PL_HIGH)
        options->printLevel = qpOASES::PL_HIGH;       
#endif
    }

    if(hasOptionsValue(optionsStruct, "enableRamping", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      options->enableRamping = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "enableFarBounds", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      options->enableFarBounds = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "enableFlippingBounds", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      options->enableFlippingBounds = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "enableRegularisation", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      options->enableRegularisation = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "enableFullLITests", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      options->enableFullLITests = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "enableNZCTests", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      options->enableNZCTests = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "enableDriftCorrection", &optionValue))
      options->enableDriftCorrection = (qpOASES::int_t)optionValue;

    if(hasOptionsValue(optionsStruct, "enableCholeskyRefactorisation", &optionValue))
      options->enableCholeskyRefactorisation = (qpOASES::int_t)optionValue;

    if(hasOptionsValue(optionsStruct, "enableEqualities", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      options->enableEqualities = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "terminationTolerance", &optionValue))
      options->terminationTolerance = optionValue;

    if(hasOptionsValue(optionsStruct, "boundTolerance", &optionValue))
      options->boundTolerance = optionValue;

    if(hasOptionsValue(optionsStruct, "boundRelaxation", &optionValue))
      options->boundRelaxation = optionValue;

    if(hasOptionsValue(optionsStruct, "epsNum", &optionValue))
      options->epsNum = optionValue;

    if(hasOptionsValue(optionsStruct, "epsDen", &optionValue))
      options->epsDen = optionValue;

    if(hasOptionsValue(optionsStruct, "maxPrimalJump", &optionValue))
      options->maxPrimalJump = optionValue;

    if(hasOptionsValue(optionsStruct, "maxDualJump", &optionValue))
      options->maxDualJump = optionValue;

    if(hasOptionsValue(optionsStruct, "initialRamping", &optionValue))
      options->initialRamping = optionValue;

    if(hasOptionsValue(optionsStruct, "finalRamping", &optionValue))
      options->finalRamping = optionValue;

    if(hasOptionsValue(optionsStruct, "initialFarBounds", &optionValue))
      options->initialFarBounds = optionValue;

    if(hasOptionsValue(optionsStruct, "growFarBounds", &optionValue))
      options->growFarBounds = optionValue;

    if(hasOptionsValue(optionsStruct, "initialStatusBounds", &optionValue))
    {
      optionValueInt = (qpOASES::int_t)optionValue;
      if(optionValueInt < -1) 
        optionValueInt = -1;
      if(optionValueInt > 1) 
        optionValueInt = 1;
      options->initialStatusBounds = (REFER_NAMESPACE_QPOASES SubjectToStatus)optionValueInt;
    }

    if(hasOptionsValue(optionsStruct, "epsFlipping", &optionValue))
      options->epsFlipping = optionValue;

    if(hasOptionsValue(optionsStruct, "numRegularisationSteps", &optionValue))
      options->numRegularisationSteps = (qpOASES::int_t)optionValue;

    if(hasOptionsValue(optionsStruct, "epsRegularisation", &optionValue))
      options->epsRegularisation = optionValue;

    if(hasOptionsValue(optionsStruct, "numRefinementSteps", &optionValue))
      options->numRefinementSteps = (qpOASES::int_t)optionValue;

    if(hasOptionsValue(optionsStruct, "epsIterRef", &optionValue))
      options->epsIterRef = optionValue;

    if(hasOptionsValue(optionsStruct, "epsLITests", &optionValue))
      options->epsLITests = optionValue;

    if(hasOptionsValue(optionsStruct, "epsNZCTests", &optionValue))
      options->epsNZCTests = optionValue;
  }


/*
 *  setup OSQP options
 */
  void setupOSQPOptions(OSQPSettings* settings, const matlab::data::StructArray& optionsStructArr)
  {
    std::shared_ptr<matlab::engine::MATLABEngine> matlab = getEngine();
    matlab::data::ArrayFactory factory;
    matlab::data::Struct optionsStruct(std::move(optionsStructArr[0]));
    /* Check for correct number of option entries;
     * may occur, e.g., if user types options.<misspelledName> = <someValue>; */
    if(optionsStructArr.getNumberOfFields() != 22)
    {
      std::ostringstream msg;
      msg << "OSQP options might be set incorrectly as struct has wrong number of entries!" << std::endl;
      matlab->feval(u"warning", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
    }

    double optionValue;
    
    if(hasOptionsValue(optionsStruct, "rho", &optionValue))
      settings->rho = (c_float) optionValue;

    if(hasOptionsValue(optionsStruct, "sigma", &optionValue))
      settings->sigma = (c_float) optionValue;
        
    if(hasOptionsValue(optionsStruct, "scaling", &optionValue))
      settings->scaling = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "adaptive_rho", &optionValue))
      settings->adaptive_rho = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "adaptive_rho_interval", &optionValue))
      settings->adaptive_rho_interval = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "adaptive_rho_tolerance", &optionValue))
      settings->adaptive_rho_tolerance = (c_float) optionValue;

    if(hasOptionsValue(optionsStruct, "adaptive_rho_fraction", &optionValue))
      settings->adaptive_rho_fraction = (c_float) optionValue;

    if(hasOptionsValue(optionsStruct, "max_iter", &optionValue))
      settings->max_iter = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "eps_abs", &optionValue))
      settings->eps_abs = (c_float) optionValue;

    if(hasOptionsValue(optionsStruct, "eps_rel", &optionValue))
      settings->eps_rel = (c_float) optionValue;

    if(hasOptionsValue(optionsStruct, "eps_prim_inf", &optionValue))
      settings->eps_prim_inf = (c_float) optionValue;

    if(hasOptionsValue(optionsStruct, "eps_dual_inf", &optionValue))
      settings->eps_dual_inf = (c_float) optionValue;
        
    if(hasOptionsValue(optionsStruct, "alpha", &optionValue))
      settings->alpha = (c_float) optionValue;
        
    char* optionStr;
    if(hasOptionsValue(optionsStruct, "linsys_solver", &optionStr))
    {
      std::ostringstream msg;
      msg << "Setting OSQP solver through matlab interface currently not supported (will use qdldl)." << std::endl;
      matlab->feval(u"warning", 0, std::vector<matlab::data::Array>({factory.createScalar(msg.str())}));
    }
        
    if(hasOptionsValue(optionsStruct, "delta", &optionValue))
      settings->delta = (c_float) optionValue;
        
    if(hasOptionsValue(optionsStruct, "polish", &optionValue))
      settings->polish = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "polish_refine_iter", &optionValue))
      settings->polish_refine_iter = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "verbose", &optionValue))
      settings->verbose = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "scaled_termination", &optionValue))
      settings->scaled_termination = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "check_termination", &optionValue))
      settings->check_termination = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "warm_start", &optionValue))
      settings->warm_start = (c_int) optionValue;

    if(hasOptionsValue(optionsStruct, "time_limit", &optionValue))
      settings->time_limit = (c_float) optionValue;
  }
};
