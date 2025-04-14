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


#ifndef LCQPOW_UTILITIES_HPP
#define LCQPOW_UTILITIES_HPP

extern "C" {
    #include <osqp.h>
}

#include <qpOASES.hpp>

namespace LCQPow {

    /**
     *  Various return values passed throughout all classes.
     */
    enum ReturnValue {
        // Special values
        NOT_YET_IMPLEMENTED = -1,                       /**< Not yet implemented (internal use only). */
        SUCCESSFUL_RETURN = 0,						    /**< Successful return. */

        // Invalid arguments
        INVALID_ARGUMENT = 100,                         /**< Generic invalid argument. */
        INVALID_PENALTY_UPDATE_VALUE = 101,             /**< Invalid penalty update value. Must be > 1. */
        INVALID_COMPLEMENTARITY_TOLERANCE = 102,        /**< Invalid complementarity tolerance. Must be no smaller than machine precision. */
        INVALID_INITIAL_PENALTY_VALUE = 103,            /**< Invalid initial penalty parameter. Must be positive. */
        INVALID_MAX_ITERATIONS_VALUE = 104,             /**< Invalid number of maximal outer iterations. Must be a positive integer. */
        INVALID_STATIONARITY_TOLERANCE = 105,           /**< Invalid stationarity tolerance. Must be no smaller than machine precision. */
        INVALID_NUMBER_OF_OPTIM_VARS = 106,             /**< Invalid number of optimization variables. Must be a positive integer. */
        INVALID_NUMBER_OF_COMP_VARS = 107,              /**< Invalid number of complementarity constraints. Must be a positive integer. */
        INVALID_NUMBER_OF_CONSTRAINT_VARS = 108,        /**< Invalid number of linear constraints. Must be a non-negative integer. */
        INVALID_QPSOLVER = 109,                         /**< Invalid QPSolver passed. */
        INVALID_OSQP_BOX_CONSTRAINTS = 110,             /**< Invalid constraints passed to OSQP solver: This solver does not handle box constraints, please pass them through linear constraints. */
        INVALID_TOTAL_ITER_COUNT = 111,                 /**< Invalid total number of iterations delta passed to output statistics (must be non-negative integer). */
        INVALID_TOTAL_OUTER_ITER = 112,                 /**< Invalid total number of outer iterations delta passed to output statistics (must be non-negative integer). */
        IVALID_SUBPROBLEM_ITER = 113,                   /**< Invalid total number of subproblem solver iterates delta passed to output statistics (must be non-negative integer). */
        INVALID_RHO_OPT = 114,                          /**< Invalid rho value at solution passed to output statistics. (must be positive double). */
        INVALID_PRINT_LEVEL_VALUE = 115,                /**< Invalid integer to be parsed to print level passed (must be in range of enum). */
        INVALID_OBJECTIVE_LINEAR_TERM = 116,            /**< Invalid objective linear term passed (must be a double array of length n). */
        INVALID_CONSTRAINT_MATRIX = 117,                /**< Invalid constraint matrix passed (matrix was null pointer but number of constraints is positive). */
        INVALID_COMPLEMENTARITY_MATRIX = 118,           /**< Invalid complementarity matrix passed (can not be null pointer). */
        INVALID_ETA_VALUE = 119,                        /**< Invalid etaDynamicPenalty value, which describes the fraction of loss required for complementarity progress (must be in (0,1)). */
        INVALID_LOWER_COMPLEMENTARITY_BOUND = 120,      /**< Lower complementarity bound must be bounded below. */
        INVALID_MAX_RHO_VALUE = 121,                    /**< Invalid maximal penalty value. Must be a positive double. */

        // Algorithmic errors
        MAX_ITERATIONS_REACHED = 200,                   /**< Maximum number of iterations reached. */
        MAX_PENALTY_REACHED = 201,                      /**< Maxium penalty value reached. */
        INITIAL_SUBPROBLEM_FAILED = 202,                /**< Failed to solve the initial QP. */
        SUBPROBLEM_SOLVER_ERROR = 203,                  /**< An error occured in the subproblem solver. */
        FAILED_SYM_COMPLEMENTARITY_MATRIX = 204,        /**< Failed to compute the symmetric complementarity matrix C. */
        FAILED_SWITCH_TO_SPARSE = 205,                  /**< Failed to switch to sparse mode (a to be created sparse matrix was nullpointer). */
        FAILED_SWITCH_TO_DENSE = 206,                   /**< Failed to switch to dense mode (an array to be created was nullpointer). */
        OSQP_WORKSPACE_NOT_SET_UP = 207,                /**< OSQP Workspace is not set up. */
        OSQP_INITIAL_PRIMAL_GUESS_FAILED = 208,         /**< OSQP failed to use the primal initial guess. */
        OSQP_INITIAL_DUAL_GUESS_FAILED = 209,           /**< OSQP failed to use the dual initial guess. */

        // Generic errors
        LCQPOBJECT_NOT_SETUP = 300,                     /**< Constructor has not been called. */
        INDEX_OUT_OF_BOUNDS = 301,                      /**< Index out of bounds. */
        UNABLE_TO_READ_FILE = 302,                      /**< Unable to read a file. */

        // Sparse matrices
        INVALID_INDEX_POINTER = 400,                    /**< Invalid index pointer for a csc matrix. */
        INVALID_INDEX_ARRAY = 401,                      /**< Invalid index array for a csc matrix. */
        DENSE_SPARSE_MISSMATCH = 402,                   /**< Solver was initialized with dense (sparse) matrices but a sparse (dense) method was chosen. */
    };


    /**
     *  Various message types.
     */
    enum MessageType {
        MESSAGE = 0,                                    /**< Standard output. */
        WARNING = 1,                                    /**< Program will continue. */
        ERROR = 2                                       /**< Program should abort. */
    };


    /**
     *  Various algorithm status.
     */
    enum AlgorithmStatus {
        PROBLEM_NOT_SOLVED = 0,                         /**< The problem was not solved. */
        W_STATIONARY_SOLUTION = 1,                      /**< The solution corresponds to a weakly stationary point. */
        C_STATIONARY_SOLUTION = 2,                      /**< The solution corresponds to a Clarke stationary point. */
        M_STATIONARY_SOLUTION = 3,                      /**< The solution corresponds to a Mordukhovich stationary point. */
        S_STATIONARY_SOLUTION = 4                       /**< The solution corresponds to a strongly stationary point. */
    };


    /**
     *  Various print levels.
     */
    enum PrintLevel {
        NONE = 0,                                       /**< No Output. */
        OUTER_LOOP_ITERATES = 1,                        /**< Print stats for each outer loop iterate. */
        INNER_LOOP_ITERATES = 2,                        /**< Print stats for each inner loop iterate. */
    };


    /**
     *  Various qp solvers.
     */
    enum QPSolver {
        QPOASES_DENSE = 0,                              /**< QP solver qpOASES in dense mode. */
        QPOASES_SPARSE = 1,                             /**< QP solver qpOASES in sparse mode. */
        OSQP_SPARSE = 2                                 /**< QP solver OSQP. */
    };


    /**
     *  The utilities class
     */
    class Utilities {

        public:

            /** C = A*B **/
            static void MatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p);


            /** C = A*B **/
            static void MatrixMultiplication(const csc* const A, const double* const b, double* c);


            /** C = A'*B **/
            static void TransponsedMatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p);


            /** c = A'*b **/
            static void TransponsedMatrixMultiplication(const csc* const A, const double* const b, double* c);


            /** C += A'*B **/
            static void AddTransponsedMatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p);


            /** c += A'*b **/
            static void AddTransponsedMatrixMultiplication(const csc* const A, const double* const b, double* c);


            /** C = A'*B + B'*A **/
            static void MatrixSymmetrizationProduct(const double* const A, const double* const B, double* C, int m, int n);


            /** C = A'*B + B'*A **/
            static csc* MatrixSymmetrizationProduct(double* L_x, int* L_i, int* L_p, double* R_x, int* R_i, int* R_p, int n);


            /** C = A'*B + B'*A **/
            static csc* MatrixSymmetrizationProduct(csc* L, csc* R);


            /** d = A*b + c **/
            static void AffineLinearTransformation(const double alpha, const double* const A, const double* const b, const double* const c, double* d, int m, int n);


            /** d = A*b + c **/
            static void AffineLinearTransformation(const double alpha, const csc* const S, const double* const b, const double* const c, double* d, int m);


            /** C = alpha*A + beta*B **/
            static void WeightedMatrixAdd(const double alpha, const double* const A, const double beta, const double* const B, double* C, int m, int n);


            /** c = alpha*a + beta*b **/
            static void WeightedVectorAdd(const double alpha, const double* const a, const double beta, const double* const b, double* c, int m);


            /** @returns p' * Q * p **/
            static double QuadraticFormProduct(const double* const Q, const double* const p, int m);


            /** @return p' * Q * p **/
            static double QuadraticFormProduct(const csc* const S, const double* const p, int m);


            /** @returns a'*b **/
            static double DotProduct(const double* const a, const double* const b, int m);


            /** @returns 1-norm **/
            static double MaxAbs(const double* const a, int m);


            /** Clear sparse matrix **/
            static void ClearSparseMat(csc* M);
      static void ClearSparseMatCPP(csc* M);


            /** Clear sparse matrix **/
            static void ClearSparseMat(csc** M);


            /** Read integral data from file **/
            static ReturnValue readFromFile(int* data, int n, const char* datafilename);


            /** Read float data from file **/
            static ReturnValue readFromFile(double* data, int n, const char* datafilename );


            /** Read float data from file **/
            static ReturnValue writeToFile(double* data, int n, const char* datafilename );


            /** Print a double valued matrix **/
            static void printMatrix(const double* const A, int m, int n, const char* const name);


            /** Print an integer valued matrix **/
            static void printMatrix(const int* const A, int m, int n, const char* const name);


            /** Print dense representation of csc matrix **/
            static void printMatrix(const csc* A, const char* const name);


            /** Printing bounds **/
            static void printStep(double* xk, double* pk, double* xk_new, double alpha, int nV);


            /** Printing bounds **/
            static void printBounds(double* lb, double* xk, double* ub, int m);


            /** Construct a csc matrix (like csc_matrix in OSQP) **/
            static csc* createCSC(int m, int n, int nnz, double* x, int* i, int* p);


            /** Copy a csc matrix (like create CSC but deep copy is made) **/
            static csc* copyCSC(int m, int n, int nnz, double* x, int* i, int* p);


            /** Copy a csc matrix (override) **/
            static csc* copyCSC(const csc* const M, bool toUpperTriangular = false);


            /** Copy an integer array **/
            static void copyIntToIntT(int* dest, const int* const src, int n);


            /** Transform a csc matrix to dense.
             *
             * @param sparse A sparse matrix.
             * @param full A target pointer for the full matrix (expected to have size m*n).
             * @param m Number of rows of `Q_sparse` (in its dense representation).
             * @param n Number of columns of `Q_sparse` (in its dense representation).
             *
             * @returns An dense array representing sparse (or null pointer if failed).
             */
            static double* csc_to_dns(const csc* const sparse);


            /** Transform a dense matrix to csc.
             *
             * @param full A dense double array.
             * @param m Number of rows of `full`.
             * @param n Number of columns of `full`.
             *
             * @returns A csc pointer to the sparse matrix.
             */
            static csc* dns_to_csc(const double* const full, int m, int n);


            // Methods below this line where taken from qpOASES implementation
            /** Returns the absolute value of x. */
            static double getAbs(double x);


            /** Checks if the absolute difference between x and y is less than TOL. */
            static bool isEqual(double x, double y, double TOL);


            /** Checks if the absolute difference between x and y is less than the constant Utilities::ZERO. */
            static bool isEqual(double x, double y);


            /** Checks if the absolute value of x is less than TOL. */
            static bool isZero(double x, double TOL);


            /** Checks if the absolute value of x is less than the constant Utilities::ZERO. */
            static bool isZero(double x);


            /** Returns the sign of a double x. */
            static double getSign(double x);


            /** Returns the maximum of two integer values x and y. */
            static int getMax(int x, int y);


            /** Returns the minimum of two integer values x and y. */
            static int getMin(int x, int y);


            /** Returns the maximum of two real values x and y. */
            static double getMax(double x, double y);


            /** Returns the minimum of two real values x and y. */
            static double getMin(double x, double y);


            /** Get the index of an integer value within a sorted list. */
            static int getIndexOfIn(int val, int* sorted_lst, int beg, int end);


            /** Checks if the ptr is null. */
			template <typename PtrType>
		    static bool isNullPtr(PtrType ptr) { 
				return (ptr == 0 || ptr == nullptr);
			}


			/** Checks if the ptr is not null. */
			template <typename PtrType>
		    static bool isNotNullPtr(PtrType ptr) { 
                return (ptr != NULL && ptr != nullptr);
			}


            /** Numerical value of machine precision (min eps, s.t. 1+eps > 1).
             *	Note: this value has to be positive! */
            #ifdef __USE_SINGLE_PRECISION__
            const double EPS = 1.193e-07f;
            #else
            constexpr static double EPS = 2.221e-16;
            #endif /* __USE_SINGLE_PRECISION__ */


            /** Numerical value of zero (for situations in which it would be
             *	unreasonable to compare with 0.0).
            *	Note: this value has to be positive! */
            constexpr static double ZERO = 1.0e-25;


            /** Numerical value of infinity (e.g. for non-existing bounds).
                Note: this value has to be positive! */
            constexpr static double INFTY = 1.0e20;


            /** Maximum number of characters within a string.
             *	Note: this value should be at least 41! */
            constexpr static uint MAX_STRING_LENGTH = 160;
    };
}

#endif  // LCQPOW_UTILITIES_HPP
