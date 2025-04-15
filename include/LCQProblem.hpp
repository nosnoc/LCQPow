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


#ifndef LCQPOW_LCQPROBLEM_HPP
#define LCQPOW_LCQPROBLEM_HPP

#include "Utilities.hpp"
#include "Subsolver.hpp"
#include "OutputStatistics.hpp"
#include "Options.hpp"

#include <qpOASES.hpp>
#include <vector>
#include <deque>

using qpOASES::QProblem;

namespace LCQPow {
	class LCQProblem
	{

		/**
		 *	PUBLIC METHODS
		 */
		public:

			/** Default constructor. */
			LCQProblem( );


			/** Construct LCQP solver with dimensions.
			 *
			 * @param _nV Number of optimization variables.
			 * @param _nC Number of linear constraints.
			 * @param _nComp Number of complementarity pairs.
			 */
			LCQProblem(
				int _nV,
				int _nC,
				int _nComp
			);


			/** Destructor. */
			~LCQProblem( );


			/** Run solver passing the desired LCQP in dense format (qpOASES is used on subsolver level).
			 *
			 * @param _Q The objective's hessian matrix.
			 * @param _g The obective's linear term.
			 * @param _L The matrix selecting the left hand side of the complementarity pairs.
			 * @param _R The matrix selecting the right hand side of the complementarity pairs.
			 * @param _lbL The lower bounds associated to the complementarity matrix `_L`. A `NULL` leads to zero bounds.
			 * @param _ubL The upper bounds associated to the complementarity matrix `_L`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _lbR The lower bounds associated to the complementarity matrix `_R`. A `NULL` leads to zero boudns.
			 * @param _ubR The upper bounds associated to the complementarity matrix `_R`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _A The constraint matrix. A `NULL` pointer can be passed if no linear constraints exist.
			 * @param _lbA The lower bounds associated to the constraint matrix `_A`. A `NULL` pointer can be passed if no lower bounds exist.
			 * @param _ubA The upper bounds associated to the constraint matrix `_A`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _lb The box constraint's lower bounds. A `NULL` pointer can be passed if no lower bounds exist.
			 * @param _ub The box constraint's upper bounds. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _x0 The initial guess for the optimal primal solution vector. If a `NULL` pointer is passed the zero vector is used.
			 * @param _y0 The initial guess for the optimal dual solution vector. If a `NULL` pointer is passed, then the initialization depends on the subsolver and its options.
			 *
			 * @returns SUCCESSFUL_RETURN if a solution is found. Otherwise the return value will indicate an occured error.
			*/
			ReturnValue loadLCQP(
				const double* const _Q,
				const double* const _g,
				const double* const _L,
				const double* const _R,
				const double* const lbL = 0,
				const double* const ubL = 0,
				const double* const lbR = 0,
				const double* const ubR = 0,
				const double* const _A = 0,
				const double* const _lbA = 0,
				const double* const _ubA = 0,
				const double* const _lb = 0,
				const double* const _ub = 0,
				const double* const _x0 = 0,
				const double* const _y0 = 0
			 );


			/** Run solver passing the desired LCQP in (file) dense format (qpOASES is used on subsolver level).
			 *  All matrices are assumed to be stored row-wise.
			 *
			 * @param Q_file The objective's hessian matrix.
			 * @param g_file The obective's linear term.
			 * @param L_file The matrix selecting the left hand side of the complementarity pairs.
			 * @param R_file The matrix selecting the right hand side of the complementarity pairs.
			 * @param lbL_file The lower bounds associated to the complementarity matrix `_L`. A `NULL` leads to zero bounds.
			 * @param ubL_file The upper bounds associated to the complementarity matrix `_L`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param lbR_file The lower bounds associated to the complementarity matrix `_R`. A `NULL` leads to zero boudns.
			 * @param ubR_file The upper bounds associated to the complementarity matrix `_R`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param A_file The constraint matrix. A `NULL` pointer can be passed if no linear constraints exist.
			 * @param lbA_file The lower bounds associated to the constraint matrix `_A`. A `NULL` pointer can be passed if no lower bounds exist.
			 * @param ubA_file The upper bounds associated to the constraint matrix `_A`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param lb_file The box constraint's lower bounds. A `NULL` pointer can be passed if no lower bounds exist.
			 * @param ub_file The box constraint's upper bounds. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param x0_file The initial guess for the optimal primal solution vector. If a `NULL` pointer is passed the zero vector is used.
			 * @param y0_file The initial guess for the optimal dual solution vector. If a `NULL` pointer is passed, then the initialization depends on the subsolver and its options.
			 *
			 * @returns SUCCESSFUL_RETURN if a solution is found. Otherwise the return value will indicate an occured error.
			*/
			ReturnValue loadLCQP(
				const char* const Q_file,
				const char* const g_file,
				const char* const L_file,
				const char* const R_file,
				const char* const lbL = 0,
				const char* const ubL = 0,
				const char* const lbR = 0,
				const char* const ubR = 0,
				const char* const A_file = 0,
				const char* const lbA_file = 0,
				const char* const ubA_file = 0,
				const char* const lb_file = 0,
				const char* const ub_file = 0,
				const char* const x0_file = 0,
				const char* const y0_file = 0
			);


			/** Run solver passing the desired LCQP in sparse format (OSQP is used on subsolver level).
			 *
			 * @param _Q Hessian matrix in csc sparse format.
			 * @param _g The objective's linear term.
			 * @param _L LHS of complementarity product in csc sparse format.
			 * @param _R RHS of complementarity product in csc sparse format.
			 * @param _lbL The lower bounds associated to the complementarity matrix `_L`. A `NULL` leads to zero bounds.
			 * @param _ubL The upper bounds associated to the complementarity matrix `_L`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _lbR The lower bounds associated to the complementarity matrix `_R`. A `NULL` leads to zero boudns.
			 * @param _ubR The upper bounds associated to the complementarity matrix `_R`. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _A Constraint matrix in csc sparse format.
			 * @param _lbA The constraints lower bounds. A `NULL` pointer can be passed if no lower bounds exist.
			 * @param _ubA The constraints upper bounds. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _lb The box constraints lower bounds. A `NULL` pointer can be passed if no lower bounds exist.
			 * @param _ub The box constraints upper bounds. A `NULL` pointer can be passed if no upper bounds exist.
			 * @param _x0 The initial guess for the optimal primal solution vector. If a `NULL` pointer is passed the zero vector is used.
			 * @param _y0 The initial guess for the optimal dual solution vector. If a `NULL` pointer is passed, then the initialization depends on the subsolver and its options.
			 *
			 * @returns SUCCESSFUL_RETURN if a solution is found. Otherwise the return value will indicate an occured error.
			*/
			ReturnValue loadLCQP(
				const csc* const _Q,
				const double* const _g,
				const csc* const _L,
				const csc* const _R,
				const double* const lbL = 0,
				const double* const ubL = 0,
				const double* const lbR = 0,
				const double* const ubR = 0,
				const csc* const _A = 0,
				const double* const _lbA = 0,
				const double* const _ubA = 0,
				const double* const _lb = 0,
				const double* const _ub = 0,
				const double* const _x0 = 0,
				const double* const _y0 = 0
			);


			/** Switch to sparse mode (if initialized with dense data but want to use sparse solver). */
			ReturnValue switchToSparseMode( );


			/** Switch to dense mode (if initialized with sparse data but want to use dense solver). */
			ReturnValue switchToDenseMode( );


			/** After problem is set up, call this function and solve the LCQP.
			 *
			 * @returns SUCCESSFUL_RETURN if a solution is found. Otherwise the return value will indicate an occured error. */
			ReturnValue runSolver( );


			/** Writes the primal solution vector.
			 *
			 * @param xOpt A pointer to the desired primal solution storage vector.
			 *
			 * @returns If the problem was solved successfully the stationarity type is passed. Else PROBLEM_NOT_SOLVED is returned.
			 */
			virtual AlgorithmStatus getPrimalSolution( double* const xOpt ) const;


			/** Writes the dual solution vector.
			 *
			 * @param yOpt A pointer to the desired dual solution storage vector (assumed to be allocated).
			 *
			 * @returns If the problem was solved successfully the stationarity type is passed. Else PROBLEM_NOT_SOLVED is returned.
			 */
			virtual AlgorithmStatus getDualSolution( double* const yOpt ) const;


			/** Get the number of primal variables. 
			 *
			 * @return Returns the number of primal variables.
			 */
			int getNumberOfPrimals( ) const;


			/** Get the number of dual variables. This depends on the utilized subproblem solver.
			 *
			 * @return Returns the number of dual variables.
			 */
			virtual int getNumberOfDuals( ) const;


			/** Get the output statistics.
			 *
			 * @param _stats The output statistics pointer to copy the stats to.
			 */
			virtual void getOutputStatistics( OutputStatistics& stats) const;


			/** Pass options for the LCQP.
			 *
			 * @param _options Options to be used.
			 */
			inline void setOptions(	const Options& _options	);

    /** Get nV
     *
     * @return Returns nV.
     */
    int getNV() const;

    /** Get nV
     *
     * @return Returns nC.
     */
    int getNC() const;

    /** Get nV
     *
     * @return Returns nComp.
     */
    int getNComp() const;

    /** Get (a copy of) Options object
     *
     * @return Returns deep copy of current Options object.
     */
     Options getOptions() const;


		/**
    	 *	PROTECTED METHODS
		 */
		protected:
			/** Clears all memory. */
			void clear( );


			/** Prints concise information on the current iteration. */
			void printIteration( );


			/** Print header every once in a while. */
			void printHeader();


			/** Print line (mainly for printing header). */
			void printLine();


			/** Set the complementarity bounds. */
			ReturnValue setComplementarityBounds(const double* const lbL_new, const double* const ubL_new, const double* const lbR_new, const double* const ubR_new);


			/** Store the (dense) Hessian matrix Q internally.
			 *
			 * @param Q_new New dense Hessian matrix (with correct dimension!), a shallow copy is made.
			 */
			inline ReturnValue setQ( const double* const Q_new );


			/** Store the (sparse) Hessian matrix Q internally.
			 *
			 * @param Q_new Hessian matrix in csc sparse format.
			 */
			inline ReturnValue setQ( const csc* const Q_new );


			/** Store the objective linear component internally.
			 *
			 * @param g_new Objective linear component.
			 */
			inline ReturnValue setG( const double* const g_new );


			/** Store the lower (box) bounds internally.
			 *
			 * @param lb_new New lower bound vector (with correct dimension!).
			 */
			inline ReturnValue setLB( const double* const lb_new );


			/** Store a specific lower (box) bound internally.
			 *
			 * @param number Number of entry to be changed.
			 * @param value New value for entry of lower bound vector.
			 */
			inline ReturnValue setLB(
				int number,
				double value
			);


			/** Store the upper (box) bounds internally.
			 *
			 * @param ub_new New upper bound vector (with correct dimension!).
			 */
			inline ReturnValue setUB( const double* const ub_new );


			/** Store a specific upper (box) bound internally.
			 *
			 * @param number Number of entry to be changed.
			 * @param value New value for entry of lower bound vector.
			 */
			inline ReturnValue setUB(
				int number,
				double value
			);


			/** Set the new linear constraint consisting of (dense) complementarity pairs and regular linear constraints.
			 *
			 * @param L_new New lhs complementarity matrix.
			 * @param R_new New rhs complementarity matrix.
			 * @param A_new New constraint matrix.
			 * @param lbA New lower bounds for A.
			 * @param ubA New upper bounds for A.
			 */
			ReturnValue setConstraints(
				const double* const L_new,
				const double* const R_new,
				const double* const A_new,
				const double* const lbA,
				const double* const ubA
			);


			/** Set the new linear constraint consisting of (dense) complementarity pairs and regular linear constraints.
			 *
			 * @param L_new LHS of complementarity product in csc sparse format.
			 * @param R_new RHS of complementarity product in csc sparse format.
			 * @param A_new Constraint matrix in csc sparse format.
			 * @param lbA The constraint's lower bounds. A `NULL` pointer can be passed if no lower bounds exist.
			 * @param ubA The constraint's upper bounds. A `NULL` pointer can be passed if no upper bounds exist.
			 */
			ReturnValue setConstraints(
				const csc* const L_new,
				const csc* const R_new,
				const csc* const A_new,
				const double* const lbA,
				const double* const ubA
			);


			/** Set the complementarity matrix (requires the constraints to be set) as the symmetrization product of L and R.
			 *  C = L'*R + R'*L
			 */
			ReturnValue setC( );


			/** Sets the initial guess x0 and y0.
			 *
			 * @param _x0 The primal initial guess.
			 * @param _y0 The dual initial guess.
			 */
			inline ReturnValue setInitialGuess(
				const double* const _x0,
				const double* const _y0
			);


			/** Set Qk. */
			void setQk( );


		/**
		 *	PROTECTED MEMBER VARIABLES
		 */
		protected:
			Options options;						/**< Class for algorithmic options. */


		/**
	     *	PRIVATE STUFF
		 */
		private:

			/** Called in runSolver to initialize variables. */
			ReturnValue initializeSolver( );

			/** Update the penalty linearization. */
			void updateLinearization( );

			/** Solves the qp subproblem wrt Q, gk, A, L, R .
			 *
			 * @param initialSolve Pass true on first solve of sequence, false on subsequent calls (initialization vs hotstart).
			 */
			ReturnValue solveQPSubproblem( bool initialSolve );

			/** Check outer stationarity at current iterate xk. */
			bool stationarityCheck( );

			/** Check satisfaction of complementarity value. */
			bool complementarityCheck( );

			/** Evaluate objective function at current iterate. */
			double getObj( );

			/** Evaluate penalty function at current iterate. */
			double getPhi( );

			/** Evaluate merit function at current iterate. */
			double getMerit( );

			/** Perform penalty update. */
			void updatePenalty( );

			/** Get optimal step length. */
			void getOptimalStepLength( );

			/** Update xk and gk. */
			void updateStep( );

			/** Update gradient of Lagrangian. */
			void updateStationarity( );

			/** Check the dynamic penalty update strategy by Leyffer. */
			bool leyfferCheckPositive( );

			/** Update Qk. */
			void updateQk( );

			/** Update outer iteration counter. */
			void updateOuterIter( );

			/** Update outer iteration counter. */
			void updateTotalIter( );

			/** Gradient perturbation method. */
			void perturbGradient( );

			/** Step perturbation method. */
			void perturbStep( );

			/** Store detailed steps to output stats. */
			void storeSteps( );

			/** Transform the dual variables from penalty form to LCQP form. */
			void transformDuals( );

			/** Determine stationarity type of optimal solution. */
			void determineStationarityType( );

			/** Get indices of weak complementarites. */
			std::vector<int> getWeakComplementarities( );

			int nV;									/**< Number of variables. */
			int nC;									/**< Number of constraints. */
			int nComp;								/**< Number of complementarity constraints. */
			int nDuals; 							/**< Number of duals variables. */
			int boxDualOffset;						/**< Offset for linear constraint duals (i.e. 0 if no BC (Box Constraints) exist nV if BC exist). */

			double* Q = NULL;						/**< Objective Hessian term. */

			double* g = NULL;						/**< Objective linear term. */

			double* lb = NULL;						/**< Lower bound vector (on variables). */
			double* ub = NULL;						/**< Upper bound vector (on variables). */
			double* lb_tmp = NULL;					/**< Temporary box constraints. */
			double* ub_tmp = NULL;					/**< Temporary box constraints. */

			double* A = NULL;						/**< Constraint matrix. */
			double* lbA = NULL;						/**< Lower bound vector (on constraints). */
			double* ubA = NULL;						/**< Upper bound vector (on constraints). */

			double* L = NULL;						/**< LHS of complementarity product. */
			double* R = NULL;						/**< RHS of complementarity product. */
			double* C = NULL;						/**< Complementarity matrix (L'*R + R'*L). */
			double* lbL = NULL;						/**< LHS Complementarity lower bounds. */
			double* ubL = NULL;						/**< LHS Complementarity upper bounds. */
			double* lbR = NULL;						/**< RHS Complementarity lower bounds. */
			double* ubR = NULL;						/**< RHS Complementarity upper bounds. */
			double* g_phi = NULL;					/**< Linear Term of phi -(l_L'*R + l_R'*L). */
			double phi_const = 0;					/**< Constant phi expression (l_L'*l_R). */

			double rho; 							/**< Current penalty value. */

			double* g_tilde = NULL;					/**< Current linear terms (g + rhok*g_phi). Updated once per inner loop. */
			double* gk = NULL;						/**< Current objective linear term. */

			double* xk = NULL;						/**< Current primal iterate. */
			double* yk = NULL;						/**< Current dual vector. */
			double* yk_A = NULL;					/**< Current dual vector w.r.t A. */
			double* xnew = NULL;					/**< Current qpSubproblem solution. */
			double* pk = NULL;						/**< xnew - xk. */

			double alphak; 							/**< Optimal step length. */
			double* lk_tmp = NULL;					/**< An auxiliar vector to help compute lkj. */

			double* Qk = NULL;						/**< Q + rho*C, required for stationarity and optimal step length. */
			double* statk = NULL;					/**< Stationarity of current iterate. */
			double* constr_statk = NULL;			/**< Constraint contribution to stationarity equation. */
			double* box_statk = NULL;				/**< Box Constraint contribution to stationarity equation. */

			int outerIter;							/**< Outer iterate counter. */
			int innerIter;							/**< Inner iterate counter. */
			int totalIter;							/**< Total iterate counter. */

			int qpIterk;							/**< Iterations taken by qpSolver to solve subproblem. */
			int qpSolverExitFlag;					/**< Most recent exit flag of QP solver. */
			AlgorithmStatus algoStat;				/**< Status of algorithm. */

			bool sparseSolver = false;				/**< Whether to use sparse algebra or dense. */

			csc* Q_sparse = NULL;					/**< Sparse objective Hessian matrix. */
			csc* A_sparse = NULL;					/**< Sparse constraint matrix. */
			csc* L_sparse = NULL;					/**< Sparse L. */
			csc* R_sparse = NULL;					/**< Sparse R. */
			csc* C_sparse = NULL;					/**< Sparse C. */
			csc* Qk_sparse = NULL;					/**< Sparse Qk. */
			std::vector<int> Qk_indices_of_C;		/**< Remember the indices of Qk corresponding to C (for fast Qk update). */

			std::deque<double> complHistory; 		/**< Vector containing the previous complementarity values. */

			Subsolver subsolver;					/**< Subsolver class for solving the QP subproblems. */

			OutputStatistics stats;					/**< Output statistics. */
	};
}


#include "LCQProblem.ipp"

#endif	/* LCQPOW_LCQPROBLEM_HPP */
