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


#include "LCQProblem.hpp"
#include "Utilities.hpp"
#include "MessageHandler.hpp"
#include "SubsolverQPOASES.hpp"
#include "SubsolverOSQP.hpp"

#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>

#include <qpOASES.hpp>
  
using qpOASES::QProblem;

namespace LCQPow {

	LCQProblem::LCQProblem() {}


	LCQProblem::LCQProblem(int _nV, int _nC, int _nComp)
	{
		/* consistency checks */
		if ( _nV <= 0 )
		{
			MessageHandler::PrintMessage( INVALID_NUMBER_OF_OPTIM_VARS, ERROR );
			return;
		}

		if ( _nComp <= 0 )
		{
			MessageHandler::PrintMessage( INVALID_NUMBER_OF_OPTIM_VARS, ERROR );
			return;
		}

		if ( _nC < 0 )
		{
			_nC = 0;
			MessageHandler::PrintMessage( INVALID_NUMBER_OF_CONSTRAINT_VARS, ERROR );
			return;
		}

		nV = _nV;
		nC = _nC;
		nComp = _nComp;

		// Allocate auxiliar vectors
		Qk = new double[nV*nV]();
		gk = new double[nV]();
		xnew = new double[nV]();
		yk_A = new double[nC + 2*nComp]();
		pk = new double[nV]();
		statk = new double[nV]();
		constr_statk = new double[nV]();
		box_statk = new double[nV]();
		lk_tmp = new double[nV]();
	}


	LCQProblem::~LCQProblem()
  {
		clear();
	}


	ReturnValue LCQProblem::loadLCQP(const double* const _Q, const double* const _g,
                                   const double* const _L, const double* const _R,
                                   const double* const _lbL, const double* const _ubL,
                                   const double* const _lbR, const double* const _ubR,
                                   const double* const _A, const double* const _lbA, const double* const _ubA,
                                   const double* const _lb, const double* const _ub,
                                   const double* const _x0, const double* const _y0)
	{
		ReturnValue ret;

		if ( nV <= 0 || nComp <= 0 )
            return( MessageHandler::PrintMessage(ReturnValue::LCQPOBJECT_NOT_SETUP, ERROR) );

		ret = setQ( _Q );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		// Only copy lb and ub to temporary variables
		// Build them once we know what solver is used
		if (Utilities::isNotNullPtr(_lb)) {
			lb_tmp = new double[nV];
			memcpy(lb_tmp, _lb, (size_t) nV*sizeof(double));
		}

		if (Utilities::isNotNullPtr(_ub)) {
			ub_tmp = new double[nV];
			memcpy(ub_tmp, _ub, (size_t) nV*sizeof(double));
		}

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setConstraints( _L, _R, _A, _lbA, _ubA );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setComplementarityBounds( _lbL, _ubL, _lbR, _ubR );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		sparseSolver = false;

		return ReturnValue::SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::loadLCQP(const char* const Q_file, const char* const g_file,
                                   const char* const L_file, const char* const R_file,
                                   const char* const lbL_file, const char* const ubL_file,
                                   const char* const lbR_file, const char* const ubR_file,
                                   const char* const A_file, const char* const lbA_file, const char* const ubA_file,
                                   const char* const lb_file, const char* const ub_file,
                                   const char* const x0_file, const char* const y0_file)
	{
		ReturnValue ret;

		double* _Q = new double[nV*nV];
		ret = Utilities::readFromFile( _Q, nV*nV, Q_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _Q;
			return MessageHandler::PrintMessage( ret, ERROR );
		}

		double* _g = new double[nV];
		ret = Utilities::readFromFile( _g, nV, g_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _g;
			return MessageHandler::PrintMessage( ret, ERROR );
		}

		double* _L = new double[nComp*nV];
		ret = Utilities::readFromFile( _L, nComp*nV, L_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _L;
			return MessageHandler::PrintMessage( ret, ERROR );
		}

		double* _R = new double[nComp*nV];
		ret = Utilities::readFromFile( _R, nComp*nV, R_file );
		if ( ret != SUCCESSFUL_RETURN ) {
			delete[] _R;
			return MessageHandler::PrintMessage( ret, ERROR );
		}

		double* _lbL = NULL;
		if (Utilities::isNotNullPtr(lbL_file)) {
			_lbL = new double[nComp];
			ret = Utilities::readFromFile( _lbL, nComp, lbL_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _lbL;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _ubL = NULL;
		if (Utilities::isNotNullPtr(ubL_file)) {
			_ubL = new double[nComp];

			ret = Utilities::readFromFile( _ubL, nComp, ubL_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _ubL;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _lbR = NULL;
		if (Utilities::isNotNullPtr(lbR_file)) {
			_lbR = new double[nComp];

			ret = Utilities::readFromFile( _lbR, nComp, lbR_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _lbR;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _ubR = NULL;
		if (Utilities::isNotNullPtr(ubR_file)) {
			_ubR = new double[nComp];

			ret = Utilities::readFromFile( _ubR, nComp, ubR_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _ubR;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _A = NULL;
		if (Utilities::isNotNullPtr(A_file)) {
			_A = new double[nC*nV];
			ret = Utilities::readFromFile( _A, nC*nV, A_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _A;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _lbA = NULL;
		if (Utilities::isNotNullPtr(lbA_file)) {
			_lbA = new double[nC];
			ret = Utilities::readFromFile( _lbA, nC, lbA_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _lbA;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _ubA = NULL;
		if (Utilities::isNotNullPtr(ubA_file)) {
			_ubA = new double[nC];
			ret = Utilities::readFromFile( _ubA, nC, ubA_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _ubA;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _lb = NULL;
		if (Utilities::isNotNullPtr(lb_file)) {
			_lb = new double[nV];
			ret = Utilities::readFromFile( _lb, nV, lb_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _lb;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _ub = NULL;
		if (Utilities::isNotNullPtr(ub_file)) {
			_ub = new double[nV];
			ret = Utilities::readFromFile( _ub, nV, ub_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _ub;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _x0 = NULL;
		if (Utilities::isNotNullPtr(x0_file)) {
			_x0 = new double[nV];
			ret = Utilities::readFromFile( _x0, nV, x0_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _x0;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		double* _y0 = NULL;
		if (Utilities::isNotNullPtr(y0_file)) {
			_y0 = new double[nC + 2*nComp];
			ret = Utilities::readFromFile( _y0, nC + 2*nComp, y0_file );
			if ( ret != SUCCESSFUL_RETURN ) {
				delete[] _y0;
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		// Fill vaues
		ret = setQ( _Q );
		delete[] _Q;

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setG( _g );
		delete[] _g;

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		// Only copy lb and ub to temporary variables
		// Build them once we know what solver is used
		if (Utilities::isNotNullPtr(_lb)) {
			lb_tmp = new double[nV];
			memcpy(lb_tmp, _lb, (size_t) nV*sizeof(double));
			delete[] _lb;
		}

		if (Utilities::isNotNullPtr(_ub)) {
			ub_tmp = new double[nV];
			memcpy(ub_tmp, _ub, (size_t) nV*sizeof(double));
			delete[] _ub;
		}

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR);

		ret = setConstraints( _L, _R, _A, _lbA, _ubA );

		delete[] _L; delete[] _R;

		if (Utilities::isNotNullPtr(_A)) {
			delete[] _A;
			_A = NULL;
		}

		if (Utilities::isNotNullPtr(_lbA)) {
			delete[] _lbA;
			_lbA = NULL;
		}

		if (Utilities::isNotNullPtr(_ubA)) {
			delete[] _ubA;
			_ubA = NULL;
		}

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setComplementarityBounds( _lbL, _ubL, _lbR, _ubR );

		if (Utilities::isNotNullPtr(_lbL)) {
			delete[] _lbL;
			_lbL = NULL;
		}

		if (Utilities::isNotNullPtr(_ubL)) {
			delete[] _ubL;
			_ubL = NULL;
		}

		if (Utilities::isNotNullPtr(_lbR)) {
			delete[] _lbR;
			_lbR = NULL;
		}

		if (Utilities::isNotNullPtr(_ubR)) {
			delete[] _ubR;
			_ubR = NULL;
		}

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setInitialGuess( _x0, _y0 );

		if (Utilities::isNotNullPtr(_x0)) delete[] _x0;
		if (Utilities::isNotNullPtr(_y0)) delete[] _y0;

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		sparseSolver = false;

		return ReturnValue::SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::loadLCQP(const csc* const _Q, const double* const _g,
                                   const csc* const _L, const csc* const _R,
                                   const double* const _lbL, const double* const _ubL,
                                   const double* const _lbR, const double* const _ubR,
                                   const csc* const _A, const double* const _lbA, const double* const _ubA,
                                   const double* const _lb, const double* const _ub,
                                   const double* const _x0, const double* const _y0)
	{
		ReturnValue ret;

		ret = setQ( _Q );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setG( _g );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setConstraints( _L, _R, _A, _lbA, _ubA );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setComplementarityBounds( _lbL, _ubL, _lbR, _ubR );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		ret = setInitialGuess( _x0, _y0 );

		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		// Only copy lb and ub to temporary variables
		// Build them once we know what solver is used
		if (Utilities::isNotNullPtr(_lb)) {
			lb_tmp = new double[nV];
			memcpy(lb_tmp, _lb, (size_t) nV*sizeof(double));
		}

		if (Utilities::isNotNullPtr(_ub)) {
			ub_tmp = new double[nV];
			memcpy(ub_tmp, _ub, (size_t) nV*sizeof(double));
		}

		sparseSolver = true;

		return ReturnValue::SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::runSolver()
	{
		// Initialize variables
		ReturnValue ret = initializeSolver();
		if (ret != SUCCESSFUL_RETURN)
			return MessageHandler::PrintMessage( ret, ERROR );

		// Initialization strategy
		if (options.getSolveZeroPenaltyFirst()) {

			// Zero penalty, i.e. pen-linearization = 0, i.e. gk = g
			memcpy(gk, g, (size_t)nV*sizeof(double));
			ret = solveQPSubproblem( true );
			if (ret != SUCCESSFUL_RETURN) {
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		} else {
			// Linearize penalty function at initial guess
			updateLinearization();
			ret = solveQPSubproblem( true );
			if (ret != SUCCESSFUL_RETURN) {
				return MessageHandler::PrintMessage( ret, ERROR );
			}
		}

		// Initialize Qk = Q + rhok*C
		setQk();

		// Initialize stats.rho_opt
		stats.updateRhoOpt( rho );

		// Outer and inner loop in one
		while ( true ) {

			// Update xk, Qk, stationarity
			updateStep( );

			// Update gradient of Lagrangian
			updateStationarity( );

			// Print iteration
			printIteration( );

			// Store steps if desired
			if (options.getStoreSteps()) {
				storeSteps( );
			}

			// Update the total iteration counter
			updateTotalIter();

			// Update inner iterate counter
			innerIter++;

			// Perform Dynamic Leyffer Strategy
			if (leyfferCheckPositive( )) {
				updatePenalty( );

				// Update iterate counters
				updateOuterIter();
				innerIter = 0;
			}

			// gk = new linearization + g
			updateLinearization();

			// Terminate, update pen, or continue inner loop
			if (stationarityCheck()) {
				if (complementarityCheck()) {
					// Switch from penalized to LCQP duals
					transformDuals();

					// Determine C-,M-,S-Stationarity
					determineStationarityType();

					// Update output statistics
					stats.updateSolutionStatus( algoStat );

					// Print solution type
					if (options.getPrintLevel() > PrintLevel::NONE)
						MessageHandler::PrintSolution( algoStat );

					return SUCCESSFUL_RETURN;
				} else {
					updatePenalty();

					// Update iterate counters
					updateOuterIter();
					innerIter = 0;
				}
			}

			// (Failed) termination condition due to number of iterations
			if ( totalIter > options.getMaxIterations() )
				return MAX_ITERATIONS_REACHED;

			// (Failed) termination condition due to penalty value
			if ( rho > options.getMaxPenaltyParameter() )
				return MAX_PENALTY_REACHED;

			// gk = new linearization + g
			updateLinearization();

			// Step computation
			ret = solveQPSubproblem( false );
			if (ret != SUCCESSFUL_RETURN) {
				return MessageHandler::PrintMessage( ret, ERROR );
			}

			// Add some +/- EPS to each coordinate
			if (options.getPerturbStep())
				perturbStep();

			// Step length computation
			getOptimalStepLength( );
		}
	}


	ReturnValue LCQProblem::setConstraints(const double* const L_new, const double* const R_new,
                                         const double* const A_new, const double* const lbA_new, const double* const ubA_new)
	{
		if ( nV == 0 || nComp == 0 )
			return LCQPOBJECT_NOT_SETUP;

		if ( Utilities::isNullPtr(A_new) && nC > 0)
			return INVALID_CONSTRAINT_MATRIX;

		// Set up new constraint matrix (A; L; R)
		A = new double[(nC + 2*nComp)*nV];

		for (int i = 0; i < nC*nV; i++)
			A[i] = A_new[i];

		for (int i = 0; i < nComp*nV; i++)
			A[i + nC*nV] = L_new[i];

		for (int i = 0; i < nComp*nV; i++)
			A[i + nC*nV + nComp*nV] = R_new[i];

		// Set up new constraint bounds (lbA; 0; 0) & (ubA; INFINITY; INFINITY)
		lbA = new double[nC + 2*nComp];
		ubA = new double[nC + 2*nComp];

		if (Utilities::isNotNullPtr(lbA_new)) 
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = lbA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = -INFINITY;
		}

		if (Utilities::isNotNullPtr(ubA_new)) 
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = ubA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = INFINITY;
		}

		// Set complementarities
		if ( Utilities::isNullPtr(L_new) || Utilities::isNullPtr(R_new) )
			return INVALID_COMPLEMENTARITY_MATRIX;

		L = new double[nComp*nV];
		R = new double[nComp*nV];

		for (int i = 0; i < nComp*nV; i++) {
			L[i] = L_new[i];
			R[i] = R_new[i];
		}

		C = new double[nV*nV];
		Utilities::MatrixSymmetrizationProduct(L, R, C, nComp, nV);

		return SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::setConstraints(const csc* const L_new, const csc* const R_new,
                                         const csc* const A_new, const double* const lbA_new, const double* const ubA_new)
	{
		// Create sparse matrices
		L_sparse = Utilities::copyCSC(L_new);
		R_sparse = Utilities::copyCSC(R_new);

		// Get number of elements
    for(int ii=0; ii<nV+1; ii++)
      std::cout << ii <<" L_sparse->p[ii]: " << L_sparse->p[ii] << std::endl;
    std::cout << "L_sparse->p[nV]: " << L_sparse->p[nV] << std::endl;
		int tmpA_nnx = L_sparse->p[nV] + R_sparse->p[nV];

		if (Utilities::isNotNullPtr(A_new)) {
			tmpA_nnx += Utilities::isNotNullPtr(A_new->p) ? A_new->p[nV] : 0;
		}

		// Data array
		double* tmpA_data = (double*)malloc((size_t)tmpA_nnx*sizeof(double));

		// Row indices
		int* tmpA_i = (int*)malloc((size_t)tmpA_nnx*sizeof(int));

		// Column pointers
		int* tmpA_p = (int*)malloc((size_t)(nV+1)*sizeof(int));

		int index_data = 0;
		tmpA_p[0] = 0;

		// Iterate over columns
		for (int i = 0; i < nV; i++) {
			tmpA_p[i+1] = tmpA_p[i];

			// First handle rows of A
			if (Utilities::isNotNullPtr(A_new)) {
				for (int j = A_new->p[i]; j < A_new->p[i+1]; j++) {
					tmpA_data[index_data] = A_new->x[j];
					tmpA_i[index_data] = A_new->i[j];
					index_data++;
					tmpA_p[i+1]++;
				}
			}

			// Then rows of L
			for (int j = L_sparse->p[i]; j < L_sparse->p[i+1]; j++) {
				tmpA_data[index_data] = L_sparse->x[j];
				tmpA_i[index_data] = nC + L_sparse->i[j];
				index_data++;
				tmpA_p[i+1]++;
			}

			// Then rows of R
			for (int j = R_sparse->p[i]; j < R_sparse->p[i+1]; j++) {
				tmpA_data[index_data] = R_sparse->x[j];
				tmpA_i[index_data] = nC + nComp + R_sparse->i[j];
				index_data++;
				tmpA_p[i+1]++;
			}
		}

		// Create sparse matrix
		A_sparse = Utilities::createCSC(nC + 2*nComp, nV, tmpA_nnx, tmpA_data, tmpA_i, tmpA_p);

		// Set up new constraint bounds (lbA; 0; 0) & (ubA; INFINITY; INFINITY)
		lbA = new double[nC + 2*nComp];
		ubA = new double[nC + 2*nComp];

		if (Utilities::isNotNullPtr(lbA_new)) 
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = lbA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				lbA[i] = -INFINITY;
		}

		if (Utilities::isNotNullPtr(ubA_new)) 
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = ubA_new[i];
		}
		else
		{
			for (int i = 0; i < nC; i++)
				ubA[i] = INFINITY;
		}

		C_sparse = Utilities::MatrixSymmetrizationProduct(L_sparse, R_sparse);

		if (Utilities::isNullPtr(C_sparse)) {
			return FAILED_SYM_COMPLEMENTARITY_MATRIX;
		}

		return SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::setComplementarityBounds(const double* const lbL_new, const double* const ubL_new, const double* const lbR_new, const double* const ubR_new)
  {
		if (Utilities::isNotNullPtr(lbL_new )) {
			lbL = new double[nComp];
		}

		if (Utilities::isNotNullPtr(ubL_new )) {
			ubL = new double[nComp];
		}

		if (Utilities::isNotNullPtr(lbR_new )) {
			lbR = new double[nComp];
		}

		if (Utilities::isNotNullPtr(ubR_new )) {
			ubR = new double[nComp];
		}

		// Bounds on Lx
		for (int i = 0; i < nComp; i++) {
			if (Utilities::isNotNullPtr(lbL_new )) { 
				if (lbL_new[i] <= -INFINITY)
					return INVALID_LOWER_COMPLEMENTARITY_BOUND;

				lbL[i] = lbL_new[i];
				lbA[nC + i] = lbL_new[i];
			} else {
				lbA[nC + i] = 0;
			}

			if (Utilities::isNotNullPtr(ubL_new )) { 
				ubL[i] = ubL_new[i];
				ubA[nC + i] = ubL_new[i];
			} else {
				ubA[nC + i] = INFINITY;
			}
		}

		// Bounds on Rx
		for (int i = 0; i < nComp; i++) {
			if (Utilities::isNotNullPtr(lbR_new )) { 
				if (lbR_new[i] <= -INFINITY)
					return INVALID_LOWER_COMPLEMENTARITY_BOUND;

				lbR[i] = lbR_new[i];
				lbA[nC + nComp + i] = lbR_new[i];
			} else {
				lbA[nC + nComp + i] = 0;
			}

			if (Utilities::isNotNullPtr(ubR_new )) { 
				ubR[i] = ubR_new[i];
				ubA[nC + nComp + i] = ubR_new[i];
			} else {
				ubA[nC + nComp + i] = INFINITY;
			}
		}

		return SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::setQ(const csc* const Q_new)
	{
		if (nV <= 0)
			return LCQPOBJECT_NOT_SETUP;

		Q_sparse = Utilities::copyCSC(Q_new);
    for(int ii=0; ii<nV+1; ii++)
      std::cout << ii <<" Q_sparse->p[ii]: " << Q_sparse->p[ii] << std::endl;

		return ReturnValue::SUCCESSFUL_RETURN;
	}


	void LCQProblem::setQk()
	{
		if (sparseSolver) {
			std::vector<double> Qk_data;
			std::vector<int> Qk_row;
			int* Qk_p = (int*)malloc((size_t)(nV+1)*sizeof(int));
			Qk_p[0] = 0;

			// Iterate over columns
			for (int j = 0; j < nV; j++) {
				Qk_p[j+1] = Qk_p[j];

				int idx_Q = Q_sparse->p[j];
				int idx_C = C_sparse->p[j];

				while (idx_Q < Q_sparse->p[j+1] || idx_C < C_sparse->p[j+1]) {
					// Only elements of Q left in this column
					if (idx_Q < Q_sparse->p[j+1] && idx_C >= C_sparse->p[j+1]) {
						Qk_data.push_back(Q_sparse->x[idx_Q]);
						Qk_row.push_back(Q_sparse->i[idx_Q]);

						idx_Q++;
					}

					// Only elements of C left in this column
					else if (idx_Q >= Q_sparse->p[j+1] && idx_C < C_sparse->p[j+1]) {
						Qk_data.push_back(rho*C_sparse->x[idx_C]);
						Qk_row.push_back(C_sparse->i[idx_C]);

						Qk_indices_of_C.push_back(Qk_p[j+1]);

						idx_C++;
					}

					// Element of Q is in higher row than C
					else if (Q_sparse->i[idx_Q] < C_sparse->i[idx_C]) {
						Qk_data.push_back(Q_sparse->x[idx_Q]);
						Qk_row.push_back(Q_sparse->i[idx_Q]);

						idx_Q++;
					}

					// Element of C are in higher row than Q
					else if (Q_sparse->i[idx_Q] > C_sparse->i[idx_C]) {
						Qk_data.push_back(rho*C_sparse->x[idx_C]);
						Qk_row.push_back(C_sparse->i[idx_C]);

						Qk_indices_of_C.push_back(Qk_p[j+1]);

						idx_C++;
					}

					// Add both and update index only once!
					else {
						Qk_data.push_back(Q_sparse->x[idx_Q] + rho*C_sparse->x[idx_C]);
						Qk_row.push_back(Q_sparse->i[idx_Q]);

						idx_Q++;
						idx_C++;

						Qk_indices_of_C.push_back(Qk_p[j+1]);
					}

					Qk_p[j+1]++;
				}
			}

			int Qk_nnx = (int) Qk_data.size();
			double* Qk_x = (double*)malloc((size_t)Qk_nnx*sizeof(double));
			int* Qk_i =  (int*)malloc((size_t)Qk_nnx*sizeof(int));

			for (size_t i = 0; i < (size_t) Qk_nnx; i++) {
				Qk_x[i] = Qk_data[i];
				Qk_i[i] = Qk_row[i];
			}

			Qk_sparse = Utilities::createCSC(nV, nV, Qk_nnx, Qk_x, Qk_i, Qk_p);

			Qk_data.clear();
			Qk_row.clear();
		} else {
			Utilities::WeightedMatrixAdd(1, Q, rho, C, Qk, nV, nV);
		}
	}


	ReturnValue LCQProblem::initializeSolver()
	{
		ReturnValue ret = SUCCESSFUL_RETURN;
		if (options.getQPSolver() == QPSolver::QPOASES_DENSE) {
			nDuals = nV + nC + 2*nComp;
			boxDualOffset = nV;

			if (sparseSolver) {
				return DENSE_SPARSE_MISSMATCH;
			}

			ret = setLB( lb_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;

			ret = setUB ( ub_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;

			Subsolver tmp(nV, nC + 2*nComp, Q, A);
			subsolver = tmp;
		} else if (options.getQPSolver() == QPSolver::QPOASES_SPARSE) {
			nDuals = nV + nC + 2*nComp;
			boxDualOffset = nV;

			if (!sparseSolver) {
				return DENSE_SPARSE_MISSMATCH;
			}

			ret = setLB( lb_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;

			ret = setUB( ub_tmp );

			if (ret != SUCCESSFUL_RETURN)
				return ret;

			Subsolver tmp(nV, nC + 2*nComp, Q_sparse, A_sparse, options.getQPSolver());
			subsolver = tmp;

		} else if (options.getQPSolver() == QPSolver::OSQP_SPARSE) {
			if (Utilities::isNotNullPtr(lb) || Utilities::isNotNullPtr(ub)) {
				return INVALID_OSQP_BOX_CONSTRAINTS;
			}

			nDuals = nC + 2*nComp;
			boxDualOffset = 0;

			// If solving with OSQP we ignore the dual guess on the box constraints
			if (Utilities::isNotNullPtr(yk)) {
				double* yk_tmp = new double[nDuals];

				// Shift the duals
				for (int i = 0; i < nDuals; i++)
					yk_tmp[i] = yk[nV + i];

				delete[] yk;
				yk = new double[nDuals];
				memcpy(yk, yk_tmp, (size_t)nDuals);
				delete[] yk_tmp;
			}

			if (!sparseSolver) {
				return DENSE_SPARSE_MISSMATCH;
			}

			if (Utilities::isNotNullPtr(lb_tmp) || Utilities::isNotNullPtr(ub_tmp)) {
				return ReturnValue::INVALID_OSQP_BOX_CONSTRAINTS;
			}

			Subsolver tmp(nV, nDuals, Q_sparse, A_sparse, options.getQPSolver());
			subsolver = tmp;
		} else {
			return ReturnValue::NOT_YET_IMPLEMENTED;
		}

		// Linear objective component
		g_tilde = new double[nV];
		memcpy(g_tilde, g, (size_t)nV*sizeof(double));

		// Phi expressions
		if (Utilities::isNotNullPtr(lbL) || Utilities::isNotNullPtr(lbR)) 
			phi_const = Utilities::DotProduct(lbL, lbR, nComp);

		// g_phi
		if (Utilities::isNotNullPtr(lbL) || Utilities::isNotNullPtr(lbR)) {
			g_phi = new double[nV]();

			// (R'*lb_L contribution)
			if (Utilities::isNotNullPtr(lbL)) {
				if (sparseSolver)
					Utilities::AddTransponsedMatrixMultiplication(R_sparse, lbL, g_phi);
				else
					Utilities::AddTransponsedMatrixMultiplication(R, lbL, g_phi, nComp, nV, 1);
			}

			// (L'*lb_R contribution)
			if (Utilities::isNotNullPtr(lbR)) {
				if (sparseSolver)
					Utilities::AddTransponsedMatrixMultiplication(L_sparse, lbR, g_phi);
				else
					Utilities::AddTransponsedMatrixMultiplication(L, lbR, g_phi, nComp, nV, 1);
			}

			// Sign must be negative (really have 0 <= Lx - lbL and 0 <= Rx - lbR)
			for (int i = 0; i < nV; i++)
				g_phi[i] = -g_phi[i];
		}

		// Initialize variables and counters
		alphak = 1;
		rho = options.getInitialPenaltyParameter( );
		outerIter = 0;
		innerIter = 0;
		totalIter = 0;
		algoStat = AlgorithmStatus::PROBLEM_NOT_SOLVED;

		// Set solver options
		if (options.getQPSolver() < QPSolver::OSQP_SPARSE)
			subsolver.setOptions(options.getqpOASESOptions());
		else
			subsolver.setOptions(options.getOSQPOptions());

		// Reset output statistics
		stats.reset();

		// Set seed
		srand( (unsigned int)time( NULL ) );

		// Clean up temporary box constraints
		if (Utilities::isNotNullPtr(lb_tmp)) {
			delete[] lb_tmp;
			lb_tmp = NULL;
		}

		if (Utilities::isNotNullPtr(ub_tmp)) {
			delete[] ub_tmp;
			ub_tmp = NULL;
		}

		// Print new line before printing anything else (might not have been printed by other users...)
		if (options.getPrintLevel() > PrintLevel::NONE)
			printf("\n");

		return ret;
	}


	ReturnValue LCQProblem::switchToSparseMode()
	{

		// Only perform this action if required
		if (sparseSolver){
			return SUCCESSFUL_RETURN;
		}

		Q_sparse = Utilities::dns_to_csc(Q, nV, nV);
		A_sparse = Utilities::dns_to_csc(A, nC + 2*nComp, nV);

		L_sparse = Utilities::dns_to_csc(L, nComp, nV);
		R_sparse = Utilities::dns_to_csc(R, nComp, nV);
		C_sparse = Utilities::dns_to_csc(C, nV, nV);

		// Make sure that all sparse matrices are not null pointer
		if (Utilities::isNullPtr(Q_sparse) || Utilities::isNullPtr(A_sparse) || Utilities::isNullPtr(L_sparse) || Utilities::isNullPtr(R_sparse) || Utilities::isNullPtr(C_sparse)) {
			return FAILED_SWITCH_TO_SPARSE;
		}

		// Clean up dense data (only if succeeded)
		delete[] Q; Q = NULL;
		delete[] A; A = NULL;
		delete[] L; L = NULL;
		delete[] R; R = NULL;
		delete[] C; C = NULL;

		// Toggle sparsity flag
		sparseSolver = true;

		return SUCCESSFUL_RETURN;
	}


	ReturnValue LCQProblem::switchToDenseMode()
	{

		// Only perform this action if required
		if (!sparseSolver){
			return SUCCESSFUL_RETURN;
		}

		Q = Utilities::csc_to_dns(Q_sparse);
		A = Utilities::csc_to_dns(A_sparse);

		L = Utilities::csc_to_dns(L_sparse);
		R = Utilities::csc_to_dns(R_sparse);
		C = Utilities::csc_to_dns(C_sparse);

		// Make sure that all sparse matrices are not null pointer
		if (Utilities::isNullPtr(Q) || Utilities::isNullPtr(A) || Utilities::isNullPtr(L) || Utilities::isNullPtr(R) || Utilities::isNullPtr(C)) {
			return FAILED_SWITCH_TO_DENSE;
		}

		// Clean up sparse data (only if succeeded)
		Utilities::ClearSparseMat(&C_sparse);
		Utilities::ClearSparseMat(&A_sparse);
		Utilities::ClearSparseMat(&Q_sparse);
		Utilities::ClearSparseMat(&L_sparse);
		Utilities::ClearSparseMat(&R_sparse);

		// Toggle sparsity flag
		sparseSolver = false;

		return SUCCESSFUL_RETURN;
	}


	void LCQProblem::updateLinearization()
	{
		if (sparseSolver) {
			Utilities::AffineLinearTransformation(rho, C_sparse, xk, g_tilde, gk, nV);
		} else {
			Utilities::AffineLinearTransformation(rho, C, xk, g_tilde, gk, nV, nV);
		}
	}


	ReturnValue LCQProblem::solveQPSubproblem(bool initialSolve)
	{
		// First solve convex subproblem
		ReturnValue ret = subsolver.solve( initialSolve, qpIterk, qpSolverExitFlag, gk, lbA, ubA, xk, yk, lb, ub );

		// Update stats
		stats.updateSubproblemIter(qpIterk);
		stats.updateQPSolverExitFlag(qpSolverExitFlag);

		// If no initial guess was passed, then need to allocate memory
		if (Utilities::isNullPtr(xk)) {
			xk = new double[nV];
		}

		if (Utilities::isNullPtr(yk)) {
			yk = new double[nDuals];
		}

		// Return on error
		if (ret != SUCCESSFUL_RETURN)
			return ret;

		// Update xnew, yk
		subsolver.getSolution(xnew, yk);

		// Update yk_A
		for (int i = 0; i < nC + 2*nComp; i++)
			yk_A[i] = yk[boxDualOffset + i];

		// Update pk
		Utilities::WeightedVectorAdd(1, xnew, -1, xk, pk, nV);

		return SUCCESSFUL_RETURN;
	}


	bool LCQProblem::stationarityCheck()
  {
		return Utilities::MaxAbs(statk, nV) < options.getStationarityTolerance();
	}


	bool LCQProblem::complementarityCheck()
  {
		return getPhi() < options.getComplementarityTolerance();
	}


	double LCQProblem::getObj()
  {
		double lin = Utilities::DotProduct(g, xk, nV);

		if (sparseSolver) {
			return lin + Utilities::QuadraticFormProduct(Q_sparse, xk, nV)/2.0;
		} else {
			return lin + Utilities::QuadraticFormProduct(Q, xk, nV)/2.0;
		}
	}


	double LCQProblem::getPhi()
  {
		double phi_lin = 0;

		// Linear term
		if (Utilities::isNotNullPtr(g_phi))
			phi_lin += Utilities::DotProduct(g_phi, xk, nV);

		// Quadratic term
		if (sparseSolver) {
			return phi_const + phi_lin + Utilities::QuadraticFormProduct(C_sparse, xk, nV)/2.0;
		} else {
			return phi_const + phi_lin + Utilities::QuadraticFormProduct(C, xk, nV)/2.0;
		}
	}


	double LCQProblem::getMerit()
  {
		double lin = Utilities::DotProduct(g, xk, nV);

		if (sparseSolver) {
			return lin + Utilities::QuadraticFormProduct(Qk_sparse, xk, nV)/2.0;
		} else {
			return lin + Utilities::QuadraticFormProduct(Qk, xk, nV)/2.0;
		}
	}


	void LCQProblem::updatePenalty()
  {
		// Clear Leyffer history
		if (options.getNDynamicPenalty() > 0)
			complHistory.clear();

		rho *= options.getPenaltyUpdateFactor();
		stats.updateRhoOpt( rho );

		// On penalty update also update Qk = Q + rhok C
		updateQk();

		// Update g_tilde = g + rho*g_phi
		if (Utilities::isNotNullPtr(g_phi)) {
			Utilities::WeightedVectorAdd(1.0, g, rho, g_phi, g_tilde, nV);
		}
	}

	void LCQProblem::getOptimalStepLength()
  {

		double qk;

		if (sparseSolver) {
			qk = Utilities::QuadraticFormProduct(Qk_sparse, pk, nV);
			Utilities::AffineLinearTransformation(1, Qk_sparse, xk, g_tilde, lk_tmp, nV);
		} else {
			qk = Utilities::QuadraticFormProduct(Qk, pk, nV);
			Utilities::AffineLinearTransformation(1, Qk, xk, g_tilde, lk_tmp, nV, nV);
		}

		double lk = Utilities::DotProduct(pk, lk_tmp, nV);

		alphak = 1;

		// Convex Descent Case
		if (qk > 0 && lk < 0) {
			alphak = std::min(-lk/qk, 1.0);
		}
	}


	void LCQProblem::updateStep()
  {
		// xk = xk + alphak*pk
		Utilities::WeightedVectorAdd(1, xk, alphak, pk, xk, nV);
	}


	void LCQProblem::updateStationarity()
  {
		// stat = Qk*xk + g - A'*yk_A - yk_x
		// 1) Objective contribution: Qk*xk + g
		if (sparseSolver) {
			Utilities::AffineLinearTransformation(1, Qk_sparse, xk, g_tilde, statk, nV);
		} else {
			Utilities::AffineLinearTransformation(1, Qk, xk, g_tilde, statk, nV, nV);
		}

		// 2) Constraint contribution: A'*yk
		// Utilities::TransponsedMatrixMultiplication(A, yk_A, constr_statk, nC + 2*nComp, nV, 1);
		if (sparseSolver) {
			Utilities::TransponsedMatrixMultiplication(A_sparse, yk_A, constr_statk);
		} else {
			Utilities::TransponsedMatrixMultiplication(A, yk_A, constr_statk, nC + 2*nComp, nV, 1);
		}

		Utilities::WeightedVectorAdd(1, statk, -1, constr_statk, statk, nV);

		// 3) Box constraint contribution
		if (Utilities::isNotNullPtr(lb) || Utilities::isNotNullPtr(ub)) {
			for (int i = 0; i < nV; i++)
				box_statk[i] = yk[i];

			Utilities::WeightedVectorAdd(1, statk, -1, box_statk, statk, nV);
		}
	}


	bool LCQProblem::leyfferCheckPositive()
  {

		size_t n = (size_t)options.getNDynamicPenalty();

		// Only perform Leyffer check if desired
		if (n <= 0)
			return false;

		// Evaluate current complementarity satisfaction
		double complCur = getPhi();

		// Don't perform in first getNDynamicPenalty steps
		if (complHistory.size() < n) {
			complHistory.push_back(complCur);
			return false;
		}

		// Don't increase penalty if already at satisfactory level
		if (complementarityCheck()) {
			complHistory.pop_front();
			complHistory.push_back(complCur);
			return false;
		}

		bool retFlag = true;
		for (size_t i = 0; i < n; i++) {
			if (complCur < options.getEtaDynamicPenalty()*complHistory[i]) {
				// In this case phi(xkj) < eta*max{phi(xkj-1),...,phi(xkj-n)}
				retFlag = false;
				break;
			}
		}

		// Update history vector
		complHistory.pop_front();
		complHistory.push_back(complCur);

		return retFlag;
	}


	void LCQProblem::updateQk()
  {
		// Smart update in sparse case
		if (sparseSolver) {
			double factor = rho*(1 - 1.0/options.getPenaltyUpdateFactor());
			for (size_t j = 0; j < Qk_indices_of_C.size(); j++) {
				Qk_sparse->x[Qk_indices_of_C[j]] += factor*C_sparse->x[j];
			}
		} else {
			Utilities::WeightedMatrixAdd(1, Q, rho, C, Qk, nV, nV);
		}
	}


	void LCQProblem::updateOuterIter()
  {
		outerIter++;
		stats.updateIterOuter(1);
	}


	void LCQProblem::updateTotalIter()
  {
		totalIter++;
		stats.updateIterTotal(1);
	}


	void LCQProblem::perturbGradient()
  {

		int randNum;
		for (int i = 0; i < nV; i++) {
			// Random number -1, 0, 1
			randNum = (rand() % 3) - 1;

			gk[i] += randNum*Utilities::EPS;
		}
	}


	void LCQProblem::perturbStep()
  {

		int randNum;
		for (int i = 0; i < nV; i++) {
			// Random number -1, 0, 1
			randNum = (rand() % 3) - 1;

			xk[i] += randNum*Utilities::EPS;
		}
	}


	void LCQProblem::storeSteps()
  {
		stats.updateTrackingVectors(
			xk,
			innerIter,
			qpIterk,
			alphak,
			Utilities::MaxAbs(pk, nV),
			Utilities::MaxAbs(statk, nV),
			getObj(),
			getPhi(),
			getMerit(),
			nV
		);
	}


	void LCQProblem::transformDuals()
  {

		double* tmp = new double[nComp];

		// y_L = y - rho*R*xk
		if (sparseSolver) {
			Utilities::MatrixMultiplication(R_sparse, xk, tmp);
		} else {
			Utilities::MatrixMultiplication(R, xk, tmp, nComp, nV, 1);
		}

		for (int i = 0; i < nComp; i++) {
			yk[boxDualOffset + nC + i] = yk[boxDualOffset + nC + i] - rho*tmp[i];
		}

		// y_R = y - rho*L*xk
		if (sparseSolver) {
			Utilities::MatrixMultiplication(L_sparse, xk, tmp);
		} else {
			Utilities::MatrixMultiplication(L, xk, tmp, nComp, nV, 1);
		}

		for (int i = 0; i < nComp; i++) {
			yk[boxDualOffset + nC + nComp + i] = yk[boxDualOffset + nC + nComp + i] - rho*tmp[i];
		}

		// clear memory
		delete[] tmp;
	}


	void LCQProblem::determineStationarityType()
  {

		std::vector<int> weakComp = getWeakComplementarities( );

		bool s_stationary = true;
		bool m_stationary = true;

		for (int i : weakComp) {
			double dualProd = yk_A[nC + i]*yk_A[nC + nComp + i];
			double dualMin = std::min(yk_A[nC + i], yk_A[nC + nComp + i]);

			// Check failure of s-stationarity
			if (dualMin < 0)
				s_stationary = false;

			// Check failure of m-/c-stationarity
			if (std::abs(dualProd) >= options.getComplementarityTolerance() && dualMin <= 0) {

				// Check failure of c-stationarity
				if (dualProd <= options.getComplementarityTolerance()) {
					algoStat = AlgorithmStatus::W_STATIONARY_SOLUTION;
					return;
				}

				m_stationary = false;
			}
		}

		if (s_stationary) {
			algoStat = AlgorithmStatus::S_STATIONARY_SOLUTION;
			return;
		}


		if (m_stationary) {
			algoStat = AlgorithmStatus::M_STATIONARY_SOLUTION;
			return;
		}

		algoStat = AlgorithmStatus::C_STATIONARY_SOLUTION;
		return;
	}


	std::vector<int> LCQProblem::getWeakComplementarities()
	{
		double* Lx = new double[nComp];
		double* Rx = new double[nComp];

		if (sparseSolver) {
			Utilities::MatrixMultiplication(L_sparse, xk, Lx);
			Utilities::MatrixMultiplication(R_sparse, xk, Rx);
		} else {
			Utilities::MatrixMultiplication(L, xk, Lx, nComp, nV, 1);
			Utilities::MatrixMultiplication(R, xk, Rx, nComp, nV, 1);
		}

		std::vector<int> indices;

		for (int i = 0; i < nComp; i++) {
			if (Lx[i] <= options.getComplementarityTolerance())
				if (Rx[i] <= options.getComplementarityTolerance())
					indices.push_back(i);
		}

		// Free memory
		delete[] Lx;
		delete[] Rx;

		return indices;
	}


	AlgorithmStatus LCQProblem::getPrimalSolution(double* const xOpt) const
	{
		if (Utilities::isNotNullPtr(xOpt) && Utilities::isNotNullPtr(xk)) {
			for (int i = 0; i < nV; i++)
				xOpt[i] = xk[i];
		}

		return algoStat;
	}


	AlgorithmStatus LCQProblem::getDualSolution(double* const yOpt) const
	{
		if (Utilities::isNotNullPtr(yOpt) && Utilities::isNotNullPtr(yk)) {
			for (int i = 0; i < nDuals; i++)
				yOpt[i] = yk[i];
		}

		return algoStat;
	}


	int LCQProblem::getNumberOfPrimals() const
	{
		return nV;
	}


	int LCQProblem::getNumberOfDuals() const
	{
		return nDuals;
	}

int LCQProblem::getNV() const
{
  return nV;
}

int LCQProblem::getNC() const
{
  return nC;
}
int LCQProblem::getNComp() const
{
  return nComp;
}


	void LCQProblem::getOutputStatistics(OutputStatistics& _stats) const
	{
		_stats = stats;
	}


	/*
	 *	 p r i n t I t e r a t i o n
	 */
	void LCQProblem::printIteration()
	{
		if (options.getPrintLevel() == PrintLevel::NONE)
			return;

		if (options.getPrintLevel() == PrintLevel::OUTER_LOOP_ITERATES && innerIter > 0)
			return;

		// Print header every 10 iters
		bool headerInner = (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES && innerIter % 10 == 0);
		bool headerOuter = (options.getPrintLevel() == PrintLevel::OUTER_LOOP_ITERATES && outerIter % 10 == 0);

		if (headerInner || headerOuter)
			printHeader();

		const char* sep = " | ";

		// Print outer iterate
		printf("%6d", outerIter);

		// Print innter iterate
		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES)
			printf("%s%*d", sep, 6, innerIter);

		// Print stationarity violation
		double tmpdbl = Utilities::MaxAbs(statk, nV);
		printf("%s%10.3g", sep, tmpdbl);

		// Print complementarity violation
		printf("%s%10.3g", sep, getPhi());

		// Print current penalty parameter
		printf("%s%10.3g", sep, rho);

		// Print infinity norm of computed full step
		tmpdbl = Utilities::MaxAbs(pk, nV);
		printf("%s%10.3g", sep, tmpdbl);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES) {
			// Print optimal step length
			printf("%s%10.3g", sep, alphak);

			// Print number of qpOASES iterations
			printf("%s%6d", sep, qpIterk);
		}

		// Print new line
		printf(" \n");
	}


	void LCQProblem::printHeader()
	{
		printLine();

		const char* sep = " | ";
		const char* outer = " outer";
		const char* inner = " inner";
		const char* stat = "  station ";
		const char* comp = "  complem ";
		const char* pen = "    rho   ";
		const char* np = "  norm p  ";
		const char* sl = "   alpha  ";
		const char* subIt = "sub it";

		printf("%s",outer);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES)
			printf("%s%s", sep, inner);

		printf("%s%s", sep, stat);
		printf("%s%s", sep, comp);
		printf("%s%s", sep, pen);
		printf("%s%s", sep, np);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES) {
			printf("%s%s%s%s", sep, sl, sep, subIt);
		}

		printf(" \n");

		printLine();
	}

	/*
	 * 	 p r i n t L i n e
	 */
	void LCQProblem::printLine()
	{
		const char* iSep = "------";
		const char* dSep = "----------";
		const char* node = "-+-";

		printf("%s", iSep);

		// Print innter iterate
		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES)
			printf("%s%s", node, iSep);

		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);
		printf("%s%s", node, dSep);

		if (options.getPrintLevel() >= PrintLevel::INNER_LOOP_ITERATES) {
			printf("%s%s%s%s", node, dSep, node, iSep);
		}

		printf("-\n");
	}


	/// Clear allocated memory
	void LCQProblem::clear()
	{
		if (Utilities::isNotNullPtr(Q)) {
			delete[] Q;
			Q = NULL;
		}

		if (Utilities::isNotNullPtr(g)) {
			delete[] g;
			g = NULL;
		}

		if (Utilities::isNotNullPtr(lb)) {
			delete[] lb;
			lb = NULL;
		}

		if (Utilities::isNotNullPtr(ub)) {
			delete[] ub;
			ub = NULL;
		}

		if (Utilities::isNotNullPtr(lb_tmp)) {
			delete[] lb_tmp;
			lb_tmp = NULL;
		}

		if (Utilities::isNotNullPtr(ub_tmp)) {
			delete[] ub_tmp;
			ub_tmp = NULL;
		}

		if (Utilities::isNotNullPtr(A)) {
			delete[] A;
			A = NULL;
		}

		if (Utilities::isNotNullPtr(lbA)) {
			delete[] lbA;
			lbA = NULL;
		}

		if (Utilities::isNotNullPtr(ubA)) {
			delete[] ubA;
			ubA = NULL;
		}

		if (Utilities::isNotNullPtr(L)) {
			delete[] L;
			L = NULL;
		}

		if (Utilities::isNotNullPtr(R)) {
			delete[] R;
			R = NULL;
		}

		if (Utilities::isNotNullPtr(C)) {
			delete[] C;
			C = NULL;
		}

		if (Utilities::isNotNullPtr(lbL)) {
			delete[] lbL;
			lbL = NULL;
		}

		if (Utilities::isNotNullPtr(ubL)) {
			delete[] ubL;
			ubL = NULL;
		}

		if (Utilities::isNotNullPtr(lbR)) {
			delete[] lbR;
			lbR = NULL;
		}

		if (Utilities::isNotNullPtr(ubR)) {
			delete[] ubR;
			ubR = NULL;
		}

		if (Utilities::isNotNullPtr(g_phi)) {
			delete[] g_phi;
			g_phi = NULL;
		}

		if (Utilities::isNotNullPtr(gk)) {
			delete[] gk;
			gk = NULL;
		}

		if (Utilities::isNotNullPtr(g_tilde)) {
			delete[] g_tilde;
			g_tilde = NULL;
		}

		if (Utilities::isNotNullPtr(xk)) {
			delete[] xk;
			xk = NULL;
		}

		if (Utilities::isNotNullPtr(yk)) {
			delete[] yk;
			yk = NULL;
		}

		if (Utilities::isNotNullPtr(yk_A)) {
			delete[] yk_A;
			yk_A = NULL;
		}

		if (Utilities::isNotNullPtr(xnew)) {
			delete[] xnew;
			xnew = NULL;
		}

		if (Utilities::isNotNullPtr(pk)) {
			delete[] pk;
			pk = NULL;
		}

		if (Utilities::isNotNullPtr(Qk)) {
			delete[] Qk;
			Qk = NULL;
		}

		if (Utilities::isNotNullPtr(statk)) {
			delete[] statk;
			statk = NULL;
		}

		if (Utilities::isNotNullPtr(constr_statk)) {
			delete[] constr_statk;
			constr_statk = NULL;
		}

		if (Utilities::isNotNullPtr(box_statk)) {
			delete[] box_statk;
			box_statk = NULL;
		}

		if (Utilities::isNotNullPtr(lk_tmp)) {
			delete[] lk_tmp;
			lk_tmp = NULL;
		}

		Utilities::ClearSparseMat(&C_sparse);
		Utilities::ClearSparseMat(&A_sparse);
		Utilities::ClearSparseMat(&Q_sparse);
		Utilities::ClearSparseMat(&Qk_sparse);
		Utilities::ClearSparseMat(&L_sparse);
		Utilities::ClearSparseMat(&R_sparse);
	}
}





/*
 *	end of file
 */
