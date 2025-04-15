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
 *	but WITQOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "SubsolverQPOASES.hpp"

#include <qpOASES.hpp>

namespace LCQPow {


    SubsolverQPOASES::SubsolverQPOASES( ) { }


    SubsolverQPOASES::SubsolverQPOASES( int _nV, int _nC,
                                        double* _Q, double* _A)
    {
        nV = _nV;
        nC = _nC;

        isSparse = false;
        qp = qpOASES::QProblem(nV, nC);

        Q = new double[nV*nV];
        A = new double[nC*nV];

        memcpy(Q, _Q, (size_t)(nV*nV)*sizeof(double));
        memcpy(A, _A, (size_t)(nC*nV)*sizeof(double));
    }


    SubsolverQPOASES::SubsolverQPOASES( int _nV, int _nC,
                                        csc* _Q, csc* _A)
    {
        nV = _nV;
        nC = _nC;

        isSparse = true;

        // Determine whether to use Schur complement method
        useSchur = false;
        #ifdef SOLVER_MA57
            useSchur = true;
        #endif

        if (useSchur) {
            qpSchur = qpOASES::SQProblemSchur(nV, nC);
        } else {
            qp = qpOASES::QProblem(nV, nC);
        }

        if (Utilities::isNotNullPtr(Q_sparse)) {
            free(Q_sparse);
            Q_sparse = NULL;
        }

        if (Utilities::isNotNullPtr(A_sparse)) {
            free(A_sparse);
            A_sparse = NULL;
        }

        Q_i = new int[_Q->p[_nV]];
        Q_x = new double[_Q->p[_nV]];
        Q_p = new int[_nV+1];

        A_i = new int[_A->p[_nV]];
        A_x = new double[_A->p[_nV]];
        A_p = new int[_nV+1];

        Utilities::copyIntToIntT(Q_p, _Q->p, nV+1);
        Utilities::copyIntToIntT(Q_i, _Q->i, _Q->p[nV]);
        Utilities::copyIntToIntT(A_p, _A->p, nV+1);
        Utilities::copyIntToIntT(A_i, _A->i, _A->p[nV]);

        memcpy(Q_x, _Q->x, (size_t)_Q->p[nV]*sizeof(double));
        memcpy(A_x, _A->x, (size_t)_A->p[nV]*sizeof(double));

        Q_sparse = new qpOASES::SymSparseMat(nV, nV, Q_i, Q_p, Q_x);
        A_sparse = new qpOASES::SparseMatrix(nC, nV, A_i, A_p, A_x);
        Q_sparse->createDiagInfo();
        A_sparse->createDiagInfo();
    }


    SubsolverQPOASES::SubsolverQPOASES(const SubsolverQPOASES& rhs)
    {
        copy( rhs );
    }


    SubsolverQPOASES::~SubsolverQPOASES()
    {
        clear();
    }


    SubsolverQPOASES& SubsolverQPOASES::operator=(const SubsolverQPOASES& rhs)
    {
        if (this != &rhs) {
            copy( rhs );
        }

        return *this;
    }


    void SubsolverQPOASES::setOptions( qpOASES::Options options )
    {
        if (useSchur) {
            qpSchur.setOptions( options );
        } else {
            qp.setOptions( options );
        }
    }


    ReturnValue SubsolverQPOASES::solve(    bool initialSolve, int& iterations, int& exit_flag,
                                            const double* const g,
                                            const double* const lbA, const double* const ubA,
                                            const double* const x0, const double* const y0,
                                            const double* const lb, const double* const ub )
    {
        qpOASES::returnValue ret;

        int nwsr = 1000000;

        if (initialSolve) {
            if (isSparse) {
                if (useSchur) {
                    ret = qpSchur.init(Q_sparse, g, A_sparse, lb, ub, lbA, ubA, nwsr, (double*)0, x0, y0);
                } else {
                    ret = qp.init(Q_sparse, g, A_sparse, lb, ub, lbA, ubA, nwsr, (double*)0, x0, y0);
                }
            } else {
                ret = qp.init(Q, g, A, lb, ub, lbA, ubA, nwsr, (double*)0, x0, y0);
            }
        } else {
            if (useSchur) {
                ret = qpSchur.hotstart(g, lb, ub, lbA, ubA, nwsr);
            } else {
                ret = qp.hotstart(g, lb, ub, lbA, ubA, nwsr);
            }
        }

        iterations = (int)(nwsr);
        exit_flag = (int)(ret);

        if (ret != qpOASES::returnValue::SUCCESSFUL_RETURN)
            return ReturnValue::SUBPROBLEM_SOLVER_ERROR;

        return ReturnValue::SUCCESSFUL_RETURN;
    }


    void SubsolverQPOASES::getSolution( double* x, double* y )
    {
        if (useSchur) {
            qpSchur.getPrimalSolution( x );
            qpSchur.getDualSolution( y );
        } else {
            qp.getPrimalSolution( x );
            qp.getDualSolution( y );
        }
    }


    void SubsolverQPOASES::copy(const SubsolverQPOASES& rhs)
    {

        clear();

        nV = rhs.nV;
        nC = rhs.nC;

        isSparse = rhs.isSparse;
        useSchur = rhs.useSchur;

        if (isSparse) {
            Q_i = new int[rhs.Q_p[nV]];
            Q_x = new double[rhs.Q_p[nV]];
            Q_p = new int[nV+1];

            A_i = new int[rhs.A_p[nV]];
            A_x = new double[rhs.A_p[nV]];
            A_p = new int[nV+1];


            memcpy(Q_p, rhs.Q_p, (size_t)(nV+1)*sizeof(int));
            memcpy(Q_i, rhs.Q_i, (size_t)(rhs.Q_p[nV])*sizeof(int));
            memcpy(Q_x, rhs.Q_x, (size_t)(rhs.Q_p[nV])*sizeof(double));
            memcpy(A_p, rhs.A_p, (size_t)(nV+1)*sizeof(int));
            memcpy(A_i, rhs.A_i, (size_t)(rhs.A_p[nV])*sizeof(int));
            memcpy(A_x, rhs.A_x, (size_t)(rhs.A_p[nV])*sizeof(double));

            Q_sparse = new qpOASES::SymSparseMat(nV, nV, Q_i, Q_p, Q_x);
            A_sparse = new qpOASES::SparseMatrix(nC, nV, A_i, A_p, A_x);

            Q_sparse->createDiagInfo();
            A_sparse->createDiagInfo();
        } else {
            Q = new double[nV*nV];
            A = new double[nC*nV];

            memcpy(Q, rhs.Q, (size_t)(nV*nV)*sizeof(double));
            memcpy(A, rhs.A, (size_t)(nC*nV)*sizeof(double));
        }

        if (useSchur) {
            qpSchur = rhs.qpSchur;
        } else {
            qp = rhs.qp;
        }
    }


    void SubsolverQPOASES::clear() 
    {
        if (Utilities::isNotNullPtr(Q)) {
            delete[] Q;
            Q = NULL;
        }

        if (Utilities::isNotNullPtr(A)) {
            delete[] A;
            A = NULL;
        }

        if (Utilities::isNotNullPtr(Q_i)) {
            delete[] Q_i;
            Q_i = NULL;
        }

        if (Utilities::isNotNullPtr(Q_p)) {
            delete[] Q_p;;
            Q_p = NULL;
        }

        if (Utilities::isNotNullPtr(Q_x)) {
            delete[] Q_x;
            Q_x = NULL;
        }

        if (Utilities::isNotNullPtr(A_x)) {
            delete[] A_x;
            A_x = NULL;
        }

        if (Utilities::isNotNullPtr(A_p)) {
            delete[] A_p;
            A_p = NULL;
        }

        if (Utilities::isNotNullPtr(A_x)) {
            delete[] A_x;
            A_x = NULL;
        }

        if (Utilities::isNotNullPtr(Q_sparse)) {
            delete Q_sparse;
            Q_sparse = 0;
        }

        if (Utilities::isNotNullPtr(A_sparse)) {
            delete A_sparse;
            A_sparse = 0;
        }
    }
}
