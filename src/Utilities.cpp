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

#include "Utilities.hpp"
#include "MessageHandler.hpp"

#include <iostream>
#include <vector>

#include <qpOASES.hpp>

extern "C"
{
#include <osqp.h>
}

namespace LCQPow
{

void Utilities::MatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p)
{
  for(int i = 0; i < m; i++)
  {
    for(int j = 0; j < p; j++)
    {
      C[i * p + j] = 0;
      for(int k = 0; k < n; k++)
      {
        C[i * p + j] += A[i * n + k] * B[k * p + j];
      }
    }
  }
}

void Utilities::MatrixMultiplication(const csc* const A, const double* const b, double* c)
{
  for(int i = 0; i < A->m; i++)
    c[i] = 0;

  for(int j = 0; j < A->n; j++)
  {
    for(int i = A->p[j]; i < A->p[j + 1]; i++)
    {
      c[A->i[i]] += A->x[i] * b[j];
    }
  }
}

void Utilities::TransponsedMatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n,
                                                int p)
{
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < p; j++)
    {
      C[i * p + j] = 0;

      for(int k = 0; k < m; k++)
      {
        C[i * p + j] += A[k * n + i] * B[k * p + j];
      }
    }
  }
}

void Utilities::TransponsedMatrixMultiplication(const csc* const A, const double* const b, double* c)
{
  for(int j = 0; j < A->n; j++)
  {
    c[j] = 0;
    for(int k = A->p[j]; k < A->p[j + 1]; k++)
    {
      c[j] += b[A->i[k]] * A->x[k];
    }
  }
}

void Utilities::AddTransponsedMatrixMultiplication(const double* const A, const double* const B, double* C, int m,
                                                   int n, int p)
{
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < p; j++)
    {
      for(int k = 0; k < m; k++)
      {
        C[i * p + j] += A[k * n + i] * B[k * p + j];
      }
    }
  }
}

void Utilities::AddTransponsedMatrixMultiplication(const csc* const A, const double* const b, double* c)
{
  for(int j = 0; j < A->n; j++)
  {
    for(int k = A->p[j]; k < A->p[j + 1]; k++)
    {
      c[j] += b[A->i[k]] * A->x[k];
    }
  }
}

void Utilities::MatrixSymmetrizationProduct(const double* const A, const double* const B, double* C, int m, int n)
{
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j <= i; j++)
    {
      C[i * n + j] = 0;
      for(int k = 0; k < m; k++)
      {
        C[i * n + j] += A[k * n + i] * B[k * n + j] + B[k * n + i] * A[k * n + j];
      }

      // Make symmetric
      C[j * n + i] = C[i * n + j];
    }
  }
}

csc* Utilities::MatrixSymmetrizationProduct(double* L_x, int* L_i, int* L_p, double* R_x, int* R_i, int* R_p, int n)
{
  std::vector<int> C_rows;
  std::vector<double> C_data;
  int* C_p = (int*) malloc((size_t) (n + 1) * sizeof(int));
  C_p[0] = 0;

  for(int j = 0; j < n; j++)
  {
    C_p[j + 1] = C_p[j];
    for(int i = 0; i < n; i++)
    {
      double tmp = 0;

      // (L'*R)_ij
      for(int k = L_p[i]; k < L_p[i + 1]; k++)
      {
        int ind_s2 = getIndexOfIn(L_i[k], R_i, R_p[j], R_p[j + 1]);
        if(ind_s2 != -1)
        {
          tmp += L_x[k] * R_x[ind_s2];
        }
      }

      // (R'*L)_ij
      for(int k = R_p[i]; k < R_p[i + 1]; k++)
      {
        int ind_s1 = getIndexOfIn(R_i[k], L_i, L_p[j], L_p[j + 1]);

        if(ind_s1 != -1)
        {
          tmp += R_x[k] * L_x[ind_s1];
        }
      }

      // If the entry is non-zero append it to the data
      if(!isZero(tmp))
      {
        C_rows.push_back(i);
        C_data.push_back(tmp);
        C_p[j + 1]++;
      }
    }
  }

  if(C_p[n] == 0)
    return 0;

  int* C_i = (int*) malloc((size_t) C_p[n] * sizeof(int));
  double* C_x = (double*) malloc((size_t) C_p[n] * sizeof(double));

  for(int i = 0; i < C_p[n]; i++)
  {
    C_i[i] = C_rows[(size_t) i];
    C_x[i] = C_data[(size_t) i];
  }

  csc* M = createCSC(n, n, C_p[n], C_x, C_i, C_p);
  return M;
}

csc* Utilities::MatrixSymmetrizationProduct(csc* L, csc* R)
{
  return MatrixSymmetrizationProduct(L->x, L->i, L->p, R->x, R->i, R->p, L->n);
}

void Utilities::AffineLinearTransformation(const double alpha, const double* const A, const double* const b,
                                           const double* const c, double* d, int m, int n)
{
  for(int i = 0; i < m; i++)
  {

    double tmp = 0;
    for(int k = 0; k < n; k++)
    {
      tmp += A[i * n + k] * b[k];
    }

    d[i] = alpha * tmp + c[i];
  }
}

void Utilities::AffineLinearTransformation(const double alpha, const csc* const S, const double* const b,
                                           const double* const c, double* d, int m)
{
  for(int j = 0; j < m; j++)
  {

    double tmp = 0;
    for(int k = S->p[j]; k < S->p[j + 1]; k++)
    {
      tmp += S->x[k] * b[S->i[k]];
    }

    d[j] = alpha * tmp + c[j];
  }
}

void Utilities::WeightedMatrixAdd(const double alpha, const double* const A, const double beta, const double* const B,
                                  double* C, int m, int n)
{
  for(int i = 0; i < m; i++)
    for(int j = 0; j < n; j++)
      C[i * n + j] = alpha * A[i * n + j] + beta * B[i * n + j];
}

void Utilities::WeightedVectorAdd(const double alpha, const double* const a, const double beta, const double* const b,
                                  double* c, int m)
{
  WeightedMatrixAdd(alpha, a, beta, b, c, m, 1);
}

double Utilities::QuadraticFormProduct(const double* const Q, const double* const p, int m)
{
  double ret = 0;
  for(int i = 0; i < m; i++)
  {
    double tmp = 0;
    for(int j = 0; j < m; j++)
      tmp += Q[i * m + j] * p[j];

    ret += tmp * p[i];
  }

  return ret;
}

double Utilities::QuadraticFormProduct(const csc* const S, const double* const p, int m)
{
  double ret = 0;
  for(int j = 0; j < m; j++)
  {

    double tmp = 0;
    for(int k = S->p[j]; k < S->p[j + 1]; k++)
    {
      tmp += S->x[k] * p[S->i[k]];
    }

    ret += p[j] * tmp;
  }

  return ret;
}

double Utilities::DotProduct(const double* const a, const double* const b, int m)
{
  double ret = 0;
  for(int i = 0; i < m; i++)
    ret += a[i] * b[i];

  return ret;
}

double Utilities::MaxAbs(const double* const a, int m)
{
  double mx = 0;
  double mn = 0;

  for(int i = 0; i < m; i++)
  {
    if(a[i] > mx)
      mx = a[i];
    else if(a[i] < mn)
      mn = a[i];
  }

  return Utilities::getMax(mx, -mn);
}

void Utilities::ClearSparseMat(csc* M)
{
  if(isNotNullPtr(M))
  {
    if(isNotNullPtr(M->p))
    {
      free(M->p);
      M->p = NULL;
    }
    if(isNotNullPtr(M->i))
    {
      free(M->i);
      M->i = NULL;
    }
    if(isNotNullPtr(M->x))
    {
      free(M->x);
      M->x = NULL;
    }

    free(M);
    M = NULL;
  }
}

void Utilities::ClearSparseMatCPP(csc* M)
{
  if(M != nullptr)
  {
    if(M->p != nullptr)
    {
      delete M->p;
      M->p = nullptr;
    }
    if(M->i != nullptr)
    {
      delete M->i;
      M->i = nullptr;
    }
    if(M->x != nullptr)
    {
      delete M->x;
      M->x = nullptr;
    }

    delete M;
    M = nullptr;
  }
}

void Utilities::ClearSparseMat(csc** M)
{
  if(isNotNullPtr(*M))
  {
    if(isNotNullPtr((*M)->p))
    {
      free((*M)->p);
      (*M)->p = NULL;
    }
    if(isNotNullPtr((*M)->i))
    {
      free((*M)->i);
      (*M)->i = NULL;
    }
    if(isNotNullPtr((*M)->x))
    {
      free((*M)->x);
      (*M)->x = NULL;
    }

    free(*M);
    *M = NULL;
  }
}

ReturnValue Utilities::readFromFile(int* data, int n, const char* datafilename)
{
  int i;
  FILE* datafile;

  /* 1) Open file. */
  if((datafile = fopen(datafilename, "r")) == 0)
  {
    fclose(datafile);
    return UNABLE_TO_READ_FILE;
  }

  /* 2) Read data from file. */
  for(i = 0; i < n; ++i)
  {
    if(fscanf(datafile, "%d\n", &(data[i])) == 0)
    {
      fclose(datafile);
      return UNABLE_TO_READ_FILE;
    }
  }

  /* 3) Close file. */
  fclose(datafile);

  return SUCCESSFUL_RETURN;
}

ReturnValue Utilities::readFromFile(double* data, int n, const char* datafilename)
{
  int i;
  FILE* datafile;

  /* 1) Open file. */
  if((datafile = fopen(datafilename, "r")) == 0)
  {
    return UNABLE_TO_READ_FILE;
  }

  /* 2) Read data from file. */
  for(i = 0; i < n; ++i)
  {
    if(fscanf(datafile, "%lf\n", &(data[i])) == 0)
    {
      fclose(datafile);
      return UNABLE_TO_READ_FILE;
    }
  }

  /* 3) Close file. */
  fclose(datafile);

  return SUCCESSFUL_RETURN;
}

ReturnValue Utilities::writeToFile(double* data, int n, const char* datafilename)
{
  int i;
  FILE* datafile;

  /* 1) Open file. */
  if((datafile = fopen(datafilename, "w")) == 0)
  {
    fclose(datafile);
    return UNABLE_TO_READ_FILE;
  }

  /* 2) Read data from file. */
  for(i = 0; i < n; ++i)
  {
    if(fprintf(datafile, "%f\n", data[i]) == 0)
    {
      fclose(datafile);
      return UNABLE_TO_READ_FILE;
    }
  }

  /* 3) Close file. */
  fclose(datafile);

  return SUCCESSFUL_RETURN;
}

void Utilities::printMatrix(const double* const A, int m, int n, const char* const name)
{
  printf("Printing matrix %s:\n", name);

  for(int i = 0; i < m; i++)
  {
    for(int j = 0; j < n; j++)
      printf("%.5f ", A[i * n + j]);

    printf("\n");
  }

  printf("\n");
}

void Utilities::printMatrix(const int* const A, int m, int n, const char* const name)
{
  printf("Printing matrix %s:\n", name);

  for(int i = 0; i < m; i++)
  {
    for(int j = 0; j < n; j++)
      printf("%d ", A[i * n + j]);

    printf("\n");
  }

  printf("\n");
  fflush(stdout);
}

void Utilities::printMatrix(const csc* A, const char* const name)
{
  if(A->m * A->n <= 0)
    return;

  // Get dense representation
  double* dense = Utilities::csc_to_dns(A);

  // Print the dense matrix
  Utilities::printMatrix(dense, A->m, A->n, name);

  // Clear memory
  delete[] dense;
}

void Utilities::printStep(double* xk, double* pk, double* xk_new, double alpha, int nV)
{
  printf("Printing Step:\n");

  for(int i = 0; i < nV; i++)
    printf("%.2f + %.2f * %.2f = %.2f \n", xk[i], alpha, pk[i], xk_new[i]);

  printf("\n");
}

void Utilities::printBounds(double* lb, double* xk, double* ub, int m)
{
  printf("Printing box constraints:\n");

  for(int i = 0; i < m; i++)
    printf("%.2f <= %.2f <= %.2f \n", lb[i], xk[i], ub[i]);

  printf("\n");
}

csc* Utilities::createCSC(int m, int n, int nnx, double* x, int* i, int* p)
{
  csc* M = (csc*) malloc(sizeof(csc));

  if(isNullPtr(M))
    return 0;

  M->m = m;
  M->n = n;
  M->p = p;
  M->i = i;
  M->x = x;
  M->nz = -1;
  M->nzmax = nnx;

  return M;
}

csc* Utilities::copyCSC(int m, int n, int nnx, double* x, int* i, int* p)
{
  csc* M = (csc*) malloc(sizeof(csc));

  if(isNullPtr(M))
    return 0;

  // Allocate space
  int* rows = (int*) malloc((size_t) nnx * sizeof(int));
  double* data = (double*) malloc((size_t) nnx * sizeof(double));
  int* cols = (int*) malloc((size_t) (n + 1) * sizeof(int));

  // Copy sparse matrix data
  memcpy(rows, i, (size_t) nnx * sizeof(int));
  memcpy(data, x, (size_t) nnx * sizeof(double));
  memcpy(cols, p, (size_t) (n + 1) * sizeof(int));

  // Assign copied data
  M->m = m;
  M->n = n;
  M->p = cols;
  M->i = rows;
  M->x = data;
  M->nz = -1;
  M->nzmax = nnx;

  return M;
}

csc* Utilities::copyCSC(const csc* const _M, bool toUpperTriangular)
{
  csc* M = (csc*) malloc(sizeof(csc));

  if(isNullPtr(M))
    return 0;

  if(toUpperTriangular)
  {
    // Allocate space
    std::vector<int> rows;
    std::vector<double> data;
    int* p = (int*) malloc((size_t) (_M->n + 1) * sizeof(int));

    p[0] = 0;

    for(int j = 0; j < _M->n; j++)
    {
      p[j + 1] = p[j];

      for(int i = _M->p[j]; i < _M->p[j + 1]; i++)
      {
        // Ignore entries below diagonal
        if(_M->i[i] > j)
          continue;

        // Push back elements on or above diagonal
        rows.push_back(_M->i[i]);
        data.push_back(_M->x[i]);
        p[j + 1]++;
      }
    }

    // copy std vector to arrays
    int* i = (int*) malloc((size_t) p[_M->n] * sizeof(int));
    double* x = (double*) malloc((size_t) p[_M->n] * sizeof(double));
    for(size_t k = 0; k < (size_t) p[_M->n]; k++)
    {
      i[k] = rows[k];
      x[k] = data[k];
    }

    // Assign copied data
    M->m = _M->m;
    M->n = _M->n;
    M->p = p;
    M->i = i;
    M->x = x;
    M->nz = -1;
    M->nzmax = p[_M->n];
  }
  else
  {
    // Allocate space
    int* rows = (int*) malloc((size_t) _M->nzmax * sizeof(int));
    double* data = (double*) malloc((size_t) _M->nzmax * sizeof(double));
    int* cols = (int*) malloc((size_t) (_M->n + 1) * sizeof(int));

    // Copy sparse matrix data
    memcpy(rows, _M->i, (size_t) _M->nzmax * sizeof(int));
    memcpy(data, _M->x, (size_t) _M->nzmax * sizeof(double));
    memcpy(cols, _M->p, (size_t) (_M->n + 1) * sizeof(int));

    // Assign copied data
    M->m = _M->m;
    M->n = _M->n;
    M->p = cols;
    M->i = rows;
    M->x = data;
    M->nz = -1;
    M->nzmax = _M->nzmax;
  }

  return M;
}

void Utilities::copyIntToIntT(int* dest, const int* const src, int n)
{
  for(int i = 0; i < n; i++)
    dest[i] = src[i];
}

double* Utilities::csc_to_dns(const csc* const sparse)
{
  int m = sparse->m;
  int n = sparse->n;
  double* full = new double[m * n]();

  for(int j = 0; j < n; j++)
  {
    for(int i = sparse->p[j]; i < sparse->p[j + 1]; i++)
    {
      // Reached final element
      if(i == sparse->nzmax)
      {
        return full;
      }

      // Ensure validity of index
      if(sparse->i[i] * n + j >= m * n || sparse->i[i] * n + j < 0)
      {
        MessageHandler::PrintMessage(INDEX_OUT_OF_BOUNDS, ERROR);
        return 0;
      }

      full[sparse->i[i] * n + j] = sparse->x[i];
    }
  }

  return full;
}

csc* Utilities::dns_to_csc(const double* const full, int m, int n)
{
  std::vector<double> Q_data;
  std::vector<int> Q_rows;
  int* Q_p = (int*) malloc((size_t) (n + 1) * sizeof(int));
  Q_p[0] = 0;

  for(int i = 0; i < n; i++)
  {
    // Begin column pointer with previous value
    Q_p[i + 1] = Q_p[i];

    for(int j = 0; j < m; j++)
    {
      if(full[j * n + i] > 0 || full[j * n + i] < 0)
      {
        Q_data.push_back(full[j * n + i]);
        Q_rows.push_back(j);
        Q_p[i + 1]++;
      }
    }
  }

  int* Q_i = (int*) malloc((size_t) Q_p[n] * sizeof(int));
  double* Q_x = (double*) malloc((size_t) Q_p[n] * sizeof(double));

  for(int i = 0; i < Q_p[n]; i++)
  {
    Q_i[i] = Q_rows[(size_t) i];
    Q_x[i] = Q_data[(size_t) i];
  }

  csc* sparse = createCSC(m, n, Q_p[n], Q_x, Q_i, Q_p);
  return sparse;
}

double Utilities::getAbs(double x)
{
#ifdef __NO_FMATH__
  return (x >= 0.0) ? x : -x;
#else
  return fabs(x);
#endif
}

bool Utilities::isEqual(double x, double y, double TOL)
{
  if(getAbs(x - y) <= TOL)
    return true;
  else
    return false;
}

bool Utilities::isEqual(double x, double y) { return isEqual(x, y, Utilities::ZERO); }

bool Utilities::isZero(double x, double TOL)
{
  if(getAbs(x) <= TOL)
    return true;
  else
    return false;
}

bool Utilities::isZero(double x) { return isZero(x, Utilities::ZERO); }

double Utilities::getSign(double arg)
{
  if(arg >= 0.0)
    return 1.0;
  else
    return -1.0;
}

int Utilities::getMax(int x, int y) { return (y < x) ? x : y; }

int Utilities::getMin(int x, int y) { return (y > x) ? x : y; }

double Utilities::getMax(double x, double y) { return (y < x) ? x : y; }

double Utilities::getMin(double x, double y) { return (y > x) ? x : y; }

int Utilities::getIndexOfIn(int val, int* sorted_lst, int beg, int end)
{
  for(int i = beg; i < end; i++)
  {
    if(sorted_lst[i] == val)
      return i;

    if(sorted_lst[i] > val)
      break;
  }

  return -1;
}
} // namespace LCQPow
