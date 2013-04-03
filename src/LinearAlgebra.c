// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Copyright - Maurice Berk (maurice@mauriceberk.com) (http://www2.imperial.ac.uk/~mab201)

#include <R.h>
#include <R_ext/BLAS.h>

#include "LinearAlgebra.h"
#include "Matrix.h"
#include "Vector.h"

extern void F77_NAME(dgesv)(int*, int*, double*, int*, int*, double*, int*, int*);
extern void F77_NAME(dpotf2)(char*, int*, double*, int*, int*);
extern void F77_NAME(dgels)(char*, int*, int*, int*, double*, int*, double*, int*, double*, int*, int*);
extern void F77_NAME(dgesvd)(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
extern void F77_NAME(dgetrf)(int*, int*, double*, int*, int*, int*);

void invertMatrix(Matrix* matrix, Matrix* result)
{
  int row, column, info;
  double* matrixMemoryCopy = calloc(matrix->rows * matrix->columns, sizeof(double));
  int* ipiv = calloc(matrix->rows, sizeof(int));
  memcpy(matrixMemoryCopy, matrix->pointer, sizeof(double) * matrix->rows * matrix->columns);

  for(row = 0; row < result->rows; row++)
  {
    for(column = 0; column < result->columns; column++)
    {
      result->pointer[column * result->rows + row] = column == row;
    }
  }

  F77_NAME(dgesv)(&matrix->rows,
                  &matrix->rows,
                  matrixMemoryCopy,
                  &matrix->rows,
                  ipiv,
                  result->pointer,
                  &result->rows,
                  &info);

  free(ipiv);
  free(matrixMemoryCopy);
}

void multiplyMatrices(Matrix* A, Matrix* B, Matrix* result, char isATransposed, char isBTransposed)
{
  double one = 1.0;
  double zero = 0.0;
  
  char transA = isATransposed ? 'T' : 'N';
  char transB = isBTransposed ? 'T' : 'N';

  //In case the result is A or B
  double* resultCopy = calloc(result->rows * result->columns, sizeof(double));

  F77_NAME(dgemm)(&transA,
                  &transB,
                  isATransposed ? &A->columns : &A->rows,
                  isBTransposed ? &B->rows : &B->columns,
                  isATransposed ? &A->rows : &A->columns,
                  &one,
                  A->pointer,
                  &A->rows,
                  B->pointer,
                  &B->rows,
                  &zero,
                  resultCopy,
                  &result->rows);

  memcpy(result->pointer, resultCopy, sizeof(double) * result->rows * result->columns);
  free(resultCopy);
}

void matrixVectorMultiply(Matrix* A, Vector* x, Vector* y, Vector* result, double alpha, double beta, char isATransposed)
{
  char trans = isATransposed ? 'T' : 'N';
  int one = 1;
  double* resultCopy = calloc(result->length, sizeof(double));

  if(beta)
  {
    memcpy(resultCopy, y->pointer, y->length * sizeof(double));
  }
  
  F77_NAME(dgemv)(&trans,
                  &A->rows,
                  &A->columns,
                  &alpha,
                  A->pointer,
                  &A->rows,
                  x->pointer,
                  &one,
                  &beta,
                  resultCopy,
                  &one);

  memcpy(result->pointer, resultCopy, sizeof(double) * result->length);
  free(resultCopy);
}

void symmetricRank1Update(Matrix* A, Vector* x, double alpha)
{
  char uplo = 'U';
  int one = 1;
  int i, j;

  F77_NAME(dsyr)(&uplo,
                 &A->rows,
                 &alpha,
                 x->pointer,
                 &one,
                 A->pointer,
                 &A->rows);

  for(i = 0; i < A->rows; i++)
  {
    for(j = i; j < A->columns; j++)
    {
      A->pointer[j + i * A->rows] = A->pointer[i + j * A->rows];
    }
  }
}

void symmetricRankKUpdate(Matrix* A, Matrix* C, double alpha, double beta)
{
  char uplo = 'U';
  char trans = 'N';
  int i, j;

  F77_NAME(dsyrk)(&uplo,
                  &trans,
                  &C->rows,
                  &A->columns,
                  &alpha,
                  A->pointer,
                  &A->rows,
                  &beta,
                  C->pointer,
                  &C->rows);

  //Copy the upper triangle of C into the lower triangle
  for(i = 0; i < C->rows; i++)
  {
    for(j = i; j < C->columns; j++)
    {
      C->pointer[j + i * C->rows] = C->pointer[i + j * C->rows];
    }
  }
}

void choleskyFactorization(Matrix* matrix, Matrix* result)
{
  char uplo = 'U';
  int info;
  int i, j;

  memcpy(result->pointer, matrix->pointer, sizeof(double) * matrix->rows * matrix->columns);

  F77_NAME(dpotf2)(&uplo,
                   &result->rows,
                   result->pointer,
                   &result->rows,
                   &info);

  //Not sure if it's necessary to zero the lower triangle but do it just in case
  for(i = 0; i < result->rows; i++)
  {
    for(j = 0; j < result->columns; j++)
    {
      if(i > j) result->pointer[i + j * result->rows] = 0.0;
    }
  }
}

void leastSquares(Matrix* X, Vector* y, Vector* beta)
{
  char trans = 'N';
  int one = 1;

  double* workspace = calloc(1, sizeof(double));
  int workspaceLength = -1;
  int info;

  double* XMemoryCopy = calloc(X->rows * X->columns, sizeof(double));
  double* yMemoryCopy = calloc(y->length, sizeof(double));

  memcpy(XMemoryCopy, X->pointer, sizeof(double) * X->rows * X->columns);
  memcpy(yMemoryCopy, y->pointer, sizeof(double) * y->length);

  // First make a workspace query by setting workspaceLength to -1
  F77_NAME(dgels)(&trans,
                  &X->rows,
                  &X->columns,
                  &one,
                  XMemoryCopy,
                  &X->rows,
                  yMemoryCopy,
                  &y->length,
                  workspace,
                  &workspaceLength,
                  &info);

  // Now free up the original workspace and reallocate the right amount
  workspaceLength = workspace[0];

  free(workspace);
  workspace = calloc(workspaceLength, sizeof(double));
  F77_NAME(dgels)(&trans,
                  &X->rows,
                  &X->columns,
                  &one,
                  XMemoryCopy,
                  &X->rows,
                  yMemoryCopy,
                  &y->length,
                  workspace,
                  &workspaceLength,
                  &info);

  //Finally, copy the result into beta->pointer
  memcpy(beta->pointer, yMemoryCopy, sizeof(double) * X->columns);

  free(workspace);
  free(yMemoryCopy);
  free(XMemoryCopy);
}

void penalisedLeastSquares(Matrix* X, Vector* y, Vector* beta, Matrix* roughnessMatrix, double* lambda)
{
  /*Matrix factorization;
  double* roughnessMatrixMemoryCopy = calloc(roughnessMatrix->rows * roughnessMatrix->columns, double);
  

  factorization.rows = roughnessMatrix.rows;
  factorization.*/
}

// Basically turn it into a standard least squares problem by using data augmentation
void penalisedLeastSquaresUsingFactorization(Matrix* X, Vector* y, Vector* beta, Matrix* factorization)
{
  int i, j;

  // Augmented y is just y with some trailing zeros
  // Augmented X is trickier to construct
  Vector augmentedY;
  Matrix augmentedX;
  
  augmentedY.length = y->length + factorization->rows;

  augmentedX.rows = X->rows + factorization->rows;
  augmentedX.columns = X->columns;
  
  augmentedX.pointer = calloc(augmentedX.rows * augmentedX.columns, sizeof(double));

  augmentedY.pointer = calloc(augmentedY.length, sizeof(double));

  memcpy(augmentedY.pointer, y->pointer, sizeof(double) * y->length);

  for(i = 0; i < X->rows; i++)
  {
    for(j = 0; j < X->columns; j++)
    {
      augmentedX.pointer[i + j * augmentedX.rows] = X->pointer[i + j * X->rows];
    }
  }
  for(i = X->rows; i < augmentedX.rows; i++)
  {
    for(j = 0; j < X->columns; j++)
    {
      augmentedX.pointer[i + j * augmentedX.rows] = factorization->pointer[(i - X->rows) + j * factorization->rows];
    }
  }

  leastSquares(&augmentedX, &augmentedY, beta);

  free(augmentedX.pointer);
  free(augmentedY.pointer);
}

void logDeterminant(Matrix* A, double* logDeterminant)
{
  int info;
  double* matrixMemoryCopy = calloc(A->rows * A->columns, sizeof(double));
  int* ipiv = calloc(A->rows, sizeof(int));
  double sign = 1.0;
  int i;

  memcpy(matrixMemoryCopy, A->pointer, sizeof(double) * A->rows * A->columns);

  F77_NAME(dgetrf)(&A->rows,
                   &A->rows,
                   matrixMemoryCopy,
                   &A->rows,
                   ipiv,
                   &info);

  *logDeterminant = 0.0;
  for(i = 0; i < A->rows; i++)
  {
    if(ipiv[i] != (i+1))
    {
      sign = -sign;
    }
  }
  for(i = 0; i < A->rows; i++)
  {
    if(matrixMemoryCopy[i + i * A->rows] < 0)
    {
      *logDeterminant += log(-matrixMemoryCopy[i + i * A->rows]);
      sign = -sign;
    }
    else
    {
      *logDeterminant += log(matrixMemoryCopy[i + i * A->rows]);
    }
  }

  *logDeterminant *= sign;

  free(ipiv);
  free(matrixMemoryCopy);
}

void singularValueDecomposition(Matrix* A, Matrix* U, Matrix* Vtransposed, Matrix* Sigma)
{
  double* matrixMemoryCopy = calloc(A->rows * A->columns, sizeof(double));
  double* singularValues = calloc(Sigma->rows, sizeof(double));
  double optimalWorkLength;
  double* work;
  int workLength = -1;
  int info;

  int i;

  char jobu = 'A';
  char jobvt = 'A';

  memcpy(matrixMemoryCopy, A->pointer, sizeof(double) * A->rows * A->columns);

  //First workspace query
  F77_NAME(dgesvd)(&jobu,
                   &jobvt,
                   &A->rows,
                   &A->columns,
                   matrixMemoryCopy,
                   &A->rows,
                   singularValues,
                   U->pointer,
                   &U->rows,
                   Vtransposed->pointer,
                   &Vtransposed->rows,
                   &optimalWorkLength,
                   &workLength,
                   &info);

  workLength = (int) optimalWorkLength;
  work = calloc(workLength, sizeof(double));

  //Now do it for real
  F77_NAME(dgesvd)(&jobu,
                   &jobvt,
                   &A->rows,
                   &A->columns,
                   matrixMemoryCopy,
                   &A->rows,
                   singularValues,
                   U->pointer,
                   &U->rows,
                   Vtransposed->pointer,
                   &Vtransposed->rows,
                   work,
                   &workLength,
                   &info);

  for(i = 0; i < Sigma->rows; i++)
  {
    Sigma->pointer[i + i * Sigma->rows] = singularValues[i];
  }

  free(work);
  free(singularValues);
  free(matrixMemoryCopy);
}

void matrixSquareRoot(Matrix* A, Matrix* S)
{
  //Approach is to do an SVD of A, then set S to V * sqrt(Sigma) * Vtransposed
  int i, j;

  Matrix* U = allocateMatrix(A->rows, A->columns);
  Matrix* Vtransposed = allocateMatrix(A->rows, A->columns);
  Matrix* Sigma = allocateMatrix(A->rows, A->columns);

  //First do an SVD
  singularValueDecomposition(A, U, Vtransposed, Sigma);

  //Now calculate sqrt(Sigma) * Vtransposed and stick it in U
  for(i = 0; i < Sigma->rows; i++)
  {
    for(j = 0; j < Sigma->columns; j++)
    {
      U->pointer[i + j * U->rows] = Vtransposed->pointer[i + j * Vtransposed->rows] * sqrt(Sigma->pointer[i + i * Sigma->rows]);
    }
  }

  multiplyMatrices(Vtransposed, U, S, 1, 0);

  freeMatrix(U);
  freeMatrix(Vtransposed);
  freeMatrix(Sigma);
}