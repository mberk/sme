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

#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include "Matrix.h"
#include "Vector.h"

void invertMatrix(Matrix* matrix, Matrix* result);
void multiplyMatrices(Matrix* A, Matrix* B, Matrix* result, char isATransposed, char isBTransposed);

void matrixVectorMultiply(Matrix* A, Vector* x, Vector* y, Vector* result, double alpha, double beta, char isATransposed);

void symmetricRankKUpdate(Matrix* A, Matrix* C, double alpha, double beta);

void symmetricRank1Update(Matrix* A, Vector* x, double alpha);

void leastSquares(Matrix* X, Vector* y, Vector* beta);
void penalisedLeastSquares(Matrix* X, Vector* y, Vector* beta, Matrix* roughnessMatrix, double* lambda);

void choleskyFactorization(Matrix* matrix, Matrix* result);

void penalisedLeastSquaresUsingFactorization(Matrix* X, Vector* y, Vector* beta, Matrix* factorization);

void logDeterminant(Matrix* A, double* logDeterminant);

void singularValueDecomposition(Matrix* A, Matrix* U, Matrix* Vtransposed, Matrix* Sigma);

void matrixSquareRoot(Matrix* A, Matrix* S);

#endif