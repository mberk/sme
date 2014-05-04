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

#ifndef SME_H
#define SME_H

#include "Matrix.h"
#include "Vector.h"

void SMEOptimization(
         double* y,
         double* timePoints,
         int* individual,
         double* X, //N * p basis matrix
         int* N, //Total number of observations
         int* n, //Total number of individuals
         int* Ni, //n-length vector of number of observations per individual
         int* p, //Number of unique time points
         double* lambdaMu,
         double* lambdaV,
         double* G, //p * p roughness matrix
         double* mu, //Output: p-length vector of fitted mean curve values
         double* sigmaSquared, //Output: error variance
         double* D, //Output: p * p covariance matrix of the random-effects
         double* v, // p * n length vector of random effects for all subjects
         double* likelihood,
         double* dfMu,
         double* dfV,
         int* iterations,
         int* maxIterations,
         double* deltaEM,
         double* deltaNM,
         int* criteria,
         int* verbose,
         int* info);

void SME(double* y,
         double* timePoints,
         int* individual,
         double* X, //N * p basis matrix
         int* N, //Total number of observations
         int* n, //Total number of individuals
         int* Ni, //n-length vector of number of observations per individual
         int* p, //Number of unique time points
         double* lambdaMu,
         double* lambdaV,
         double* G, //p * p roughness matrix
         double* mu, //Output: p-length vector of fitted mean curve values
         double* sigmaSquared, //Output: error variance
         double* D, //Output: p * p covariance matrix of the random-effects
         double* v, // p * n length vector of random effects for all subjects
         double* likelihood,
         double* dfMu,
         double* dfV,
         int* iterations,
         int* maxIterations,
         double* deltaEM,
         int* verbose,
         int* info);

void calculateYiPrecision(Matrix** Xi, Matrix* Dv, double* sigmaSquared, int n, Matrix** inverseVi);

void EStep(Vector** yi, Matrix** Xi, Matrix** inverseVi, int N, int n, Matrix* Dv, Vector* mu, Vector** vi, Vector** epsiloni);
void MStep(Vector** yi,
           Matrix* X,
           Matrix** Xi,
           Vector** vi,
           Vector** epsiloni,
           Matrix** inverseVi,
           Vector* mu,
           Matrix* muPenaltyFactorization,
           Matrix* Dv,
           Matrix* D,
           Matrix* G,
           double lambdaV,
           double* sigmaSquared,
           int N,
           int n);

void calculateLikelihood(Vector** yi, Matrix** Xi, Matrix** inverseVi, Vector* mu, int n, double *likelihood);
void calculateDegreesOfFreedom(Matrix* X, Matrix** Xi, Matrix** inverseVi, Matrix* G, Matrix* Dv, double lambdaMu, int n, double* dfMu, double* dfV);

#endif