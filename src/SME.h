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

void SMEMultiple(
  int* M,  //Number of variables being fit
  double* y,  //All observations across all variables
  double* timePoints, //All time points across all variables
  int* individual,  //All individual identifiers across all variables
  double* X,  //M Ni*p basis matrices, stored one after the other
  int* Ni,  //M-length vector Observations per variable
  int* ni,  //M-length vector Individuals per variable
  int* Nij, //Observations per individual per variable
  int* p, //M-length vector Number of unique time points per variable
  double* lambdaMu, //Output - optimal lambda mu per variable
  double* lambdaV,  //Output - optimal lambda V per variable
  double* G,  //M p*p matrices, stored one after the other
  double* mu, //Output - M p length vectors of fitted coefficients
  double* sigmaSquared, //Output - M error variances
  double* D,  //Output - M p*p covariance matrices
  double* v,  //Output - M ni*p vectors of random effects for all subjects
  double* likelihood, //Output - M likelihoods
  double* dfMu, //Output - M degrees of freedom for mu
  double* dfV,  //Output - M degrees of freedom for V
  int* iterations,  //Output - number of iterations it took for final EM fit per variable
  int* maxIterations, //1 control parameters indicating maximum iterations for EM per variable
  double* deltaEM, //1 control parameters indicating convergence criteria for EM
  double* deltaNM, //1 control parameters indicating convergence criteria for NM
  int* verbose, //1 integer indicating if the entire procedure should be verbose
  int* info, //M length integers returning status of fits
  int* numberOfThreads);  //Number of threads to use

void SMEOptimizationMultiple(
  int* M,  //Number of variables being fit
  double* y,  //All observations across all variables
  double* timePoints, //All time points across all variables
  int* individual,  //All individual identifiers across all variables
  double* X,  //M Ni*p basis matrices, stored one after the other
  int* Ni,  //M-length vector Observations per variable
  int* ni,  //M-length vector Individuals per variable
  int* Nij, //Observations per individual per variable
  int* p, //M-length vector Number of unique time points per variable
  double* lambdaMu, //Output - optimal lambda mu per variable
  double* lambdaV,  //Output - optimal lambda V per variable
  double* G,  //M p*p matrices, stored one after the other
  double* mu, //Output - M p length vectors of fitted coefficients
  double* sigmaSquared, //Output - M error variances
  double* D,  //Output - M p*p covariance matrices
  double* v,  //Output - M ni*p vectors of random effects for all subjects
  double* likelihood, //Output - M likelihoods
  double* dfMu, //Output - M degrees of freedom for mu
  double* dfV,  //Output - M degrees of freedom for V
  int* iterations,  //Output - number of iterations it took for final EM fit per variable
  int* maxIterations, //1 control parameters indicating maximum iterations for EM per variable
  double* deltaEM, //1 control parameters indicating convergence criteria for EM
  double* deltaNM, //1 control parameters indicating convergence criteria for NM
  int* criteria,  //1 control parameter indicating which smoothing parameter selection criteria to use
  int* verbose, //1 integer indicating if the entire procedure should be verbose
  int* info, //M length integers returning status of fits
  int* numberOfThreads);  //Number of threads to use

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
