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

// Copyright - Maurice Berk (maurice.berk01@imperial.ac.uk) (http://www2.imperial.ac.uk/~mab201)

#include <R.h>
#include "NelderMead.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "SME.h"
#include "Vector.h"
#include "Utility.h"
#include "LinearAlgebra.h"
#include "BlockDiagonalMatrix.h"

#define INFO_LIKELIHOOD_DECREASED -1

#define CRITERIA_AIC  1
#define CRITERIA_AICc 2
#define CRITERIA_BICN 3
#define CRITERIA_BICn 4

typedef struct
{
  double* y;
  double* timePoints;
  int* individual;
  double* X;
  int* N;
  int* n;
  int* Ni;
  int* p;
  double* G;
  double* mu;
  double* sigmaSquared;
  double* D;
  double* v;
  double* likelihood;
  double* dfMu;
  double* dfV;
  int* zeroIntercept;
  int* iterations;
  int* maxIterations;
  double* deltaEM;
  int* verbose;
  int* info;
  int* criteria;
} SMEParameters;

double SMEWrapper(
  int n,
  double* par,
  void *ex)
{
  SMEParameters* smeParameters = (SMEParameters*) ex;

  double lambdaMu = abs(par[0]);//exp(par[0]);
  double lambdaV = abs(par[1]);//exp(par[1]);

  SME(smeParameters->y,
      smeParameters->timePoints,
      smeParameters->individual,
      smeParameters->X,
      smeParameters->N,
      smeParameters->n,
      smeParameters->Ni,
      smeParameters->p,
      &lambdaMu,
      &lambdaV,
      smeParameters->G,
      smeParameters->mu,
      smeParameters->sigmaSquared,
      smeParameters->D,
      smeParameters->v,
      smeParameters->likelihood,
      smeParameters->dfMu,
      smeParameters->dfV,
      smeParameters->zeroIntercept,
      smeParameters->iterations,
      smeParameters->maxIterations,
      smeParameters->deltaEM,
      smeParameters->verbose,
      smeParameters->info);

  double score;
  if(*smeParameters->info)
  {
    //EM algorithm failed for these values of smoothing parameters
    score = 1e35;
  }
  else if(*smeParameters->dfMu < 2.0 || *smeParameters->dfV < 0.0)
  {
    //Smoothing parameters were too large and resulted in numerical instability
    score = 1e35;
  }
  else
  {
    if(*smeParameters->criteria == CRITERIA_AIC)
    {
      score = -2.0 * *smeParameters->likelihood + 2.0 * (*smeParameters->dfMu + *smeParameters->dfV);
    }
    else if(*smeParameters->criteria == CRITERIA_BICN)
    {
      score = -2.0 * *smeParameters->likelihood + log(*smeParameters->N) + 2.0 * (*smeParameters->dfMu + *smeParameters->dfV);
    }
    else if(*smeParameters->criteria == CRITERIA_BICn)
    {
      score = -2.0 * *smeParameters->likelihood + log(*smeParameters->n) + 2.0 * (*smeParameters->dfMu + *smeParameters->dfV);
    }
    else// if(*smeParameters->criteria == CRITERIA_AICc)
    {
      //For now we will treat any other value as corrected AIC. The R code should be responsible for capturing invalid arguments
      score = -2.0 * *smeParameters->likelihood + 2.0 * *smeParameters->N * (*smeParameters->dfMu + *smeParameters->dfV) / (*smeParameters->N - (*smeParameters->dfMu + *smeParameters->dfV) - 1.0);
    }
  }
  
  if(*smeParameters->verbose)
  {
    Rprintf("%f,%f (%f,%f) (%f,%f) gave a score of %f in %d iterations\n", par[0], par[1], lambdaMu, lambdaV, *smeParameters->dfMu, *smeParameters->dfV, score, *smeParameters->iterations);
  }
  
  return score;
}

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
  int* zeroIntercept, //M control parameters indicating if variable should be fit with zero intercept
  int* iterations,  //Output - number of iterations it took for final EM fit per variable
  int* maxIterations, //1 control parameters indicating maximum iterations for EM per variable
  double* deltaEM, //1 control parameters indicating convergence criteria for EM
  double* deltaNM, //1 control parameters indicating convergence criteria for NM
  int* criteria,  //1 control parameter indicating which smoothing parameter selection criteria to use
  int* verbose, //1 integer indicating if the entire procedure should be verbose
  int* info, //M length integers returning status of fits
  int* numberOfThreads)  //Number of threads to use
{
  int i;

#ifdef _OPENMP
  if(*numberOfThreads == -1)
  {
    omp_set_dynamic(1);
  }
  else
  {
    omp_set_dynamic(0);
    omp_set_num_threads(*numberOfThreads);
  }
#endif

  SMEParameters* allParameters = calloc(*M, sizeof(SMEParameters));
  
  double* currentY = y;
  double* currentTimePoints = timePoints;
  int* currentIndividual = individual;
  double* currentX = X;
  int* currentNij = Nij;
  double* currentG = G;
  double* currentMu = mu;
  double* currentD = D;
  double* currentV = v;
  for(i = 0; i < *M; i++)
  {
    allParameters[i].y = currentY;
    allParameters[i].timePoints = currentTimePoints;
    allParameters[i].individual = currentIndividual;
    allParameters[i].X = currentX;
    allParameters[i].N = &Ni[i];
    allParameters[i].n = &ni[i];
    allParameters[i].Ni = currentNij;
    allParameters[i].p = &p[i];
    allParameters[i].G = currentG;
    allParameters[i].mu = currentMu;
    allParameters[i].sigmaSquared = &sigmaSquared[i];
    allParameters[i].D = currentD;
    allParameters[i].v = currentV;
    allParameters[i].likelihood = &likelihood[i];
    allParameters[i].dfMu = &dfMu[i];
    allParameters[i].dfV = &dfV[i];
    allParameters[i].zeroIntercept = &zeroIntercept[i];
    allParameters[i].iterations = &iterations[i];
    allParameters[i].maxIterations = maxIterations;
    allParameters[i].deltaEM = deltaEM;
    allParameters[i].verbose = verbose;
    allParameters[i].info = &info[i];
    allParameters[i].criteria = criteria;

    currentY += Ni[i];
    currentTimePoints += Ni[i];
    currentIndividual += Ni[i];
    currentX += (Ni[i] * p[i]);
    currentNij += ni[i];
    currentG += (p[i] * p[i]);
    currentMu += p[i];
    currentD += (p[i] * p[i]);
    currentV += (ni[i] * p[i]);
  }

  #pragma omp parallel for
  for(i = 0; i < *M; i++)
  {
    SMEOptimization(
         allParameters[i].y, 
         allParameters[i].timePoints,
         allParameters[i].individual,
         allParameters[i].X, //N * p basis matrix
         allParameters[i].N, //Total number of observations
         allParameters[i].n, //Total number of individuals
         allParameters[i].Ni, //n-length vector of number of observations per individual
         allParameters[i].p, //Number of unique time points
         &lambdaMu[i],
         &lambdaV[i],
         allParameters[i].G, //p * p roughness matrix
         allParameters[i].mu, //Output: p-length vector of fitted mean curve values
         allParameters[i].sigmaSquared, //Output: error variance
         allParameters[i].D, //Output: p * p covariance matrix of the random-effects
         allParameters[i].v, // p * n length vector of random effects for all subjects
         allParameters[i].likelihood,
         allParameters[i].dfMu,
         allParameters[i].dfV,
         allParameters[i].zeroIntercept,
         allParameters[i].iterations,
         allParameters[i].maxIterations,
         deltaEM,
         deltaNM,
         criteria,
         allParameters[i].verbose,
         allParameters[i].info);
  }

  free(allParameters);
}

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
  int* zeroIntercept, //M control parameters indicating if variable should be fit with zero intercept
  int* iterations,  //Output - number of iterations it took for final EM fit per variable
  int* maxIterations, //1 control parameters indicating maximum iterations for EM per variable
  double* deltaEM, //1 control parameters indicating convergence criteria for EM
  double* deltaNM, //1 control parameters indicating convergence criteria for NM
  int* verbose, //1 integer indicating if the entire procedure should be verbose
  int* info, //M length integers returning status of fits
  int* numberOfThreads)  //Number of threads to use
{
  int i;

  if(*numberOfThreads == -1)
  {
    omp_set_dynamic(1);
  }
  else
  {
    omp_set_dynamic(0);
    omp_set_num_threads(*numberOfThreads);
  }

  SMEParameters* allParameters = calloc(*M, sizeof(SMEParameters));
  
  double* currentY = y;
  double* currentTimePoints = timePoints;
  int* currentIndividual = individual;
  double* currentX = X;
  int* currentNij = Nij;
  double* currentG = G;
  double* currentMu = mu;
  double* currentD = D;
  double* currentV = v;
  for(i = 0; i < *M; i++)
  {
    allParameters[i].y = currentY;
    allParameters[i].timePoints = currentTimePoints;
    allParameters[i].individual = currentIndividual;
    allParameters[i].X = currentX;
    allParameters[i].N = &Ni[i];
    allParameters[i].n = &ni[i];
    allParameters[i].Ni = currentNij;
    allParameters[i].p = &p[i];
    allParameters[i].G = currentG;
    allParameters[i].mu = currentMu;
    allParameters[i].sigmaSquared = &sigmaSquared[i];
    allParameters[i].D = currentD;
    allParameters[i].v = currentV;
    allParameters[i].likelihood = &likelihood[i];
    allParameters[i].dfMu = &dfMu[i];
    allParameters[i].dfV = &dfV[i];
    allParameters[i].zeroIntercept = &zeroIntercept[i];
    allParameters[i].iterations = &iterations[i];
    allParameters[i].maxIterations = maxIterations;
    allParameters[i].deltaEM = deltaEM;
    allParameters[i].verbose = verbose;
    allParameters[i].info = &info[i];

    currentY += Ni[i];
    currentTimePoints += Ni[i];
    currentIndividual += Ni[i];
    currentX += (Ni[i] * p[i]);
    currentNij += ni[i];
    currentG += (p[i] * p[i]);
    currentMu += p[i];
    currentD += (p[i] * p[i]);
    currentV += (ni[i] * p[i]);
  }

  #pragma omp parallel for
  for(i = 0; i < *M; i++)
  {
    SME(allParameters[i].y,
        allParameters[i].timePoints,
        allParameters[i].individual,
        allParameters[i].X, //N * p basis matrix
        allParameters[i].N, //Total number of observations
        allParameters[i].n, //Total number of individuals
        allParameters[i].Ni, //n-length vector of number of observations per individual
        allParameters[i].p, //Number of unique time points
        &lambdaMu[i],
        &lambdaV[i],
        allParameters[i].G, //p * p roughness matrix
        allParameters[i].mu, //Output: p-length vector of fitted mean curve values
        allParameters[i].sigmaSquared, //Output: error variance
        allParameters[i].D, //Output: p * p covariance matrix of the random-effects
        allParameters[i].v, // p * n length vector of random effects for all subjects
        allParameters[i].likelihood,
        allParameters[i].dfMu,
        allParameters[i].dfV,
        allParameters[i].zeroIntercept,
        allParameters[i].iterations,
        allParameters[i].maxIterations,
        allParameters[i].deltaEM,
        allParameters[i].verbose,
        allParameters[i].info);
  }

  free(allParameters);
}


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
         int* zeroIntercept,
         int* iterations,
         int* maxIterations,
         double* deltaEM,
         double* deltaNM,
         int* criteria,
         int* verbose,
         int* info)
{
  double lambdas[2] = {10000.0, 10000.0};
  double minimum;
  int numberOfParameters = 2;
  int fail = 0;
  int maxRetries = 3;
  int numberOfRetries;

  double absoluteTolerance = -INFINITY;
  double relativeTolerance = *deltaNM;

  double alpha = 1.0;
  double beta = 0.5;
  double gamma = 2.0;

  int maxNMIterations = 500;
  int numberOfFunctionCalls = 0;

  SMEParameters smeParameters;
  
  smeParameters.y = y;
  smeParameters.timePoints = timePoints;
  smeParameters.individual = individual;
  smeParameters.X = X;
  smeParameters.N = N;
  smeParameters.n = n;
  smeParameters.Ni = Ni;
  smeParameters.p = p;
  smeParameters.G = G;
  smeParameters.mu = mu;
  smeParameters.sigmaSquared = sigmaSquared;
  smeParameters.D = D;
  smeParameters.v = v;
  smeParameters.likelihood = likelihood;
  smeParameters.dfMu = dfMu;
  smeParameters.dfV = dfV;
  smeParameters.zeroIntercept = zeroIntercept;
  smeParameters.iterations = iterations;
  smeParameters.maxIterations = maxIterations;
  smeParameters.deltaEM = deltaEM;
  smeParameters.verbose = verbose;
  smeParameters.info = info;
  smeParameters.criteria = criteria;

  for(numberOfRetries = 0, *smeParameters.info = 1; *smeParameters.info && numberOfRetries < maxRetries; lambdas[0] *= 10.0, lambdas[1] *= 10.0, numberOfRetries++)
  {
    NelderMead(numberOfParameters,
          lambdas,
          lambdas,
          &minimum,
          SMEWrapper,
          &fail,
          absoluteTolerance,
          relativeTolerance,
          &smeParameters,
          alpha,
          beta,
          gamma,
          0,/**verbose ? 10 : 0,*/
          &numberOfFunctionCalls,
          maxNMIterations);

    *lambdaMu = abs(lambdas[0]);//exp(lambdas[0]);
    *lambdaV = abs(lambdas[1]);//exp(lambdas[1]);
    
    if(*verbose) Rprintf("NM chose %f,%f (%f,%f) which gave minimum %f with %03d runs of the EM aglorithm\n", lambdas[0], lambdas[1], *lambdaMu, *lambdaV, minimum, numberOfFunctionCalls);

    //Refit with optimal lambdas
    SME(y, 
        timePoints,
        individual,
        X, //N * p basis matrix
        N, //Total number of observations
        n, //Total number of individuals
        Ni, //n-length vector of number of observations per individual
        p, //Number of unique time points
        lambdaMu,
        lambdaV,
        G, //p * p roughness matrix
        mu, //Output: p-length vector of fitted mean curve values
        sigmaSquared, //Output: error variance
        D, //Output: p * p covariance matrix of the random-effects
        v, // p * n length vector of random effects for all subjects
        likelihood,
        dfMu,
        dfV,
        zeroIntercept,
        iterations,
        maxIterations,
        deltaEM,
        verbose,
        info);
  }

  double AICc;
  if(*smeParameters.info)
  {
    //EM algorithm failed for these values of smoothing parameters
    AICc = 1e35;
  }
  else
  {
    AICc = -2.0 * *smeParameters.likelihood + 2.0 * *smeParameters.N * (*smeParameters.dfMu + *smeParameters.dfV) / (*smeParameters.N - (*smeParameters.dfMu + *smeParameters.dfV) - 1.0);
  }

  if(*verbose)
  {
    Rprintf("Final info = %d, final AICc = %f\n", *info, AICc);
  }
}

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
         int* zeroIntercept,
         int* iterations,
         int* maxIterations,
         double* deltaEM,
         int* verbose,
         int* info)
{
  int i, j;

  double oldLikelihood = -INFINITY;

  *info = 0;

  Vector yAsVector;
  Vector timePointsAsVector;
  Vector muAsVector;

  Matrix XAsMatrix;
  Matrix muPenalty;
  Matrix muPenaltyFactorization;
  Matrix DAsMatrix;
  Matrix Dv;
  Matrix GAsMatrix;
  Matrix inverseV;

  int numberOfVectors, numberOfMatrices;

  yAsVector.pointer = y;
  yAsVector.length = *N;

  timePointsAsVector.pointer = timePoints;
  timePointsAsVector.length = *N;

  muAsVector.pointer = mu;
  muAsVector.length = *p;
  
  XAsMatrix.pointer = X;
  XAsMatrix.rows = *N;
  XAsMatrix.columns = *p;

  GAsMatrix.pointer = G;
  GAsMatrix.rows= *p;
  GAsMatrix.columns = *p;
  
  DAsMatrix.pointer = D;
  DAsMatrix.rows = *p;
  DAsMatrix.columns = *p;
  
  Dv.pointer = calloc(*p * *p, sizeof(double));
  Dv.rows = *p;
  Dv.columns = *p;

  inverseV.pointer = calloc(*N * *N, sizeof(double));
  inverseV.rows = *N;
  inverseV.columns = *N;

  Vector** yi;
  Vector** timei;
  Matrix** Xi;

  yi = splitVector(&yAsVector, individual, &numberOfVectors);

  timei = splitVector(&timePointsAsVector, individual, &numberOfVectors);

  Xi = splitMatrix(&XAsMatrix, individual, &numberOfMatrices);
  //Matrix* Z = constructBlockDiagonalMatrix(Zi, numberOfMatrices);
  Matrix** inverseVi = calloc(*n, sizeof(Matrix*));

  Vector** epsiloni = calloc(*n, sizeof(Vector*));
  Vector** vi = calloc(*n, sizeof(Vector*));

  for(i = 0; i < *n; i++)
  {
    inverseVi[i] = calloc(1, sizeof(Matrix));
    inverseVi[i]->pointer = calloc(Ni[i] * Ni[i], sizeof(double));
    inverseVi[i]->rows = Ni[i];
    inverseVi[i]->columns = Ni[i];
    
    epsiloni[i] = calloc(1, sizeof(Vector));
    epsiloni[i]->length = Ni[i];
    epsiloni[i]->pointer = calloc(Ni[i], sizeof(double));

    vi[i] = calloc(1, sizeof(Vector));
    vi[i]->pointer = &v[i * *p];
    vi[i]->length = *p;
  }
  
  muPenalty.pointer = calloc(*p * *p, sizeof(double));
  muPenalty.rows = *p;
  muPenalty.columns = *p;
  memcpy(muPenalty.pointer, G, *p * *p * sizeof(double));

  for(i = 0; i < *p * *p; i++)
  {
    muPenalty.pointer[i] *= *lambdaMu;
  }

  muPenaltyFactorization.pointer = calloc(*p * *p, sizeof(double));
  muPenaltyFactorization.rows = *p;
  muPenaltyFactorization.columns = *p;

  matrixSquareRoot(&muPenalty, &muPenaltyFactorization);

  //Initialisation
  *sigmaSquared = 1.0;
  identityMatrix(&DAsMatrix);
  
  for(i = 0; i < *p; i++)
  {
    for(j = 0; j < *p; j++)
    {
      //This expression is valid because D is initialised to the identity matrix
      //It will have to change if the initialisation is improved
      Dv.pointer[i + j * *p] = *lambdaV * G[i + j * *p] + (i == j);
    }
  }

  invertMatrix(&Dv, &Dv);

  calculateYiPrecision(Xi, &Dv, sigmaSquared, *n, inverseVi);

  penalisedLeastSquaresUsingFactorization(&XAsMatrix, &yAsVector, &muAsVector, &muPenaltyFactorization);

  if(*zeroIntercept)
  {
    muAsVector.pointer[0] = 0.0;
  }
  calculateLikelihood(yi, Xi, inverseVi, &muAsVector, *n, likelihood);
  
  for(*iterations = 0; *iterations < *maxIterations && (*likelihood - oldLikelihood) > *deltaEM; (*iterations)++)
  {
    EStep(yi, Xi, inverseVi, *N, *n, &Dv, &muAsVector, vi, epsiloni, *zeroIntercept);
    MStep(yi,
          &XAsMatrix,
          Xi,
          vi,
          epsiloni,
          inverseVi,
          &muAsVector,
          &muPenaltyFactorization,
          &Dv,
          &DAsMatrix,
          &GAsMatrix,
          *lambdaV,
          sigmaSquared,
          *N,
          *n,
          *zeroIntercept);

    oldLikelihood = *likelihood;
    calculateLikelihood(yi, Xi, inverseVi, &muAsVector, *n, likelihood);
  }

  calculateDegreesOfFreedom(&XAsMatrix, Xi, inverseVi, &GAsMatrix, &Dv, *lambdaMu, *n, dfMu, dfV);
  EStep(yi, Xi, inverseVi, *N, *n, &Dv, &muAsVector, vi, epsiloni, *zeroIntercept);

  if((*likelihood - oldLikelihood) < 0)
  {
    *info = INFO_LIKELIHOOD_DECREASED;
  }
  if(*verbose) Rprintf("EM converged in %03d iterations\n", *iterations + 1);

  for(i = 0; i < *n; i++)
  {
    freeVector(yi[i]);
    freeVector(timei[i]);
    freeMatrix(Xi[i]);
    freeMatrix(inverseVi[i]);
    freeVector(epsiloni[i]);

    //Don't call freeVector here because the memory doesn't belong to us
    free(vi[i]);
  }

  //freeMatrix(Z);
  free(vi);
  free(epsiloni);
  free(inverseV.pointer);
  free(inverseVi);
  free(Dv.pointer);
  free(muPenalty.pointer);
  free(muPenaltyFactorization.pointer);
  free(Xi);
  free(timei);
  free(yi);
}

void EStep(Vector** yi, Matrix** Xi, Matrix** inverseVi, int N, int n, Matrix* Dv, Vector* mu, Vector** vi, Vector** epsiloni, int zeroIntercept)
{
  Vector yiCentered;
  Vector inverseViYiCentered;
  int i;

  double* yiCenteredBuffer = calloc(N, sizeof(double));
  double* currentYiCenteredPointer = yiCenteredBuffer;

  double* inverseViYiCenteredBuffer = calloc(N, sizeof(double));
  double* currentInverseViYiCenteredPointer = inverseViYiCenteredBuffer;

  for(i = 0; i < n; i++)
  {
    //yiCentered.pointer = calloc(yi[i]->length, sizeof(double));
    yiCentered.length = yi[i]->length;
    yiCentered.pointer = currentYiCenteredPointer;
    
    //inverseViYiCentered.pointer = calloc(yi[i]->length, sizeof(double));
    inverseViYiCentered.length = yi[i]->length;
    inverseViYiCentered.pointer = currentInverseViYiCenteredPointer;

    //yiCentered <- (yi - Xi %*% eta)
    matrixVectorMultiply(Xi[i], mu, yi[i], &yiCentered, -1.0, 1.0, 0);
    //inverseViYiCentered <- solve(Vi) %*% (yi - Xi %*% eta)
    matrixVectorMultiply(inverseVi[i], &yiCentered, 0, &inverseViYiCentered, 1.0, 0.0, 0);
    //vi <- t(Xi) %*% solve(Vi) %*% (yi - Xi %*% eta)
    matrixVectorMultiply(Xi[i], &inverseViYiCentered, 0, vi[i], 1.0, 0.0, 1);
    //vi <- Dv %*% t(Xi) %*% solve(Vi) %*% (yi - Xi %*% eta)
    matrixVectorMultiply(Dv, vi[i], 0, vi[i], 1.0, 0.0, 0);

    matrixVectorMultiply(Xi[i], vi[i], &yiCentered, epsiloni[i], -1.0, 1.0, 0);
    if(zeroIntercept)
    {
      vi[i]->pointer[0] = 0.0;
    }

    //free(inverseViYiCentered.pointer);
    //free(yiCentered.pointer);
    currentYiCenteredPointer += yi[i]->length;
    currentInverseViYiCenteredPointer += yi[i]->length;
  }

  free(inverseViYiCenteredBuffer);
  free(yiCenteredBuffer);
}

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
           int n,
           int zeroIntercept)
{
  Matrix U, A, conditionalCovariance, Dinverse; // A <- Dv %*% t(X) %*% t(chol(solve(Vi[[1]])))
  Vector yCentered, yiCentered; // y - Z %*% v
  int i, j;
  
  A.rows = D->rows;

  conditionalCovariance.rows = Dv->rows;
  conditionalCovariance.columns = Dv->rows;
  conditionalCovariance.pointer = calloc(conditionalCovariance.rows * conditionalCovariance.columns, sizeof(double));
  
  yCentered.length = N;
  yCentered.pointer = calloc(N, sizeof(double));

  Dinverse.rows = Dv->rows;
  Dinverse.columns = Dv->columns;
  Dinverse.pointer = calloc(Dinverse.rows * Dinverse.columns, sizeof(double));

  double* currentYCenteredPointer = yCentered.pointer;

  double sumOfSquares = 0.0, trace = 0.0;
  zeroMatrix(D);
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < epsiloni[i]->length; j++)
    {
      sumOfSquares += epsiloni[i]->pointer[j] * epsiloni[i]->pointer[j];
      trace += inverseVi[i]->pointer[j + j * inverseVi[i]->rows];
    }

    U.pointer = calloc(inverseVi[i]->rows * inverseVi[i]->columns, sizeof(double));
    U.rows = inverseVi[i]->rows;
    U.columns = inverseVi[i]->columns;
    
    A.pointer = calloc(A.rows * inverseVi[i]->rows, sizeof(double));
    A.columns = inverseVi[i]->rows;

    choleskyFactorization(inverseVi[i], &U);
    // A <- Dv %*% t(X)
    //multiplyMatrices(Dv, Xi[i], &A, 0, 1);
    // A <- Dv %*% t(X) %*% t(chol(solve(Vi[[1]])))
    //multiplyMatrices(&A, &U, &A, 0, 1);
    multiplyMatrices(Xi[i], &U, &A, 1, 1);

    symmetricRankKUpdate(&A, &conditionalCovariance, 1.0, 1.0);
    symmetricRank1Update(D, vi[i], 1.0);

    yiCentered.length = yi[i]->length;
    yiCentered.pointer = currentYCenteredPointer;

    matrixVectorMultiply(Xi[i], vi[i], yi[i], &yiCentered, -1.0, 1.0, 0);

    currentYCenteredPointer += yi[i]->length;

    free(A.pointer);
    free(U.pointer);
  }
  U.pointer = calloc(Dv->columns * Dv->rows, sizeof(double));
  U.rows = Dv->columns;
  U.columns = Dv->rows;
  choleskyFactorization(&conditionalCovariance, &U);
  A.pointer = calloc(Dv->columns * U.columns, sizeof(double));
  A.rows = Dv->rows;
  A.columns = Dv->columns;
  multiplyMatrices(Dv, &U, &A, 0, 1);
  symmetricRankKUpdate(&A, &conditionalCovariance, 1.0, 0.0);
  free(U.pointer);
  free(A.pointer);

  for(i = 0; i < conditionalCovariance.rows * conditionalCovariance.columns; i++)
  {
    D->pointer[i] += n * Dv->pointer[i] - conditionalCovariance.pointer[i];
    D->pointer[i] /= (double) n;
  }

  invertMatrix(D, &Dinverse);

  for(i = 0; i < Dv->rows; i++)
  {
    for(j = 0; j < Dv->columns; j++)
    {
      Dv->pointer[i + j * Dv->rows] = lambdaV * G->pointer[i + j * G->rows] + Dinverse.pointer[i + j * Dinverse.rows];
    }
  }

  invertMatrix(Dv, Dv);

  *sigmaSquared = (1.0) / N * (sumOfSquares + *sigmaSquared * N - *sigmaSquared * *sigmaSquared * trace);

  penalisedLeastSquaresUsingFactorization(X, &yCentered, mu, muPenaltyFactorization);


  if(zeroIntercept)
  {
    mu->pointer[0] = 0.0;
  }
  calculateYiPrecision(Xi, Dv, sigmaSquared, n, inverseVi);

  free(Dinverse.pointer);
  free(yCentered.pointer);
  free(conditionalCovariance.pointer);
}

void calculateYiPrecision(Matrix** Xi, Matrix* Dv, double* sigmaSquared, int n, Matrix** inverseVi)
{
  int i;
  Matrix U, XiUTransposed;

  XiUTransposed.columns = Dv->rows;

  U.pointer = calloc(Dv->rows * Dv->rows, sizeof(double));
  U.rows = Dv->rows;
  U.columns = Dv->columns;

  choleskyFactorization(Dv, &U);

  for(i = 0; i < n; i++)
  {
    XiUTransposed.pointer = calloc(Xi[i]->rows * U.rows, sizeof(double));
    XiUTransposed.rows = Xi[i]->rows;

    identityMatrix(inverseVi[i]);

    multiplyMatrices(Xi[i], &U, &XiUTransposed, 0, 1);

    symmetricRankKUpdate(&XiUTransposed, inverseVi[i], 1.0, *sigmaSquared);
    invertMatrix(inverseVi[i], inverseVi[i]);

    free(XiUTransposed.pointer);
  }

  free(U.pointer);
}

void calculateLikelihood(Vector** yi, Matrix** Xi, Matrix** inverseVi, Vector* mu, int n, double *likelihood)
{
  Matrix U;
  Vector Uyi;
  Matrix UXi;
  int i, j;

  //matrixVectorMultiply(Xi[i], vi[i], 0, vi[i], 1.0, 0.0, 1);

  double sumOfSquares = 0.0;
  double determinant;

  *likelihood = 0.0;

  for(i = 0; i < n; i++)
  {
    UXi.rows = Xi[i]->rows;
    UXi.columns = Xi[i]->columns;
    UXi.pointer = calloc(UXi.rows * UXi.columns, sizeof(double));

    U.rows = inverseVi[i]->rows;
    U.columns = inverseVi[i]->columns;
    U.pointer = calloc(U.rows * U.columns, sizeof(double));

    Uyi.pointer = calloc(yi[i]->length, sizeof(double));
    Uyi.length = yi[i]->length;

    choleskyFactorization(inverseVi[i], &U);
    matrixVectorMultiply(&U, yi[i], 0, &Uyi, 1.0, 0.0, 0);
    multiplyMatrices(&U, Xi[i], &UXi, 0, 0);

    matrixVectorMultiply(&UXi, mu, &Uyi, &Uyi, -1.0, 1.0, 0);
    for(j = 0; j < Uyi.length; j++)
    {
      sumOfSquares += Uyi.pointer[j] * Uyi.pointer[j];
    }

    logDeterminant(inverseVi[i], &determinant);

    *likelihood += determinant;

    free(UXi.pointer);
    free(Uyi.pointer);
    free(U.pointer);
  }

  *likelihood -= sumOfSquares;
  *likelihood *= 0.5;
}

void calculateDegreesOfFreedom(Matrix* X, Matrix** Xi, Matrix** inverseVi, Matrix* G, Matrix* Dv, double lambdaMu, int n, double* dfMu, double* dfV)
{
  int i, j;
  
  Matrix** tXiinverseViXis = calloc(n, sizeof(Matrix*));
  Matrix** tXiinverseViXiDvs = calloc(n, sizeof(Matrix*));
  Matrix U, A;
  Matrix Hessian, inverseHessian;
  
  Hessian.rows = Dv->rows;
  Hessian.columns = Dv->columns;
  Hessian.pointer = calloc(Dv->rows * Dv->columns, sizeof(double));
  
  inverseHessian.rows = Dv->rows;
  inverseHessian.columns = Dv->columns;
  inverseHessian.pointer = calloc(Dv->rows * Dv->columns, sizeof(double));
  
  *dfMu = 0.0;
  *dfV = 0.0;
  
  for(i = 0; i < n; i++)
  {
    tXiinverseViXis[i] = calloc(1, sizeof(Matrix));
    tXiinverseViXiDvs[i] = calloc(1, sizeof(Matrix));
  
    tXiinverseViXis[i]->rows = Dv->rows;
    tXiinverseViXis[i]->columns = Dv->rows;
    tXiinverseViXis[i]->pointer = calloc(Dv->rows * Dv->columns, sizeof(double));
    
    tXiinverseViXiDvs[i]->rows = Dv->rows;
    tXiinverseViXiDvs[i]->columns = Dv->rows;
    tXiinverseViXiDvs[i]->pointer = calloc(Dv->rows * Dv->columns, sizeof(double));

    U.rows = inverseVi[i]->rows;
    U.columns = inverseVi[i]->columns;
    U.pointer = calloc(U.rows * U.columns, sizeof(double));
    
    A.rows = Xi[i]->columns;
    A.columns = U.rows;
    A.pointer = calloc(A.rows * A.columns, sizeof(double));
    
    choleskyFactorization(inverseVi[i], &U);
    multiplyMatrices(Xi[i], &U, &A, 1, 1);
    symmetricRankKUpdate(&A, tXiinverseViXis[i], 1.0, 1.0);
    for(j = 0; j < (Hessian.rows * Hessian.columns); j++)
    {
      Hessian.pointer[j] += tXiinverseViXis[i]->pointer[j];
    }

    //symmetricRankKUpdate(&A, &Hessian, 1.0, 1.0);

    multiplyMatrices(tXiinverseViXis[i], Dv, tXiinverseViXiDvs[i], 0, 0);
    for(j = 0; j < tXiinverseViXiDvs[i]->rows; j++)
    {
      *dfV += tXiinverseViXiDvs[i]->pointer[j + j * tXiinverseViXiDvs[i]->rows];
    }
    
    free(A.pointer);
    free(U.pointer);
  }
  
  invertMatrix(&Hessian, &inverseHessian);
  
  A.rows = Dv->rows;
  A.columns = Dv->columns;
  A.pointer = calloc(A.rows * A.columns, sizeof(double));

  U.rows = Dv->rows;
  U.columns = Dv->columns;
  U.pointer = calloc(U.rows * U.columns, sizeof(double));

  multiplyMatrices(&inverseHessian, G, &A, 0, 0);  
  for(i = 0; i < A.rows; i++)
  {
    for(j = 0; j < A.rows; j++)
    {
      A.pointer[i + j * A.rows] = (i==j) + A.pointer[i + j * A.rows] * lambdaMu;
    }
  }
  invertMatrix(&A, &U);

  for(i = 0; i < U.rows; i++)
  {
    *dfMu += U.pointer[i + i * U.rows];
  }

  for(i = 0; i < Hessian.rows; i++)
  {
    for(j = 0; j < Hessian.columns; j++)
    {
      Hessian.pointer[i + j * Hessian.rows] += lambdaMu * G->pointer[i + j * G->rows];
    }
  }

  invertMatrix(&Hessian, &inverseHessian);
  
  for(i = 0; i < n; i++)
  {
    multiplyMatrices(tXiinverseViXiDvs[i], tXiinverseViXis[i], &A, 0, 0);
    multiplyMatrices(&A, &inverseHessian, &U, 0, 0);
    for(j = 0; j < U.rows; j++)
    {
      *dfV -= U.pointer[j + j * U.rows];
    }

    freeMatrix(tXiinverseViXis[i]);
    freeMatrix(tXiinverseViXiDvs[i]);
  }

  free(tXiinverseViXis);
  free(tXiinverseViXiDvs);
  free(U.pointer);
  free(A.pointer);
  free(inverseHessian.pointer);
  free(Hessian.pointer);
}