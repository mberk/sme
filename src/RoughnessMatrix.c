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

#include "Matrix.h"
#include "RoughnessMatrix.h"
#include "LinearAlgebra.h"

Matrix* createRoughnessMatrix(double* timePoints, int numberOfTimePoints)
{
  int r;

  // h is the vector of differences between successive time points
  double* h = calloc(numberOfTimePoints - 1, sizeof(double));

  Matrix* roughnessMatrix = allocateMatrix(numberOfTimePoints, numberOfTimePoints);
  Matrix* A = allocateMatrix(numberOfTimePoints, numberOfTimePoints - 2);
  Matrix* B = allocateMatrix(numberOfTimePoints - 2, numberOfTimePoints - 2);
  Matrix* BInverse = allocateMatrix(numberOfTimePoints - 2, numberOfTimePoints - 2);
  Matrix* ATimesBInverse = allocateMatrix(numberOfTimePoints, numberOfTimePoints - 2);

  for(r = 0; r < numberOfTimePoints - 1; r++)
  {
    h[r] = timePoints[r+1] - timePoints[r];
  }

  for(r = 0; r < numberOfTimePoints - 2; r++)
  {
    A->pointer[r * A->rows + r] = 1.0 / h[r];
    A->pointer[r * A->rows + r + 1] = -((1.0 / h[r]) + (1.0 / h[r+1]));
    A->pointer[r * A->rows + r + 2] = 1.0 / h[r+1];
  }

  B->pointer[0] = (h[0] + h[1]) / 3.0;
  B->pointer[1] = h[1] / 6.0;

  if(numberOfTimePoints > 4)
  {
    for(r = 0; r < numberOfTimePoints - 4; r++)
    {
      B->pointer[(r + 1) * B->rows + r] = h[r+1] / 6.0;
      B->pointer[(r + 1) * B->rows + r + 1] = (h[r+1] + h[r+2]) / 3.0;
      B->pointer[(r + 1) * B->rows + r + 2] = h[r+2] / 6.0;
    }
  }
  
  B->pointer[(numberOfTimePoints - 3) * B->rows + (numberOfTimePoints - 4)] =
    h[numberOfTimePoints - 3] / 6.0;

  B->pointer[(numberOfTimePoints - 3) * B->rows + (numberOfTimePoints - 3)] =
    (h[numberOfTimePoints - 3] + h[numberOfTimePoints - 2]) / 3.0;

  invertMatrix(B, BInverse);
  multiplyMatrices(A, BInverse, ATimesBInverse, 0, 0);
  multiplyMatrices(ATimesBInverse, A, roughnessMatrix, 0, 1);

  free(h);

  freeMatrix(BInverse);
  freeMatrix(B);
  freeMatrix(A);
  freeMatrix(ATimesBInverse);

  return roughnessMatrix;
}