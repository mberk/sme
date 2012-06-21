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
#include <R_ext/Utils.h>

#include "IncidenceMatrix.h"
#include "Matrix.h"

void incidenceMatrix(double* timePoints, int* numberOfTimePoints, double* uniqueTimePoints, int* numberOfUniqueTimePoints, double* result)
{
  Matrix* X = createIncidenceMatrix(timePoints, *numberOfTimePoints, uniqueTimePoints, *numberOfUniqueTimePoints);
  Rprintf("Copying %d doubles\n", *numberOfTimePoints * *numberOfUniqueTimePoints);
  memcpy(result, X->pointer, sizeof(double) * *numberOfTimePoints * *numberOfUniqueTimePoints);
  freeMatrix(X);
}

Matrix* createIncidenceMatrix(double* timePoints, int numberOfTimePoints, double* uniqueTimePoints, int numberOfUniqueTimePoints)
{
  int i, j;

  // Just make sure that the unique time points are sorted first
  //R_rsort(timePoints, numberOfTimePoints);
  R_rsort(uniqueTimePoints, numberOfUniqueTimePoints);

  Matrix* incidenceMatrix = allocateMatrix(numberOfTimePoints, numberOfUniqueTimePoints);

  for(i = 0; i < numberOfTimePoints; i++)
  {
    for(j = 0; j < numberOfUniqueTimePoints; j++)
    {
      if(timePoints[i] == uniqueTimePoints[j])
      {
        incidenceMatrix->pointer[j * incidenceMatrix->rows + i] = 1.0;
        break;
      }
    }
  }

  return incidenceMatrix;
}