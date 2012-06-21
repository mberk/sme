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

#include "BlockDiagonalMatrix.h"
#include "Matrix.h"


// Calling function is responsible for freeing result!
Matrix* constructBlockDiagonalMatrix(Matrix** matrices, int numberOfMatrices)
{
  int i, j, k;
  double* currentMatrixStart;

  Matrix* result = calloc(1, sizeof(Matrix));

  // First determine dimensions of the result
  for(i = 0; i < numberOfMatrices; i++)
  {
    result->rows += matrices[i]->rows;
    result->columns += matrices[i]->columns;
  }

  result->pointer = calloc(result->rows * result->columns, sizeof(double));

  currentMatrixStart = result->pointer;

  for(i = 0; i < numberOfMatrices; i++)
  {
    for(j = 0; j < matrices[i]->rows; j++)
    {
      for(k = 0; k < matrices[i]->columns; k++)
      {
        currentMatrixStart[j + k * result->rows] = matrices[i]->pointer[j + k * matrices[i]->rows];
      }
    }

    currentMatrixStart += matrices[i]->columns * result->rows + matrices[i]->rows;
  }

  return result;
}