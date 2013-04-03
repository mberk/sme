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

#include "Matrix.h"

Matrix* allocateMatrix(int rows, int columns)
{
  Matrix* matrix = calloc(1, sizeof(Matrix));
  matrix->pointer = calloc(rows * columns, sizeof(double));
  matrix->rows = rows;
  matrix->columns = columns;

  return matrix;
}

void identityMatrix(Matrix* matrix)
{
  int i, j;

  for(i = 0; i < matrix->rows; i++)
  {
    for(j = 0; j < matrix->columns; j++)
    {
      matrix->pointer[i + j * matrix->rows] = (double)(i == j);
    }
  }
}

void zeroMatrix(Matrix* matrix)
{
  int i;
  
  for(i = 0; i < matrix->rows * matrix->columns; i++)
  {
    matrix->pointer[i] = 0.0;
  }
}

void printMatrix(Matrix* matrix)
{
  int i, j;

  Rprintf("A %d x %d matrix:\n", matrix->rows, matrix->columns);

  for(i = 0; i < matrix->rows; i++)
  {
    for(j = 0; j < matrix->columns; j++)
    {
      Rprintf("%f ", matrix->pointer[i + j * matrix->rows]);
    }
    Rprintf("\n");
  }
}

void freeMatrix(Matrix* matrix)
{
  free(matrix->pointer);
  free(matrix);
}