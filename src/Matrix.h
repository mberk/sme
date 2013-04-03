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

#ifndef MATRIX_H
#define MATRIX_H

#define setElement(MatrixPointer,Row,Column,Value) MatrixPointer->pointer[MatrixPointer->rows * Column + Row] = Value

typedef struct
{
  double* pointer;
  int rows;
  int columns;
} Matrix;

Matrix* allocateMatrix(int rows, int columns);

void zeroMatrix(Matrix* matrix);
void identityMatrix(Matrix* matrix);

void printMatrix(Matrix* matrix);
void freeMatrix(Matrix* matrix);

#endif