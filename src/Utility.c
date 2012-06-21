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
#include "Utility.h"
#include "Vector.h"
#include "Matrix.h"

// Returns an array of the unique elements in the array argument
// The number of unique elements is returned in the numberOfUniqueElements pointer
// Assume that arraySize >= 1 !!
double* getUniqueDoubleArray(double* array, int arraySize, int* numberOfUniqueElements)
{
  int i, uniqueIndex;
  double* uniqueElements;

  // First sort the array
  R_rsort(array, arraySize);

  // Count the number of unique elements
  *numberOfUniqueElements = 1;

  for(i = 1; i < arraySize; i++)
  {
    if(array[i] != array[i-1])
    {
      (*numberOfUniqueElements)++;
    }
  }

  uniqueElements = calloc(*numberOfUniqueElements, sizeof(double));
  uniqueElements[0] = array[0];
  uniqueIndex = 1;

  for(i = 1; i < arraySize; i++)
  {
    if(array[i] != array[i-1])
    {
      uniqueElements[uniqueIndex++] = array[i];
    }
  }

  return uniqueElements;
}

int* getUniqueIntegerArray(int* array, int arraySize, int* numberOfUniqueElements)
{
  int i, uniqueIndex;
  int* uniqueElements;

  // First sort the array
  R_isort(array, arraySize);

  // Count the number of unique elements
  *numberOfUniqueElements = 1;

  for(i = 1; i < arraySize; i++)
  {
    if(array[i] != array[i-1])
    {
      (*numberOfUniqueElements)++;
    }
  }

  uniqueElements = calloc(*numberOfUniqueElements, sizeof(int));
  uniqueElements[0] = array[0];
  uniqueIndex = 1;

  for(i = 1; i < arraySize; i++)
  {
    if(array[i] != array[i-1])
    {
      uniqueElements[uniqueIndex++] = array[i];
    }
  }

  return uniqueElements;
}

// Calling function is responsible for freeing the return value
Vector** splitVector(Vector* vector, int* vectorIdentifiers, int* numberOfVectors)
{
  int* uniqueVectorIdentifiers;
  int* indices;
  Vector** vectors;
  int i, j;

  // First work out how many vectors there are
  uniqueVectorIdentifiers =
    getUniqueIntegerArray(vectorIdentifiers, vector->length, numberOfVectors);

  vectors = calloc(*numberOfVectors, sizeof(Vector*));
  indices = calloc(*numberOfVectors, sizeof(int));
  
  for(i = 0; i < *numberOfVectors; i++)
  {
    vectors[i] = calloc(1, sizeof(Vector));
  }

  for(i = 0; i < vector->length; i++)
  {
    for(j = 0; j < *numberOfVectors; j++)
    {
      if(vectorIdentifiers[i] == uniqueVectorIdentifiers[j])
      {
        vectors[j]->length++;
        break;
      }
    }
  }

  for(i = 0; i < *numberOfVectors; i++)
  {
    vectors[i]->pointer = calloc(vectors[i]->length, sizeof(double));
  }

  for(i = 0; i < vector->length; i++)
  {
    for(j = 0; j < *numberOfVectors; j++)
    {
      if(vectorIdentifiers[i] == uniqueVectorIdentifiers[j])
      {
        vectors[j]->pointer[indices[j]++] = vector->pointer[i];
        break;
      }
    }
  }

  free(indices);
  free(uniqueVectorIdentifiers);

  return vectors;
}

//Splits a matrix into a number of sub matrices, where each row is associated with a different
//matrix
Matrix** splitMatrix(Matrix* matrix, int* matrixIdentifiers, int* numberOfMatrices)
{
  int* uniqueMatrixIdentifiers;
  int* indices;
  Matrix** matrices;
  int i, j, k;

  //First work out how many matrices there are
  uniqueMatrixIdentifiers =
    getUniqueIntegerArray(matrixIdentifiers, matrix->rows, numberOfMatrices);

  matrices = calloc(*numberOfMatrices, sizeof(Matrix*));
  indices = calloc(*numberOfMatrices, sizeof(int));
  
  for(i = 0; i < *numberOfMatrices; i++)
  {
    matrices[i] = calloc(1, sizeof(Matrix));
    matrices[i]->columns = matrix->columns;
  }

  for(i = 0; i < matrix->rows; i++)
  {
    for(j = 0; j < *numberOfMatrices; j++)
    {
      if(matrixIdentifiers[i] == uniqueMatrixIdentifiers[j])
      {
        (matrices[j]->rows)++;
        break;
      }
    }
  }

  for(i = 0; i < *numberOfMatrices; i++)
  {
    matrices[i]->pointer = calloc(matrices[i]->rows * matrices[i]->columns, sizeof(double));
  }
  
  for(i = 0; i < matrix->rows; i++)
  {
    for(j = 0; j < *numberOfMatrices; j++)
    {
      if(matrixIdentifiers[i] == uniqueMatrixIdentifiers[j])
      {
        //Copy all k elements of the row
        for(k = 0; k < matrix->columns; k++)
        {
          matrices[j]->pointer[indices[j] + k * matrices[j]->rows] = matrix->pointer[i + k * matrix->rows];
        }
        indices[j]++;
      }
    }
  }

  free(indices);
  free(uniqueMatrixIdentifiers);

  return(matrices);
}

void split(double* matrixBuffer, int* rows, int* columns, int* matrixIdentifiers)
{
  Matrix matrix;
  int numberOfMatrices;
  int i;

  matrix.pointer = matrixBuffer;
  matrix.rows = *rows;
  matrix.columns = *columns;

  Matrix** matrices = splitMatrix(&matrix, matrixIdentifiers, &numberOfMatrices);

  Rprintf("There are %d submatrices\n", numberOfMatrices);

  for(i = 0; i < numberOfMatrices; i++)
  {
    printMatrix(matrices[i]);
    freeMatrix(matrices[i]);
  }

  free(matrices);
}