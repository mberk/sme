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

#include "Vector.h"

void freeVector(Vector* vector)
{
  free(vector->pointer);
  free(vector);
}

void printVector(Vector* vector)
{
  int i;

  Rprintf("A %d length vector:\n", vector->length);

  for(i = 0; i < vector->length; i++)
  {
    Rprintf("%f\n" , vector->pointer[i]);
  }
}