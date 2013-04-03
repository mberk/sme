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

#ifndef NELDERMEAD_H
#define NELDERMEAD_H

#include <R.h>

typedef double optimfn(int, double *, void *);

typedef struct
{
  double* pointer;
  double value;
  int length;
} Vertex;

#define big             1.0e+35

void NelderMead(int n, double *Bvec, double *X, double *Fmin, optimfn fminfn,
	   int *fail, double abstol, double intol, void *ex,
	   double alpha, double bet, double gamm, int trace,
	   int *fncount, int maxit);

#endif