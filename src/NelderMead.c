#include <string.h>

#include "NelderMead.h"

void NelderMead(int n, double *Bvec, double *X, double *Fmin, optimfn fminfn,
	   int *fail, double abstol, double intol, void *ex,
	   double alpha, double bet, double gamm, int trace,
	   int *fncount, int maxit)
{
    char action[50];
    int C;

    // Changed from Rboolean to int
    int calcvert;

    double convtol, f;
    int funcount=0, H, i, j, L=0;
    int n1=0;
    double oldsize;
    double size, step, temp, trystep;
    char tstr[12];
    double VH, VL, VR;

    if (maxit <= 0)
    {
      *Fmin = fminfn(n, Bvec, ex);
      *fncount = 0;
      *fail = 0;
      return;
    }

    if (trace) Rprintf("  Nelder-Mead direct search function minimizer\n");

    Vertex* vertices = calloc(n + 2, sizeof(Vertex));
    for(i = 0; i < (n+2); i++)
    {
      vertices[i].length = n;
      vertices[i].pointer = calloc(n, sizeof(double));
    }

    *fail = FALSE;
    f = fminfn(n, Bvec, ex);
    if (!R_FINITE(f))
    {
      //error(_("function cannot be evaluated at initial parameters"));
      *fail = TRUE;
    }
    else
    {
      if (trace) Rprintf("function value for initial parameters = %f\n", f);
      funcount = 1;
      convtol = intol * (fabs(f) + intol);
      if (trace) Rprintf("  Scaled convergence tolerance is %g\n", convtol);
      n1 = n + 1;
      C = n + 2;
      //P[n1 - 1][0] = f;
      vertices[0].value = f;
      for (i = 0; i < n; i++)
        //P[i] = Bvec[i];
        vertices[0].pointer[i] = Bvec[i];

      L = 1;
      size = 0.0;

      step = 0.0;
      for (i = 0; i < n; i++)
      {
        if (0.1 * fabs(Bvec[i]) > step) step = 0.1 * fabs(Bvec[i]);
      }
      if (step == 0.0) step = 0.1;
      if (trace) Rprintf("Stepsize computed as %f\n", step);
      for (j = 2; j <= n1; j++)
      {
        strcpy(action, "BUILD          ");
        for (i = 0; i < n; i++)
          //P[i][j - 1] = Bvec[i];
          //i-th coordinate of the j-1th vertex
          vertices[j-1].pointer[i] = Bvec[i];

        trystep = step;
        //j-2th coordinate of the j-1th vertex
        //while (P[j - 2][j - 1] == Bvec[j - 2])
        while(vertices[j - 1].pointer[j - 2] == Bvec[j - 2])
        {
          //P[j - 2][j - 1] = Bvec[j - 2] + trystep;
          vertices[j - 1].pointer[j - 2] = Bvec[j - 2] + trystep;
          trystep *= 10;
        }
        size += trystep;
      }
      oldsize = size;
      calcvert = 1;
      do
      {
        if (calcvert)
        {
          for (j = 0; j < n1; j++)
          {
            if (j + 1 != L)
            {
              for (i = 0; i < n; i++)
                //Bvec[i] = P[i][j];
                Bvec[i] = vertices[j].pointer[i];//P[i][j];
              f = fminfn(n, Bvec, ex);
              if (!R_FINITE(f)) f = big;
              funcount++;
              vertices[j].value = f;
              //P[n1 - 1][j] = f;
            }
          }
          calcvert = 0;
        }

        //VL = P[n1 - 1][L - 1];
        VL = vertices[L - 1].value;
        VH = VL;
        H = L;

        for (j = 1; j <= n1; j++)
        {
          if (j != L)
          {
            //f = P[n1 - 1][j - 1];
            f = vertices[j - 1].value;
            if (f < VL)
            {
              L = j;
              VL = f;
            }
            if (f > VH)
            {
              H = j;
              VH = f;
            }
          }
        }

        if (VH <= VL + convtol || VL <= abstol) break;

        if (trace)
        {
          snprintf(tstr, 12, "%5d", funcount);
          Rprintf("%s%s %f %f\n", action, tstr, VH, VL);
        }

        for (i = 0; i < n; i++)
        {
          //temp = -P[i][H - 1];
          temp = -vertices[H - 1].pointer[i];
          for (j = 0; j < n1; j++)
            //temp += P[i][j];
            temp += vertices[j].pointer[i];
            //P[i][C - 1] = temp / n;
          vertices[C - 1].pointer[i] = temp / n;
        }
        for (i = 0; i < n; i++)
          //Bvec[i] = (1.0 + alpha) * P[i][C - 1] - alpha * P[i][H - 1];
          Bvec[i] = (1.0 + alpha) * vertices[C - 1].pointer[i] - alpha * vertices[H - 1].pointer[i];
        f = fminfn(n, Bvec, ex);
        if (!R_FINITE(f)) f = big;
        funcount++;
        strcpy(action, "REFLECTION     ");
        VR = f;
        if (VR < VL)
        {
          //P[n1 - 1][C - 1] = f;
          vertices[C - 1].value = f;
          for (i = 0; i < n; i++)
          {
            //f = gamm * Bvec[i] + (1 - gamm) * P[i][C - 1];
            f = gamm * Bvec[i] + (1 - gamm) * vertices[C - 1].pointer[i];
            //P[i][C - 1] = Bvec[i];
            vertices[C - 1].pointer[i] = Bvec[i];
            Bvec[i] = f;
          }
          f = fminfn(n, Bvec, ex);
          if (!R_FINITE(f)) f = big;
          funcount++;
          if (f < VR)
          {
            for (i = 0; i < n; i++)
              //P[i][H - 1] = Bvec[i];
              vertices[H - 1].pointer[i] = Bvec[i];
            //P[n1 - 1][H - 1] = f;
            vertices[H - 1].value = f;
            strcpy(action, "EXTENSION      ");
          }
          else
          {
            for (i = 0; i < n; i++)
              //P[i][H - 1] = P[i][C - 1];
              vertices[H - 1].pointer[i] = vertices[C - 1].pointer[i];
            //P[n1 - 1][H - 1] = VR;
            vertices[H - 1].value = VR;
          }
        }
        else
        {
          strcpy(action, "HI-REDUCTION   ");
          if (VR < VH)
          {
            for (i = 0; i < n; i++)
              //P[i][H - 1] = Bvec[i];
              vertices[H - 1].pointer[i] = Bvec[i];
            //P[n1 - 1][H - 1] = VR;
            vertices[H - 1].value = VR;
            strcpy(action, "LO-REDUCTION   ");
          }

          for (i = 0; i < n; i++)
            //Bvec[i] = (1 - bet) * P[i][H - 1] + bet * P[i][C - 1];
            Bvec[i] = (1 - bet) * vertices[H - 1].pointer[i] + bet * vertices[C - 1].pointer[i];
          f = fminfn(n, Bvec, ex);
          if (!R_FINITE(f)) f = big;
          funcount++;

          //if (f < P[n1 - 1][H - 1])
          if (f < vertices[H - 1].value)
          {
            for (i = 0; i < n; i++)
              //P[i][H - 1] = Bvec[i];
              vertices[H - 1].pointer[i] = Bvec[i];
            vertices[H - 1].value = f;
          }
          else
          {
            if (VR >= VH)
            {
              strcpy(action, "SHRINK         ");
              calcvert = 1;
              size = 0.0;
              for (j = 0; j < n1; j++)
              {
                if (j + 1 != L)
                {
                  for (i = 0; i < n; i++)
                  {
                    //P[i][j] = bet * (P[i][j] - P[i][L - 1])
                    //  + P[i][L - 1];
                    vertices[j].pointer[i] = bet * (vertices[j].pointer[i] - vertices[L - 1].pointer[i])
                      + vertices[L - 1].pointer[i];
                    //size += fabs(P[i][j] - P[i][L - 1]);
                    size += fabs(vertices[j].pointer[i] - vertices[L - 1].pointer[i]);
                  }
                }
              }
              if (size < oldsize)
              {
                oldsize = size;
              }
              else
              {
                if (trace) Rprintf("Polytope size measure not decreased in shrink\n");
                *fail = 10;
                break;
              }
            }
          }
        }
      } while (funcount <= maxit);
    }

    if (trace)
    {
      Rprintf("Exiting from Nelder Mead minimizer\n");
      Rprintf("    %d function evaluations used\n", funcount);
    }
    //*Fmin = P[n1 - 1][L - 1];
    *Fmin = vertices[L - 1].value;
    for (i = 0; i < n; i++)
      //X[i] = P[i][L - 1];
      X[i] = vertices[L - 1].pointer[i];
    if (funcount > maxit) *fail = 1;
    *fncount = funcount;
    
    for(i = 0; i < (n+2); i++)
    {
      free(vertices[i].pointer);
    }
    
    free(vertices);
}
