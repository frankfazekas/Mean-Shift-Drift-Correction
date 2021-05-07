#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define OK 0
#define ERR_OVERFLOW 1

double sqrt();

void Fcrosspairs(n1, x1, y1, t1, n2, x2, y2, t2, rmax, taumin, taumax,
        noutmax, nout, dxout, dyout, dtout, status)
     /* inputs */
     int n1, n2, noutmax;
     double *x1, *y1, *x2, *y2, rmax, *t1, *t2, taumin, taumax;
     /* outputs */
     int *nout, *status;
     double *dxout, *dyout, *dtout;
     
{
  int k, kmax, i, j, jleft;
  double x1i, y1i, r2max, xleft, dx, dy, dx2, d2, dt;
 
  r2max = rmax * rmax;
 
  *status = OK;
  *nout = 0;
  k = 0;   /* k is the next available storage location 
              and also the current length of the list */
  kmax = noutmax;
  
  if(n1 == 0 || n2 == 0) 
    return;

  i = 0;
  
  while (i < n1) {

      /* Possibly check for matlab interrupt. */

      for (; i < n1; i++) {

          x1i = x1[i];
          y1i = y1[i];

          /* adjust starting position jleft */
          xleft = x1i - rmax;
          while ((x2[jleft] < xleft) && (jleft + 1 < n2))
              ++jleft;

          /* process from j=jleft until dx > rmax */
          for (j = jleft; j < n2; j++) {
              dx = x2[j] - x1i;
              dx2 = dx * dx;
              if (dx2 > r2max)
                  break;
              dy = y2[j] - y1i;
              d2 = dx2 + dy * dy;
              dt = t2[j] - t1[i];
              if (d2 <= r2max && taumin <= dt && dt <= taumax) {
                  /* add this (i, j) pair to output */
                  if (k >= kmax) {
                      *nout = k;
                      *status = ERR_OVERFLOW;
                      return;
                  }
                  dxout[k] = dx;
                  dyout[k] = dy;
                  dtout[k] = dt;
                  ++k;
              }
          }
      }
  }
  *nout = k;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int n1, n2;
    double *x1, *y1, *t1, *x2, *y2, *t2;
    double rmax, taumin, taumax;
    double *dxout, *dyout, *dtout;
    int *nout, noutmax, status;
    int i;

    /* RHS args are:
     *  0   x1 // coords for first pp
     *  1   y1
     *  2   t1
     *  3   x2 // coords for second pp
     *  4   y2
     *  5   t2
     *  6   rmax
     *  7   taumin
     *  8   taumax
     *  9   noutmax
     */

    /* LHS args are:
     * 0,1,2    dxout, dyout, dtout
     * 3    optional error
     */

    
    n1 = mxGetNumberOfElements(prhs[0]);
    n2 = mxGetNumberOfElements(prhs[3]);

    x1 = (double *) mxGetPr(prhs[0]);
    y1 = (double *) mxGetPr(prhs[1]);
    t1 = (double *) mxGetPr(prhs[2]);
    x2 = (double *) mxGetPr(prhs[3]);
    y2 = (double *) mxGetPr(prhs[4]);
    t2 = (double *) mxGetPr(prhs[5]);

    rmax = mxGetScalar(prhs[6]);
    taumin = mxGetScalar(prhs[7]);
    taumax = mxGetScalar(prhs[8]);

    nout = (int *) mxGetPr(prhs[9]);
    noutmax = *nout;

    for (i=0; i<3; i++)
        plhs[i] = mxCreateNumericMatrix(1, noutmax, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    dxout = (double *) mxGetPr(plhs[0]);
    dyout = (double *) mxGetPr(plhs[1]);
    dtout = (double *) mxGetPr(plhs[2]);
    status = (int) mxGetScalar(plhs[3]);

    Fcrosspairs(n1, x1, y1, t1, n2, x2, y2, t2,
            rmax, taumin, taumax, noutmax,
            nout, dxout, dyout, dtout, &status);

   dxout = mxRealloc(dxout, (*nout) * sizeof(double));
   dyout = mxRealloc(dyout, (*nout) * sizeof(double));
   dtout = mxRealloc(dtout, (*nout) * sizeof(double));

   mxSetPr(plhs[0], dxout);
   mxSetPr(plhs[1], dyout);
   mxSetPr(plhs[2], dtout);

   for (i=0; i<3; i++)
       mxSetN(plhs[i], *nout);
   mxSetN(plhs[3], 1);
}
  
