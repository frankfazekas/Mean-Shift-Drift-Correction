#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define OK 0
#define ERR_OVERFLOW 1

double sqrt();

void Fcrosspairs_3D(n1, x1, y1, z1, t1, n2, x2, y2, z2, t2, rmax, taumin, taumax,
        noutmax, nout, dxout, dyout, dzout, dtout, status)
     /* inputs */
     int n1, n2, noutmax;
     double *x1, *y1, *z1, *x2, *y2, *z2, rmax, *t1, *t2, taumin, taumax;
     /* outputs */
     int *nout, *status;
     double *dxout, *dyout, *dzout, *dtout;
     
{
  int k, kmax, i, j, jleft;
  double x1i, y1i, z1i, r2max, xleft, dx, dy, dz, dx2, d2, dt;
 
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

      for (; i < n1; i++) {

          x1i = x1[i];
          y1i = y1[i];
          z1i = z1[i];

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
              dz = z2[j] - z1i;
              d2 = dx2 + dy * dy + dz * dz;
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
                  dzout[k] = dz;
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
    double *x1, *y1, *z1, *t1, *x2, *y2, *z2, *t2;
    double rmax, taumin, taumax;
    double *dxout, *dyout, *dzout, *dtout;
    int *nout, noutmax, status;
    int i;

    /* RHS args are:
     *  0   x1 // coords for first pp
     *  1   y1
     *  2   z1
     *  3   t1
     *  4   x2 // coords for second pp
     *  5   y2
     *  6   z2
     *  7   t2
     *  8   rmax
     *  9   taumin
     *  10   taumax
     *  11   noutmax
     */

    /* LHS args are:
     * 0,1,2,3    dxout, dyout, dzoutdtout
     * 4    optional error
     */

     /* Check for proper number of arguments. */
    if (nrhs != 12) {
     mexErrMsgIdAndTxt("MATLAB:timestwo:invalidNumInputs",
      "Twelve inputs required.");
    }
    else if (nlhs > 5) {
     mexErrMsgIdAndTxt("MATLAB:timestwo:maxlhs",
      "Too many output arguments.");
    }

    n1 = mxGetNumberOfElements(prhs[0]);
    n2 = mxGetNumberOfElements(prhs[4]);

    x1 = (double *) mxGetPr(prhs[0]);
    y1 = (double *) mxGetPr(prhs[1]); 
    z1 = (double *) mxGetPr(prhs[2]);
    t1 = (double *) mxGetPr(prhs[3]);
    x2 = (double *) mxGetPr(prhs[4]);
    y2 = (double *) mxGetPr(prhs[5]);
    z2 = (double *) mxGetPr(prhs[6]);
    t2 = (double *) mxGetPr(prhs[7]);

    rmax = mxGetScalar(prhs[8]);
    taumin = mxGetScalar(prhs[9]);
    taumax = mxGetScalar(prhs[10]);

    nout = (int *) mxGetPr(prhs[11]);
    noutmax = *nout;

    for (i=0; i<4; i++)
        plhs[i] = mxCreateNumericMatrix(1, noutmax, mxDOUBLE_CLASS, mxREAL);
    plhs[4] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    dxout = (double *) mxGetData(plhs[0]);
    dyout = (double *) mxGetData(plhs[1]);
    dzout = (double *) mxGetData(plhs[2]);
    dtout = (double *) mxGetData(plhs[3]);
    status = (int) mxGetScalar(plhs[4]);

    Fcrosspairs_3D(n1, x1, y1, z1, t1, n2, x2, y2, z2, t2,
            rmax, taumin, taumax, noutmax,
            nout, dxout, dyout, dzout, dtout, &status);

    for (i = 0; i < 4; i++)
        mxSetN(plhs[i], *nout);
    mxSetN(plhs[4], 1);

}
  
