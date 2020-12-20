/* LTIITR */
#include <mex.h>
#include "f2c.h"
#include "matrix.h"

void CheckMxValid(const mxArray*theArray,char*ArrayName)
{
    if(!mxIsNumeric(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be numeric.");
    }
    if(!mxIsDouble(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be double precision.");
    }
    if(mxIsSparse(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must not be sparse.");
    }
}

void CheckMxReal(const mxArray*theArray,char*ArrayName)
{
    CheckMxValid(theArray,ArrayName);
    if(mxIsComplex(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be real.");
    }
}

void ltiitr(double*A, double*B, double*u, double*w, double*xzero, mwSize n, 
            mwSize m, mwSize N, double*xt, double*tt, double*x)
{
    double t1,t2;
    mwSize i,j,k,l;

    for(l= 0;l<n;l++){
        xt[l]= xzero[l];
    }
    for(k= 0;k<N;k++){

        for(j= 0;j<n;j++){
            x[k+j*N]= xt[j];
        }
        for(i= 0;i<n;i++){
            t1= 0;
            t2= 0;
            for(l= 0;l<n;l++){
                t1+= (A[l*n+i]*xt[l]);
            }
            for(l= 0;l<m;l++){
                t2+= (B[l*n+i]*u[l*N+k]);
            }
            tt[i]= t1+t2;
        }

        for(l= 0;l<n;l++){
            xt[l]= tt[l]+w[l*N+k];
        }
    }
}

void pltiitr(double*A, double*B, double*u, double*xzero, mwSize n, mwSize m, mwSize N,
            double*xt, double*tt, double*x)
{
    double t1,t2;
    mwSize i,j,k,l;

    for(l= 0;l<n;l++){
        xt[l]= xzero[l];
    }

    for(k= 0;k<N;k++){

        for(j= 0;j<n;j++){
            x[k+j*N]= xt[j];
        }
        for(i= 0;i<n;i++){
            t1= 0;
            t2= 0;
            for(l= 0;l<n;l++){
                t1+= (A[l*n+i]*xt[l]);
            }
            for(l= 0;l<m;l++){
                t2+= (B[l*n+i]*u[l*N+k]);
            }
            tt[i]= t1+t2;
        }

        for(l= 0;l<n;l++){
            xt[l]= tt[l];
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
  /* Declare pointers used for int arguments of lowess */
  double *A, *B, *u, *w, *xzero, *x;
  double *xt, *tt;
  mwSize n, N, m ,i;  

  if(nrhs != 5) {
      mexErrMsgTxt("LTIITR requires 5 input arguments");
  }
  
  CheckMxReal(prhs[0],"A");
  CheckMxReal(prhs[1],"B");
  CheckMxReal(prhs[2],"u");
  CheckMxReal(prhs[3],"w");
  CheckMxReal(prhs[4],"x0");
  
  /* assign pointers to the input arguments */
  A     = mxGetPr(prhs[0]);
  B     = mxGetPr(prhs[1]);
  u     = mxGetPr(prhs[2]);
  xzero = mxGetPr(prhs[4]);
   
  /* Get dimensions */
  n = mxGetM(prhs[0]);
  N = mxGetM(prhs[2]);
  m = mxGetN(prhs[2]);
  
  /* Check input arguments */
  if (mxGetN(prhs[0])!= n) {
    mexErrMsgTxt("A matrix has wrong number of columns.");
  }
  if (mxGetM(prhs[1])!= n) {
    mexErrMsgTxt("B matrix has wrong number of rows.");
  }
  if (mxGetN(prhs[1])!= m) {
    mexErrMsgTxt("B matrix has wrong number of columns.");
  }
  if ((!mxIsEmpty(prhs[3])) && (mxGetN(prhs[3])!= n)) {
    mexErrMsgTxt("W has wrong number of columns.");
  } 
  if (!mxIsEmpty(prhs[4])){
      if (mxGetM(prhs[4])!= n){
          mexErrMsgTxt("X0 and A must have the same number of rows.");
      }
      if (mxGetN(prhs[4])!= 1){
          mexErrMsgTxt("X0 can only have one column.");
      }
  }
  else {
      xzero = mxCalloc(n,sizeof(double));
      for(i=0; i<n; i++){
          xzero[i] = 0;
      }
  }

  /* Create matrix for the output arguments and assign a pointer */
  plhs[0] = mxCreateDoubleMatrix(N,n,mxREAL);
  x = mxGetPr(plhs[0]);
   
  /* Allocate memmory */ 
  xt = mxCalloc(n,sizeof(double));
  tt = mxCalloc(n,sizeof(double));
  
  /* Call the function */
  if (!mxIsEmpty(prhs[3])) {
    w = mxGetPr(prhs[3]);
    ltiitr(A,B,u,w,xzero,n,m,N,xt,tt,x);
  }
  else {
    pltiitr(A,B,u,xzero,n,m,N,xt,tt,x);
  }

  /* Free allocated memory */
  mxFree(xt);
  mxFree(tt);
}

