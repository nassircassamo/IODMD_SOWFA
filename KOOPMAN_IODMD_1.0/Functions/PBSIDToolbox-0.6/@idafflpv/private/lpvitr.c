/*
 * lpvitr.c 
 *
 * Simulation of the LPV state equation
 *
 * Vincent Verdult, July 2000
 * copyright (c) 2000 Vincent Verdult
 */

#include "mex.h"
#define P_mxCalloc(N,T) (T*)mxCalloc((mwSize)N,sizeof(T))

/* Function: lpvitr
 *
 * Compute for k = 1...N
 *   x(k+1) = A(:,1:n)*x(k) + A(:,n+1:n*(s+1))*kron(p(k),x(k))
 *          + B(:,1:m)*u(k) + B(:,m+1:m*(s+1))*kron(p(k),u(k)) + w(k)
 *
 *
 * This is a fast implementation of the following MATLAB code
 *
 * for k=1:N
 *   x(:,k)=x0;
 *   for i=1:s
 *     px((i-1)*n+1:i*n,1)=p(k,i)*x0;
 *     pu((i-1)*n+1:i*n,1)=p(k,i)*u(k,:)';
 *   end
 *   x0=A(:,1:n)*x0+A(:,n+1:n*(s+1))*px+B(:,1:m)*u(k,:)'+B(:,m+1:m*(s+1))*pu+w(k,:)';
 * end
 *  x=x';
 * 
 * The matrices are stored in an array, columnwise.
 * For an m by n matrix C, the (i,j) entry is given by C(j*m+i).
 *
 * auxiliary variables: px             stores kron(p(k),x(k))
 *                      pu             stores kron(p(k),u(k))
 *                      xt             state at a certain time instant
 *                      tt             state at a certain time instant
 */

void lpvitr(double *A, double *B, double *p, double *u, 
            double *w, double *xzero, mwSize n, mwSize m, mwSize s, mwSize N, 
            double *px, double *pu, double *xt, double *tt, double *x)
{ 
  mwSize i,j,k,l;
  
  /* Copy contents of xzero to xt */
  for (l=0; l<n; l++) {
    xt[l]=xzero[l];
  }
 
  for (k=0; k<N; k++) {
    /* Store the state at time instant k */
    for (j=0; j<n; j++) { 
      x[k+j*N] = xt[j]; 
    }     
    /* Compute px = kron(p(k),x(k)) */
    for (i=0; i<s; i++) {
      for (j=0; j<n; j++) {
       px[i*n+j] = p[i*N+k] * xt[j];
      }   
    } 
    /* Compute pu = kron(p(k),u(k)) */  
    for (i=0; i<s; i++) {
      for (j=0; j<m; j++) {
       pu[i*m+j] = p[i*N+k] * u[j*N+k];
      }   
    } 
    /* Compute xt */
    for (i=0; i<n; i++) {
      tt[i] = w[i*N+k];
      for (l=0; l<n; l++) {
        tt[i] += (A[l*n+i] * xt[l]);
      }        
      for (l=0; l<s*n; l++) {
        tt[i] += (A[n*n+l*n+i] * px[l]);
      }       
      for (l=0; l<m; l++) {
        tt[i] += (B[l*n+i] * u[l*N+k]);
      }  
      for (l=0; l<s*m; l++) {
        tt[i] += (B[n*m+l*n+i] * pu[l]);
      }        
    }
    for (l=0; l<n; l++) {
      xt[l]=tt[l];
    }
  }
}   

/* Function: plpvitr
 *
 * Same as lpvitr, but without w.
 */

void plpvitr(double *A, double *B, double *p, double *u, 
            double *xzero, mwSize n, mwSize m, mwSize s, mwSize N, 
            double *px, double *pu, double *xt, double *tt, double *x)
{ 
  mwSize i,j,k,l;
  
  /* Copy contents of xzero to xt */
  for (l=0; l<n; l++) {
    xt[l]=xzero[l];
  }
 
  for (k=0; k<N; k++) {
    /* Store the state at time instant k */
    for (j=0; j<n; j++) { 
      x[k+j*N] = xt[j]; 
    }     
    /* Compute px = kron(p(k),x(k)) */
    for (i=0; i<s; i++) {
      for (j=0; j<n; j++) {
       px[i*n+j] = p[i*N+k] * xt[j];
      }   
    } 
    /* Compute pu = kron(p(k),u(k)) */  
    for (i=0; i<s; i++) {
      for (j=0; j<m; j++) {
       pu[i*m+j] = p[i*N+k] * u[j*N+k];
      }   
    } 
    /* Compute xt */
    for (i=0; i<n; i++) {
      tt[i] = 0;
      for (l=0; l<n; l++) {
        tt[i] += (A[l*n+i] * xt[l]);
      }        
      for (l=0; l<s*n; l++) {
        tt[i] += (A[n*n+l*n+i] * px[l]);
      }       
      for (l=0; l<m; l++) {
        tt[i] += (B[l*n+i] * u[l*N+k]);
      }  
      for (l=0; l<s*m; l++) {
        tt[i] += (B[n*m+l*n+i] * pu[l]);
      }        
    }
    for (l=0; l<n; l++) {
      xt[l]=tt[l];
    }
  }
}   


/************************************ 
 * MATLAB Interface routine 
 ************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
  /* Declare pointers used for integer arguments of lowess */
  double *A, *B, *p, *u, *w, *xzero, *x;
  double *px, *pu, *xt, *tt;
  mwSize n, N, s, m;  

  /* assign pointers to the input arguments */
  A          = mxGetPr(prhs[0]);
  B          = mxGetPr(prhs[1]);
  p          = mxGetPr(prhs[2]); 
  u          = mxGetPr(prhs[3]);
  xzero      = mxGetPr(prhs[5]);
  
  /* Get dimensions */
  n = mxGetM(prhs[0]);
  N = mxGetM(prhs[2]);
  s = mxGetN(prhs[2]);
  m = mxGetN(prhs[3]);
  
  /* Check input arguments */
  if (mxGetN(prhs[0])!=n*(s+1)) {
    mexErrMsgTxt("A matrix has wrong number of columns.");
  }
  if (mxGetM(prhs[1])!=n) {
    mexErrMsgTxt("B matrix has wrong number of rows.");
  }
  if (mxGetN(prhs[1])!=m*(s+1)) {
    mexErrMsgTxt("B matrix has wrong number of columns.");
  }
  if (mxGetM(prhs[3])!=N) {
    mexErrMsgTxt("P and U must have the same number of rows.");
  }
  if ((!mxIsEmpty(prhs[4])) && (mxGetN(prhs[4])!=n)) {
    mexErrMsgTxt("W has wrong number of columns.");
  }
  if ((!mxIsEmpty(prhs[4])) && (mxGetM(prhs[4])!=N)) {
    mexErrMsgTxt("P and W must have the same number of rows.");
  }
  if (mxGetM(prhs[5])!=n) {
    mexErrMsgTxt("X0 and A must have the same number of rows.");
  }
  if (mxGetN(prhs[5])!=1) {
    mexErrMsgTxt("X0 can only have one column.");
  }


  /* Create matrix for the output arguments and assign a pointer */
  plhs[0]=mxCreateDoubleMatrix(N,n,mxREAL);
  x = mxGetPr(plhs[0]);
   
  /* Allocate memmory */ 
  px = P_mxCalloc(n*s,double);
  pu = P_mxCalloc(m*s,double);
  xt = P_mxCalloc(n,double);
  tt = P_mxCalloc(n,double);
  
  /* Call the function */
  if (!mxIsEmpty(prhs[4])) {
     w = mxGetPr(prhs[4]);
     lpvitr(A,B,p,u,w,xzero,n,m,s,N,px,pu,xt,tt,x);
  }
  else {
    plpvitr(A,B,p,u,xzero,n,m,s,N,px,pu,xt,tt,x);
  }

  /* Free allocated memory */
  mxFree(px);
  mxFree(pu);
  mxFree(xt);
  mxFree(tt);
}








