#include <mex.h>
#include "f2c.h"
#include "matrix.h"
#include "slicot.h"

/*
  FASTSVD is a mex-file interface to the SLICOT routime MB03UD. The
    objective of FASTSVD is to compute all, or part, of the singular
    value decomposition of a real upper triangular matrix.

   Karel Hinnen, 2004
*/

/* Declaration of max_i subroutine */
mwSize max_i(i,j) { if(i>j) return i; else return j; }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Declaration of workspaces, warning and info params */
  double *A;
  mwSize LDA;
  double *Q, *S;
  double *DWORK;
  mwSize LDWORK;
  mwSize N;
  integer INFO  = 0;

  /* Check is the nuber of input arguments is right */
  if ( (nrhs < 1) || (nrhs > 1) ) {
    mexErrMsgTxt("Error calling fastsvd: wrong number of input arguments.\n");
    return;
  }

  A   = mxGetPr(prhs[0]);
  LDA = mxGetM(prhs[0]);
  N   = mxGetN(prhs[0]);

  /* Check of the input matrix is square */
  if (LDA != N){
    mexErrMsgTxt("Error calling fastsvd: the matrix should be square. \n");
  }

  /* Free memory for left singular matrix and singular values */
  plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
  Q       = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
  S       = mxGetPr(plhs[1]);

  /* Allocating of workspace */
  LDWORK = max_i(1,10*N);
  DWORK  = (double *)mxMalloc(LDWORK*sizeof(double));

  /*  Calculation of left singular matrix and singular values: */
  mb03ud("V","N",&N,A,&LDA,Q,&N,S,DWORK,&LDWORK,&INFO);

  if (INFO > 0){
    printf("the SVD algorithm has failed to converge \n");
  }

  mxFree(DWORK);
}

