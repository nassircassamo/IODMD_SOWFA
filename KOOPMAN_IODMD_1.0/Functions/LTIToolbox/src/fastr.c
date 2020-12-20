#include <mex.h>
#include "f2c.h"
#include "matrix.h"
#include "slicot.h"

/*
   FASTR is a mex-function that is used as a subroutine in an efficient
   implementation of the SSARX algorithm. The routine uses standard functions
   of the SLICOT library to built a Hankel matrix from the measurement data
   and to compute the R factor.

   Matlab SYNTAX
   R = fastr(u,y,s,method);

   Karel Hinnen, October 2005
 */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Declaration of required local variables */
    char *string;
    mwSize strlen;
    char meth = 'N';
    char alg  = 'C';
    mwSize s, M, L, NSMP;
    double *U; mwSize LDU;
    double *Y; mwSize LDY;
    double *R; mwSize LDR;
    double *d_info;

    integer *IWORK;
    double *DWORK;
    mwSize LDWORK, LIWORK;
    mwSize NRMAX;
    integer IWARN = 0;
    integer INFO  = 0;

    /* Check is the number of input arguments is right*/
    if ( (nrhs < 3) || (nrhs > 4) ) {
        mexErrMsgTxt("Error calling fastr: wrong number of arguments.\n");
        return;
    }

    /* Obtain the input and output signal u and y from workspace */
    U    = mxGetPr(prhs[0]);
    LDU  = mxGetM(prhs[0]);
    M    = mxGetN(prhs[0]);
    Y    = mxGetPr(prhs[1]);
    LDY  = mxGetM(prhs[1]);
    L    = mxGetN(prhs[1]);
    if (M == 0) {
        NSMP = LDY;
        LDU  = 1;
    }
    else {
        NSMP = min(LDU,LDY);
    }
    s = (integer)*mxGetPr(prhs[2]);


    /* Check is there are sufficient samples to build the Hankel matrix */
    if ( NSMP < 2*(M+L+1)*s-1 ) {
        mexErrMsgTxt("Number of samples too small for the chosen number of"
        "blockrows.\n");
        return;
    }

    /* Determine the type of algorithm that is used in computing the RQ factor */
    if (nrhs == 4) {
        strlen = mxGetM(prhs[3])*mxGetN(prhs[3])*sizeof(mxChar)+1;
        string = (char *)mxMalloc(strlen);
        if (mxGetString(prhs[3],string,strlen) > 0) {
            mexErrMsgTxt("Invalid argument. Fourth argument should be of type char.\n");
            return;
        }
        alg = string[0];
        mxFree(string);
    }
    if ( (alg != 'C') && (alg != 'F') && (alg != 'Q') ) {
        mexErrMsgTxt("Algorithm type should be 'C','F' or 'Q'.\n");
        return;
    }
    mexPrintf("Specified algorithm: %c \n",alg);

    /* Free memory for R matrix: */
    LDR = 2*(M+L)*s;
    plhs[0] = mxCreateDoubleMatrix(LDR,LDR,mxREAL);
    R = mxGetPr(plhs[0]);

    if (nlhs > 1){
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        d_info  = mxGetPr(plhs[1]);
    }

    /* Allocating workspace: integer workspace */
    if (alg == 'F') {
        LIWORK = M+L;
    }
    else {
        LIWORK = 0;
    }
    IWORK = (integer*)mxMalloc(LIWORK*sizeof(integer));

    /* Free memory as workspace */
    NRMAX  = (M+L)*s;

    if ( alg == 'C')
        LDWORK = 1;
    else if (alg == 'F')
        LDWORK = 4*NRMAX*(M+L+1)+2*NRMAX;
    else
        LDWORK = 6*NRMAX;

    /*  printf("Number of doubles in LDWORK: %d \n",LDWORK); */
    DWORK = (double *)mxMalloc(LDWORK*sizeof(double));

    /*  Calculation: */
    ib01md(&meth,&alg,"O","C",&s,&M,&L,&NSMP,U,&LDU,Y,&LDY,R,&LDR,IWORK,DWORK,&LDWORK,&IWARN,&INFO);
    
    if (IWARN == 2){
        printf("WARNING: a fast algorithm was requested (ALG = 'C' or 'F') "
        "but failed; \n the QR algorithm has been used instead. \n");
    }
    if (nlhs > 1){
        d_info[0] = (double)INFO;
    }

    mxFree(IWORK);
    mxFree(DWORK);
}
