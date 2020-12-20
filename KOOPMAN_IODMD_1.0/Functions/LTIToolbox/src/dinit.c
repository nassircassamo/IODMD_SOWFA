/* DINIT */
#include <mex.h>
#include "matrix.h"
#include "f2c.h"
#include "slicot.h"

void CheckMxValid(const mxArray*theArray,char*ArrayName)
{
    if (!mxIsNumeric(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be numeric.");
    }
    if (!mxIsDouble(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be double precision.");
    }
    if (mxIsSparse(theArray))
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

void CalculateX0(double*u,double*y,double*A,double*B,double*C,
                 double*D,double*x,mwSize N,mwSize n,mwSize l,mwSize m)
{
    mwSize LDWORK,LIWORK;
    integer *IWORK,IWARN,INFO;
    double *V,TOL,*DWORK;
    mwSize LDWa,LDWb,LDWc,LDWd,LDWf,LDWr,LDWt;
    TOL = 0; IWARN = 0; INFO = 0;

    LIWORK = n;

    LDWt = N;
    LDWa = n;
    LDWr = n*m+LDWa;
    LDWb = l*m;
    LDWc = l*n;
    LDWd = 2*n*n+n;
    LDWf = m+max(2*LDWr,m);
    LDWORK = LDWt*l*(LDWr+1)+(LDWb+LDWr)*(LDWr+1)+n*n*m+LDWc+max(LDWd,LDWf);

    IWORK = mxMalloc(LIWORK*sizeof(integer));
    DWORK = mxMalloc(LDWORK*sizeof(double));

    V = mxMalloc(n*n*sizeof(double));

    ib01cd("X","U","D",&n,&m,&l,&N,A,&n,B,&n,C,&l,D,&l,u,&N,y,&N,
            x,V,&n,&TOL,IWORK,DWORK,&LDWORK,&IWARN,&INFO);
    switch(IWARN)
    {
    case 6:
        mexPrintf("DINIT: A is unstable, x0 results might be inaccurate.\n");
    }

    mxFree(V);
    mxFree(DWORK);
    mxFree(IWORK);
}

void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[])
{
    double *u,*y,*x;
    double *A,*B,*C,*D;
    mwSize N,n,l,m;

    if(nrhs!=6)
    {
        mexErrMsgTxt("DINIT requires 6 input arguments.");
    }
    if(nlhs!=1)
    {
        mexErrMsgTxt("DINIT requires 1 output argument.");
    }

    CheckMxReal(prhs[0],"A");
    CheckMxReal(prhs[1],"B");
    CheckMxReal(prhs[2],"C");
    CheckMxReal(prhs[3],"D");
    CheckMxReal(prhs[4],"u");
    CheckMxReal(prhs[5],"y");

    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    D = mxGetPr(prhs[3]);
    u = mxGetPr(prhs[4]);
    y = mxGetPr(prhs[5]);

    N = mxGetM(prhs[4]);
    m = mxGetN(prhs[4]);
    l = mxGetN(prhs[5]);
    n = mxGetM(prhs[0]);

    if (mxGetN(prhs[0])!=n)
    {
        mexErrMsgTxt("A must be square.");
    }
    if (mxGetM(prhs[1])!=n)
    {
        mexErrMsgTxt("B has a wrong number of rows.");
    }
    if (mxGetN(prhs[1])!=m)
    {
        mexErrMsgTxt("B has a wrong number of columns.");
    }
    if (mxGetM(prhs[2])!=l)
    {
        mexErrMsgTxt("C has a wrong number of rows.");
    }
    if (mxGetN(prhs[2])!=n)
    {
        mexErrMsgTxt("C has a wrong number of columns.");
    }
    if (mxGetM(prhs[3])!=l)
    {
        mexErrMsgTxt("D has a wrong number of rows.");
    }
    if (mxGetN(prhs[3])!=m)
    {
        mexErrMsgTxt("D has a wrong number of columns.");
    }
    if (mxGetM(prhs[5])!=N)
    {
        mexErrMsgTxt("The input and output must have an equal number of samples.");
    }

    if (l<1)
    {

        mexErrMsgTxt("DINIT requires an output.");
    }
    if (m<1)
    {
        mexErrMsgTxt("DINIT requires an input.");
    }
    if (N<n)
    {
        mexErrMsgTxt("The number of samples is too small.");
    }

    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    x = mxGetPr(plhs[0]);

    if (n>0) {
		CalculateX0(u,y,A,B,C,D,x,N,n,l,m);
	}

}
