/* SIMLNS */
#include <math.h>
#include <mex.h>
#include "f2c.h"
#include "matrix.h"
#include "slicot.h"

void CheckMxReal(const mxArray*theArray,char*ArrayName)
{
    if (!mxIsNumeric(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be numeric.");
    }
    if (!mxIsDouble(theArray)) {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be double precision.");
    }    
    if (mxIsSparse(theArray)) {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must not be sparse.");
    }
    if (mxIsComplex(theArray)) {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be real.");
    }
}

void simmaplns(double *A, double *B, double *C, double *K, mwSize n, mwSize l, mwSize m,
               double *U2, mwSize fD, mwSize fx)
{
    mwSize LDU,LWORK;
    double *M;
    double *TAU, *WORK, OPTWORK;
    integer INFO;
    mwSize SizeU2;
    mwSize fK;
    mwSize SizeR;
    mwSize SizeC;
    mwSize StartRowK;
    mwSize i,j,k;

    fK = (K==NULL)?0:1;
    SizeR = n*n + n*m + l*n + fD*l*m + fK*n*l;
    SizeC = n*n;
    M = mxCalloc(SizeR*SizeC,sizeof(double));

    for (i= 0;i<n;i++)
    {
        for (j= 0;j<n;j++)
        {
            for (k= 0;k<n;k++)
            {
                M[i*(n+l)+k+(i*n+j)*SizeR]= A[k+n*j];
            }
            for (k= 0;k<l;k++)
            {
                M[i*(n+l)+n+k+(i*n+j)*SizeR]= C[k+l*j];
            }
        }
    }

    for (i= 0;i<n;i++)
    {
        for (j= 0;j<n;j++)
        {
            for (k= 0;k<n;k++)
            {
                M[j*(n+l)+k+(i*n+k)*SizeR]-= A[i+n*j];
            }
        }
        for (j= 0;j<m;j++)
        {
            for (k= 0;k<n;k++)
            {
                if(fD)
                {
                    M[(n+j)*(n+l)+k+(i*n+k)*SizeR]-= B[i+n*j];
                }
                else
                {
                    M[n*(n+l)+j*n+k+(i*n+k)*SizeR]-= B[i+n*j];
                }
            }
        }
    }

    if (fK == 1)
    {
        StartRowK = n*n+n*m+l*n+fD*l*m;
        for (i=0;i<n;i++)
        {
            for (j=0;j<l;j++)
            {
                for (k=0;k<n;k++)
                {
                    M[j*n+StartRowK+k+(i*n+k)*SizeR] = -K[i+n*j];
                }
            }
        }
    }

    TAU =  mxMalloc(SizeR*sizeof(double));
    LWORK = -1;
    dgeqrf(&SizeR,&SizeC,M,&SizeR,TAU,&OPTWORK,&LWORK,&INFO);
    
    LWORK = ceil(OPTWORK);
    WORK = mxMalloc(LWORK*sizeof(double));
    dgeqrf(&SizeR,&SizeC,M,&SizeR,TAU,WORK,&LWORK,&INFO);

    LDU = SizeR + fx*n;
    for (i=0; i<SizeR-SizeC; i++)
    {
        U2[SizeC+i+i*LDU] = 1;
    }
    SizeU2 = SizeR-SizeC;
    LWORK = -1;
    dormqr("L","N",&SizeR,&SizeU2,&SizeC,M,&SizeR,TAU,U2,&LDU,&OPTWORK,&LWORK,&INFO);
    
    LWORK = ceil(OPTWORK);
    mxFree(WORK);
    WORK = mxMalloc(LWORK*sizeof(double));
    dormqr("L","N",&SizeR,&SizeU2,&SizeC,M,&SizeR,TAU,U2,&LDU,WORK,&LWORK,&INFO);

    mxFree(WORK);
    mxFree(TAU);
    mxFree(M);

    if (fx == 1)
    {
        for (i=0; i<n; i++)
        {
            U2[SizeR+i+(SizeR-SizeC+i)*LDU] = 1;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C, *K, *U2, fDind, fxind;
    mwSize n,l,m;
    mwSize fD, fK, fx;
    mwSize i;

    if(nrhs!=6) {
        mexErrMsgTxt("SIMLNS requires 6 input arguments");
    }

    CheckMxReal(prhs[0],"A");
    CheckMxReal(prhs[1],"B");
    CheckMxReal(prhs[2],"C");
    if (!mxIsEmpty(prhs[3])) 
    {
        CheckMxReal(prhs[3],"K");
    }

    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    K = mxGetPr(prhs[3]);
    fDind = mxGetScalar(prhs[4]);
    fxind = mxGetScalar(prhs[5]);

    n = mxGetM(prhs[0]);
    m = mxGetN(prhs[1]);
    l = mxGetM(prhs[2]);

    if (mxGetN(prhs[0])!=n){
        mexErrMsgTxt("A matrix must be square.");
    }
    if (mxGetM(prhs[1])!=n){
        mexErrMsgTxt("B matrix has wrong number of rows.");
    }
    if (mxGetN(prhs[2])!=n){
        mexErrMsgTxt("C matrix has wrong number of columns.");
    }
    if (!mxIsEmpty(prhs[3]))
    {
        if (mxGetM(prhs[3])!=n)
        {
            mexErrMsgTxt("K matrix has wrong number of rows.");
        }
        if (mxGetN(prhs[3])!=l)
        {
            mexErrMsgTxt("K matrix has wrong number of columns.");
        }
        fK = 1;
    }
    else
    {
        fK = 0;
    }
    if (!mxIsEmpty(prhs[4]))
    {
        if ((fDind > 0.99)&(fDind < 1.01))
        {
            fD = 1;
        }
        else if ((fDind > -0.01)&(fDind < 0.01))
        {
            fD = 0;
        }
        else
        {
            mexErrMsgTxt("fD should be empty, 0 or 1.");
        }
    }
    else
    {
        fD = 0;
    }

    if (!mxIsEmpty(prhs[5]))
    {
        if ((fxind > 0.99)&(fxind < 1.01))
        {
            fx = 1;
        }
        else if ((fxind > -0.01)&(fxind < 0.01))
        {
            fx = 0;
        }
        else
        {
            mexErrMsgTxt("fx should be empty, 0 or 1.");
        }
    }
    else
    {
        fx = 0;
    }

    plhs[0] = mxCreateDoubleMatrix(n*n + n*m + l*n + fD*l*m + fK*l*n + fx*n,
                                  n*m + l*n + fD*l*m + fK*l*n + fx*n, mxREAL);

    U2 = mxGetPr(plhs[0]);

    if (n == 0)
    {
        if (fD == 1)
        {
            for (i=0; i<l*m; i++)
            {
                U2[i+i*l*m] = 1;
            }
        }
    }
    else
    {
        if (fK == 1)
        {
            simmaplns(A,B,C,K,n,l,m,U2,fD,fx);
        }
        else
        {
            simmaplns(A,B,C,NULL,n,l,m,U2,fD,fx);
        }
    }
}

