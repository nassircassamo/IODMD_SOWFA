/* DMODPO */
#include <mex.h>
#include "matrix.h"
#include "f2c.h"
#include "slicot.h"

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

void CalculateMatrices(double*R,mwSize n,mwSize l,mwSize m,mwSize s,mwSize N,
                       double*A,double*B,double*C,double*D,
                       double*K,mwSize fBD,mwSize fK)
{
    mwSize LDR,LIWORK,LDWORK;
    integer *IWORK,*BWORK,IWARN,INFO;
    mwSize LDW1,LDW2,LDW3;
    mwSize LIW1,LIW2;
    double *Q,*RY,*S,TOL,*DWORK;

    TOL= 0;
    LDR= max(max(2*(m+l)*s,3*m*s),5);

    if(fBD)
    {
        LIW1= max(l*s,max(m,1)*s+fK*n);
    }
    else
    {
        LIW1= max(m,1)*s+fK*n;
    }
    LIW2= fK*n*n;
    LIWORK= max(LIW1,LIW2);

    if(fBD)
    {
        LDW1= max(2*(l*s-l)*n+n*n+7*n,(l*s-l)*n+n+6*m*s);
        LDW1= max(LDW1,(l*s-l)*n+n+max(l+m*s,l*s+max(3*l*s,m)));
    }
    else
    {
        LDW1= max(2*(l*s-l)*n+2*n,(l*s-l)*n+n*n+7*n);
    }
    if(fK)
    {
        LDW2= max((l*s-l)*n+n+n*n+2*n+max(5*n,(2*m+l)*s+l),4*(max(m,1)*s+n));
        LDW2= max(LDW2,max(m,1)*s+2*n+l);
        LDW2= LDW2+l*s*n;
        LDW2= max(LDW2,l*s*n+max(m,1)*s*(n+l)*(max(m,1)*(n+l)+1)+max((n+l)*(n+l),4*max(m,1)*(n+l)+1));
        LDW3= max(4*n*n+2*n*l+l*l+max(3*l,n*l),14*n*n+12*n+5);
    }
    else
    {
        LDW2= 0;
        LDW3= 0;
    }
    LDWORK= max(LDW1,max(LDW2,LDW3));
    IWORK= mxMalloc(LIWORK*sizeof(integer));
    DWORK= mxMalloc(LDWORK*sizeof(double));
    BWORK= NULL;
    Q= NULL;
    RY= NULL;
    S= NULL;
    if(fK)
    {
        Q= mxMalloc(n*n*sizeof(double));
        RY= mxMalloc(l*l*sizeof(double));
        S= mxMalloc(n*l*sizeof(double));
        BWORK= mxMalloc(2*n*sizeof(integer));
    }

    if(fK)
    {
        if(fBD)
        {
            ib01bd("M","A","K",&s,&n,&m,&l,&N,R,&LDR,A,&n,C,&l,B,&n,
                   D,&l,Q,&n,RY,&l,S,&n,K,&n,&TOL,IWORK,DWORK,&LDWORK,
                   BWORK,&IWARN,&INFO);
        }
        else
        {
            ib01bd("M","C","K",&s,&n,&m,&l,&N,R,&LDR,A,&n,C,&l,B,&n,
                   D,&l,Q,&n,RY,&l,S,&n,K,&n,&TOL,IWORK,DWORK,&LDWORK,
                   BWORK,&IWARN,&INFO);
        }
        if(IWARN==5)
        {
            mexWarnMsgTxt("Covariance matrices are too small: returning K=0");
        }
    }
    else
    {
        if(fBD)
        {
            ib01bd("M","A","N",&s,&n,&m,&l,&N,R,&LDR,A,&n,C,&l,B,&n,
                   D,&l,Q,&n,RY,&l,S,&n,K,&n,&TOL,IWORK,DWORK,&LDWORK,
                   BWORK,&IWARN,&INFO);
        }
        else
        {
            ib01bd("M","C","N",&s,&n,&m,&l,&N,R,&LDR,A,&n,C,&l,B,&n,
                   D,&l,Q,&n,RY,&l,S,&n,K,&n,&TOL,IWORK,DWORK,&LDWORK,
                   BWORK,&IWARN,&INFO);
        }
    }
    if(IWARN==4)
    {
        mexWarnMsgTxt("Regression problem is rank-deficient");
    }

    if(fK)
    {
        mxFree(BWORK);
        mxFree(S);
        mxFree(RY);
        mxFree(Q);
    }
    mxFree(DWORK);
    mxFree(IWORK);
}


void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[])
{
    mwSize n,N,l,m,s,LDR;
    mwSize fBD,fK,count;
    double*A,*B,*C,*D,*K;
    double*R,*Rinternal,*nmat;

    if(nrhs!=2)
    {
        mexErrMsgTxt("DMODPO requires 2 input arguments");
    }
    if((nlhs<2)|(nlhs> 5))
    {
        mexErrMsgTxt("DMODPO requires 2, 3, 4 or 5 output arguments");
    }

    CheckMxReal(prhs[0],"R");
    CheckMxReal(prhs[1],"n");


    R= mxGetPr(prhs[0]);
    nmat= mxGetPr(prhs[1]);


    if(nmat==NULL)
    {
        mexErrMsgTxt("n must be specified.");
    }
    n= nmat[0];


    LDR= mxGetM(prhs[0]);
    if(LDR<5)
    {
        mexErrMsgTxt("R is too small.");
    }


    N= R[LDR-1];
    l= R[2*LDR-1];
    m= R[3*LDR-1];
    s= R[4*LDR-1];

    if(mxGetN(prhs[0])!=max((4*(m+l)*s),9))
    {
        mexErrMsgTxt("R has an incorrect size.");
    }

    if(!(s> 1))
    {
        mexErrMsgTxt("s must be at least 2.");
    }
    if((n<0)|(n>=s))
    {
        mexErrMsgTxt("n must be nonnegative and smaller than s.");
    }
    /*if(!(m> 0))
    {
        mexErrMsgTxt("DMODPO requires an input.");
    }*/
    if(!(l> 0))
    {
        mexErrMsgTxt("DMODPO requires an output.");
    }
    if(N<(2*(m+l+1)*abs(s)-1))
    {
        mexErrMsgTxt("The number of samples must be at least 2*(m+l+1)*s-1.");
    }
    if(LDR!=max(max(2*(m+l)*s,3*m*s),5))
    {
        mexErrMsgTxt("The size of R is incorrect.");
    }

    fBD= 0;fK= 0;
    switch(nlhs)
    {
    case 2:fBD= 0;fK= 0;break;
    case 3:fBD= 0;fK= 1;break;
    case 4:fBD= 1;fK= 0;break;
    case 5:fBD= 1;fK= 1;break;
    }

    Rinternal= mxCalloc(LDR*(2*(m+l)*s),sizeof(double));
    for(count= 0;count<(2*(m+l)*s)*LDR;count++)
    {
        Rinternal[count]= R[2*(m+l)*s*LDR+count];
    }
    A= NULL;B= NULL;C= NULL;D= NULL;K= NULL;
    switch(nlhs)
    {
    case 2:
        plhs[0]= mxCreateDoubleMatrix(n,n,mxREAL);
        plhs[1]= mxCreateDoubleMatrix(l,n,mxREAL);
        A= mxGetPr(plhs[0]);
        C= mxGetPr(plhs[1]);
        break;
    case 3:
        plhs[0]= mxCreateDoubleMatrix(n,n,mxREAL);
        plhs[1]= mxCreateDoubleMatrix(l,n,mxREAL);
        plhs[2]= mxCreateDoubleMatrix(n,l,mxREAL);
        A= mxGetPr(plhs[0]);
        C= mxGetPr(plhs[1]);
        K= mxGetPr(plhs[2]);
        break;
    case 4:
        plhs[0]= mxCreateDoubleMatrix(n,n,mxREAL);
        plhs[1]= mxCreateDoubleMatrix(n,m,mxREAL);
        plhs[2]= mxCreateDoubleMatrix(l,n,mxREAL);
        plhs[3]= mxCreateDoubleMatrix(l,m,mxREAL);
        A= mxGetPr(plhs[0]);
        B= mxGetPr(plhs[1]);
        C= mxGetPr(plhs[2]);
        D= mxGetPr(plhs[3]);
        break;
    case 5:
        plhs[0]= mxCreateDoubleMatrix(n,n,mxREAL);
        plhs[1]= mxCreateDoubleMatrix(n,m,mxREAL);
        plhs[2]= mxCreateDoubleMatrix(l,n,mxREAL);
        plhs[3]= mxCreateDoubleMatrix(l,m,mxREAL);
        plhs[4]= mxCreateDoubleMatrix(n,l,mxREAL);
        A= mxGetPr(plhs[0]);
        B= mxGetPr(plhs[1]);
        C= mxGetPr(plhs[2]);
        D= mxGetPr(plhs[3]);
        K= mxGetPr(plhs[4]);
        break;
    }

    if(n> 0)
    {
        CalculateMatrices(Rinternal,n,l,m,s,N,A,B,C,D,K,fBD,fK);
    }
    else
    {
        if(fBD)
        {
            mexPrintf("WARNING: Estimating D in 0th order system is not functional.\n");
        }
    }

    mxFree(Rinternal);
}
