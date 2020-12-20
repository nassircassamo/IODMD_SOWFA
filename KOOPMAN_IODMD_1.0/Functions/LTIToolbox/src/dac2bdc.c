/* DAC2BDC */
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
    if (mxIsComplex(theArray))
    {
        mexPrintf("Input parameter %s is invalid:\n",ArrayName);
        mexErrMsgTxt("Array must be real.");
    }
}

void EstimateD(mwSize l, mwSize m, mwSize N, double*u, double*y, double*D)
{
    mwSize count,count2;
    mwSize LDWORK;
    integer INFO;
    double *utemp, *ytemp;
    double *DWORK;

    utemp = mxMalloc(N*m*sizeof(double));
    ytemp = mxMalloc(N*l*sizeof(double));

    for (count= 0;count<N*m;count++) {
		utemp[count] = u[count];
	}
    for (count= 0;count<N*l;count++) {
		ytemp[count] = y[count];
	}

    LDWORK = max(1,min(N,m)+max(min(N,m),l));
    DWORK = mxMalloc(LDWORK*sizeof(double));

    dgels("N",&N,&m,&l,utemp,&N,ytemp,&N,DWORK,&LDWORK,&INFO);

    for(count= 0;count<l;count++) {
        for(count2= 0;count2<m;count2++) {
            D[count+l*count2] = ytemp[count2+N*count];
		}
	}

    mxFree(DWORK);
    mxFree(ytemp);
    mxFree(utemp);
}

void CalculateBD(double*u,double*y,double*A,double*B,
                 double*C,double*D,mwSize N,mwSize l,mwSize m,
                 mwSize n,mwSize fD)
{
    mwSize count;
    mwSize LDWORK, LIWORK;
    integer *IWORK, IWARN, INFO;
    double *V, TOL, *DWORK, *utemp;
    double *x;
    TOL = 0; IWARN = 0; INFO = 0;

    LIWORK = max(n*m+n,m);
    LDWORK = N*l*(n*m+n+l*m+1)+max(n+max(2*n*n+n,m+max(2*(n*m+n),m)),6*(n*m+n));
    LDWORK = max(LDWORK,N*l*(n*m+n+1)+2*m*m+6*m);
    IWORK = mxMalloc(LIWORK*sizeof(integer));
    DWORK = mxMalloc(LDWORK*sizeof(double));
    V = mxMalloc(n*n*sizeof(double));
    x = mxMalloc(n*sizeof(double));
    utemp = mxMalloc(N*m*sizeof(double));
    for (count= 0;count<N*m;count++) {
        utemp[count]= u[count];
	}

    if(fD)
    {
        ib01cd("X","C","D",&n,&m,&l,&N,A,&n,B,&n,C,&l,D,&l,utemp,&N,y,&N,
                x,V,&n,&TOL,IWORK,DWORK,&LDWORK,&IWARN,&INFO);
    }
    else
    {
        ib01cd("X","C","B",&n,&m,&l,&N,A,&n,B,&n,C,&l,D,&l,utemp,&N,y,&N,
                x,V,&n,&TOL,IWORK,DWORK,&LDWORK,&IWARN,&INFO);
    }

    mxFree(utemp);
    mxFree(x);
    mxFree(V);
    mxFree(DWORK);
    mxFree(IWORK);
}

void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[])
{
    double *u, *y;
    double *A, *B, *C, *D;
    mwSize N,l,m,n;
    mwSize fD;

    if (nrhs!=4)
    {
        mexErrMsgTxt("DAC2BDC requires 4 input arguments.");
    }
    if ((nlhs<1)||(nlhs> 2))
    {
        mexErrMsgTxt("DAC2BDC requires 1 or 2 output arguments.");
    }

    CheckMxReal(prhs[0],"A");
    CheckMxReal(prhs[1],"C");
    CheckMxReal(prhs[2],"u");
    CheckMxReal(prhs[3],"y");

    A = mxGetPr(prhs[0]);
    C = mxGetPr(prhs[1]);
    u = mxGetPr(prhs[2]);
    y = mxGetPr(prhs[3]);

    N = mxGetM(prhs[2]);
    m = mxGetN(prhs[2]);
    l = mxGetN(prhs[3]);
    n = mxGetM(prhs[0]);

    fD = 0;
    if (nlhs==2)
    {
        fD= 1;
    }

    if (mxGetN(prhs[0])!=n)
    {
        mexErrMsgTxt("A must be square.");
    }
    if (mxGetM(prhs[1])!=l)
    {
        mexErrMsgTxt("C has a wrong number of rows.");
    }
    if (mxGetN(prhs[1])!=n)
    {
        mexErrMsgTxt("C has a wrong number of columns.");
    }
    if (mxGetM(prhs[3])!=N)
    {
        mexErrMsgTxt("The input and output must have an equal number of samples.");
    }

    if (l<1)
    {
        mexErrMsgTxt("DAC2BDC requires an output");
    }
    if (m<1)
    {
        mexErrMsgTxt("DAC2BDC requires an input");
    }
    if (N<((n*m)+(fD?l*m+1:1)+n))
    {
        mexErrMsgTxt("The number of samples is too small.");
    }

    plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL);
    if(fD)plhs[1] = mxCreateDoubleMatrix(l,m,mxREAL);
    B = mxGetPr(plhs[0]);
    D = fD?mxGetPr(plhs[1]):NULL;

    if (n> 0)
    {
        CalculateBD(u,y,A,B,C,D,N,l,m,n,fD);
    }
    else
    {
        if(fD)
        {
            EstimateD(l,m,N,u,y,D);
        }
    }
}
