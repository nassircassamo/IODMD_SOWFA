/* LTIFRF */
#include <mex.h>
#include "f2c.h"
#include "matrix.h"
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

void TransformToUpperHessenberg(double*A,double*dA,double*B,double*C,
                                double*AArray,double*dAArray,double*BArray,double*CArray,
                                mwSize n,mwSize l,mwSize m)
{
    double *TauArray, *WorkArray;
    mwSize i;
    integer c1,info;
    mwSize LDWORK;

    for (i= 0;i<(n*n);i++) {
		AArray[i]= A[i];
	}
    for (i= 0;i<(n*m);i++) {
		BArray[i]= B[i];
	}
    for (i= 0;i<(l*n);i++) {
		CArray[i]= C[i];
	}
    if((dA!=NULL)&&(dAArray!=NULL)) {
        for(i= 0;i<(n*n);i++) {
			dAArray[i]= dA[i];
		}
	}

    if (n==1) {
		return;
	}

    TauArray = (double*)mxMalloc(n*sizeof(double));
    LDWORK = max(n,max(l,m));
    WorkArray = (double*)mxMalloc(LDWORK*sizeof(double));
    c1 = 1;

    dgehrd(&n,&c1,&n,AArray,&n,TauArray,WorkArray,&n,&info);

    dormhr("L","T",&n,&m,&c1,&n,AArray,&n,TauArray,BArray,&n,WorkArray,&LDWORK,&info);


    dormhr("R","N",&l,&n,&c1,&n,AArray,&n,TauArray,CArray,&l,WorkArray,&LDWORK,&info);


    if((dA!=NULL)&&(dAArray!=NULL))
    {
        dormhr("L","T",&n,&n,&c1,&n,AArray,&n,TauArray,dAArray,&n,WorkArray,&LDWORK,&info);

        dormhr("R","N",&n,&n,&c1,&n,AArray,&n,TauArray,dAArray,&n,WorkArray,&LDWORK,&info);
    }
    mxFree(TauArray);
    mxFree(WorkArray);
}

void ltifrN(double*A,double*B,double*C,double*D,
            double*wr,double*wi,double*Hr,double*Hi,
            mwSize n,mwSize l,mwSize m,mwSize N)
{
    mwSize i,j;
    integer info,c1;
    double *AArray,*BArray,*CArray;
    doublecomplex *GArray;

    doublecomplex *freq,*hinvb;
    double rcond,*evre,*evim;
    integer *iwork;
    double *dwork;
    doublecomplex *zwork;
    mwSize ldwork,lzwork;

    if ((n>=(m-1))&&(n>=(l-1)))
    {
        ldwork= 2*n;
    }
    else
    {
        ldwork = (m> l)?n+m-1:n+l-1;
    }
    lzwork = n*n+2*n;
    c1 = 1;
    freq = (doublecomplex*)mxMalloc(sizeof(freq));
    AArray = (double*)mxMalloc(n*n*sizeof(double));
    BArray = (double*)mxMalloc(n*m*sizeof(double));
    CArray = (double*)mxMalloc(l*n*sizeof(double));
    GArray = (doublecomplex*)mxMalloc(l*m*sizeof(doublecomplex));
    evre = (double*)mxMalloc(n*sizeof(double));
    evim = (double*)mxMalloc(n*sizeof(double));
    hinvb = (doublecomplex*)mxMalloc(n*m*sizeof(doublecomplex));
    iwork = (integer*)mxMalloc(n*sizeof(integer));
    dwork = (double*)mxMalloc(ldwork*sizeof(double));
    zwork = (doublecomplex*)mxMalloc(lzwork*sizeof(doublecomplex));

    for (i= 0;i<(n*n);i++) {
		AArray[i]= A[i];
	}
    for (i= 0;i<(n*m);i++) {
		BArray[i]= B[i];
	}
    for (i= 0;i<(l*n);i++) {
		CArray[i]= C[i];
	}

    for(i= 0;i<N;i++)
    {
        (freq[0]).r = wr[i];
        (freq[0]).i = wi[i];

        if(i==0)
        {
            tb05ad("A","G",&n,&m,&l,freq,AArray,&n,BArray,&n,CArray,&l,
                    &rcond,GArray,&l,evre,evim,hinvb,&n,iwork,dwork,&ldwork,
                    zwork,&lzwork,&info);
        }
        else
        {
            tb05ad("N","H",&n,&m,&l,freq,AArray,&n,BArray,&n,CArray,&l,
                    &rcond,GArray,&l,evre,evim,hinvb,&n,iwork,dwork,&ldwork,
                    zwork,&lzwork,&info);
        }

        for(j= 0;j<l*m;j++)
        {
            Hr[j+i*(l*m)]= (GArray[j]).r;
            Hi[j+i*(l*m)]= (GArray[j]).i;
        }

        if(D!=NULL)
        {
            for(j= 0;j<l*m;j++)
            {
                Hr[j+i*(l*m)]+= D[j];
            }
        }
    }
    mxFree(freq);
    mxFree(zwork);
    mxFree(dwork);
    mxFree(iwork);
    mxFree(hinvb);
    mxFree(evim);
    mxFree(evre);
    mxFree(GArray);
    mxFree(CArray);
    mxFree(BArray);
    mxFree(AArray);
}

void ltifrD(double*D,double*Hr,double*Hi,
            mwSize l,mwSize m,mwSize N)
{
    mwSize i,j;

    for(i= 0;i<N;i++)
    {
        for(j= 0;j<l*m;j++)
        {
            Hr[j+i*(l*m)]= D[j];
            Hi[j+i*(l*m)]= 0;
        }
    }
}

void ltifrdA(double*A,double*B,double*C,double*dA,
             double*wr,double*wi,double*Hr,double*Hi,
             mwSize n,mwSize l,mwSize m,mwSize N)
{
    mwSize i,j,k,p;
    integer *ipiv,info,c1;
    doublecomplex *G, *GB, *GdA;
    double *AArray,*dAArray,*BArray,*CArray;
    double *FirstTerm;

    AArray = (double*)mxMalloc(n*n*sizeof(double));
    dAArray = (double*)mxMalloc(n*n*sizeof(double));
    BArray = (double*)mxMalloc(n*m*sizeof(double));
    CArray = (double*)mxMalloc(l*n*sizeof(double));
    ipiv = (integer*)mxMalloc((n+2)*sizeof(integer));
    G = (doublecomplex*)mxMalloc(n*n*sizeof(doublecomplex));
    GB = (doublecomplex*)mxMalloc(n*m*sizeof(doublecomplex));
    GdA = (doublecomplex*)mxMalloc(n*n*sizeof(doublecomplex));
    FirstTerm = (double*)mxMalloc(2*l*n*sizeof(double));
    c1 = 1;

    TransformToUpperHessenberg(A,dA,B,C,AArray,dAArray,BArray,CArray,n,l,m);

    for (i= 0;i<N;i++)
    {
        for(j= 0;j<n*n;j++)
        {
            (G[j]).r= -AArray[j];
            (G[j]).i= 0;
        }
        for(j= 0;j<n;j++)
        {
            (G[(j+j*n)]).r+= wr[i];
            (G[(j+j*n)]).i= wi[i];
        }

        for(j= 0;j<(n*n);j++)
        {
            (GdA[j]).r= dAArray[j];
            (GdA[j]).i= 0;
        }
        for(j= 0;j<(n*m);j++)
        {
            (GB[j]).r= BArray[j];
            (GB[j]).i= 0;
        }

        mb02sz(&n,G,&n,ipiv,&info);
        mb02rz("N",&n,&n,G,&n,ipiv,GdA,&n,&info);
        mb02rz("N",&n,&m,G,&n,ipiv,GB,&n,&info);

        for(j= 0;j<2*n*l;j++)
        {
            FirstTerm[j]= 0;
        }
        for(k= 0;k<n;k++)
        {
            for(p= 0;p<n;p++)
            {
                for(j= 0;j<l;j++)
                {
                    FirstTerm[2*(j+k*l)]+= CArray[j+p*l]*((GdA[(p+k*n)]).r);
                    FirstTerm[2*(j+k*l)+1]+= CArray[j+p*l]*((GdA[(p+k*n)]).i);
                }
            }
        }

        for(k= 0;k<m;k++)
        {
            for(p= 0;p<n;p++)
            {
                for(j= 0;j<l;j++)
                {
                    Hr[j+k*l+i*(l*m)]+= FirstTerm[2*(j+p*l)]*((GB[(p+k*n)]).r);
                    Hr[j+k*l+i*(l*m)]-= FirstTerm[2*(j+p*l)+1]*((GB[(p+k*n)]).i);
                    Hi[j+k*l+i*(l*m)]+= FirstTerm[2*(j+p*l)]*((GB[(p+k*n)]).i);
                    Hi[j+k*l+i*(l*m)]+= FirstTerm[2*(j+p*l)+1]*((GB[(p+k*n)]).r);
                }
            }
        }
    }

    mxFree(FirstTerm);
    mxFree(GdA);
    mxFree(GB);
    mxFree(G);
    mxFree(ipiv);
    mxFree(CArray);
    mxFree(BArray);
    mxFree(dAArray);
    mxFree(AArray);
}

void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[])
{
    double*A,*B,*C,*D,*dA,*wr,*wi,*Hr,*Hi,*ColFlag;
    mwSize n,N,l,m;
    mwSize WRealFlag;
    mwSize Hdims[3];

    if(nrhs!=7)
    {
        mexErrMsgTxt("LTIFRF requires 7 input arguments");
    }

    CheckMxReal(prhs[0],"A");
    CheckMxReal(prhs[1],"B");
    CheckMxReal(prhs[2],"C");
    CheckMxReal(prhs[3],"D");
    CheckMxReal(prhs[4],"dA");
    CheckMxValid(prhs[5],"w");
    CheckMxReal(prhs[6],"OutOpt");

    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    D = mxGetPr(prhs[3]);
    dA = mxGetPr(prhs[4]);
    wr = mxGetPr(prhs[5]);
    wi = mxGetPi(prhs[5]);
    ColFlag = mxGetPr(prhs[6]);
    WRealFlag = 0;

    n = mxGetM(prhs[0]);
    m = mxGetN(prhs[1]);
    l = mxGetM(prhs[2]);
    N = mxGetM(prhs[5]);

    if (mxGetN(prhs[0])!=n){
        mexErrMsgTxt("A matrix must be square.");
    }
    if (mxGetM(prhs[1])!=n){
        mexErrMsgTxt("B matrix has wrong number of rows.");
    }
    if (mxGetN(prhs[2])!=n){
        mexErrMsgTxt("C matrix has wrong number of columns.");
    }
    if (!mxIsEmpty(prhs[4])){
        if ((mxGetM(prhs[4])!=n)||(mxGetN(prhs[4])!=n)){
            mexErrMsgTxt("dA and A matrix must be the same size.");
        }
    }

    if (n> 0)
    {
        if (!mxIsEmpty(prhs[3])){

            if (mxGetM(prhs[3])!=l){
                mexErrMsgTxt("D matrix has wrong number of rows.");
            }
            if (mxGetN(prhs[3])!=m){
                mexErrMsgTxt("D matrix has wrong number of columns.");
            }
        }
        else
        {
            D = NULL;
        }
    }
    else
    {
        A = NULL;
        B = NULL;
        C = NULL;
        if(!mxIsEmpty(prhs[3])){
            l= mxGetM(prhs[3]);
            m= mxGetN(prhs[3]);
        }
    }

    if(wi==NULL)
    {
        wi= mxCalloc(N,sizeof(double));
        WRealFlag= 1;
    }

    Hdims[0]= l;Hdims[1]= m;Hdims[2]= N;
    if(!mxIsEmpty(prhs[6]))
    {
        if(ColFlag[0]==1)
            plhs[0] = mxCreateDoubleMatrix(l*m*N,1,mxCOMPLEX);
        else
            if(ColFlag[0]==2)
                plhs[0] = mxCreateDoubleMatrix(l,m*N,mxCOMPLEX);
            else
                plhs[0] = mxCreateNumericArray(3,Hdims,mxDOUBLE_CLASS,mxCOMPLEX);
    }
    else
    {
        plhs[0] = mxCreateNumericArray(3,Hdims,mxDOUBLE_CLASS,mxCOMPLEX);
    }
    Hr = mxGetPr(plhs[0]);
    Hi = mxGetPi(plhs[0]);

    if (dA==NULL){
        if (A!=NULL) {
            if((l!=0)&(m!=0)) {
                ltifrN(A,B,C,D,wr,wi,Hr,Hi,n,l,m,N);
            }
        }
        else
        {
            if (D!=NULL)
            {
                ltifrD(D,Hr,Hi,l,m,N);
            }
        }
    }
    else
    {
        if (D!=NULL){
            mexWarnMsgTxt("D is ignored for the two-component dA calculation");
        }
        if ((l!=0)&(m!=0))
        {
            ltifrdA(A,B,C,dA,wr,wi,Hr,Hi,n,l,m,N);
        }
    }

    if (WRealFlag)
    {
        mxFree(wi);
    }
}
