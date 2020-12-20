/* FCORDOM */
#include <math.h>
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

void ForsytheRecursion(double*Br,double*Bi,double*wr,double*wi,
                       double*M,double*Z,double*wscale,mwSize s,mwSize N,mwSize l,mwSize m)
{
    mwSize scount;
    mwSize FreqCount,count,count2;

    double*PrevCol1,*PrevCol2,*CurrCol;

    double wmin,wmax,*wtemp;

    double*VecPointer,VecNorm;
    mwSize VecCount;

    wtemp= mxMalloc(2*N*sizeof(double));

    wmin= mxGetInf();
    wmax= 0;
    for(count= 0;count<N;count++)
    {
        *wscale= sqrt(wr[count]*wr[count]+wi[count]*wi[count]);
        if(*wscale<wmin)wmin= *wscale;
        if(*wscale> wmax)wmax= *wscale;
    }
    *wscale= (wmin+wmax)/2;
    for(count= 0;count<N;count++)
    {
        wtemp[2*count]= wr[count]/ *wscale;
        wtemp[2*count+1]= wi[count]/ *wscale;
    }

    for(scount= 0;scount<s;scount++)
    {

        switch(scount)
        {

        case 0:
            {

                for(count= 0;count<l;count++)
                {
                    for(FreqCount= 0;FreqCount<N;FreqCount++)
                    {
                        if((Br!=NULL)&&(Bi!=NULL))
                        {
                            for(count2= 0;count2<m;count2++)
                            {
                                M[2*(FreqCount*m+count2+count*N*m)]= Br[FreqCount*l*m+count2*l+count];
                                M[2*(FreqCount*m+count2+count*N*m)+1]= Bi[FreqCount*l*m+count2*l+count];
                            }
                        }
                        else
                        {
                            if(count<=min(l,m))
                            {
                                M[2*(FreqCount*m+count+count*N*m)]= 1;
                            }
                        }
                    }

                    VecPointer= M+2*count*N*m;

                    VecNorm= 0;
                    for(VecCount= 0;VecCount<2*N*m;VecCount++)
                    {
                        VecNorm+= VecPointer[VecCount]*VecPointer[VecCount];
                    }
                    VecNorm= (VecNorm==0)?0:1/sqrt(VecNorm);

                    for(VecCount= 0;VecCount<2*N*m;VecCount++)
                    {
                        VecPointer[VecCount]*= VecNorm;
                    }
                    VecNorm= (VecNorm==0)?0:1/VecNorm;
                    Z[s*(m+l)*count]= VecNorm;
                }
                break;
            }

        case 1:
            {

                PrevCol1= M;
                CurrCol= M+2*N*m*l;
                for(count= 0;count<l;count++)
                {
                    for(FreqCount= 0;FreqCount<N;FreqCount++)
                    {

                        for(count2= 0;count2<m;count2++)
                        {

                            CurrCol[2*(FreqCount*m+count2+count*N*m)]=
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)]*wtemp[2*FreqCount]-
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)+1]*wtemp[2*FreqCount+1];
                            CurrCol[2*(FreqCount*m+count2+count*N*m)+1]=
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)]*wtemp[2*FreqCount+1]+
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)+1]*wtemp[2*FreqCount];
                        }
                    }

                    VecPointer= CurrCol+2*count*N*m;

                    VecNorm= 0;
                    for(VecCount= 0;VecCount<2*N*m;VecCount++)
                    {
                        VecNorm+= VecPointer[VecCount]*VecPointer[VecCount];
                    }
                    VecNorm= (VecNorm==0)?0:1/sqrt(VecNorm);

                    for(VecCount= 0;VecCount<2*N*m;VecCount++)
                    {
                        VecPointer[VecCount]*= VecNorm;
                    }
                    VecNorm= (VecNorm==0)?0:1/VecNorm;

                    Z[scount+s*(m+l)*count]= VecNorm;
                }
                break;
            }

        default:
            {

                PrevCol2= M+2*N*m*l*(scount-2);
                PrevCol1= M+2*N*m*l*(scount-1);
                CurrCol= M+2*N*m*l*scount;
                for(count= 0;count<l;count++)
                {
                    for(FreqCount= 0;FreqCount<N;FreqCount++)
                    {
                        for(count2= 0;count2<m;count2++)
                        {

                            CurrCol[2*(FreqCount*m+count2+count*N*m)]=
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)]*wtemp[2*FreqCount]-
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)+1]*wtemp[2*FreqCount+1];
                            CurrCol[2*(FreqCount*m+count2+count*N*m)+1]=
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)]*wtemp[2*FreqCount+1]+
                                PrevCol1[2*(FreqCount*m+count2+count*N*m)+1]*wtemp[2*FreqCount];

                            CurrCol[2*(FreqCount*m+count2+count*N*m)]+=
                                Z[(scount-1)+count*(m+l)*s]*PrevCol2[2*(FreqCount*m+count2+count*N*m)];
                            CurrCol[2*(FreqCount*m+count2+count*N*m)+1]+=
                                Z[(scount-1)+count*(m+l)*s]*PrevCol2[2*(FreqCount*m+count2+count*N*m)+1];
                        }
                    }

                    VecPointer= CurrCol+2*count*N*m;
                    VecNorm= 0;
                    for(VecCount= 0;VecCount<2*N*m;VecCount++)
                    {
                        VecNorm+= VecPointer[VecCount]*VecPointer[VecCount];
                    }
                    VecNorm= (VecNorm==0)?0:1/sqrt(VecNorm);

                    for(VecCount= 0;VecCount<2*N*m;VecCount++)
                    {
                        VecPointer[VecCount]*= VecNorm;
                    }
                    VecNorm= (VecNorm==0)?0:1/VecNorm;
                    Z[scount+s*(m+l)*count]= VecNorm;
                }
                break;
            }

        }
    }

    mxFree(wtemp);
}

void CompressData(double*Hr,double*Hi,double*wr,double*wi,
                  double*R,double*Sout,mwSize N,mwSize l,mwSize m,mwSize s)
{
    mwSize LDR,LDU,count,count2,p;
    integer INFO,LWORK; 
    double wscale,*TAU,OPTWORK,*WORK,*S;
    double*R22;
    double c0,c1;
    double*Z,*Wtilde;
    double*WG;
    mwSize sml,Nm;

    WG= mxCalloc(2*s*(m+l)*N*m,sizeof(double));
    S= mxMalloc(s*l*sizeof(double));
    Z= R+s*(m+l)*(m+l)*s;
    R22= R+s*(m+l)*s*(m+l)+s*m;
    LDU= s*l;
    LDR= s*(m+l);

    Wtilde= mxCalloc(2*s*N,sizeof(double));
    ForsytheRecursion(NULL,NULL,wr,wi,Wtilde,Z,&wscale,s,N,1,1);
    for(count= 0;count<s;count++)
    {

        for(count2= 0;count2<N;count2++)
        {
            for(p= 0;p<m;p++)
            {
                WG[2*(m*count2+p+(m*count+p)*N*m)]= Wtilde[2*(count2+count*N)];
                WG[2*(m*count2+p+(m*count+p)*N*m)+1]= Wtilde[2*(count2+count*N)+1];
            }
        }
    }
    mxFree(Wtilde);

    ForsytheRecursion(Hr,Hi,wr,wi,WG+2*N*s*m*m,Z,&wscale,s,N,l,m);

    sml= s*(m+l);
    Nm= 2*N*m;
    c0= 0;
    c1= 1;
    dgemm("T","N",&sml,&sml,&Nm,&c1,WG,&Nm,WG,&Nm,&c0,R,&sml);
    dpotrf("L",&sml,R,&sml,&INFO);
    for(count= 0;count<s*(m+l);count++)
    {
        for(count2= 0;count2<count;count2++)
        {
            R[count2+count*s*(m+l)]= 0;
        }
    }
    if(INFO!=0)
    {
        mexWarnMsgTxt("Cholesky-factorization failed; falling back on QR-factorization.");

        Nm= 2*N*m;
        sml= s*(m+l);
        TAU= mxMalloc(sml*sizeof(double));
        LWORK= -1;
        dgeqrf(&Nm,&sml,WG,&Nm,TAU,&OPTWORK,&LWORK,&INFO);
        LWORK= ceil(OPTWORK);
        WORK= mxMalloc(LWORK*sizeof(double));
        dgeqrf(&Nm,&sml,WG,&Nm,TAU,WORK,&LWORK,&INFO);
        mxFree(WORK);
        mxFree(TAU);

        for(count= 0;count<s*(m+l);count++)
        {
            for(count2= count;count2<s*(m+l);count2++)
            {
                R[count2+count*s*(m+l)]= WG[count+count2*2*N*m];
            }
        }
    }

    for(count= 0;count<s*l;count++)
    {
        for(count2= count;count2<s*l;count2++)
        {
            R[((m+l)*s*s*(m+l)+m*s)+count2+count*s*(m+l)]= R[(m*s*s*(m+l)+m*s)+count2+count*s*(m+l)];
        }
    }

    LWORK= -1;
    dgesvd("O","N",&LDU,&LDU,R22,&LDR,S,NULL,&LDU,NULL,&LDU,&OPTWORK,&LWORK,&INFO);
    LWORK= ceil(OPTWORK);
    WORK= mxMalloc(LWORK*sizeof(double));
    dgesvd("O","N",&LDU,&LDU,R22,&LDR,S,NULL,&LDU,NULL,&LDU,WORK,&LWORK,&INFO);
    mxFree(WORK);

    for(count= 0;count<s;count++)Sout[count]= S[count];
    R[1*s*(m+l)]= m;
    R[2*s*(m+l)]= l;
    R[3*s*(m+l)]= s;
    R[2*s*(m+l)+1]= wscale;

    mxFree(S);
    mxFree(WG);
}

void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[])
{
    double*Hr,*Hi,*wr,*wi,*smat;
    double*R,*Sout,*Rold;
    mwSize N,l,m,s,*Hdims;
    mwSize HNoReal,HNoImag,WNoReal,WNoImag;
    mwSize count,FreqsOK;

    if((nrhs!=3)&&(nrhs!=4))
    {
        mexErrMsgTxt("FORDOMC requires 3 or 4 input arguments.");
    }
    if(nrhs==4)
    {
        mexErrMsgTxt("FORDOMC does not support data batching yet.");
    }

    CheckMxValid(prhs[0],"H");
    CheckMxValid(prhs[1],"w");
    CheckMxReal(prhs[2],"s");

    Hr= mxGetPr(prhs[0]);
    Hi= mxGetPi(prhs[0]);
    wr= mxGetPr(prhs[1]);
    wi= mxGetPi(prhs[1]);
    smat= mxGetPr(prhs[2]);

    if(smat==NULL)
    {
        mexErrMsgTxt("s must be specified.");
    }
    s= smat[0];

    if(mxGetNumberOfDimensions(prhs[0])!=3)
    {
        mexErrMsgTxt("H must be a 3D array.");
    }
    Hdims= (integer*)mxGetDimensions(prhs[0]);
    l= Hdims[0];
    m= Hdims[1];
    N= Hdims[2];

    Rold= NULL;
    if(nrhs==4)
    {
        CheckMxReal(prhs[3],"Previous R-factor");
        Rold= mxGetPr(prhs[3]);
    }


    HNoReal= 0;HNoImag= 0;
    WNoReal= 0;WNoImag= 0;
    if(Hr==NULL)HNoReal= 1;
    if(Hi==NULL)HNoImag= 1;
    if(wr==NULL)WNoReal= 1;
    if(wi==NULL)WNoImag= 1;

    if(mxGetM(prhs[1])!=N)
    {
        mexErrMsgTxt("w must contain the same number of frequencies as H.");
    }
    if(mxGetN(prhs[1])!=1)
    {
        mexErrMsgTxt("w must have one column.");
    }
    if(s<2)
    {
        mexErrMsgTxt("s must be at least 2.");
    }
    if(l<1)
    {
        mexErrMsgTxt("FORDOMC requires an output.");
    }
    if(m<1)
    {
        mexErrMsgTxt("FORDOMC requires an input.");
    }
    if(2*N*m<l*s)
    {
        mexErrMsgTxt("The number of samples is too small.");
    }
    if(nrhs==4)
    {
        if((mxGetM(prhs[3])!=s*(m+l))|(mxGetN(prhs[3])!=s*(m+2*l)))
        {
            mexErrMsgTxt("The size of the old R-factor is incorrect.");
        }
        if(Rold[s*(m+l)]!=m)
        {
            mexErrMsgTxt("The number of inputs does not correspond to the old R-factor.");
        }
        if(Rold[2*s*(m+l)]!=l)
        {
            mexErrMsgTxt("The number of outputs does not correspond to the old R-factor.");
        }
        if(Rold[3*s*(m+l)]!=s)
        {
            mexErrMsgTxt("The blocksize does not correspond to the old R-factor.");
        }
        if(Rold[2*s*(m+l)+1]!=0)
        {
            mexErrMsgTxt("The previous R-factor seems to belong to a continuous-time model.");
        }
    }

    FreqsOK= 1;
    if(!WNoReal)
    {
        if(!WNoImag)
        {
            for(count= 0;count<N;count++)
            {
                if(fabs(wr[count])> (10*mxGetEps()))FreqsOK= 0;
            }
        }
    }
    if(!FreqsOK)
    {
        mexErrMsgTxt("At least one complex frequency is not imaginary.\n");
    }

    plhs[0]= mxCreateDoubleMatrix(s,1,mxREAL);
    Sout= mxGetPr(plhs[0]);
    plhs[1]= mxCreateDoubleMatrix(s*(m+l),s*(m+2*l),mxREAL);
    R= mxGetPr(plhs[1]);

    if(HNoReal)Hr= mxCalloc(N*l*m,sizeof(double));
    if(HNoImag)Hi= mxCalloc(N*l*m,sizeof(double));
    if(WNoReal)wr= mxCalloc(N,sizeof(double));
    if(WNoImag)wi= mxCalloc(N,sizeof(double));

    CompressData(Hr,Hi,wr,wi,R,Sout,N,l,m,s);

    if(HNoReal)mxFree(Hr);
    if(HNoImag)mxFree(Hi);
    if(WNoReal)mxFree(wr);
    if(WNoImag)mxFree(wi);
}
