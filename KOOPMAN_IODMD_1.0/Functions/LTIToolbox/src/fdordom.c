/* FDORDOM */
#include <math.h>
#include <time.h>
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

void FillVandermondeProduct(double*Hr,double*Hi,double*wr,double*wi,
                            double*X,mwSize N,mwSize l,mwSize m,mwSize s,double*KKT)
{
    double*wtemp;
    double*GwNC;
    double*Gw;
    double*GGw;
    double zrtemp,zitemp;
    double X11Value;
    double*X21Block1,*X21Block2,*X22Block;
    mwSize FreqCount,BlockCount,CopyCount,count,count2,p;
    mwSize Fml,Fll;

    wtemp= mxMalloc(2*N*sizeof(double));
    GwNC= mxMalloc(2*l*m*N*sizeof(double));
    Gw= mxMalloc(2*l*m*N*sizeof(double));
    GGw= mxMalloc(2*l*l*N*sizeof(double));
    X21Block1= mxMalloc(m*l*sizeof(double));
    X21Block2= mxMalloc(m*l*sizeof(double));
    X22Block= mxMalloc(l*l*sizeof(double));

    for(FreqCount= 0;FreqCount<N;FreqCount++)
    {

        wtemp[2*FreqCount]= 1;
        wtemp[2*FreqCount+1]= 0;

        for(count= 0;count<l;count++)
        {
            for(count2= 0;count2<m;count2++)
            {
                GwNC[2*(l*m*FreqCount+count2*l+count)]= Hr[l*m*FreqCount+count+count2*l];
                GwNC[2*(l*m*FreqCount+count2*l+count)+1]= Hi[l*m*FreqCount+count+count2*l];
                Gw[2*(l*m*FreqCount+count2*l+count)]= Hr[l*m*FreqCount+count+count2*l];
                Gw[2*(l*m*FreqCount+count2*l+count)+1]= -Hi[l*m*FreqCount+count+count2*l];
            }
        }

        for(count= 0;count<l;count++)
        {
            for(count2= 0;count2<l;count2++)
            {
                zrtemp= 0;zitemp= 0;
                for(p= 0;p<m;p++)
                {
                    zrtemp+= Hr[l*m*FreqCount+count+p*l]*Hr[l*m*FreqCount+count2+p*l]
                             +Hi[l*m*FreqCount+count+p*l]*Hi[l*m*FreqCount+count2+p*l];
                    zitemp+= -Hr[l*m*FreqCount+count+p*l]*Hi[l*m*FreqCount+count2+p*l]
                             +Hi[l*m*FreqCount+count+p*l]*Hr[l*m*FreqCount+count2+p*l];
                }
                GGw[2*(l*l*FreqCount+count+count2*l)]= zrtemp;
                GGw[2*(l*l*FreqCount+count+count2*l)+1]= zitemp;
            }
        }
    }

    for(BlockCount= 0;BlockCount<s;BlockCount++)
    {

        X11Value= 0;
        for(count= 0;count<l*m;count++)X21Block1[count]= 0;
        for(count= 0;count<l*m;count++)X21Block2[count]= 0;
        for(count= 0;count<l*l;count++)X22Block[count]= 0;
        for(FreqCount= 0;FreqCount<N;FreqCount++)
        {
            X11Value+= wtemp[2*FreqCount];
            for(count= 0;count<l*m;count++)
            {
                X21Block1[count]+= GwNC[2*(FreqCount*m*l+count)];
                X21Block2[count]+= Gw[2*(FreqCount*m*l+count)];
            }

            for(count= 0;count<l*l;count++)
            {
                X22Block[count]+= GGw[2*(FreqCount*l*l+count)];
            }
        }

        for(CopyCount= 0;CopyCount<(s-BlockCount);CopyCount++)
        {


            KKT[CopyCount+(CopyCount+BlockCount)*s]+= X11Value;
            for(count= 0;count<m;count++)
            {
                X[(CopyCount+BlockCount)*m+count+(CopyCount*m+count)*(m+l)*s]+= X11Value;
            }

            for(count= 0;count<m;count++)
            {
                for(count2= 0;count2<l;count2++)
                {
                    X[(CopyCount+BlockCount)*l+m*s+count2+(CopyCount*m+count)*(m+l)*s]
                    += X21Block1[count2+count*l];
                }
            }
            if(BlockCount> 0)
            {


                for(count= 0;count<m;count++)
                {
                    for(count2= 0;count2<l;count2++)
                    {
                        X[CopyCount*l+count2+m*s+((CopyCount+BlockCount)*m+count)*(m+l)*s]
                        += X21Block2[count*l+count2];
                    }
                }
            }

            for(count= 0;count<l;count++)
            {
                for(count2= (BlockCount==0)?count:0;count2<l;count2++)
                {
                    X[m*s+(CopyCount+BlockCount)*l+count2+(m*s+CopyCount*l+count)*(m+l)*s]
                    += X22Block[count2+count*l];
                }
            }
        }

        for(FreqCount= 0;FreqCount<N;FreqCount++)
        {
            Fml= FreqCount*m*l;
            Fll= FreqCount*l*l;

            zrtemp= wtemp[2*FreqCount]*wr[FreqCount]-wtemp[2*FreqCount+1]*wi[FreqCount];
            zitemp= wtemp[2*FreqCount]*wi[FreqCount]+wtemp[2*FreqCount+1]*wr[FreqCount];
            wtemp[2*FreqCount]= zrtemp;
            wtemp[2*FreqCount+1]= zitemp;

            for(count= 0;count<l*m;count++)
            {
                zrtemp= GwNC[2*(Fml+count)]*wr[FreqCount]
                        -GwNC[2*(Fml+count)+1]*wi[FreqCount];
                zitemp= GwNC[2*(Fml+count)]*wi[FreqCount]
                        +GwNC[2*(Fml+count)+1]*wr[FreqCount];
                GwNC[2*(Fml+count)]= zrtemp;
                GwNC[2*(Fml+count)+1]= zitemp;

                zrtemp= Gw[2*(Fml+count)]*wr[FreqCount]
                        -Gw[2*(Fml+count)+1]*wi[FreqCount];
                zitemp= Gw[2*(Fml+count)]*wi[FreqCount]
                        +Gw[2*(Fml+count)+1]*wr[FreqCount];
                Gw[2*(Fml+count)]= zrtemp;
                Gw[2*(Fml+count)+1]= zitemp;
            }

            for(count= 0;count<l*l;count++)
            {
                zrtemp= GGw[2*(Fll+count)]*wr[FreqCount]
                        -GGw[2*(Fll+count)+1]*wi[FreqCount];
                zitemp= GGw[2*(Fll+count)]*wi[FreqCount]
                        +GGw[2*(Fll+count)+1]*wr[FreqCount];
                GGw[2*(Fll+count)]= zrtemp;
                GGw[2*(Fll+count)+1]= zitemp;
            }
        }
    }

    mxFree(X22Block);
    mxFree(X21Block2);
    mxFree(X21Block1);
    mxFree(GwNC);
    mxFree(Gw);
    mxFree(GGw);
    mxFree(wtemp);
}

void ObtainNormalRFactor(double*Hr,double*Hi,double*wr,double*wi,
                         double*R,mwSize N,mwSize l,mwSize m,mwSize s,double*Rold)
{
    double*WG;
    mwSize BlockCol,Col,Row,Freq;
    integer LWORK,INFO,AM,AN;
    double*W,*G,*TAU,*WORK,OPTWORK;
    mwSize LDWG;

    if(Rold==NULL)
    {
        WG= mxCalloc(2*N*m*s*(m+l),sizeof(double));
        LDWG= 2*N*m;
        W= WG;
        G= W+LDWG*m*s;
    }
    else
    {
        WG= mxCalloc(2*(N*m+s*(m+l))*s*(m+l),sizeof(double));
        LDWG= 2*N*m+s*(m+l);
        W= WG+s*(m+l);
        G= W+LDWG*m*s;
        for(Col= 0;Col<s*(m+l);Col++)
        {
            for(Row= 0;Row<=Col;Row++)
            {
                WG[Row+Col*LDWG]= Rold[Col+Row*s*(m+l)];
            }
        }
    }

    for(BlockCol= 0;BlockCol<s;BlockCol++)
    {
        if(BlockCol==0)
        {

            for(Col= 0;Col<m;Col++)
            {
                for(Freq= 0;Freq<N;Freq++)
                {
                    W[2*(Freq*m+Col)+Col*LDWG]= 1;
                }
            }

        }
        else
        {

            for(Col= 0;Col<m;Col++)
            {
                for(Freq= 0;Freq<N;Freq++)
                {
                    if(Col==0)
                    {

                        W[m*(2*Freq+BlockCol*LDWG)]= W[m*(2*Freq+(BlockCol-1)*LDWG)]*wr[Freq]
                                                     -W[m*(2*Freq+(BlockCol-1)*LDWG)+1]*wi[Freq];
                        W[m*(2*Freq+BlockCol*LDWG)+1]= W[m*(2*Freq+(BlockCol-1)*LDWG)]*wi[Freq]
                                                       +W[m*(2*Freq+(BlockCol-1)*LDWG)+1]*wr[Freq];
                    }
                    else
                    {

                        W[(2*(Freq*m+Col)+(BlockCol*m+Col)*LDWG)]=
                            W[(2*(Freq*m+(Col-1))+(BlockCol*m+(Col-1))*LDWG)];
                        W[(2*(Freq*m+Col)+(BlockCol*m+Col)*LDWG)+1]=
                            W[(2*(Freq*m+(Col-1))+(BlockCol*m+(Col-1))*LDWG)+1];
                    }
                }
            }
        }
    }

    for(BlockCol= 0;BlockCol<s;BlockCol++)
    {
        if(BlockCol==0)
        {

            for(Freq= 0;Freq<N;Freq++)
            {
                for(Col= 0;Col<l;Col++)
                {
                    for(Row= 0;Row<m;Row++)
                    {
                        G[(2*(Freq*m+Row)+Col*LDWG)]= Hr[l*m*Freq+Row*l+Col];
                        G[(2*(Freq*m+Row)+Col*LDWG)+1]= Hi[l*m*Freq+Row*l+Col];
                    }
                }
            }
        }
        else
        {

            for(Col= 0;Col<l;Col++)
            {
                for(Freq= 0;Freq<N;Freq++)
                {
                    for(Row= 0;Row<m;Row++)
                    {

                        G[(2*(Row+Freq*m)+(Col+BlockCol*l)*LDWG)]=
                            G[(2*(Row+Freq*m)+(Col+(BlockCol-1)*l)*LDWG)]*wr[Freq]
                            -G[(2*(Row+Freq*m)+(Col+(BlockCol-1)*l)*LDWG)+1]*wi[Freq];
                        G[(2*(Row+Freq*m)+(Col+BlockCol*l)*LDWG)+1]=
                            G[(2*(Row+Freq*m)+(Col+(BlockCol-1)*l)*LDWG)]*wi[Freq]
                            +G[(2*(Row+Freq*m)+(Col+(BlockCol-1)*l)*LDWG)+1]*wr[Freq];
                    }
                }
            }

        }
    }

    AM= LDWG;
    AN= s*(m+l);
    TAU= mxMalloc(AN*sizeof(double));
    LWORK= -1;
    dgeqrf(&AM,&AN,WG,&AM,TAU,&OPTWORK,&LWORK,&INFO);
    LWORK= ceil(OPTWORK);
    WORK= mxMalloc(LWORK*sizeof(double));
    dgeqrf(&AM,&AN,WG,&AM,TAU,WORK,&LWORK,&INFO);

    mxFree(WORK);
    mxFree(TAU);

    for(Row= 0;Row<s*(m+l);Row++)
    {
        for(Col= 0;Col<=Row;Col++)
        {
            R[Row+Col*s*(m+l)]= WG[Col+Row*LDWG];
        }
    }

    mxFree(WG);
}

void CompressData(double*Hr,double*Hi,double*wr,double*wi,
                  double*R,double*Sout,mwSize N,mwSize l,mwSize m,mwSize s,double*Rold)
{
    mwSize count,count2,p;
    integer INFO,LWORK;
    mwSize LDR,LDU,LDK;
    double*S,*WORK,OPTWORK;
    double*R22;
    double*KKT,*K;
    double c1;

    LDU= l*s;
    LDR= s*(m+l);
    LDK= s;
    KKT= mxCalloc(LDK*LDK,sizeof(double));
    R22= R+(m+l)*s*s*(m+l)+m*s;
    c1= 1;

    if(Rold!=NULL)
    {
        for(count= 0;count<LDR;count++)
        {

            for(count2= count;count2<LDR;count2++)
            {
                R[count2+count*LDR]= Rold[count2+count*LDR];
            }
        }

        dtrmm("R","L","T","N",&LDR,&LDR,&c1,Rold,&LDR,R,&LDR);
        for(count= 0;count<LDR;count++)
        {

            for(count2= 0;count2<count;count2++)
            {
                R[count2+count*LDR]= 0;
            }
        }
    }

    for(count= 0;count<LDK;count++)
    {
        for(count2= count;count2<LDK;count2++)
        {
            KKT[count+count2*LDK]= R[m*count2+m*count*LDR];
        }
    }

    FillVandermondeProduct(Hr,Hi,wr,wi,R,N,l,m,s,KKT);

    dpotrf("U",&LDK,KKT,&LDK,&INFO);
    if(INFO!=0)
    {
        mexErrMsgTxt("FORDOM: Frequency Vandermonde matrix is non positive definite.");
    }

    LDR= s*(m+l);
    dpotrf("L",&LDR,R,&LDR,&INFO);

    if(INFO!=0)
    {
        mexWarnMsgTxt("Cholesky-factorization failed; falling back on QR-factorization.");
        ObtainNormalRFactor(Hr,Hi,wr,wi,R,N,l,m,s,Rold);
    }

    for(count= 0;count<s*l;count++)
    {
        for(count2= 0;count2<s*l;count2++)
        {
            R[((m+l)*s*s*(m+l)+m*s)+count2+count*s*(m+l)]= R[(m*s*s*(m+l)+m*s)+count2+count*s*(m+l)];
        }
    }

    K= mxCalloc(s*l*s*l,sizeof(double));
    for(count= 0;count<s;count++)
    {
        for(count2= count;count2<s;count2++)
        {
            for(p= 0;p<l;p++)
            {
                K[l*count+p+(l*count2+p)*s*l]= KKT[count+s*count2];
            }
        }
    }
    mxFree(KKT);

    dtrtrs("U","N","N",&LDU,&LDU,K,&LDU,R22,&LDR,&INFO);

    S= mxMalloc(s*l*sizeof(double));
    LWORK= -1;
    dgesvd("O","N",&LDU,&LDU,R22,&LDR,S,NULL,&LDU,NULL,&LDU,&OPTWORK,&LWORK,&INFO);
    LWORK= ceil(OPTWORK);
    WORK= mxMalloc(LWORK*sizeof(double));
    dgesvd("O","N",&LDU,&LDU,R22,&LDR,S,NULL,&LDU,NULL,&LDU,WORK,&LWORK,&INFO);
    mxFree(WORK);

    dtrmm("L","U","N","N",&LDU,&LDU,&c1,K,&LDU,R22,&LDR);

    for(count= 0;count<s;count++)Sout[count]= S[count];
    R[1*s*(m+l)]= m;
    R[2*s*(m+l)]= l;
    R[3*s*(m+l)]= s;

    mxFree(K);
    mxFree(S);
}

void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[])
{
    double*Hr,*Hi,*wr,*wi,*smat;
    double*R,*Sout,*Rold;
    mwSize N,l,m,s,*Hdims;
    mwSize HNoReal,HNoImag,WNoReal,WNoImag;
    mwSize count,FreqMagnitudeOK;

    if((nrhs!=3)&&(nrhs!=4))
    {
        mexErrMsgTxt("FORDOM requires 3 or 4 input arguments.");
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
        mexErrMsgTxt("FORDOM requires an output.");
    }
    if(m<1)
    {
        mexErrMsgTxt("FORDOM requires an input.");
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

    FreqMagnitudeOK= 1;
    for(count= 0;count<N;count++)
    {
        if(WNoReal)
        {
            if(fabs(wi[count]-1)> 10*mxGetEps())FreqMagnitudeOK= 0;
        }
        else
        {
            if(WNoImag)
            {
                if(fabs(wr[count]-1)> 10*mxGetEps())FreqMagnitudeOK= 0;
            }
            else
            {
                if(fabs(sqrt(wr[count]*wr[count]+wi[count]*wi[count])-1)> 10*mxGetEps())
                    FreqMagnitudeOK= 0;
            }
        }
    }
    if(!FreqMagnitudeOK)
    {
        mexErrMsgTxt("At least one complex frequency is not unit-magnitude.\n");
    }

    plhs[0]= mxCreateDoubleMatrix(s,1,mxREAL);
    Sout= mxGetPr(plhs[0]);
    plhs[1]= mxCreateDoubleMatrix(s*(m+l),s*(m+2*l),mxREAL);
    R= mxGetPr(plhs[1]);

    if(HNoReal)Hr= mxCalloc(N*l*m,sizeof(double));
    if(HNoImag)Hi= mxCalloc(N*l*m,sizeof(double));
    if(WNoReal)wr= mxCalloc(N,sizeof(double));
    if(WNoImag)wi= mxCalloc(N,sizeof(double));

    CompressData(Hr,Hi,wr,wi,R,Sout,N,l,m,s,Rold);

    if(HNoReal)mxFree(Hr);
    if(HNoImag)mxFree(Hi);
    if(WNoReal)mxFree(wr);
    if(WNoImag)mxFree(wi);
}

