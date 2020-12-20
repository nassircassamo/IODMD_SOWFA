/* DORDPO */
#include <mex.h>
#include "f2c.h"
#include "matrix.h"
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

void CompressData(double*u,double*y,mwSize s,mwSize N,mwSize l,mwSize m,
                  double*R,double*S,double*Sout,mwSize jobtype)
{
    mwSize LDR,twomls,LDWORK;
    integer *IWORK,IWARN,INFO;
    double RCOND,TOL,*DWORK;
    double *Rright,*Rbackup,*Rswap,TEMP;
    mwSize i1;

    mwSize count,count2,countm1;

    mwSize fUsedQR;

    LDR= max(max(2*(m+l)*s,3*m*s),5);
    twomls= 2*(m+l)*s;
    Rright= R+2*(m+l)*s*LDR;
    RCOND= 0;
    TOL= 0;

    LDWORK= N*2*(m+l)*s+4*(m+l)*s;

    IWORK= mxMalloc((m+l)*sizeof(integer));
    DWORK= mxMalloc(LDWORK*sizeof(double));

    i1= 1;
    fUsedQR= 0;

    if(jobtype==0)
    {
        ib01md("M","C","F","N",&s,&m,&l,&N,u,&N,y,&N,
                R,&LDR,IWORK,DWORK,&LDWORK,&IWARN,&INFO);
        Rright= R+2*(m+l)*s*LDR;
        for(count= 0;count<2*(m+l)*s*LDR;count++)
        {
            Rright[count]= R[count];
        }
        Rswap= Rright;
        TEMP= Rswap[0];
        Rswap[0]= Rswap[m*s+(m*s)*LDR];
        Rswap[m*s+(m*s)*LDR]= TEMP;
        for(count= 2;count<=m*s;count++)
        {
            dswap(&count,Rswap+(count-1)*LDR,&i1,Rswap+m*s+(m*s+(count-1))*LDR,&i1);
            countm1= count-1;
            dswap(&countm1,Rswap+(m*s+(count-1))*LDR,&i1,Rswap+(count-1)+m*s*LDR,&LDR);
        }

        for(count= 2*m*s;count<2*(m+l)*s;count++)
        {
            countm1= m*s;
            dswap(&countm1,Rswap+count*LDR,&i1,Rswap+count*LDR+m*s,&i1);
        }
        dpotrf("U",&twomls,Rright,&LDR,&INFO);
        if(INFO!=0)
        {
            fUsedQR= 1;
            mexWarnMsgTxt("Cholesky failed: using QR for this and any subsequent batches");
            ib01md("M","Q","F","N",&s,&m,&l,&N,u,&N,y,&N,
                    R,&LDR,IWORK,DWORK,&LDWORK,&IWARN,&INFO);
            Rright= R+2*(m+l)*s*LDR;
            for(count= 0;count<2*(m+l)*s*LDR;count++)
            {
                Rright[count]= R[count];
            }
        }
    }
    else if(jobtype==1)
    {

        Rbackup= mxMalloc(LDR*2*(m+l)*s*sizeof(double));
        for(count= 0;count<2*(m+l)*s;count++)
        {
            for(count2= 0;count2<=count;count2++)
            {
                Rbackup[count*LDR+count2]= R[count*LDR+count2];
            }
        }

        ib01md("M","C","I","N",&s,&m,&l,&N,u,&N,y,&N,
                R,&LDR,IWORK,DWORK,&LDWORK,&IWARN,&INFO);
        Rright= R+2*(m+l)*s*LDR;
        for(count= 0;count<2*(m+l)*s*LDR;count++)
        {
            Rright[count]= R[count];
        }
        Rswap= Rright;
        TEMP= Rswap[0];
        Rswap[0]= Rswap[m*s+(m*s)*LDR];
        Rswap[m*s+(m*s)*LDR]= TEMP;
        for(count= 2;count<=m*s;count++)
        {
            dswap(&count,Rswap+(count-1)*LDR,&i1,Rswap+m*s+(m*s+(count-1))*LDR,&i1);
            countm1= count-1;
            dswap(&countm1,Rswap+(m*s+(count-1))*LDR,&i1,Rswap+(count-1)+m*s*LDR,&LDR);
        }

        for(count= 2*m*s;count<2*(m+l)*s;count++)
        {
            countm1= m*s;
            dswap(&countm1,Rswap+count*LDR,&i1,Rswap+count*LDR+m*s,&i1);
        }
        dpotrf("U",&twomls,Rright,&LDR,&INFO);

        if(INFO!=0)
        {
            for(count= 0;count<2*(m+l)*s;count++)
            {
                for(count2= 0;count2<=count;count2++)
                {
                    R[count*LDR+count2]= Rbackup[count*LDR+count2];
                }
            }
            mxFree(Rbackup);

            fUsedQR= 1;
            mexWarnMsgTxt("Cholesky failed: using QR for this and any subsequent batches");

            Rswap= R;
            TEMP= Rswap[0];
            Rswap[0]= Rswap[m*s+(m*s)*LDR];
            Rswap[m*s+(m*s)*LDR]= TEMP;
            for(count= 2;count<=m*s;count++)
            {
                dswap(&count,Rswap+(count-1)*LDR,&i1,Rswap+m*s+(m*s+(count-1))*LDR,&i1);
                countm1= count-1;
                dswap(&countm1,Rswap+(m*s+(count-1))*LDR,&i1,Rswap+(count-1)+m*s*LDR,&LDR);
            }

            for(count= 2*m*s;count<2*(m+l)*s;count++)
            {
                countm1= m*s;
                dswap(&countm1,Rswap+count*LDR,&i1,Rswap+count*LDR+m*s,&i1);
            }

            dpotrf("U",&twomls,R,&LDR,&INFO);
            if(INFO!=0)
            {
                mexErrMsgTxt("Previous R-factor is corrupt!");
            }
            ib01md("M","Q","I","N",&s,&m,&l,&N,u,&N,y,&N,
                    R,&LDR,IWORK,DWORK,&LDWORK,&IWARN,&INFO);

            Rright= R+2*(m+l)*s*LDR;
            for(count= 0;count<2*(m+l)*s*LDR;count++)
            {
                Rright[count]= R[count];
            }
        }
    }
    else if(jobtype==2)
    {

        fUsedQR= 1;
        ib01md("M","Q","I","N",&s,&m,&l,&N,u,&N,y,&N,
                R,&LDR,IWORK,DWORK,&LDWORK,&IWARN,&INFO);

        Rright= R+2*(m+l)*s*LDR;
        for(count= 0;count<2*(m+l)*s*LDR;count++)
        {
            Rright[count]= R[count];
        }
    }
    ib01nd("M","M",&s,&m,&l,Rright,&LDR,S,&RCOND,IWORK,DWORK,&LDWORK,&IWARN,&INFO);

    for(count= 0;count<s;count++)
    {
        Sout[count]= S[count];
    }
    R[LDR-1]= N;
    R[LDR-2]= (double)fUsedQR;
    R[2*LDR-1]= l;
    R[3*LDR-1]= m;
    R[4*LDR-1]= s;

    mxFree(DWORK);
    mxFree(IWORK);
}

void mexFunction(int nlhs,mxArray*plhs[],int nrhs,const mxArray*prhs[])
{
    double*u,*y,*smat;
    double*R,*Rold,*S,*Sout;
    mwSize N,l,m,s,LDR;
    mwSize count,count2,jobtype;

    if((nrhs!=3)&(nrhs!=4))
    {
        mexErrMsgTxt("DORDPO requires 3 or 4 input arguments");
    }

    CheckMxReal(prhs[0],"u");
    CheckMxReal(prhs[1],"y");
    CheckMxReal(prhs[2],"s");
    if(nrhs==4)CheckMxReal(prhs[3],"Rold");
    u= mxGetPr(prhs[0]);
    y= mxGetPr(prhs[1]);
    smat= mxGetPr(prhs[2]);
    Rold= NULL;
    if(nrhs==4)Rold= mxGetPr(prhs[3]);

    if(smat==NULL)
    {
        mexErrMsgTxt("s must be specified.");
    }

    s = (integer)smat[0];
    N = mxGetM(prhs[1]);
    m = mxGetN(prhs[0]);
    l = mxGetN(prhs[1]);

    if((mxGetM(prhs[0])!=N)&(m!=0))
    {
        mexErrMsgTxt("u and y must have an equal number of samples.");
    }
    if(!(s> 1))
    {
        mexErrMsgTxt("s must be at least 2.");
    }
    /*if(!(m> 0))
    {
        mexErrMsgTxt("DORDPO requires an input.");
    } */
    if(!(l> 0))
    {
        mexErrMsgTxt("DORDPO requires an output.");
    }
    if(N<=(2*(m+l+1)*s-1))
    {
        mexErrMsgTxt("The number of samples must be at least 2*(m+l+1)*s-1");
    }

    if(nrhs==3)
    {
        jobtype= 0;
    }
    else
    {

        LDR= max(max(2*(m+l)*s,3*m*s),5);
        if((mxGetM(prhs[3])!=LDR)|(mxGetN(prhs[3])!=max((4*(m+l)*s),9)))
        {
            mexErrMsgTxt("The size of the old R-factor is incorrect.");
        }
        if(Rold[2*LDR-1]!=l)
        {
            mexErrMsgTxt("The current number of outputs does not match that in the old R-factor.");
        }
        if(Rold[3*LDR-1]!=m)
        {
            mexErrMsgTxt("The current number of inputs does not match that in the old R-factor.");
        }
        if(Rold[4*LDR-1]!=s)
        {
            mexErrMsgTxt("The current value of s does not match that in the old R-factor.");
        }
        jobtype= (Rold[LDR-2]> 0)?2:1;
    }

    if(nlhs!=2)
    {
        mexErrMsgTxt("DORDPO requires 2 output arguments");
    }

    LDR= max(max(2*(m+l)*s,3*m*s),5);
    plhs[0]= mxCreateDoubleMatrix(s,1,mxREAL);
    plhs[1]= mxCreateDoubleMatrix(LDR,max(4*(m+l)*s,9),mxREAL);
    Sout= mxGetPr(plhs[0]);
    R= mxGetPr(plhs[1]);
    S= mxMalloc(s*l*sizeof(double));

    if(jobtype> 0)
    {

        for(count= 0;count<2*(m+l)*s;count++)
        {
            for(count2= 0;count2<=count;count2++)
            {
                R[count*LDR+count2]= Rold[count*LDR+count2];
            }
        }
    }
    CompressData(u,y,s,N,l,m,R,S,Sout,jobtype);
}
