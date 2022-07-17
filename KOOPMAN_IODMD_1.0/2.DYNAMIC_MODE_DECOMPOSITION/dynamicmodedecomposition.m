function [sys_red,FITje,U,S,V,method,X,X_p,Xd,dirdmd,x]=dynamicmodedecomposition(states, Inputs, Outputs, Deterministic ,method,r,maindir,f,dt)

%Def: This function aims to build a reduced order model from the states,
% input/output information and deterministic states gathered in the
% simulation and resampled

%Input arguments:
    %states: full snapshot matrix to be further splitted into its snapshot
    %time shifted version
    %Inputs: Input data matrix [U1 U2 U3 ... Un];
    %Outputs: Output data matrix [Y1 Y3 .... Yn];
    %Deterministic states, if any [Xd1 Xd2 Xd3 Xd4 ... Xdn];
    %method: algortihm to be used to build ROM
        %0: original Dynamic Mode Decomposition algortihm, where the A
        %state space matrix is computed (there is no contorl action)
        %1: Dyanmic Mode Decomposition with Control 
        %2. Input Output Dynamic Mode Decomposition
        %3. Extended Input Output Dynamic Mode DECOMPOSITION
        %4. Initial approach with systems identification
    %r: (maximum) truncation level (number of singular values to retain)
    
%Output arguments:
    % sys_red: the state space systems, according to the number of singular
    % values used
    % FITje: the Fit of the model for the given data
    % U, S, V: the matrices resulting from SVD to be used for DMD states reconstruction
    %method: method is later used, as reconstruction depends on the method
    %used, even though methodology is similar
    % X, X_p: DMD matrices to be used later for reconstruction 

%log:
    %0. first commit October 2020
    %1. function revised and comments added
    

%% DMD - Dynamic Mode Decomposition

% x(k+1) ~ A x(k)
% X' ~ AX

% X = [   |  |        |  ]
%     [   x1 x2 ... xm-1 ]
%     [   |  |        |  ]

% X'= [   |  |        |  ]
%     [   x2 x3 ...  xm  ]
%     [   |  |        |  ]

%define necessary matrices for DMD
X      = states(:,1:end-1);
X_p    = states(:,2:end);

out      = Outputs(:,1:end-1);
out_p    = Outputs(:,2:end);

inp      =Inputs(:,1:end-1);
inp_p    =Inputs(:,2:end);

Xd     =Deterministic(:,1:end-1);
Xd_p   =Deterministic(:,2:end);

%% (0) DMD - there is no control action
if method==0 %dmd algortihm
    
     dirdmd='DMDresults_DMD';
     dirdmd=strcat(maindir,dirdmd);
    if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
    end
    
    [U, S, V]=svds(X,r);

    for si=1:1:r
        
        part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si,r);
        
        Util=U(:,1:si);
        Sigtil=S(1:si,1:si);
        Vtil=V(:,1:si);
        
        approxA{si} = Util'*(X_p)*Vtil*inv(Sigtil);
        approxB{si} = zeros(si, 1);
        approxC{si}=zeros(2,si);
        approxD{si}=zeros(2,1);
        sys_red{si}=ss(approxA{si},approxB{si},approxC{si},approxD{si},2);
       
        FITje=0;
    end

%%  (1) DMD - there is a extrnal forcing term
elseif method==1 %dmdc algortihm

     
     dirdmd='DMDresults_DMDc';
     dirdmd=strcat(maindir,dirdmd);
    if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
    end
    
    %Goal of DMDc is to characterize the relationship between three
    %measurments: current measurment x(k), the future state x(k+1) and the
    %current control u(k). Relationship is approximated by the canonical
    %discrete linear dynamical system. 
    %Assumptions: all states are observable
    %are C will be identity.

    % x(k+1) ~  Ax(k) + Bu(k)
    % y(k)   ~  Cx(k) + Du(k)

    % X' ~ AX + BU 
    % Y  ~ CX + DU

    % X' ~ AX + BU ~ [A B][ X ] = G?
    %                     [ U ]

    %2 Perform SVD on the augmented data matrix SVD(?)=USV such that
    %G=X'VS^(-1)U*

    %3. Separate A and B by splitting left singular vectors into 2 separate
    %components

    % [A B] ~ [X'VS^(-1)U*1, X'VS^(-1)U*2 ]

    % U*1:
    % U*2:

    %4. Construct reduced order subspace from the output measurments X'. U
    %cannot be used to find low rank model of the dynamics and input matrixes
    %since it is defined for the input space which now includes both the state
    %measurments and the exogeneous inputs

    %SVD(X')=�S^V*

    % �  = �*A� = �*X'VS^(-1)U*1 �
    % B~ = �*B  = �*X'VS^(-1)U*2
    
    Omega = [X;inp];
    
    [U, Sig, V]=svds(Omega,r);
    [Uf, Sigf, Vf]=svds(X_p,r);
    FITje=zeros(2, r);
    OMEGA={};
    DAMPING={};
    x=cell(r,1);
    
    for si=1:1:r
        
        part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si,r);
        
        Util=U(:,1:si);
        Sigtil=Sig(1:si,1:si);
        Vtil=V(:,1:si);
        
        Uhat=Uf(:,1:si);
       % Sighat=Sigf(1:si,1:si);
       % Vbar=Vf(:,1:si);
        
        n=size(X,1);
        q=size(inp,1);
        U_1=Util(1:n,:);
        U_2=Util(n+1:n+q,:);
        
        approxA{si} = Uhat'*(X_p)*Vtil*inv(Sigtil)*U_1'*Uhat;
        approxB{si} = Uhat'*(X_p)*Vtil*inv(Sigtil)*U_2';
        approxC{si}=eye(si,si);
        approxD{si}=zeros(si,1);
        sys_red{si}=ss(approxA{si},approxB{si},approxC{si},approxD{si},2);
        
        
        close all
    end
    [xo]=dinit(sys_red{r}.A,sys_red{r}.B,sys_red{r}.C,sys_red{r}.D,[Inputs]',[Uhat'*states]'); 
    [ysim, t, xout]=lsim(sys_red{r}, [Inputs]',[],xo);   
    x=ysim;
    Xd={};    
    FITje=0;
    
           
  %% (2) ioDMD: Input Output Dynamic Mode Decomposition 
  
elseif method==2 %ioDMD
    
    %The goal of ioDMD is to capitalize on DMDc and extend it so a full
    %state space system may be obtained via usual subspace system
    %identification methods, estimating matrices A,B,C,D via lieast squares
    
    dirdmd='DMDresults_IODMD';
    dirdmd=strcat(maindir,dirdmd);
    if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
    end
    
    dirdmdident='DMDresults_IODMD/ident';
    dirdmdident=strcat(maindir,dirdmdident);
    if ~exist(dirdmdident,'dir') 
        mkdir(dirdmdident);
    end
    
    [U,S,V]=svds(X,r);
    FITje=zeros(2, r);
    FITje_fil=zeros(2, r);
    OMEGA={};
    DAMPING={};
    x=cell(r,1);

    for si=1:1:r
        
        part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si,r);
        
        Util=U(:,1:si);
        Sigtil=S(1:si,1:si);
        Vtil=V(:,1:si);
       
        all=[ Util'*X_p;out]*pinv([Sigtil*Vtil';inp]);
            
        A{si}=all(1:si,1:si);
        B{si}=all(1:si,si+1:end);
        C{si}=all(si+1:end, 1:si);
        D{si}=all(si+1:end, si+1:end);

        sys_red{si}=ss(A{si},B{si},C{si},D{si},2);
        
        %%REGULARISATION
%         kapa=[Sigtil*Vtil';inp]*[Sigtil*Vtil';inp]';
%         W=eye(size(sys_red{si}.A));
%         WW=[W, zeros(size(sys_red{si}.A,1),size(sys_red{si}.B,2));...
%              zeros(size(sys_red{si}.B,2),size(sys_red{si}.A,1)), ...
%              zeros(size(sys_red{si}.B,2) )];
% 
%         nrs=[sys_red{si}.A sys_red{si}.B;...
%             sys_red{si}.C sys_red{si}.D]*kapa*(kapa+1000000*WW)^(-1);
% 
%         An=nrs(1:si,1:si);
%         Bn=nrs(1:si,si+1:end);
%         Cn=nrs(si+1:end,1:si);
%         Dn=nrs(si+1:end,si+1:end);
%         nss=ss(An,Bn,Cn,Dn,2);
% 
%         sys_red{si}=nss;

        [FITje,OMEGA,DAMPING,fig1,x]=evaluatemodel(sys_red,si,Inputs,Outputs,FITje,OMEGA,DAMPING,'identification',x,states,U,Deterministic,method);
        warning off
        %export_fig(fig1,strcat(dirdmdident,'/image',num2str(10000+si)),'-nocrop','-m2')
        print2eps(strcat(dirdmdident,'/image',num2str(10000+si)),fig1)
        warning on
        close all
    end
        
    [fig200]=VAFpermodes(FITje,r,{});
    warning off
    export_fig(fig200,strcat(dirdmdident,'/image',num2str(1000+length(sys_red)+1)),'-nocrop','-m2')
    warning on
    Xd={};
    
%% (3) extioDMD: Extended Input Output Dynamic Mode Decomposition
    
elseif method==3
    
    % The goal of extioDMD is to obtain a state space system, as ioDMD
    % performs, but etending the existing state to others which are not
    % necessarily related with the previous. This cpaitalizes on the
    % convergence of DMD results and the Koopman Operator, where there is
    % significant evidence that by incuding non linear observables (which
    % may be functions of the pre exisitng states, or not) DMD provides
    % better results.
    %These observables used to extend the current state space are referred
    %to as determinisitc states, as they are measurable and known for the
    %current scenario
    
    dirdmd='DMDresults_EIODMD';
    dirdmd=strcat(maindir,dirdmd);
        if ~exist(dirdmd,'dir') 
            mkdir(dirdmd);
        end
        
    dirdmdident='DMDresults_EIODMD/ident';
    dirdmdident=strcat(maindir,dirdmdident);
        if ~exist(dirdmdident,'dir') 
            mkdir(dirdmdident);
        end    
    
    [U,S,V]=svds(X,r);
    FITje=zeros(2, r);
    FITje_fil=zeros(2, r);
    OMEGA={};
    DAMPING={};
    x=cell(r,1);
    
    for si=1:1:r
        part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si,r);
        
        Util=U(:,1:si);
        Sigtil=S(1:si,1:si);
        Vtil=V(:,1:si);
       
        all=[ Xd_p; Util'*X_p;out]*pinv([Xd; Sigtil*Vtil';inp]);
            
        A{si}=all(1:size(Xd,1)+si,1:size(Xd,1)+si);
        B{si}=all(1:size(Xd,1)+si,size(Xd,1)+si+1:end);
        C{si}=all(size(Xd,1)+si+1:end, 1:size(Xd,1)+si);
        D{si}=all(size(Xd,1)+si+1:end, size(Xd,1)+si+1:end);
        
        sys_red{si}=ss(A{si},B{si},C{si},D{si},2);
         
        %same as before   
        [FITje,OMEGA,DAMPING,fig1,x]=evaluatemodel(sys_red,si,Inputs,Outputs,FITje,OMEGA,DAMPING,'identification',x,states,U,Deterministic,method);
        warning off
        export_fig(fig1,strcat(dirdmdident,'/image',num2str(10000+si)),'-nocrop','-m2')
        warning on
        close all
    end
        
    [fig200]=VAFpermodes(FITje,r,{});
    warning off
    export_fig(fig200,strcat(dirdmdident,'/image',num2str(1000+length(sys_red)+1)),'-nocrop','-m2')
    warning on
           
    
elseif method==4
    %% Professor Wingerden Least Square Solution for state space problem

    dirdmd='DMDresults_Wing';
    dirdmd=strcat(maindir,dirdmd);
    if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
    end
    
    % Take singular value decomposition of X with rank r
    % X ~ USV*
    % U belongs to set C with size nxr
    % S belongs to set C with sixe rxr
    % V belongs to set C iwth size mxr, and * denotes the conjugate transpose 
    % r is the rank of the reduced SVD approximation to X

    % left singular vectors U are POD modes
    % Columns of U are orthonormal, so U*U=I and V*V=I
   

    %Projection of full state space onto POD modes 
    %Make use of SVD decompositoin in this phase

    % U*X' ~ U*A(USV*) + U*BU 
    % Y  ~ C(USV*) + DU

    %A new state X^= U*X=U*USV*=SV*
    % �  = U*A*U
    % B^ = U*B

    % U*X' ~ �SV* + B^U 

    % Matrixes � and B^can now be found via least squares 

    % || U*X' - [ � B^][ SV* ] ||
    %                  [  U  ]

    % [ � B^ ] = U*X [ SV* ]* x  [ SV*VS    SV*U_*  ] ^-1
    %                [ U_   ]    [ U_SV     U_U_*   ]
  
    %%%---%%%---%%%
    
    %including deterministic states, the matrix problem formulation will
    %be:
    
    % [ I 0  ] [ Xd ] ~ [ I 0  ] A [ Xd   ] + [ I 0  ] BU
    % [ 0 U* ] [ X' ] ~ [ 0 U* ]   [ USV* ]   [ 0 U* ]
    
    % with 
    %
    % � = [ I 0  ] A
    %     [ 0 U* ]
    
    % || [ I 0  ] [ Xd ]  -  �[ Xd   ] - B^ U    || 
    % || [ 0 U* ] [ X' ]      [ USV* ]           || 
    
    % [ � B^ ]
    
    %truncation/number of singular values for SVD decomposition
    
    Xd=[X1; X2; X3;X4];
    a=size(Xd);
    states=[states];
    [Uo,So,Vo]=svds(X,r);
    U=blkdiag(eye(size(Xd,1)),Uo); 
    QQ=[Xd; states];  
    X=QQ(:,1:end-1);
    X_p=QQ(:,2:end);
    S=blkdiag(eye(size(Xd,1)),So);
    V=[Xd(:,1:end-1)' Vo];
    FITje=zeros(2, r+a(1));
    OMEGA={};
    DAMPING={};
    for si=1:1:(r+a(1))
        part=2; subpart=5; [f]= MPC_progress(part,subpart,f,si,r);
        
        Attt{si}=U(:,1:si)'*QQ(:,2:end)*[S(1:si,1:si)*V(:,1:si)';...
            [U1(:,1:end-1)]]'*inv([S(1:si,1:si)*V(:,1:si)'*V(:,1:si)*S(1:si,1:si) S(1:si,1:si)*V(:,1:si)'* [U1(:,1:end-1)]'; [U1(:,1:end-1)]*V(:,1:si)*S(1:si,1:si) [U1(:,1:end-1)]*[U1(:,1:end-1)]' ] );
        
        At{si}=Attt{si}(:,1:si);
        Bt{si}=Attt{si}(:,si+1:end);
        Ctt{si}=[Y1(:,1:end-1);Y2(:,1:end-1)]*pinv(([S(1:si,1:si)*V(:,1:si)'; U1(:,1:end-1)]));
        
        Ct{si} = Ctt{si}(:,1:si);
        Dt{si} = Ctt{si}(:,si+1:end);
        
        sys_red{si}=ss(At{si},Bt{si},Ct{si}, Dt{si},2);
        
        [FITje,OMEGA,DAMPING,fig1]=evaluatemodel(sys_red,si,Inputs,Outputs,FITje,OMEGA,DAMPING,'identification');
        
        warning off
        export_fig(fig1,strcat(dirdmd,'/image',num2str(10000+si)),'-nocrop','-m2')
        warning on
        close all
    end
        [fig200]=VAFpermodes(FITje,r,{});
        warning off
        export_fig(fig200,strcat(dirdmd,'/image',num2str(10000+si+1)),'-nocrop','-m2')
        warning on
           
end    


    
end

