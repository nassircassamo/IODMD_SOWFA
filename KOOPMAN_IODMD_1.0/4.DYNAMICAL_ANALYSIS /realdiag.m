function [modalmodels,total,cmodalssd,xmodal,in]=realdiag(model,Inputs,Outputs,Diam,Uups,dt,dir)

%% Generate directory 
dirmodal='/ModalAnalysesVAF';
dir=strcat(dir,dirmodal);
    if ~exist(dir,'dir') 
        mkdir(dir);
    end

%% Initial conditions
%[xo]=dinit(model.A,model.B,model.C,model.D,[Inputs]',[Outputs]');

%% Convert to continuous time
sysc=d2c(model);

A=sysc.A;
B=sysc.B;
C=sysc.C;
D=sysc.D;

%% Find modal state space representation
%The diagonal blocks of Am are in increasing order of imaginary part of the
%iegenvalues of A
%Vm is the transformation matrix from (A,B,C) to (Am,Bm,Cm)
[V,D]=eig(A);
cev=abs(imag(diag(D)))/2/pi; %continuous time frequency in Herz
cevstr=cev*Diam/Uups; %continuous time frequency in Herz

[x ,in]=sort(diag(abs(D)));
D=D(in,in);
V=V(in,in);
B=B(in,:);
C=C(:,in);
[Vm, Am]=cdf2rdf(V,D);
Bm=inv(Vm)*B;
Cm=C*Vm;
%xo=xo(in);
cev=cev(in);
cevstr=cevstr(in);

%vector of eigenvalues order
veceigen=diag(D);

%vector of imaginary parts of eigenvlaues
eigenimag=imag(diag(D));

%% Generate continuos SS and respective discrete representation
cmodalss=ss(Am,Bm,Cm,sysc.D);
cmodalssd=c2d(cmodalss,dt);
[xo]=dinit(cmodalssd.A,cmodalssd.B,cmodalssd.C,cmodalssd.D,[Inputs]',[Outputs]');

%% Determine number of states, inputs and outputs
n = size(cmodalssd.A,1);
q = size(cmodalssd.B,2);
l = size(cmodalssd.C,1);

%% Fetch modal state space representatios
alpha=1;
mm=1;
modalmodels=cell(1,1);
% for eigenvalue=1:length(veceigen)
while alpha <= length(veceigen)
    
    if all(eigenimag(alpha))
        
        %Fetch correspondent mode state space matrices
        At{alpha}=cmodalssd.A(alpha:alpha+1,alpha:alpha+1);
        Bt{alpha}=cmodalssd.B(alpha:alpha+1,1:q);
        Ct{alpha}=cmodalssd.C(1:l,alpha:alpha+1);
        %Dt{alpha}=cmodalssd.D;
        Dt{alpha}=zeros(l,q);
        
        %make modal discrete state space system 
        modalmodels{mm}=ss(At{alpha},Bt{alpha},...
            Ct{alpha},Dt{alpha},2);
        modalic{mm}=xo(alpha:alpha+1);
        freq{mm}=cev(alpha);
        freqstr(mm)=cevstr(alpha);
        
        mm=mm+1;
        alpha=alpha+2;
        
    else
        At{alpha}=cmodalssd.A(alpha:alpha,alpha:alpha);
        Bt{alpha}=cmodalssd.B(alpha:alpha,1:q);
        Ct{alpha}=cmodalssd.C(1:l,alpha:alpha);
        %Dt{alpha}=cmodalssd.D;
        Dt{alpha}=zeros(l,q);
        
        modalmodels{mm}=ss(At{alpha},Bt{alpha},...
            Ct{alpha},Dt{alpha},2);
        modalic{mm}=xo(alpha:alpha);
        freq{mm}=cev(alpha:alpha);
        freqstr(mm)=cevstr(alpha:alpha);
        
        mm=mm+1;
        alpha=alpha+1;
        
    end    
end

modalmodels{mm}=ss(0,zeros(1,q),zeros(l,1),cmodalssd.D,dt);
modalic{mm}=0;

%% Calculate all outputs 
total=0;
for modal=1:length(modalmodels)
     %[xo]=dinit(modalmodels{modal}.A,modalmodels{modal}.B,modalmodels{modal}.C,modalmodels{modal}.D,[Inputs]',[Outputs]');
     [ysim, t, xout]=lsim(modalmodels{modal}, [Inputs]',[],modalic{modal});  
     total=ysim+total;
end

%% Calculate modal contribution: SEE EFFECT OF INCLUDING JUST MODAL MODEL

% %Establish baseline scenario
% [ysimtotal, t, xsimtotal]=lsim(cmodalssd, [Inputs]',[],xo);  
% xsimtotal=xsimtotal';
% FFF0=sqrt(sum(xsimtotal.^2));
% FF0=sqrt(sum(sum(xsimtotal.^2)));
% 
 alpha=1;
 mm=1;
 Att=zeros(n,n);
 Btt=zeros(n,q);
 Ctt=zeros(l,n);
 Dtt=zeros(q,l);
% 
 while alpha <= length(veceigen)
     if all(eigenimag(alpha))
         
         %Build state space where only modal entries are present
         Att(alpha:alpha+1,alpha:alpha+1)=modalmodels{mm}.A;
         Btt(alpha:alpha+1,1:q)=modalmodels{mm}.B;
         Ctt(1:l,alpha:alpha+1)=modalmodels{mm}.C;
         Dtt=zeros(l,q);
         %Dtt=modalmodels{mm}.D;
         
         %simualte modal with only entries of mode mm 
         bigssmodal{mm}=ss(Att,Btt,Ctt,Dtt,dt);
         [ysim, t, xout]=lsim(bigssmodal{mm}, [Inputs]',[],xo);  
         xmodal{mm}=xout';
         
          %advance
          mm=mm+1;
          alpha=alpha+2;
         
         %put back to zero 
         Att=zeros(n,n);
         Btt=zeros(n,q);
         Ctt=zeros(l,n);
         Dtt=zeros(q,l);
%         
     else
         Att(alpha:alpha,alpha:alpha)=modalmodels{mm}.A;
         Btt(alpha:alpha,1:q)=modalmodels{mm}.B;
         Ctt(1:l,alpha:alpha)=modalmodels{mm}.C;
         Dtt=zeros(l,q);
         %Dtt=modalmodels{mm}.D;
          
         %simualte modal with only entries of mode mm 
         bigssmodal{mm}=ss(Att,Btt,Ctt,Dtt,dt);
         [ysim, t, xout]=lsim(bigssmodal{mm}, [Inputs]',[],xo);  
         xmodal{mm}=xout';

         %make graph showing how the modes explain output and states
         %modalanalysisnew(ysimtotal,ysim,FFF0,FFF,FITmode,mm,dir);
%         
         %advance
         mm=mm+1;
         alpha=alpha+1;
%         
%         %put back to zero 
         Att=zeros(n,n);
         Btt=zeros(n,q);
         Ctt=zeros(l,n);
         Dtt=zeros(q,l);
%         
     end
%      
 end


%% Calculate modal contribution: SEE EFFECT OF INCLUDING JUST MODAL MODEL
%Establish baseline scenario
[ysimtotal, t, xsimtotal]=lsim(cmodalssd, [Inputs]',[],xo);  
xsimtotal=xsimtotal';
FFF0=sqrt(sum(xsimtotal.^2));
FF0=sqrt(sum(sum(xsimtotal.^2)));

alpha=1;
mm=1;
Att=cmodalssd.A;
Btt=cmodalssd.B;
Ctt=cmodalssd.C;
%
while alpha <= length(veceigen)
    if all(eigenimag(alpha))
        
        %Build state space where only modal entries are present
        Att(alpha:alpha+1,alpha:alpha+1)=zeros(2,2);
        Btt(alpha:alpha+1,1:q)=zeros(2,q);
        Ctt(1:l,alpha:alpha+1)=zeros(l,2);
        Dtt=zeros(l,q);
        
        %simualte modal with only entries of mode mm 
        bigssmodal{mm}=ss(Att,Btt,Ctt,Dtt,dt);
        [ysim, t, xout]=lsim(bigssmodal{mm}, [Inputs]',[],xo);  
        xout=xout';
        yD=cmodalssd.D*Inputs;
        yD=yD';
        ysim=yD+ysim;
        %save modal contribution
        xmode{mm}=xout;
        FFF{mm}=sqrt(sum(xout.^2));
        FF{mm}=sqrt(sum(sum(xout.^2)));
        FFFr(mm)=100*(FF{mm})./(FF0);
        
        %VAF of mode
        FITmode(1:2,mm)=vaf(ysimtotal,ysim);
        FITmode(3,mm)=vaf(FFF0,FFF{mm});
        
        %make graph showing how the modes explain output and states
        modalanalysisnew(ysimtotal,ysim,FFF0,FFF,FITmode,mm,dir,freqstr);
        
        %advance
        mm=mm+1;
        alpha=alpha+2;
        
        %put back to zero 
        Att=cmodalssd.A;
        Btt=cmodalssd.B;
        Ctt=cmodalssd.C;
        
    else
        Att(alpha:alpha,alpha:alpha)=zeros(1,1);
        Btt(alpha:alpha,1:q)=zeros(1,q);
        Ctt(1:l,alpha:alpha)=zeros(l,1);
        Dtt=modalmodels{mm}.D;
        
        %simualte modal with only entries of mode mm 
        bigssmodal{mm}=ss(Att,Btt,Ctt,Dtt,dt);
        [ysim, t, xout]=lsim(bigssmodal{mm}, [Inputs]',[],xo);  
        xout=xout';
        yD=cmodalssd.D*Inputs;
        yD=yD';
        ysim=yD+ysim;
        
        %save modal contribution
        xmode{mm}=xout;
        FFF{mm}=sqrt(sum(xout.^2));
        FF{mm}=sqrt(sum(sum(xout.^2)));
        FFFr(mm)=100*(FF{mm})./(FF0);

        %VAF of mode
        FITmode(1:2,mm)=vaf(ysimtotal,ysim);
        FITmode(3,mm)=vaf(FFF0,FFF{mm});
        
        %make graph showing how the modes explain output and states
        modalanalysisnew(ysimtotal,ysim,FFF0,FFF,FITmode,mm,dir,freqstr);
        
        %advance
        mm=mm+1;
        alpha=alpha+1;
        
        %put back to original
        Att=cmodalssd.A;
        Btt=cmodalssd.B;
        Ctt=cmodalssd.C;
        
    end
end

%analyse matrix D influence
 teste=ss(Att,Btt,Ctt,zeros(l,q),2);
 [yteste, t, xoutD]=lsim(teste, [Inputs]',[],xo);
 xoutD=xoutD';
 VAFD=vaf(ysimtotal,yteste);
 
 FFF{mm}=sqrt(sum(xoutD.^2));
 
 VAFD(3)=vaf(FFF0,FFF{mm});      
  
 VAFD=100-VAFD;
 

%% Make final graph summarising everything
FITmode2=100-FITmode;

figure999=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');
subplot(3,10,1:6)
yaxis([0 100])
s=stem(freqstr,FITmode2(1,:),'k','filled');
grid on
xlabel('Adimensionalised frequency St=f*D/U')
ylabel('1-VAF_{mode} [%]')
title('Modal contribution to turbine 1 generator power output at each frequency')
set(gca, 'FontSize', 12)
ax=gca;
ax.XTick=0:0.2:3;
subplot(3,10,7:9)
s=stem(freqstr,FITmode2(1,:),'k','filled');
grid on
xlabel('Adimensionalised frequency St=f*D/U')
title('Modal contribution at each frequency (detail)')
set(gca, 'FontSize', 12)
axis([0 0.4 0 16])
ax=gca;
ax.XTick=0:0.05:0.4;
subplot(3,10,10)
s=stem(1,VAFD(1),'k');
grid on
title('Feedthrough')
set(gca, 'FontSize', 12)

subplot(3,10,11:16)
s=stem(freqstr,FITmode2(2,:),'k','filled');
yaxis([0 100])
grid on
xlabel('Adimensionalised frequency St=f*D/U')
ylabel('1-VAF_{mode} [%]')
title('Modal contribution to turbine 2 generator power output at each frequency')
set(gca, 'FontSize', 12)
ax=gca;
ax.XTick=0:0.2:3;
subplot(3,10,17:19)
s=stem(freqstr,FITmode2(2,:),'k','filled');
grid on
xlabel('Adimensionalised frequency St=f*D/U')
title('Modal contribution at each frequency (detail)')
set(gca, 'FontSize', 12)
axis([0 0.4 0 60])
ax=gca;
ax.XTick=0:0.05:0.4;
subplot(3,10,20)
s=stem(1,VAFD(2),'k');
grid on
title('Feedthrough')
set(gca, 'FontSize', 12)
yaxis([0 60])

subplot(3,10,21:26)
s=stem(freqstr,FITmode2(3,:),'k','filled');
grid on
xlabel('Adimensionalised frequency St=f*D/U')
ylabel('1-VAF_{mode} [%]')
title('Modal contribution to flow evolution (states) at each frequency')
set(gca, 'FontSize', 12)
ax=gca;
ax.XTick=0:0.2:3;
subplot(3,10,27:29)
s=stem(freqstr,FITmode2(3,:),'k','filled');
grid on
xlabel('Adimensionalised frequency St=f*D/U')
title('Modal contribution at each frequency (detail)')
set(gca, 'FontSize', 12)
axis([0 0.4 0 20])
ax=gca;
ax.XTick=0:0.05:0.4;
subplot(3,10,30)
s=stem(1,VAFD(3),'k');
grid on
title('Feedthrough')
set(gca, 'FontSize', 12)
yaxis([0 60])

export_fig(figure999,strcat(dir,'/modalcontribution'),'-nocrop','-m2');

end

