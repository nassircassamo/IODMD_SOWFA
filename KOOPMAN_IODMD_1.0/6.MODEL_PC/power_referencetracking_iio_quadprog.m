function [u]=power_referencetracking_iio_quadprog(mpcmodel,Hp,Hc,Inputs,Outputs,dirdmd,qq,rr)

dirdmdmpc=strcat(dirdmd,'/iioMPC');
if ~exist(dirdmdmpc,'dir') 
        mkdir(dirdmdmpc);
    end
%---------------------------------------------------------------------
% This file simulates the Model Predictive Control of a discrete-time
% MIMO system 
%---------------------------------------------------------------------

%% Take matrices from state spcae model to be used for MPC
A=mpcmodel.A;
B=mpcmodel.B;
C=mpcmodel.C;
D=mpcmodel.D;

%% Determine number of states, inputs and outputs
n = size(A,1);
q = size(B,2);
l = size(C,1);

%% Make extended matrices for increment input output (IIO) MPC formulation 
Ae=[eye(l) C; zeros(n,l) A];
Be=[zeros(l,q);B];
Ce=[eye(l) C];

%% Construct matrices to compute global optimum solution
% Construct a vector with the diagonal elements necessary to build other
% matrices
RRR = [];

for i= 1:Hp
  RRR = [RRR;Ce*(Ae^(i-1))*Be];
end

%==================
% Construct the H bar matrix
%==================
if q==1 & l==1
    Ru=zeros(Hp,Hp);
    for i = 1:Hp
        Ru=Ru+diag(RRR(i,1)*ones(Hp-i+1,1),-i+1);
    end
else
    Ru=[];
    for i = 1:Hp
        Ru = [Ru [zeros((i-1)*l,q);RRR(1:l*(Hp+1-i),:)]];
    end
end

% Reduce dimension of Ru, in case of Hc < Hp:
Ru = Ru(:,1:Hc*q);

%Fetch H bar
Hbar=Ru;

%==================
% Construct the Tau matrix
%==================
Tau=[];
for i = 1:Hp
  Tau = [Tau;Ce*(Ae^(i))];
end

%==================
% Construct the Imc matrix which contains Hc identity-matrices stacked
% on each other.
%==================
Imc=[];
for i = 1:Hc
  Imc = [Imc ; eye(q)];
end 


%==================
% Construct the Lmc matrix which contains Hc identity-matrices in
% the lower triangular part of the matrix
%==================

Lmc=[];
for i = 1:Hc
  Lmc = [Lmc [zeros((i-1)*q,q);Imc(1:q*(Hc+1-i),:)]];
end 
  

%% Generate time axis
ti=0;       		% starting time
tstep=2;    		% sample interval
tf = 2000+Hp*tstep;	   	% final instant
time=ti:tstep:tf;

%% Generate reference for tracking
 Pref_wt(:,1) = [0*ones(100,1); 0*ones(200,1);   0*ones(200,1); 0*ones(200,1);  0*ones(200+Hp+1,1)];
 Pref_wt(:,2) = [0*ones(100,1); 0.25*ones(200,1); 0.25*ones(200,1); 0.5*ones(200,1); 0.5*ones(200+Hp+1,1)];

%generate reference for total power in wind farm
for b=1:length(Pref_wt)
    Pref(b,1)=Pref_wt(b,1)+Pref_wt(b,2);
end

F=zeros(Hp, Hp*l);
for o=1:Hp
    for t=1:l
        alpha=o*2-1;
        F(o,alpha+t-1)=1;
    end
end

% Set initial states

xe(:,1) = [zeros(l,1); zeros(n,1)];  
u_k1=0;
%% =========================================================================
                          % SIMULATION LOOP
%  =========================================================================
t=0;
Q=qq*eye(Hp);
%Q(150,150)=qq*10;
R=rr*eye(Hc);

HH=2*(Hbar'*F'*Q*F*Hbar+R);

%Arest=-Imc;
%Arest=[-Imc;Imc];
Arest=[-Lmc;Lmc];

umin=0;
umax=4;

%brest = [umin*ones(q*Hc,1)];
%brest = [umin*ones(q*Hc,1)
        % umax*ones(q*Hc,1)];
        
 brest = [-umin*ones(q*Hc,1)+Imc*u_k1
         umax*ones(q*Hc,1)-Imc*u_k1];

cca=0;

for k=1:length(Pref)-Hp
    
            
    brest = [-umin*ones(q*Hc,1)+Imc*cca
         umax*ones(q*Hc,1)-Imc*cca];
     
    Ptilde=[];
    
    for v=1:Hp
        Ptilde=[Ptilde;Pref(k+v,1)];
    end
    
    % Calculate constant term to be used in MPC solution
    cterm{k}=F*Tau * xe(:,k) - Ptilde;
    cc=2*(Hbar'*F'*Q*cterm{k})';
    
%      if k==1
%      brest = [-umin*ones(q*Hc,1)+Imc*u_k1
%           umax*ones(q*Hc,1)-Imc*u_k1];
%      else
%          brest = [-umin*ones(q*Hc,1)+Imc*u(k-1)
%           umax*ones(q*Hc,1)-Imc*u(k-1)];
%          
%     end
     
    %------
  	% Calculate global optimum solution for non constrained optimization
  	% problem. Mathemtical solution is known
  	%------
    %deltau_opt=-pinv(Hbar'*F'*Q*F*Hbar+R)*Hbar'*F'*Q*cterm{k}; %global optimal solution to MPC problem under no constraints
    
    
    deltau_opt=quadprog(HH,cc,Arest,brest,[],[],[],[],[]);
    deltau(k)=deltau_opt(1); %receeding horizon
    cca=deltau(k)+cca;
    u(k)=cca;
    %------
  	% Simualte model and save information
  	%------ 
    xe(:,k+1)=Ae*xe(:,k)+Be*deltau(k);
    predictedpower(:,k)=Ce*xe(:,k);
end

%% Compute tracking error

%% Graphical visualisation
maxd=length(predictedpower);
 
%% Evaluate Turbine power tracking performance
figure600=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');
subplot(4,1,[1 2])
baseline=0*ones(1,maxd);
p5=plot(time(1:maxd),baseline,'LineWidth',1.6'); 
hold on
p5.Color=[0.25, 0.25, 0.25];
p5.LineStyle='--';

%static yaw deflection
p6=plot(time(1:maxd),baseline+0.6,'LineWidth',1.6'); 
p6.Color=[0.6350, 0.0780, 0.1840];
p6.LineStyle='--';

%reference tracking for dynamic yaw deflection
p7=plot(time(1:maxd),Pref(1:maxd),'LineWidth',1.6');
hold on
p7.Color='b';
p7.LineStyle='--';

%wind farm generator power output
p8=plot(time(1:maxd),predictedpower(1,1:maxd)'+predictedpower(2,1:maxd)','LineWidth',1.6');
hold on
p8.Color='g';

grid on
ylabel(' Power [MW] ')
title('Generator power tracking and comparison with different scenarios');
legend({'Baseline scenario - No wake control',...
    'Wake Control: static yaw deflection',...
    'Wake Control: dynamic yaw deflection reference'...
    'Model Predictive Contorl: total generated power'},...
    'Location','best','Orientation','horizontal');
legend('boxoff')
set(gca,'fontsize', 12)

subplot(4,1,3)
% p1=plot(time(1:maxd),Pref(1:maxd,1)','LineWidth',1.6'); %1
% p1.Color='b';
% hold on
p1m=plot(time(1:maxd),predictedpower(1,1:maxd)','LineWidth',1.6');
p1m.Color=[0.2 0.2 0.8];
ylabel(' Power [MW]')
title('Turbine 1 generator power tracking');
%legend({'Model generator power output'},'Location','bestoutside','Orientation','horizontal')
%legend('boxoff')
grid on
set(gca,'fontsize', 12)

subplot(4,1,4)
% p2=plot(time(1:maxd),Pref(1:maxd,2)','LineWidth',1.6'); %1
% p2.Color='b';
% hold on
p2m=plot(time(1:maxd),predictedpower(2,1:maxd)','LineWidth',1.6');
p2m.Color=[0.2 0.2 0.8];
ylabel(' Power [MW]')
title('Turbine 2 generator power tracking');
%legend({'Model generator power output'},'Location','bestoutside','Orientation','horizontal')
%legend('boxoff')
grid on
set(gca,'fontsize', 12)

export_fig(figure600,strcat(dirdmdmpc,'/image','TurbinesOutput'),'-nocrop','-m2'); 
shg
%% Evaluate control action through time
figure601=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');
subplot(2,1,1)
p3=plot(time(1:maxd),u(1:end),'LineWidth',1.6'); %1
p3.Color=[0.8500, 0.3250, 0.0980];
ylabel('\gamma [º]')
xlabel('Sampling instants k')
title('Optimal control action for generator power tracking');
grid on
set(gca,'fontsize', 12)

subplot(2,1,2)
deltau=diff(u);
p4=plot(time(2:maxd),deltau,'LineWidth',1.6'); %1
p4.Color=[0.8500, 0.3250, 0.0980];
ylabel('\Delta \gamma [º s^{-1}]');
xlabel('Sampling instants k')
title('Optimal control action  rate for generator power tracking');
grid on
set(gca,'fontsize', 12)
shg
export_fig(figure601,strcat(dirdmdmpc,'/image','ControlAction'),'-nocrop','-m2'); 


end


