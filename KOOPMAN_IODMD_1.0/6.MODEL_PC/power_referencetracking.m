function [u,predictedpower,Pref]=power_referencetracking(mpcmodel,Hp,Hc,Inputs,Outputs,scalingfactors,dirdmd,qq,rr)
tic
 dirdmdmpc=strcat(dirdmd,'/MPC');
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

%% Construct matrices to compute global optimum solution
% Construct a vector with the diagonal elements necessary to build other
% matrices
RRR = [];

for i= 1:Hp
  RRR = [RRR;C*(A^(i-1))*B];
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
Ht=Hbar(Hp*l-1:end,:);

%==================
% Construct the Tau matrix
%==================
Tau=[];
for i = 1:Hp
  Tau = [Tau;C*(A^(i))];
end
Taut=Tau(Hp*l-1:end,:);

%% Generate time axis
ti=0;       		% starting time
tstep=2;    		% sample interval
tf = 4000+Hp*tstep;	   	% final instant
time=ti:tstep:tf;

%% Generate reference for tracking
 Pref_wt(:,1) = [0*ones(100,1); 0*ones(200,1);   0*ones(200,1); 0*ones(200,1);  0*ones(200+Hp+1,1); 0*ones(200+Hp+1,1); 0*ones(200+Hp+1,1)];
 Pref_wt(:,2) = [0*ones(100,1); 0.25*ones(200,1); 0.5*ones(200,1); 1*ones(200,1); 1.5*ones(200+Hp+1,1); 2*ones(200+Hp+1,1); 2.5*ones(200+Hp+1,1)];
 
 %generate reference for total power in wind farm
 for b=1:length(Pref_wt)
     Pref(b,1)=Pref_wt(b,1)+Pref_wt(b,2);
 end

%generate total windfarm power reference directly
%Pref(:,1) = [0; -0.1*ones(300,1); 0.8*ones(400+Hp+1,1)]/scalingfactors(1);

F=zeros(Hp, Hp*l);
for o=1:Hp
    for t=1:l
        alpha=o*2-1;
        F(o,alpha+t-1)=1;
    end
end
Ft=[1 1];

%% Compute initial states
[xo]=dinit(A,B,C,D,[Inputs]',[Outputs]');

x(:,1)=xo;
%% =========================================================================
                          % SIMULATION LOOP
%  =========================================================================
t=0;
Q=qq*eye(Hp);
R=rr*eye(Hc);

Qp=5000000000;

for k=1:length(Pref)-Hp
    
    Ptilde=[];
    for v=1:Hp
        Ptilde=[Ptilde;Pref(k+v,1)];
    end
    
    % Calculate constant term to be used in MPC solution
    cterm{k}=F*Tau * x(:,k) - Ptilde;
    ctermt{k}=Ft*Taut * x(:,k) - Ptilde(end);
    %------
  	% Calculate global optimum solution for non constrained optimization
  	% problem. Mathemtical solution is known
  	%------
    %u_opt=-pinv(Hbar'*F'*Q*F*Hbar+R)*Hbar'*F'*Q*cterm{k}; %global optimal solution to MPC problem under no constraints
    u_opt=-pinv(Hbar'*F'*Q*F*Hbar+R+Ht'*Ft'*Qp*Ft*Ht)*(Hbar'*F'*Q*cterm{k}+Ht'*Ft'*Qp*ctermt{k});
    u(k)=u_opt(1); %receeding horizon
    
    %------
  	% Simualte model and save information
  	%------ 
    x(:,k+1)=A*x(:,k)+B*u(k);
    predictedpower(:,k)=C*x(:,k);
    
%     %scaling outputs
%     predictedpower(1,k)=predictedpower(1,k)*scalingfactors(1);
%     predictedpower(1,k)=predictedpower(1,k)*scalingfactors(2);
    
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


%% Evaluate cost function and 
figure
set(gcf,'color','w','Position', get(0, 'Screensize'));
rectangle('Position', [-1 -1 2 2], 'Curvature', 1,'EdgeColor', 'k', 'LineStyle', '--');
axis(1.2*[-1 1 -1 1])
axis square
xlabel('Real axis \Re')
ylabel('Imaginary axis \Im')
title('Eigenvalue \lambda visualisation on the complex plane')
set(gca, 'FontSize', 14)
grid on
grid minor
hold on
plot(real(eig(A)),imag(eig(A)),'b+'); 
hold on; 
plot(real(eig(A+B*ones(1,size(F*Hbar,1))*inv(Hbar'*F'*Q*F*Hbar+eye(size(F,1))*R)*Hbar'*F'*(eye(size(F,1))*Q)*F*Tau)),...
    imag(eig((A-B*ones(1,size(F*Hbar,1))*inv(Hbar'*F'*(eye(size(F,1))*Q)*F*Hbar+eye(size(F,1))*R)*Hbar'*F'*(eye(size(F,1))*Q)*F*Tau))),'r+')
toc

end


