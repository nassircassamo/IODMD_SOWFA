%% Example: Example from Reynders
close all; clear; clc;

%% LTI model

% LTI system matrices
h = 1; % Sample time
n = 4; % number of states
A22 = [0.5321 0.8349; -0.8349 0.5231];
A33 = [-0.9165 0.1313; -0.1313 -0.9615];
B2 = [0.0012; 0.0007]; 
C2 = [-0.0267 2.6723];
C3 = [-0.0792 7.9236];
D = -0.008;
K2 = [0.0060; 0.0009];
K3 = [0.0242; -0.1218];
A = [A22 K2*C3; zeros(2) A33];
B = [B2; zeros(2,1)];
K = [K2; K3];
C = [C2 C3];

% open-loop system
OL = ss(A,[B K],C,[D eye(1)],1);


%% Open-loop identification experiment
% Simulation of the model in open loop

% input signals
N = 10000; % number of samples
t = (0:N-1)';   % time samples
u = 1000*randn(N,1); % excitation signal

% noise
e = randn(N,1); % noise signal

% simulation
y0 = lsim(OL,[u zeros(N,1)],t);
y = lsim(OL,[u e],t);
disp('Signal to noise ratio (SNR) (open-loop)')
snr(y,y0)

% Defining a number of constants
p = 15;     % past window size
f = 15;     % future window size
[u,Du,y,Dy] = sigscale(u,y); % signal scaling

% PBSID-opt
[S,x] = dordvarx(u,y,f,p,'tikh','gcv');
figure, semilogy(S,'x');
title('Singular values')
x = dmodx(x,n);
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,u,y,f,p);
Dat = iddata(y',u',h);
Mi = abcdk2idss(Dat,Ai,Bi,Ci,Di,Ki);
Mi

% Variance-accounted-for (by Kalman filter)
yest = predict(Mi,Dat);
x0 = findstates(Mi,Dat);
disp('VAF of identified system')
vaf(y,yest.y)

% PBSID-opt (greybox)
[S,x,xd,xs] = dordvarx_grey(u,y,f,p,'tikh','gcv');
figure, semilogy(S,'x');
title('Singular values');
[x,R1,R2,R3] = dmodx_grey(x,xd,xs,n,2,4);
figure, plot(R1,'*');
title('Canonical values between full state and deterministic only');
figure, plot(R2,'*');
title('Canonical values between full state and stochastic only');
figure, plot(R3,'*');
title('Canonical values between deterministic and stochastic only');
[Ap,Bp,Cp,Dp,Kp] = dx2abcdk(x,u,y,f,p);
Dat = iddata(y',u',h);
Mp = abcdk2idss(Dat,Ap,Bp,Cp,Dp,Kp);
Mp

% Variance-accounted-for (by Kalman filter)
yest = predict(Mp,Dat);
x0 = findstates(Mp,Dat);
disp('VAF of identified system')
vaf(y,yest.y)


%% Identification results
%

% Plot eigenvalues
figure
hold on
title('poles of identified LTI systems')
[cx,cy] = pol2cart(linspace(0,2*pi),ones(1,100));
plot(cx,cy,'k');
plot(real(eig(OL.a)),imag(eig(OL.a)),'k+','LineWidth',0.1,'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',10);
plot(real(eig(Mi.a)),imag(eig(Mi.a)),'rx');
plot(real(eig(Mp.a)),imag(eig(Mp.a)),'gx');
axis([-1 1 -1 1]);
axis square
legend('STABBND','TRUE','PBSID-opt','GREY','Location','East');
hold off

% Bodediagram (open loop)
OLi = ss(Mi);
OLp = ss(Mp);
OLi = Dy*OLi*inv([Du 0; 0 1]);
OLp = Dy*OLp*inv([Du 0; 0 1]);
figure, bodemag(OL(1,1),'k',OLi(1,1),'r',OLp(1,1),'g');
figure, bodemag(OL(1,2),'k',OLi(1,2),'r',OLp(1,2),'g');


