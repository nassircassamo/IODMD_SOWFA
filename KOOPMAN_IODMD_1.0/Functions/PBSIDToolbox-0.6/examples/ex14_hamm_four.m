%% Example 14: Fourth-order Hammerstein model
close all; clear; clc;

%% The fourth-order LTI model with coloured process noise

% state-space matrices
A = [0.67 0.67 0 0; -0.67 0.67 0 0; 0 0 -0.67 -0.67; 0 0 0.67 -0.67];
B = [0.6598 -0.5256; 1.9698 0.4845; 4.3171 -0.4879; -2.6436 -0.3416];
K = [-0.6968 -0.1474; 0.1722 0.5646; 0.6484 -0.4660; -0.9400 0.1032];
C = [-0.3749 0.0751 -0.5225 0.5830; -0.8977 0.7543 0.1159 0.0982];
D = zeros(2);

% open-loop system
OL = ss(A,[B K],C,[D eye(2)],1);

%% Open-loop identification experiment
% Simulation of the model in open loop

% input signals
N = 1000; % number of samples
t = (0:N-1)';   % time samples
r = randn(N,2); % excitation signal
fu = [sinc(r(:,1)).*r(:,1).^2 sin(r(:,2))];

% noise
e = 0.04.*randn(N,2); % noise signal

% simulation of open loop
y0 = lsim(OL,[fu zeros(N,2)],t);
y = lsim(OL,[fu e],t);
disp('Signal to noise ratio (SNR) (open-loop)')
snr(y,y0)

%%
% Identification of the model in open loop

% parameters
n = 4;    % order of system
f = 10;   % future window size
p = 10;   % past window size

% PBSID-varx
[us,Du,ys,Dy] = sigscale(r,y);
[S,x,fui] = hordvarx(us,ys,f,p,'tikh','gcv');
figure, semilogy(S,'*');
x = hmodx(x,n);
[Ai,Bi,Ci,Di] = hx2abcdk(x,fui,ys,f,p);

%%
% Verification results

% verification using variance accounted for (VAF) (open loop)
Q = fu(p+1:p+size(fui,2),:)'*pinv(fui);
OLi = ss(Ai,(Bi/Q)/Du,Dy*Ci,Dy*((Di/Q)/Du),1);
figure, bodemag(OL(1:2,1:2),'k',OLi,'r');
y = lsim(OL(1:2,1:2),fu,t);
yi = lsim(OLi,fu,t);
disp('VAF (open-loop)')
vaf(y,yi)



% simulation (open loop)
figure, bodemag(OL(1:2,1:2),'k',OLi,'r');
legend('REAL','IDENT');

% plot the non-linear function
fui = Q*Du*fui;
figure, subplot(1,2,1), plot(r(:,1),fu(:,1),'k.',r(p+1:p+size(fui,2),1),fui(1,:)','r.')
subplot(1,2,2), plot(r(:,2),fu(:,2),'k.',r(p+1:p+size(fui,2),2),fui(2,:)','r.')

