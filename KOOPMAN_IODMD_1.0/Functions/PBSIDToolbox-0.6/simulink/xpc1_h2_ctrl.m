% Identification of 2D Airfoil
clc; clear; close all;

% High frequency roll-off filter
omega=10*pi*2;
filter = tf([omega^2],[1 omega*2*0.5 omega^2]);

% Load Data
load model.mat
P = ss(Mi);
P = d2c(P);
P = series([filter 0 0; 0 1 0; 0 0 1],P);
figure, bodemag(P);
figure, sigma(P);

% H2 design
GP = [P; 0.01 0 0];
GP = minreal(GP([1:3 1:2],[2:3 1]));
[C,CL,gam,info] = h2syn(GP,2,1);
gam
figure, bodemag(C);
figure, sigma(C);

% disturbance
t = (0:0.01:10)';
u1 = [5*sin(t.*4.7) 5*sin(t.*4.7)];
u2 = [5*sin(t.*17.1) 5*sin(t.*17.1+pi)];

% Rescale;
H = P(1:2,2:3);
P = P(1:2,1);

% Closed-loop
CL = feedback(P,C,+1);
figure, pzmap(CL);
figure, bodemag(CL);
figure, sigma(CL);

% Sensitivity function
S = minreal(inv(eye(2) - P*C));
figure, bodemag(S);
figure, sigma(S);
figure, lsim(S,u1,t);
figure, lsim(S,u2,t);

% Noise Sensitivity function
N = minreal(inv(eye(2) - P*C))*H;
figure, bodemag(N,H);
figure, sigma(N,H);

% Input Sensitivity function
I = minreal(inv(1 - C*P))*C;
figure, bodemag(I);
figure, sigma(I);
figure, lsim(I,u1,t);
figure, lsim(I,u2,t);

% resample
Ts = 0.005;
Fs = c2d(filter,Ts,'tustin');
Cs = c2d(C,Ts,'tustin');
figure, bodemag(C*filter,Cs*Fs);

% save controller
save h2ctrl.mat Fs Cs

% State-feedback regulator design for antiwindup
O = ss(Cs.a',Cs.c',Cs.b',Cs.d',Cs.ts);
Q = eye(2);
R = 100;
[L,S,E] = lqry(O,Q,R);
L = L';
Cs = ss(Cs.a,[Cs.b L],Cs.c,[Cs.d 0],Cs.ts);

% save controller
save h2ctrl_aw.mat Fs Cs
