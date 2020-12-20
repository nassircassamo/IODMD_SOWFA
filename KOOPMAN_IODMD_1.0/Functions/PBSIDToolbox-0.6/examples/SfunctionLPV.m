function [Sys,x0]=SfunctionLPV(t,x,u,flag,A,B,C,D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define topology model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumContStates  = 0;
NumDiscStates  = 7;
NumOutputs     = 6;
NumInputs      = 10;
DirFeedthrough = 0;
NumSampleTimes = 1;   % at least one sample time is needed

n = 7;
r = 7;

A0 = A(:,1:n);
A1 = A(:,n+1:2*n);
A2 = A(:,2*n+1:3*n);
A3 = A(:,3*n+1:4*n);
B0 = B(:,1:r);
B1 = B(:,r+1:2*r);
B2 = B(:,2*r+1:3*r);
B3 = B(:,3*r+1:4*r);
C0 = C(:,1:n);
C1 = C(:,n+1:2*n);
C2 = C(:,2*n+1:3*n);
C3 = C(:,3*n+1:4*n);
D0 = D(:,1:r);
D1 = D(:,r+1:2*r);
D2 = D(:,2*r+1:3*r);
D3 = D(:,3*r+1:4*r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if abs(flag) == 0,	% return sizes of parameters and initial conditions
   
Sys = [NumContStates NumDiscStates NumOutputs NumInputs DirFeedthrough NumSampleTimes];

x0 = [0 0 0 0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update of discrete states  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif abs(flag) == 2,

% rename inputs
mu1 = u(8);
mu2 = u(9);
mu3 = u(10);
u = u(1:7);

% define model equations
Sys = A0*x + (A1.*mu1)*x + (A2.*mu2)*x + (A3.*mu3)*x +...
      B0*u + (B1.*mu1)*u + (B2.*mu2)*u + (B3.*mu3)*u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif abs(flag) == 3,	% return output vector

% rename inputs
mu1 = u(8);
mu2 = u(9);
mu3 = u(10);
u = u(1:7);    
    
Sys = C0*x + (C1.*mu1)*x + (C2.*mu2)*x + (C3.*mu3)*x +...
      D0*u + (D1.*mu1)*u + (D2.*mu2)*u + (D3.*mu3)*u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
else	% no need to return anything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sys=[];
end 