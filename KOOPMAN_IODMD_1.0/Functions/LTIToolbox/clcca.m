function [A,B,C,D,K,Sn] = clcca(u,y,sysd,k,n,estd)
%CLCCA     This routine implements the second stage of the 2CCA (a
%          subspace algorithm based on Canonical Correlation Analysis)
%          algorithm developed in [1]. The model structure used in the
%          this identification algorithm is the standard linear state
%          space innovations model: 
%          
%          x(k+1)  = A x(k) + B u(k) + K e(k)
%          y(k)    = C x(k) + D u(k) + e(k) 
% 
%          Given input data u, output data y and an initial estimate
%          of the deterministic part (meaning quadruple (A,B,C,D)),
%          this routine re-estimates (A,B,C,D) and gives and estimate
%          of K. 
%          
%          The data can be generated in open- or closed- loop. [1,2]
% 
% Syntax: [A,B,C,D,K,Sn] = clcva(u,y,sysd,k,n,estd)
% 
% Inputs: u      Input data
%         y      Output data 
%         sysd   Estimate of the deterministic part in MATLAB 5 LTI
%                model structure (can be estimated by using other
%                subspace routines) 
%         k      The truncation index (It can be computed by using ARX
%                modelling together with AIC, see ARXSTRUC)
%         N      The model order 
%         estd   Estimate 'D' matrix or not. (1 = yes or 0 = no)
%
% Outputs: A,B,C,D,K   state space matrices 
%          Sn          singular values 
% 
% See also: dordpo, dmodpo (subapace routines for open-loop data) 
%           dordeiv, dmodeiv, dac2bd_eiv, dac2b_eiv, dac2bd_cl, dac2b_cl
%           (subspace routines for closed-loop data) arxstruc (for
%           selection of truncation index) 
% 
% References: [1] K. Peternell, W. Scherrer and  M. Deistler.
%                 Statistical analysis of novel subspace
%                 identification method. Signal Processing, vol. 52,
%                 pp. 161-177, 1996.
%             [2] C.T. Chou and M. Verhaegen. Closed-loop
%                 identification using canonical correlation
%                 analysis. To appear in Proc. ECC '99. 

% C.T. Chou, Sept 1998, May 1999. 
% Revised by Ivo Houtzager, 2009
% Copyright (c) 1998-2009, Delft Center of Systems and Control 

% get dimension parameters
[Nd,m] = size(u);
p = size(y,2);

% construct Hankel and Toeplitz matrices
% number of columns in the Hankel and Toeplitz matrices
Nh = Nd - 2*k + 1;
% Form Toeplitz matrices using past input and output
Nr = (m+p)*k;  % number of rows
Zp = zeros((m+p)*k,Nh);
skip = (m+p);
for i = 1:p    % Yp
    Zp(  i:skip:Nr,:) = toeplitz(y(k:-1:1,i),y(k:Nd-k,i));
end
for i = 1:m    % Up
    Zp(i+p:skip:Nr,:) = toeplitz(u(k:-1:1,i),u(k:Nd-k,i));
end
% Form Hankel matrices using future input and output
Nr = m*k; % number of rows
Uf = zeros(Nr,Nh);
for i = 1:m
    Uf(i:m:Nr,:) = hankel(u(k+1:2*k,i),u(2*k:Nd,i));  % Uf
end
Nr = p*k;
Yf = zeros(Nr,Nh);
for i = 1:p
    Yf(i:p:Nr,:) = hankel(y(k+1:2*k,i),y(2*k:Nd,i));  % Yf
end

% the Markov parameters
sysd.Ts = -1;
impd = impulse(sysd,k); % this gives k impusle response
buhat = zeros(k*p,k*m);
% fill in the first block column, this depends on ...
if (m == 1) && (p == 1)
    buhat(:,1) = impd;
elseif (m == 1) && (p > 1)
    buhat(:,1) = reshape(impd',k*p,1);
else
    impd = shiftdim(impd,1);
    for i = 1:k
        buhat((i-1)*p+1:i*p,1:m) = squeeze(impd(:,:,i));
    end
end
% fill in the rest of the block columns
for i = 2:k
    buhat((i-1)*p+1:k*p,(i-1)*m+1:i*m) = buhat(1:(k-i+1)*p,1:m);
end

% the canonical variables are Yf - buhat * Uf and Zp
Ym = Yf - buhat * Uf;
bzhat = Ym/Zp;
% weighting W1
R_Ym = triu(qr(Ym'));
R_Ym = R_Ym(1:k*p,1:k*p);
[U_R_Ym,S_R_Ym,V_R_Ym] = svd(R_Ym);
W1inv = V_R_Ym*S_R_Ym*V_R_Ym';
W1 = inv(W1inv);

% weighting W2
R_Zp = triu(qr(Zp'));
R_Zp = R_Zp(1:(m+p)*k,1:(m+p)*k);
[U_R_Zp,S_R_Zp,V_R_Zp] = svd(R_Zp);
W2 = V_R_Zp*S_R_Zp*V_R_Zp';

% svd of weighted estimate of bzhat
[U,S,V] = svd(W1*bzhat*W2);
Sn = diag(S);

% estimation of Kk
Kkhat = sqrt(S(1:n,1:n))*V(:,1:n)'/W2;

% estimation of states
Xhat = Kkhat*Zp;

% estimation of state space matrices
% C and D matrices and innovations
Y1 = Yf(1:p,:);
U1 = Uf(1:m,:);
if estd
    CD = Y1/[Xhat ; U1];
    C = CD(:,1:n);
    D = CD(:,n+1:n+m);
    E = Y1-C*Xhat-D*U1;
else
    C = Y1/Xhat;
    D = zeros(p,m);
    E = Y1-C*Xhat;
end

% estimation of (A,B,K)
Xhat_shifted = zeros(size(Xhat));
Xhat_shifted(:,1:end-1) = Xhat(:,2:end);
ABK = Xhat_shifted/[Xhat ; U1 ; E];
A = ABK(:,1:n);
B = ABK(:,n+1:n+m);
K = ABK(:,n+m+1:n+m+p);







