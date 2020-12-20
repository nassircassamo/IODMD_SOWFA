%% Example 19c: Using sparse estimation in the identification of a second-order LPV model
% Using the LPV system of example 19, we demonstate the regularization
% options of |lordvarx|, more in particular the BPDN approach for
% regularization, which can reduce the sensitivity of the performance of
% the algorithm to the choice of the past window parameter _p_.

close all; clear; clc;

%% Flapping dynamics of a wind turbine

% System matrices
A1 = [0 0.0734; -6.5229 -0.4997];
A2 = [-0.0021 0; -0.0138 0.5196];
A12 = [A1 A2];
B12 = [-0.7221 0; -9.6277 0];
C12 = [1 0 0 0];
D12 = [0 0];
n = size(A12,1);    % The order of the system
m = size(A12,2)/n;  % The number of scheduling parameters
r = size(B12,2)/m;  % The number of inputs
l = size(C12,1);    % The number of outputs

%% Obtaining the identification data
% Simulation of the model in open loop

% defining a number of constants
j = 10;    % period
np = 7;    % number of periods
N = np*j;  % number of data points

% measured data from the scheduling parameters
mu3 = cos(2*pi*(1:N)'./j)+0.2;

% make affine LPV system
M = idafflpv(A12,B12,C12,D12,eye(2),zeros(2,1),1);

% simulation of the system with noise, to obtain identification data
t = (0:N-1)';
u = randn(N,r);
e = 0.05.*randn(N,l);
y0 = sim(M,u,t,mu3);
y = sim(M,u,t,mu3,e);
disp('Signal to noise ratio (SNR) (open-loop)')
snr(y,y0)

% simulation of the system without noise, to obtain validation data
uval = randn(N,r);
yval = sim(M,uval,t,mu3);

%% LPV identification with PBSID using Tikhonov regularization
% We study the effect of the past window parameter _p_ to the performance of
% the identification algorithm when we use Tikhonov regularization. We use
% Generalized Cross-Validation (GCV) to select the trade-off parameter in
% the regularization.

disp('----------------------------------------------------------------------------------------------');
disp(' past window p     VAF on identification data     VAF on validation data     Calculation Time');

c = warning('query','regress:RankDefDataMat');
warning('off','regress:RankDefDataMat')
prange = 3:16;
for count = 1:length(prange)
    p = prange(count);  % past window size
    f = p;              % future window size
    % LPV identification with noise
    mu = [ones(N,1) mu3];
    tic
    [S,x] = lordvarx(u,y,mu,f,p,'tikh','gcv',[0 1 0]);
    time = toc;
    x = lmodx(x,n);
    [A,B,C,D,K] = lx2abcdk(x,u,y,mu,f,p,[0 1 0]);
    Mk = idafflpv(A,B,C,D,K,zeros(2,1),1);
    % Simulation of identified LPV system
    yidk = sim(Mk,u,t,mu3);   vafid = vaf(y0,yidk); % with identification data
    yvalk = sim(Mk,uval,t,mu3);  vafval = vaf(yval,yvalk); % with validation data
    disp([sprintf('%8.u',p),sprintf('%25.1f',vafid),'%',sprintf('%27.1f',vafval),'%',sprintf('%28.2e',time),'s']);
end
%%
% We see that the performance of the identification algorithm, expressed in
% the Variance Accounted For on validation data, decreases with the past
% window if we use a fixed ammount of data points, as the parameter 
% estimation problem becomes ill-conditioned, since the number of
% parameters to be estimated increases.

%% LPV identification with PBSID using BPDN regularization
% By using Basis Pirsuit DeNoising (BPDN) in the parameter estimation 
% problem in the PBSID algorithm, we improve its conditioning. BPDN sets 
% insignificant parameters to zero, making the estimated model less 
% sensitive to noise, possibly at the cost of a bias error. We see 
% that increasing the past window beyond its optimal value, will not 
% degrade the quality of the model as much as in the previous case.

c = warning('query','lordvarx:BpdnThenNoKernel');
warning('off','lordvarx:BpdnThenNoKernel')
disp('----------------------------------------------------------------------------------------------');
disp(' past window p     VAF on identification data     VAF on validation data     Calculation Time');
for count = 1:length(prange)
    p = prange(count);  % past window size
    f = p;              % future window size
    % LPV identification with noise
    mu = [ones(N,1) mu3];
    tic
    [S,x] = lordvarx(u,y,mu,f,p,'bpdn','sv',[0 1 0]);
    time = toc;
    x = lmodx(x,n);
    [A,B,C,D,K] = lx2abcdk(x,u,y,mu,f,p,[0 1 0]);
    Mk = idafflpv(A,B,C,D,K,zeros(2,1),1);
    % Simulation of identified LPV system
    yidk = sim(Mk,u,t,mu3);   vafid = vaf(y0,yidk); % with identification data
    yvalk = sim(Mk,uval,t,mu3);  vafval = vaf(yval,yvalk); % with validation data
    disp([sprintf('%8.u',p),sprintf('%25.1f',vafid),'%',sprintf('%27.1f',vafval),'%',sprintf('%28.2e',time),'s']);
end
warning(c.state,'lordvarx:BpdnThenNoKernel');
%%
% The option 'sv' means that in the BPDN solver, a trade-off between 
% sparsity of the solution and the residual error is made automatically on 
% the basis of a part of the data not used for the regression (i.e. 
% validation data). By default, the last quarter of the data sequences 
% entered to LORDVARX is used as validation data.
%% Conclusion
% Using BPDN regularization in |lordvarx| will eliminate the need to 'tune'
% the past window parameter _p_ , possibly at the cost of a bias error. It is 
% especially useful if we want to identify models from a relatively small
% ammount of data.
% A downside of the BPDN approach is the increased calculation time.