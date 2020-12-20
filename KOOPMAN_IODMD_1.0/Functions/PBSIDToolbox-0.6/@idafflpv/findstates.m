function x0 = findstates(sys,u,y,t,p,type,wn)
%FINDSTATES  Estimate initial states of the model for a given data set.
%  X0 = FINDSTATES(M,U,Y,T,MU) returns the residue response and initial
%  state of the IDAFFLPV model M to the input and scheduling signal
%  described by U, Y, MU and T. The time vector T consists of regularly
%  spaced time samples, U, Y, and MU is are matrices with as many columns
%  as inputs and scheduling variables and whose i-th row specifies the
%  input value at time T(i). For discrete-time models, U, Y, and MU should
%  be sampled at the same rate as M.
%
%  X0 = FINDSTATES(M,U,Y,T,MU,'Type') specifies the type of LPV predictor.
%  'K' specifies that K is not dependent on scheduling, or 'CD' specifies
%  that C and D is not dependent on scheduling. Default is Type = 'CD'.
%
%  X0 = FINDSTATES(M,U,Y,T,MU,'Type',P) specifies the past window P. The
%  default is P = 5*size(M.a,1).
%
%  NOTE: This is only possible if K or C and D are not dependent on the
%  scheduling sequence.

% Define sizes
N = length(t);
if size(u,1) < size(u,2);
    u = u';
end
if size(y,1) < size(y,2);
    y = y';
end
if size(t,1) < size(t,2);
    t = t';
end
if size(p,1) < size(p,2);
    p = p';
end
[Ns,nu] = size(u);
[Nn,np] = size(p);
[Ni,ny] = size(y);
[Ny,Nu,Nx,Np] = size(sys);

% Computability and consistency checks
if ~isequal(N,Ns,Nn,Ni)
    error('Number of samples in vector T, U and P must be equal.')
end
if ny ~= Ny
    error('Input data Y must have as many columns as system outputs.')
end
if nu ~= Nu
    error('Input data U must have as many columns as system inputs.')
end
if np ~= Np
    error('Input data P must have as many columns as scheduling parameters.')
end

% Assign values to unspecified parameters
if nargin < 7 || isempty(wn)
    wn = 5*Nx;
end
if nargin < 6 || isempty(type)
    type = 'CD';
end

% Get the system matrices
[a b c d k] = getABCDK(sys);
if strcmpi(type,'CD')
    for i = 1:Np+1
        a(:,(i-1)*Nx+1:i*Nx) = a(:,(i-1)*Nx+1:i*Nx) - k(:,(i-1)*Ny+1:i*Ny)*c(:,1:Nx);
        b(:,(i-1)*Nu+1:i*Nu) = b(:,(i-1)*Nu+1:i*Nu) - k(:,(i-1)*Ny+1:i*Ny)*d(:,1:Nu);
    end
elseif strcmpi(type,'K')
    for i = 1:Np+1
        a(:,(i-1)*Nx+1:i*Nx) = a(:,(i-1)*Nx+1:i*Nx) - k(:,1:Ny)*c(:,(i-1)*Nx+1:i*Nx);
        b(:,(i-1)*Nu+1:i*Nu) = b(:,(i-1)*Nu+1:i*Nu) - k(:,1:Ny)*d(:,(i-1)*Nu+1:i*Nu);
    end
else
    error('Type not recognized!')
end
bt = zeros(Nx,(Np+1)*(Nu+Ny));
dt = zeros(Ny,(Np+1)*(Nu+Ny));
for i = 1:Np+1
    bt(:,(i-1)*(Nu+Ny)+1:i*(Nu+Ny)) = [b(:,(i-1)*Nu+1:i*Nu) k(:,(i-1)*Ny+1:i*Ny)];
    dt(:,(i-1)*(Nu+Ny)+1:i*(Nu+Ny)) = [d(:,(i-1)*Nu+1:i*Nu) zeros(ny)];
end
b = bt;
d = dt;
u = [u y];

if isct(sys) % Continuous models
    yr = y;   
    options = optimset('Display','off');
    x0 = lsqnonlin(@residue,ones(Nx,1),[],[],options,u(1:wn,:),yr(1:wn,:),t(1:wn,:),p(1:wn,:),a,b,c,d,type);
else
    % Get the system matrices
    [Altv,Bltv,Cltv,Dltv] = lpv2ltv(a,b,c,d,p(1:wn,:));
    [HU,Gamma] = lift(Altv,Bltv,Cltv,Dltv,wn);
    Y = y(1:wn,:)';
    Y = Y(:);
    U = u(1:wn,:)';
    U = U(:);
    Y = Y - HU*U;
    x0 = pinv(Gamma)*Y;
end
end

% State-derivative function used for the simulation of continuous models
function cost = residue(x0,u,yr,t,p,a,b,c,d,type)
Ny = size(yr,2);

if strcmpi(type,'K')
    error('Continuous models with varying C and D is not yet implemented.')
end

% Determine the states
tspan = [t(1) t(end)];
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[tc,x] = ode15s(@statesim,tspan,x0,options,t,u,p,a,b);

% Determine the output
y = zeros(length(tc),Ny);
uc = interp1q(t,u,tc)';
pc = interp1q(t,p,tc)';
for i = 1:length(tc)
    y(i,:) = c*(kron([1; pc(:,i)],x(i,:)')) + d*(kron([1; pc(:,i)],uc(:,i)));
end

% Format output arrays
y = y.';
t = tc;
yc = interp1q(t,yr,tc)';
cost = pec(yc,y);
end

% State-derivative function used for the simulation of continuous models
function dx = statesim(ts,xs,t,u,p,a,b)
    ps = interp1q(t,p,ts)';
    us = interp1q(t,u,ts)';
    dx = a*(kron([1; ps],xs)) + b*(kron([1; ps],us));
end

function [HU,Gamma] = lift(A,B,C,D,p)
%LIFT Lift state-spave matrices
% written by, I. Houtzager [2007]
% Delft Center of Systems and Control

% determine lpv system sizes
r = size(B{1},2); % The number of inputs
l = size(C{1},1); % The number of outputs
n = size(A{1},1); % The number of states

% build lifted impulse matrix
HU = zeros(p*l,p*r);
Gamma = zeros(p*l,n);
for i = 1:p
    for j = 1:p
        if i == j
            if isempty(D)
                HU((i-1)*l+1:i*l,(j-1)*r+1:j*r) = zeros(l,r);
            else
                HU((i-1)*l+1:i*l,(j-1)*r+1:j*r) = D{i};
            end
        elseif j < i
            if j == i-1
                HU((i-1)*l+1:i*l,(j-1)*r+1:j*r) = C{i}*B{j};
            else
                T = C{i};
                for k = i-1:-1:j+1
                    T = T*A{k};
                end
                HU((i-1)*l+1:i*l,(j-1)*r+1:j*r) = T*B{j};
            end
        end
    end
    T = C{i};
    for k = i-1:-1:1
        T = T*A{k};
    end
    Gamma((i-1)*l+1:i*l,:) = T;
end
end
