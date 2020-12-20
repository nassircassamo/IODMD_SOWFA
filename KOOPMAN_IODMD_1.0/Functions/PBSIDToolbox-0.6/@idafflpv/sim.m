function [y,t,x] = sim(sys,u,t,p,e,x0)
%SIM Linear response simulation of affine LPV state-space model.
%  [Y,T,X] = SIM(M,U,T,MU) returns the output response of the IDAFFLPV
%  model M to the input and scheduling signal described by U, MU and T.
%  The time vector T consists of regularly spaced time samples, U and MU is
%  are matrices with as many columns as inputs and scheduling variables and
%  whose i-th row specifies the input value at time T(i). For discrete-time
%  models, U should be sampled at the same rate as M.
%
%  [Y,T,X] = SIM(M,U,T,MU,E) adds the innovation noise to the simulation.
%
%  [Y,T,X] = SIM(M,U,T,MU,E,X0) specifies the initial state vector X0 at
%  time T(1). X0 is set to zero when omitted.

% Get the system matrices
[a b c d k] = getABCDK(sys);

% Define sizes
N = length(t);
Ny = size(c,1);
[Ns,nu] = size(u);
[Nn,np] = size(p);
if size(u,1) < size(u,2);
    u = u';
end
if size(t,1) < size(t,2);
    t = t';
end
if size(p,1) < size(p,2);
    p = p';
end
if ~(nargin < 5 || isempty(e))
    if size(e,1) < size(e,2);
        e = e';
    end
else
    e = zeros(N,Ny);
end
[Ny,Nu,Nx,Np] = size(sys);

% Computability and consistency checks
if ~isequal(N,Ns,Nn)
    error('Number of samples in vector T, U and P must be equal.')
end
if nu ~= Nu
    error('Input data U must have as many columns as system inputs.')
end
if np ~= Np
    error('Input data P must have as many columns as scheduling parameters.')
end

% Assign values to unspecified parameters
if nargin < 6 || isempty(x0)
    x0 = zeros(Nx,1);
elseif length(x0)~=Nx
    error('Length of initial condition X0 must match number of states.')
end

if isct(sys) % Continuous models             
    % Determine the states
    tspan = [t(1) t(end)];
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [tc,x] = ode15s(@statesim,tspan,x0,options,t,u,p,e,a,b,k);   
    
    % Determine the output 
    y = zeros(length(tc),Ny);
    uc = interp1q(t,u,tc)';
    pc = interp1q(t,p,tc)';
    ec = interp1q(t,e,tc)';
    for i = 1:length(tc)
        y(i,:) = c*(kron([1; pc(:,i)],x(i,:)')) + d*(kron([1; pc(:,i)],uc(:,i))) + ec(:,i);
    end
   
    % Format output arrays
    y = y.';
    t = tc;
    
else % Discrete models
    Ts = sys.Ts;
    if length(t)>1 && Ts>0 && abs(t(2)-t(1)-Ts)>1e-4*Ts
        error('Time step must match sample time of discrete-time models.')
    end
    
    bt = zeros(Nx,(Np+1)*(Nu+Ny));
    dt = zeros(Ny,(Np+1)*(Nu+Ny));
    for i = 1:Np+1
        bt(:,(i-1)*(Nu+Ny)+1:i*(Nu+Ny)) = [b(:,(i-1)*Nu+1:i*Nu) k(:,(i-1)*Ny+1:i*Ny)];
        if i == 1
            dt(:,(i-1)*(Nu+Ny)+1:i*(Nu+Ny)) = [d(:,(i-1)*Nu+1:i*Nu) eye(Ny)];
        else
            dt(:,(i-1)*(Nu+Ny)+1:i*(Nu+Ny)) = [d(:,(i-1)*Nu+1:i*Nu) zeros(Ny)];
        end
    end
    [y,x] = lpvsim(a,bt,c,dt,p,[u e],[],x0);
end

end

% State-derivative function used for the simulation of continuous models
function dx = statesim(ts,xs,t,u,p,e,a,b,k)
    ps = interp1q(t,p,ts)';
    us = interp1q(t,u,ts)';
    es = interp1q(t,e,ts)';
    dx = a*(kron([1; ps],xs)) + b*(kron([1; ps],us))+ k*(kron([1; ps],es));
end

