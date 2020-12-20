function [e,t] = resid(sys,u,y,t,p,x0,type)
%RESID Compute the residuals associated with an IDAFFLPV.
%  E = RESID(M,U,Y,T,MU) returns the residue response of the IDAFFLPV model
%  M to the input and scheduling signal described by U, Y, MU and T. The
%  time vector T consists of regularly spaced time samples, U, Y, and MU is
%  are matrices with as many columns as inputs and scheduling variables and
%  whose i-th row specifies the input value at time T(i). For discrete-time
%  models, U, Y, and MU should be sampled at the same rate as M.
%
%  E = RESID(M,U,Y,T,MU,X0) specifies the initial state vector X0 at time
%  T(1). X0 is taken from model.
%
%  E = RESID(M,U,Y,T,MU,X0,'Type') specifies the type of LPV predictor. 'K'
%  specifies that K does not dependent on scheduling, or 'CD' specifies
%  that C and D does not dependent on scheduling. Default is Type = 'CD'.
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
if nargout < 7 || isempty(type)
    type = 'CD';
end
if nargout < 5 || isempty(x0)
    x0 = sys.x0;
elseif length(x0)~=Nx
    error('Length of initial condition X0 must match number of states.')
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
    if i == 1
        dt(:,(i-1)*(Nu+Ny)+1:i*(Nu+Ny)) = [d(:,(i-1)*Nu+1:i*Nu) -eye(Ny)];
    else
        dt(:,(i-1)*(Nu+Ny)+1:i*(Nu+Ny)) = [d(:,(i-1)*Nu+1:i*Nu) zeros(ny)];
    end
end
b = bt;
d = dt;
u = [u y];

if isct(sys) % Continuous models  
    if strcmpi(type,'K')
        error('Continuous models with varying C and D is not yet implemented.')
    end
        
    % Determine the states
    tspan = [t(1) t(end)];
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [tc,x] = ode15s(@statesim,tspan,x0,options,t,u,p,a,b);   
    
    % Determine the output 
    e = zeros(length(tc),Ny);
    uc = interp1q(t,u,tc)';
    pc = interp1q(t,p,tc)';
    for i = 1:length(tc)
        e(i,:) = c*(kron([1; pc(:,i)],x(i,:)')) + d*(kron([1; pc(:,i)],uc(:,i)));
    end
   
    % Format output arrays
    e = e.';
    t = tc;
else % Discrete models
    Ts = sys.Ts;
    if length(t)>1 && Ts>0 && abs(t(2)-t(1)-Ts)>1e-4*Ts
        error('Time step must match sample time of discrete-time models.')
    end
    
    % Discrete simulation of LPV systems
    e = lpvsim(a,b,c,d,p,u,[],x0);
end

end

% State-derivative function used for the simulation of continuous models
function dx = statesim(ts,xs,t,u,p,a,b)
    ps = interp1q(t,p,ts)';
    us = interp1q(t,u,ts)';
    dx = a*(kron([1; ps],xs)) + b*(kron([1; ps],us));
end