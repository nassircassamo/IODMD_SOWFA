function sys = pem(sys,u,y,t,p,type,options)
%PEM Computes the prediction error estimate of an IDAFFLPV.
%  M = PEM(M,U,Y,T,MU) returns optimized model of the IDAFFLPV model M to
%  the input and scheduling signal described by U, Y, MU and T. The time
%  vector T consists of regularly spaced time samples, U, Y, and MU is are
%  matrices with as many columns as inputs and scheduling variables and
%  whose i-th row specifies the input value at time T(i). For discrete-time
%  models, U, Y, and MU should be sampled at the same rate as M.
%
%  NOTE: This is only possible if K or C and D are not dependent on the
%  scheduling sequence.
%
%  NOTE: Note that zero elements in M are NOT optimized (are constrained).
%
%  M = PEM(M,U,Y,T,MU,'Type') specifies the type of LPV predictor. 'K'
%  specifies that K does not dependent on scheduling, or 'CD' specifies
%  that C and D does not dependent on scheduling. Default is Type = 'CD'.
%
%  M = PEM(M,U,Y,T,MU,'Type',options) specifies options for the non-linear
%  least sqaures solver. The following options need to be set:
%
%     options.Display = {0,1}
%     options.TolX = 1e-4
%     options.TolFun = 1e-4
%     options.Jacobian = {0,1}
%     options.MaxFunEval = 200

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
if nargin < 7 || isempty(options)
    fopts = zeros(1,18);
else
	fopts = zeros(1,18);
	fopts(1) = options.Display;
	fopts(2) = options.TolX;
	fopts(3) = options.TolFun;
	fopts(9) = options.Jacobian;
	fopts(14) = options.MaxFunEval;	
end
if nargin < 6 || isempty(type)
    type = 'CD';
end

if isct(sys) % Continuous models
    error('Continuous systems not implemented yet.')
else
    % Get the system matrices
    [a b c d k] = getABCDK(sys);
    if strcmpi(type,'CD')
        for i = 1:Np+1
            a(:,(i-1)*Nx+1:i*Nx) = a(:,(i-1)*Nx+1:i*Nx) - k(:,(i-1)*Ny+1:i*Ny)*c(:,1:Nx);
            b(:,(i-1)*Nu+1:i*Nu) = b(:,(i-1)*Nu+1:i*Nu) - k(:,(i-1)*Ny+1:i*Ny)*d(:,1:Nu);
        end
        c(:,Nx+1:(Np+1)*Nx) = zeros(Ny,Np*Nx);
        d(:,Nu+1:(Np+1)*Nu) = zeros(Ny,Np*Nu);
    elseif strcmpi(type,'K')
        for i = 1:Np+1
            a(:,(i-1)*Nx+1:i*Nx) = a(:,(i-1)*Nx+1:i*Nx) - k(:,1:Ny)*c(:,(i-1)*Nx+1:i*Nx);
            b(:,(i-1)*Nu+1:i*Nu) = b(:,(i-1)*Nu+1:i*Nu) - k(:,1:Ny)*d(:,(i-1)*Nu+1:i*Nu);
        end
        k(:,Ny+1:(Np+1)*Ny) = zeros(Ny,Np*Ny);
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
    
    % Get parameters
    [th,dim]=lpv2par(a,b,c,d);
    ac = (a == 0);
    bc = (b == 0);
    cc = (c == 0);
    dc = (d == 0);
    thc= lpv2par(ac,bc,cc,dc);
    th = lpvopt(th,thc,dim,u,y,p,fopts);
    
    % Get the system matrices
    [a,b,c,d]=par2lpv(th,dim);
    bt = zeros(Nx,(Np+1)*Nu);
    dt = zeros(Ny,(Np+1)*Nu);
    k = zeros(Nx,(Np+1)*Ny);
    for i = 1:Np+1
        bt(:,(i-1)*Nu+1:i*Nu) = b(:,(i-1)*(Nu+Ny)+1:(i-1)*(Nu+Ny)+Nu);
        dt(:,(i-1)*Nu+1:i*Nu) = d(:,(i-1)*(Nu+Ny)+1:(i-1)*(Nu+Ny)+Nu);
        k(:,(i-1)*Ny+1:i*Ny) = b(:,(i-1)*(Nu+Ny)+Nu+1:i*(Nu+Ny));
    end
    d = dt;
    b = bt;
    if strcmpi(type,'CD')
        for i = 1:Np+1
            a(:,(i-1)*Nx+1:i*Nx) = a(:,(i-1)*Nx+1:i*Nx) + k(:,(i-1)*Ny+1:i*Ny)*c(:,1:Nx);
            b(:,(i-1)*Nu+1:i*Nu) = b(:,(i-1)*Nu+1:i*Nu) + k(:,(i-1)*Ny+1:i*Ny)*d(:,1:Nu);
        end
    elseif strcmpi(type,'K')
        for i = 1:Np+1
            a(:,(i-1)*Nx+1:i*Nx) = a(:,(i-1)*Nx+1:i*Nx) + k(:,1:Ny)*c(:,(i-1)*Nx+1:i*Nx);
            b(:,(i-1)*Nu+1:i*Nu) = b(:,(i-1)*Nu+1:i*Nu) + k(:,1:Ny)*d(:,(i-1)*Nu+1:i*Nu);
        end
    else
        error('Type not recognized!')
    end
    sys.a = a;
    sys.b = b;
    sys.c = c;
    sys.d = d;
    sys.k = k;
    [e,x0] = pe(sys,u(:,1:Nu),y,t,p,type);
    sys.x0 = x0;
    sys.NoiseVariance = cov(e);
end

