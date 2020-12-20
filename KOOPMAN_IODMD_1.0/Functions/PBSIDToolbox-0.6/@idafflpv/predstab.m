function tau = predstab(sys,mumin,mumax,type,resize,showeigs)
%PREDSTAB Assesses the quadratic stability of the predictor form of affine
% LPV system.
% tau = predstab(sys,mumin,mumax,type) Assesses the quadratic 
% stability of the predictor form of the LPV model:
% 
%     x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K kron(mu(k),e(k))
%     y(k)   = C x(k) + D u(k) + e(k)
%
% given by sys and returns tau = 1 if quadratic stability holds for mu(k) 
% varying in [mumin,mumax].
%
% type = 'K' specifies that K is not dependent on the scheduling, or 'CD' 
% specifies that C and D is not dependent on the scheduling. 
% Default is type = 'CD'.
% 
% tau = predstab(sys,mumin,mumax,type,1), returns the largest tau > 0 
% such that quadratic stability is guaranteed for the scheduling varying 
% in the range:
%
%    mu0-dmu*tau <= mu(k) <= mu0+dmu*tau
%
% where mu0 = mumin+mumax./2, dmu = mumax-mumin./2.
%
% Hence, the system
%
%  xhat(k+1) = (A-CK) kron(mu(k),xhat(k)) + (B-DK) kron(mu(k),u(k)) + K kron(mu(k),y(k))
%  yhat(k) = C x(k) + D u(k)
%
% is quadratically stable for mu in [mumin,mumax] if tau>=1.
% 
% When showeigs = 1 (default showeigs = 0), it shows the spectral radius of
% the predictor state matrix at each cornerpoint of the scheduling.

%  Pieter Gebraad
%  Delft Center of Systems and Control
%  Delft University of Technology
%  The Netherlands, 2011
    
    if nargin < 4 || isempty(type)
        type = 'CD';
    end

    if nargin<6 || isempty(showeigs)
        showeigs = 0;
    end
    
    if nargin<5 || isempty(resize)
        resize = 0;
    end
    
    % check if YALMIP is installed
    if exist('yalmip','file') ~= 2
        error('predstab:yalmipnotfound','The function PREDSTAB uses YALMIP and an SDP solver (e.g. SEDUMI) to check the stability of the LPV predictor form. \n YALMIP is not found in your MATLAB path. \n Please install YALMIP and an SDP solver and add it to your MATLAB path \n YALMIP is available at: http://users.isy.liu.se/johanl/yalmip \n A free SDP solver called SEDUMI is available at: http://sedumi.ie.lehigh.edu');
    end
    
    % Get the system matrices of the predictor
    [A,~,C,~,K] = getABCDK(sys);
    [Ny,Nu,Nx,Np] = size(sys);

    if strcmpi(type,'CD')
        for i = 1:Np
            A(:,(i-1)*Nx+1:i*Nx) = A(:,(i-1)*Nx+1:i*Nx) - K(:,(i-1)*Ny+1:i*Ny)*C(:,1:Nx);
        end
    elseif strcmpi(type,'K')
        for i = 1:Np
            A(:,(i-1)*Nx+1:i*Nx) = A(:,(i-1)*Nx+1:i*Nx) - K(:,1:Ny)*C(:,(i-1)*Nx+1:i*Nx);
        end
    else
        error('Type not recognized!')
    end
    
    if (length(mumin)~=Np)
        error 'Incorrect size of mumin';
    end
    if (length(mumin)~=Np)
        error 'Incorrect size of mumax';
    end
    
    if size(mumin,2)~=1
        mumin = mumin';
    end
    if size(mumax,2)~=1
        mumax = mumax';
    end
    
    % check stability for corner points of convex polytope
    much = [mumin,mumax];
    much = allcomb(much);
    tau = quadstab(A,much,Nx,Np,showeigs);
       
    % bisection algorithm to find largest possible tau
    
    ind = ones(Np,1);
    ind = allcomb([ind,2*ind]);
    ind = sub2ind(size(ind),(1:Np)'*ones(1,size(ind,2)),ind);
    
    mu0 = mumin+mumax./2; dmu = mumax-mumin./2;
    if resize
        if ~tau % if unstable, halve tau until stable
            tau=1;
            st = 0;
            while ~st && (tau>1e-8)
                tau2 = tau;
                tau = tau/2;
                range = [mu0-dmu*tau,mu0+dmu*tau];
                much = range(ind);
                st = quadstab(A,much,Nx,Np,0);
            end
            tau1 = tau;
            if ~st
               disp 'Unable to stabilize for very small tau';
               resize = 0;
            end
        else % if stable, double tau until unstable
            st = 1;
            while st && (tau<1e5)
                tau1 = tau;
                tau = tau*2;
                much = [mu0-dmu*tau,mu0+dmu*tau];
                much = allcomb(much);
                st = quadstab(A,much,Nx,Np,0);
            end
            tau2 = tau;
            if st
               disp 'Unable to stabilize for very large tau';
               resize = 0;
            end
        end
    end
    if resize
        % bisection
        while (tau2-tau1) > 1e-5
            tau = (tau1+tau2)/2;
            range = [mu0-dmu*tau,mu0+dmu*tau];
            much = range(ind);
            if(quadstab(A,much,Nx,Np,0))
                tau1 = tau;
            else
                tau2 = tau;
            end
        end
        tau = tau1;
    end
    
    function tau = quadstab(A,much,Nx,Np,showeigs)
        tau = 1;
        %check stability on cornerpoints of polytope
        %yalmip clear; options = sdpsettings('verbose',0);
        %P = sdpvar(Nx,Nx);
        %LMI = set(P>0);
        for j = 1:size(much,2)
            Ap = A(1:Nx,1:Nx);
            for c = 1:Np
              Ap = Ap + much(c,j)*A(:,c*Nx+1:(c+1)*Nx);
              sr = max(abs(eig(Ap)));
            end
            if showeigs
                disp(['cornerpoint:     mu = ',num2str(much(:,j)','%9.2g')]);
                disp(['spectral radius A-KC: ',num2str(sr','%9.2g')]);
            end
            %LMI =LMI + set(Ap'*P*Ap-P<0);
            if sr>1
                tau = 0;
            end
        end
        %solvesdp(LMI,[],options);
        %primal=checkset(LMI);
        %tau = all(primal>0);
        
    end
end

