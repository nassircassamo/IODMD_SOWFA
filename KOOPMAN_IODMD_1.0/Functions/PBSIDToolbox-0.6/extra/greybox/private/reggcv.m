function reg_min=reggcv(Y,Vn,Sn,method,show)
%REGGCV      Compute regularization using generalized cross validation.
%            Determine the regularization parameter for ordkernel
%            using Generalized Cross-Validation (GCV). It plots the
%            GCV function as a function of the regularization
%            parameter and finds its minimum.
%
% Syntax:
%            reg=reggcv(Y,V,S)
%            reg=reggcv(Y,V,S,method,show)
%
% Input:
%            Y,V,S    Data matrices from lpvkernel or bilkernel.
%            method   Regularization method to be used.
%                     'Tikh' - Tikhonov regularization (default).
%                     'tsvd' - Truncated singular value decomposition.
%            show     Display intermediate steps of the algorithm.
%
% Output:
%            reg      Regularization paramater for the kernel
%                     subspace identification method of ordkernel.

% Written by Vincent Verdult, May 2004.
% Based on Regularization Tools by P. C. Hansen

% default method
if nargin<4
    method='Tikh';
end
if nargin<5
    show=0;
end
if size(Y,1)~=size(Vn,1)
    error('The number of rows in Y must equal the number of rows in V.')
end
if size(Vn,1)~=size(Vn,2)
    error('V must be a square matrix.')
end
if size(Sn,2)~=1
    error('S must be a column vector.')
end
if size(Vn,1)~=size(Sn,1)
    error('The number of rows in S must equal the number of rows in V.')
end


% Initialization.
N = size(Sn,1);
beta = Vn'*Y;


% Tikhonov regularization
if (strncmp(method,'Tikh',4) || strncmp(method,'tikh',4))
    
    npoints = 200;         % Number of points on the curve.
    smin_ratio = 16*eps;   % Smallest regularization parameter.
    reg_param = zeros(npoints,1);
    G = zeros(npoints,1);
    s = sqrt(Sn);
    reg_param(npoints) = max([s(N),s(1)*smin_ratio]);
    ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
    for i=npoints-1:-1:1
        reg_param(i) = ratio*reg_param(i+1);
    end
    
    if show==1
        disp('Calculating GCV curve.')
    end
    % Vector of GCV-function values.
    for i=1:npoints
        G(i) = reggcvfun(reg_param(i),Sn,beta);
    end
    
    % Plot GCV function.
    if show==1
        loglog(reg_param,G,'-'), xlabel('\lambda'), ylabel('G(\lambda)')
        title('GCV function')
    end
    
    % Find minimum
    if show==1
        disp('Searching GCV minimum')
        OPT=optimset('Display','iter');
    else
        OPT=optimset('Display','off');
    end
    [minG,minGi] = min(G); % Initial guess.
    reg_min = fminbnd(@reggcvfun,...
        reg_param(min(minGi+1,npoints)),...
        reg_param(max(minGi-1,1)),OPT,Sn,beta); % Minimizer.
    minG = reggcvfun(reg_min,Sn,beta); % Minimum of GCV function.
    if show==1
        ax = axis;
        HoldState = ishold; hold on;
        loglog(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
        title(['GCV function, minimum at \lambda = ',num2str(reg_min)])
        axis(ax)
        if (~HoldState)
            hold off
        end
    end
    % Truncated SVD
elseif (strncmp(method,'tsvd',4) || strncmp(method,'TSVD',4))
    
    rho=zeros(1,N-1);
    G=zeros(1,N-1);
    rho(N-1) = sum(beta(N,:).^2);
    G(N-1) = rho(N-1);
    for k=N-2:-1:1
        rho(k) = rho(k+1) + sum(beta(k+1,:).^2);
        G(k) = rho(k)/((N - k)^2);
    end
    reg_param = (1:N-1)';
    
    % Plot GCV function.
    if show==1
        semilogy(reg_param,G,'o'), xlabel('k'), ylabel('G(k)')
        title('GCV function')
    end
    
    % Find minimum
    [minG,reg_min] = min(G);
    if show==1
        ax = axis;
        HoldState = ishold; hold on;
        semilogy(reg_min,minG,'*r',[reg_min,reg_min],[minG/1000,minG],':r')
        title(['GCV function, minimum at k = ',num2str(reg_min)])
        axis(ax);
        if (~HoldState)
            hold off
        end
    end
    
else
    error('Illegal method.')
end
end

function G=reggcvfun(lam,s2,beta)
% reggcvfun  Computes GCV function for reggcv.
%            Auxiliary function for reggcv.
%
% Written by Vincent Verdult, May 2004.
% Based on Regularization Tools by P. C. Hansen

f=lam^2./(s2+lam^2);
G=norm((f*ones(1,size(beta,2))).*beta,'fro')^2/(sum(f)^2);
end