function [Ak,Bk,Ck,Dk,Kk,err,eigA,regA] = rpbsid(u,y,f,p,n,W,idopts,rlsopts,A,B,C,D,K,s)
%RPBSID Recursive Predictor-based Subspace IDentification 
%  [A,B,C,D,K]=rpbsid(u,y,f,p,n) reursively estimates the matrices A, B, C,
%  D and K of the state space model:
%
%       x(N) = A x(N-1) + B u(N-1) + K e(N-1)
%       y(N-1) = C x(N-1) + D u(N-1) + e(N-1)
%
%  where N is the number of observations. The input matrix u and output
%  matrix y must have the same number of observations but can have
%  different numbers of variables. The past and future window size p and f
%  must be higher then the expected order n. 
%
%  [A,B,C,D,K]=rpbsid(u,y,f,p,n,S) specifies the n times f*l permutation
%  matrix S. The default is S=[eye(n) zeros(n,n-f)].
%
%  [A,B,C,D,K]=rpbsid(u,y,f,p,n,S,idopts) specifies the identification
%  options. The default is idopts =
%  struct('method','varx','weight',0,'ltv',0,'noD',0,'past',0,'Kalm',0);
%
%  [A,B,C,D,K]=rpbsid(u,y,f,p,n,S,idopts,rlsopts) specified the recursive
%  least squares options. The default is rlsopts = struct('lambda',[0.999
%  0.999 0.999],'ireg',[1e-6 1e-6 1e6],'reg',0);
%
%  [A,B,C,D,K]=rpbsid(u,y,f,p,n,S,idopts,rlsopts,S,A0,B0,C0,D0,K0,Ts)
%  specifies the initial state-space matrices and sampling time.
%
%  [A,B,C,D,K,err]=rpbsid(u,y,f,p,n,S) resturns the prediction error of the
%  recursive least squares solvers.
%
%  [A,B,C,D,K,err,eigA]=rpbsid(u,y,f,p,n,S) returns the eigenvalue and
%  damping vectors over time.
%
%  [A,B,C,D,K,err,eigA,regA]=rpbsid(u,y,f,p,n,S) returns the regularisation 
%  over time.
%
%  [A,B,C,D,K,err,eigA,regA]=rpbsid(u,y,f,p,n,S,idopts,rlsopts,A0,B0,C0,D0,K0,Ts)
%  specifies the sampling of returned vector err, eigA. The default is Ts=0
%  (off).
%
%  References:
%  [1] Ali H. Sayed, "Adaptive Filters", Wiley and Sons, 2008

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  The Netherlands, 2010

% check number if input arguments
if nargin < 5
    error('RPBSID requires at least five input arguments.')
end

% Determine sizes
if size(u,1) > size(u,2)
    u = u';
end
if size(y,1) > size(y,2)
    y = y';
end
r = size(u,1);
l = size(y,1);
N = size(y,2);

% Assign known values to the parameters
if nargin < 14 || isempty(s)
    s = 0;
end
if nargin < 13 || isempty(K)
    K = zeros(n,l);
end
if nargin < 12 || isempty(D)
    D = zeros(l,r);
end
if nargin < 11 || isempty(C)
    C = zeros(l,n);
end
if nargin < 10 || isempty(B)
    B = zeros(n,r);
end
if nargin < 9 || isempty(A)
    A = zeros(n,n);
end
if nargin < 8 || isempty(idopts)
   rlsopts = struct('ireg',[1e-6 1e-6 1e-6],'lambda',[0.999 0.999 0.999],'reg',0);
end
if nargin < 7 || isempty(idopts)
   idopts = struct('method','varx','weight',0,'ltv',0,'noD',0,'past',0,'Kalm',0);
end
if nargin < 6 || isempty(W)
   W = [eye(n) zeros(n,f*l-n)];
end
switch lower(idopts.method)
    case 'fir'
        m = r;
    case 'varx'
        m = r+l;
    case 'varmax'
        m = r+2*l;
    otherwise
        error('Unknown type.')
end
       
% Initialisation of recursive least squares iterations

% Estimation of the Markov parameters
Plk = rlsopts.ireg(1);                  % Initial inverse of sample covariance matrix
CKD = zeros(l,p*m+~idopts.noD*r);       % Initial solution
CKDS = zeros(f*l,p*m+~idopts.noD*r);    % Initial solutions
P = zeros((p+f-1)*m,1);                 % Initial regression vector
reg_min = rlsopts.reg(1);

% Estimation of the output matrices
Pcd = rlsopts.ireg(2);
if idopts.noD
    CD = C;
else
    CD = [C D];
end

% Estimation of the state matrices
Pabk = rlsopts.ireg(3);
if strcmpi(idopts.method,'fir')
    ABK = [A B];
else
    ABK = [A B K];
end

% Initialisation of forward Ricatti iterations (if selected)
if idopts.Kalm
    Px = rlsopts.ireg(3).*eye(n);
    Q = rlsopts.ireg(3).*eye(n);
    R = rlsopts.ireg(3).*eye(l);
    S = zeros(n,l);
end

% Allocate state-space matrices for return
if s > 0
    Ak = zeros(n,n,floor(N/s));
    Bk = zeros(n,r,floor(N/s));
    Ck = zeros(l,n,floor(N/s));
    Dk = zeros(l,r,floor(N/s));
    Ak(:,:,1) = A;
    Bk(:,:,1) = B;
    Ck(:,:,1) = C;
    Dk(:,:,1) = D;
    if strcmpi(idopts.method,'varx') || strcmpi(idopts.method,'varmax')
        Kk = zeros(n,l,floor(N/s));
        Kk(:,:,1) = K;
    end
end
if nargout > 5
    err = zeros(2*l+n,N);
end
if nargout > 6
    eigA = zeros(n,N);
end
if nargout > 7
    regA = zeros(1,N);
end
if idopts.past
    if norm(W) > 1;
        Pw = 1/rlsopts.ireg(1);
    else
        Pw = rlsopts.ireg(1);
    end
    W = W';
end

% Store vectors for next iteration
start = 2;
U1 = u(:,start-1);
Y1 = y(:,start-1);
E1 = zeros(l,1);
Xf1 = zeros(n,1);

% Start recursive identification
h = 1;
startA = 3;
for k = start:1:N   
    % New signal vector
    switch lower(idopts.method)
        case 'fir'
            P = [P(m+1:end,:); U1];
        case 'varx'
            P = [P(m+1:end,:); U1; Y1];
        case 'varmax'
            P = [P(m+1:end,:); U1; Y1; E1];
        otherwise
            error('Unknown type.')
    end
    Y = y(:,k);
    U = u(:,k);
    if idopts.noD
        Z = P((f-1)*m+1:end,:);
    else
        Z = [P((f-1)*m+1:end,:); U];
    end

    if k >= p   
        % Solve Regression problem recursively
        if reg_min ~= 0
            [CKD,Plk] = rls_ew_track_reg(Z,Y,CKD,Plk,rlsopts.lambda(1),reg_min);
        else
            [CKD,Plk] = rls_ew_track(Z,Y,CKD,Plk,rlsopts.lambda(1));
        end
        if nargout > 5
            err(1:l,k-f+1) = Y - CKD*Z;
        end
        if nargout > 7
            regA(:,k) = reg_min;
        end
        CKDS(1:(f-1)*l,:) = CKDS(l+1:f*l,:);
        CKDS((f-1)*l+1:f*l,:) = CKD;
    end
    
    if k >= startA*p
        % Construction of observability times controllability
        LK = zeros(l*f,m*p);
        if idopts.ltv
            if idopts.weight
                for i = 0:f-1
                    LK(i*l+1:(i+1)*l,i*m+1:p*m) = CKDS(i*l+1:(i+1)*l,1:(p-i)*m);
                    if i ~= 0
                        for j = 0:i-1
                            LK(i*l+1:(i+1)*l,:) = LK(i*l+1:(i+1)*l,:) + CKDS(i*l+1:(i+1)*l,(p-i+j)*m+r+(1:l))*LK(j*l+1:(j+1)*l,:);
                        end
                    end
                end
            else
                for i = 1:f
                    LK((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = CKDS((i-1)*l+1:i*l,1:(p-i+1)*m);
                end
            end
        else
            if idopts.weight
                for i = 0:f-1
                    LK(i*l+1:(i+1)*l,i*m+1:p*m) = CKD(:,1:(p-i)*m);
                    if i ~= 0
                        for j = 0:i-1
                            LK(i*l+1:(i+1)*l,:) = LK(i*l+1:(i+1)*l,:) + CKD(:,(p-i+j)*m+r+(1:l))*LK(j*l+1:(j+1)*l,:);
                        end
                    end
                end
            else
                for i = 1:f
                    LK((i-1)*l+1:i*l,p*m-(p-i+1)*m+1:p*m) = CKD(:,1:(p-i+1)*m);
                end
            end
        end

        % Predicte future signal vector (= state estimate X)
        if idopts.ltv
            Xf = LK*P(1:p*m,:);
            Uf1 = P((p-1)*m+1:(p-1)*m+r,1);
            Yf1 = P((p-1)*m+r+1:(p-1)*m+r+l,1);
        else
            Xf = LK*P((f-1)*m+1:(p+f-1)*m,:);
            Uf1 = P((p+f-2)*m+1:(p+f-2)*m+r,1);
            Yf1 = P((p+f-2)*m+r+1:(p+f-2)*m+r+l,1);
        end
        if idopts.past
            [W,Pw] = rls_ew_track(W'*Xf,Xf,W,Pw,rlsopts.lambda(1));
            Xf = W'*Xf;
        else
            Xf = W*Xf;
        end
    end
    
    if k >= startA*p+1
        % The estimation of the system matrices 
        if idopts.noD
            [CD,Pcd] = rls_ew_track(Xf1,Yf1,CD,Pcd,rlsopts.lambda(2));
            if nargout > 5
                err(l+1:2*l,k) = Yf1 - CD*Xf1;
            end
            Ef1 = Yf1 - CD*Xf1;
        else
            [CD,Pcd] = rls_ew_track([Xf1; Uf1],Yf1,CD,Pcd,rlsopts.lambda(2));
            if nargout > 5
                err(l+1:2*l,k-f+1) = Yf1 - CD*[Xf1; Uf1];
            end
            Ef1 = Yf1 - CD*[Xf1; Uf1];
        end
        if strcmpi(idopts.method,'fir')
            [ABK,Pabk] = rls_ew_track([Xf1; Uf1],Xf,ABK,Pabk,rlsopts.lambda(3));
            if nargout > 5
                err(2*l+1:2*l+n,k) = Xf - ABK*[Xf1; Uf1];
            end
            if nargout > 6
                eigA(:,k) = sort(eig(ABK(:,1:n)));
            end
        else
            [ABK,Pabk] = rls_ew_track([Xf1; Uf1; Ef1],Xf,ABK,Pabk,rlsopts.lambda(3));
            if nargout > 5
                err(2*l+1:2*l+n,k) = Xf - ABK*[Xf1; Uf1; Ef1];
            end
            if nargout > 6
                eigA(:,k) = sort(eig(ABK(:,1:n)));
            end
        end  
           
        % Estimate stable Kalman gain by the forward Riccati iteration
        if idopts.Kalm
            VW = [Xf; Yf1] - [ABK(:,1:n+r); CD zeros(l,idopts.noD*r)]*[Xf1; Uf1];
            VW = VW*VW';
            Q = 0.5.*(VW(1:n,1:n) + rlsopts.lambda(3).*Q);
            R = 0.5.*(VW(n+1:n+l,n+1:n+l) + rlsopts.lambda(3).*R);
            S = 0.5.*(VW(1:n,n+1:n+l) + rlsopts.lambda(3).*S);
            K = (ABK(:,1:n)*Px*CD(:,1:n)' + S)/(R + CD(:,1:n)*Px*CD(:,1:n)');
            Px = ABK(:,1:n)*Px*ABK(:,1:n)' + Q - ABK(:,(n+r+1):(n+r+l))*(ABK(:,1:n)*Px*CD(:,1:n)' + S)';
        end
    end
       
    % Store state-space matrices (if selected)
    if k-f == h*s && s ~= 0
        if k >= startA*p+1
            Ak(:,:,h) = ABK(:,1:n);
            Ck(:,:,h) = CD(:,1:n);
            if strcmpi(idopts.method,'fir')
                Bk(:,:,h) = ABK(:,n+1:n+r);
                if idopts.noD
                    Dk(:,:,h) = zeros(l,r);
                else
                    Dk(:,:,h) = CD(:,n+1:n+r);
                end
            else
                if idopts.Kalm
                    Bk(:,:,h) = [ABK(:,n+1:n+r) K];
                else
                    Bk(:,:,h) = ABK(:,n+1:end);
                end
                if idopts.noD
                    Dk(:,:,h) = [zeros(l,r) eye(l)];
                else
                    Dk(:,:,h) = [CD(:,n+1:n+r) eye(l)];
                end
            end
        end
        h = h + 1;
    end
    
    % Store vectors and matrices for next iteration
    if k >= startA*p
        Xf1 = Xf;
    end
    if k < startA*p+1
        E = E1;
    end
    Y1 = Y;
    U1 = U;
    E1 = E;
end

if nargout >= 1 && s == 0
    % Store state-space matrices
    Ak = ABK(:,1:n);
    Ck = CD(:,1:n);
    if strcmpi(idopts.method,'fir')
        Bk = ABK(:,n+1:n+r);
        if idopts.noD
            Dk = zeros(l,r);
        else
            Dk = CD(:,n+1:n+r);
        end
    else
        Bk = ABK(:,n+1:n+r);
        if idopts.Kalm
            Kk = K;
        else
            Kk = ABK(:,n+r+1:end);
        end
        if idopts.noD
            Dk = zeros(l,r);
        else
            Dk = CD(:,n+1:n+r);
        end
    end
end

end % end of function RPBSID


function [theta,P] = rls_ew_track(z,y,theta,P,lambda)
%RLS_EW_TRACK Exponentially Weighted RLS iteration
%  [THETA,P]=RLS_EW_TRACK(Z,Y,THETA,P,LAMBDA) applies one iteration of
%  exponentially weighted regularized least-squares problem. In recursive
%  least-squares, we deal with the issue of an increasing amount of date Z
%  and Y. At each iteration, THETA is the solution. The scalar LAMBDA is
%  called the forgetting factor since past data are exponentially weighted
%  less heavily than more recent data.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  The Netherlands, 2010

% Assign default values to unspecified parameters
mz = size(z,1);
if (nargin < 5) || isempty(lambda)
    lambda = 1;
end
if (nargin < 4) || isempty(P)
    P = zeros(mz);
elseif isscalar(P)
    P = (1/P).*eye(mz); 
end
if (nargin < 3) || isempty(theta)
    theta = zeros(size(y,1),mz);
end

% Exponentially-Weighted, Tracking, Regularized, Least-Squares Iteration
P = (1/lambda).*(P - P*(z*z')*P./(lambda + z'*P*z));
P = 0.5.*(P+P'); % force symmetric
e = y - theta*z;
theta = theta + e*z'*P;

end % end of function RLS_EW_TRACK


function [theta,P] = rls_ew_track_reg(z,y,theta,P,lambda,reg_min)
%RLS_EW_TRACK_REG Exponentially Weighted and Regularized RLS iteration
%  [THETA,P]=RLS_EW_TRACK_REG(Z,Y,THETA,P,LAMBDA,REG) applies one iteration
%  of exponentially weighted regularized least-squares problem. In
%  recursive least-squares, we deal with the issue of an inceasing amount
%  of date Z and Y. At each iteration, THETA is the solution. The scalar
%  LAMBDA is called the forgetting factor since past data are exponentially
%  weighted less heavily than more recent data.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  The Netherlands, 2010

% Assign default values to unspecified parameters
mz = size(z,1);
if (nargin < 6) || isempty(reg_min)
    reg_min = 0;
end
if (nargin < 5) || isempty(lambda)
    lambda = 1;
end
if (nargin < 4) || isempty(P)
    P = zeros(mz);
elseif isscalar(P)
    P = (1/P).*eye(mz); 
end
if (nargin < 3) || isempty(theta)
    theta = zeros(size(y,1),mz);
end

% Exponentially-Weighted, Tracking, Regularized, Least-Squares Iteration
P = (1/lambda).*(P - P*(z*z')*P./(lambda + z'*P*z));
P = 0.5.*(P+P'); % force symmetric
if isscalar(reg_min)
    opts.SYM = true;
    opts.POSDEF = true;
    P1 = linsolve((eye(size(P)) + reg_min^2.*P),P,opts);
    e = y - theta*z;
    theta = theta + e*z'*P1;
elseif strcmpi(reg_min,'tikh')
    [U,S,V] = svd(pinv(P));
    s = diag(S);
    YP = (y-theta*z)';
    if isscalar(opt)
        reg_min = opt;
    elseif strcmpi(opt,'lcurve')
        reg_min = reglcurve(YP,U,s);
    elseif strcmpi(opt,'gcv')
        reg_min = reggcv(YP,U,s);
    end
    theta = theta + (V*(diag(s./(s.^2 + reg_min^2)))*U'*YP)';
elseif strcmpi(reg_min,'tsvd')
    [U,S,V] = svd(pinv(P));
    s = diag(S);
    YP = (y-theta*z)';
    if isscalar(opt)
        k_min = opt;
    elseif strcmpi(opt,'lcurve')
        k_min = reglcurve(YP,U,s,'tsvd');
    elseif strcmpi(opt,'gcv')
        k_min = reggcv(YP,U,s,'tsvd');
    end
    theta = theta + (V(:,1:k_min)*diag(1./s(1:k_min))*U(:,1:k_min)'*YP)';
end
end % end of function RLS_EW_TRACK

