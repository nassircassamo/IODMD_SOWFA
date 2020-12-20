function [K,gamma] = lobshinfsyn(A,C,Q,R,S,murange,c)
%LOBSHINFSYN H-infinity design of the observer gain for an LPV model:
%
%     x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + w(k)
%     y(k)   = C x(k) + D u(k) + v(k)
%
% The observer gives an estimate of the state x_est and an estimate of the
% output y_est, through:
%
%     x_est(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K kron(mu(k),e(k))
%     y_est(k)   = C x(k) + D u(k) + e(k)
%     e(k) = y(k) - y_est(k)
%
% [K,gamma] = lobshinfsyn(A,C,Q,R,S,murange) returns the Kalman gain K that 
% yields a stable observer resulting in minimal error x_est - x in an Hinf 
% sense for the scheduling in the range [mumin, mumax].
% A, C are the state-space matrices in the model structure above, Q,R,S are 
% the (estimates) of the noise covariance matrices:
%
%   E([v(k);w(k)][v(j)',w(j)']) = [R,S',S,Q] delta(k-j)>=0
%
% where delta(k) is the unit pulse. murange is an mx2 matrix of the form
% containing the minimum and maximum values of of the scheduling, i.e. 
% murange=[mumin,mumax].
% gamma is the hinf norm of the observer error system (transfer between 
% scaled noise and x_est-x) i.e. a smaller gamma indicates better 
% estimation of the states, given the noise properties Q,R,S.
%
% If c=1, the Kalman gain is kept scheduling independent, i.e. we have that
% K(:,l+1:s*l) = 0.
% 
% Uses YALMIP and an SDP solver (e.g. SEDUMI), which can be downloaded at:
%  http://users.isy.liu.se/johanl/yalmip/
%  http://sedumi.ie.lehigh.edu

%  Pieter Gebraad
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2011

% check if YALMIP is installed
if exist('yalmip','file') ~= 2
    error('lobshinfsyn:yalmipnotfound','The function LOBSHINFSYN uses YALMIP and an SDP solver (e.g. SEDUMI) to find optimal stabilizing LPV observer gains. \n YALMIP is not found in your MATLAB path. \n Please install YALMIP and an SDP solver and add it to your MATLAB path \n YALMIP is available at: http://users.isy.liu.se/johanl/yalmip \n A free SDP solver called SEDUMI is available at: http://sedumi.ie.lehigh.edu');
end

% check number of input arguments
if nargin<7
    c = 0;
end

% find required sizes
n = size(A,1);
s = size(murange,1);
l = size(C,1);

% check sizes
if ~isequal(size(Q),[n,n])
    error 'Incorrect size of Q'
end    
if ~isequal(size(R),[l,l])
    error 'Incorrect size of R'
end
if ~isequal(size(S),[n,l])
    error 'Incorrect size of S'
end
if ~isequal(size(A),[n,(s+1)*n])
    error 'Incorrect size of A'
end
if ~isequal(size(murange,2),2)
    error 'murange should contain two columns: murange=[mumin,mumax]'
end

% corner points of convex polytope around scheduling data
much = allcomb(murange);

% Solve the LMI's to find the Hinf observer stable for scheduling
% inside polytope
Q = sqrtm(Q - S*(R\(S'))); % = (Q_x)^(1/2)

yalmip clear;
options = sdpsettings('verbose',0);
P = sdpvar(n,n); gamma = sdpvar(1,1);
if c
    R = sqrtm(R);
    X = S/(R');
    M = sdpvar(n,l);
    LMI = set(P>zeros(n));
    for j = 1:(1/2)*size(much,2)
        Ap = A(1:n,1:n);
        Cp = C(1:l,1:n);
        for i = 1:s
            Ap = Ap + much(i,j)*A(:,(i-1)*n+1:i*n);
            Cp = Cp + much(i,j)*C(:,(i-1)*n+1:i*n);
        end
        LMI = LMI + set([P,            zeros(n,n),     P*Ap-M*Cp,  P*X-M*R,       P*Q;
                        zeros(n,n),    eye(n),         eye(n),     zeros(n,n+l);          
                        Ap'*P-Cp'*M',  eye(n),         P,          zeros(n,n+l);     
                        X'*P-R'*M',    zeros(l,2*n),               gamma*eye(l),  zeros(l,n);                 
                        Q'*P,          zeros(n,2*n+l),                            gamma*eye(n)]>0);
    end
    diagnostic = solvesdp(LMI,gamma,options);

    if diagnostic.problem ~= 0
        if diagnostic.problem == -3 || diagnostic.problem == -2
            error('lobshinfsyn:nosolver','The function LOBSHINFSYN uses YALMIP and an SDP solver (e.g. SEDUMI) to find optimal stabilizing LPV observer gains. \n No SDP solver is found in your MATLAB path. Please add an SDP solver to YALMIP. \n A free SDP solver called SEDUMI is available at: http://sedumi.ie.lehigh.edu')    
        else
            error('lobshinfsyn:sdpsolveerror',['An error/warning occured in solving the LMIs in LOBSHINFSYN: ',yalmiperror(diagnostic.problem)])    
        end
    end
    
    gamma = double(gamma);
    P = double(P);
    M = double(M);
    K = P\M;
    K = [K,zeros(n,s*l)];
else
    RC = R\C(:,1:n);
%    LMI = set(P-eye(n)>0);
    for j = 1:size(much,2)
        Ap = A(1:n,1:n);
        for i = 1:s
            Ap = Ap + much(i,j)*A(:,(i-1)*n+1:i*n);
        end
        Ap = Ap-S*RC;
        LMI = LMI + set([P,             zeros(n,n),     P*Ap,                   P*Q;
                         zeros(n,n),    eye(n),         eye(n),                 zeros(n,n);
                         Ap'*P,         eye(n),         P+gamma*C(:,1:n)'*RC,   zeros(n,n);
                         Q'*P,          zeros(n,2*n),                           gamma*eye(n)]>0);
    end
    diagnostic = solvesdp(LMI,gamma,options);
    
    if diagnostic.problem ~= 0
        if diagnostic.problem == -3 || diagnostic.problem == -2
            error('lobshinfsyn:nosolver','The function LOBSHINFSYN uses YALMIP and an SDP solver (e.g. SEDUMI) to find optimal stabilizing LPV observer gains. \n No SDP solver is found in your MATLAB path. Please add an SDP solver to YALMIP. \n A free SDP solver called SEDUMI is available at: http://sedumi.ie.lehigh.edu')    
        else
            warning('lobshinfsyn:sdpsolveerror',['A warning occured in solving the LMIs in LOBSHINFSYN: ',yalmiperror(diagnostic.problem)])    
        end
    end

    gamma = double(gamma);
    P = double(P);
    Cp = (P-eye(n))\(C(:,1:n)');
    for i = 1:s+1
      K(:,(i-1)*l+1:i*l) = (gamma*A(:,(i-1)*n+1:i*n)*Cp+S)/(gamma*C(:,1:n)*Cp+R);
    end
end
gamma = sqrt(gamma);

end
