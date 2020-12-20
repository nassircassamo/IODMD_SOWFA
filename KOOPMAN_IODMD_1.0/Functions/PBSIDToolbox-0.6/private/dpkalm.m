function [K,S,E] = dpkalm(A,C,Q,R,S,F,N)
%DPKALM Dynamic programming of Kalman design for a discrete state-space system.
%   [K,P,E] = DPKALM(A,C,Q,R,S,F,N) calculates the optimal gain matrices
%   K[n] such that:
%
%   For a discrete-time state-space model SYS, u[n] = -K[n]x[n] minimizes
%
%     J = 0.5*Sum {x[n]'Qx[n] + u[n]'Ru[n] + x[n]'Su[n]} + 0.5*x[-1]'Fx[-1]
%
%   subject to  x[n+1] = Ax[n] + Bu[n] with n = 0...N-1.
%
%   Also returned are the the solution S of the associated algebraic
%   Riccati equation and the closed-loop eigenvalues E = EIG(A-K*C).

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number of input arguments
if nargin ~=7
    error('DPKALM requires at least five or six input arguments')
end
if ndims(A) > 2;
    array = 2;
else
    array = 1;
end

% check dimensions and symmetry
[nax ns nna] = size(A);
[ncx nc nnc] = size(C);
[nqx nq nnq] = size(Q);
[nfx nf nnf] = size(F);
[nrx nr nnr] = size(R);
if ~isequal(nna,nnc,nnq,nnr,N);
    if array == 2
        error('The number of arrays have to be equal to N.');
    end
end
if ~isequal(nc,ns);
    error('The A and C matrices must have the same number of columns.')
end
if ~isequal(nqx,nax,nfx,ns,nq,nf);
    error('The A, Q and F matrices must be the same size.')
end
if ~isequal(nrx,nr,ncx);
    error('The R matrix must be square with as many rows as C.')
end
if ~isreal(Q) || ~isreal(R)
   error('The weight matrices Q, R must be real valued.')
end

% backwards iteration
K = zeros(nb,ns,N);
E = zeros(ns,N);
X = zeros(ns,ns,N+1);
X(:,:,1) = F;
switch array
    case 1
        for t = 1:N
            K(:,:,t) = (R + C*X(:,:,t)*C')\(C*X(:,:,t)*A' + S');
            X(:,:,t+1) = (A - K(:,:,t)*C(:,:,t))*X(:,:,t)*(A - K(:,:,t)*C(:,:,t))' + K(:,:,t)'*R*K(:,:,t) + Q;
            if nargout >= 2
                E(:,t) = eig(A - K(:,:,t)*C);
            end
        end
    case 2
        for t = 1:N
            K(:,:,t) = (R(:,:,t) + C(:,:,t)*X(:,:,t)*C(:,:,t)')\(C(:,:,t)*X(:,:,t+1)*A(:,:,t)' + S(:,:,t)');
            X(:,:,t+1) = (A(:,:,t) - K(:,:,t)*C(:,:,t))*X(:,:,t)*(A(:,:,t) - K(:,:,t)*C(:,:,t))' + K(:,:,t)'*R(:,:,t)*K(:,:,t) + Q(:,:,t);
            if nargout >= 2
                E(:,t) = eig(A(:,:,t) - K(:,:,t)*C(:,:,t));
            end
        end
end

