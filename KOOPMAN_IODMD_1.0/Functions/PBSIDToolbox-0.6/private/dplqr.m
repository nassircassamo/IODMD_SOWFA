function [K,S,E] = dplqr(A,B,Q,R,F,N)
%DPLQR Dynamic programming of LQR design for a discrete state-space system.
%   [K,S,E] = DPLQR(A,B,Q,R,F,N) calculates the optimal gain matrices K[n]
%   such that:
%
%   For a discrete-time state-space model SYS, u[n] = -K[n]x[n] minimizes
%
%     J = 0.5*Sum {x[n]'Qx[n] + u[n]'Ru[n]} + 0.5*x[N]'Fx[N]
%
%   subject to  x[n+1] = Ax[n] + Bu[n] with n = 0...N-1.
%
%   Also returned are the the solution S of the associated algebraic
%   Riccati equation and the closed-loop eigenvalues E = EIG(A-B*K).

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number of input arguments
if nargin < 5
    error('DPLQR requires at least five or six input arguments')
end
if ndims(A) > 2;
    array = 2;
else
    array = 1;
end

% check dimensions and symmetry
[nax ns nna] = size(A);
[nbx nb nnb] = size(B);
[nqx nq nnq] = size(Q);
[nfx nf nnf] = size(F);
[nrx nr nnr] = size(R);
if ~isequal(nna,nnb,nnq,nnr,N);
    if array == 2
        error('The number of arrays have to be equal to N.');
    end
end
if ~isequal(nbx,nax);
    error('The A and B matrices must have the same number of rows.')
end
if ~isequal(nqx,nax,nfx,ns,nq,nf);
    error('The A, Q and F matrices must be the same size.')
end
if ~isequal(nrx,nr,nb);
    error('The R matrix must be square with as many columns as B.')
end
if ~isreal(Q) || ~isreal(R)
   error('The weight matrices Q, R must be real valued.')
end

% backwards iteration
K = zeros(nb,ns,N);
E = zeros(ns,N);
S = zeros(ns,ns,N);
S(:,:,N+1) = F;
switch array
    case 1
        for t = N:-1:1
            K(:,:,t) = (R + B'*S(:,:,t+1)*B)\B'*S(:,:,t+1)*A;
            S(:,:,t) = (A - B*K(:,:,t))'*S(:,:,t+1)*(A - B*K(:,:,t)) + K(:,:,t)'*R*K(:,:,t) + Q;
            if nargout >= 2
                E(:,t) = eig(A - B*K(:,:,t));
            end
        end
    case 2
        for t = N:-1:1
            K(:,:,t) = (R(:,:,t) + B(:,:,t)'*S(:,:,t+1)*B(:,:,t))\B(:,:,t)'*S(:,:,t+1)*A(:,:,t);
            S(:,:,t) = (A(:,:,t) - B(:,:,t)*K(:,:,t))'*S(:,:,t+1)*(A(:,:,t) - B(:,:,t)*K(:,:,t)) + K(:,:,t)'*R(:,:,t)*K(:,:,t) + Q(:,:,t);
            if nargout >= 2
                E(:,t) = eig(A(:,:,t) - B(:,:,t)*K(:,:,t));
            end
        end
end

