function [B,D,R] = dac2bd(A,C,varargin)
%DAC2DB     Estimates the matrices B and D of the state space model
%
%                      x(k+1) = A x(k) + B u(k)
%                      y(k)   = C x(k) + D u(k) + v(k)
%           using the knowledge of the pair A, C. This function enables
%           the concatenation of different input-output data batches.
%
%           B and D are calculated by solving a linear least squares
%           problem. The influence of the initial state x0 in the data
%           batches is compensated for in the solution. After B and D
%           are found, x0 can be calculated for a single data batch
%           using the function dinit.
%
% Syntax:
%           [B,D]=dac2bd(A,C,u,y);
%           [B,D]=dac2bd(A,C,u1,y1,...,up,yp);
%
% Input:
%   A, C    The estimated (A,C) pair of the quadruple of system matrices
%   u,y     The input respectively output data of the system to be
%           identified. Multiple data batches (ui,yi) may be specified on
%           the command line.
%
% Output:
%   B, D    The estimated (B,D) pair.
%
% See also: dac2b dinit

% Based on dac2bd.m in the SMI-1.0 toolbox, by Michel Verhaegen
% Revised by Niek Bergboer, 2001
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1994-2007, Delft Center of Systems and Control 

% Check input arguments
if nargin < 4
    error('DAC2BD requires at least four input arguments.')
end
if rem(nargin-2,2) ~= 0
    error('An integer number of data batches must be specified.');
end

% Get number of batches
Nbatch = (nargin-2)/2;
N = zeros(Nbatch,1);

% Get dimensions
n = size(A,1);

% Get signal dimensions from the first batch
l = size(varargin{2},2);
m = size(varargin{1},2);
if l == 0
    error('DAC2BD requires an output')
end
if m == 0
    error('DAC2BD requires an input')
end
N(1) = size(varargin{2},1);

% Check consistency in the first batch.
if size(varargin{1},1) ~= N(1)
    error('Input and output should have same length.')
end

% Process the other batches and check for consistency
for Batch = 2:Nbatch,
    if size(varargin{2*(Batch-1)+1},2)~=m
        error('All batches should have the same number of inputs.');
    end
    if size(varargin{2*(Batch-1)+2},2)~=l
        error('All batches should have the same number of outputs.');
    end
    N(Batch) = size(varargin{2*(Batch-1)+2},1);
    if size(varargin{2*(Batch-1)+2},1)~=N(Batch),
        error('Input and output should have same length.')
    end
end

% Check input arguments
if size(A,2) ~= n
    error('A should be square.')
end
if size(C,2) ~= n
    error('C matrix must have the same number of colums as A.')
end
if (max(abs(eig(A)))>1)
    disp('A matrix has an instable pole. The estimate of B and D might be very bad.')
end
for Batch = 1:Nbatch,
    if (N(Batch)<n+n*m+m*l+1)
        error('Not enough data points in batch %d to find an estimate.',Batch);
    end
end;

if Nbatch == 1
    if exist('dac2bdc','file') == 3
        try
            [B,D] = dac2bdc(A,C,varargin{1},varargin{2});
            return
        catch
            warning('LTI:mexFailed','Loading of mex failed, use m-file instead.');
        end
    end
end

[Q,A1] = schur(A);
C1 = C*Q;

for Batch = 1:Nbatch,
    %Define signals
    u = varargin{2*(Batch-1)+1};
    y = varargin{2*(Batch-1)+2};
    Y = zeros(l*N(Batch),1);
    Y(:) = y;

    temp = zeros(n*l,N(Batch));
    Gamma = zeros(N(Batch)*l,n);
    Gamma(1:l,:) = C1;
    An = A1;
    for i = 1:floor(log(N(Batch))/log(2)),
        Gamma(2^(i-1)*l+1:2^i*l,:) = Gamma(1:2^(i-1)*l,:)*An;
        An = An*An;
    end
    Gamma(2^i*l+1:N(Batch)*l,:) = Gamma(1:N(Batch)*l-2^i*l,:) * An;
    temp(:) = Gamma';
    for j = 1:l
        Gamma(N(Batch)*(j-1)+1:N(Batch)*j,:) = temp(n*(j-1)+1:n*j,:)';
    end

    Yij = zeros(N(Batch)*l,n*m);
    e = eye(n);
    temp  =  zeros(N(Batch)*l,1);
    for i = 1:n
        if i==n
            s = i;
        elseif A(i+1,i)==0
            s = i;
        else
            s = i+1;
        end
        Aout = A1(1:s,1:s);
        Bout = e(1:s,i);
        Cout = C1(:,1:s);
        for j = 1:m
            x  = ltiitr(Aout,Bout,u(:,j),[],zeros(length(Aout),1));
            yij =  x*Cout.';
            temp(:) = yij;
            Yij(:,(j-1)*n+i) = temp;
        end
    end

    Uij = zeros(N(Batch) * l,m * l);
    for j = 1:l,
        Uij((j-1) * N(Batch)+1:j * N(Batch),(j-1) * m+1:j * m) = u;
    end

    % Perform initial QR
    PhiTh = [Gamma Yij Uij Y];
    R = triu(qr(PhiTh));

    % Batch concatenation
    if Batch == 1,
        R = R(n+1:n+n*m+m*l,n+1:n+n*m+m*l+1);
        Rold = R;
    else
        R = triu(qr([R(n+1:n+n*m+m*l,n+1:n+n*m+m*l+1); Rold]));
        R = R(1:n*m+m*l,:);
        Rold = R;
    end
end

BD = R(1:n*m+m*l,1:n*m+m*l)\R(1:n*m+m*l,n*m+m*l+1);
B = zeros(n,m);B(:) = BD(1:n * m);
B = Q*B;
DT = zeros(m,l);DT(:) = BD(n * m+1:n * m+l * m);
D = DT';
