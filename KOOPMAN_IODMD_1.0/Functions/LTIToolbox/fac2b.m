function B = fac2b(A,C,varargin) 
%FAC2B      Frequency domain estimation of B based on A, C 
%           and measured FRF. D will be assumed to be zero. 
%           General model structure: 
%                          -1 
%           H = C (w I - A)   B 
% 
% Syntax: 
%           B=fac2b(A,C,H,w) 
%           B=fac2b(A,C,H1,w1,...,Hp,wp) 
% 
% Inputs: 
% A,C       Known or estimated A and C matrices of an LTI state-space 
%           model. 
% H         Measured frequency response function (FRF) 
% w         Complex frequencies at which the FRF is measured 
%           Multiple data batches can be specified in the parameter 
%           list, but only for discrete-time models. 
% 
% Outputs: 
% B         The state space model's B matrix 
% 
% Uses functions: 
%     ltifrf 
% 
% See Also : fac2bd, dac2b 
 
% Based partly on code by Thomas McKelvey, 1995 
% Revised by Niek Bergboer, September 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 1995-2007, Delft Center of Systems and Control

% Check input arguments
if nargin < 4
    error('FAC2B requires at least four input arguments.')
end
if rem(nargin-2,2) ~= 0
    error('An integer number of data batches must be specified.');
end
Nbatch = (nargin-2)/2;
N = zeros(Nbatch,1);

if ndims(varargin{1}) ~= 3,
    error('H must be a 3D array.');
end
l = size(varargin{1},1);
m = size(varargin{1},2);
N(1) = size(varargin{1},3);
n = size(A,1);

if size(A,2) ~= n
    error('A must be square.');
end
if size(C,2) ~= n
    error('The size of C does not correspond to the number of states');
end
if size(C,1) ~= l
    error('The size of C does not correspond to the number of outputs')
end
if size(varargin{2},1) ~= N(1)
    error('The number of samples in H and w does not correspond.');
end
if size(varargin{2},2) ~= 1
    error('w must have 1 column.');
end
if N(1)*l<n+l+m,
    error('The number of samples is too small.');
end

for Batch = 2:Nbatch,
    if ndims(varargin{(2*(Batch-1))+1})~=3,
        error('H for batch %d must be a 3D-array.',Batch);
    end
    if size(varargin{(2*(Batch-1))+1},1)~=l,
        error('The number of outputs must be the same among all batches.');
    end
    if size(varargin{(2*(Batch-1))+1},2)~=m,
        error('The number of inputs must be the same among all batches.');
    end
    N(Batch) = size(varargin{(2*(Batch-1))+1},3);
end

if max(abs(real(varargin{2})))<=10*eps
    wscale = (max(imag(varargin{2}))+min(imag(varargin{2})))/2;
    if Nbatch>1,
        error('Concatenating data batches is not supported for continuous-time models.');
    end
    %timing = 'cont';
else
    wscale = 1;
    %timing = 'disc';
end

% Note: if Nbatch>1, the first batch was discrete-time
for i = 2:Nbatch
    if max(abs(real(varargin{2*(Batch-1)+2})))<=10*eps,
        error('Discrete-time and continuous-time batches cannot be mixed.');
    end
end

Rold = [];
for Batch = 1:Nbatch,
    H = varargin{2*(Batch-1)+1};
    w = varargin{2*(Batch-1)+2};
    Theta = zeros(N(Batch) * l,m);
    for i = 1:m
        Theta(:,i) = reshape(H(:,i,:),N(Batch) * l,1);
    end

    Phi = zeros(N(Batch) * l,n);
    Phi(:,1:n) = ltifrf((A/wscale)',(C/wscale)',eye(n),[],[],w/wscale,2).';

    M = zeros(size(Rold,1)+2*N(Batch)*l,size(Phi,2)+size(Theta,2));
    if ~isempty(Rold)
        M(1:size(Rold,1),:) = Rold;
    end
    M(size(Rold,1)+1:size(Rold,1)+N(Batch)*l,1:size(Phi,2)) = real(Phi);
    Phi = imag(Phi);
    M(size(Rold,1)+N(Batch)*l+1:size(Rold,1)+2*N(Batch)*l,1:size(Phi,2)) = Phi;
    M(size(Rold,1)+1:size(Rold,1)+N(Batch)*l,size(Phi,2)+1:size(Phi,2)+size(Theta,2)) = real(Theta);
    Theta = imag(Theta);
    M(size(Rold,1)+N(Batch)*l+1:size(Rold,1)+2*N(Batch)*l,size(Phi,2)+1:size(Phi,2)+size(Theta,2)) = Theta;
    clear Phi;
    clear Theta;

    R = qr(M,0);
    Rold = triu(R(1:size(R,2)-m,:));
    clear M;
end
R = triu(R(1:n,:));
B = R(:,1:n)\R(:,n+1:n+m);




