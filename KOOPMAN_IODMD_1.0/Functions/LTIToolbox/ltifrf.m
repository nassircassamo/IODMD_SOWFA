function H = ltifrf(A,B,C,D,dA,w,outopt) 
%LTIFRF    Linear Time Invariant Frequency Response Function 
% 
% Syntax: 
%          H = ltifrf(A,B,C,[],[],w,outopt) 
%          H = ltifrf(A,B,C,D,[],w,outopt) 
%          H = ltifrf([],[],[],D,[],w,outopt) 
%          H = ltifrf(A,B,C,[],dA,w,outopt) 
% 
% Description: 
%          ltifrf returns the Frequency Response Function (FRF) 
%          of a linear time-invariant state-space model, evaluated 
%          at the complex frequencies provided in w: 
% 
%                                    -1 
%                   H = C (w I  -  A)   B 
%                             n 
% 
% Input: 
%  A       State-space model matrix A 
%  B       State-space model matrix B 
%  C       State-space model matrix C 
%  D       (optional) State-space model matrix D 
%  dA      (optional) Calculates the change in FRF 
%          given the deviation dA in A. D and dA 
%          are mutually exclusive. 
%  w       Vector of complex frequencies. exp(j*omega) 
%          for discrete-time systems and j*omega for 
%          continuous-time systems. 
%  outopt  Controls how H will be returned (see below) 
% 
% Outputs: 
%  H       The FRF. Usually a 3D-array of size 
%          (#outputs x #inputs x #frequencies). 
% 
%          However, if outopt is non-empty and 1, 
%          H will be a vector of size 
%          (#outputs #inputs #frequencies x 1). 
%          If outopt is 2, H will be a matrix of size 
%          (#outputs x #inputs #frequencies). 
% 
 
% Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 
 
if nargin < 7
    error('LTIFRF requires seven input arguments');
end
if isempty(outopt)
    outopt = 0;
end

% Get dimensions
n = size(A,1);
m = size(B,2);
l = size(C,1);
N = size(w,1);

% Check input arguments: see whether n is compatible among A,B,C and dA
if size(A,2) ~= n
    error('A matrix must be square.');
end
if size(B,1) ~= n
    error('B matrix has wrong number of rows.');
end
if size(C,2) ~= n
    error('C matrix has wrong number of columns.');
end
if ~isempty(dA)
    if size(dA,1)~=n || size(dA,2)~=n
        error('dA and A matrix must be the same size.');
    end
end

if n > 0
    % n>0: case 1 or 2 (see MWEB documentation)
    if ~isempty(D),
        % Case 1: check dimensions of D
        if size(D,1)~=l,
            error('D matrix has wrong number of rows.');
        end
        if size(D,2)~=m,
            error('D matrix has wrong number of columns.');
        end
    else
        % Case 2: set D empty
        D = [];
    end
else
    % n==0: case 3
    A = [];
    B = [];
    C = [];
    if ~isempty(D),
        % Case 3: Since A,B and C are empty, set l and m from D
        l = size(D,1);
        m = size(D,2);
    end
end

VersionString = version;
matlabver = str2double(VersionString(1,1));
if isempty(dA),
    % Cases 1, 4 or 5: dA is empty
    if ~isempty(A),
        % Cases 1 or 4: call ltifrN
        if matlabver > 5,
            H = mimofr(A,B,C,[],w);
        else
            H = zeros(l,m,N);
            for i = 1:m
                H(:,i,:) = reshape(C*ltifr(A,B(:,i),w),l,1,N);
            end
        end
        if ~isempty(D),
            % Add the D-component
            for i = 1:l
                for j = 1:m
                    H(i,j,:) = H(i,j,:) + D(i,j)*ones(1,1,N);
                end
            end
        end
    else
        % Case 5: call ltifrD (static gain)
        H = zeros(l,m,N);

        % Every plane of the 3D-array is D
        for i = 1:l
            for j = 1:m
                H(i,j,:) = D(i,j)*ones(1,1,N);
            end
        end
    end
else
    % Cases 2 or 3: dA is non-empty
    if ~isempty(D)
        warning('LTI:ignoreD','D is ignored for the two-component dA calculation');
    end;
    H = zeros(l,m,N);

    % Calculate two FRFs and multiply them
    if matlabver > 5,
        H1 = mimofr(A,dA,C,[],w);
        H2 = mimofr(A,B,[],[],w);
    else
        H1 = zeros(l,n,N);
        H2 = zeros(n,m,N);
        for i = 1:n
            H1(:,i,:) = reshape(C*ltifr(A,dA(:,i),w),l,1,N);
        end
        for i = 1:m
            H2(:,i,:) = reshape(ltifr(A,B(:,i),w),n,1,N);
        end
    end

    % Multiply them
    for i = 1:N
        H(:,:,i) = H1(:,:,i)*H2(:,:,i);
    end
end

% Reshape output if so desired
if outopt == 1
    H = H(:);
elseif outopt == 2
    H = reshape(H,l,m*N);
end
 
 
 
 
 

