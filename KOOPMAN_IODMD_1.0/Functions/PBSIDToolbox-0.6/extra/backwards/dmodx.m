function X = dmodx(X,n)
%DMODX  Closed-loop LTI system identification using the PBSIDopt method.
%  x=dmodx(X,n) estimates the state sequence x of the identifiable system
%  with order n. The order n can be determined from the singular values
%  given by dordvarx or dordfir. The matrix X is also calculated by
%  dordvarx or dordfir.
%
%  See also: dordfir, dordvarx.
%
%  References:
%    [1] A. Chiuso, G. Picci, ``Consistency Analysis of Certain Closed-Loop
%    Subspace Identification Methods'', Automatica, Special Issue on System
%    Identification,  41(3), pp.377--391, 2005.
%
%    [2] A. Chiuso, ``The role of vector auto regressive modeling in
%    predictor based subspace identification'', Automatica, 43(6), 
%    pp.1034–-1048, 2007.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number of arguments
if nargin < 2
    error('DMODX requires at least two input arguments.');
end

% check for batches
if iscell(X)
    batch = length(X);
else
    batch = 1;
end

% do for all batches
for k = 1:batch
    if batch == 1
        x = X;
    else
        x = X{k};
    end
    
    % check the dimensions of the inputs
    mx = size(x,1);
    if (n < 1) || isempty(n)
        error('System order of zero or lower does not make sense!')
    end
    if mx < n
        error('The number of rows of matrix X must be equal or higher then the order n.')
    end
    
    % use only n rows
    x = x(1:n,:);
    
    if batch == 1
        X = x;
    else
        X{k} = x;
    end
end






