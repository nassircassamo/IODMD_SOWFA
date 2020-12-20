function options = optim5to6(fopts) 
%OPTIM5TO6   Translates a foptions-vector to an optimset-structure 
% 
% Syntax 
%            options=optim5to6(fopts) 
% 
% Description 
%            This function translates a MATLAB 5 optimization options 
%            vector -- as generated using foptions -- into a MATLAB 6 
%            compatible optimset-structure. 
%            Translated fields are: 
% 
%                 Field    Description 
%                 -------------------------- 
%                    1     Display 
%                    2     TolX 
%                    3     TolFun 
%                    9     Jacobian 
%                   14     MaxFunEval 
% 
% Inputs 
%   fopts    A MATLAB 5 compatible foptions-vector 
% 
% Outputs 
%   options  A MATLAB 6 compatible optimset-structure 
% 
% Remarks 
%            MATLAB 5 uses a default parameter and function tolerance 
%            of 1e-4. This is indicated by the second and third element 
%            if fopts that are 1e-4 in the default case. 
%            MATLAB 6 uses a default value of 1e-6 for both tolerances, 
%            but setting the tolerances to 1e-6 if the fopts vector 
%            contains the default values is impossible: there is no 
%            way of telling whether the user just used the default values 
%            or that he actually specified 1e-4 as tolerance. 
%            Consequently, the tolerances are just copied, AND THERE WILL 
%            THUS BE DIFFERENT RESULTS IN AN OPTIMIZATION WHEN USING 
%            A DEFAULT FOPTIONS-VECTOR OR A DEFAULT OPTIMSET-STRUCTURE. 
% 
% Uses functions 
%            mkoptstruc 
% 
% See also 
%            foptions, optimset 
 
% Niek Bergboer, 2001 
% Revised by Ivo Houtzager, 2007
% Copyright (c) 2001-2007, Delft Center of Systems and Control 
 
if ndims(fopts)~=2 || ~all(size(fopts)==[1 18]),
    error('The size of fopts should be 1-by-18');
end

options = mkoptstruc;
% Field 1 determines Display verbosity
if fopts(1) == 0,
    options.Display = 'final';
elseif fopts(1) == 1,
    options.Display = 'iter';
else
    error('Illegal value in options(1)');
end

% Field 2 is TolX
options.TolX = fopts(2);

% Field 3 is TolFun
options.TolFun = fopts(3);

% Field 9 is Jacobian
if fopts(9) == 0,
    options.Jacobian = 'off';
elseif fopts(9) == 1,
    options.Jacobian = 'on';
else
    error('Illegal value in options(9)');
end

% Options 14 is MaxFunEval: leave default when 0
if fopts(14) ~= 0,
    options.MaxFunEval = fopts(14);
end
 
 
 
 

