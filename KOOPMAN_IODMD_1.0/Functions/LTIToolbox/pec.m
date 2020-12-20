function p = pec(y,y_est)
%PEC         Compute the Prediction Error Cost (VAF) between two signals. 
% 
% Syntax: 
%            p = pec(y,y_estimate) 
% 
% Description: 
%            The PEC is calculated as: 
% 
%               p = 1/N ||y - y_est||^2_F 
% 
%            The PEC is often used to verify the 
%            correctness of a model, by comparing the real 
%            output with the estimated output of the model. 
% 
% Inputs: 
%  y         Signal 1, often the real output. 
%  y_est     Signal 2, often the estimated output of a model. 
% 
% Output: 
%  p         PEC, computed for the two signals 

% Ivo Houtzager, 2008
% Copyright (c) 1996-2008, Delft Center of Systems and Control 

if nargin < 2
    error('PEC requires two input arguments.');
end
if size(y,2) > size(y,1)
    y = y';
end
if size(y_est,2) > size(y_est,1)
    y_est = y_est';
end

N = size(y,1);
if size(y_est,1) ~= N
    error('Both signals should have an equal number of samples.');
end
if size(y,2) ~= size(y_est,2)
    error('Both signals should have an equal number of components.');
end

e = y - y_est;
p = norm(e'*e,'fro')./N;