function v = fit(y,y_est)
%FIT         Compute the best FIT percentage (FIT) between two signals. 
% 
% Syntax: 
%            v = fit(y,y_estimate) 
% 
% Description: 
%            The FIT is calculated as: 
% 
%                            ||y-y_est||^2_F 
%               v = ( 1 -  -------------------  ) * 100% 
%                          ||y - mean(y)||^2_F 
% 
%            The FIT of two signals that are the same is 
%            100%. If they differ, the  FIT will be lower. 
%            When y and y_est have multiple columns, the FIT 
%            is calculated for every column in y and y_est. 
%            The FIT is often used to verify the 
%            correctness of a model, by comparing the real 
%            output with the estimated output of the model. 
% 
% Inputs: 
%  y         Signal 1, often the real output. 
%  y_est     Signal 2, often the estimated output of a model. 
% 
% Output: 
%  v         FIT, computed for the two signals 

% Ivo Houtzager, 2008
% Copyright (c) 1996-2008, Delft Center of Systems and Control 

if nargin<2
  error('Not enough input variables')
end
if size(y,2)>size(y,1)
  y=y';
end
if size(y_est,2)>size(y_est,1)
  y_est=y_est';
end

N=size(y,1);

if ~(size(y_est,1)==N)
  error('Both signals should have same length')
end

v = zeros(size(y,2),1);
for i = 1:size(y,2)
    v(i) = max(diag(100*(1-norm(y(:,i)-y_est(:,i))./norm(y(:,i)-ones(size(y(:,i),1),1)*mean(y(:,i))))),0);
end





