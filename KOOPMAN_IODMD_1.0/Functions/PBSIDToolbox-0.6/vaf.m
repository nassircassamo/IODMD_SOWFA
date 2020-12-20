function v = vaf(y,y_est) 
%VAF Variance Accounted For [%]
%  V = VAF(Y,Yest) computes the percentage of the Variance Accounted For
%  (VAF) between the two signals Y and Yest. The VAF is calculated as: 
% 
%                             variance(y-y_est) 
%               v = max( 1 -  ----------------- , 0 ) * 100% 
%                                variance(y) 
%
%  The VAF of two signals that are the same is 100%. If they differ, the
%  VAF will be lower. When Y and Yest have multiple columns, the VAF is
%  calculated for every column in Y and Yest. The VAF is often used to
%  verify the correctness of a stable model, by comparing the real output
%  with the estimated output of the stable model.
%
%  See also PEC, FIT.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check input arguments
if nargin < 2
    error('VAF requires two input arguments!');
end

% check dimensions of inputs
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

% calculate the variance accounted for
v = max(diag(100*(eye(size(y,2))-cov(y-y_est)./cov(y))),0);
 
 
 

