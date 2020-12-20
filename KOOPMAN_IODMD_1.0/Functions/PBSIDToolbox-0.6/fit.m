function v = fit(y,y_est)
%FIT Model fit [%]
%  V = FIT(Y,Yest) computes the percentage of the model fit (FIT) between
%  the two signals Y and Yest. The FIT is calculated as: 
% 
%                               ||y-y_est||_2 
%               v = max( 1 -  ----------------- , 0 ) * 100% 
%                             ||y - mean(y)||_2 
%
%  The FIT of two signals that are the same is 100%. If they differ, the
%  FIT will be lower. When Y and Yest have multiple columns, the FIT is
%  calculated for every column in Y and Yest. The FIT is often used to
%  verify the correctness of a stable model, by comparing the real output
%  with the estimated output of the stable model.
%
%  See also VAF, PEC.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check input arguments
if nargin < 2
    error('FIT requires two input arguments!');
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

% calculate the model fit
v = zeros(size(y,2),1);
for i = 1:size(y,2)
    v(i) = max(diag(100*(1-norm(y(:,i)-y_est(:,i))./norm(y(:,i)-ones(size(y(:,i),1),1)*mean(y(:,i))))),0);
end





