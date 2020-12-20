function v = pec(y,y_est)
%PEC Prediction error cost
%  V = PEC(Y,Yest) computes the percentage of the model fit (FIT) between
%  the two signals Y and Yest. The FIT is calculated as: 
% 
%               v = 1/sqrt(N) * ||y - y_est||_2
%
%  The PEC of two signals that are the same is 0. If they differ, the
%  PEC will be higher. When Y and Yest have multiple columns, the PEC is
%  calculated for every column in Y and Yest. The PEC is often used to
%  verify the correctness of a model, by comparing the real output
%  with the estimated output of the model.
% 
%  See also VAF, FIT.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check input arguments
if nargin < 2
    error('PEC requires two input arguments!');
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

% compute prediction error cost
v = zeros(size(y,2),1);
for i = 1:size(y,2)
    v(i) = norm(y(:,i) - y_est(:,i))./sqrt(N);
end
