
%% VAF
% Calculates the Variance Accounted For between two signals.

%% Syntax
% |v = vaf(y,yhat)|

%% Description
% The function |vaf| calculates the Variance Accounted For between two
% signals. The VAF between |y| and |yhat| for the _i_ th component is
% defined as
%   
% $$ \mathrm{VAF}_i = \left(1 - \frac{\mathrm{var}(y_i -
% \hat{y}_i)}{\mathrm{var}(y_i)} \right) \cdot 100 \% $$
% 
% The VAF of two signals that are the same is 100%. If they differ, the
% VAF will be lower. If |y| and |yhat| have multiple columns, the VAF is
% calculated for every column in |y| and |yhat| seperately.
% 
% The VAF is often used to verify the correctness of a model, by comparing
% the real output with the estimated output of the model.

%% Inputs
% |y| is the measured output _y(k)_.
%%
% |yhat| is the estimated output _yhat(k)_.
          
%% Outputs
% |v| is the Variance Accounted For as defined above.

%% Used By
% This is a top-level function that is used directly by the user.


