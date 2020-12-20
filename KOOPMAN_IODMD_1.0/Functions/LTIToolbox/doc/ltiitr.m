
%% LTIITR
% Calculates an LTI state-trajectory

%% Syntax
% |x = ltiitr(A,B,u,w,x0)|

%% Description
% In its most general setting, this function iterates the state equation of
% an linear time-inavriant (LTI) system. It computes the state _x(k)_ for
% _k_ = _1,2,...,N_ satisfying the LTI state equation:
% 
% $$ x(k+1) = Ax(k) + Bx(k) + w(k) $$
% 
% This function is used internally by <dfunlti.html |dfunlti|>, <dac2b.html
% |dac2b|>, <dac2bd.html |dac2bd|>, <dinit.html |dinit|> and <dltisim.html
% |dltisim|>. It is not meant for stand-alone use.

%% Inputs
% |A| is an LTI state-transition matrix of size _n_ x _n_
%% 
% |B| is an LTI input matrix of size _n_ x _m_.
%% 
% |u| is a _N_ x _m_ matrix containing _N_ samples of the _m_ inputs.
%% 
% |w| is a (optional) _N_ x _n_ matrix containing the process noise.
%% 
% |x0| is the (optional) initial state, an _n_ x _1_ vector.
          
%% Outputs
% |x| is the computed state, an _N_ x _n_ matrix.

%% Algorithm
% A direct iteration of the system's state-transition equation is used to
% obtain the state-trajectory for all time-instants.

%% Used By
% <dfunlti.html |dfunlti|>, <dac2b.html |dac2b|>, <dac2bd.html |dac2bd|>,
% <dinit.html |dinit|>, <dltisim.html |dltisim|>

%% See Also
% <dfunlti.html |dfunlti|>, <dac2b.html |dac2b|>, <dac2bd.html |dac2bd|>,
% <dinit.html |dinit|>, <dltisim.html |dltisim|>


