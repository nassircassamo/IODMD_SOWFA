function [us,Du,ys,Dy,zs,Dz] = sigscale(u,y,z)
%SIGSCALE Scaling of Identifications Signals
%  [US,DU,YS,DY]=SIGSCALE(U,Y) scales the input and output signals U and
%  Y, such that the returned signal have the variance VAR(US)=1 and
%  VAR(YS)=1. The scaling matrices DU and DY can be used to convert the
%  amplitudes back to original size Y=DY*YS and U=DU*US. The original
%  state-space system is obtained from the identified system as follows:
%   
%   x(k+1) = A x(k) + B/Du u(k)
%   y(k) = Dy*C x(k) + Dy*D/Du u(k)

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number if input arguments
if nargin < 2
    error('SIGSCALE requires two input arguments.')
end

% check dimensions of inputs
if size(y,2) < size(y,1)
    y = y';
end
if size(u,2) < size(u,1)
    u = u';
end
% N = size(y,2);
% if ~isequal(N,length(u))
%     error('The number of rows of vectors/matrices u and y must be the same.')
% end

Dy = diag(sqrt(var(y,0,2))); % output scaling matrix
ys = diag(1./diag(Dy))*y; % scaling

Du = diag(sqrt(var(u,0,2))); % input scaling matrix
us = diag(1./diag(Du))*u; % scaling

if nargin > 2
    if size(z,2) < size(z,1)
        z = z';
    end
    Dz = diag(sqrt(var(z,0,2))); % output scaling matrix
    zs = diag(1./diag(Dz))*z; % scaling
end
    