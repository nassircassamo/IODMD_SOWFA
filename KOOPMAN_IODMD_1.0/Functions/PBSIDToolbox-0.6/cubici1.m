function r=cubici1(fnew,fold,graddnew,graddold,stepsize)
%CUBICI1 Cubicly interpolates 2 points and gradients to estimate minimum.
%
%	This function uses cubic interpolation and the values of two 
%	points and their gradients in order to estimate the minimum of a 
%	a function along a line.

%	Copyright (c) 1990 by the MathWorks, Inc.
%	Andy Grace 7-9-90.

if fnew==Inf, fnew=1/eps; end
z=3*(fold-fnew)/stepsize+graddold+graddnew;
w=real(sqrt(z*z-graddold*graddnew));
r=stepsize*((z+w-graddold)/(graddnew-graddold+2*w));
         
