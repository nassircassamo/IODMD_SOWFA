function [stepmin,fbest]=cubici3(fnew,fold,graddnew,graddold,stepsize)
%CUBICI3  Cubicly interpolates 2 points and gradients to find step and min.
%
%	This function uses cubic interpolation and the values of 
%	two points and their gradients in order estimate the minimum of a 
%	a function along a line.

%	Copyright (c) 1990 by the MathWorks, Inc.
%	Andy Grace 7-9-90.
if fnew==Inf, fnew=1/eps; end
amat=[1/3*stepsize^3 , 0.5*stepsize^2; stepsize^2     stepsize];
bmat=[fnew-graddold*stepsize-fold; graddnew-graddold];
abd=amat\bmat;
root=real(sqrt(abd(2)^2-4*abd(1)*graddold));
x1=(-abd(2)+root)/(2*abd(1));
if 2*abd(1)*x1+abd(2)>0
    stepmin=x1;
   else
    stepmin=(-abd(2)-root)/(2*abd(1));;
end
if stepmin<0,  stepmin=-stepmin; end
fbest=1/3*abd(1)*stepmin^3+0.5*abd(2)*stepmin^2+graddold*stepmin+fold;
