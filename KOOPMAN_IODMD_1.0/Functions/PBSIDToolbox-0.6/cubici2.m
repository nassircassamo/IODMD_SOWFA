function [stepmin]=cubici2(graddold,matl,matx)
%CUBICI2 Cubicly interpolates 3 points and 1 gradient.
%
%	This function uses cubic interpolation and
% 	the values of 3  points and one gradient.

%	Copyright (c) 1990 by the MathWorks, Inc.
%	Andy Grace 7-9-90.
abd=[1/3*matx.^3, 0.5*matx.^2, ones(3,1)]\(matl-graddold*matx);
root=real(sqrt(abd(2)^2-4*abd(1)*graddold));
x1=(-abd(2)+root)/(2*abd(1));
if 2*abd(1)*x1+abd(2)>0
    stepmin=x1;
   else
    stepmin=(-abd(2)-root)/(2*abd(1));
end
if stepmin<0|isnan(stepmin)|stepmin==Inf , stepmin=abs(quadi(matx,matl));end
if isnan(stepmin),stepmin=matx(2)/2; end

% fbest=1/3*abd(1)*stepmin^3+0.5*abd(2)*stepmin^2+graddold*stepmin+matl(1);
