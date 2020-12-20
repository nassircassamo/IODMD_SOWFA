function [a,b,c,d,k,x0,Ts] = ssdata(sys)
%IDAFFLPV/SSDATA  Returns state-space matrices for IDAFFLPV models.
%   [A,B,C,D,K,x0,Ts] = SSDATA(SYS) returns the A,B,C,D matrices of the
%   state-space model SYS.  If SYS is not a state-space model, 
%   it is first converted to the state-space representation.

if nargout < 4
    % Extract data
    [a,b,c,d] = getABCD(sys);
else
    % Extract data
    [a,b,c,d,k] = getABCDK(sys); 
end
   
% Sample time
if nsys==0
   Ts = 0;
else
   Ts = sys.Ts;
end
x0 = sys.x0;
