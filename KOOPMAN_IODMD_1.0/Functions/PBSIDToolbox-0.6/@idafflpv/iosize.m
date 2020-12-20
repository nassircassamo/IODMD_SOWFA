function [ny,nu] = iosize(sys) 
%IDAFFLPV/IOSIZE Returns I/O size of dynamical systems.
%   [NY,NU] = IOSIZE(SYS)
%   S = IOSIZE(SYS) returns S = [NY NU].

[ny,nu] = size(sys);
if nargout==1
   ny = [ny nu];
end
