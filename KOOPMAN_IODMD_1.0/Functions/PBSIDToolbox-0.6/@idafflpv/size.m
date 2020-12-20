function [Ny,Nu,Ns,Np] = size(sys)
%IDAFFLPV/SIZE Size of affine LPV state-space system

% Number of states
Ns = size(sys.a,1);

% Determine number scheduling parameters
Np = size(sys.a,2)/Ns;

% Number of inputs
Nu = size(sys.b,2)/Np;

% Number of outputs
Ny = size(sys.c,1);

Np = Np - 1;