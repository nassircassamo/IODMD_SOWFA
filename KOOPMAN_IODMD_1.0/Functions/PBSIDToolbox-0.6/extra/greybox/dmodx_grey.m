function [X,R1,R2,R3] = dmodx(X,Xd,Xs,n,nd,ns)
%DMODX  Closed-loop LTI system identification using the PBSIDopt method.
%  x=dmodx(X,n) estimates the state sequence x of the identifiable system
%  with order n. The order n can be determined from the singular values
%  given by dordvarx or dordfir. The matrix X is also calculated by
%  dordvarx or dordfir.
%
%  See also: dordfir, dordvarx.
%
%  References:
%    [1] A. Chiuso, G. Picci, ``Consistency Analysis of Certain Closed-Loop
%    Subspace Identification Methods'', Automatica, Special Issue on System
%    Identification,  41(3), pp.377--391, 2005.
%
%    [2] A. Chiuso, ``The role of vector auto regressive modeling in
%    predictor based subspace identification'', Automatica, 43(6), 
%    pp.1034–-1048, 2007.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% check number of arguments
if nargin < 2
    error('DMODX requires at least two input arguments.');
end

% use only n rows
X = X(1:n,:);
Xd = Xd(1:n,:);
Xs = Xs(1:n,:);
[At,Bt,R3] = canoncorr(Xd',Xs');
[Ad,Bd,R1] = canoncorr(X',Xd');
Xd = Bd'*Xd;
Xd = (Ad(:,1:nd)')\Xd(1:nd,:);

[As,Bs,R2] = canoncorr(X',Xs');
Xs = Bs'*Xs;
Xs = (As(:,1:ns)')\Xs(1:ns,:);
[At,Bt,R4] = canoncorr(X',[Xd; Xs]');
X = (At')\(Bt'*[Xd; Xs]);


% Xd = At*Xd;
% Xs = Bt*Xs;
%X = Ad'*flipud(As')*X;
% [U,S,V] = svd([Xd(1:nd,:); Xs(1:ns,:)],'econ');
% figure, semilogy(diag(S),'*');
% X = diag(sqrt(diag(S(1:n,1:n))))*V(:,1:n)';
%Xd = [zeros(nd,length(Xd)); Xd(nd+1:n,:)];
%X = (eye(n) - Xd*pinv(Xd))*X; 
%Xs = [zeros(ns,length(Xs)); Xs(ns+1:n,:)];
%X = (eye(n) - Xs*pinv(Xs))*X; 






