function [statesrebuildfinal]=rebuild(phi,b,LambdaDiag,mn,X,Xd)

mm1=size(X,2);

dt=2;
omega=log(LambdaDiag)/dt;
%time_dynamics=zeros(mn,mm1);
t=(0:mm1-1)*dt;  %timevetor

%% Reconstruction of the dynamics for wake sterring happens 

for iter=1:length(t)
    time_dynamics(:,iter)=b.*exp(omega*t(iter));
end

statesrebuild=phi*time_dynamics;

%% Get only relevant states 
% (if extension of basic states has been done, it should be discarded)

if isempty(Xd)
    statesrebuildfinal=statesrebuild;
else
    statesrebuildfinal=statesrebuild(size(Xd,1)+1:end,:);
end

