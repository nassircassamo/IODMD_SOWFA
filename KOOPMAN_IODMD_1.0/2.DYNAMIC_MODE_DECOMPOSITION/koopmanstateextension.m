function [nonlobs]=koopmanstateextension(u, v, w,rho)

%% Possible extensions
absvel=u.^2+v.^2+w.^2;

steadyu=u(:,100:204);
steadyv=v(:,100:204);
steadyw=w(:,100:204);

 for i=1:size(steadyu,1)
    meanu(i)=mean(steadyu(i,:));
    meanv(i)=mean(steadyv(i,:));
    meanw(i)=mean(steadyw(i,:));
 end
     
    meanu=meanu';
    meanv=meanv';
    meanw=meanw';
    
for l=1:size(u,2)
    detu(:,l)=u(:,l)-meanu;
    detv(:,l)=v(:,l)-meanv;
    detw(:,l)=w(:,l)-meanw;
end

% RANS
    uv=-detu.*detv*rho;
    
    absveldet=detu.^2+detv.^2+detw.^2;
    absveldet=detu.^2+detv.^2+detw.^2;

%% Declare Non Linear Observables for states extension
nonlobs=absvel;