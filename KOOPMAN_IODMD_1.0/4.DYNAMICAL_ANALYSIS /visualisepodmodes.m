function []=visualisepodmodes(phi,f, P,x,y,z,Decimate,D,LambdaDiag,damping,method,Xd,dirdmd)

%POD modes are ordered by energy amount, so let us first order them and
%then visualise them 

%% TAKE ONLY MODES RELATIVE TO RELEVANT STATES
[sortedP,I]=sort(P,'descend');
sortedphi=phi(:,I);
sortedf=f(I);
sorteddamping=damping(I);

%% RENAMING FOR PLOTTING
P=sortedP;
phi=sortedphi;
f=sortedf;
damping=sorteddamping;

%% MODE VISUALISAITON 
sos=96600; %size of one state, meaning sie of whole velocity component x y or z

statesused=size(phi,1);
u=phi(1:sos, :);
v=phi(sos+1:end,:);

if statesused < 100000
    %only streamwise u component has been used 
    
    %visualise certain modes
    fig5= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    
    mode=1;
    %mode=34;
    subplot(5,2,1)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=2;
    %mode=36;
    subplot(5,2,2)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=4;
    %mode=38;
    subplot(5,2,3)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=6;
    %mode=40;
    subplot(5,2,4)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=10;
    %mode=42;
    subplot(5,2,5)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=16;
    %mode=44;
    subplot(5,2,6)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=32;
    %mode=46;
    subplot(5,2,7)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=40;
    %mode=48;
    subplot(5,2,8)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=80;
    %mode=50;
    subplot(5,2,9)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=length(f);
    %mode=52;
    subplot(5,2,10)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    export_fig(fig5,strcat(dirdmd,'/image','podmodes'),'-nocrop','-m2');   
   
elseif statesused > 100000
    
    % both streamwise and spanwise velocity components have been used to
    % perform DMD
    %% FOR U
    fig5= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    
    mode=1;
    subplot(5,2,1)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=2;
    subplot(5,2,2)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=4;
    subplot(5,2,3)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=6;
    subplot(5,2,4)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=10;
    subplot(5,2,5)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=20;
    subplot(5,2,6)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=40;
    subplot(5,2,7)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=80;
    subplot(5,2,8)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=120;
    subplot(5,2,9)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=150;
    subplot(5,2,10)
    plotmode(x,y,z,u,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    export_fig(fig5,strcat(dirdmd,'/image','podmodesu'),'-nocrop','-m2'); 
    
    
    %%
    
    fig6= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    
    mode=1;
    subplot(5,2,1)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=2;
    subplot(5,2,2)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=4;
    subplot(5,2,3)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=6;
    subplot(5,2,4)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=10;
    subplot(5,2,5)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=20;
    subplot(5,2,6)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=40;
    subplot(5,2,7)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=80;
    subplot(5,2,8)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=120;
    subplot(5,2,9)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    mode=150;
    subplot(5,2,10)
    plotmode(x,y,z,v,mode,Decimate,D,f,phi,P,LambdaDiag,damping)
    
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    export_fig(fig6,strcat(dirdmd,'/image','podmodesv'),'-nocrop','-m2'); 
end



end
