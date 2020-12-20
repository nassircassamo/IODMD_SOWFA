function []=comparereconstruction(states, statesrebuild,D,dirdmd,x,y,z,Decimate,dirName,filename,initialtime,Inputs)

    n=1;
    cases = dirName;
% 
%     [nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
%     [nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
     [nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');
%     
%     azirs(:,:)          =resample(azi(end-750*10:1:end,1:end),1,10);
%     rotorPowerrs{n}(:,:)=resample(rotorPower{n}(end-750*10:1:end,1:end),1,10);
    yawanglers(:,:)     =resample(yawangle(end-750*10:1:end,1:end),1,10);

    Uups=9; %[m/s]
    [xx,yy,zz,X,Y,Z]=retakepoints([],x,y,z,Decimate);
   % [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);
   % X = length(xx);
   % Y = length(yy);
   % Z = length(zz);
    
    %% First figure
    fig500= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    set(gcf,'color','w','Position', get(0, 'Screensize'));   
    Inputs(2,:)=Inputs(1,:);
    Inputs(2,:)=0;
    Inputs=Inputs'+270;
    yawanglers=Inputs;
    %yawanglers=270*ones(size(Inputs,1),size(Inputs,2));%pitch
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    
    i=initialtime;
    subplot(5,2,1)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,2)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+2;
    subplot(5,2,3)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,4)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+13;
    subplot(5,2,5)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,6)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+10;
    subplot(5,2,7)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,8)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+10;
    subplot(5,2,9)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,10)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    [ax4,h3]=suplabel('DMD flow field reconstruction                                                                                               SOWFA flow field'  ,'t');
    set(h3,'FontSize',16)
    export_fig(fig500,strcat(dirdmd,'/image',strcat(filename,'1')),'-nocrop','-m2'); 
 
   
    %% SECOND FIGURE
    fig502= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    set(gcf,'color','w','Position', get(0, 'Screensize'));  
    
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    
    i=i+10;
    subplot(5,2,1)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,2)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+10;
    subplot(5,2,3)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,4)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+10;
    subplot(5,2,5)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,6)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+10;
    subplot(5,2,7)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,8)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=i+10;
    subplot(5,2,9)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,10)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    [ax4,h3]=suplabel('DMD flow field reconstruction                                                                                               SOWFA flow field'  ,'t');
    set(h3,'FontSize',16)
    
    export_fig(fig502,strcat(dirdmd,'/image',strcat(filename,'2')),'-nocrop','-m2'); 
    close all
    
    