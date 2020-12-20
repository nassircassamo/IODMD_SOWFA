function []=comparereconstruction_pitch(states, statesrebuild,D,dirdmd,x,y,z,Decimate,dirName,filename,initialtime,Inputs)

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
   
   [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    X = length(xx);
    Y = length(yy);
    Z = length(zz);

    Xm_d = Xm_sh/D;
    Ym_d = Ym_sh/D;
    Zm_d = Zm_sh/D;

    
    %% First figure
    fig500= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    set(gcf,'color','w','Position', get(0, 'Screensize'));   
    
    i=initialtime;
    subplot(3,3,1) %DMD RECONSTRUCTION
    flow=statesrebuild(:,i);
    UmeanAbs_sh_u = reshape(flow,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{DMD}  / U_\infty ';
    tt=title(['DMD  \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,2) %TRUE FLOW FIELD FROM SOWFA
    flowsowfa=states(:,i);
    UmeanAbs_sh_u = reshape(flowsowfa,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{SOWFA}  / U_\infty ';
    tt=title([' SOWFA \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,3) %RELATIVE ERROR
    error=abs( abs((states(:,i)-statesrebuild(:,i))) ./ abs((states(:,i))))*100;
    error3d=reshape(error,Y,X,Z);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    err=error3d(:,k,:);
    Usq=squeeze(err);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 100])
    c=colorbar('eastoutside');
    c.Label.String = '[%]';
    set(gca,'fontsize', 14)
    tt=title('Relative deviation');
    tt.FontWeight='normal';
    
    %%
    i=i+20;
    subplot(3,3,4) %DMD RECONSTRUCTION
    flow=statesrebuild(:,i);
    UmeanAbs_sh_u = reshape(flow,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{DMD}  / U_\infty ';
     tt=title(['DMD  \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,5) %TRUE FLOW FIELD FROM SOWFA
    flowsowfa=states(:,i);
    UmeanAbs_sh_u = reshape(flowsowfa,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{SOWFA}  / U_\infty ';
     tt=title([' SOWFA \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,6) %RELATIVE ERROR
    error=abs( abs((states(:,i)-statesrebuild(:,i))) ./ abs((states(:,i))))*100;
    error3d=reshape(error,Y,X,Z);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    err=error3d(:,k,:);
    Usq=squeeze(err);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 100])
    c=colorbar('eastoutside');
    c.Label.String = '[%]';
    set(gca,'fontsize', 14)
    tt=title('Relative deviation');
    tt.FontWeight='normal';
   
   i=i+20;
    subplot(3,3,7) %DMD RECONSTRUCTION
    flow=statesrebuild(:,i);
    UmeanAbs_sh_u = reshape(flow,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{DMD}  / U_\infty ';
     tt=title(['DMD  \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,8) %TRUE FLOW FIELD FROM SOWFA
    flowsowfa=states(:,i);
    UmeanAbs_sh_u = reshape(flowsowfa,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{SOWFA}  / U_\infty ';
     tt=title([' SOWFA \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,9) %RELATIVE ERROR
    error=abs( abs((states(:,i)-statesrebuild(:,i))) ./ abs((states(:,i))))*100;
    error3d=reshape(error,Y,X,Z);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    err=error3d(:,k,:);
    Usq=squeeze(err);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 100])
    c=colorbar('eastoutside');
    c.Label.String = '[%]';
    set(gca,'fontsize', 14)
    tt=title('Relative deviation');
    tt.FontWeight='normal';
    

    shg
    export_fig(fig500,strcat(dirdmd,'/image',strcat(filename,'2')),'-nocrop','-m2'); 
     
    
     %% SECOND FIGURE 
    fig502= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    set(gcf,'color','w','Position', get(0, 'Screensize'));   
    
   i=i+10;
    subplot(3,3,1) %DMD RECONSTRUCTION
    flow=statesrebuild(:,i);
    UmeanAbs_sh_u = reshape(flow,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{DMD}  / U_\infty ';
     tt=title(['DMD  \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,2) %TRUE FLOW FIELD FROM SOWFA
    flowsowfa=states(:,i);
    UmeanAbs_sh_u = reshape(flowsowfa,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{SOWFA}  / U_\infty ';
     tt=title([' SOWFA \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,3) %RELATIVE ERROR
    error=abs( abs((states(:,i)-statesrebuild(:,i))) ./ abs((states(:,i))))*100;
    error3d=reshape(error,Y,X,Z);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    err=error3d(:,k,:);
    Usq=squeeze(err);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 100])
    c=colorbar('eastoutside');
    c.Label.String = '[%]';
    set(gca,'fontsize', 14)
    tt=title('Relative deviation');
    tt.FontWeight='normal';
    
    %%
   i=i+10;
    subplot(3,3,4) %DMD RECONSTRUCTION
    flow=statesrebuild(:,i);
    UmeanAbs_sh_u = reshape(flow,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{DMD}  / U_\infty ';
     tt=title(['DMD  \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,5) %TRUE FLOW FIELD FROM SOWFA
    flowsowfa=states(:,i);
    UmeanAbs_sh_u = reshape(flowsowfa,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{SOWFA}  / U_\infty ';
     tt=title([' SOWFA \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,6) %RELATIVE ERROR
    error=abs( abs((states(:,i)-statesrebuild(:,i))) ./ abs((states(:,i))))*100;
    error3d=reshape(error,Y,X,Z);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    err=error3d(:,k,:);
    Usq=squeeze(err);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 100])
    c=colorbar('eastoutside');
    c.Label.String = '[%]';
    set(gca,'fontsize', 14)
    tt=title('Relative deviation');
    tt.FontWeight='normal';
   
   i=i+10;
    subplot(3,3,7) %DMD RECONSTRUCTION
    flow=statesrebuild(:,i);
    UmeanAbs_sh_u = reshape(flow,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{DMD}  / U_\infty ';
     tt=title(['DMD  \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,8) %TRUE FLOW FIELD FROM SOWFA
    flowsowfa=states(:,i);
    UmeanAbs_sh_u = reshape(flowsowfa,Y,X,Z);
    k=72;
    Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.2])
    set(gca,'fontsize', 14)
    c=colorbar('eastoutside');
    c.Label.String = 'u_{SOWFA}  / U_\infty ';
     tt=title([' SOWFA \theta: ',num2str(Inputs(i)),'º   \deltat: ',num2str((i-initialtime)*2),' s']);
    tt.FontWeight='normal';
    
    subplot(3,3,9) %RELATIVE ERROR
    error=abs( abs((states(:,i)-statesrebuild(:,i))) ./ abs((states(:,i))))*100;
    error3d=reshape(error,Y,X,Z);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    err=error3d(:,k,:);
    Usq=squeeze(err);
    a=pcolor(Ym_shs',Zm_shs',Usq);
    set(a,'ZData',-1+zeros(size(Usq)));
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    hold on
    pbaspect([1 1 1])
    hold on
    [p]=plotrotor2D(0,D,0,115);
    daspect([ 1 1 1])
    caxis([0 100])
    c=colorbar('eastoutside');
    c.Label.String = '[%]';
    set(gca,'fontsize', 14)
    tt=title('Relative deviation');
    tt.FontWeight='normal';
    
    shg
    export_fig(fig502,strcat(dirdmd,'/image',strcat(filename,'2')),'-nocrop','-m2'); 
 

    