    function []=hubheightcut(dirName,filename,D,DCAT)

    warning off
    F=cbrewer('seq', 'Blues', 256); 
    F=flipud(F);
    F(:,1)=F(:,1)./F(end,1);
    F(:,2)=F(:,2)./F(end,2);
    F(:,3)=F(:,3)./F(end,3);
    warning on

     n=1;
    %cases = {'steps_yaw'};
    cases=dirName;

    [nTurbine,time4,dt,nVal,thrust]        = readTurbineOutputGlobal(cases{n},'rotorAxialForce');
    [nTurbine,time4,dt,nVal,thrustv]       = readTurbineOutputGlobal(cases{n},'rotorVerticalForce');
    [nTurbine,time4,dt,nVal,thrusth]       = readTurbineOutputGlobal(cases{n},'rotorHorizontalForce');
    [nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
    [nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
    [nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');

   thrustrs(:,:)       =resample(thrust,1,10);
    thrustvrs(:,:)      =resample(thrustv,1,10);
    thrusthrs(:,:)      =resample(thrusth,1,10);
    azirs(:,:)          =resample(azi,1,10);
    rotorPowerrs{n}(:,:)=resample(rotorPower{n},1,10);
    yawanglers(:,:)     =resample(yawangle,1,10);

    s=tf('s'); 
    omega=2*pi*2;
    filter=c2d(omega^2/(s^2+2*0.7*omega*s+omega^2),0.1); 
    [a,b]=tfdata(filter);
    [b,a] = butter(12,0.08,'low');          

    thrustVec{n}=[thrustrs,thrusthrs,thrustvrs]./(sqrt(thrustrs(1)^2+thrustvrs(1)^2+thrusthrs(1)^2));
    thrustVec{n}=filtfilt(b,a,thrustVec{n});

    thrustVec{n}(:,1)=(thrustVec{n}(:,1)-mean(thrustVec{n}(:,1))).*cosd(5)+(thrustVec{n}(:,5)-mean(thrustVec{n}(:,5)))*sind(5);
    thrustVec{n}(:,2)=(thrustVec{n}(:,2)-mean(thrustVec{n}(:,2))).*cosd(5)+(thrustVec{n}(:,6)-mean(thrustVec{n}(:,6)))*sind(5);
    thrustVec{n}(:,3)=thrustVec{n}(:,3)-mean(thrustVec{n}(:,3));
    thrustVec{n}(:,4)=thrustVec{n}(:,4)-mean(thrustVec{n}(:,4));
    thrustVec{n}(:,5)=(thrustVec{n}(:,1)-mean(thrustVec{n}(:,1))).*sind(5)+(thrustVec{n}(:,5)-mean(thrustVec{n}(:,5)))*cosd(5);
    thrustVec{n}(:,6)=(thrustVec{n}(:,2)-mean(thrustVec{n}(:,2))).*sind(5)+(thrustVec{n}(:,6)-mean(thrustVec{n}(:,6)))*cosd(5);
    dT{n}=mean(diff(time4))*10;
    Azi{n}=azirs;

    load(filename);
    Uups=9; %[m/s]
    
    [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);
    X = length(xx);
    Y = length(yy);
    Z = length(zz);
    Xm_d = Xm_sh/D;
    Ym_d = Ym_sh/D;
    Zm_d = Zm_sh/D;
    
    fig1= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    
    %% I=1
    subplot(5,2,1)
    i=DCAT;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
    
    subplot(5,2,2)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
     %% I=2
    subplot(5,2,3)
    i=i+3;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,4)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
     %% I=3
    subplot(5,2,5)
    i=i+17;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,6)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
    
     %% I=4
    subplot(5,2,7)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,8)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
     %% I=1
    subplot(5,2,9)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,10)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
    shg
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig1= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    
    %% I=1
    subplot(5,2,1)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,2)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
      
     %% I=2
    subplot(5,2,3)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,4)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
     %% I=3
    subplot(5,2,5)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,6)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
     %% I=4
    subplot(5,2,7)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,8)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
     %% I=1
    subplot(5,2,9)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')
      
    subplot(5,2,10)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['First turbine yaw angle of: -', ...
    num2str(round(270-yawanglers(i,1))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'v / U_\infty';
    set(gca,'fontsize', 14) 
    caxis([-0.6 0.5])
    hold off 
    set(gca,'Ydir','reverse')
    shg
    
    
    