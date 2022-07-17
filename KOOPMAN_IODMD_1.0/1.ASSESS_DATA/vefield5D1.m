    function []=vefield5D1(dirName,filename,D,DCAT)
    
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
    [time2,pitch] = readPitchData(strcat(cases{n},'/1000/bladePitch')); %pitch=pitch{:};
    [nTurbine,time4,dt,nVal,thrustv]       = readTurbineOutputGlobal(cases{n},'rotorVerticalForce');
    [nTurbine,time4,dt,nVal,thrusth]       = readTurbineOutputGlobal(cases{n},'rotorHorizontalForce');
    [nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
    [nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
    [nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');
    
    pitch=resampleedgeeffect(pitch{1}(:,1),10);
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
    
    
    subplot(3,4,1)
    k=72;
    i=DCAT;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 2
    subplot(3,4,2)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 3
    subplot(3,4,3)
    k=72;
    i=i+10;;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 4
    subplot(3,4,4)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 5
    subplot(3,4,5)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 6
    subplot(3,4,6)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 7
    subplot(3,4,7)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 8
    subplot(3,4,8)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 9
    subplot(3,4,9)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 10
    subplot(3,4,10)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 11
    subplot(3,4,11)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% 12
    hp4 = get(subplot(3,4,12),'Position');
    subplot(3,4,12)
    k=72;
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)))
    a.FaceAlpha=0.9;
    shading interp
    colormap(jet(4096))
    view(2)
    hold on
    quiver(Ym_shs',Zm_shs',Vsq, Wsq,'color',[0 0 0]);
    axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
    xlabel('Distance [m]');
    ylabel('Distance [m]');
    pbaspect([1 1 1])
    title('x/D = 2');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.6])
    set(gca,'fontsize', 12)
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title([num2str(minutos),' min. and ',num2str(segundos),' sec.',' Power ', num2str(round(rotorPowerrs{1}(i,2)/1e6,2)),' MW' ]);
    titlee.FontSize=16;
    titlee.FontWeight='normal';
    set(gca,'fontname','times','fontsize',16)  % Set it to times
    
    %% color bar
    
    c=colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.025  hp4(2)+hp4(3)*4.5]);
    c.Label.String = '|| (u,v,w) || / U_h ';
    
    