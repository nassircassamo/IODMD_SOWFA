function []=wake_vorticity_pitch(dirName,filename,D,DCAT)

%% WAKE PROPAGATION DURING FIRST YAW

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
    
    [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx,yy,zz);
    X = length(xx);
    Y = length(yy);
    Z = length(zz);
    Xm_d = Xm_sh/D;
    Ym_d = Ym_sh/D;
    Zm_d = Zm_sh/D;
    
    fig1= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
   
    
    %% 1
    subplot(5,1,1)
    i=DCAT;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
    STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
    
    %% 2
    subplot(5,1,2)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
    %% 3
    subplot(5,1,3)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
    %% 4
    
    subplot(5,1,4)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
   titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
    %% 5
    subplot(5,1,5)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
    set(gcf, 'Position', get(0, 'Screensize'));
    
    
    shg
    
    %% NEW FIGURE
    fig2= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    
    %% 6
    subplot(5,1,1)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    %% 7
    subplot(5,1,2)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
    %% 8
    subplot(5,1,3)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    %% 9
    
    subplot(5,1,4)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
     %% 10
    
    subplot(5,1,5)
    i=i+10;
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,5.65));
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on
    p.FaceColor = 'yellow';
    p.EdgeColor = 'none';
    p.FaceAlpha =0.1;
    isonormals(Xm_d,Ym_d,Zm_d,absVor,p);  
    hold on;
    p.FaceColor = [88./255 162/255 206/255];
    p.EdgeColor = 'none';
    p.FaceAlpha =0.6;
    p.FaceLighting='gouraud';
    daspect([0.8 1 1])
    view(0,0)
    grid on
    axis tight
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.985:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.85,1.75,1.5,STring1,'Fontsize',12,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(3.0,1.75,1.5,STring2,'Fontsize',12,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(12.0,2.5,1.5,STring3,'Fontsize',12,'FontName','Tahoma')
     STring4={sprintf('Collective pitch angle: %2.0f �',round(abs(pitch(1,i))))};
    text(3.0,1.75,1.7,STring4,'Fontsize',12,'FontName','Tahoma')
    titlee=title({
        [' Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']});
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    titlee.FontName='Tahoma';
    titlee.Position=[ titlee.Position(1)*1.1  titlee.Position(2)-titlee.Position(2)*0.4  titlee.Position(3)];
    
    shg
    
