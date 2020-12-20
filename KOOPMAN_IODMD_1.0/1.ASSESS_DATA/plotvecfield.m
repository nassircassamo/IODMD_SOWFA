function [dirpathvelfield]=plotvecfield(D,dirpathvelfield,cases,pitchmode,flow)

if ~exist(dirpathvelfield,'dir') 
    mkdir(dirpathvelfield);
end

%% Graphical interface characteristics
%blue scale interface where higher deficit correponds to stronger colour
warning off
F=cbrewer('seq', 'Blues', 256); 
F=flipud(F);
F(:,1)=F(:,1)./F(end,1);
F(:,2)=F(:,2)./F(end,2);
F(:,3)=F(:,3)./F(end,3);
warning on

%% Read direct simulation data for turbine operation features
n=1;

[nTurbine,time4,dt,nVal,thrust]        = readTurbineOutputGlobal(cases{n},'rotorAxialForce');
[nTurbine,time3,dt,nVal,powerGenerator] = readTurbineOutputGlobal(cases{n},'generatorPower');
[time2,pitch]                          = readPitchData(strcat(cases{n},'/1000/bladePitch')); %pitch=pitch{:};
[nTurbine,time4,dt,nVal,thrustv]       = readTurbineOutputGlobal(cases{n},'rotorVerticalForce');
[nTurbine,time4,dt,nVal,thrusth]       = readTurbineOutputGlobal(cases{n},'rotorHorizontalForce');
[nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
[nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
[nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');

%% Resample data in order to have according dimensions
dev=1;
thrustrs(:,:)        =resample(thrust(1+dev:1:end,1:end),1,10);
thrustvrs(:,:)       =resample(thrustv(1+dev:1:end,1:end),1,10);
thrusthrs(:,:)       =resample(thrusth(1+dev:1:end,1:end),1,10);
azirs(:,:)           =resample(azi(1+dev:1:end,1:end),1,10);
rotorPowerrs{n}(:,:) =resample(rotorPower{n}(1+dev:1:end,1:end),1,10);
yawanglers(:,:)      =resample(yawangle(1+dev:1:end,1:end),1,10);
pitch                =resampleedgeeffect(pitch{1}(1+dev:1:end,1),10);
powergen             =resample(powerGenerator(1+dev:1:end,1:end),1,10);

% IIR filter design
s=tf('s'); 
omega=2*pi*2;
filter=c2d(omega^2/(s^2+2*0.7*omega*s+omega^2),0.1); 
[a,b]=tfdata(filter);
[b,a] = butter(12,0.08,'low');          

thrustVec{n}=[thrustrs,thrusthrs,thrustvrs]./(sqrt(thrustrs(1)^2+thrustvrs(1)^2+thrusthrs(1)^2));
thrustVec{n}=filtfilt(b,a,thrustVec{n});

%Tranformação que não percebo completamente bem
thrustVec{n}(:,1)=(thrustVec{n}(:,1)-mean(thrustVec{n}(:,1))).*cosd(5)+(thrustVec{n}(:,5)-mean(thrustVec{n}(:,5)))*sind(5);
thrustVec{n}(:,2)=(thrustVec{n}(:,2)-mean(thrustVec{n}(:,2))).*cosd(5)+(thrustVec{n}(:,6)-mean(thrustVec{n}(:,6)))*sind(5);
thrustVec{n}(:,3)=thrustVec{n}(:,3)-mean(thrustVec{n}(:,3));
thrustVec{n}(:,4)=thrustVec{n}(:,4)-mean(thrustVec{n}(:,4));
thrustVec{n}(:,5)=(thrustVec{n}(:,1)-mean(thrustVec{n}(:,1))).*sind(5)+(thrustVec{n}(:,5)-mean(thrustVec{n}(:,5)))*cosd(5);
thrustVec{n}(:,6)=(thrustVec{n}(:,2)-mean(thrustVec{n}(:,2))).*sind(5)+(thrustVec{n}(:,6)-mean(thrustVec{n}(:,6)))*cosd(5);
dT{n}=mean(diff(time4))*10;
Azi{n}=azirs;

%% Load data realtive to (1) grid spacing and (2) velocity field
load(flow);

%data from velocity field was pre processed, and only points located every 4
%in the grid were chosen. Therefore, resample x, y and z to get grid were
%velocity field was sampled
[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
[Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx,yy-500,zz);
X = length(xx);
Y = length(yy);
Z = length(zz);
    
[l,c]=size(QQ_u);

Uups=9; %[m/s]

for i=10:1:(c-8)

    %velocity fields QQ are processed such as when they reach this phase
    %velocity fiel QQ is, for each time instant (column) a staking of Y, Z
    %and X has been done
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    
    %curl_ circulation density of the fluid and provides mathematical
    %insights into fluid rotation based on fluid velocity field
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    
    figure1= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
   
     %% Write time
  
    %% D=0
%     w=1;
%     subplot(2,2,1)
      k=2;
%     Usecu=UmeanAbs_sh_u(:,k,:);
%     Vsecv=UmeanAbs_sh_v(:,k,:);
%     Wsecv=UmeanAbs_sh_w(:,k,:);
%     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
%     Usq=squeeze(Usecu);
%     Vsq=squeeze(Vsecv);
%     Wsq=squeeze(Wsecv);
%     absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2);
%     a=surf(Ym_shs',Zm_shs',absvel);
%     a.FaceAlpha=0.6;
%     shading interp
%     colormap(jet(4096))
%     view(2)
%     hold on
%     quiver(Ym_shs',Zm_shs',Vsq, Wsq);
%     axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
%     xlabel('Distance [m]');
%     ylabel('Distance [m]');
    
    %% D=1
%     w=1;
%     subplot(2,5,2)
      k=k+14;
%     Usecu=UmeanAbs_sh_u(:,k,:);
%     Vsecv=UmeanAbs_sh_v(:,k,:);
%     Wsecv=UmeanAbs_sh_w(:,k,:);
%     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
%     Usq=squeeze(Usecu);
%     Vsq=squeeze(Vsecv);
%     Wsq=squeeze(Wsecv);
%     absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2);
%     a=surf(Ym_shs',Zm_shs',absvel);
%     a.FaceAlpha=0.6;
%     shading interp
%     colormap(jet(4096))
%     view(2)
%     hold on
%     quiver(Ym_shs',Zm_shs',Vsq, Wsq);
%     axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
%     xlabel('Distance [m]');
%     ylabel('Distance [m]');
    
    %% D=2
    subplot(2,3,1)
    k=k+14;
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
    caxis([0 1.4])
    set(gca,'fontsize', 14)

       

    %% D=3
    subplot(2,3,2)
    k=k+14;
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)));
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
    title('x/D = 3');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.4])
    set(gca,'fontsize', 14)
    
    %% D=4
    subplot(2,3,3)
    k=k+14;
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)));
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
    title('x/D = 4');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.4])
    set(gca,'fontsize', 14)
    
    %% D=5
    subplot(2,3,4)
    k=k+14;
    Usecu=UmeanAbs_sh_u(:,k,:);
    Vsecv=UmeanAbs_sh_v(:,k,:);
    Wsecv=UmeanAbs_sh_w(:,k,:);
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
    absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
    a=pcolor(Ym_shs',Zm_shs',absvel);
    set(a,'ZData',-1+zeros(size(absvel)));
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
    title('x/D = 5');
    [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
    daspect([ 1 1 1])
    caxis([0 1.4])
    set(gca,'fontsize', 14)
    
    %% D=6
     subplot(2,3,5)
     k=k+14;
     Usecu=UmeanAbs_sh_u(:,k,:);
     Vsecv=UmeanAbs_sh_v(:,k,:);
     Wsecv=UmeanAbs_sh_w(:,k,:);
     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
     Usq=squeeze(Usecu);
     Vsq=squeeze(Vsecv);
     Wsq=squeeze(Wsecv);
     absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
     a=pcolor(Ym_shs',Zm_shs',absvel);
     set(a,'ZData',-1+zeros(size(absvel)));
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
     title('x/D = 6');
     [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
     daspect([ 1 1 1])
     caxis([0 1.4])
     set(gca,'fontsize', 14)
     
    %% D=7
     hp4 = get(subplot(2,3,6),'Position');
     
    
     w=1;
     subplot(2,3,6)
     k=k+14;
     Usecu=UmeanAbs_sh_u(:,k,:);
     Vsecv=UmeanAbs_sh_v(:,k,:);
     Wsecv=UmeanAbs_sh_w(:,k,:);
     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
     Usq=squeeze(Usecu);
     Vsq=squeeze(Vsecv);
     Wsq=squeeze(Wsecv);
     absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
     a=pcolor(Ym_shs',Zm_shs',absvel);
     set(a,'ZData',-1+zeros(size(absvel)));
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
     title('x/D = 7');
     [p]=plotrotor2D(yawanglers(i,1)+90,D,0,115);
     daspect([ 1 1 1])
     caxis([0 1.4])
     set(gca,'fontsize', 14)
     hold off
    
    %% D
%     subplot(2,5,9)
%     k=k+14;
%     Usecu=UmeanAbs_sh_u(:,k,:);
%     Vsecv=UmeanAbs_sh_v(:,k,:);
%     Wsecv=UmeanAbs_sh_w(:,k,:);
%     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
%     Usq=squeeze(Usecu);
%     Vsq=squeeze(Vsecv);
%     Wsq=squeeze(Wsecv);
%     absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2);
%     a=surf(Ym_shs',Zm_shs',absvel);
%     a.FaceAlpha=0.6;
%     shading interp
%     colormap(jet(4096))
%     view(2)
%     hold on
%     quiver(Ym_shs',Zm_shs',Vsq, Wsq);
%     axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
%     xlabel('Distance [m]');
%     ylabel('Distance [m]');
    
    %% D
%     w=1;
%     subplot(2,5,10)
%     k=k+14;
%     Usecu=UmeanAbs_sh_u(:,k,:);
%     Vsecv=UmeanAbs_sh_v(:,k,:);
%     Wsecv=UmeanAbs_sh_w(:,k,:);
%     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
%     Usq=squeeze(Usecu);
%     Vsq=squeeze(Vsecv);
%     Wsq=squeeze(Wsecv);
%     absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2);
%     a=surf(Ym_shs',Zm_shs',absvel);
%     a.FaceAlpha=0.6;
%     shading interp
%     colormap(jet(4096))
%     view(2)
%     hold on
%     quiver(Ym_shs',Zm_shs',Vsq, Wsq);
%     axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
%     xlabel('Distance [m]');
%     ylabel('Distance [m]');

%%  
     k=i;
     number=k/30;
     integ=floor(number);
     fract=number-integ;
     
     minutos=integ;
     segundos=60*fract;
    
     %% Write plot title
     if pitchmode==0
           [ax4,h3]=suplabel(['First turbine yaw angle of: ', ...
                num2str(round(abs(270-yawanglers(i,1)))), ' degrees. Time: ',num2str(minutos)...
                ,' minutes and ',num2str(segundos),' seconds'],'t');
     elseif pitchmode==1
         [ax4,h3]=suplabel(['Collective pitch angle of first turbine: ', ...
                num2str(round(abs(pitch(i)))), ' degrees. Time: ',num2str(minutos)...
                ,' minutes and ',num2str(segundos),' seconds'],'t');
     end
     set(h3,'FontSize',18)
     set(h3,'FontWeight','normal')
     
     %% COLOR BAR
    
    c=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.025  hp4(2)+hp4(3)*3.25]);
    c.Label.String = '|| (u,v,w) || / U_h ';
  
 %%
    warning off
    export_fig(figure1,strcat(dirpathvelfield,'/image',num2str(10000+i)),'-nocrop','-m2')
    warning on
    close all
    
end

    