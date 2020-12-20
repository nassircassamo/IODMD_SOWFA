function [dirpathfinal]=makefinalframes(D,dirpathfinal,cases)

%% Generate directory in external hard drive

if ~exist(dirpathfinal,'dir') 
    mkdir(dirpathfinal);
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
[nTurbine,time4,dt,nVal,thrustv]       = readTurbineOutputGlobal(cases{n},'rotorVerticalForce');
[nTurbine,time4,dt,nVal,thrusth]       = readTurbineOutputGlobal(cases{n},'rotorHorizontalForce');
[nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
[nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
[nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');

%% Resample data in order to have according dimensions
thrustrs(:,:)       =resample(thrust(end-750*10:1:end,1:end),1,10);
thrustvrs(:,:)      =resample(thrustv(end-750*10:1:end,1:end),1,10);
thrusthrs(:,:)      =resample(thrusth(end-750*10:1:end,1:end),1,10);
azirs(:,:)          =resample(azi(end-750*10:1:end,1:end),1,10);
rotorPowerrs{n}(:,:)=resample(rotorPower{n}(end-750*10:1:end,1:end),1,10);
yawanglers(:,:)     =resample(yawangle(end-750*10:1:end,1:end),1,10);

 % IIR filter design
s=tf('s'); 
omega=2*pi*2;
filter=c2d(omega^2/(s^2+2*0.7*omega*s+omega^2),0.1); 
[a,b]=tfdata(filter);
[b,a] = butter(12,0.08,'low');          

thrustVec{n}=[thrustrs,thrusthrs,thrustvrs]./(sqrt(thrustrs(1)^2+thrustvrs(1)^2+thrusthrs(1)^2));
thrustVec{n}=filtfilt(b,a,thrustVec{n});

%Tranforma��o que n�o percebo completamente bem
thrustVec{n}(:,1)=(thrustVec{n}(:,1)-mean(thrustVec{n}(:,1))).*cosd(5)+(thrustVec{n}(:,5)-mean(thrustVec{n}(:,5)))*sind(5);
thrustVec{n}(:,2)=(thrustVec{n}(:,2)-mean(thrustVec{n}(:,2))).*cosd(5)+(thrustVec{n}(:,6)-mean(thrustVec{n}(:,6)))*sind(5);
thrustVec{n}(:,3)=thrustVec{n}(:,3)-mean(thrustVec{n}(:,3));
thrustVec{n}(:,4)=thrustVec{n}(:,4)-mean(thrustVec{n}(:,4));
thrustVec{n}(:,5)=(thrustVec{n}(:,1)-mean(thrustVec{n}(:,1))).*sind(5)+(thrustVec{n}(:,5)-mean(thrustVec{n}(:,5)))*cosd(5);
thrustVec{n}(:,6)=(thrustVec{n}(:,2)-mean(thrustVec{n}(:,2))).*sind(5)+(thrustVec{n}(:,6)-mean(thrustVec{n}(:,6)))*cosd(5);
dT{n}=mean(diff(time4))*10;
Azi{n}=azirs;

%% Load data realtive to (1) grid spacing and (2) velocity field
load('U_data_complete_vec');

%data from velocity field was pre processed, and only points located every 4
%in the grid were chosen. Therefore, resample x, y and z to get grid were
%velocity field was sampled
[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
[Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx,yy,zz);
X = length(xx);
Y = length(yy);
Z = length(zz);
    
[l,c]=size(QQ_u);

for i=10:1:(c-8)

 %velocity fields QQ are processed such as when they reach this phase
    %velocity fiel QQ is, for each time instant (column) a staking of Y, Z
    %and X has been done
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    
    %curl_ circulation density of the fluid and provides mathematical
    %insights into fluid rotation based on fluid velocity field
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    
    figure1= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    handle = tight_subplot(2,1,[-0.1 -0.1],[0.05 0.1],0.05);
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2);
    set(figure1,'CurrentAxes',handle(1))
    
    %% Plot yaw angle
    rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    hold on
    %%
    
    %normalize meshgrid dimensions by rotor diameter
    Xm_d = Xm_sh/D;
    Ym_d = Ym_sh/D;
    Zm_d = Zm_sh/D;
       
    p = patch(isosurface(Xm_d,Ym_d,Zm_d,UmeanAbs_sh_u,4.75));
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
    
    daspect([1 1 1])
    view(3);
    view(-25,20)
    grid on
    axis tight
    %axis([450/D 2000/D 300/D 700/D 0 300/D]) %defining x, y and z limits based on a rotor diameter sacle
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;
    ax.XTick = 0:0.99:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    grid off
    colormap(F)
    
    %% 3D Arrow plotting for thrust variation (normalized wrt first value)
    hh=mArrow3([500/D 500/D 115/D],[1-1.0*thrustVec{n}((i),1) 500/D-15*thrustVec{n}((i),3) 115/D-15*thrustVec{n}((i),5)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hh.FaceLighting='gouraud';
    % hh=mArrow3([500/D+5D 500/D 115/D],[1-1.0*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),1) 500/D-15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),3) 115/D+15*thrustVec{n}(floor(1/dT{n})*i-floor(1000/dT{n}),5)],'color','yellow','stemWidth',0.05,'facealpha',0.7)
    hhh=mArrow3([500/D+5 500/D 115/D],[6-1.0*thrustVec{n}((i),2) 500/D-15*thrustVec{n}((i),4) 115/D-15*thrustVec{n}((i),6)],'color',[255 196 0]./255,'stemWidth',0.05,'facealpha',0.7);
    hhh.FaceLighting='gouraud';
   
    %% Write wind farm power output
    STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
    text(7.0,4,1.5,STring1,'Fontsize',16,'FontName','Tahoma')
    STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
    text(2.0,3.1,1.5,STring2,'Fontsize',16,'FontName','Tahoma')
    STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
    text(10.0,3.4,1.75,STring3,'Fontsize',16,'FontName','Tahoma')
    
     k=i+250;
     number=k/30;
     integ=floor(number);
     fract=number-integ;
     
     minutos=integ;
     segundos=60*fract;
     

%% LOWER PLOT
    
    set(figure1,'CurrentAxes',handle(2))
    rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
    rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);
    
    %% Define axis aspect and length
    daspect([1 1 1])
    view(3);
    view(-25,20)
    grid off
    axis tight
    
    axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
    ax = gca;

    ax.XTick = 0:0.99:11;
    ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ax.YTickLabel = {'','',''};
    ax.ZTick = 0:0.5:2; 
    ax.ZTickLabel = {''};
    
   
     %% Make slices and plot directly in figure
     
    Uups=9;
    
    VEL=sqrt(UmeanAbs_sh_u.^2+UmeanAbs_sh_v.^2+UmeanAbs_sh_w.^2)/Uups;
    yslice=500/D;
    hold on
    %s3=slice(Xm_d,Ym_d,Zm_d,VEL,[],yslice,[]);
    warning off
    hold on;
    VEL=sqrt(UmeanAbs_sh_u.^2+UmeanAbs_sh_v.^2+UmeanAbs_sh_w.^2)/Uups;
    xslicetb1=[500/D ;500/D+1;500/D+2;500/D+3;500/D+4];
    xslicetb2=[500/D+5 ;500/D+6;500/D+7;500/D+8;500/D+9];
    yslice=500/D;% [500]/D;
    zslice=[115/D];

    s1=slice(Xm_d,Ym_d,Zm_d,VEL,xslicetb1,[],[]);
    s2=slice(Xm_d,Ym_d,Zm_d,VEL,xslicetb2,[],[]);

    shading interp
    colormap(jet(4096))
    

 %% Write plot title
    titlee=suptitle(['First turbine yaw angle of: ', ...
        num2str(round(abs(270-yawanglers(i,1)))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=18;
    
    %%
    warning off
    export_fig(figure1,strcat(dirpathfinal,'/image',num2str(10000+i)),'-nocrop','-m2')
    warning on
    close all

end

